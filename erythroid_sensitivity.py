"""Is the erythroid tracking claim driven by the globins alone?

The manuscript scores erythroid content from a small marker panel and reports that
every one of the top 100 differentially expressed genes in GSE33615 tracks it. This
script rebuilds that score from subsets of the panel, so a reader can see whether the
result depends on any single marker or on the globins as a group.

Reuses the loader, seed and background of top100_erythroid.py so the numbers are
directly comparable.
"""
import os
os.makedirs('results', exist_ok=True)
import gzip, itertools
import numpy as np, pandas as pd
from scipy import stats

SM = 'GSE33615_series_matrix.txt.gz'
SOFT = 'GSE33615_family.soft.gz'
RNG_SEED = 20260817

print('loading sample metadata ...', flush=True)
meta = {}
with gzip.open(SM, 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!Sample_'):
            k = line.split('\t')[0][1:]
            v = [x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            meta.setdefault(k, []).append(v)
gsm = meta['Sample_geo_accession'][0]
title = meta['Sample_title'][0]
chars = meta['Sample_characteristics_ch1']

def field(prefix):
    for block in chars:
        if any(str(x).lower().startswith(prefix) for x in block):
            return [str(x).split(':', 1)[1].strip() if ':' in str(x) else '' for x in block]
    return [''] * len(gsm)

samples = pd.DataFrame({'gsm': gsm, 'title': title, 'disease': field('disease state')})
samples['group'] = np.where(samples['disease'].str.contains('ATL', case=False, na=False), 'ATL',
                    np.where(samples['disease'].str.contains('health|normal', case=False, na=False),
                             'control', 'UNASSIGNED'))
m = samples['group'] == 'UNASSIGNED'
samples.loc[m, 'group'] = np.where(samples.loc[m, 'title'].str.contains('ATL', case=False),
                                   'ATL', 'control')
print(samples['group'].value_counts().to_string(), flush=True)

print('loading expression matrix ...', flush=True)
with gzip.open(SM, 'rt', errors='replace') as f:
    for i, line in enumerate(f):
        if line.startswith('!series_matrix_table_begin'):
            skip = i + 1
            break
X = pd.read_csv(SM, sep='\t', skiprows=skip, compression='gzip',
                index_col=0, comment='!', low_memory=False)
X = X.apply(pd.to_numeric, errors='coerce').loc[:, samples['gsm']]

print('loading platform annotation ...', flush=True)
rows = []
with gzip.open(SOFT, 'rt', errors='replace') as f:
    inplat = False
    for line in f:
        if line.startswith('!platform_table_begin'):
            cols = next(f).rstrip('\n').split('\t'); inplat = True; continue
        if line.startswith('!platform_table_end'):
            break
        if inplat:
            rows.append(line.rstrip('\n').split('\t'))
ann = pd.DataFrame(rows, columns=cols)[['ID', 'GENE_SYMBOL', 'CONTROL_TYPE']]
ann['ID'] = pd.to_numeric(ann['ID'], errors='coerce')
ann = ann.dropna(subset=['ID']).set_index('ID')
ann.index = ann.index.astype(int)

X = np.log2(X.clip(lower=1))
ranks = X.rank(method='first')
mean_sorted = np.sort(X.values, axis=0).mean(axis=1)
Xq = pd.DataFrame(np.interp(ranks.values, np.arange(1, X.shape[0] + 1), mean_sorted),
                  index=X.index, columns=X.columns)
sym = ann.reindex(Xq.index)['GENE_SYMBOL'].fillna('')
ctrl = ann.reindex(Xq.index)['CONTROL_TYPE'].fillna('')
keep = (sym != '') & (ctrl.str.lower().isin(['', 'false']))
Xq, sym = Xq[keep], sym[keep]
order = Xq.mean(axis=1).sort_values(ascending=False).index
Xq, sym = Xq.loc[order], sym.loc[order]
dup = ~sym.duplicated()
G = Xq[dup]; G.index = sym[dup]
grp = samples.set_index('gsm').loc[G.columns, 'group']
is_atl = (grp == 'ATL').values
print(f'gene-level matrix: {G.shape[0]} genes x {G.shape[1]} samples\n', flush=True)

def cohens_d(row):
    a, b = row[is_atl], row[~is_atl]
    sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / sd if sd > 0 else 0.0

print('ranking genes by effect size ...', flush=True)
d = G.apply(cohens_d, axis=1)
ranked = d.abs().sort_values(ascending=False)
top100 = list(ranked.index[:100]); top500 = list(ranked.index[:500])

FULL = [g for g in ['HBB', 'HBA1', 'HBA2', 'ALAS2', 'SLC4A1', 'AHSP'] if g in G.index]
GLOBIN = [g for g in FULL if g.startswith('HB')]
NONGLOBIN = [g for g in FULL if not g.startswith('HB')]
print('panel present on this platform:', FULL)
print('globins:', GLOBIN, ' non-globins:', NONGLOBIN, '\n', flush=True)

rng = np.random.default_rng(RNG_SEED)
bg_idx = rng.choice(G.shape[0], 4000, replace=False)
BG = G.iloc[bg_idx]

def score(panel):
    return np.asarray(pd.DataFrame(stats.zscore(G.loc[panel].values, axis=1),
                                   columns=G.columns).mean())

def evaluate(label, panel, out):
    ery = score(panel)
    bg = np.array([stats.spearmanr(BG.iloc[i].values, ery)[0] for i in range(len(BG))])
    bgf = (np.abs(bg) > 0.3).mean()
    row = {'panel': label, 'n_markers': len(panel), 'markers': ' '.join(panel),
           'background_frac': round(bgf * 100, 1)}
    for name, lst in [('top100', top100), ('top500', top500)]:
        rho = np.array([stats.spearmanr(G.loc[g].values, ery)[0] for g in lst])
        row[f'{name}_frac_above_0.3'] = round((np.abs(rho) > 0.3).mean() * 100, 1)
        row[f'{name}_median_abs_rho'] = round(float(np.median(np.abs(rho))), 3)
    print(f'  {label:26s} n={len(panel)}  top100 {row["top100_frac_above_0.3"]:5.1f}%  '
          f'top500 {row["top500_frac_above_0.3"]:5.1f}%  background {row["background_frac"]:5.1f}%  '
          f'median |rho| top100 {row["top100_median_abs_rho"]:.3f}', flush=True)
    out.append(row)

out = []
print('PANEL SENSITIVITY: percent of top differentially expressed genes with |rho| > 0.3')
print('-' * 74, flush=True)
evaluate('full panel (as published)', FULL, out)
for g in FULL:
    evaluate(f'leave out {g}', [x for x in FULL if x != g], out)
evaluate('globins only', GLOBIN, out)
evaluate('non-globins only', NONGLOBIN, out)
for g in FULL:
    evaluate(f'{g} alone', [g], out)

R = pd.DataFrame(out)
R.to_csv('results/results_erythroid_panel_sensitivity.csv', index=False)
print('\nwrote results_erythroid_panel_sensitivity.csv')

# correlation among the markers themselves
print('\nMARKER CROSS-CORRELATION (Spearman, across samples)')
C = pd.DataFrame(index=FULL, columns=FULL, dtype=float)
for a, b in itertools.product(FULL, FULL):
    C.loc[a, b] = stats.spearmanr(G.loc[a].values, G.loc[b].values)[0]
print(C.round(2).to_string())
C.round(4).to_csv('results/results_erythroid_marker_correlation.csv')
print('\nwrote results_erythroid_marker_correlation.csv')
