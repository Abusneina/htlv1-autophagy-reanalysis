"""Do the top differentially expressed genes in GSE33615, ranked by effect size
alone and without regard to pathway membership, track erythroid content the same
way the published Wang & Iha gene lists do (see test_wang_iha.py)?

Post hoc: this script was written after the frozen protocol, to confirm the
"every one of the top 100 differentially expressed genes tracks erythroid
content" claim made in Sections 3.3 and 4.2 of the manuscript. It reuses the
same GSE33615 loader as run_gse33615.py, the same erythroid marker panel and
the same 4,000-gene random background under the same seed as test_wang_iha.py,
so its output is directly comparable to that script's.
"""
import gzip
import numpy as np, pandas as pd
from scipy import stats

SM = 'data/GSE33615_series_matrix.txt.gz'
SOFT = 'data/GSE33615_family.soft.gz'
RNG_SEED = 20260817

# ---------- sample metadata (matches run_gse33615.py) ----------
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


disease = field('disease state')
samples = pd.DataFrame({'gsm': gsm, 'title': title, 'disease': disease})
samples['group'] = np.where(samples['disease'].str.contains('ATL', case=False, na=False), 'ATL',
                    np.where(samples['disease'].str.contains('health|normal', case=False, na=False),
                             'control', 'UNASSIGNED'))
mask = samples['group'] == 'UNASSIGNED'
samples.loc[mask, 'group'] = np.where(
    samples.loc[mask, 'title'].str.contains('ATL', case=False), 'ATL', 'control')
print('== sample assignment ==')
print(samples['group'].value_counts().to_string())

# ---------- expression matrix (matches run_gse33615.py) ----------
with gzip.open(SM, 'rt', errors='replace') as f:
    for i, line in enumerate(f):
        if line.startswith('!series_matrix_table_begin'):
            skip = i + 1
            break
X = pd.read_csv(SM, sep='\t', skiprows=skip, compression='gzip',
                index_col=0, comment='!', low_memory=False)
X = X.apply(pd.to_numeric, errors='coerce')
X = X.loc[:, samples['gsm']]
print(f'raw matrix: {X.shape[0]} probes x {X.shape[1]} samples')

rows = []
with gzip.open(SOFT, 'rt', errors='replace') as f:
    inplat = False
    for line in f:
        if line.startswith('!platform_table_begin'):
            cols = next(f).rstrip('\n').split('\t')
            inplat = True
            continue
        if line.startswith('!platform_table_end'):
            break
        if inplat:
            rows.append(line.rstrip('\n').split('\t'))
ann = pd.DataFrame(rows, columns=cols)
ann = ann[['ID', 'GENE_SYMBOL', 'CONTROL_TYPE']]
ann['ID'] = pd.to_numeric(ann['ID'], errors='coerce')
ann = ann.dropna(subset=['ID']).set_index('ID')
ann.index = ann.index.astype(int)

X = np.log2(X.clip(lower=1))
ranks = X.rank(method='first')
mean_sorted = np.sort(X.values, axis=0).mean(axis=1)
Xq = pd.DataFrame(
    np.interp(ranks.values, np.arange(1, X.shape[0] + 1), mean_sorted),
    index=X.index, columns=X.columns)

sym = ann.reindex(Xq.index)['GENE_SYMBOL'].fillna('')
ctrl = ann.reindex(Xq.index)['CONTROL_TYPE'].fillna('')
keep = (sym != '') & (ctrl.str.lower().isin(['', 'false']))
Xq, sym = Xq[keep], sym[keep]
order = Xq.mean(axis=1).sort_values(ascending=False).index
Xq, sym = Xq.loc[order], sym.loc[order]
dup = ~sym.duplicated()
G = Xq[dup]
G.index = sym[dup]
grp = samples.set_index('gsm').loc[G.columns, 'group']
is_atl = (grp == 'ATL').values
print(f'gene-level matrix: {G.shape[0]} genes x {G.shape[1]} samples\n')

# ---------- rank genes by effect size ----------
def cohens_d(row):
    a, b = row[is_atl], row[~is_atl]
    sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / sd if sd > 0 else 0.0


d = G.apply(cohens_d, axis=1)
ranked = d.abs().sort_values(ascending=False)
top100 = ranked.index[:100]
top500 = ranked.index[:500]
print(f"|Cohen's d| at rank 100: {ranked.iloc[99]:.2f}   at rank 500: {ranked.iloc[499]:.2f}\n")

# ---------- erythroid tracking, against the test_wang_iha.py background ----------
ery_genes = [g for g in ['HBB', 'HBA1', 'HBA2', 'ALAS2', 'SLC4A1', 'AHSP'] if g in G.index]
print('erythroid markers used:', ery_genes)
ery = np.asarray(pd.DataFrame(stats.zscore(G.loc[ery_genes].values, axis=1),
                              columns=G.columns).mean())

rng = np.random.default_rng(RNG_SEED)
bg_idx = rng.choice(G.shape[0], 4000, replace=False)
bg_ery = np.array([stats.spearmanr(G.iloc[i].values, ery)[0] for i in bg_idx])
bg_frac = (np.abs(bg_ery) > 0.3).mean()


def report(name, gene_list):
    rho = np.array([stats.spearmanr(G.loc[g].values, ery)[0] for g in gene_list])
    frac = (np.abs(rho) > 0.3).mean()
    u, p = stats.mannwhitneyu(np.abs(rho), np.abs(bg_ery))
    print(f'{name} ({len(gene_list)} genes): median |rho| with erythroid content '
          f'{np.median(np.abs(rho)):.3f} vs background {np.median(np.abs(bg_ery)):.3f}; '
          f'|rho|>0.3 in {frac*100:.0f}% vs {bg_frac*100:.0f}% of background; P={p:.2e}')
    return rho


print()
rho100 = report('Top 100 by |Cohen\'s d|', list(top100))
rho500 = report('Top 500 by |Cohen\'s d|', list(top500))

# ---------- write results ----------
pd.DataFrame({'gene': list(top100), 'cohens_d': d.loc[top100].values,
             'rho_erythroid': rho100}).sort_values(
    'cohens_d', key=abs, ascending=False).to_csv(
    'results_top100_erythroid.csv', index=False)
pd.DataFrame({'gene': list(top500), 'cohens_d': d.loc[top500].values,
             'rho_erythroid': rho500}).sort_values(
    'cohens_d', key=abs, ascending=False).to_csv(
    'results_top500_erythroid.csv', index=False)
print('\nwrote results_top100_erythroid.csv and results_top500_erythroid.csv')

lines = ['# Top differentially expressed genes and erythroid content (post hoc)', '',
         f"Ranking statistic: Cohen's d between ATL and control, {G.shape[0]} genes total.",
         f'Erythroid markers: {", ".join(ery_genes)}.',
         f'Background: 4,000 genes, seed {RNG_SEED} (same as test_wang_iha.py).', '']
for name, gene_list, rho in [('Top 100', list(top100), rho100), ('Top 500', list(top500), rho500)]:
    frac = (np.abs(rho) > 0.3).mean()
    lines.append(f'- {name}: {frac*100:.1f}% of genes have |rho| with erythroid content above 0.3 '
                 f'(median |rho| {np.median(np.abs(rho)):.3f})')
lines += ['', 'Per-gene correlations: results_top100_erythroid.csv and results_top500_erythroid.csv.']
open('top100_erythroid_summary.md', 'w').write('\n'.join(lines) + '\n')
print('wrote top100_erythroid_summary.md')
