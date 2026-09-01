#!/usr/bin/env python3
"""
HTLV-1 / autophagy-lysosome program in ATL
Arm 1: expression cohort, module-level test of Model A vs Model B
Frozen parameters from the 2026-08-17 protocol.
"""
import gzip, re, sys
import numpy as np, pandas as pd
from scipy import stats

FDR = 0.05
MIN_COVERAGE = 0.60
RNG = np.random.default_rng(20260817)

SM = '/mnt/user-data/uploads/GSE33615_series_matrix_txt.gz'
SOFT = '/mnt/user-data/uploads/GSE33615_family_soft.gz'
MODS = '/home/claude/htlv/code/modules_frozen_v1.tsv'

# ---------- sample metadata ----------
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
subtype = field('atl subtype')
sex = field('gender')

samples = pd.DataFrame({'gsm': gsm, 'title': title, 'disease': disease,
                        'subtype': subtype, 'sex': sex})
samples['group'] = np.where(samples['disease'].str.contains('ATL', case=False, na=False), 'ATL',
                    np.where(samples['disease'].str.contains('health|normal', case=False, na=False),
                             'control', 'UNASSIGNED'))
# fall back to the title when the disease field is blank
mask = samples['group'] == 'UNASSIGNED'
samples.loc[mask, 'group'] = np.where(
    samples.loc[mask, 'title'].str.contains('ATL', case=False), 'ATL', 'control')

print('== sample assignment ==')
print(samples['group'].value_counts().to_string())
print(samples.groupby(['group', 'subtype']).size().to_string())
print()

# ---------- expression matrix ----------
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

# ---------- probe annotation ----------
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
ann = pd.DataFrame(rows, columns=cols)
ann = ann[['ID', 'GENE_SYMBOL', 'CONTROL_TYPE']]
ann['ID'] = pd.to_numeric(ann['ID'], errors='coerce')
ann = ann.dropna(subset=['ID']).set_index('ID')
ann.index = ann.index.astype(int)
print(f'platform annotation: {ann.shape[0]} probes, '
      f'{ann["GENE_SYMBOL"].replace("", np.nan).notna().sum()} with a symbol')

# ---------- normalize, collapse ----------
X = np.log2(X.clip(lower=1))
# quantile normalization between arrays
ranks = X.rank(method='first')
mean_sorted = np.sort(X.values, axis=0).mean(axis=1)
Xq = pd.DataFrame(
    np.interp(ranks.values, np.arange(1, X.shape[0] + 1), mean_sorted),
    index=X.index, columns=X.columns)

sym = ann.reindex(Xq.index)['GENE_SYMBOL'].fillna('')
ctrl = ann.reindex(Xq.index)['CONTROL_TYPE'].fillna('')
keep = (sym != '') & (ctrl.str.lower().isin(['', 'false']))
Xq, sym = Xq[keep], sym[keep]
order = Xq.mean(axis=1).sort_values(ascending=False).index      # max mean intensity
Xq, sym = Xq.loc[order], sym.loc[order]
dup = ~sym.duplicated()
G = Xq[dup]; G.index = sym[dup]
print(f'gene-level matrix: {G.shape[0]} genes x {G.shape[1]} samples\n')

# ---------- modules ----------
mt = pd.read_csv(MODS, sep='\t')
mods = {m: sorted(set(g['gene'])) for m, g in mt.groupby('module')}
axis_genes = mods.pop('axis')

cov = {}
for m, genes in list(mods.items()) + [('axis', axis_genes)]:
    present = [g for g in genes if g in G.index]
    cov[m] = (len(present), len(genes), len(present) / len(genes))
print('== module coverage on GPL4133 ==')
for m, (p, t, f) in cov.items():
    flag = '' if f >= MIN_COVERAGE else '   << below floor'
    print(f'  {m:11s} {p:2d}/{t:2d}  {f:.2f}{flag}')
missing = {m: [g for g in genes if g not in G.index]
           for m, genes in list(mods.items()) + [('axis', axis_genes)]}
print('  missing:', {m: v for m, v in missing.items() if v})
print()

# ---------- ssGSEA (Barbie et al.) ----------
def ssgsea(expr, gene_sets, alpha=0.25):
    genes = expr.index.to_numpy()
    out = {}
    R = expr.rank(axis=0, method='average').values
    N = expr.shape[0]
    for name, gs in gene_sets.items():
        idx = np.isin(genes, gs)
        if idx.sum() < 5:
            continue
        scores = []
        for j in range(expr.shape[1]):
            order = np.argsort(-R[:, j])
            inset = idx[order]
            w = (R[order, j] ** alpha) * inset
            cs_in = np.cumsum(w) / w.sum()
            cs_out = np.cumsum(~inset) / (N - inset.sum())
            scores.append(np.sum(cs_in - cs_out))
        out[name] = scores
    S = pd.DataFrame(out, index=expr.columns).T
    return (S - S.values.min()) / (S.values.max() - S.values.min())

S = ssgsea(G, mods)
grp = samples.set_index('gsm').loc[S.columns, 'group']
a, b = S.loc[:, (grp == 'ATL').values], S.loc[:, (grp == 'control').values]

res = []
for m in S.index:
    t, p = stats.mannwhitneyu(a.loc[m], b.loc[m], alternative='two-sided')
    d = (a.loc[m].mean() - b.loc[m].mean()) / np.sqrt(
        (a.loc[m].var(ddof=1) + b.loc[m].var(ddof=1)) / 2)
    res.append([m, a.loc[m].mean(), b.loc[m].mean(), d, p, cov[m][2]])
res = pd.DataFrame(res, columns=['module', 'mean_ATL', 'mean_control',
                                 'cohens_d', 'p', 'coverage'])
res['q'] = stats.false_discovery_control(res['p'])
res['verdict'] = np.where(res['coverage'] < MIN_COVERAGE, 'not testable',
                  np.where(res['q'] < FDR,
                           np.where(res['cohens_d'] > 0, 'UP in ATL', 'DOWN in ATL'),
                           'no shift'))
print('== module-level result, ATL vs healthy CD4+ ==')
print(res.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
print()

# ---------- acetylation axis ----------
ax = []
for g in axis_genes:
    if g not in G.index:
        continue
    va, vb = G.loc[g, (grp == 'ATL').values], G.loc[g, (grp == 'control').values]
    t, p = stats.mannwhitneyu(va, vb, alternative='two-sided')
    ax.append([g, va.mean(), vb.mean(), va.mean() - vb.mean(), p])
ax = pd.DataFrame(ax, columns=['gene', 'mean_ATL', 'mean_control', 'log2FC', 'p'])
ax['q'] = stats.false_discovery_control(ax['p'])
print('== acetylation axis ==')
print(ax.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
print()

# ---------- axis ratio vs module scores ----------
if all(g in G.index for g in ['EP300', 'SIRT1']):
    ratio = G.loc['EP300'] - G.loc['SIRT1']
    print('== Spearman: EP300-SIRT1 ratio vs module score ==')
    for m in S.index:
        r, p = stats.spearmanr(S.loc[m], ratio)
        print(f'  {m:11s} rho={r:+.3f}  p={p:.3g}')

S.T.join(samples.set_index('gsm')).to_csv('/home/claude/htlv/results_module_scores.csv')
res.to_csv('/home/claude/htlv/results_modules.csv', index=False)
ax.to_csv('/home/claude/htlv/results_axis.csv', index=False)
samples.to_csv('/home/claude/htlv/results_sample_assignment.csv', index=False)
