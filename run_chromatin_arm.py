"""Chromatin arm: H3K27ac signal at autophagy-module promoters, Tax vs empty.
Single library per condition, so this is descriptive; significance is judged by a
competitive permutation against size-matched random gene sets."""
import sys, glob, numpy as np, pandas as pd, pyBigWig

PROMOTER = 2000          # +/- bp around TSS, fixed before analysis
rng = np.random.default_rng(20260817)

reg = pd.read_csv('regions_hg38.tsv', sep='\t')
mods_tbl = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mods_tbl.groupby('module')}
axis = mods.pop('axis')

def quantify(path):
    bw = pyBigWig.open(path)
    chroms = bw.chroms()
    vals = {}
    for r in reg.itertuples():
        if r.chrom not in chroms: continue
        s = max(0, r.tss - PROMOTER); e = min(chroms[r.chrom], r.tss + PROMOTER)
        try:
            v = bw.stats(r.chrom, s, e, type='mean')[0]
        except Exception:
            v = None
        vals[r.symbol] = 0.0 if v is None else float(v)
    bw.close()
    return pd.Series(vals)

files = sorted(glob.glob('bw2/*.bw'))
print('tracks found:'); [print('  ', f.split('/')[-1]) for f in files]
sig = {}
for f in files:
    name = f.split('/')[-1].replace('GSM9453007_','').replace('GSM9453008_','').replace('GSM9453009_','').replace('GSM9453011_','').replace('.bw','')
    sig[name] = quantify(f)
    print(f'  quantified {name}: {len(sig[name])} promoters, median {sig[name].median():.3f}')
S = pd.DataFrame(sig)
S.to_csv('results_h3k27ac_promoters.csv')

def competitive(lfc, genes, n=20000):
    idx = [g for g in genes if g in lfc.index]
    if len(idx) < 5: return np.nan, np.nan, np.nan, 0
    obs = lfc.loc[idx].mean(); pool = lfc.values
    null = np.array([rng.choice(pool, len(idx), replace=False).mean() for _ in range(n)])
    return obs, max(2*min((null>=obs).mean(), (null<=obs).mean()), 1/n), (obs-null.mean())/null.std(), len(idx)

pairs = []
cols = list(S.columns)
for a in cols:
    for b in cols:
        if 'Empty' in a and 'Tax' in b and a.split('_')[1] == b.split('_')[1]:
            pairs.append((a, b))
if not pairs:
    print('\nOnly one condition present. Reporting promoter signal ranks only.')
    for m, genes in mods.items():
        idx = [g for g in genes if g in S.index]
        pct = (S[cols[0]].rank(pct=True)[idx]).mean()
        print(f'  {m:11s} mean promoter-signal percentile = {pct:.3f}  (n={len(idx)})')
    for g in axis:
        if g in S.index:
            print(f'  {g:8s} percentile = {S[cols[0]].rank(pct=True)[g]:.3f}')
    sys.exit()

# quantile-normalize the two tracks before differencing (depth differs)
for empty, tax in pairs:
    x, y = S[empty].copy(), S[tax].copy()
    keep = (x + y) > 0
    x, y = x[keep], y[keep]
    ry = y.rank(method='first')
    y_qn = pd.Series(np.interp(ry, np.arange(1, len(y)+1), np.sort(x.values)), index=y.index)
    lfc = np.log2(y_qn + 0.1) - np.log2(x + 0.1)
    print(f'\n== {tax} vs {empty} ==  ({len(lfc)} promoters, quantile-matched)')
    print(f'   global median log2 change {lfc.median():+.3f}')
    for m, genes in mods.items():
        obs, p, z, n = competitive(lfc, genes)
        print(f'   {m:11s} mean log2FC={obs:+.3f}  z={z:+.2f}  p_competitive={p:.4f}  (n={n})')
    for g in axis + ['RELB','NFKB2','BATF3','SQSTM1','MAP1LC3B','TFEB','CTSD','BECN1','LAMP1']:
        if g in lfc.index:
            print(f'   {g:9s} log2FC={lfc[g]:+.2f}')
    lfc.to_csv(f'results_h3k27ac_lfc_{tax}.csv')
