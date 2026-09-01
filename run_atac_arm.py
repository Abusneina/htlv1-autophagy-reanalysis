"""ATAC arm: chromatin accessibility at autophagy-module promoters, Tax vs no Tax,
JPX-9 system, on the empty and IRF4 backgrounds. Single library per condition."""
import glob, numpy as np, pandas as pd, pyBigWig
PROMOTER = 1000                     # ATAC is sharper than H3K27ac; TSS +/- 1 kb
rng = np.random.default_rng(20260817)

reg = pd.read_csv('regions_hg38.tsv', sep='\t')
mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module')}
axis = mods.pop('axis')

def quantify(path):
    bw = pyBigWig.open(path); ch = bw.chroms(); out = {}
    for r in reg.itertuples():
        if r.chrom not in ch: continue
        s, e = max(0, r.tss - PROMOTER), min(ch[r.chrom], r.tss + PROMOTER)
        try: v = bw.stats(r.chrom, s, e, type='mean')[0]
        except Exception: v = None
        out[r.symbol] = 0.0 if v is None else float(v)
    bw.close(); return pd.Series(out)

S = {}
for f in sorted(glob.glob('atac/*.bw')):
    n = f.split('/')[-1].split('_ATAC_JPX-9_')[1].replace('.bw','')
    S[n] = quantify(f); print(f'  {n}: median {S[n].median():.3f}')
S = pd.DataFrame(S); S.to_csv('results_atac_promoters.csv')

def competitive(lfc, genes, n=20000):
    idx = [g for g in genes if g in lfc.index]
    if len(idx) < 5: return np.nan, np.nan, np.nan, 0
    obs = lfc.loc[idx].mean(); pool = lfc.values
    null = np.array([rng.choice(pool, len(idx), replace=False).mean() for _ in range(n)])
    return obs, max(2*min((null>=obs).mean(), (null<=obs).mean()), 1/n), (obs-null.mean())/null.std(), len(idx)

for base, tax, label in [('Empty','Empty_Tax','Tax on empty background'),
                         ('IRF4','IRF4_Tax','Tax on IRF4 background')]:
    x, y = S[base].copy(), S[tax].copy()
    keep = (x + y) > 0; x, y = x[keep], y[keep]
    ry = y.rank(method='first')
    yq = pd.Series(np.interp(ry, np.arange(1, len(y)+1), np.sort(x.values)), index=y.index)
    lfc = np.log2(yq + 0.1) - np.log2(x + 0.1)
    print(f'\n== {label}: {tax} vs {base} ==  ({len(lfc)} promoters, quantile-matched)')
    print(f'   dynamic range sd={lfc.std():.3f}')
    for m, genes in mods.items():
        o, p, z, n = competitive(lfc, genes)
        print(f'   {m:11s} mean log2FC={o:+.3f}  z={z:+.2f}  p_competitive={p:.4f}  (n={n})')
    pct = lfc.rank(pct=True)
    for g in axis + ['RELB','NFKB2','BATF3','DUSP10','TFEB','TFE3','MAP1LC3B','SQSTM1','CTSD','LAMP1','BECN1']:
        if g in lfc.index:
            print(f'   {g:9s} log2FC={lfc[g]:+.2f}  pct={pct[g]:.3f}')
    lfc.to_csv(f'results_atac_lfc_{tax}.csv')
