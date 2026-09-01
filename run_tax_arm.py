import numpy as np, pandas as pd, pyreadr
from scipy import stats

ann = pyreadr.read_r('/tmp/grch38.rda')['grch38'][['ensgene','symbol']].dropna().drop_duplicates('ensgene')
sym = ann.set_index('ensgene')['symbol']

C = pd.read_csv('/mnt/user-data/uploads/GSE316417_RNA-seq_raw_counts_all_samples_txt.gz', sep='\t', index_col=0)
C.columns = [c.replace('_count','') for c in C.columns]
C = C[C.sum(axis=1) > 0]
cpm = C / C.sum() * 1e6
keep = (cpm > 1).sum(axis=1) >= 3
cpm = cpm[keep]
g = sym.reindex(cpm.index)
cpm = cpm[g.notna().values]; g = g.dropna()
L = np.log2(cpm + 1)
L.index = g.values
L = L.groupby(level=0).max()
print(f'expressed genes with a symbol: {L.shape[0]}  samples: {L.shape[1]}')

mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module')}
axis = mods.pop('axis')
print('coverage:', {m: round(len([x for x in gs if x in L.index])/len(gs), 2)
                    for m, gs in list(mods.items())+[('axis', axis)]})

CONTRASTS = [
    ('JPX-9, Tax induction (no IRF4)',      'JPX9_Empty',      'JPX9_Empty_Tax'),
    ('JPX-9, Tax induction (with IRF4)',    'JPX9_IRF4',       'JPX9_IRF4_Tax'),
    ('Jurkat, Tax on IRF4-WT background',   'Jurkat_IRF4_WT',  'Jurkat_IRF4_WT_Tax'),
]

rng = np.random.default_rng(20260817)
def competitive(lfc, genes, n=20000):
    """Is the module's mean log-ratio larger than random gene sets of the same size?"""
    idx = [g for g in genes if g in lfc.index]
    if len(idx) < 5: return np.nan, np.nan, 0
    obs = lfc.loc[idx].mean()
    pool = lfc.values
    null = np.array([rng.choice(pool, len(idx), replace=False).mean() for _ in range(n)])
    p = 2 * min((null >= obs).mean(), (null <= obs).mean())
    z = (obs - null.mean()) / null.std()
    return obs, max(p, 1/n), z, len(idx)

for label, base, tax in CONTRASTS:
    if base not in L.columns or tax not in L.columns: continue
    lfc = L[tax] - L[base]
    print(f'\n== {label} ==   (single library per condition, descriptive)')
    print(f'   global median shift {lfc.median():+.3f}')
    for m, genes in mods.items():
        obs, p, z, n = competitive(lfc, genes)
        print(f'   {m:11s} mean log2FC={obs:+.3f}  z={z:+.2f}  p_competitive={p:.4f}  (n={n})')
    for gname in axis + ['SQSTM1','MAP1LC3B','BECN1','LAMP1','CTSD','TFEB','NFKB1','RELA','RELB','NFKB2']:
        if gname in lfc.index:
            print(f'   {gname:9s} log2FC={lfc[gname]:+.2f}')

# tumor-line perturbations with replicates: JQ1 in C91/PL (n=3 vs 3)
dm = [c for c in L.columns if 'DMSO' in c]; jq = [c for c in L.columns if 'JQ1' in c]
if len(dm) == 3 and len(jq) == 3:
    print('\n== C91/PL, JQ1 vs DMSO (n=3 each, the only replicated contrast) ==')
    lfc = L[jq].mean(axis=1) - L[dm].mean(axis=1)
    for m, genes in mods.items():
        idx = [g for g in genes if g in L.index]
        t, p = stats.mannwhitneyu(L.loc[idx, jq].mean(axis=1), L.loc[idx, dm].mean(axis=1))
        obs, pc, z, n = competitive(lfc, genes)
        print(f'   {m:11s} mean log2FC={obs:+.3f}  z={z:+.2f}  p_competitive={pc:.4f}')
