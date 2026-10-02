"""Gene-level composition sensitivity.

Rebuilds the synthetic specimens exactly as csi_metric.py does, confirms that the
module-level index reproduces the published values, then computes the same slope for
every expressed gene: change in standardized effect (synthetic specimen versus pure
CD4) per 10 percentage points of non-T content. A per-gene table lets any gene set,
including licensed collections not in the deposit, be decomposed by lookup.
"""
import os
os.makedirs('results', exist_ok=True)
import sys
# Input: the GEO file GSE107011_Processed_data_TPM.txt.gz placed in this folder (not committed),
# or pass its path as the first argument.
GSE107011 = sys.argv[1] if len(sys.argv) > 1 else 'GSE107011_Processed_data_TPM.txt.gz'
import numpy as np, pandas as pd, pyreadr, time

t0 = time.time()
SEED = 20260817
rng = np.random.default_rng(SEED)

print('reading GSE107011 TPM matrix ...', flush=True)
d = pd.read_csv(GSE107011, sep='\t', index_col=0)
d.index = [i.split('.')[0] for i in d.index]
ann = pyreadr.read_r('grch38.rda')['grch38'][['ensgene', 'symbol']].dropna().drop_duplicates('ensgene')
g = ann.set_index('ensgene')['symbol'].reindex(d.index)
d = d[g.notna().values]; d.index = g.dropna().values
d = d.groupby(level=0).max()
print(f'  {d.shape[0]} genes x {d.shape[1]} samples  ({time.time()-t0:.0f}s)', flush=True)

T = ['CD4_naive', 'Th1', 'Th2', 'Th17', 'Treg']
O = {'C_mono': 0.40, 'I_mono': 0.08, 'NC_mono': 0.07, 'B_naive': 0.25, 'NK': 0.12, 'Neutrophils': 0.08}
donors = sorted(dn for dn in {c.split('_', 1)[0] for c in d.columns}
                if all(f'{dn}_{ct}' in d.columns for ct in T + list(O)))
FR = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
print(f'  donors with every population: {donors}', flush=True)

cols = {}
for dn in donors:
    cd4 = d[[f'{dn}_{ct}' for ct in T]].mean(axis=1)
    oth = sum(w * d[f'{dn}_{ct}'] for ct, w in O.items())
    for f in FR:
        cols[f'{dn}_f{int(f*100)}'] = (1 - f) * cd4 + f * oth
X = np.log2(pd.DataFrame(cols) + 1)
X = X[X.var(axis=1) > 0]
print(f'  {X.shape[1]} synthetic specimens, {X.shape[0]} genes with variance', flush=True)

# ---------- module-level reproduction, same scorer as csi_metric.py ----------
genes_arr = X.index.to_numpy(); N = X.shape[0]
R = X.rank(axis=0, method='average').values
ORDERS = [np.argsort(-R[:, j]) for j in range(X.shape[1])]
RS = [R[o, j] for j, o in enumerate(ORDERS)]
POS = {gg: i for i, gg in enumerate(genes_arr)}
RANKPOS = [np.empty(N, dtype=np.int64) for _ in ORDERS]
for j_, o in enumerate(ORDERS):
    RANKPOS[j_][o] = np.arange(N)
WPOW = [RS[j_] ** 0.25 for j_ in range(len(ORDERS))]

def score(gs):
    ii = np.fromiter((POS[x] for x in gs if x in POS), dtype=np.int64)
    k = ii.size
    out = []
    for j_ in range(len(ORDERS)):
        p = np.sort(RANKPOS[j_][ii]); w = WPOW[j_][p]; tail = N - p
        s_in = float(np.sum(w * tail) / w.sum()); s_hits = float(np.sum(tail))
        s_out = (N * (N + 1) / 2 - s_hits) / (N - k)
        out.append(s_in - s_out)
    return pd.Series(out, index=X.columns)

def slope(vals):
    base = np.array([vals[f'{dn}_f0'] for dn in donors])
    xs, ys = [], []
    for f in FR[1:]:
        cur = np.array([vals[f'{dn}_f{int(f*100)}'] for dn in donors])
        sd = np.sqrt((cur.var(ddof=1) + base.var(ddof=1)) / 2)
        ys.append((cur.mean() - base.mean()) / sd if sd > 0 else 0.0)
        xs.append(f)
    return np.polyfit(xs, ys, 1)[0] / 10

mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
pub = pd.read_csv('results/results_CSI_index.csv').set_index('set')['CSI']
print('\nmodule-level reproduction against the published index:', flush=True)
for m, v in mt.groupby('module'):
    if m == 'axis':
        continue
    gs = [x for x in set(v['gene']) if x in X.index]
    val = slope(score(gs))
    ref = pub.get('MODULE_' + m.upper(), np.nan)
    print(f'  MODULE_{m.upper():11s} recomputed {val:+.4f}   published {ref:+.4f}', flush=True)

# ---------- per-gene slope ----------
print('\nper-gene composition sensitivity ...', flush=True)
res = {}
for i, gene in enumerate(X.index):
    res[gene] = slope(X.loc[gene])
    if (i + 1) % 5000 == 0:
        print(f'  {i+1}/{X.shape[0]} genes  ({time.time()-t0:.0f}s)', flush=True)
G = pd.Series(res, name='gene_CSI').to_frame()
mean_cd4 = X[[c for c in X.columns if c.endswith('_f0')]].mean(axis=1)
oth_mean = pd.DataFrame({dn: np.log2(sum(w * d[f'{dn}_{ct}'] for ct, w in O.items()) + 1) for dn in donors}).mean(axis=1)
G['log2_TPM_CD4'] = mean_cd4.reindex(G.index)
G['log2_TPM_nonT_pool'] = oth_mean.reindex(G.index)
G['nonT_minus_CD4'] = G.log2_TPM_nonT_pool - G.log2_TPM_CD4
G['percentile'] = G.gene_CSI.rank(pct=True)
G = G.sort_values('gene_CSI', ascending=False)
G.to_csv('results/results_gene_level_CSI.csv')
print(f'  wrote results_gene_level_CSI.csv ({len(G)} genes)', flush=True)

# ---------- CLEAR module decomposition ----------
clear = [x for x in set(mt[mt.module == 'clear'].gene)]
C = G.reindex(clear).dropna().sort_values('gene_CSI', ascending=False)
coh = pd.read_csv('results/results_composition_per_gene.csv')
coh = coh[coh.module == 'clear'].set_index('gene')   # output of composition_analysis.py
C['rho_monocyte_GSE33615'] = coh['rho_monocyte'].reindex(C.index)
C['cohens_d_ATL_GSE33615'] = coh['cohens_d_ATL'].reindex(C.index)
C.to_csv('results/results_clear_gene_level_CSI.csv')
print('\nCLEAR module, gene by gene (mixtures vs the unsorted cohort):', flush=True)
print(C[['gene_CSI', 'percentile', 'nonT_minus_CD4', 'rho_monocyte_GSE33615', 'cohens_d_ATL_GSE33615']].round(3).to_string(), flush=True)
from scipy import stats
both = C.dropna(subset=['rho_monocyte_GSE33615'])
r, p = stats.spearmanr(both.gene_CSI, both.rho_monocyte_GSE33615)
print(f'\nagreement between mixture sensitivity and cohort monocyte tracking across CLEAR genes: '
      f'Spearman rho = {r:+.3f}, P = {p:.3g}, n = {len(both)}', flush=True)
print(f'done in {time.time()-t0:.0f}s')
