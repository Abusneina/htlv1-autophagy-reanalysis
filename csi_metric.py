"""Composition Sensitivity Index (CSI).

Definition. For a gene set S, CSI(S) is the change in standardized effect size
(Cohen d, synthetic specimen versus pure CD4) per 10 percentage points of
non-T-cell content, estimated by least squares across mixing levels of 2 to 30
percent. Mixtures are built per donor from purified populations, so between-sample
variation is real donor variation. CSI is positive when contamination inflates a
set's score.

Null model. For each set size, 300 random gene sets of the same size are drawn
from the expressed transcriptome and their CSI computed. Each set's z score and
percentile are reported against that size-matched null, so the index is not a
function of set size.

Interpretation. |z| > 3 high risk, |z| > 2 caution, otherwise robust.
"""
import glob
import numpy as np, pandas as pd, pyreadr
rng = np.random.default_rng(20260817)

d = pd.read_csv('/mnt/user-data/uploads/GSE107011_Processed_data_TPM_txt.gz', sep='\t', index_col=0)
d.index = [i.split('.')[0] for i in d.index]
ann = pyreadr.read_r('/tmp/grch38.rda')['grch38'][['ensgene', 'symbol']].dropna().drop_duplicates('ensgene')
g = ann.set_index('ensgene')['symbol'].reindex(d.index)
d = d[g.notna().values]; d.index = g.dropna().values
d = d.groupby(level=0).max()

T = ['CD4_naive', 'Th1', 'Th2', 'Th17', 'Treg']
O = {'C_mono': 0.40, 'I_mono': 0.08, 'NC_mono': 0.07, 'B_naive': 0.25, 'NK': 0.12, 'Neutrophils': 0.08}
donors = sorted(dn for dn in {c.split('_', 1)[0] for c in d.columns}
                if all(f'{dn}_{ct}' in d.columns for ct in T + list(O)))
FR = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]

cols = {}
for dn in donors:
    cd4 = d[[f'{dn}_{ct}' for ct in T]].mean(axis=1)
    oth = sum(w * d[f'{dn}_{ct}'] for ct, w in O.items())
    for f in FR:
        cols[f'{dn}_f{int(f*100)}'] = (1 - f) * cd4 + f * oth
X = np.log2(pd.DataFrame(cols) + 1)
X = X[X.var(axis=1) > 0]
print(f'{len(donors)} donors, {X.shape[1]} synthetic specimens, {X.shape[0]} genes')

genes_arr = X.index.to_numpy()
R = X.rank(axis=0, method='average').values
N = X.shape[0]
ORDERS = [np.argsort(-R[:, j]) for j in range(X.shape[1])]
RS = [R[o, j] for j, o in enumerate(ORDERS)]

POS = {g: i for i, g in enumerate(genes_arr)}
RANKPOS = [np.empty(N, dtype=np.int64) for _ in ORDERS]
for j_, o in enumerate(ORDERS):
    RANKPOS[j_][o] = np.arange(N)
WPOW = [RS[j_] ** 0.25 for j_ in range(len(ORDERS))]

def score(gs, alpha=0.25):
    """ssGSEA enrichment score, summed over all N ranked positions.

    Between consecutive hits the in-set cumulative fraction is constant while the
    out-of-set fraction keeps rising, so the full sum is obtained in closed form
    from the hit positions alone rather than by walking the whole list.
    """
    ii = np.fromiter((POS[g] for g in gs if g in POS), dtype=np.int64)
    k = ii.size
    if k < 5:
        return None
    out = []
    for j_ in range(len(ORDERS)):
        p = np.sort(RANKPOS[j_][ii])              # 0-based positions of hits
        w = WPOW[j_][p]
        tail = N - p                              # positions from each hit to the end
        s_in = float(np.sum(w * tail) / w.sum())  # sum over i of cumulative in-set fraction
        s_hits = float(np.sum(tail))              # sum over i of hit count up to i
        s_out = (N * (N + 1) / 2 - s_hits) / (N - k)
        out.append(s_in - s_out)
    return pd.Series(out, index=X.columns)

def score_bruteforce(gs, alpha=0.25):
    idx = np.isin(genes_arr, gs); k = idx.sum()
    out = []
    for j_, o in enumerate(ORDERS):
        ins = idx[o]; w = (RS[j_] ** alpha) * ins
        out.append(np.sum(np.cumsum(w) / w.sum() - np.cumsum(~ins) / (N - k)))
    return pd.Series(out, index=X.columns)

_test = list(rng.choice(genes_arr, 120, replace=False))
_a, _b = score(_test), score_bruteforce(_test)
print(f'fast vs brute-force scorer: max abs difference {np.abs(_a - _b).max():.2e}')
assert np.abs(_a - _b).max() < 1e-6, 'fast scorer does not reproduce the reference'

def csi(sc):
    base = np.array([sc[f'{dn}_f0'] for dn in donors])
    xs, ys = [], []
    for f in FR[1:]:
        cur = np.array([sc[f'{dn}_f{int(f*100)}'] for dn in donors])
        sd = np.sqrt((cur.var(ddof=1) + base.var(ddof=1)) / 2)
        ys.append((cur.mean() - base.mean()) / sd if sd > 0 else 0.0)
        xs.append(f)
    return np.polyfit(xs, ys, 1)[0] / 10

sets = {}
mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
for m, v in mt.groupby('module'):
    if m != 'axis':
        gs = [x for x in set(v['gene']) if x in X.index]
        if len(gs) >= 10:
            sets['MODULE_' + m.upper()] = gs
COLLECTION={'h.all.v7.0.symbols.gmt':'Hallmark','kegg.gmt':'KEGG','kegg_medicus.gmt':'KEGG_MEDICUS','reactome.gmt':'Reactome'}
origin={}
for path in sorted(glob.glob('code/*.gmt')):
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        gs = [x for x in p[2:] if x in X.index]
        if len(gs) >= 15:
            sets[p[0]] = gs; origin[p[0]] = COLLECTION.get(path.split('/')[-1], path)
print(f'{len(sets)} gene sets scored')

rows = []
for name, gs in sets.items():
    sc = score(gs)
    if sc is not None:
        rows.append([name, len(gs), csi(sc)])
res = pd.DataFrame(rows, columns=['set', 'n_genes', 'CSI'])
res['collection'] = [origin.get(s_, 'module') for s_ in res.set]

print('building size-matched null...')
sizes = sorted({int(b) for b in np.clip(np.round(res.n_genes/25)*25, 15, 500)})
null = {}
for sz in sizes:
    vals = []
    for _ in range(120):
        sc = score(list(rng.choice(X.index, sz, replace=False)))
        vals.append(csi(sc))
    null[sz] = np.array(vals)
print(f'  {len(sizes)} size bins, null CSI sd {np.mean([v.std() for v in null.values()]):.3f}')

def zp(row):
    sz = min(sizes, key=lambda s: abs(s - row.n_genes))
    n = null[sz]
    return pd.Series({'z': (row.CSI - n.mean()) / n.std(), 'pct': float((n < row.CSI).mean())})

res = pd.concat([res, res.apply(zp, axis=1)], axis=1)
res['risk'] = np.where(res.z.abs() > 3, 'high', np.where(res.z.abs() > 2, 'caution', 'robust'))
res = res.sort_values('CSI', ascending=False)
res.to_csv('results_CSI_index.csv', index=False)

print('\nby collection:'); print(res.groupby('collection').risk.value_counts().unstack(fill_value=0).to_string())
print('\n== Composition Sensitivity Index, most sensitive ==')
print(res.head(10).to_string(index=False, float_format=lambda x: f'{x:.3f}'))
print('\n== least sensitive ==')
print(res.tail(6).to_string(index=False, float_format=lambda x: f'{x:.3f}'))
print('\nrisk classes:', res.risk.value_counts().to_dict())
print('\nautophagy modules:')
print(res[res.set.str.startswith('MODULE_')].to_string(index=False, float_format=lambda x: f'{x:.3f}'))
