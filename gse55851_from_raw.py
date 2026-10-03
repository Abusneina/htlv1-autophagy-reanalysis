"""Rebuild the sorted-cohort (GSE55851) expression matrix from the per-sample Agilent files in
GSE55851_RAW.tar and test whether it reproduces the deposited module scores.

The deposited pipeline (run_gse55851.py) starts from GEO's processed series matrix. This script
starts one step earlier, from the Agilent Feature Extraction files, and applies the identical
steps afterwards: log2 with a floor of 1, quantile normalization across samples, probe to gene
collapse by highest mean, ssGSEA with exponent 0.25, and min-max scaling across all scores.
Each candidate signal column is tested against results/results_gse55851_scores.csv.
"""
import gzip, io, os, re, sys, tarfile, numpy as np, pandas as pd
from scipy import stats
TAR = sys.argv[1] if len(sys.argv) > 1 else 'GSE55851_RAW.tar'
SIGNALS = ['gProcessedSignal', 'gMedianSignal', 'gBGSubSignal']

dep = pd.read_csv('results/results_gse55851_scores.csv', index_col=0)       # deposited scores
order = list(dep.gsm)                                                       # sample order of the deposited run
mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(g['gene'])) for m, g in mt.groupby('module')}; axis = mods.pop('axis')

print('reading', TAR, flush=True)
raw, ann_txt = {}, None
with tarfile.open(TAR) as tf:
    for mem in tf.getmembers():
        data = gzip.decompress(tf.extractfile(mem).read()).decode('utf-8', 'replace')
        if mem.name.startswith('GPL'): ann_txt = data
        else: raw[mem.name.split('_')[0]] = data
print(f'  {len(raw)} sample files', flush=True)

def parse(txt):
    lines = txt.split('\n')
    h = [i for i, l in enumerate(lines) if l.startswith('FEATURES')][0]
    df = pd.read_csv(io.StringIO('\n'.join(lines[h:])), sep='\t', low_memory=False)
    df = df[df['FEATURES'] == 'DATA'] if 'FEATURES' in df.columns else df
    return df.set_index('FeatureNum')

feat = {g: parse(t) for g, t in raw.items()}
first = feat[order[0]]
print(f'  {len(first)} features per array; first sample: {order[0]}', flush=True)

# annotation
al = ann_txt.split('\n'); h = [i for i, l in enumerate(al) if l.startswith('ID\t')][0]
ann = pd.read_csv(io.StringIO('\n'.join(al[h:])), sep='\t', low_memory=False)
ann['ID'] = pd.to_numeric(ann['ID'], errors='coerce'); ann = ann.dropna(subset=['ID']).set_index('ID'); ann.index = ann.index.astype(int)
print(f'  annotation rows {len(ann)}, with gene symbol {(ann.GENE_SYMBOL.fillna("")!="").sum()}', flush=True)

def ssgsea(expr, gene_sets, alpha=0.25):
    genes = expr.index.to_numpy(); R = expr.rank(axis=0, method='average').values; N = expr.shape[0]; out = {}
    for name, gs in gene_sets.items():
        idx = np.isin(genes, gs)
        if idx.sum() < 5: continue
        sc = []
        for j in range(expr.shape[1]):
            o = np.argsort(-R[:, j]); ins = idx[o]; w = (R[o, j] ** alpha) * ins
            sc.append(np.sum(np.cumsum(w) / w.sum() - np.cumsum(~ins) / (N - ins.sum())))
        out[name] = sc
    S = pd.DataFrame(out, index=expr.columns).T
    return (S - S.values.min()) / (S.values.max() - S.values.min())

def gene_matrix(sig):
    X = pd.DataFrame({g: feat[g][sig] for g in order}).apply(pd.to_numeric, errors='coerce')
    X = np.log2(X.clip(lower=1))
    ranks = X.rank(method='first'); ms = np.sort(X.values, axis=0).mean(axis=1)
    X = pd.DataFrame(np.interp(ranks.values, np.arange(1, X.shape[0] + 1), ms), index=X.index, columns=X.columns)
    sym = ann.reindex(X.index)['GENE_SYMBOL'].fillna('')
    ctrl = ann.reindex(X.index).get('CONTROL_TYPE', pd.Series('', index=X.index)).fillna('')
    keep = (sym != '') & (ctrl.astype(str).str.lower().isin(['', 'false']))
    X, sym = X[keep], sym[keep]
    o = X.mean(axis=1).sort_values(ascending=False).index; X, sym = X.loc[o], sym.loc[o]
    d = ~sym.duplicated(); G = X[d]; G.index = sym[d]
    return G

if __name__ == '__main__':
    best = None
    for sig in SIGNALS:
        print(f'\ncandidate signal: {sig}', flush=True)
        G = gene_matrix(sig); print(f'  gene matrix {G.shape}', flush=True)
        S = ssgsea(G, mods).T
        S.index = order
        cmp_ = dep.set_index('gsm')[S.columns]
        diff = (S - cmp_).abs().max().max()
        r = np.corrcoef(S.values.ravel(), cmp_.values.ravel())[0, 1]
        print(f'  against the deposited scores: max |difference| = {diff:.4f}, correlation over all 84 values = {r:.4f}', flush=True)
        if best is None or diff < best[1]: best = (sig, diff, r)
    print(f'\nbest candidate: {best[0]} (max |difference| {best[1]:.4f}, correlation {best[2]:.4f})')
