"""Spatial dependence in the chromatin arm.

Three questions, in order:
 1. Does the promoter-level change correlate between nearby genes, and over what distance?
 2. Do the module genes actually cluster, or are they scattered as Section 2.8 asserts?
 3. If sets are drawn as contiguous genomic blocks instead of scattered genes, how much
    wider does the null get, and therefore how much larger is the true detection bound?

Also revalidates the normal approximation at the smallest module size in every contrast,
not only at n = 13 in Jurkat.
"""
import os
os.makedirs('results', exist_ok=True)
import numpy as np, pandas as pd, sys

SEED = 20260817
rng = np.random.default_rng(SEED)

reg = pd.read_csv('regions_hg38.tsv', sep='\t')
mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module') if m != 'axis'}

ARMS = {
 'H3K27ac Jurkat': 'results/results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv',
 'H3K27ac TL-Om1': 'results/results_h3k27ac_lfc_CT_TL-Om1_Tax_H3K27ac.csv',
 'ATAC empty':     'results/results_atac_lfc_Empty_Tax.csv',
 'ATAC IRF4':      'results/results_atac_lfc_IRF4_Tax.csv',
}

# genomic order: chrom, then TSS
reg = reg.sort_values(['chrom', 'tss']).reset_index(drop=True)
pos = reg.set_index('symbol')[['chrom', 'tss']]

# ---------------------------------------------------------------- 1. autocorrelation
BINS = [(0, 1e5), (1e5, 5e5), (5e5, 1e6), (1e6, 5e6), (5e6, 2e7), (2e7, 1e9)]
print('=' * 74)
print('1. SPATIAL AUTOCORRELATION OF PROMOTER CHANGE (same chromosome)')
print('=' * 74)
auto_rows = []
for arm, path in ARMS.items():
    lfc = pd.read_csv(path, index_col=0).iloc[:, 0]
    d = reg[reg.symbol.isin(lfc.index)].copy()
    d['lfc'] = lfc.reindex(d.symbol).values
    print(f'  {arm}: {len(d)} genes with coordinates and signal', flush=True)
    out = []
    for lo, hi in BINS:
        xs, ys = [], []
        for ch, g in d.groupby('chrom'):
            t = g.tss.values; v = g.lfc.values
            n = len(t)
            # sample pairs rather than enumerate all, for speed
            for _ in range(min(4000, n * 6)):
                i = rng.integers(0, n)
                j = rng.integers(0, n)
                if i == j:
                    continue
                dist = abs(t[i] - t[j])
                if lo <= dist < hi:
                    xs.append(v[i]); ys.append(v[j])
        r = np.corrcoef(xs, ys)[0, 1] if len(xs) > 200 else np.nan
        out.append((lo, hi, len(xs), r))
        auto_rows.append([arm, f'{lo/1e3:.0f}-{hi/1e3:.0f} kb', len(xs), round(r, 4) if r == r else None])
    for lo, hi, n, r in out:
        lab = f'{lo/1e6:>6.2f} - {hi/1e6:<6.2f} Mb'
        print(f'     {lab}  pairs {n:>6}  r = {r:+.4f}' if r == r else f'     {lab}  too few pairs')
    print(flush=True)

pd.DataFrame(auto_rows, columns=['contrast', 'distance_bin', 'n_pairs', 'r']).to_csv(
    'results/results_spatial_autocorrelation.csv', index=False)

# ---------------------------------------------------------------- 2. module clustering
print('=' * 74)
print('2. DO MODULE GENES CLUSTER?  (nearest same-module neighbor on the same chromosome)')
print('=' * 74)
clus_rows = []
for m, genes in mods.items():
    p = pos.reindex([g for g in genes if g in pos.index]).dropna()
    n = len(p)
    close_1mb = 0
    nearest = []
    for ch, g in p.groupby('chrom'):
        t = np.sort(g.tss.values)
        if len(t) < 2:
            continue
        gaps = np.diff(t)
        for k in range(len(t)):
            cand = []
            if k > 0: cand.append(gaps[k-1])
            if k < len(t)-1: cand.append(gaps[k])
            nearest.append(min(cand))
        close_1mb += int((gaps < 1e6).sum())
    same_chrom_pairs = sum(len(g)*(len(g)-1)//2 for _, g in p.groupby('chrom'))
    med = np.median(nearest)/1e6 if nearest else np.nan
    print(f'  {m:11s} n={n:2d}  genes sharing a chromosome with another module gene: '
          f'{sum(len(g) for _, g in p.groupby("chrom") if len(g) > 1):2d}  '
          f'pairs <1 Mb apart: {close_1mb}  median nearest neighbor: '
          + (f'{med:.1f} Mb' if med == med else 'n/a'), flush=True)
    clus_rows.append([m, n, same_chrom_pairs, close_1mb, round(med, 2) if med == med else None])
pd.DataFrame(clus_rows, columns=['module', 'n_genes', 'same_chrom_pairs',
                                 'pairs_within_1Mb', 'median_nn_Mb']).to_csv(
    'results/results_module_clustering.csv', index=False)
print()

# ---------------------------------------------------------------- 3. block null
print('=' * 74)
print('3. BLOCK-PERMUTATION NULL: set of n genes drawn as k blocks of contiguous genes')
print('   k = n is the scattered null used in the paper; k = 1 is fully contiguous.')
print('=' * 74)
NDRAW = 20000
block_rows = []
for arm, path in ARMS.items():
    lfc = pd.read_csv(path, index_col=0).iloc[:, 0]
    d = reg[reg.symbol.isin(lfc.index)].copy()
    d['lfc'] = lfc.reindex(d.symbol).values
    v = d.lfc.values                      # already in genomic order
    N = len(v)
    s = lfc.std()
    print(f'\n  {arm}   genome-wide SD = {s:.3f}   {N} ordered genes', flush=True)
    for n in (12, 13, 15):
        line = [f'    n={n:2d}']
        base = None
        for k in (n, 6, 4, 2, 1):
            if k > n:
                continue
            sizes = [n // k + (1 if i < n % k else 0) for i in range(k)]
            draws = np.empty(NDRAW)
            for it in range(NDRAW):
                tot = 0.0
                for bs in sizes:
                    st = rng.integers(0, N - bs)
                    tot += v[st:st + bs].sum()
                draws[it] = tot / n
            sd = draws.std()
            if base is None:
                base = sd
            mde = 2.80 * sd
            line.append(f'k={k:<2d} sd={sd:.4f} MDE={mde:.3f} ({2**mde:.2f}x) infl={sd/base:.2f}')
            block_rows.append([arm, n, k, round(sd, 5), round(mde, 4),
                               round(2 ** mde, 3), round(sd / base, 3)])
            print(f'      n={n:2d} k={k:<2d}  null sd={sd:.4f}  MDE={mde:.3f} log2 '
                  f'({2**mde:.2f}-fold)  inflation vs scattered = {sd/base:.2f}x', flush=True)
pd.DataFrame(block_rows, columns=['contrast', 'n_genes', 'n_blocks', 'null_sd',
                                  'MDE_log2', 'MDE_fold', 'inflation_vs_scattered']).to_csv(
    'results/results_block_permutation.csv', index=False)

# ------------------------------------------- 4. normal approximation, every contrast
print()
print('=' * 74)
print('4. NORMAL APPROXIMATION CHECKED AT THE SMALLEST MODULE SIZE IN EVERY CONTRAST')
print('=' * 74)
norm_rows = []
for arm, path in ARMS.items():
    x = pd.read_csv(path, index_col=0).iloc[:, 0].values
    s = x.std()
    for n in (12, 13, 15):
        perm = np.array([rng.choice(x, n, replace=False).mean() for _ in range(20000)])
        approx = s / np.sqrt(n)
        norm_rows.append([arm, n, round(perm.std(), 5), round(approx, 5),
                          round(perm.std() / approx, 4), round(perm.mean(), 5)])
        print(f'  {arm:16s} n={n:2d}  permutation sd={perm.std():.4f}  '
              f's/sqrt(n)={approx:.4f}  ratio={perm.std()/approx:.4f}  '
              f'null mean={perm.mean():+.5f}', flush=True)
pd.DataFrame(norm_rows, columns=['contrast', 'n_genes', 'permutation_sd', 'normal_approx_sd',
                                 'ratio', 'null_mean']).to_csv(
    'results/results_normal_approximation_check.csv', index=False)
print('\nDONE. Four result tables written.')
