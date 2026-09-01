"""Minimum detectable effect for the competitive gene-set test in each arm.
The null distribution of a size-n random-set mean is approximately normal with
sd = s/sqrt(n) (s = sd of the genome-wide statistic). MDE at 80% power,
two-sided alpha=0.05, is therefore (1.96+0.84)*s/sqrt(n). The normal
approximation is checked against the permutation null before use."""
import numpy as np, pandas as pd
rng = np.random.default_rng(20260817)

mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module') if m != 'axis'}

ARMS = {
 'H3K27ac Jurkat Tax':   'results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv',
 'H3K27ac TL-Om1 Tax':   'results_h3k27ac_lfc_CT_TL-Om1_Tax_H3K27ac.csv',
 'ATAC Tax (empty bg)':  'results_atac_lfc_Empty_Tax.csv',
 'ATAC Tax (IRF4 bg)':   'results_atac_lfc_IRF4_Tax.csv',
}

# validate the normal approximation once
v = pd.read_csv(ARMS['H3K27ac Jurkat Tax'], index_col=0).iloc[:,0].values
n = 13
perm = np.array([rng.choice(v, n, replace=False).mean() for _ in range(20000)])
print(f'permutation null sd = {perm.std():.4f}   normal approx s/sqrt(n) = {v.std()/np.sqrt(n):.4f}')
print(f'permutation null mean = {perm.mean():+.4f}\n')

rows = []
for arm, path in ARMS.items():
    x = pd.read_csv(path, index_col=0).iloc[:,0]
    s = x.std()
    for m, genes in mods.items():
        idx = [g for g in genes if g in x.index]
        n = len(idx)
        mde = 2.80 * s / np.sqrt(n)
        obs = x.loc[idx].mean()
        rows.append([arm, m, n, round(s,3), round(mde,3), round(obs,3),
                     round(abs(obs)/mde, 2)])
R = pd.DataFrame(rows, columns=['arm','module','n_genes','genome_sd',
                                'MDE_log2_80pct','observed','obs_over_MDE'])
print(R.to_string(index=False))
R.to_csv('results_power_analysis.csv', index=False)

print('\nInterpretation aid: fold-change equivalent of each MDE')
for arm in ARMS:
    sub = R[R.arm == arm]
    lo, hi = sub.MDE_log2_80pct.min(), sub.MDE_log2_80pct.max()
    print(f'  {arm:22s} MDE {lo:.2f} to {hi:.2f} log2  '
          f'({2**lo:.2f}x to {2**hi:.2f}x)')
