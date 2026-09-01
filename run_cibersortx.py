"""CIBERSORTx LM22 deconvolution of GSE33615: group tests and composition-adjusted
module effects. Independent confirmation of the MCP-counter analysis."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

C = pd.read_csv('/mnt/user-data/uploads/CIBERSORTx_Job2_Results.csv').set_index('Mixture')
qc = C[['P-value','Correlation','RMSE']]
F = C.drop(columns=['P-value','Correlation','RMSE'])
print(f'LM22 fractions: {F.shape[0]} samples x {F.shape[1]} populations')
print(f'fit quality: all P < 0.001; correlation {qc.Correlation.min():.2f} to {qc.Correlation.max():.2f} '
      f'(median {qc.Correlation.median():.2f}); RMSE median {qc.RMSE.median():.2f}')
print(f'fraction sums: {F.sum(axis=1).min():.3f} to {F.sum(axis=1).max():.3f}')

key = pd.read_csv('results_sample_assignment.csv').set_index('gsm')
grp = key.loc[F.index, 'group']
is_atl = (grp == 'ATL').values
print(f'groups: {int(is_atl.sum())} ATL, {int((~is_atl).sum())} control')

rows = []
for p in F.columns:
    a, b = F.loc[is_atl, p], F.loc[~is_atl, p]
    if F[p].sum() == 0:
        continue
    sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    d = (a.mean() - b.mean()) / sd if sd > 0 else np.nan
    rows.append([p, a.mean(), b.mean(), a.mean()-b.mean(), d, stats.mannwhitneyu(a, b)[1]])
R = pd.DataFrame(rows, columns=['population','mean_ATL','mean_control','difference','cohens_d','p'])
R['q'] = stats.false_discovery_control(R['p'])
R = R.sort_values('q')
print('\n== LM22 population fractions, ATL PBMC vs sorted CD4+ controls ==')
print(R.to_string(index=False, float_format=lambda x: f'{x:.4g}'))
R.to_csv('results_cibersortx_groups.csv', index=False)

# module effects adjusted for LM22 composition principal components
S = pd.read_csv('results_module_scores.csv', index_col=0)
use = [p for p in F.columns if F[p].std() > 0]
Z = (F[use] - F[use].mean()) / F[use].std()
u, s, vt = np.linalg.svd(Z.values - Z.values.mean(0), full_matrices=False)
var = s**2 / (s**2).sum()
npc = max(2, int(np.argmax(np.cumsum(var) > 0.80) + 1))
PC = u[:, :npc] * s[:npc]
print(f'\ncomposition PCs retained: {npc} (cumulative variance {np.cumsum(var)[npc-1]:.2f})')

out = []
for m in ['initiation','elongation','fusion','clear']:
    y = S.loc[F.index, m].values
    raw = sm.OLS(y, sm.add_constant(is_atl.astype(float))).fit()
    adj = sm.OLS(y, sm.add_constant(np.column_stack([is_atl.astype(float), PC]))).fit()
    out.append([m, raw.params[1], raw.pvalues[1], adj.params[1], adj.pvalues[1]])
O = pd.DataFrame(out, columns=['module','beta_unadj','p_unadj','beta_adj_LM22','p_adj_LM22'])
print('\n== module effect, unadjusted vs adjusted for LM22 composition PCs ==')
print(O.to_string(index=False, float_format=lambda x: f'{x:.4g}'))
O.to_csv('results_module_adjusted_cibersortx.csv', index=False)

# per-gene: does monocyte fraction predict the apparent ATL effect?
exec(open('code/run_gse33615.py').read().split('# ---------- modules ----------')[0].replace('print(','_p=('))
mono = F.loc[G.columns, 'Monocytes'].values
mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
rows = []
for _, r in mt[mt.module != 'axis'].iterrows():
    if r.gene not in G.index: continue
    x = G.loc[r.gene]
    rho = stats.spearmanr(x.values, mono)[0]
    d = (x[is_atl].mean() - x[~is_atl].mean()) / np.sqrt(
        (x[is_atl].var(ddof=1) + x[~is_atl].var(ddof=1)) / 2)
    rows.append([r.module, r.gene, rho, d])
P = pd.DataFrame(rows, columns=['module','gene','rho_monocyte_LM22','cohens_d'])
rr, pp = stats.pearsonr(P.rho_monocyte_LM22, P.cohens_d)
print(f'\nLM22 monocyte fraction vs apparent ATL effect, {len(P)} genes: r={rr:+.3f} p={pp:.3g}')
P.to_csv('results_cibersortx_per_gene.csv', index=False)

# agreement between the two methods
E = pd.read_csv('results_mcpcounter_estimates.csv', index_col=0)
pairs = [('Monocytes','Monocytic lineage'), ('T cells CD4 memory resting','T cells'),
         ('B cells naive','B lineage'), ('NK cells activated','NK cells'), ('Neutrophils','Neutrophils')]
print('\n== agreement between LM22 fractions and MCP-counter estimates ==')
for lm, mc in pairs:
    if lm in F.columns and mc in E.columns and F[lm].std() > 0:
        r_, p_ = stats.spearmanr(F.loc[E.index, lm], E[mc])
        print(f'  {lm:28s} vs {mc:20s} rho={r_:+.3f} (p={p_:.2g})')
