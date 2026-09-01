"""MCP-counter deconvolution of GSE33615 using the published signature set.
Becht et al., Genome Biology 2016. Population abundance is the arithmetic mean of
log2 expression across the transcriptomic markers of that population, which is the
published estimator."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

exec(open('code/run_gse33615.py').read().split('# ---------- modules ----------')[0].replace('print(', '_p=('))

sig = pd.read_csv('code/mcpcounter_signatures.txt', sep='\t')
sig.columns = [c.strip('"') for c in sig.columns]
sig = sig.rename(columns={'HUGO symbols': 'gene', 'Cell population': 'pop'})
sig['gene'] = sig['gene'].str.strip('"'); sig['pop'] = sig['pop'].str.strip('"')

grp = samples.set_index('gsm').loc[G.columns, 'group']
is_atl = (grp == 'ATL').values

est, cover = {}, {}
for pop, blk in sig.groupby('pop'):
    idx = [g for g in blk['gene'].unique() if g in G.index]
    cover[pop] = (len(idx), blk['gene'].nunique())
    if len(idx) >= 3:
        est[pop] = G.loc[idx].mean()
E = pd.DataFrame(est)
print('MCP-counter populations estimated (markers found / total):')
for p, (a, b) in sorted(cover.items()):
    print(f'  {p:28s} {a:3d}/{b:3d}' + ('' if a >= 3 else '   excluded'))

rows = []
for p in E.columns:
    a, b = E.loc[is_atl, p], E.loc[~is_atl, p]
    d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    rows.append([p, a.mean(), b.mean(), a.mean() - b.mean(), d,
                 stats.mannwhitneyu(a, b)[1]])
R = pd.DataFrame(rows, columns=['population', 'mean_ATL', 'mean_control',
                                'difference', 'cohens_d', 'p'])
R['q'] = stats.false_discovery_control(R['p'])
R = R.sort_values('q')
print('\n== population abundance, ATL PBMC vs sorted CD4+ controls ==')
print(R.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
R.to_csv('results_mcpcounter_groups.csv', index=False)
E.to_csv('results_mcpcounter_estimates.csv')

# module effect before and after adjusting for composition principal components
S = pd.read_csv('results_module_scores.csv', index_col=0)
Z = (E - E.mean()) / E.std()
u, s, vt = np.linalg.svd(Z.values - Z.values.mean(0), full_matrices=False)
var = s**2 / (s**2).sum()
npc = max(2, int(np.argmax(np.cumsum(var) > 0.80) + 1))
PC = u[:, :npc] * s[:npc]
print(f'\ncomposition PCs retained: {npc} (cumulative variance {np.cumsum(var)[npc-1]:.2f})')

out = []
for m in ['initiation', 'elongation', 'fusion', 'clear']:
    y = S.loc[G.columns, m].values
    raw = sm.OLS(y, sm.add_constant(is_atl.astype(float))).fit()
    adj = sm.OLS(y, sm.add_constant(np.column_stack([is_atl.astype(float), PC]))).fit()
    out.append([m, raw.params[1], raw.pvalues[1], adj.params[1], adj.pvalues[1]])
O = pd.DataFrame(out, columns=['module', 'beta_unadj', 'p_unadj', 'beta_adj_PC', 'p_adj_PC'])
print('\n== module effect, unadjusted vs adjusted for composition PCs ==')
print(O.to_string(index=False, float_format=lambda x: f'{x:.4g}'))
O.to_csv('results_module_adjusted_mcp.csv', index=False)

# does each module gene's association with monocyte abundance predict its ATL effect?
mono = E['Monocytic lineage'] if 'Monocytic lineage' in E.columns else None
if mono is not None:
    mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
    rows = []
    for _, r in mt[mt.module != 'axis'].iterrows():
        if r.gene not in G.index: continue
        x = G.loc[r.gene]
        rho = stats.spearmanr(x.values, mono.values)[0]
        d = (x[is_atl].mean() - x[~is_atl].mean()) / np.sqrt(
            (x[is_atl].var(ddof=1) + x[~is_atl].var(ddof=1)) / 2)
        rows.append([r.module, r.gene, rho, d])
    P = pd.DataFrame(rows, columns=['module', 'gene', 'rho_mono_MCP', 'cohens_d'])
    r, p = stats.pearsonr(P.rho_mono_MCP, P.cohens_d)
    print(f'\nMCP monocytic abundance vs apparent ATL effect, {len(P)} genes: r={r:+.3f} p={p:.3g}')
    P.to_csv('results_mcp_per_gene.csv', index=False)
