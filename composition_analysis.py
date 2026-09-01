"""Is the CLEAR/lysosomal signal in GSE33615 an artifact of cell composition?
No external reference matrix was reachable from this environment, so the test is
internal: genes whose expression tracks a monocyte marker score are identified
within the cohort itself, and the module effect is then re-estimated with lineage
scores as covariates."""
import gzip, numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

exec(open('code/run_gse33615.py').read().split('# ---------- modules ----------')[0]
     .replace('print(', '_ = ('))

mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module') if m != 'axis'}
grp = samples.set_index('gsm').loc[G.columns, 'group']
is_atl = (grp == 'ATL').values

LINEAGE = {
 'monocyte':   ['CD14','LYZ','CSF1R','FCN1','S100A8','S100A9','VCAN','ITGAM','CD68','MNDA'],
 'Bcell':      ['CD19','MS4A1','CD79A','CD79B','PAX5'],
 'NK':         ['NCAM1','KLRD1','GNLY','NKG7','PRF1'],
 'neutrophil': ['ELANE','MPO','CEACAM8','DEFA4','FCGR3B'],
 'platelet':   ['PF4','PPBP','ITGA2B','GP9'],
 'erythroid':  ['HBB','HBA1','ALAS2','SLC4A1'],
 'Tcell':      ['CD3D','CD3E','CD2','LCK','IL7R'],
}
scores = {}
for k, gs in LINEAGE.items():
    idx = [g for g in gs if g in G.index]
    scores[k] = G.loc[idx].apply(stats.zscore, axis=1).mean()
    print(f'{k:11s} n={len(idx):2d}  ATL mean {scores[k][is_atl].mean():+.2f}  '
          f'control {scores[k][~is_atl].mean():+.2f}  p={stats.mannwhitneyu(scores[k][is_atl], scores[k][~is_atl])[1]:.3g}')
Sc = pd.DataFrame(scores)

print('\n=== per-gene test: does the ATL effect track monocyte content? ===')
mono = Sc['monocyte'].values
rows = []
for m, genes in mods.items():
    for g in genes:
        if g not in G.index: continue
        x = G.loc[g]
        rho = stats.spearmanr(x.values, mono)[0]
        d = (x[is_atl].mean() - x[~is_atl].mean()) / np.sqrt(
            (x[is_atl].var(ddof=1) + x[~is_atl].var(ddof=1)) / 2)
        rows.append([m, g, rho, d])
P = pd.DataFrame(rows, columns=['module','gene','rho_monocyte','cohens_d_ATL'])
P.to_csv('results_composition_per_gene.csv', index=False)
for m in mods:
    s = P[P.module == m]
    r, p = stats.pearsonr(s.rho_monocyte, s.cohens_d_ATL)
    print(f'  {m:11s} correlation between monocyte-tracking and ATL effect: r={r:+.3f} (p={p:.3g}, n={len(s)})')
r, p = stats.pearsonr(P.rho_monocyte, P.cohens_d_ATL)
print(f'  ALL {len(P)} module genes: r={r:+.3f} (p={p:.3g})')

print('\n=== module score, adjusted for lineage composition ===')
S = pd.read_csv('results_module_scores.csv', index_col=0)
for m in mods:
    y = S.loc[G.columns, m].values
    X = sm.add_constant(np.column_stack([is_atl.astype(float), Sc.values]))
    fit = sm.OLS(y, X).fit()
    raw = sm.OLS(y, sm.add_constant(is_atl.astype(float))).fit()
    print(f'  {m:11s} unadjusted ATL beta={raw.params[1]:+.4f} (p={raw.pvalues[1]:.3g})  ->  '
          f'adjusted beta={fit.params[1]:+.4f} (p={fit.pvalues[1]:.3g})')
