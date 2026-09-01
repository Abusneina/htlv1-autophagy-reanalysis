import pandas as pd, numpy as np
from scipy import stats
S = pd.read_csv('results_module_scores.csv', index_col=0)
exec(open('code/run_gse33615.py').read().split("# ---------- modules ----------")[0].replace("print(","_p=(" ))
mono = ['CD14','LYZ','CSF1R','FCN1','S100A8','S100A9','VCAN','ITGAM']
tcell = ['CD3D','CD3E','CD2','IL7R','LCK']
prol = ['MKI67','TOP2A','CCNB1','PCNA']
grp = samples.set_index('gsm').loc[G.columns,'group']
for name,genes in [('monocyte',mono),('T cell',tcell),('proliferation',prol)]:
    g = [x for x in genes if x in G.index]
    sc = G.loc[g].apply(stats.zscore, axis=1).mean()
    a,b = sc[(grp=='ATL').values], sc[(grp=='control').values]
    print(f'{name:14s} n={len(g)}  ATL={a.mean():+.2f} ctrl={b.mean():+.2f}  p={stats.mannwhitneyu(a,b)[1]:.2g}')
    r,p = stats.spearmanr(np.asarray(sc), S.loc[G.columns,'clear'].values)
    print(f'{"":14s} rho with CLEAR score = {r:+.3f} (p={p:.2g})')
    globals()[name.split()[0]+"_sc"] = pd.Series(np.asarray(sc), index=G.columns)
# partial: does CLEAR still separate after regressing out monocyte score?
import numpy.linalg as la
y = S.loc[G.columns,'clear'].values
x = np.asarray(monocyte_sc)
resid = y - np.polyval(np.polyfit(x,y,1), x)
a,b = resid[(grp=='ATL').values], resid[(grp=='control').values]
print(f'\nCLEAR residual after removing monocyte signal: ATL={a.mean():+.3f} ctrl={b.mean():+.3f} p={stats.mannwhitneyu(a,b)[1]:.2g}')
ys = S.loc[G.columns,'initiation'].values
resid2 = ys - np.polyval(np.polyfit(x,ys,1), x)
a2,b2 = resid2[(grp=='ATL').values], resid2[(grp=='control').values]
print(f'initiation residual: ATL={a2.mean():+.3f} ctrl={b2.mean():+.3f} p={stats.mannwhitneyu(a2,b2)[1]:.2g}')
s1 = G.loc["SIRT1"]
r,p = stats.spearmanr(np.asarray(s1), x); print(f'\nSIRT1 vs monocyte score rho={r:+.3f} p={p:.2g}')
