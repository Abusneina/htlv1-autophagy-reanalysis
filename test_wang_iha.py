"""Do the published ATL 'ferroptosis' and 'autophagy' DEG lists (Wang & Iha,
Genes 2023;14:2005) track specimen composition, and do they replicate in
FACS-sorted cells?"""
import json, numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

D = json.load(open('wang_iha_degs.json'))
exec(open('code/run_gse33615.py').read().split('# ---------- modules ----------')[0].replace('print(','_p=('))
grp = samples.set_index('gsm').loc[G.columns,'group']; is_atl=(grp=='ATL').values

C = pd.read_csv('/mnt/user-data/uploads/CIBERSORTx_Job2_Results.csv').set_index('Mixture')
F = C.drop(columns=['P-value','Correlation','RMSE']).loc[G.columns]
mono = F['Monocytes'].values
_e = [g for g in ['HBB','HBA1','HBA2','ALAS2','SLC4A1','AHSP'] if g in G.index]
ery = np.asarray(pd.DataFrame(stats.zscore(G.loc[_e].values, axis=1), columns=G.columns).mean())
print('erythroid markers used:', _e)

def profile(name, genes):
    genes=[g for g in genes if g in G.index]
    rows=[]
    for g in genes:
        x=G.loc[g].values
        rows.append([g, stats.spearmanr(x,mono)[0], stats.spearmanr(x,ery)[0]])
    P=pd.DataFrame(rows,columns=['gene','rho_mono','rho_ery'])
    # background: all genes on the array
    bg_m=[]; bg_e=[]
    rng=np.random.default_rng(20260817)
    idx=rng.choice(G.shape[0], 4000, replace=False)
    for i in idx:
        x=G.iloc[i].values
        bg_m.append(stats.spearmanr(x,mono)[0]); bg_e.append(stats.spearmanr(x,ery)[0])
    bg_m=np.array(bg_m); bg_e=np.array(bg_e)
    print(f'\n== {name} ({len(genes)} genes recovered on the array) ==')
    for lab,col,bg in [('monocyte','rho_mono',bg_m),('erythroid','rho_ery',bg_e)]:
        frac=(P[col].abs()>0.3).mean(); bgfrac=(np.abs(bg)>0.3).mean()
        u,p=stats.mannwhitneyu(P[col].abs(), np.abs(bg))
        print(f'   |rho| with {lab:9s} content: median {P[col].abs().median():.3f} vs background {np.median(np.abs(bg)):.3f}'
              f'   |rho|>0.3 in {frac*100:.0f}% vs {bgfrac*100:.0f}%   P={p:.2g}')
    return P

Pf = profile('Ferroptosis-related DEGs', D['ferroptosis'])
Pa = profile('Autophagy-related DEGs',  D['autophagy'])

# --- how many survive composition adjustment ---
use=[c for c in F.columns if F[c].std()>0]
Z=(F[use]-F[use].mean())/F[use].std()
u_,s_,vt_=np.linalg.svd(Z.values-Z.values.mean(0),full_matrices=False)
PC=u_[:,:11]*s_[:11]
def survive(name, genes):
    genes=[g for g in genes if g in G.index]
    raw_p=[]; adj_p=[]
    for g in genes:
        y=G.loc[g].values
        raw=sm.OLS(y,sm.add_constant(is_atl.astype(float))).fit()
        adj=sm.OLS(y,sm.add_constant(np.column_stack([is_atl.astype(float),PC]))).fit()
        raw_p.append(raw.pvalues[1]); adj_p.append(adj.pvalues[1])
    raw_p=np.array(raw_p); adj_p=np.array(adj_p)
    rq=stats.false_discovery_control(raw_p); aq=stats.false_discovery_control(adj_p)
    print(f'\n{name}: significant unadjusted {int((rq<0.05).sum())}/{len(genes)}, '
          f'after composition adjustment {int((aq<0.05).sum())}/{len(genes)} '
          f'({100*(1-(aq<0.05).sum()/max((rq<0.05).sum(),1)):.0f}% lost)')
    return pd.DataFrame({'gene':genes,'q_unadj':rq,'q_adj':aq})
Sf=survive('Ferroptosis DEGs', D['ferroptosis'])
Sa=survive('Autophagy DEGs',  D['autophagy'])

# --- replication in the sorted cohort ---

# --- LM22 has no erythroid category; adjust on the erythroid score explicitly ---
def survive_ery(name, genes):
    genes=[g for g in genes if g in G.index]
    raw=[];adj=[]
    for g in genes:
        y=G.loc[g].values
        r=sm.OLS(y,sm.add_constant(is_atl.astype(float))).fit()
        a=sm.OLS(y,sm.add_constant(np.column_stack([is_atl.astype(float),ery,mono]))).fit()
        raw.append(r.pvalues[1]); adj.append(a.pvalues[1])
    rq=stats.false_discovery_control(np.array(raw)); aq=stats.false_discovery_control(np.array(adj))
    n0=int((rq<0.05).sum()); n1=int((aq<0.05).sum())
    print(f'{name}: adjusting on erythroid + monocyte content -> {n1}/{len(genes)} remain '
          f'(was {n0}/{len(genes)}; {100*(n0-n1)/max(n0,1):.0f}% lost)')
    return pd.DataFrame({'gene':genes,'q_unadj':rq,'q_adj_ery_mono':aq})
print()
Ef=survive_ery('Ferroptosis DEGs', D['ferroptosis'])
Ea=survive_ery('Autophagy DEGs',  D['autophagy'])
Ef.to_csv('results_wangiha_ferroptosis.csv',index=False); Ea.to_csv('results_wangiha_autophagy.csv',index=False)
print('\ngenes losing significance (autophagy):', list(Ea[(Ea.q_unadj<0.05)&(Ea.q_adj_ery_mono>=0.05)].gene))
print('genes losing significance (ferroptosis):', list(Ef[(Ef.q_unadj<0.05)&(Ef.q_adj_ery_mono>=0.05)].gene))
