"""Which pathway signatures in GSE33615 are driven by specimen composition?
Every Hallmark gene set is scored, tested for ATL versus control, then re-tested
with LM22 composition principal components as covariates."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

exec(open('code/run_gse33615.py').read().split('# ---------- modules ----------')[0].replace('print(','_p=('))
grp = samples.set_index('gsm').loc[G.columns,'group']; is_atl=(grp=='ATL').values

C=pd.read_csv('/mnt/user-data/uploads/CIBERSORTx_Job2_Results.csv').set_index('Mixture')
F=C.drop(columns=['P-value','Correlation','RMSE']).loc[G.columns]
mono=F['Monocytes'].values
_e=[g for g in ['HBB','HBA1','HBA2','ALAS2','SLC4A1','AHSP'] if g in G.index]
ery=np.asarray(pd.DataFrame(stats.zscore(G.loc[_e].values,axis=1),columns=G.columns).mean())
use=[c for c in F.columns if F[c].std()>0]
Z=(F[use]-F[use].mean())/F[use].std()
u_,s_,vt_=np.linalg.svd(Z.values-Z.values.mean(0),full_matrices=False)
PC=u_[:,:11]*s_[:11]

sets={}
for line in open('code/h.all.v7.0.symbols.gmt'):
    p=line.rstrip('\n').split('\t'); sets[p[0].replace('HALLMARK_','')]=[g for g in p[2:] if g in G.index]
sets={k:v for k,v in sets.items() if len(v)>=15}
print(f'{len(sets)} Hallmark sets with at least 15 genes on this array')

def ssgsea(expr, gs, alpha=0.25):
    genes=expr.index.to_numpy(); R=expr.rank(axis=0,method='average').values; N=expr.shape[0]; out={}
    for name,g in gs.items():
        idx=np.isin(genes,g)
        sc=[]
        for j in range(expr.shape[1]):
            o=np.argsort(-R[:,j]); ins=idx[o]; w=(R[o,j]**alpha)*ins
            sc.append(np.sum(np.cumsum(w)/w.sum()-np.cumsum(~ins)/(N-ins.sum())))
        out[name]=sc
    S=pd.DataFrame(out,index=expr.columns).T
    return (S-S.values.min())/(S.values.max()-S.values.min())

S=ssgsea(G,sets)
rows=[]
for name,genes in sets.items():
    y=S.loc[name,G.columns].values
    raw=sm.OLS(y,sm.add_constant(is_atl.astype(float))).fit()
    adj=sm.OLS(y,sm.add_constant(np.column_stack([is_atl.astype(float),PC]))).fit()
    ery_adj=sm.OLS(y,sm.add_constant(np.column_stack([is_atl.astype(float),ery,mono]))).fit()
    rho=[abs(stats.spearmanr(G.loc[g].values,mono)[0]) for g in genes]
    rhoe=[abs(stats.spearmanr(G.loc[g].values,ery)[0]) for g in genes]
    rows.append([name,len(genes),raw.params[1],raw.pvalues[1],adj.pvalues[1],ery_adj.pvalues[1],
                 np.mean(np.array(rho)>0.3),np.mean(np.array(rhoe)>0.3)])
R=pd.DataFrame(rows,columns=['pathway','n_genes','beta','p_unadj','p_adj_LM22','p_adj_ery_mono',
                             'frac_mono_linked','frac_ery_linked'])
R['q_unadj']=stats.false_discovery_control(R.p_unadj)
R['q_adj_LM22']=stats.false_discovery_control(R.p_adj_LM22)
sig=R[R.q_unadj<0.05]
lost=sig[stats.false_discovery_control(sig.p_adj_LM22)>=0.05]
print(f'\nHallmark sets significant before adjustment: {len(sig)}/{len(R)}')
print(f'losing significance after composition adjustment: {len(lost)}/{len(sig)} ({100*len(lost)/len(sig):.0f}%)')
R=R.sort_values('frac_ery_linked',ascending=False)
R.to_csv('results_pathway_survey.csv',index=False)
print('\nmost composition-linked pathways (fraction of genes tracking erythroid content):')
print(R.head(12)[['pathway','n_genes','frac_ery_linked','frac_mono_linked','q_unadj','q_adj_LM22']].to_string(index=False,float_format=lambda x:f'{x:.3g}'))
print('\nleast composition-linked:')
print(R.tail(8)[['pathway','n_genes','frac_ery_linked','frac_mono_linked','q_unadj','q_adj_LM22']].to_string(index=False,float_format=lambda x:f'{x:.3g}'))
print(f'\nbackground fraction of all array genes tracking erythroid content at |rho|>0.3: '
      f'{np.mean([abs(stats.spearmanr(G.iloc[i].values,ery)[0])>0.3 for i in np.random.default_rng(1).choice(G.shape[0],2000,replace=False)]):.2f}')

# ---- Supplementary Figure S1: the full Hallmark survey ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'savefig.dpi':300})
NPG=['#E64B35','#4DBBD5','#00A087','#3C5488','#8A5A00']
D=R.sort_values('frac_ery_linked')
fig,ax=plt.subplots(1,2,figsize=(7.5,6.4),gridspec_kw={'width_ratios':[2.6,1]})
y=np.arange(len(D))
col=[NPG[0] if q<0.05 else '0.75' for q in D.q_unadj]
ax[0].barh(y,D.frac_ery_linked,color=col,height=0.72)
ax[0].axvline(0.29,color='0.25',ls='--',lw=1)
ax[0].set_yticks(y); ax[0].set_yticklabels([p.replace('_',' ').title() for p in D.pathway],fontsize=6.4)
ax[0].set_xlabel('Fraction of genes tracking erythroid content ($|\\rho|$ > 0.3)')
ax[0].text(0.295,len(D)-0.3,'background',fontsize=6,color='0.25',va='bottom')
ax[0].set_ylim(-1,len(D)+1.2); ax[0].set_xlim(0,0.60)
from matplotlib.patches import Patch
ax[0].legend(handles=[Patch(color=NPG[0],label='significant before adjustment'),
                      Patch(color='0.75',label='not significant')],frameon=False,fontsize=6,loc='lower right',bbox_to_anchor=(1.0,0.02))
lost=[(q1<0.05) and (q2>=0.05) for q1,q2 in zip(D.q_unadj,D.q_adj_LM22)]  # FDR-based, matches the main text
ax[1].scatter(-np.log10(D.q_unadj+1e-12),-np.log10(D.q_adj_LM22+1e-12),
              s=18,c=[NPG[4] if l else '0.6' for l in lost],edgecolor='white',linewidth=0.4)
lim=max(-np.log10(D.q_unadj+1e-12).max(),-np.log10(D.q_adj_LM22+1e-12).max())*1.05
ax[1].plot([0,lim],[0,lim],color='0.7',lw=0.8,ls=':')
ax[1].axhline(-np.log10(0.05),color='0.4',lw=0.8); ax[1].axvline(-np.log10(0.05),color='0.4',lw=0.8)
ax[1].set_xlabel('$-\\log_{10}q$ unadjusted'); ax[1].set_ylabel('$-\\log_{10}q$ after composition adjustment')
ax[1].text(0.97,0.04,f'{sum(lost)} sets lose\nsignificance',transform=ax[1].transAxes,va='bottom',ha='right',fontsize=6.8,color=NPG[4])
for a,l in zip(ax,['a','b']): a.text(-0.12,1.02,l,transform=a.transAxes,fontweight='bold',fontsize=10)
plt.tight_layout(); plt.savefig('fig/FigureS1.png',bbox_inches='tight'); plt.close()
print('Supplementary Figure S1 written')
