import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':0.8,'savefig.dpi':300})
NPG=['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#8A5A00']
R=pd.read_csv('results_insilico_mixing.csv')
C=pd.read_csv('results_CSI_index.csv')

fig,ax=plt.subplots(1,3,figsize=(7.5,3.0))

# (a) dose-response for the four modules
order=['initiation','elongation','fusion','clear']
for i,m in enumerate(order):
    s=R[R.set==m].sort_values('frac')
    ax[0].plot(s.frac*100,s.cohens_d,marker='o',ms=3.5,lw=1.4,color=NPG[i],label=m)
ax[0].axhline(1.197,color='0.35',ls='--',lw=1)
ax[0].text(29.5,1.30,'effect observed in GSE33615',fontsize=6.5,color='0.35',ha='right')
ax[0].axhline(0,color='0.85',lw=0.7)
ax[0].set_xlabel('Non-T content of specimen (%)')
ax[0].set_ylabel("Difference from pure CD4 (Cohen's $d$)")
ax[0].legend(frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(0.0,1.0))
ax[0].set_xlim(0,31)

# (b) distribution of CSI across 1635 sets
bins=np.linspace(C.CSI.min(),C.CSI.max(),70)
ax[1].hist(C.CSI,bins=bins,color='0.8',edgecolor='none')
for i,m in enumerate(order):
    v=C.loc[C.set=='MODULE_'+m.upper(),'CSI']
    if len(v):
        ax[1].axvline(v.iloc[0],color=NPG[i],lw=1.3)
        if m in ('initiation','clear'):
            ax[1].text(v.iloc[0],230,m,rotation=90,fontsize=6,color=NPG[i],ha='right',va='top')
ax[1].set_yscale('log'); ax[1].set_xlabel('Composition Sensitivity Index')
ax[1].set_ylabel('Number of gene sets')
ax[1].text(0.97,0.95,f'{len(C)} sets',transform=ax[1].transAxes,ha='right',va='top',fontsize=7)

# (c) risk class by collection
piv=C.groupby(['collection','risk']).size().unstack(fill_value=0)
piv=piv.loc[[c for c in ['Hallmark','KEGG','KEGG_MEDICUS','Reactome'] if c in piv.index]]
pct=piv.div(piv.sum(axis=1),axis=0)*100
b=np.zeros(len(pct))
for cls,col in [('robust',NPG[6]),('caution',NPG[4]),('high',NPG[0])]:
    if cls in pct: ax[2].barh(range(len(pct)),pct[cls],left=b,color=col,height=0.62,label=cls); b=b+pct[cls].values
ax[2].set_yticks(range(len(pct))); ax[2].set_yticklabels(pct.index,fontsize=7); ax[2].invert_yaxis()
ax[2].set_xlabel('Percentage of gene sets'); ax[2].set_xlim(0,100)
ax[2].legend(frameon=False,fontsize=6.5,ncol=3,loc='upper center',bbox_to_anchor=(0.5,1.16))
for a,l in zip(ax,['a','b','c']): a.text(-0.24,1.06,l,transform=a.transAxes,fontweight='bold',fontsize=10)
plt.tight_layout(); plt.savefig('fig/Figure6.png',bbox_inches='tight'); plt.close()
print('Figure 6 written')
