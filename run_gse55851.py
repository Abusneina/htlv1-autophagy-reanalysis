import gzip, numpy as np, pandas as pd
from scipy import stats
SM='/mnt/user-data/uploads/GSE55851_series_matrix_txt.gz'
SOFT='/mnt/user-data/uploads/GSE55851_family_soft.gz'
MIN_COVERAGE=0.60; FDR=0.05

meta={}
with gzip.open(SM,'rt',errors='replace') as f:
    for i,line in enumerate(f):
        if line.startswith('!series_matrix_table_begin'): skip=i+1; break
        if line.startswith('!Sample_'):
            k=line.split('\t')[0][1:]
            meta.setdefault(k,[]).append([x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]])
gsm=meta['Sample_geo_accession'][0]; title=meta['Sample_title'][0]
sub=[b for b in meta['Sample_characteristics_ch1'] if 'atl subtype' in b[0].lower()][0]
sub=[s.split(':',1)[1].strip() for s in sub]
sm=pd.DataFrame({'gsm':gsm,'title':title,'subtype':sub})
sm['patient']=sm['title'].str.rsplit('-',n=1).str[0]
sm['fraction']=sm['title'].str.rsplit('-',n=1).str[1]
sm['class']=np.where(sm['subtype'].str.contains('Normal'),'Normal',
             np.where(sm['subtype'].str.contains('carrier'),'AC','ATL'))
print(sm.to_string(index=False)); print()

X=pd.read_csv(SM,sep='\t',skiprows=skip,compression='gzip',index_col=0,comment='!',low_memory=False)
X=X.apply(pd.to_numeric,errors='coerce').loc[:,sm['gsm']]
rows=[]
with gzip.open(SOFT,'rt',errors='replace') as f:
    inp=False
    for line in f:
        if line.startswith('!platform_table_begin'):
            cols=next(f).rstrip('\n').split('\t'); inp=True; continue
        if line.startswith('!platform_table_end'): break
        if inp: rows.append(line.rstrip('\n').split('\t'))
ann=pd.DataFrame(rows,columns=cols)
ann['ID']=pd.to_numeric(ann['ID'],errors='coerce')
ann=ann.dropna(subset=['ID']).set_index('ID'); ann.index=ann.index.astype(int)
X=np.log2(X.clip(lower=1))
ranks=X.rank(method='first'); ms=np.sort(X.values,axis=0).mean(axis=1)
X=pd.DataFrame(np.interp(ranks.values,np.arange(1,X.shape[0]+1),ms),index=X.index,columns=X.columns)
sym=ann.reindex(X.index)['GENE_SYMBOL'].fillna('')
ctrl=ann.reindex(X.index).get('CONTROL_TYPE',pd.Series('',index=X.index)).fillna('')
keep=(sym!='')&(ctrl.astype(str).str.lower().isin(['','false']))
X,sym=X[keep],sym[keep]
o=X.mean(axis=1).sort_values(ascending=False).index; X,sym=X.loc[o],sym.loc[o]
d=~sym.duplicated(); G=X[d]; G.index=sym[d]
print(f'gene matrix: {G.shape}')

mt=pd.read_csv('code/modules_frozen_v1.tsv',sep='\t')
mods={m:sorted(set(g['gene'])) for m,g in mt.groupby('module')}
axis=mods.pop('axis')
print('coverage:', {m:round(len([g for g in gs if g in G.index])/len(gs),2) for m,gs in list(mods.items())+[('axis',axis)]})

def ssgsea(expr, gene_sets, alpha=0.25):
    genes=expr.index.to_numpy(); R=expr.rank(axis=0,method='average').values; N=expr.shape[0]; out={}
    for name,gs in gene_sets.items():
        idx=np.isin(genes,gs)
        if idx.sum()<5: continue
        sc=[]
        for j in range(expr.shape[1]):
            order=np.argsort(-R[:,j]); ins=idx[order]; w=(R[order,j]**alpha)*ins
            sc.append(np.sum(np.cumsum(w)/w.sum()-np.cumsum(~ins)/(N-ins.sum())))
        out[name]=sc
    S=pd.DataFrame(out,index=expr.columns).T
    return (S-S.values.min())/(S.values.max()-S.values.min())

S=ssgsea(G,mods); S.columns=sm['title'].values; G.columns=sm['title'].values
print('\n== module scores by sample ==')
print(S.T.join(sm.set_index('title')[['class','fraction']]).to_string(float_format=lambda x:f'{x:.3f}'))

atlN=sm[(sm['class']=='ATL')&(sm['fraction']=='N')]['title']
norm=sm[sm['class']=='Normal']['title']
print('\n== ATL tumor fraction (N, n=%d) vs Normal CD4 (P, n=%d), all FACS-sorted ==' % (len(atlN),len(norm)))
for m in S.index:
    a,b=S.loc[m,atlN],S.loc[m,norm]
    p=stats.mannwhitneyu(a,b,alternative='two-sided')[1]
    print(f'  {m:11s} ATL-N={a.mean():.3f}  Normal={b.mean():.3f}  diff={a.mean()-b.mean():+.3f}  p={p:.3g}')
for g in axis+['CD14','LYZ']:
    if g in G.index:
        a,b=G.loc[g,atlN],G.loc[g,norm]
        print(f'  {g:8s} log2FC={a.mean()-b.mean():+.2f}  p={stats.mannwhitneyu(a,b)[1]:.3g}')

print('\n== within-patient paired: tumor fraction N vs normal-like fraction P ==')
pairs=[p for p in sm['patient'].unique() if {'N','P'}<=set(sm[sm.patient==p]['fraction'])]
for m in S.index:
    dif=[S.loc[m,f'{p}-N']-S.loc[m,f'{p}-P'] for p in pairs]
    print(f'  {m:11s} mean paired diff={np.mean(dif):+.3f}  ' + '  '.join(f'{p}:{d:+.3f}' for p,d in zip(pairs,dif)))
for g in ['SIRT1','EP300','CREBBP']:
    if g in G.index:
        dif=[G.loc[g,f'{p}-N']-G.loc[g,f'{p}-P'] for p in pairs]
        print(f'  {g:11s} mean paired log2FC={np.mean(dif):+.2f}  ' + '  '.join(f'{p}:{d:+.2f}' for p,d in zip(pairs,dif)))
S.T.join(sm.set_index('title')).to_csv('results_gse55851_scores.csv')
