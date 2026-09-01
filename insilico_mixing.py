"""In silico specimen mixing: how much monocyte contamination does it take to
create a false pathway signal? Purified populations from GSE107011 (Monaco et al.)
are mixed at known fractions to build synthetic specimens with ground truth."""
import numpy as np, pandas as pd, pyreadr
from scipy import stats
rng = np.random.default_rng(20260817)

# ---- reference expression, mapped to symbols ----
d = pd.read_csv('/mnt/user-data/uploads/GSE107011_Processed_data_TPM_txt.gz', sep='\t', index_col=0)
d.index = [i.split('.')[0] for i in d.index]
ann = pyreadr.read_r('/tmp/grch38.rda')['grch38'][['ensgene','symbol']].dropna().drop_duplicates('ensgene')
sym = ann.set_index('ensgene')['symbol']
g = sym.reindex(d.index)
d = d[g.notna().values]; d.index = g.dropna().values
d = d.groupby(level=0).max()
print(f'reference matrix: {d.shape[0]} genes x {d.shape[1]} samples')

def pool(names):
    cols=[c for c in d.columns if c.split('_',1)[1] in names]
    return d[cols].mean(axis=1)

CD4   = pool(['CD4_naive','CD4_TE','Th1','Th2','Th17','Th1/Th17','Treg','TFH'])
MONO  = pool(['C_mono','I_mono','NC_mono'])
BCELL = pool(['B_naive','B_SM','B_NSM','B_Ex'])
NK    = pool(['NK'])
NEUT  = pool(['Neutrophils'])
print('pools built (CD4, monocyte, B, NK, neutrophil)')

# ---- modules and Hallmark sets ----
mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
sets = {m: [x for x in set(v['gene']) if x in d.index] for m, v in mt.groupby('module') if m != 'axis'}
for line in open('code/h.all.v7.0.symbols.gmt'):
    p=line.rstrip('\n').split('\t')
    gs=[x for x in p[2:] if x in d.index]
    if len(gs)>=15: sets[p[0].replace('HALLMARK_','')]=gs

def ssgsea(expr, gs, alpha=0.25):
    genes=expr.index.to_numpy(); R=expr.rank(axis=0,method='average').values; N=expr.shape[0]; out={}
    for name,g_ in gs.items():
        idx=np.isin(genes,g_)
        sc=[]
        for j in range(expr.shape[1]):
            o=np.argsort(-R[:,j]); ins=idx[o]; w=(R[o,j]**alpha)*ins
            sc.append(np.sum(np.cumsum(w)/w.sum()-np.cumsum(~ins)/(N-ins.sum())))
        out[name]=sc
    return pd.DataFrame(out, index=expr.columns).T


# ---- build synthetic specimens from individual donors, preserving real variation ----
# Four donors have every population needed. Each synthetic specimen mixes one
# donor's own CD4 subsets with that donor's own non-T populations, so between-sample
# variance is real donor-to-donor variation rather than an assumed noise model.
T_SETS=['CD4_naive','Th1','Th2','Th17','Treg']
O_SETS={'C_mono':0.40,'I_mono':0.08,'NC_mono':0.07,'B_naive':0.25,'NK':0.12,'Neutrophils':0.08}
donors=sorted({c.split('_',1)[0] for c in d.columns
               if c.split('_',1)[1] in set(T_SETS)|set(O_SETS)})
donors=[dn for dn in donors if all(f'{dn}_{ct}' in d.columns for ct in list(T_SETS)+list(O_SETS))]
print(f'donors usable for per-donor mixing: {len(donors)} {donors}')

FRACTIONS=[0.0,0.02,0.05,0.10,0.15,0.20,0.30]
cols={}; meta=[]
for dn in donors:
    cd4=d[[f'{dn}_{ct}' for ct in T_SETS]].mean(axis=1)
    other=sum(w*d[f'{dn}_{ct}'] for ct,w in O_SETS.items())
    for f in FRACTIONS:
        n=f'case_f{int(f*100)}_{dn}'
        cols[n]=(1-f)*cd4+f*other; meta.append((n,f,'case' if f>0 else 'control'))
X=np.log2(pd.DataFrame(cols)+1)
M=pd.DataFrame(meta,columns=['sample','frac','arm']).set_index('sample')
print(f'synthetic specimens: {X.shape[1]} ({len(donors)} donors x {len(FRACTIONS)} levels)')

S = ssgsea(X, sets)
ctrl=[s for s in M.index if M.loc[s,'frac']==0.0]
rows=[]
for f in [x for x in FRACTIONS if x>0]:
    case=[s for s in M.index if M.loc[s,'frac']==f]
    for name in S.index:
        a,b=S.loc[name,case],S.loc[name,ctrl]
        sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        dd=(a.mean()-b.mean())/sd if sd>0 else np.nan
        p=stats.wilcoxon(a.values,b.values)[1] if len(a)==len(b) else stats.mannwhitneyu(a,b)[1]
        rows.append([name,f,dd,p])
R=pd.DataFrame(rows,columns=['set','frac','cohens_d','p'])
R['q']=np.nan
for f in [x for x in FRACTIONS if x>0]:
    m=R.frac==f
    R.loc[m,'q']=stats.false_discovery_control(R.loc[m,'p'])
R.to_csv('results_insilico_mixing.csv',index=False)

print('\n== magnitude of the artefactual difference, by contamination level ==')
print('   (four donors, so the paired test cannot reach P<0.05; effect sizes are the readout)')
for f in [x for x in FRACTIONS if x>0]:
    s_=R[(R.frac==f)]
    print(f'  {int(f*100):2d}% non-T content: median |d| = {s_.cohens_d.abs().median():.2f}, '
          f'{int((s_.cohens_d.abs()>0.8).sum()):2d}/{len(s_)} sets exceed a large effect (|d|>0.8), '
          f'max |d| = {s_.cohens_d.abs().max():.1f} ({s_.loc[s_.cohens_d.abs().idxmax(),"set"]})')

print('\n== the four autophagy modules ==')
piv=R[R.set.isin(['initiation','elongation','fusion','clear'])].pivot(index='set',columns='frac',values='cohens_d')
print(piv.round(2).to_string())
print('\nfirst contamination level at which each module becomes significant:')
for m in ['initiation','elongation','fusion','clear']:
    s=R[(R.set==m)&(R.q<0.05)].sort_values('frac')
    print(f'  {m:11s} ' + (f'{int(s.frac.iloc[0]*100)}%' if len(s) else 'never'))


# ---- how much contamination reproduces the effect actually observed in GSE33615? ----
OBS={'clear':1.197,'initiation':-0.636,'elongation':0.057,'fusion':-0.169}
print('\n== contamination level that reproduces the GSE33615 observation ==')
for m,obs in OBS.items():
    cur=R[R.set==m].sort_values('frac')
    x=cur.frac.values; y=cur.cohens_d.values
    if (y.max()>=abs(obs)) and obs>0:
        lvl=np.interp(obs,y,x)
        print(f'  {m:11s} observed d={obs:+.2f} -> reproduced by {lvl*100:.1f}% non-T content')
    else:
        print(f'  {m:11s} observed d={obs:+.2f} -> not reproduced by contamination at any level tested')
print('\n  measured difference in GSE33615: monocytes 16.5% vs 6.2%, plus B cell and NK excess')

# ---- which Hallmark sets are most contamination-sensitive? ----
slope=(R[R.frac==0.10].set_index('set').cohens_d).sort_values()
print('\nmost contamination-sensitive Hallmark sets at 10% non-T content:')
print(slope.tail(8).round(2).to_string())
print('\nleast affected:')
print(slope.abs().sort_values().head(8).round(2).to_string())
R.to_csv('results_insilico_mixing.csv',index=False)
