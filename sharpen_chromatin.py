"""Sharpen the single-library chromatin arms three ways:
1. a baseline-matched null (module promoters are highly marked at baseline, so an
   unmatched random draw is the wrong comparison set)
2. bootstrap confidence intervals on the module mean
3. combination across the four independent Tax contrasts, which are separate
   experiments even though none is internally replicated"""
import numpy as np, pandas as pd, pyBigWig
from scipy import stats
rng = np.random.default_rng(20260817)

mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module') if m != 'axis'}
NFKB = ['RELB','NFKB2','NFKB1','RELA','REL','NFKBIA','TNFAIP3','BIRC3','CD40','ICAM1',
        'TRAF1','CCL5','IL2RA','TNFRSF8','CFLAR','BCL2A1','PLEK','LTA']

CON = {'H3K27ac Jurkat': ('results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv',
                          'bw2/GSM9453007_CT_Jurkat_Empty_H3K27ac.bw', 2000),
       'H3K27ac TL-Om1': ('results_h3k27ac_lfc_CT_TL-Om1_Tax_H3K27ac.csv',
                          'bw2/GSM9453009_CT_TL-Om1_Empty_H3K27ac.bw', 2000),
       'ATAC empty':     ('results_atac_lfc_Empty_Tax.csv',
                          'atac/GSM9453003_ATAC_JPX-9_Empty.bw', 1000),
       'ATAC IRF4':      ('results_atac_lfc_IRF4_Tax.csv',
                          'atac/GSM9453004_ATAC_JPX-9_IRF4.bw', 1000)}
reg = pd.read_csv('regions_hg38.tsv', sep='\t')

def baseline(path, win):
    bw = pyBigWig.open(path); ch = bw.chroms(); out={}
    for r in reg.itertuples():
        if r.chrom not in ch: continue
        s,e = max(0,r.tss-win), min(ch[r.chrom], r.tss+win)
        try: v = bw.stats(r.chrom,s,e,type='mean')[0]
        except Exception: v=None
        out[r.symbol]=0.0 if v is None else float(v)
    bw.close(); return pd.Series(out)

def matched_null(lfc, base, genes, n=20000):
    """draw random sets matched on baseline-signal decile"""
    idx=[g for g in genes if g in lfc.index and g in base.index]
    dec = pd.qcut(base.rank(method='first'), 10, labels=False)
    want = dec[idx].value_counts()
    pools = {d: dec.index[dec==d].intersection(lfc.index) for d in want.index}
    obs = lfc[idx].mean()
    null=np.empty(n)
    for i in range(n):
        pick=[]
        for d,k in want.items():
            pick += list(rng.choice(pools[d], k, replace=False))
        null[i]=lfc[pick].mean()
    p=max(2*min((null>=obs).mean(),(null<=obs).mean()),1/n)
    return obs, p, (obs-null.mean())/null.std(), len(idx)

def boot_ci(lfc, genes, n=10000):
    idx=[g for g in genes if g in lfc.index]
    v=lfc[idx].values
    bs=np.array([rng.choice(v,len(v),replace=True).mean() for _ in range(n)])
    return np.percentile(bs,[2.5,97.5])

rows=[]
for name,(f,bwf,win) in CON.items():
    lfc = pd.read_csv(f, index_col=0).iloc[:,0]
    base = baseline(bwf, win)
    print(f'\n== {name} ==  baseline percentile of module promoters:')
    for m,genes in list(mods.items())+[('NF-kB benchmark',NFKB)]:
        idx=[g for g in genes if g in base.index]
        print(f'   {m:16s} median baseline percentile {base.rank(pct=True)[idx].median():.2f}')
    print('   baseline-matched competitive test:')
    for m,genes in list(mods.items())+[('NF-kB benchmark',NFKB)]:
        obs,p,z,k = matched_null(lfc, base, genes)
        lo,hi = boot_ci(lfc, genes)
        print(f'   {m:16s} mean={obs:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  z={z:+.2f}  p_matched={p:.4f}  (n={k})')
        rows.append([name,m,obs,lo,hi,z,p])
R=pd.DataFrame(rows,columns=['contrast','set','mean_log2','ci_lo','ci_hi','z_matched','p_matched'])
R.to_csv('results_chromatin_matched.csv',index=False)

print('\n== combination across the four independent Tax contrasts ==')
for m in list(mods)+['NF-kB benchmark']:
    sub=R[R.set==m]
    zs=sub.z_matched.values
    zc=zs.sum()/np.sqrt(len(zs))
    pc=2*(1-stats.norm.cdf(abs(zc)))
    print(f'   {m:16s} z per contrast {np.round(zs,2)}  combined z={zc:+.2f}  P={pc:.3f}')
