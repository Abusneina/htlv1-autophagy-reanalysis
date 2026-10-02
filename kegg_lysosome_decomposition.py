"""Why is KEGG_LYSOSOME the most composition-sensitive autophagy-adjacent set?

1. Recompute the set-level index with the KEGG legacy file supplied by the user
   (MSigDB 2026.1) and compare with the published value, to detect version drift.
2. Decompose the set gene by gene using results_gene_level_CSI.csv.
3. Score every member's association with the monocyte marker score in the unsorted
   GSE33615 cohort, using the same loader and score as composition_analysis.py, and
   test whether the genes contamination moves most in the mixtures are the ones that
   track monocyte content in the cohort.
"""
import os
os.makedirs('results', exist_ok=True)
import sys
# Input: the MSigDB KEGG legacy symbols file, downloaded by the user (licensed, not committed),
# placed in this folder, or pass its path as the first argument.
GMT = sys.argv[1] if len(sys.argv) > 1 else 'c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt'
import numpy as np, pandas as pd, re
from scipy import stats

GMT = GMT
lyso = None
for line in open(GMT):
    p = line.rstrip('\n').split('\t')
    if p[0] == 'KEGG_LYSOSOME':
        lyso = sorted(set(p[2:]))
print(f'KEGG_LYSOSOME in MSigDB 2026.1: {len(lyso)} genes', flush=True)

G = pd.read_csv('results/results_gene_level_CSI.csv', index_col=0)
L = G.reindex(lyso).dropna(subset=['gene_CSI']).sort_values('gene_CSI', ascending=False)
print(f'  expressed in the sorted reference: {len(L)}', flush=True)

# ---- 1. set-level index recomputed with this file version ----
src = open('gene_level_csi.py').read()
head = src.split('# ---------- per-gene slope ----------')[0]
ns = {'print': lambda *a, **k: None, 'sys': type('S', (), {'argv': ['x']})}
exec(head, ns)
gs = [g for g in lyso if g in ns['X'].index]
val = ns['slope'](ns['score'](gs))
pub = pd.read_csv('results/results_CSI_index.csv').set_index('set').loc['KEGG_LYSOSOME']
print(f'  set-level index recomputed {val:+.4f} on {len(gs)} genes; published {pub.CSI:+.4f} '
      f'on {int(pub.n_genes)} genes', flush=True)

# ---- 2. family breakdown ----
def family(g):
    if g.startswith('CTS'): return 'cathepsins'
    if g.startswith('ATP6'): return 'vacuolar ATPase'
    if re.match(r'^(HEX|GLA|GLB|GBA|GAA|GUSB|GNS|GALC|GALNS|ARS|IDS|IDUA|NAGA|NAGLU|NEU|MAN2B|MANBA|FUCA|ASAH|SMPD|SGSH|AGA|LIPA|PPT|TPP|PSAP|NPC|SCARB|HGSNAT|PLA2G15|ACP)', g):
        return 'other hydrolases and lipid handling'
    if re.match(r'^(LAMP|CD63|CD68|LAPTM|MCOLN|SLC|CLN|TPCN|SORT|M6PR|IGF2R|LITAF|GGA|AP[1-4])', g):
        return 'membrane, transport and sorting'
    return 'other'
L['family'] = [family(g) for g in L.index]

# ---- 3. cohort monocyte association ----
print('\nloading GSE33615 with the published loader ...', flush=True)
code = open('run_gse33615.py').read().split('# ---------- modules ----------')[0]
cns = {'print': lambda *a, **k: None, 'sys': type('S', (), {'argv': ['x']})}
exec(code, cns)
Gx, samples = cns['G'], cns['samples']
grp = samples.set_index('gsm').loc[Gx.columns, 'group']; is_atl = (grp == 'ATL').values
mono_genes = [g for g in ['CD14','LYZ','CSF1R','FCN1','S100A8','S100A9','VCAN','ITGAM','CD68','MNDA'] if g in Gx.index]
mono = stats.zscore(Gx.loc[mono_genes].values.astype(float), axis=1).mean(axis=0)
rho, dval = {}, {}
for g in L.index:
    if g in Gx.index:
        x = Gx.loc[g].values
        rho[g] = stats.spearmanr(x, mono)[0]
        a, b = x[is_atl], x[~is_atl]
        sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        dval[g] = (a.mean() - b.mean()) / sd if sd > 0 else 0.0
L['rho_monocyte_GSE33615'] = pd.Series(rho)
L['cohens_d_ATL_GSE33615'] = pd.Series(dval)
# sanity check: the loader must reproduce the deposited output of composition_analysis.py exactly
chk = pd.read_csv('results/results_composition_per_gene.csv').drop_duplicates('gene').set_index('gene')['rho_monocyte']
ov = [g for g in chk.index if g in L.index and not np.isnan(L.loc[g, 'rho_monocyte_GSE33615'])]
if ov:
    diff = (L.loc[ov, 'rho_monocyte_GSE33615'] - chk[ov]).abs().max()
    print(f'  loader check against results_composition_per_gene.csv: max |difference| in rho = {diff:.2e} over {len(ov)} genes', flush=True)

L.to_csv('results/results_kegg_lysosome_decomposition.csv')

# ---- report ----
print('\nTop 15 members by gene-level index:')
print(L.head(15)[['gene_CSI', 'percentile', 'nonT_minus_CD4', 'family', 'rho_monocyte_GSE33615',
                  'cohens_d_ATL_GSE33615']].round(3).to_string())
print('\nBy family (median gene-level index, share of members above the 90th genome percentile):')
fam = L.groupby('family').agg(n=('gene_CSI', 'size'), median_CSI=('gene_CSI', 'median'),
                              share_top10pct=('percentile', lambda s: (s > 0.9).mean()))
print(fam.round(3).sort_values('median_CSI', ascending=False).to_string())
pos = L.gene_CSI.clip(lower=0)
top10 = pos.sort_values(ascending=False).head(10)
print(f'\nshare of the summed positive gene-level index carried by the top 10 genes: '
      f'{top10.sum()/pos.sum():.2f}  ({", ".join(top10.index)})')
print(f'members above the 90th genome percentile: {(L.percentile > 0.9).sum()} of {len(L)}; '
      f'above the 50th: {(L.percentile > 0.5).sum()}')
both = L.dropna(subset=['rho_monocyte_GSE33615'])
r, p = stats.spearmanr(both.gene_CSI, both.rho_monocyte_GSE33615)
print(f'agreement, gene-level index vs cohort monocyte association: Spearman rho = {r:+.3f}, '
      f'P = {p:.2g}, n = {len(both)}')
bothx = both.drop([g for g in mono_genes if g in both.index])   # monocyte-score markers are circular here
rx, px = stats.spearmanr(bothx.gene_CSI, bothx.rho_monocyte_GSE33615)
print(f'agreement excluding the monocyte-score markers ({", ".join(g for g in mono_genes if g in both.index)}): '
      f'Spearman rho = {rx:+.3f}, P = {px:.2g}, n = {len(bothx)}   <- value cited in the manuscript', flush=True)
r2, p2 = stats.spearmanr(both.gene_CSI, both.cohens_d_ATL_GSE33615)
print(f'agreement, gene-level index vs apparent ATL effect in the cohort:  Spearman rho = {r2:+.3f}, '
      f'P = {p2:.2g}, n = {len(both)}')
fam.to_csv('results/results_kegg_lysosome_families.csv')
print('\nwrote results_kegg_lysosome_decomposition.csv and results_kegg_lysosome_families.csv')
