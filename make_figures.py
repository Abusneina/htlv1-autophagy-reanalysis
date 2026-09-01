import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':0.8,
                     'xtick.major.width':0.8,'ytick.major.width':0.8,'savefig.dpi':300})
NPG = ['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000']

# ---------------- Figure 1: composition confounding in GSE33615 ----------------
P = pd.read_csv('results_cibersortx_per_gene.csv').rename(
        columns={'rho_monocyte_LM22':'rho_monocyte','cohens_d':'cohens_d_ATL'})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
mods = ['initiation','elongation','fusion','clear']
for i, m in enumerate(mods):
    sub = P[P.module == m]
    ax[0].scatter(sub.rho_monocyte, sub.cohens_d_ATL, s=26, color=NPG[i],
                  edgecolor='white', linewidth=0.5, label=m, zorder=3)
r, p = stats.pearsonr(P.rho_monocyte, P.cohens_d_ATL)
b = np.polyfit(P.rho_monocyte, P.cohens_d_ATL, 1)
xx = np.linspace(P.rho_monocyte.min(), P.rho_monocyte.max(), 50)
ax[0].plot(xx, np.polyval(b, xx), color='0.35', lw=1.2, ls='--', zorder=2)
ax[0].axhline(0, color='0.85', lw=0.7, zorder=1); ax[0].axvline(0, color='0.85', lw=0.7, zorder=1)
ax[0].set_xlabel('Correlation of gene with monocyte fraction (Spearman $\\rho$)')
ax[0].set_ylabel("Apparent ATL effect (Cohen's $d$)")
ax[0].text(0.04, 0.95, f'$r$ = {r:+.3f}\n$P$ = {p:.1e}\n$n$ = {len(P)} genes',
           transform=ax[0].transAxes, va='top', fontsize=7.5)
ax[0].legend(frameon=False, fontsize=7, loc='lower right')

R = pd.read_csv('results_cibersortx_groups.csv')
R = R[R.q < 0.05].sort_values('cohens_d')
y = np.arange(len(R))
ax[1].barh(y-0.19, R.mean_ATL*100, height=0.36, color=NPG[0], label='ATL (PBMC)')
ax[1].barh(y+0.19, R.mean_control*100, height=0.36, color=NPG[3], label='Control (sorted CD4+)')
ax[1].set_yticks(y); ax[1].set_yticklabels(
    [n.replace('T cells ','T ').replace('Dendritic cells','DC').replace('Mast cells','Mast')
     for n in R.population], fontsize=7.5)
ax[1].set_xlabel('Estimated fraction of the specimen (%)')
ax[1].legend(frameon=False, fontsize=7, loc='center right', bbox_to_anchor=(1.0, 0.42))
ax[1].set_xlim(0, max(R.mean_ATL.max(), R.mean_control.max())*128)
for a, l in zip(ax, ['a','b']):
    a.text(-0.17, 1.06, l, transform=a.transAxes, fontweight='bold', fontsize=10)
plt.tight_layout(); plt.savefig('fig/Figure1.png', bbox_inches='tight'); plt.close()

# ---------------- Figure 2: the two cohorts disagree ----------------
fig, ax = plt.subplots(1, 3, figsize=(7.5, 3.0))
d33 = pd.read_csv('results_modules.csv')
order = ['initiation','elongation','fusion','clear']
d33 = d33.set_index('module').loc[order]
cols = [NPG[i] for i in range(4)]
ax[0].bar(range(4), d33.cohens_d, color=cols, width=0.62)
ax[0].axhline(0, color='0.4', lw=0.8)
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(order, rotation=35, ha='right')
ax[0].set_ylabel("ATL vs healthy (Cohen's $d$)")
for i, (q, dd) in enumerate(zip(d33.q, d33.cohens_d)):
    if q < 0.05: ax[0].text(i, dd + 0.07, '*', ha='center', fontsize=10)
ylo, yhi = min(d33.cohens_d.min(), 0), d33.cohens_d.max()
ax[0].set_ylim(ylo*1.30 if ylo < 0 else -0.1, yhi*1.62)
ax[0].text(0.03, 0.97, 'GSE33615\nPBMC vs sorted CD4+', transform=ax[0].transAxes,
           ha='left', va='top', fontsize=7.5)

S55 = pd.read_csv('results_gse55851_scores.csv', index_col=0)
atlN = S55[(S55['class']=='ATL') & (S55.fraction=='N')]
norm = S55[S55['class']=='Normal']
diff = [atlN[m].mean() - norm[m].mean() for m in order]
ax[1].bar(range(4), diff, color=cols, width=0.62)
ax[1].axhline(0, color='0.4', lw=0.8)
ax[1].set_xticks(range(4)); ax[1].set_xticklabels(order, rotation=35, ha='right')
ax[1].set_ylabel('Tumor fraction minus normal CD4+\n(module score)')
ax[1].set_ylim(min(diff)*1.30, max(diff)*1.85)
ax[1].text(0.03, 0.97, 'GSE55851\nall FACS-sorted CD4+', transform=ax[1].transAxes,
           ha='left', va='top', fontsize=7.5)

pats = sorted({t.rsplit('-',1)[0] for t in S55.index if t.endswith('-N')} &
              {t.rsplit('-',1)[0] for t in S55.index if t.endswith('-P')})
for i, m in enumerate(order):
    vals = [S55.loc[f'{p}-N', m] - S55.loc[f'{p}-P', m] for p in pats]
    jit = np.linspace(-0.11, 0.11, len(vals))
    ax[2].scatter(i + jit, vals, s=22, color=cols[i], edgecolor='white', linewidth=0.4, zorder=3)
    ax[2].plot([i-0.24, i+0.24], [np.mean(vals)]*2, color='0.25', lw=1.4, zorder=4)
ax[2].axhline(0, color='0.4', lw=0.8)
ax[2].set_xticks(range(4)); ax[2].set_xticklabels(order, rotation=35, ha='right')
ax[2].set_ylabel('Within-patient difference\n(N minus P fraction)')
lo2, hi2 = ax[2].get_ylim(); ax[2].set_ylim(lo2, hi2 + 0.42*(hi2-lo2))
ax[2].text(0.03, 0.97, f'Paired, $n$ = {len(pats)} donors', transform=ax[2].transAxes,
           ha='left', va='top', fontsize=7.5)
for a, l in zip(ax, ['a','b','c']):
    a.text(-0.28, 1.06, l, transform=a.transAxes, fontweight='bold', fontsize=10)
plt.tight_layout(); plt.savefig('fig/Figure2.png', bbox_inches='tight'); plt.close()

# ---------------- Figure 3: chromatin arms with detection bounds ----------------
pw = pd.read_csv('results_power_analysis.csv')
arms = ['H3K27ac Jurkat Tax','H3K27ac TL-Om1 Tax','ATAC Tax (empty bg)','ATAC Tax (IRF4 bg)']
fig, ax = plt.subplots(figsize=(7.0, 3.1))
w = 0.19
for i, m in enumerate(order):
    xs, ys, mdes = [], [], []
    for j, a in enumerate(arms):
        r = pw[(pw.arm == a) & (pw.module == m)].iloc[0]
        xs.append(j + (i - 1.5) * w); ys.append(r.observed); mdes.append(r.MDE_log2_80pct)
    ax.bar(xs, ys, width=w*0.92, color=NPG[i], label=m, zorder=3)
    for x, md in zip(xs, mdes):
        ax.plot([x-w*0.46, x+w*0.46], [md, md], color='0.25', lw=1.0, zorder=4)
        ax.plot([x-w*0.46, x+w*0.46], [-md, -md], color='0.25', lw=1.0, zorder=4)
ax.axhline(0, color='0.4', lw=0.8)
ax.set_xticks(range(4)); ax.set_xticklabels(['H3K27ac\nJurkat','H3K27ac\nTL-Om1','ATAC\nempty','ATAC\nIRF4'])
ax.set_ylabel('Mean promoter change on Tax expression (log$_2$)')
ax.legend(frameon=False, fontsize=7, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.13))
ax.text(0.985, 0.03, 'horizontal lines: minimum detectable effect (log$_2$), 80% power',
        transform=ax.transAxes, ha='right', fontsize=7, color='0.25')
ax.set_ylim(-0.78, 0.80)
plt.tight_layout(); plt.savefig('fig/Figure3.png', bbox_inches='tight'); plt.close()

# ---------------- Figure 4: benchmark, the test detects a real Tax effect ----------------
h = pd.read_csv('results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv', index_col=0).iloc[:,0]
nfkb = ['RELB','NFKB2','NFKB1','RELA','REL','NFKBIA','TNFAIP3','BIRC3','CD40','ICAM1',
        'TRAF1','CCL5','IL2RA','TNFRSF8','CFLAR','BCL2A1','PLEK','LTA']
nfkb = [g for g in nfkb if g in h.index]
mt = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
sets = {'NF-$\\kappa$B response': nfkb}
for m in order:
    sets[m] = [g for g in mt[mt.module == m].gene if g in h.index]
fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
ax[0].hist(h.values, bins=90, color='0.82', edgecolor='none')
for k, (name, gs) in enumerate(sets.items()):
    mu = h[gs].mean()
    ax[0].axvline(mu, color=(NPG[7] if k == 0 else NPG[k-1]), lw=1.6 if k == 0 else 1.1,
                  label=f'{name} ({mu:+.2f})')
ax[0].set_xlim(-2.2, 2.6); ax[0].set_yscale('log')
ax[0].set_xlabel('Promoter H3K27ac change on Tax (log$_2$)')
ax[0].set_ylabel('Number of promoters')
ax[0].legend(fontsize=7, loc='upper right', handlelength=1.4, frameon=True, facecolor='white', edgecolor='none', framealpha=0.92)

names, ratios = [], []
for name, gs in sets.items():
    s = h.std(); mde = 2.80 * s / np.sqrt(len(gs))
    names.append(name); ratios.append(h[gs].mean() / mde)
colors = [NPG[7]] + NPG[0:4]
ax[1].barh(range(len(names)), ratios, color=colors, height=0.6)
ax[1].axvline(1, color='0.25', lw=1.1, ls='--')
ax[1].set_yticks(range(len(names))); ax[1].set_yticklabels(names); ax[1].invert_yaxis()
ax[1].set_xlabel('Observed effect / minimum detectable effect')
ax[1].set_ylim(len(names)-0.4, -0.9)
ax[1].text(1.06, -0.75, 'detection threshold', fontsize=7, color='0.25', va='center')
for a, l in zip(ax, ['a','b']):
    a.text(-0.18, 1.06, l, transform=a.transAxes, fontweight='bold', fontsize=10)
plt.tight_layout(); plt.savefig('fig/Figure4.png', bbox_inches='tight'); plt.close()
print('figures written')
