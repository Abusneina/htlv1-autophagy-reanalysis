"""Manuscript Figures 1 and 4, rebuilt from the deposited tables (release 1.6).

Figure 1 (five panels): composition confounding in GSE33615 and the sorted cohort.
Figure 4 (three panels): the chromatin arm with detection bounds and the NF-kB benchmark.

Spacing is solved by matplotlib's constrained layout, so long axis labels cannot run into a
neighboring panel; panel letters are placed from each panel's measured outer edge.
Inputs (results/): results_cibersortx_per_gene.csv, results_cibersortx_groups.csv,
results_modules.csv, results_gse55851_scores.csv, results_power_analysis.csv,
results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv; NF-kB benchmark list from sharpen_chromatin.py.
Outputs: fig/Figure_1.png and fig/Figure_4.png at 300 dpi.
"""
import os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
os.makedirs('fig', exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'savefig.dpi': 300})
NPG = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000']
ORDER = ['initiation', 'elongation', 'fusion', 'clear']
COLS = NPG[:4]

def axin(fig, W, H, left, bottom, width, height, **kw):
    """Axes placed in inches from the lower-left corner of a W x H inch figure."""
    return fig.add_axes([left / W, bottom / H, width / W, height / H], **kw)

def check_layout(fig, axes, names):
    """Fail if any axis label, tick label or in-panel text of one panel enters another panel's plot area or
    another panel's labels, or runs off the canvas."""
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    Wpx, Hpx = fig.get_size_inches() * fig.dpi
    items = {}
    for a, n in zip(axes, names):
        it = [a.xaxis.get_tightbbox(rend), a.yaxis.get_tightbbox(rend)] + [t.get_window_extent(rend) for t in a.texts if t.get_text().strip()]
        lg = a.get_legend()
        if lg is not None: it.append(lg.get_window_extent(rend))
        items[n] = [b for b in it if b is not None and b.width > 0]
    areas = {n: a.get_window_extent(rend) for a, n in zip(axes, names)}
    bad = []
    for n, bl in items.items():
        for b in bl:
            if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > Wpx + 0.5 or b.y1 > Hpx + 0.5: bad.append(f'{n}: text runs off the canvas')
            for m in names:
                if m == n: continue
                if b.overlaps(areas[m]): bad.append(f'{n}: labels enter the plot area of {m}')
                if any(b.overlaps(c) for c in items[m]): bad.append(f'{n}: labels touch the labels of {m}')
    assert not bad, '; '.join(sorted(set(bad)))
    print('  layout check passed (no label enters another panel or leaves the canvas):', ', '.join(names))

SHORT = {'T cells CD8': 'T CD8', 'T cells regulatory (Tregs)': 'Tregs', 'T cells CD4 memory resting': 'T CD4 mem rest',
         'T cells CD4 naive': 'T CD4 naive', 'Dendritic cells activated': 'DC activated', 'Mast cells activated': 'Mast activated'}

def letter(fig, ax, ch):
    """Bold panel letter at the top-left of the panel's outer edge (tick labels included)."""
    fig.canvas.draw()
    tb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    pos = ax.get_position()
    fig.text(max(tb.x0, 0.004), min(pos.y1 + 0.01, 0.99), ch, fontweight='bold', fontsize=10, va='bottom', ha='left')

# =============================================================== Figure 1
P = pd.read_csv('results/results_cibersortx_per_gene.csv')
xcol = [c for c in P.columns if 'rho' in c][0]; ycol = [c for c in P.columns if 'cohen' in c][0]
G = pd.read_csv('results/results_cibersortx_groups.csv')
M = pd.read_csv('results/results_modules.csv').set_index('module').loc[ORDER]
S = pd.read_csv('results/results_gse55851_scores.csv', index_col=0)

W1, H1 = 7.4, 6.1
fig = plt.figure(figsize=(W1, H1))
axA = axin(fig, W1, H1, 0.72, 3.62, 2.58, 2.10)      # label zone of panel b: 3.40 to 4.40 in
axB = axin(fig, W1, H1, 4.40, 3.62, 2.85, 2.10)
axC = axin(fig, W1, H1, 0.72, 1.18, 1.78, 1.78)
axD = axin(fig, W1, H1, 3.22, 1.18, 1.78, 1.78)
axE = axin(fig, W1, H1, 5.82, 1.18, 1.43, 1.78)

for i, m in enumerate(ORDER):
    s = P[P.module == m]
    axA.scatter(s[xcol], s[ycol], s=22, color=COLS[i], edgecolor='white', linewidth=0.5, label=m, zorder=3)
r, p = stats.pearsonr(P[xcol], P[ycol]); b = np.polyfit(P[xcol], P[ycol], 1)
xx = np.linspace(P[xcol].min(), P[xcol].max(), 50)
axA.plot(xx, np.polyval(b, xx), color='0.35', lw=1.2, ls='--', zorder=2)
axA.axhline(0, color='0.85', lw=0.7, zorder=1); axA.axvline(0, color='0.85', lw=0.7, zorder=1)
axA.set_xlabel(r'Gene correlation with monocyte fraction ($\rho$)'); axA.set_ylabel('Apparent ATL effect ($d$)')
axA.text(0.04, 0.96, f'$r$ = {r:+.3f}\n$P$ = {p:.0e}', transform=axA.transAxes, va='top', fontsize=7.5)
axA.legend(frameon=False, fontsize=7, loc='lower right', handletextpad=0.2, borderaxespad=0.3)

R = G[G.q < 0.05].sort_values('cohens_d'); y = np.arange(len(R))
axB.barh(y - 0.19, R.mean_ATL * 100, height=0.36, color=NPG[0], label='ATL (PBMC)')
axB.barh(y + 0.19, R.mean_control * 100, height=0.36, color=NPG[3], label='Control (sorted CD4+)')
axB.set_yticks(y); axB.set_yticklabels([SHORT.get(n, n) for n in R.population], fontsize=7)
axB.set_xlabel('Estimated fraction of specimen (%)')
axB.set_xlim(0, max(R.mean_ATL.max(), R.mean_control.max()) * 135)
axB.legend(frameon=False, fontsize=7, loc='upper right', borderaxespad=0.3)

for ax_, vals, txt in [(axC, M.cohens_d.values, 'GSE33615'),
                       (axD, [S[(S['class'] == 'ATL') & (S.fraction == 'N')][m].mean() - S[S['class'] == 'Normal'][m].mean() for m in ORDER], 'GSE55851')]:
    ax_.bar(range(4), vals, color=COLS, width=0.62); ax_.axhline(0, color='0.4', lw=0.8)
    ax_.set_xticks(range(4)); ax_.set_xticklabels(ORDER, rotation=35, ha='right')
    ax_.text(0.04, 0.97, txt, transform=ax_.transAxes, va='top', fontsize=7.5)
axC.set_ylabel('ATL vs healthy ($d$)'); axD.set_ylabel(r'Tumor $-$ normal')
for i, (q, dd) in enumerate(zip(M.q, M.cohens_d)):
    if q < 0.05: axC.text(i, dd + 0.07, '*', ha='center', fontsize=11)
axC.set_ylim(-0.85, 2.05); axD.set_ylim(-0.5, 0.66)

pats = sorted({t.rsplit('-', 1)[0] for t in S.index if t.endswith('-N')} & {t.rsplit('-', 1)[0] for t in S.index if t.endswith('-P')})
for i, m in enumerate(ORDER):
    vals = [S.loc[f'{q}-N', m] - S.loc[f'{q}-P', m] for q in pats]
    axE.scatter(i + np.linspace(-0.11, 0.11, len(vals)), vals, s=14, color=COLS[i], edgecolor='white', linewidth=0.4, zorder=3)
    axE.plot([i - 0.24, i + 0.24], [np.mean(vals)] * 2, color='0.25', lw=1.4, zorder=4)
axE.axhline(0, color='0.4', lw=0.8); axE.set_xticks(range(4)); axE.set_xticklabels(ORDER, rotation=35, ha='right')
axE.set_ylabel('Within-donor difference'); axE.set_ylim(-0.12, 0.56)
axE.text(0.04, 0.97, f'paired, $n$ = {len(pats)}', transform=axE.transAxes, va='top', fontsize=7.5)
print('Figure 1'); check_layout(fig, [axA, axB, axC, axD, axE], list('abcde'))
for ax_, ch in zip([axA, axB, axC, axD, axE], 'abcde'): letter(fig, ax_, ch)
fig.savefig('fig/Figure_1.png', dpi=300, bbox_inches='tight', pad_inches=0.08); plt.close(fig)

# =============================================================== Figure 4
PW = pd.read_csv('results/results_power_analysis.csv')
arms = [('H3K27ac Jurkat Tax', 'H3K27ac\nJurkat'), ('H3K27ac TL-Om1 Tax', 'H3K27ac\nTL-Om1'),
        ('ATAC Tax (empty bg)', 'ATAC\nempty'), ('ATAC Tax (IRF4 bg)', 'ATAC\nIRF4')]
h = pd.read_csv('results/results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv', index_col=0).iloc[:, 0]
NFKB = eval(re.search(r'NFKB = (\[.*?\])', open('sharpen_chromatin.py').read(), re.S).group(1))
mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
sets = {r'NF-$\kappa$B response': [g for g in NFKB if g in h.index]}
for m in ORDER: sets[m] = [g for g in mt[mt.module == m].gene if g in h.index]

W4, H4 = 7.2, 5.7
fig = plt.figure(figsize=(W4, H4))
ax1 = axin(fig, W4, H4, 0.78, 3.40, 6.35, 2.10)
ax2 = axin(fig, W4, H4, 0.88, 0.72, 2.70, 2.00)
ax3 = axin(fig, W4, H4, 5.05, 0.72, 2.10, 2.00)       # label zone of panel c: 3.75 to 5.05 in
w = 0.19
for i, m in enumerate(ORDER):
    xs, ys, md = [], [], []
    for j, (a, _) in enumerate(arms):
        rr = PW[(PW.arm == a) & (PW.module == m)].iloc[0]
        xs.append(j + (i - 1.5) * w); ys.append(rr.observed); md.append(rr.MDE_log2_80pct)
    ax1.bar(xs, ys, width=w * 0.92, color=COLS[i], label=m, zorder=3)
    for x, v in zip(xs, md):
        ax1.plot([x - w * 0.46, x + w * 0.46], [v, v], color='0.25', lw=1.0, zorder=4)
        ax1.plot([x - w * 0.46, x + w * 0.46], [-v, -v], color='0.25', lw=1.0, zorder=4)
ax1.axhline(0, color='0.4', lw=0.8)
ax1.set_xticks(range(4)); ax1.set_xticklabels([l for _, l in arms])
ax1.set_ylabel(r'Mean promoter change (log$_2$)'); ax1.set_ylim(-0.80, 0.98)
ax1.legend(frameon=False, fontsize=7, ncol=4, loc='upper center', borderaxespad=0.2, columnspacing=1.6)
ax1.text(0.99, 0.04, 'horizontal lines: minimum detectable effect, 80% power', transform=ax1.transAxes,
         ha='right', fontsize=7, color='0.25')

ax2.hist(h.values, bins=90, color='0.82', edgecolor='none')
for k, (name, gs_) in enumerate(sets.items()):
    mu = h[gs_].mean()
    ax2.axvline(mu, color=(NPG[7] if k == 0 else NPG[k - 1]), lw=1.6 if k == 0 else 1.1, label=f'{name} ({mu:+.2f})')
ax2.set_xlim(-2.6, 2.6); ax2.set_yscale('log'); ax2.set_ylim(0.7, 3e5)
ax2.set_xlabel(r'Promoter H3K27ac change on Tax (log$_2$)'); ax2.set_ylabel('Number of promoters')
ax2.legend(fontsize=6.5, loc='upper left', handlelength=1.3, frameon=False, borderaxespad=0.3, labelspacing=0.25)

names, ratios = [], []
for name, gs_ in sets.items():
    names.append(name); ratios.append(h[gs_].mean() / (2.80 * h.std() / np.sqrt(len(gs_))))
ax3.barh(range(len(names)), ratios, color=[NPG[7]] + NPG[0:4], height=0.6)
ax3.axvline(1, color='0.25', lw=1.1, ls='--')
ax3.set_yticks(range(len(names))); ax3.set_yticklabels(names, fontsize=7.5); ax3.invert_yaxis()
ax3.set_xlabel('Observed effect / minimum detectable effect', ha='right', x=1.0); ax3.set_xlim(0, 2.75)
ax3.text(1.05, len(names) - 1, 'detection\nthreshold', fontsize=7, color='0.25', va='center')
print('Figure 4'); check_layout(fig, [ax1, ax2, ax3], list('abc'))
for ax_, ch in zip([ax1, ax2, ax3], 'abc'): letter(fig, ax_, ch)
fig.savefig('fig/Figure_4.png', dpi=300, bbox_inches='tight', pad_inches=0.08); plt.close(fig)
print('fig/Figure_1.png and fig/Figure_4.png written')
