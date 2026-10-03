"""Graphical abstract (Virology). Every plotted number is read from a deposited table.

Inputs (all in results/): results_cibersortx_groups.csv, results_cibersortx_per_gene.csv,
results_power_analysis.csv, results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv; the NF-kB
benchmark gene list is read from sharpen_chromatin.py.
Outputs: fig/Graphical_Abstract.tif, fig/Graphical_Abstract.png and
results/results_graphical_abstract_data.csv, which lists every plotted value.
"""
import os, re
os.makedirs('fig', exist_ok=True); os.makedirs('results', exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Virology asks for 531 x 1328 px (h x w) or proportionally more; readable at 5 x 13 cm.
# Built at 2x that ratio, 300 dpi.
W_CM, H_CM = 13.0, 5.2
fig = plt.figure(figsize=(W_CM / 2.54, H_CM / 2.54), dpi=300)

# palette: maroon / teal / amber / slate, no blue theme
MAROON, TEAL, AMBER, SLATE, GREY = '#7A1F2B', '#1B6F6A', '#C8862B', '#4A4A4A', '#BFBFBF'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.2})

# ---------------- plotted data, read from the deposit ----------------
Gr = pd.read_csv('results/results_cibersortx_groups.csv').set_index('population')
mono_ctrl, mono_atl = Gr.loc['Monocytes', ['mean_control', 'mean_ATL']] * 100
cd8_ctrl, cd8_atl = Gr.loc['T cells CD8', ['mean_control', 'mean_ATL']] * 100
PW = pd.read_csv('results/results_power_analysis.csv')
PWJ = PW[PW.arm.str.contains('Jurkat')].set_index('module')
ratio = (PWJ.observed.abs() / PWJ.MDE_log2_80pct)
NFKB = eval(re.search(r'NFKB = (\[.*?\])', open('sharpen_chromatin.py').read(), re.S).group(1))
lfc = pd.read_csv('results/results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv', index_col=0).iloc[:, 0]
bench = lfc.reindex(NFKB).dropna()
bench_ratio = bench.mean() / (2.80 * lfc.std() / np.sqrt(len(bench)))
assert abs(bench_ratio - 2.47) < 0.005, bench_ratio
prov = [('panel 1', 'Monocytes, sorted CD4+ controls (%)', mono_ctrl), ('panel 1', 'Monocytes, ATL blood (%)', mono_atl),
        ('panel 1', 'CD8 T cells, sorted CD4+ controls (%)', cd8_ctrl), ('panel 1', 'CD8 T cells, ATL blood (%)', cd8_atl)]

# explicit positions (fractions of the figure): panel 3 keeps a label zone of its own, 0.63 to 0.78
POS = [[0.090, 0.385, 0.195, 0.535], [0.440, 0.385, 0.175, 0.535], [0.790, 0.385, 0.190, 0.535]]

# ---- panel 1: composition ----
ax = fig.add_axes(POS[0]); axes_all = [ax]
x = np.arange(2)
ax.bar(x - 0.19, [mono_ctrl, cd8_ctrl], width=0.36, color=GREY, edgecolor='none')
ax.bar(x + 0.19, [mono_atl, cd8_atl], width=0.36, color=MAROON, edgecolor='none')
ax.set_xticks(x)
ax.set_xticklabels(['Monocytes', 'CD8 T cells'])
ax.set_ylabel('Percent of specimen')
ax.set_ylim(0, 21)
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(length=2)
ax.text(0.02, 0.96, 'Sorted CD4+ controls', color=GREY, fontsize=5.6,
        transform=ax.transAxes, va='top', fontweight='bold')
ax.text(0.02, 0.86, 'ATL patient blood', color=MAROON, fontsize=5.6,
        transform=ax.transAxes, va='top', fontweight='bold')
ax.text(0.5, -0.24, 'Cases and controls are\ndifferent cell mixtures',
        transform=ax.transAxes, ha='center', va='top', fontsize=5.6, color=SLATE)

# ---- panel 2: per-gene correlation ----
ax = fig.add_axes(POS[1]); axes_all.append(ax)
P = pd.read_csv('results/results_cibersortx_per_gene.csv')
xcol = [c for c in P.columns if 'rho' in c][0]; ycol = [c for c in P.columns if 'cohen' in c][0]
mono, eff = P[xcol].values, P[ycol].values
r_real = np.corrcoef(mono, eff)[0, 1]
assert abs(r_real - 0.841) < 0.0005, r_real
ax.scatter(mono, eff, s=9, color=TEAL, alpha=0.85, edgecolors='white', linewidths=0.3)
b = np.polyfit(mono, eff, 1)
xs = np.linspace(mono.min(), mono.max(), 10)
ax.plot(xs, np.polyval(b, xs), color=MAROON, lw=1.1)
ax.axhline(0, color=GREY, lw=0.5); ax.axvline(0, color=GREY, lw=0.5)
prov += [('panel 2', f'{g} (Spearman rho with LM22 monocyte fraction; Cohen d ATL vs control)', f'{a:.6f}; {c:.6f}') for g, a, c in zip(P.gene, mono, eff)]
prov += [('panel 2', 'Pearson r across the points', r_real)]
ax.set_xlabel('Correlation with monocyte fraction')
ax.set_ylabel('Apparent ATL effect ($d$)')
ax.tick_params(length=2)
ax.spines[['top', 'right']].set_visible(False)
ax.text(0.04, 0.95, 'r = 0.841', transform=ax.transAxes, va='top',
        fontsize=7.0, color=MAROON, fontweight='bold')
ax.text(0.5, -0.36, 'Monocyte content predicts the\napparent effect of each of\n51 autophagy genes',
        transform=ax.transAxes, ha='center', va='top', fontsize=5.6, color=SLATE)

# ---- panel 3: chromatin, effect over detection bound ----
ax = fig.add_axes(POS[2]); axes_all.append(ax)
labels = ['NF-\u03baB control', 'CLEAR', 'initiation', 'fusion', 'elongation']
vals = [bench_ratio, ratio['clear'], ratio['initiation'], ratio['fusion'], ratio['elongation']]
prov += [('panel 3', f'{l}: observed effect / detection bound, Jurkat H3K27ac', v) for l, v in zip(labels, vals)]
cols = [AMBER, TEAL, TEAL, TEAL, TEAL]
y = np.arange(len(labels))[::-1]
ax.barh(y, vals, height=0.62, color=cols, edgecolor='none')
ax.axvline(1.0, color=SLATE, ls='--', lw=0.9)
ax.set_yticks(y); ax.set_yticklabels(labels)
ax.set_xlabel('Observed effect / detection bound', ha='right', x=1.0)
ax.set_xlim(0, 2.9)
ax.set_ylim(-0.6, 4.9)
ax.spines[['top', 'right']].set_visible(False)
ax.tick_params(length=2)
ax.text(1.07, 0.55, 'detection\nthreshold', fontsize=5.2, color=SLATE, va='center')
ax.text(0.5, -0.36, 'Tax moves NF-\u03baB loci;\nautophagy loci stay below\nthe bound (Jurkat)',
        transform=ax.transAxes, ha='center', va='top', fontsize=5.2, color=SLATE)


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

check_layout(fig, axes_all, ['a', 'b', 'c'])
fig.canvas.draw()
for ax_, ch in zip(axes_all, 'abc'):
    tb = ax_.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(max(tb.x0, 0.003), ax_.get_position().y1 + 0.025, ch, fontweight='bold', fontsize=8, va='bottom', ha='left')

out = 'fig/Graphical_Abstract.tif'
fig.savefig(out, dpi=300, format='tiff', pil_kwargs={'compression': 'tiff_lzw'},
            facecolor='white')
fig.savefig('fig/Graphical_Abstract.png', dpi=300, facecolor='white')
pd.DataFrame(prov, columns=['panel', 'quantity', 'value']).to_csv('results/results_graphical_abstract_data.csv', index=False)
plt.close(fig)

from PIL import Image
for f in [out, 'fig/Graphical_Abstract.png']:
    im = Image.open(f)
    print(f.split('/')[-1], im.size, '(w x h)', 'ratio', round(im.size[0] / im.size[1], 3),
          '| required ratio', round(1328 / 531, 3))
