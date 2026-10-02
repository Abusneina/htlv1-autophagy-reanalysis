"""Three follow-up checks requested on manuscript v13.

 1. Null width at each module's ACTUAL genomic arrangement, not only at the
    worst-case single block. A 'module-shaped' null is built by sliding the module's
    real gene positions (in genomic gene order) by a common random offset, so the
    spacing between module genes, including which ones share a chromosome and how
    far apart they are, is preserved exactly while their location is randomized.
 2. Whether the ATAC discordance between the empty and IRF4 backgrounds is
    genome-wide or confined to background promoters: correlation of the two ATAC
    contrasts across all promoters, across the top Tax-responsive promoters in
    Jurkat H3K27ac, and across the NF-kB benchmark.
 3. Baseline headroom of the NF-kB benchmark promoters in every contrast.
"""
import os
os.makedirs('results', exist_ok=True)
import re, numpy as np, pandas as pd
from scipy import stats

SEED = 20260817
rng = np.random.default_rng(SEED)

src = open('sharpen_chromatin.py').read()
NFKB = eval(re.search(r'NFKB = (\[.*?\])', src, re.S).group(1))
print('NF-kB benchmark:', len(NFKB), 'genes', flush=True)

reg = pd.read_csv('regions_hg38.tsv', sep='\t').sort_values(['chrom', 'tss']).reset_index(drop=True)
mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(x['gene'])) for m, x in mt.groupby('module') if m != 'axis'}
mods['NF-kB benchmark'] = NFKB

ARMS = {
 'H3K27ac Jurkat': 'results/results_h3k27ac_lfc_CT_Jurkat_Tax_H3K27ac.csv',
 'H3K27ac TL-Om1': 'results/results_h3k27ac_lfc_CT_TL-Om1_Tax_H3K27ac.csv',
 'ATAC empty':     'results/results_atac_lfc_Empty_Tax.csv',
 'ATAC IRF4':      'results/results_atac_lfc_IRF4_Tax.csv',
}
LFC = {a: pd.read_csv(p, index_col=0).iloc[:, 0] for a, p in ARMS.items()}

# ------------------------------------------------ 1. module-shaped null
print('\n' + '=' * 74)
print('1. NULL AT EACH MODULE\'S ACTUAL GENOMIC ARRANGEMENT (20,000 slides)')
print('=' * 74, flush=True)
rows = []
NDRAW = 20000
for arm, lfc in LFC.items():
    d = reg[reg.symbol.isin(lfc.index)].reset_index(drop=True)
    v = lfc.reindex(d.symbol).values
    N = len(v)
    idx_of = {g: i for i, g in enumerate(d.symbol)}
    print(f'\n  {arm}', flush=True)
    for m, genes in mods.items():
        ii = np.array(sorted(idx_of[g] for g in genes if g in idx_of))
        n = len(ii)
        # scattered null
        scat = np.array([v[rng.choice(N, n, replace=False)].mean() for _ in range(NDRAW)])
        # module-shaped null: common offset, wrap around the ordered genome
        offs = rng.integers(0, N, NDRAW)
        shaped = np.array([v[(ii + o) % N].mean() for o in offs])
        infl = shaped.std() / scat.std()
        mde_s = 2.80 * scat.std(); mde_a = 2.80 * shaped.std()
        rows.append([arm, m, n, round(scat.std(), 5), round(shaped.std(), 5), round(infl, 3),
                     round(2 ** mde_s, 3), round(2 ** mde_a, 3)])
        print(f'    {m:16s} n={n:2d}  scattered sd {scat.std():.4f}  actual-arrangement sd '
              f'{shaped.std():.4f}  inflation {infl:.2f}x  bound {2**mde_s:.2f}x -> {2**mde_a:.2f}x',
              flush=True)
R1 = pd.DataFrame(rows, columns=['contrast', 'set', 'n_genes', 'scattered_null_sd',
                                 'arrangement_null_sd', 'inflation', 'bound_fold_scattered',
                                 'bound_fold_arrangement'])
R1.to_csv('results/results_arrangement_null.csv', index=False)

# ------------------------------------------------ 2. ATAC discordance
print('\n' + '=' * 74)
print('2. WHERE DO THE TWO ATAC CONTRASTS DISAGREE?')
print('=' * 74, flush=True)
a, b = LFC['ATAC empty'], LFC['ATAC IRF4']
common = a.index.intersection(b.index)
jur = LFC['H3K27ac Jurkat'].reindex(common).dropna()
top = jur.sort_values(ascending=False).index[:500]
out2 = []
for label, genes in [('all promoters', common),
                     ('top 500 Jurkat H3K27ac Tax gains', top),
                     ('NF-kB benchmark', [g for g in NFKB if g in common])]:
    r, p = stats.spearmanr(a[genes], b[genes])
    out2.append([label, len(genes), round(r, 4), p,
                 round(float(a[genes].mean()), 3), round(float(b[genes].mean()), 3)])
    print(f'  {label:34s} n={len(genes):5d}  rho={r:+.3f}  P={p:.2g}  '
          f'mean empty {a[genes].mean():+.3f}  mean IRF4 {b[genes].mean():+.3f}', flush=True)
# does either ATAC contrast track the Jurkat H3K27ac response?
for lab, s in [('ATAC empty', a), ('ATAC IRF4', b)]:
    r, p = stats.spearmanr(s[jur.index], jur)
    out2.append([f'{lab} vs Jurkat H3K27ac, all promoters', len(jur), round(r, 4), p, None, None])
    print(f'  {lab} vs Jurkat H3K27ac (all)     rho={r:+.3f}  P={p:.2g}', flush=True)
pd.DataFrame(out2, columns=['subset', 'n', 'spearman', 'P', 'mean_empty', 'mean_IRF4']).to_csv(
    'results/results_atac_discordance.csv', index=False)

# ------------------------------------------------ 3. baseline headroom
print('\n' + '=' * 74)
print('3. BASELINE PERCENTILE OF THE NF-kB BENCHMARK PROMOTERS')
print('=' * 74, flush=True)
BASE = {
 'H3K27ac Jurkat': ('results/results_h3k27ac_promoters.csv', 'CT_Jurkat_Empty_H3K27ac'),
 'H3K27ac TL-Om1': ('results/results_h3k27ac_promoters.csv', 'CT_TL-Om1_Empty_H3K27ac'),
 'ATAC empty':     ('results/results_atac_promoters.csv', 'Empty'),
 'ATAC IRF4':      ('results/results_atac_promoters.csv', 'IRF4'),
}
out3 = []
for arm, (path, col_hint) in BASE.items():
    t = pd.read_csv(path, index_col=0)
    cols = [c for c in t.columns if col_hint.lower() in c.lower() and 'tax' not in c.lower()]
    if not cols:
        print(f'  {arm}: baseline column not found among {list(t.columns)[:8]}')
        continue
    base = t[cols[0]]
    pct = base.rank(pct=True)
    b_nf = pct.reindex(NFKB).dropna()
    b_mod = {m: pct.reindex(g).dropna().median() for m, g in mods.items() if m != 'NF-kB benchmark'}
    out3.append([arm, cols[0], round(b_nf.median() * 100, 1)] +
                [round(b_mod[m] * 100, 1) for m in sorted(b_mod)])
    print(f'  {arm:16s} baseline column {cols[0]:28s} NF-kB median percentile '
          f'{b_nf.median()*100:5.1f}   modules ' +
          ', '.join(f'{m} {b_mod[m]*100:.0f}' for m in sorted(b_mod)), flush=True)
pd.DataFrame(out3, columns=['contrast', 'baseline_column', 'nfkb_median_pct'] +
             [f'{m}_median_pct' for m in sorted(b_mod)]).to_csv('results/results_baseline_headroom.csv', index=False)
print('\nDONE')
