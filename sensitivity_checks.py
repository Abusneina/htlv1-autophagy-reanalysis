"""Sensitivity checks requested in review (release 1.4).

A. ssGSEA weighting exponent. The unsorted cohort (GSE33615) is rescored with
   alpha = 0, 0.25 (published) and 0.5, and the ATL-versus-control result for each
   module is reported. Uses the loader and scoring function of run_gse33615.py.
B. TFEB inside the CLEAR module. CLEAR is rescored with and without TFEB in the
   unsorted cohort and in the synthetic specimens.
C. Composition of the non-T pool. The synthetic specimens of insilico_mixing.py are
   rebuilt with alternative non-T pools, and the non-T content at which the CLEAR
   effect reaches the value observed in GSE33615 (Cohen d = 1.20) is reported.

Inputs: the files used by run_gse33615.py and gene_level_csi.py.
Outputs: results/results_sensitivity_alpha.csv, results/results_sensitivity_tfeb.csv,
         results/results_sensitivity_nonT_pool.csv
"""
import os, sys, numpy as np, pandas as pd
from scipy import stats
os.makedirs('results', exist_ok=True)
SEED = 20260817

# ------------------------------------------------------------- A and B, GSE33615
print('A. loading GSE33615 with the published loader ...', flush=True)
src = open('run_gse33615.py').read()
head, tail = src.split('# ---------- modules ----------')
ns = {'print': lambda *a, **k: None}
exec(head, ns)
G, samples = ns['G'], ns['samples']
# the published ssGSEA function, taken verbatim from run_gse33615.py
fn = tail[tail.find('def ssgsea'):tail.find('S = ssgsea(G, mods)')]
exec(fn, ns)
ssgsea = ns['ssgsea']

mt = pd.read_csv('modules_frozen_v1.tsv', sep='\t')
mods = {m: sorted(set(g['gene'])) for m, g in mt.groupby('module') if m != 'axis'}
grp = samples.set_index('gsm').loc[G.columns, 'group']
atl = (grp == 'ATL').values

def contrast(S):
    rows = []
    for m in S.index:
        a, b = S.loc[m][atl], S.loc[m][~atl]
        p = stats.mannwhitneyu(a, b, alternative='two-sided')[1]
        d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        rows.append([m, d, p])
    return pd.DataFrame(rows, columns=['module', 'cohens_d', 'p'])

out = []
for alpha in (0.0, 0.25, 0.5):
    R = contrast(ssgsea(G, mods, alpha=alpha))
    R['q'] = stats.false_discovery_control(R['p'])
    R.insert(0, 'alpha', alpha)
    out.append(R)
    print(f'   alpha = {alpha:<4}  ' + '   '.join(f'{r.module} d={r.cohens_d:+.2f} q={r.q:.2g}' for r in R.itertuples()), flush=True)
A = pd.concat(out, ignore_index=True)
pub = pd.read_csv('results/results_modules.csv').set_index('module')
chk = A[A.alpha == 0.25].set_index('module')
print(f'   check: alpha 0.25 reproduces results_modules.csv, max |diff| in d = '
      f'{(chk.cohens_d - pub.cohens_d.reindex(chk.index)).abs().max():.1e}', flush=True)
A.to_csv('results/results_sensitivity_alpha.csv', index=False)

print('B. CLEAR with and without TFEB in GSE33615 ...', flush=True)
clear = mods['clear']; clear_noT = [g for g in clear if g != 'TFEB']
S = ssgsea(G, {'clear': clear, 'clear_without_TFEB': clear_noT}, alpha=0.25)
B = contrast(S)
rest = ssgsea(G, {'rest': clear_noT}, alpha=0.25).loc['rest']
rho_tfeb = stats.spearmanr(G.loc['TFEB'].values, rest.values)[0] if 'TFEB' in G.index else np.nan
B['cohort'] = 'GSE33615'
print('   ' + '   '.join(f'{r.module} d={r.cohens_d:+.2f} P={r.p:.2g}' for r in B.itertuples())
      + f'   TFEB vs rest of module, Spearman rho = {rho_tfeb:+.2f}', flush=True)

# ------------------------------------------------------------- C (and B), GSE107011 mixtures
print('C. building synthetic specimens from GSE107011 ...', flush=True)
import pyreadr
d = pd.read_csv('GSE107011_Processed_data_TPM.txt.gz', sep='\t', index_col=0)
d.index = [i.split('.')[0] for i in d.index]
ann = pyreadr.read_r('grch38.rda')['grch38'][['ensgene', 'symbol']].dropna().drop_duplicates('ensgene')
g = ann.set_index('ensgene')['symbol'].reindex(d.index)
d = d[g.notna().values]; d.index = g.dropna().values
d = d.groupby(level=0).max()

T = ['CD4_naive', 'Th1', 'Th2', 'Th17', 'Treg']
POPS = ['C_mono', 'I_mono', 'NC_mono', 'B_naive', 'NK', 'Neutrophils']
donors = sorted(dn for dn in {c.split('_', 1)[0] for c in d.columns}
                if all(f'{dn}_{ct}' in d.columns for ct in T + POPS))
FR = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]

# cohort-matched pool: weights proportional to the excess of each population in ATL
# specimens over controls in the deposited LM22 fractions (negative excess set to zero)
Gr = pd.read_csv('results/results_cibersortx_groups.csv').set_index('population')['difference']
mono_ex = max(Gr.get('Monocytes', 0), 0)
b_ex = max(Gr.get('B cells naive', 0), 0)
nk_ex = max(Gr.get('NK cells resting', 0), 0) + max(Gr.get('NK cells activated', 0), 0)
ne_ex = max(Gr.get('Neutrophils', 0), 0)
mono_split = np.array([0.40, 0.08, 0.07]) / 0.55
cm = dict(zip(['C_mono', 'I_mono', 'NC_mono'], mono_ex * mono_split))
cm.update({'B_naive': b_ex, 'NK': nk_ex, 'Neutrophils': ne_ex})
tot = sum(cm.values()); cm = {k: v / tot for k, v in cm.items()}

POOLS = {
    'published (monocytes 55%, B 25%, NK 12%, neutrophils 8%)':
        {'C_mono': 0.40, 'I_mono': 0.08, 'NC_mono': 0.07, 'B_naive': 0.25, 'NK': 0.12, 'Neutrophils': 0.08},
    'monocytes only':
        {'C_mono': 0.40 / 0.55, 'I_mono': 0.08 / 0.55, 'NC_mono': 0.07 / 0.55},
    'equal weights across the six populations':
        {k: 1 / 6 for k in POPS},
    'matched to the ATL-minus-control excess in LM22': cm,
    'no monocytes (B, NK, neutrophils only)':
        {'B_naive': 0.25 / 0.45, 'NK': 0.12 / 0.45, 'Neutrophils': 0.08 / 0.45},
}

def specimens(weights):
    cols = {}
    for dn in donors:
        cd4 = d[[f'{dn}_{ct}' for ct in T]].mean(axis=1)
        oth = sum(w * d[f'{dn}_{ct}'] for ct, w in weights.items())
        for f in FR:
            cols[f'{dn}_f{int(f*100)}'] = (1 - f) * cd4 + f * oth
    X = np.log2(pd.DataFrame(cols) + 1)
    return X[X.var(axis=1) > 0]

def d_by_fraction(scores):
    base = np.array([scores[f'{dn}_f0'] for dn in donors]); ys = []
    for f in FR[1:]:
        cur = np.array([scores[f'{dn}_f{int(f*100)}'] for dn in donors])
        sd = np.sqrt((cur.var(ddof=1) + base.var(ddof=1)) / 2)
        ys.append((cur.mean() - base.mean()) / sd if sd > 0 else 0.0)
    return np.array(ys)

def crossing(ds, target=1.20):
    x = np.r_[0, np.array(FR[1:])] * 100; y = np.r_[0, ds]
    if y.max() < target: return np.nan
    i = np.argmax(y >= target)
    return x[i - 1] + (target - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])

rowsC, rowsB = [], []
for i, (label, w) in enumerate(POOLS.items(), 1):
    print(f'   [{i}/{len(POOLS)}] {label}', flush=True)
    X = specimens(w)
    S = ssgsea(X, {**mods, 'clear_without_TFEB': clear_noT}, alpha=0.25)
    for m in ['clear', 'initiation', 'clear_without_TFEB']:
        ds = d_by_fraction(S.loc[m])
        rowsC.append([label, m, *np.round(ds, 3), round(crossing(ds), 2) if not np.isnan(crossing(ds)) else None])
    print('        CLEAR reaches d = 1.20 at ' +
          (f'{crossing(d_by_fraction(S.loc["clear"])):.1f}% non-T' if not np.isnan(crossing(d_by_fraction(S.loc['clear']))) else 'no level tested') +
          '; initiation maximum d = ' + f'{d_by_fraction(S.loc["initiation"]).max():+.2f}', flush=True)
C = pd.DataFrame(rowsC, columns=['nonT_pool', 'set'] + [f'd_at_{int(f*100)}pct' for f in FR[1:]] + ['pct_nonT_to_reach_d_1.20'])
C.to_csv('results/results_sensitivity_nonT_pool.csv', index=False)

pubrow = C[(C.nonT_pool.str.startswith('published'))]
for m in ['clear', 'clear_without_TFEB']:
    r = pubrow[pubrow.set == m].iloc[0]
    B.loc[len(B)] = [m, np.nan, np.nan, f'GSE107011 mixtures: reaches d = 1.20 at {r["pct_nonT_to_reach_d_1.20"]}% non-T']
B['TFEB_vs_rest_spearman_GSE33615'] = round(rho_tfeb, 3)
B.to_csv('results/results_sensitivity_tfeb.csv', index=False)
print('DONE: results_sensitivity_alpha.csv, results_sensitivity_tfeb.csv, results_sensitivity_nonT_pool.csv')
