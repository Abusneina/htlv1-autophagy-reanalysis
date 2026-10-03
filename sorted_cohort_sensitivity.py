"""Sorted-cohort (GSE55851) check of the enrichment weighting exponent (release 1.4).

Rebuilds the expression matrix from the Agilent files (gse55851_from_raw.py reproduces the
deposited module scores exactly with gProcessedSignal), then rescores the four autophagy modules
with exponents 0, 0.25 and 0.5 and reports the within-donor paired differences (tumor-type
fraction N minus normal-like fraction P), their signs, and the between-donor comparison.

Input: GSE55851_RAW.tar from GEO. Output: results/results_sensitivity_alpha_sorted.csv
"""
import sys, numpy as np, pandas as pd
from scipy import stats
os_makedirs = __import__('os').makedirs; os_makedirs('results', exist_ok=True)
import gse55851_from_raw as g           # reuses the verified loader, annotation and ssGSEA

TAR = sys.argv[1] if len(sys.argv) > 1 else 'GSE55851_RAW.tar'
dep = pd.read_csv('results/results_gse55851_scores.csv', index_col=0)
meta = dep[['gsm', 'patient', 'fraction', 'class']]
G = g.gene_matrix('gProcessedSignal')
print(f'gene matrix {G.shape}', flush=True)

pairs = [p for p in meta.patient.unique() if {'N', 'P'} <= set(meta[meta.patient == p].fraction)]
by = {(r.patient, r.fraction): r.gsm for r in meta.itertuples()}
atlN = meta[(meta['class'] == 'ATL') & (meta.fraction == 'N')].gsm.tolist()
norm = meta[meta['class'] == 'Normal'].gsm.tolist()
print(f'{len(pairs)} paired donors: {pairs}; between-donor groups {len(atlN)} vs {len(norm)}', flush=True)

rows = []
for alpha in (0.0, 0.25, 0.5):
    S = g.ssgsea(G, g.mods, alpha=alpha); S.columns = meta.gsm.values
    for m in ['initiation', 'elongation', 'fusion', 'clear']:
        d = [S.loc[m, by[(p, 'N')]] - S.loc[m, by[(p, 'P')]] for p in pairs]
        a, b = S.loc[m, atlN], S.loc[m, norm]
        rows.append([alpha, m, np.mean(d), int(np.sum(np.array(d) > 0)), len(d), a.mean() - b.mean(),
                     stats.mannwhitneyu(a, b, alternative='two-sided')[1]] + [round(x, 3) for x in d])
    print(f'alpha = {alpha}: scored', flush=True)
R = pd.DataFrame(rows, columns=['alpha', 'module', 'paired_mean_diff', 'donors_up', 'donors', 'between_donor_diff',
                                'between_donor_P'] + [f'diff_{p}' for p in pairs])
R.to_csv('results/results_sensitivity_alpha_sorted.csv', index=False)
pd.set_option('display.width', 250)
print(R.round(3).to_string(index=False))
