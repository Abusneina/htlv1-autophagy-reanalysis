# Release v1.5 (3 October 2026)

Adds the sorted-cohort check of the enrichment weighting exponent that was requested in review. Nothing from release 1.4 changes.

## New scripts

| Script | Purpose | Output in `results/` |
|---|---|---|
| `gse55851_from_raw.py` | Rebuilds the GSE55851 expression matrix from the per-sample Agilent files in `GSE55851_RAW.tar` and tests four candidate signal columns against the deposited module scores | prints to the console |
| `sorted_cohort_sensitivity.py` | Rescored the four modules with exponents 0, 0.25 and 0.5 and reports the within-donor and between-donor results | `results_sensitivity_alpha_sorted.csv` |

## Result

With `gProcessedSignal` the rebuilt matrix reproduces `results_gse55851_scores.csv` with a maximum difference of 0.0000 over all 84 values; the other candidates (`gMedianSignal`, `gBGSubSignal`) differ by up to 0.10 and 0.04. Exponent 0.25 reproduces the manuscript values (paired means +0.18, +0.19, +0.03, +0.02; ratios to the detection bound 1.41 and 1.35 for initiation and elongation).

- Initiation rises in all 5 donors at every exponent (paired mean +0.16, +0.18, +0.19; ratio to bound 1.39, 1.41, 1.43).
- Elongation rises in all 5 donors at every exponent (paired mean +0.21, +0.19, +0.17; ratio to bound 1.54, 1.35, 1.17).
- Fusion and CLEAR stay below their bounds at every exponent (ratios 0.32 to 0.35 and 0.40 to below 0.01), with signs mixed, so they remain uninformative.
- Between-donor comparison (6 tumor-type fractions against 3 normal CD4+ samples): elongation (+0.33) and the CLEAR reduction (-0.44 to -0.40) are completely separated (P = 0.024) at every exponent; the initiation difference is +0.17, +0.21 and +0.26 and is completely separated only at exponent 0.5 (P = 0.167, 0.048, 0.024).

## Changed files

`README.md` (two script rows, one input row, one sentence on the rebuild) and `.gitignore` (adds `GSE55851_RAW.tar`).
