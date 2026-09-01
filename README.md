# Composition confounding in the canonical ATL expression cohort, and HTLV-1 Tax at autophagy loci

Code and results for the manuscript by A. M. Abusneina and colleagues
(University of Benghazi). No new experimental data were generated. Every input is
public; this repository contains the analysis code, the frozen pre-analysis
specification, all derived result tables and the manuscript figures.

## What this repository establishes

1. The most reused adult T-cell leukemia expression cohort (GSE33615) compares
   patient peripheral blood mononuclear cells with purified healthy CD4+ T cells.
   Monocytes make up 16.5 percent of one arm and 6.2 percent of the other, and a
   gene's association with monocyte content predicts its apparent disease effect
   (r = 0.841). Every one of the top 100 differentially expressed genes, in any
   pathway, tracks erythroid content.
2. Synthetic specimens built from purified populations at known contamination
   levels reproduce the observed lysosomal effect at 6.6 percent non-T content.
3. A Composition Sensitivity Index, defined from those mixtures, is reported for
   1,635 gene sets across Hallmark, KEGG, KEGG MEDICUS and Reactome.
4. HTLV-1 Tax produces no detectable change in H3K27ac or chromatin accessibility
   at autophagy loci, at a sensitivity that comfortably detects its NF-kB effect.

## Input data (not redistributed here)

Gene Expression Omnibus: GSE33615 and GSE55851 (expression), GSE316417 (RNA-seq),
GSE316418 (ATAC), GSE316419 (CUT&Tag), GSE107011 (sorted immune populations).
Gene set collections are from MSigDB; only the Hallmark collection is
redistributed here. See THIRD_PARTY_NOTICES.md for what must be downloaded
separately and why. CIBERSORTx LM22 fractions were computed on
the CIBERSORTx web portal and the output is included as
`CIBERSORTx_LM22_fractions_GSE33615.csv`.

## Pre-specification

`modules_frozen_v1.tsv` holds the four autophagy modules and the three-gene
acetylation panel. It was frozen before any dataset was opened and was not edited
afterwards. Do not modify it if you wish to reproduce the published numbers.

## Scripts, in execution order

| Script | Arm |
| --- | --- |
| `run_gse33615.py` | Unsorted expression cohort |
| `run_gse55851.py` | FACS-sorted expression cohort |
| `composition_analysis.py`, `check_composition.py` | Marker-based composition tests |
| `run_mcpcounter.py`, `run_cibersortx.py` | Reference-based deconvolution |
| `pathway_survey.py` | Hallmark survey and Figure 3 |
| `test_wang_iha.py` | Published gene lists tested against composition |
| `build_regions.py` | GRCh38 promoter coordinates |
| `run_chromatin_arm.py`, `run_atac_arm.py`, `run_tax_arm.py` | Tax contrasts |
| `sharpen_chromatin.py` | Baseline-matched null, bootstrap, combination |
| `power_analysis.py` | Detection bounds |
| `insilico_mixing.py` | Synthetic specimens and dose response |
| `csi_metric.py` | Composition Sensitivity Index |
| `make_figures.py`, `make_fig6.py` | Figures 1, 2, 5, 6 and Figure 4 |
| `build_ms.js` | Builds the manuscript document |

`run_deconvolution_LOCAL.R` is an R alternative for the deconvolution step, for
users who prefer to run xCell and MCP-counter locally.

## Fixed analytical choices

- ssGSEA weighting exponent 0.25, implemented from the published algorithm; the
  closed-form implementation in `csi_metric.py` was validated against the direct
  calculation to 1.5e-11
- Probe collapse to gene by maximum mean intensity; quantile normalization between
  samples; no cross-platform merging
- Module coverage floor 0.60; Benjamini-Hochberg at a false discovery rate of 0.05
  within each family of tests
- Promoter windows: TSS +/- 2 kb for H3K27ac, +/- 1 kb for ATAC; paired tracks
  quantile matched before differencing
- Competitive tests use 20,000 size-matched random gene sets, drawn to match
  baseline signal decile where stated; seed 20260817 throughout
- Minimum detectable effect: 2.80 x genome SD / sqrt(set size), 80 percent power,
  two-sided alpha 0.05
- Composition Sensitivity Index: change in Cohen d per 10 percentage points of
  non-T content, referred to a size-matched null of 120 random sets per size bin

## Limitations recorded in the code

- Every chromatin and ATAC contrast has one library per condition; detection
  bounds are reported instead of confidence intervals
- Only four donors in GSE107011 carry all populations needed for mixing, so the
  paired test cannot return P below 0.125 and effect sizes are reported instead
- Ensembl identifiers were mapped to symbols with the annotables GRCh38 table

## Environment

Python 3.12 with numpy, pandas, scipy, statsmodels, pyBigWig, pyreadr and
matplotlib; Node.js with the docx package for the manuscript build.

## License

Code and derived result tables released under the MIT License. Gene set
collections and deconvolution signatures remain under the licenses of their
original providers; see THIRD_PARTY_NOTICES.md.
