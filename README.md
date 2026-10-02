# Composition confounding in the canonical ATL expression cohort, and HTLV-1 Tax at autophagy loci

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22235481.svg)](https://doi.org/10.5281/zenodo.22235481)

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
4. In the one chromatin contrast whose NF-kB positive control clears its own detection
   bound (Jurkat H3K27ac), HTLV-1 Tax produces no detectable change in H3K27ac at
   autophagy loci. That excludes an effect of the size Tax produces at NF-kB loci and
   does not exclude one below roughly 1.5-fold relative to the rest of the genome.
   Three further contrasts have flat controls and are reported as uninformative.

## Input data (not redistributed here, with one exception)

The two GSE33615 files (`GSE33615_series_matrix.txt.gz`, `GSE33615_family.soft.gz`) are
included because the erythroid and gene-level scripts read them from this folder.
Everything else is downloaded by the user.

Gene Expression Omnibus: GSE33615 and GSE55851 (expression), GSE316417 (RNA-seq),
GSE316418 (ATAC), GSE316419 (CUT&Tag), GSE107011 (sorted immune populations).
Gene set collections are from MSigDB; only the Hallmark collection is
redistributed here. See THIRD_PARTY_NOTICES.md for what must be downloaded
separately and why. CIBERSORTx LM22 fractions were computed on
the CIBERSORTx web portal and the output is included as
`CIBERSORTx_LM22_fractions_GSE33615.csv`.

### Inputs to place in this folder before running

Every script reads its inputs from this folder and writes its tables to `results/`. None of the files below is committed (see `.gitignore`).

| File or folder | Source | Used by |
| --- | --- | --- |
| `GSE55851_series_matrix.txt.gz`, `GSE55851_family.soft.gz` | GEO, series GSE55851 | `run_gse55851.py` |
| `GSE316417_RNA-seq_raw_counts_all_samples.txt.gz` | GEO, series GSE316417, supplementary file | `run_tax_arm.py` |
| `bw2/` with the H3K27ac CUT&Tag bigWig files | GEO, series GSE316419 | `run_chromatin_arm.py`, `sharpen_chromatin.py` |
| `atac/` with the ATAC bigWig files | GEO, series GSE316418 | `run_atac_arm.py` |
| `GSE107011_Processed_data_TPM.txt.gz` | GEO, series GSE107011, supplementary file | `insilico_mixing.py`, `csi_metric.py`, `gene_level_csi.py`, `sensitivity_checks.py` |
| `grch38.rda` | https://github.com/stephenturner/annotables (file `data/grch38.rda`) | `build_regions.py`, the GSE107011 scripts, `run_tax_arm.py` |
| MSigDB symbols files for KEGG legacy, KEGG MEDICUS and Reactome | MSigDB (free registration; licensed, so not redistributed) | `csi_metric.py` reads every `.gmt` file in the folder; `kegg_lysosome_decomposition.py` reads the KEGG legacy file |

The two GSE33615 files and the Hallmark collection are committed. `test_wang_iha.py` reads `results/wang_iha_degs.json`, the gene lists published by Wang and Iha (2023), which is committed.

### Reproducibility check (release 1.4)

The expression arm was re-run from a fresh copy of this repository: `run_gse33615.py`, `run_mcpcounter.py`, `run_cibersortx.py`, `composition_analysis.py`, `check_composition.py`, `pathway_survey.py`, `test_wang_iha.py`, `power_analysis.py`, `insilico_mixing.py`, `make_figures.py` and `make_fig6.py`. Every table that had been deposited before was reproduced with a maximum difference of zero, except `results_insilico_mixing.csv`: its four autophagy modules agree to within 0.01 in Cohen d and the CLEAR crossing is 6.58 percent against 6.62 percent, but its Hallmark rows, which the manuscript does not cite, differ, and the regenerated table replaces the earlier one. Six tables that the figure and adjustment scripts write had never been deposited and are now included, among them `results_cibersortx_per_gene.csv`, the source of Figure 1a (r = 0.841). The chromatin scripts were not re-run here, since they need the bigWig downloads, and `csi_metric.py` needs the licensed MSigDB collections.

## Pre-specification

`modules_frozen_v1.tsv` holds the four autophagy modules and the three-gene
acetylation panel. It was frozen before any dataset was opened and was not edited
afterwards. Do not modify it if you wish to reproduce the published numbers.
The pre-analysis protocol, 2026-08-17_1800_HTLV1-Autophagy_Competing-Models_Protocol.docx, holds the eleven discriminating observables and the decision rule.
It was written on 17 August 2026, before any dataset was opened, and is reproduced here unedited. 

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
| `spatial_null.py` | Spatial dependence of the detection bound: autocorrelation, module clustering, block-permutation null, normal-approximation check |
| `followup_checks.py` | Null at each module's actual genomic arrangement, ATAC discordance, baseline headroom of the NF-kB benchmark |
| `erythroid_sensitivity.py` | Erythroid-tracking result rebuilt from marker subsets |
| `gene_level_csi.py` | Composition sensitivity index computed gene by gene |
| `kegg_lysosome_decomposition.py` | KEGG_LYSOSOME index decomposed by gene and family |
| `sensitivity_checks.py` | Enrichment exponent, TFEB inside CLEAR, and composition of the non-T pool |
| `make_graphical_abstract.py` | Graphical abstract; reads every plotted number from the tables and writes them to `results/results_graphical_abstract_data.csv` |
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
- Only the Jurkat H3K27ac contrast carries an NF-kB benchmark that clears its own bound; the
  TL-Om1 and both ATAC contrasts are uninformative, and the two ATAC contrasts do not
  correlate with each other (rho = -0.003)
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
