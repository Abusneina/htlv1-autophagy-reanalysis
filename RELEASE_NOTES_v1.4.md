# Release v1.4 (2 October 2026)

Adds the analyses cited in the Virology submission (Sections 2.2, 2.6, 2.8, 3.5 to 3.8). Earlier releases are unchanged.
Every new script uses the project seed 20260817, prints its progress, takes its external inputs from the repository folder, and writes to `results/`.
All six scripts were run end to end from a fresh clone of the repository and reproduced the manuscript values exactly.

## New scripts

| Script | Purpose | Outputs in `results/` |
|---|---|---|
| `spatial_null.py` | Autocorrelation by distance; module clustering; block-permutation null; normal-approximation check at every module size | `results_spatial_autocorrelation.csv`, `results_module_clustering.csv`, `results_block_permutation.csv`, `results_normal_approximation_check.csv` |
| `followup_checks.py` | Null at each module's actual genomic arrangement; ATAC discordance; baseline headroom | `results_arrangement_null.csv`, `results_atac_discordance.csv`, `results_baseline_headroom.csv` |
| `erythroid_sensitivity.py` | Erythroid result from marker subsets | `results_erythroid_panel_sensitivity.csv`, `results_erythroid_marker_correlation.csv` |
| `gene_level_csi.py` | Composition sensitivity index per gene (32,398 genes); CLEAR decomposition | `results_gene_level_CSI.csv`, `results_clear_gene_level_CSI.csv` |
| `kegg_lysosome_decomposition.py` | KEGG_LYSOSOME index by gene and family (needs the user's MSigDB KEGG file) | `results_kegg_lysosome_decomposition.csv`, `results_kegg_lysosome_families.csv` |

## Changed files

- `csi_metric.py`: docstring corrected from 300 to 120 random gene sets per size-matched null. The code and every published value are unchanged.
- `build_regions.py`: paths made relative (`modules_frozen_v1.tsv`, `grch38.rda`).
- `README.md`: new script rows and input instructions; the headline statement about Tax and chromatin now carries the same scope as the manuscript (one informative contrast, three uninformative); the input-data note now mentions the two committed GSE33615 files.
- `CITATION.cff`: the DOI now points to the concept record (10.5281/zenodo.22235481, which always resolves to the newest version) instead of the version 1.0 record, and the version field is set.
- `THIRD_PARTY_NOTICES.md`: lists the three inputs added in this release.
- `.gitignore`: added, to keep downloaded inputs out of the repository.

## New data file

`regions_hg38.tsv`: 19,439 protein-coding loci with hg38 coordinates from the annotables `grch38` table.

## Removed file

`results/results_clear_gene_decomposition.csv` is deleted. No script in the repository generated it, and its monocyte correlations differ slightly from those produced by `composition_analysis.py` (CTSD 0.651 against 0.634, MCOLN1 0.432 against 0.396). The manuscript cites `results_composition_per_gene.csv`.

## Headline results

- Promoter-change correlation is at most r = 0.16 within 100 kb and about 0.01 beyond 500 kb; one same-module gene pair lies within 1 Mb across all four modules.
- At each module's actual arrangement the null widens by 0.99 to 1.03-fold; at the single-block worst case by up to 1.32-fold. The NF-kB benchmark clears its worst-case bound 1.87-fold (2.47-fold scattered).
- Normal approximation: ratio 0.991 to 1.006 across all twelve size-by-contrast combinations.
- The IRF4-background ATAC contrast correlates with nothing (rho = -0.001 with Jurkat H3K27ac); the empty-background contrast carries a weak Tax signal (rho = +0.071).
- Erythroid tracking of the top 100 genes stays at 100 percent after dropping any one marker; four of the five markers intercorrelate at 0.93 to 1.00.
- Gene-level index: CTSD 3.38, TFEB 2.36, CTSB 2.10 within CLEAR (rho = 0.85 with the cohort monocyte association). KEGG_LYSOSOME 4.97 reproduced exactly; 54 of 120 members above the 90th gene-level percentile (rho = 0.69 with the cohort, CD68 excluded).

## Original scripts made portable

Eighteen scripts from earlier releases read inputs from the environment in which they were written and could not run elsewhere: `build_regions.py`, `check_composition.py`, `composition_analysis.py`, `csi_metric.py`, `insilico_mixing.py`, `make_fig6.py`, `make_figures.py`, `pathway_survey.py`, `power_analysis.py`, `run_atac_arm.py`, `run_chromatin_arm.py`, `run_cibersortx.py`, `run_gse33615.py`, `run_gse55851.py`, `run_mcpcounter.py`, `run_tax_arm.py`, `sharpen_chromatin.py` and `test_wang_iha.py`. Each now reads from this folder and writes its tables to `results/`. No analysis changed. `pathway_survey.py` also creates its figure folder before saving. The README lists every input to download.

## Reproducibility check

The expression arm was re-run from a fresh copy of the repository. Every previously deposited table it writes was reproduced with a maximum difference of zero, with one exception described below. The chromatin scripts were not re-run (they need the bigWig downloads), and `csi_metric.py` needs the licensed MSigDB collections; its module values were reproduced to four decimal places by `gene_level_csi.py`.

## Tables added that had never been deposited

`results_cibersortx_per_gene.csv` (source of Figure 1a; r = 0.841, P = 1.2e-14, n = 51), `results_cibersortx_groups.csv` (group fractions: monocytes 16.5 against 6.2 percent, CD8 T cells 8.3 against 1.0 percent), `results_module_adjusted_cibersortx.csv` and `results_module_adjusted_mcp.csv` (Table 1), `results_mcpcounter_estimates.csv` and `results_mcpcounter_groups.csv`.

## Table replaced

`results_insilico_mixing.csv` was regenerated. Its four autophagy modules agree with the earlier version to within 0.01 in Cohen d, and the CLEAR crossing of the observed effect is 6.58 percent against 6.62 percent, so the published 6.6 percent stands. Its Hallmark rows differ, by up to 6.4 in Cohen d for some sets; the cause was not established, the manuscript does not cite them, and the regenerated table is the one that the deposited code produces.

## Sensitivity checks (`sensitivity_checks.py`)

- Enrichment exponent: with alpha = 0 and 0.5 instead of 0.25, CLEAR gives d = 1.15 and 1.23 (q below 0.001) and initiation d = -0.62 and -0.66; elongation, never significant, changes sign. Alpha = 0.25 reproduces `results_modules.csv` exactly.
- TFEB inside CLEAR: without TFEB the cohort effect is d = 1.10 (P = 1.3e-4); TFEB correlates weakly with the rest of the module (rho = 0.14).
- Non-T pool: with the monocyte share set to 40, 55 (published) and 70 percent, CLEAR reaches the observed d = 1.20 at 7.9, 6.7 and 5.7 percent non-T content; with monocytes alone, equal weights, or weights matched to the cohort's excess it is 4.5 to 5.7 percent; a pool without monocytes needs 19 percent. No pool reproduces the observed initiation decrease (largest fall d = -0.35; the equal-weight pool moves it up to d = +0.90 at 30 percent). Without TFEB the crossing is 8.1 percent in the cohort-matched pool (79 percent monocytes), 6.4 to 11.3 percent in pools with monocyte shares of 55 percent or more, and 15.2 percent in the 40 percent pool.

Outputs: `results_sensitivity_alpha.csv`, `results_sensitivity_tfeb.csv`, `results_sensitivity_nonT_pool.csv`.

## Wang and Iha gene lists

`test_wang_iha.py` read `wang_iha_degs.json` from the top folder, but the file is in `results/`; the script now reads it there. Run end to end, it reproduced `results_wangiha_ferroptosis.csv` and `results_wangiha_autophagy.csv` with a maximum difference of zero.

## Graphical abstract

`make_graphical_abstract.py` reads every plotted number from the deposited tables and writes them to `results/results_graphical_abstract_data.csv` (61 values, including the 51 gene-level points of the middle panel). An earlier draft of the middle panel had been drawn from simulated points; it was replaced before submission.

## Open point

The text of the manuscript says Wang and Iha reported 46 ferroptosis-related genes; `wang_iha_degs.json` holds 45. The difference has not been checked against the original publication.
