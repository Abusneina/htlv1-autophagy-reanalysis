# Release v1.7 (3 October 2026)

Completes the Wang and Iha (2023) ferroptosis gene list. No analysis code changes.

## Correction

`results/wang_iha_degs.json` held 45 of the 46 ferroptosis-related differentially expressed genes that Wang and Iha report (21 up, 24 down instead of 21 up, 25 down). The missing gene is GABARAPL1, downregulated, log2 fold change -1.9746. It is a node of the ferroptosis protein interaction network in their Figure 6A, it is in the autophagy list of the same file with that fold change, and adding it gives exactly their reported 46 genes (21 up, 25 down). The autophagy list (26 genes, 9 up, 17 down) was already complete. Both lists were checked against the authors' supplementary tables (Supplementary Table S1 for ferroptosis, S3 for autophagy): all 46 and all 26 genes are present, every log2 fold change is identical, and the up and down counts agree. In Table S1 GABARAPL1 has an adjusted P value of 6.59e-09 and log2 fold change -1.9746.

## Effect on results

`test_wang_iha.py` was rerun unchanged with the 46 genes. `results_wangiha_ferroptosis.csv` is regenerated; `results_wangiha_autophagy.csv` is identical to the earlier version.

| Quantity (ferroptosis list) | 45 genes | 46 genes |
|---|---|---|
| Genes with absolute rho above 0.3 with the erythroid score | 96 percent (P = 2e-17) | 96 percent (P = 9.1e-18) |
| Genes with absolute rho above 0.3 with the monocyte fraction | 51 percent | 50 percent |
| Significant after adjusting on deconvolution components | 43 of 45 | 44 of 46 |
| Significant after adjusting on erythroid and monocyte content | 39 of 45 | 40 of 46 |

## Files changed

`results/wang_iha_degs.json`, `results/results_wangiha_ferroptosis.csv` and this note.

## Resolves

The open point recorded in the v1.4 release notes, that the manuscript says 46 ferroptosis genes while the file held 45.
