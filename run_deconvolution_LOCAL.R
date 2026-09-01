# =============================================================================
# Reference-based deconvolution of GSE33615
# RUN THIS ON YOUR OWN MACHINE. It needs internet access, which the analysis
# session did not have.
#
# What it does: estimates cell-type proportions in every GSE33615 sample using
# three independent reference-based methods, writes them to CSV, tests whether
# they differ between ATL and control samples, and re-tests the four autophagy
# modules with the estimated proportions as covariates.
#
# Expected runtime: 10 to 25 minutes, most of it xCell.
# =============================================================================

## ---- 1. Install once ---------------------------------------------------------
# Run these lines the first time only, then comment them out.
# install.packages(c("BiocManager", "remotes", "data.table"))
# BiocManager::install(c("GEOquery", "limma", "GSVA", "preprocessCore"))
# remotes::install_github("dviraran/xCell")          # xCell
# remotes::install_github("ebecht/MCPcounter",
#                         subdir = "Source")          # MCP-counter
# remotes::install_github("omnideconv/immunedeconv")  # wrapper, optional

library(GEOquery); library(xCell); library(MCPcounter); library(data.table)

OUT <- "deconvolution"; dir.create(OUT, showWarnings = FALSE)

## ---- 2. Expression matrix ----------------------------------------------------
# Point this at the file you already downloaded, or let GEOquery fetch it.
gse  <- getGEO("GSE33615", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset <- gse[[1]]

X <- exprs(eset)
if (max(X, na.rm = TRUE) > 50) X <- log2(X + 1)

# collapse probes to gene symbols by maximum mean intensity, as in the paper
fd  <- fData(eset)
sym <- fd[[intersect(c("GENE_SYMBOL", "Gene symbol", "Gene.symbol"), colnames(fd))[1]]]
sym <- sub(" ?//.*$", "", sym)
keep <- !is.na(sym) & sym != ""
X <- X[keep, ]; sym <- sym[keep]
ord <- order(rowMeans(X), decreasing = TRUE)
X <- X[ord, ]; sym <- sym[ord]
X <- X[!duplicated(sym), ]
rownames(X) <- sym[!duplicated(sym)]
cat("gene-level matrix:", dim(X), "\n")

## ---- 3. Group labels ---------------------------------------------------------
pd <- pData(eset)
dz <- apply(pd[, grep("characteristics", colnames(pd)), drop = FALSE], 1,
            function(r) paste(r, collapse = " | "))
group <- ifelse(grepl("ATL", dz, ignore.case = TRUE), "ATL", "control")
cat(table(group), "\n")
stopifnot(length(group) == ncol(X))

## ---- 4. xCell ----------------------------------------------------------------
# 64 cell types, enrichment-based. Scores are relative, not fractions.
xc <- xCellAnalysis(X)
fwrite(data.table(cell_type = rownames(xc), xc), file.path(OUT, "xcell_scores.csv"))

## ---- 5. MCP-counter ----------------------------------------------------------
# 10 populations, abundance estimates on a common arbitrary scale.
mcp <- MCPcounter.estimate(X, featuresType = "HUGO_symbols")
fwrite(data.table(cell_type = rownames(mcp), mcp), file.path(OUT, "mcpcounter_scores.csv"))

## ---- 6. Optional: CIBERSORTx --------------------------------------------------
# CIBERSORTx requires a free academic account and runs on their server:
#   https://cibersortx.stanford.edu
# Upload a tab-delimited file of the matrix below, choose the LM22 signature,
# select "relative mode", 100 permutations, and disable quantile normalization
# (the data are already normalized microarray values).
write.table(cbind(Gene = rownames(X), X), file.path(OUT, "GSE33615_for_CIBERSORTx.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 7. Which populations differ between the two specimen classes? -----------
test_block <- function(M, label) {
  res <- data.table(cell_type = rownames(M))
  res$mean_ATL     <- rowMeans(M[, group == "ATL", drop = FALSE])
  res$mean_control <- rowMeans(M[, group == "control", drop = FALSE])
  res$difference   <- res$mean_ATL - res$mean_control
  res$p <- apply(M, 1, function(v)
    tryCatch(wilcox.test(v[group == "ATL"], v[group == "control"])$p.value,
             error = function(e) NA_real_))
  res$q <- p.adjust(res$p, method = "BH")
  res <- res[order(q)]
  fwrite(res, file.path(OUT, paste0(label, "_group_test.csv")))
  cat("\n== ", label, ": populations with q < 0.05 ==\n", sep = "")
  print(head(res[q < 0.05], 15))
  res
}
r_xc  <- test_block(xc,  "xcell")
r_mcp <- test_block(mcp, "mcpcounter")

## ---- 8. Re-test the autophagy modules with estimated composition as covariates
mods <- fread("modules_frozen_v1.tsv")
mods <- mods[module != "axis"]

# module scores, same specification as the paper (ssGSEA, exponent 0.25)
library(GSVA)
gs <- split(mods$gene, mods$module)
gs <- lapply(gs, function(g) intersect(g, rownames(X)))
S  <- gsva(ssgseaParam(X, gs, alpha = 0.25, normalize = TRUE))

# use the leading composition axes rather than every score, to limit collinearity
# (this addresses the overadjustment concern raised in review)
comp <- t(mcp)
pc   <- prcomp(scale(comp))
npc  <- max(2, which(cumsum(pc$sdev^2) / sum(pc$sdev^2) > 0.80)[1])
cat("\nusing", npc, "composition principal components\n")
PCs <- pc$x[, 1:npc, drop = FALSE]

is_atl <- as.numeric(group == "ATL")
out <- rbindlist(lapply(rownames(S), function(m) {
  y <- S[m, ]
  raw <- summary(lm(y ~ is_atl))$coefficients["is_atl", ]
  adj <- summary(lm(y ~ is_atl + PCs))$coefficients["is_atl", ]
  data.table(module = m,
             beta_unadjusted = raw[1], p_unadjusted = raw[4],
             beta_adjusted   = adj[1], p_adjusted   = adj[4])
}))
fwrite(out, file.path(OUT, "module_scores_composition_adjusted.csv"))
cat("\n== module effect before and after adjusting for composition PCs ==\n")
print(out)

## ---- 9. Variance inflation, to document the collinearity concern -------------
if (requireNamespace("car", quietly = TRUE)) {
  fit_all <- lm(S[1, ] ~ is_atl + t(mcp))
  cat("\nVIFs for the full covariate set:\n"); print(car::vif(fit_all))
} else {
  cat("\ninstall.packages('car') to report variance inflation factors\n")
}

cat("\nDone. Send me the contents of the", OUT, "folder.\n")
