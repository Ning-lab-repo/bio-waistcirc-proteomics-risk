## 46_imputation_medians.R
## Per-protein values used for median imputation, so that proWC can be computed in other datasets with the same
## convention. The medians are recomputed from the raw NPX values of all 53,013 participants and checked against the
## imputed matrix used to fit the model (analysis_data_WC.csv): every originally missing value must equal the median.
## Imputation used the medians over all 53,013 participants with proteomic data (before the 134 exclusions).
suppressPackageStartupMessages({ library(data.table) })
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
raw <- fread(bf, select = c(hdr[1], pcols)); setnames(raw, 1, "id")
med <- sapply(pcols, function(v) median(raw[[v]], na.rm = TRUE))
imp <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(imp, 1, "id")
setnames(imp, make.names(names(imp)))
common <- pcols[make.names(pcols) %in% names(imp)]; cat("proteins in raw:", length(pcols), " in model matrix:", length(common), "\n")
raw <- raw[id %in% imp$id]; imp <- imp[match(raw$id, imp$id)]
cat("participants matched:", nrow(raw), "\n")
res <- rbindlist(lapply(common, function(v) { x <- raw[[v]]; m <- med[[v]]; na <- is.na(x)
  data.table(protein = v, n_missing = sum(na), median_npx = m, imputed_equals_median = if (any(na)) all(abs(imp[[make.names(v)]][na] - m) < 1e-9) else NA) }))
cat("proteins with missing values:", sum(res$n_missing > 0), "; imputed value equals the recomputed median for",
    sum(res$imputed_equals_median, na.rm = TRUE), "of them\n")
print(head(res[imputed_equals_median == FALSE], 5))
fwrite(res[, .(protein, n_missing_analysis_set = n_missing, median_npx_all_53013 = median_npx)], file.path(W, "tables", "T49_imputation_medians.csv"))
cat("DONE\n")
