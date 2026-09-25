## 54_covariate_imputation_diagnostics.R
## In the covariate file complete_data_imputed.csv (Townsend index, smoking, alcohol, television, diet items, education,
## physical activity, ethnicity, energy), missing values had been filled during data preparation by a method that was
## not documented. This script describes the filled values: for each covariate, the number of values filled, whether
## every filled value is one of the observed values, whether the filled values follow the observed distribution, and how
## closely they relate to the other covariates in the file compared with observed values.
## Output: T72.
suppressPackageStartupMessages(library(data.table))
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
imp <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7/complete_data_imputed.csv"); setnames(imp, 1, "id")
bf <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8/wc_pro_Batch.csv"; hdr <- names(fread(bf, nrows = 0))
map <- c(TDI = "Townsend deprivation index at recruitment", "Smoking status" = "Smoking status", "Alcohol intake frequency" = "Alcohol intake", Time_TV = "Time spent watching television (TV)",
         "Oily fish intake" = "Oily fish intake", "Processed meat intake" = "Processed meat intake", Education = "Education", "MET Physical activity" = "MET minutes per week for Physical activity",
         "Fruit intake" = "Fruit intake", "vegetable intake" = "vegetable intake", "Red meat intake" = "Red meat intake")
raw <- fread(bf, select = c(hdr[1], names(map))); setnames(raw, c("id", paste0("raw_", seq_along(map))))
d <- merge(raw, imp, by = "id"); fcols <- setdiff(names(imp), "id")
out <- list(); for (k in seq_along(map)) { r <- suppressWarnings(as.numeric(d[[paste0("raw_", k)]])); fc <- map[[k]]; m <- as.numeric(d[[fc]]); na <- is.na(r)
  X <- as.matrix(d[, setdiff(fcols, fc), with = FALSE]); fit <- lm.fit(cbind(1, X[!na, ]), m[!na]); pa <- as.numeric(cbind(1, X) %*% ifelse(is.na(fit$coefficients), 0, fit$coefficients))
  out[[k]] <- data.table(covariate = fc, n = nrow(d), n_imputed = sum(na), pct_imputed = 100 * mean(na), imputed_values_among_observed_pct = 100 * mean(m[na] %in% m[!na]),
    mean_observed = mean(m[!na]), mean_imputed = mean(m[na]), sd_observed = sd(m[!na]), sd_imputed = sd(m[na]),
    cor_with_other_covariates_observed = cor(m[!na], pa[!na]), cor_with_other_covariates_imputed = if (sum(na) > 2) cor(m[na], pa[na]) else NA_real_) }
t72 <- rbindlist(out); print(t72[, .(covariate = substr(covariate, 1, 28), n_imputed, pct = round(pct_imputed, 2), donor_pct = imputed_values_among_observed_pct, r_obs = round(cor_with_other_covariates_observed, 3), r_imp = round(cor_with_other_covariates_imputed, 3))])
fwrite(t72, file.path(W, "tables", "T72_covariate_imputation_diagnostics.csv")); cat("DONE\n")
