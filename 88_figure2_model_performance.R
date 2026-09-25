## 88_figure2_model_performance.R
## The accuracy figures printed in Figure 2A-D (source_table_Figure2_model_performance.csv): for pWC and proWC against
## measured WC, in the full cohort (out-of-fold predictions, lasso_WC.csv) and in the geographic hold-out (model trained
## in England, applied to Scotland and Wales, lasso_test.csv; both written by figures/Figure2__proWC_construction.R):
## n, R2 as the squared correlation (as printed in the figure), RMSE, and the slope and intercept of the fitted line of
## the prediction on measured WC. The first row is the age-and-sex-only reference model of 01_oof_prowc.R (model D),
## rounded as reported. This table was first written by an exploratory step outside this package; this script
## reproduces its values (to within 1e-12).
suppressPackageStartupMessages(library(data.table))
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F2 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure2"
OUT <- if (nzchar(Sys.getenv("OUT"))) Sys.getenv("OUT") else file.path(W, "tables")
full <- fread(file.path(F12, "lasso_WC.csv")); ho <- fread(file.path(F2, "lasso_test.csv"))
perf <- function(d, panel, cohort, metric, x) {
  y <- d$Actual_WC; f <- coef(lm(x ~ y))
  data.table(panel = panel, cohort = cohort, metric = metric, n = length(y), r2 = cor(x, y)^2, rmse = sqrt(mean((y - x)^2)),
             slope = unname(f[2]), intercept = unname(f[1]), note = "Computed from retained Figure 2 prediction file.")
}
oof <- readRDS(file.path(W, "output", "oof_models.rds"))
tab <- rbind(
  data.table(panel = "Reference model", cohort = "UK Biobank full cohort", metric = "Age + sex only", n = nrow(full), r2 = round(oof$D$r2, 3),
             rmse = NA_real_, slope = NA_real_, intercept = NA_real_, note = "Reported reference model for incremental R2 comparison."),
  perf(full, "Figure 2A", "UK Biobank full cohort", "pWC", full$pred_WC),
  perf(full, "Figure 2B", "UK Biobank full cohort", "proWC", full$BioX_Adjusted),
  perf(ho, "Figure 2C", "Internal UK geographic assessment", "pWC", ho$pred_WC),
  perf(ho, "Figure 2D", "Internal UK geographic assessment", "proWC", ho$BioX_Adjusted))
print(tab)
fwrite(tab, file.path(OUT, "T120_figure2_model_performance.csv"))
cat("DONE\n")
