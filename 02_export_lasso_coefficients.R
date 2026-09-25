## Export the proWC model definition given in Supplementary Tables 1 and 2:
##   (1) the LASSO coefficients of the model fitted to the full analysis set (n = 52,879), Supplementary Table 1
##   (2) the penalty, the residual-anchoring coefficients b0 and b1 and the accuracy figures, Supplementary Table 2
## The deposited model was fitted with glmnet at a single penalty, lambda = 0.01457593 (alpha = 1,
## standardize = TRUE, age and sex in the penalised design matrix after the 2,920 proteins). This value is
## point 49 of the grid 10^seq(0.2, -4, length.out = 100). Ten-fold cross-validation over the same grid in
## the full set had its minimum at the adjacent point, lambda = 0.01321941, with the same cross-validated
## mean squared error to four significant figures (38.09 cm^2). A model fitted there retains 1,637 proteins,
## and its predictions correlate with those of the deposited model at r > 0.9999.
## Every pWC value analysed in the paper is out of fold (lasso_WC.csv, from figures/Figure2__proWC_construction.R).
## The deposited coefficients are needed only to compute pWC in a new dataset. This script reproduces
## Supplementary Table 1 to within 1e-12 (floating-point rounding).
## Set RUN_CV=1 to repeat the full-set cross-validation over the grid. The minimum can move by one grid point
## with a different fold assignment.
suppressPackageStartupMessages({library(data.table); library(glmnet); library(doParallel)})
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
SRC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
set.seed(123); registerDoParallel(cores = 24)

cat("loading...\n"); flush.console()
dat <- fread(file.path(SRC, "analysis_data_WC.csv"))
setnames(dat, names(dat)[1], "participant_id")
prot <- setdiff(names(dat), c("participant_id","WC","Age","Sex"))
y <- as.numeric(dat$WC)
X <- cbind(as.matrix(dat[, ..prot]), Age = as.numeric(dat$Age), Sex = as.integer(dat$Sex))
cat("X:", dim(X), "\n"); flush.console()

grid <- 10^seq(0.2, -4, length.out = 100)
lam  <- grid[49]
fit  <- glmnet(X, y, alpha = 1, lambda = lam, standardize = TRUE)
cf <- as.matrix(coef(fit))
co <- data.table(term = rownames(cf), coefficient = as.numeric(cf[,1]))
co <- co[coefficient != 0]
n_prot <- sum(co$term %in% prot)
cat(sprintf("lambda = %.8f ; retained proteins = %d ; age coefficient = %g\n", lam, n_prot, cf["Age", 1]))
## restore the two gene symbols that R changed when the data file was read (HLA-DRA, ERVV-1)
co[term == "HLA.DRA", term := "HLA-DRA"]; co[term == "ERVV.1", term := "ERVV-1"]
co[, abs_coef := abs(coefficient)]; setorder(co, -abs_coef); co[, abs_coef := NULL]
fwrite(co, file.path(W, "tables", "SupplementaryTable_proWC_LASSO_coefficients.csv"))

if (Sys.getenv("RUN_CV") == "1") {
  t0 <- Sys.time()
  cvm <- cv.glmnet(X, y, alpha = 1, lambda = grid, standardize = TRUE, nfolds = 10, parallel = TRUE)
  cat("cv.glmnet done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  k <- which.min(cvm$cvm)
  cat(sprintf("cross-validation minimum at grid point %d, lambda = %.8f (MSE %.4f); at point 49: MSE %.4f\n",
              k, cvm$lambda[k], cvm$cvm[k], cvm$cvm[49]))
  fwrite(data.table(lambda = cvm$lambda, log10_lambda = log10(cvm$lambda), nzero = cvm$nzero,
                    cvm = cvm$cvm, cvup = cvm$cvup, cvlo = cvm$cvlo),
         file.path(W, "tables", "SupplementaryTable_LASSO_regularisation_path.csv"))
}

## residual anchoring and accuracy, from the out-of-fold predictions analysed in the paper
lw <- fread(file.path(SRC, "lasso_WC.csv")); setnames(lw, names(lw)[1], "participant_id")
afit <- lm(pred_WC ~ Actual_WC, data = lw)
r2v <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
oof <- readRDS(file.path(W, "output", "oof_models.rds"))   ## 01_oof_prowc.R; model D: age and sex only
st2 <- data.table(
  parameter = c("lambda", "n_proteins_retained", "anchor_intercept_b0", "anchor_slope_b1", "R2_pWC_out_of_fold",
                "RMSE_pWC_cm", "R2_proWC_squared_correlation", "RMSE_proWC_cm", "R2_age_sex_only"),
  value = c(lam, n_prot, coef(afit)[1], coef(afit)[2], cor(lw$pred_WC, lw$Actual_WC)^2,
            sqrt(mean((lw$Actual_WC - lw$pred_WC)^2)), cor(lw$BioX_Adjusted, lw$Actual_WC)^2,
            sqrt(mean((lw$Actual_WC - lw$BioX_Adjusted)^2)), oof$D$r2),
  note = c("penalty of the deposited model; the full-set ten-fold cross-validation minimum 0.01321941 is the adjacent grid value with the same cross-validated MSE (38.09)",
           "plus sex; the age coefficient was shrunk to zero",
           "UK Biobank estimate (re-estimate in a new cohort); proWC_delta = pWC - (b0 + b1 * measured_WC)",
           "UK Biobank estimate (re-estimate in a new cohort); proWC = measured_WC + proWC_delta = pWC + 0.211*WC - 19.10",
           "out-of-fold prediction vs measured WC",
           "out-of-fold",
           sprintf("squared correlation of proWC with measured WC; the variance-explained R2 is %.4f. Both follow from R2_pWC_out_of_fold and anchor_slope_b1 (1/(1 + k) and 1 - k, with k = b1^2 (1 - R2) / R2), so neither is a separate measure of accuracy", r2v(lw$Actual_WC, lw$BioX_Adjusted)),
           "equals the SD of proWC_delta, because proWC - measured WC = proWC_delta; not a separate measure of accuracy",
           "reference model"))
print(st2[, .(parameter, value)])
fwrite(st2, file.path(W, "tables", "SupplementaryTable_proWC_model_parameters.csv"))
anc <- data.table(term = c("intercept_b0", "slope_b1"), estimate = as.numeric(coef(afit)),
                  note = "proWC_delta = pWC - (b0 + b1 * measured_WC); proWC = measured_WC + proWC_delta")
fwrite(anc, file.path(W, "tables", "SupplementaryTable_proWC_anchoring_coefficients.csv"))   ## model_proWC_anchoring_coefficients.csv in the source data

saveRDS(list(lambda = lam, n_protein = n_prot), file.path(W, "output", "lasso_refit_summary.rds"))
cat("DONE\n")
