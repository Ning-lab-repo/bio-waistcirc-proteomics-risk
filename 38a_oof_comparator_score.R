## 38a_oof_comparator_score.R
## Comparator proteomic scores built exactly like proWC: LASSO over the same 2,920 proteins with age and sex
## (all penalised), out-of-fold predictions from ten-fold cross-validation with the penalty chosen by an inner
## ten-fold cross-validation in each training fold (lambda.min), and a discordance measure from anchoring the
## out-of-fold prediction on the measured trait in the full sample (as for proWCdelta).
## Usage: Rscript 38a_oof_comparator_score.R BMI   |   Rscript 38a_oof_comparator_score.R BFP
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(doParallel); library(caret) })
trait <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(trait %in% c("BMI", "BFP"))
SRC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
OUT <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/output/comparator_scores"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(123); registerDoParallel(cores = 24)
num <- function(x) suppressWarnings(as.numeric(x))

dat <- fread(file.path(SRC, "analysis_data_WC.csv")); setnames(dat, 1, "id")
prot <- setdiff(names(dat), c("id", "WC", "Age", "Sex"))
if (trait == "BMI") {
  ph <- fread(file.path(SRC, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "BMI")); setnames(ph, c("id", "y"))
} else {
  ph <- fread(file.path(RAW, "bodycomp_subset.csv"), select = c("Participant ID", "Body fat percentage | Instance 0")); setnames(ph, c("id", "y"))
}
ph[, y := num(y)]
dat[, y := ph$y[match(id, ph$id)]]
dat <- dat[!is.na(y)]
cat(trait, ": n =", nrow(dat), " proteins =", length(prot), "\n"); flush.console()
yv <- dat$y
X <- cbind(as.matrix(dat[, ..prot]), Age = num(dat$Age), Sex = as.integer(dat$Sex))
pf <- rep(1, ncol(X))
folds <- createFolds(yv, k = 10)
pred <- rep(NA_real_, length(yv)); lam <- nz <- rep(NA_real_, 10)
for (k in seq_along(folds)) {
  te <- folds[[k]]; tr <- setdiff(seq_along(yv), te); t0 <- Sys.time()
  cvm <- cv.glmnet(X[tr, , drop = FALSE], yv[tr], alpha = 1, standardize = TRUE, penalty.factor = pf, nfolds = 10, parallel = TRUE)
  pred[te] <- as.numeric(predict(cvm, s = "lambda.min", newx = X[te, , drop = FALSE]))
  lam[k] <- cvm$lambda.min; nz[k] <- cvm$nzero[which(cvm$lambda == cvm$lambda.min)]
  cat(sprintf("  [%s] fold %2d lambda=%.5f nzero=%4d %.1f min\n", trait, k, lam[k], nz[k], as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
a <- lm(pred ~ yv)
res <- data.table(id = dat$id, measured = yv, pred = pred, delta = pred - fitted(a))
m0 <- lm(yv ~ num(dat$Age) + as.integer(dat$Sex))
perf <- data.table(trait = trait, n = length(yv), r2_oof = 1 - sum((yv - pred)^2) / sum((yv - mean(yv))^2),
                   r2_age_sex = summary(m0)$r.squared, anchor_b0 = coef(a)[1], anchor_b1 = coef(a)[2],
                   median_lambda = median(lam), median_nzero = median(nz))
print(perf)
fwrite(res, file.path(OUT, paste0("oof_", trait, ".csv"))); fwrite(perf, file.path(OUT, paste0("perf_", trait, ".csv")))
file.create(file.path(OUT, paste0("DONE_", trait)))
cat("DONE\n")
