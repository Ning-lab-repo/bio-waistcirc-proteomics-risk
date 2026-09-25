## ---------------------------------------------------------------------------
## 01_oof_prowc.R
## Fully out-of-fold derivation of proWC / proWCdelta, plus sex-stratified
## modelling and a "no-protein" control model to test whether the reported
## proteomic-anthropometric sexual dimorphism (PASD) is a regression-to-the-
## mean artefact.
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(doParallel); library(caret)
})

WD  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
SRC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
set.seed(123)
registerDoParallel(cores = 24)

cat("== loading protein matrix ==\n"); flush.console()
dat <- fread(file.path(SRC, "analysis_data_WC.csv"))
setnames(dat, names(dat)[1], "participant_id")
cat("dim:", dim(dat), "\n")

prot_cols <- setdiff(names(dat), c("participant_id","WC","Age","Sex"))
cat("n proteins:", length(prot_cols), "\n")

y  <- as.numeric(dat$WC)
id <- dat$participant_id
sex <- as.integer(dat$Sex)
age <- as.numeric(dat$Age)

Xp <- as.matrix(dat[, ..prot_cols])
Xfull <- cbind(Xp, Age = age, Sex = sex)
rm(dat); gc()

pf_pen   <- rep(1, ncol(Xfull))                       # age/sex penalised (as originally coded)
pf_unpen <- c(rep(1, length(prot_cols)), 0, 0)        # age/sex unpenalised (as Methods claims)

folds <- createFolds(y, k = 10)

## --- generic out-of-fold runner -------------------------------------------
oof_run <- function(Xm, yv, foldlist, pf, label, anchor_oof = TRUE) {
  n <- length(yv)
  pred <- rep(NA_real_, n); delta <- rep(NA_real_, n)
  lam <- rep(NA_real_, 10); nz <- rep(NA_integer_, 10)
  for (k in seq_along(foldlist)) {
    te <- foldlist[[k]]; tr <- setdiff(seq_len(n), te)
    t0 <- Sys.time()
    cvm <- cv.glmnet(Xm[tr, , drop = FALSE], yv[tr], alpha = 1,
                     standardize = TRUE, penalty.factor = pf,
                     nfolds = 10, parallel = TRUE)
    lam[k] <- cvm$lambda.min
    nz[k]  <- cvm$nzero[which(cvm$lambda == cvm$lambda.min)]
    pr_tr <- as.numeric(predict(cvm, s = "lambda.min", newx = Xm[tr, , drop = FALSE]))
    pr_te <- as.numeric(predict(cvm, s = "lambda.min", newx = Xm[te, , drop = FALSE]))
    fit <- lm(pr_tr ~ yv[tr])                       # anchoring fitted in TRAINING folds only
    exp_te <- coef(fit)[1] + coef(fit)[2] * yv[te]
    pred[te]  <- pr_te
    delta[te] <- pr_te - exp_te
    cat(sprintf("  [%s] fold %2d  lambda=%.5f  nzero=%4d  %.1f min\n",
                label, k, lam[k], nz[k],
                as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
  }
  if (!anchor_oof) {                                 # in-sample anchoring (original behaviour)
    fitA <- lm(pred ~ yv)
    delta <- pred - (coef(fitA)[1] + coef(fitA)[2] * yv)
  }
  list(pred = pred, delta = delta, prowc = yv + delta,
       r2 = cor(yv, pred)^2, rmse = sqrt(mean((yv - pred)^2)),
       lambda = lam, nzero = nz)
}

res <- list()

cat("\n== Model A: pooled, age/sex PENALISED, anchoring in-sample (replicates published) ==\n")
res$A <- oof_run(Xfull, y, folds, pf_pen, "A", anchor_oof = FALSE)

cat("\n== Model B: pooled, age/sex UNPENALISED, anchoring out-of-fold (corrected primary) ==\n")
res$B <- oof_run(Xfull, y, folds, pf_unpen, "B", anchor_oof = TRUE)

cat("\n== Model D: NO PROTEINS (age+sex only), same pipeline -- PASD artefact control ==\n")
Xas <- cbind(Age = age, Sex = sex, dummy = rnorm(length(y), 0, 1e-8))
res$D <- oof_run(Xas, y, folds, c(0, 0, 1), "D", anchor_oof = TRUE)

cat("\n== Model C: SEX-STRATIFIED, age unpenalised, anchoring out-of-fold within sex ==\n")
predC <- rep(NA_real_, length(y)); deltaC <- rep(NA_real_, length(y))
for (s in c(0, 1)) {
  idx <- which(sex == s)
  Xs  <- cbind(Xp[idx, , drop = FALSE], Age = age[idx])
  ys  <- y[idx]
  pfs <- c(rep(1, length(prot_cols)), 0)
  fs  <- createFolds(ys, k = 10)
  rs  <- oof_run(Xs, ys, fs, pfs, paste0("C_sex", s), anchor_oof = TRUE)
  predC[idx]  <- rs$pred
  deltaC[idx] <- rs$delta
  cat(sprintf("  sex=%d  n=%d  R2=%.4f  RMSE=%.3f\n", s, length(idx), rs$r2, rs$rmse))
}
res$C <- list(pred = predC, delta = deltaC, prowc = y + deltaC,
              r2 = cor(y, predC)^2, rmse = sqrt(mean((y - predC)^2)))

out <- data.table(
  participant_id = id, WC = y, Age = age, Sex = sex,
  pWC_A = res$A$pred, proWC_A = res$A$prowc, proWCd_A = res$A$delta,
  pWC_B = res$B$pred, proWC_B = res$B$prowc, proWCd_B = res$B$delta,
  pWC_C = res$C$pred, proWC_C = res$C$prowc, proWCd_C = res$C$delta,
  pWC_D = res$D$pred, proWC_D = res$D$prowc, proWCd_D = res$D$delta
)
fwrite(out, file.path(WD, "output", "oof_prowc_all_models.csv"))

perf <- data.table(
  model = c("A_pooled_penalised_insample_anchor","B_pooled_unpenalised_oof_anchor",
            "C_sex_stratified_oof_anchor","D_age_sex_only_control"),
  r2 = sapply(res[c("A","B","C","D")], function(z) z$r2),
  rmse = sapply(res[c("A","B","C","D")], function(z) z$rmse)
)
print(perf)
fwrite(perf, file.path(WD, "tables", "T1_oof_model_performance.csv"))

## --- PASD test -------------------------------------------------------------
pasd <- rbindlist(lapply(c("A","B","C","D"), function(m) {
  d <- res[[m]]$delta
  data.table(model = m,
             mean_delta_male   = mean(d[sex == 1]),
             mean_delta_female = mean(d[sex == 0]),
             sex_offset        = mean(d[sex == 1]) - mean(d[sex == 0]),
             sd_delta          = sd(d),
             offset_in_SD      = (mean(d[sex == 1]) - mean(d[sex == 0])) / sd(d),
             pct_male_above0   = 100 * mean(d[sex == 1] > 0),
             pct_female_above0 = 100 * mean(d[sex == 0] > 0))
}))
print(pasd)
fwrite(pasd, file.path(WD, "tables", "T2_PASD_sex_offset_by_model.csv"))

saveRDS(res, file.path(WD, "output", "oof_models.rds"))
cat("\nDONE\n")
