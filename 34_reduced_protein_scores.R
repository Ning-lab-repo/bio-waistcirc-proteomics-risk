## How many proteins are needed? Reduced proWC scores built from the top k proteins.
## Protein selection and model fitting are done inside each training fold (10-fold CV), so the
## out-of-fold predictions are free of selection leakage. For each fold the LASSO path is fitted
## (age and sex unpenalised); for each k the first k proteins to enter the path are kept (ties broken
## by standardised coefficient) and the model is refitted by ordinary least squares on age, sex and
## those k proteins.
suppressPackageStartupMessages({library(data.table); library(glmnet); library(caret)})
set.seed(20260922)
SRC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; OUT <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/output/reduced_scores"; dir.create(OUT, showWarnings=FALSE)
dat <- fread(file.path(SRC, "analysis_data_WC.csv")); setnames(dat, 1, "id")
prot <- setdiff(names(dat), c("id","WC","Age","Sex")); P <- length(prot)
y <- as.numeric(dat$WC); age <- as.numeric(dat$Age); sex <- as.integer(dat$Sex); id <- dat$id
Xp <- as.matrix(dat[, ..prot]); rm(dat); gc()
X <- cbind(Xp, Age=age, Sex=sex); pf <- c(rep(1,P),0,0)
K <- c(5,10,20,50,100,200,500)
folds <- createFolds(y, k=10)
pick <- function(fit, Xtr, k) {
  B <- as.matrix(fit$beta[1:P, , drop=FALSE]); nz <- colSums(B != 0)
  j <- which(nz >= k)[1]; if (is.na(j)) j <- ncol(B)
  b <- B[, j]; cand <- which(b != 0)
  sdx <- apply(Xtr[, cand, drop=FALSE], 2, sd)
  cand[order(-abs(b[cand]*sdx))][seq_len(min(k, length(cand)))] }
ols <- function(Ztr, ytr) { m <- lm.fit(cbind(1, Ztr), ytr); cf <- m$coefficients; cf[is.na(cf)] <- 0; cf }
pred <- matrix(NA_real_, length(y), length(K), dimnames=list(NULL, paste0("k", K)))
pred_as <- rep(NA_real_, length(y)); selected <- list()
for (f in seq_along(folds)) {
  te <- folds[[f]]; tr <- setdiff(seq_along(y), te); t0 <- Sys.time()
  fit <- glmnet(X[tr,], y[tr], alpha=1, penalty.factor=pf, nlambda=300, lambda.min.ratio=1e-4, dfmax=650)
  cf <- ols(cbind(age[tr], sex[tr]), y[tr]); pred_as[te] <- cbind(1, age[te], sex[te]) %*% cf
  for (j in seq_along(K)) {
    keep <- pick(fit, Xp[tr,], K[j])
    cf <- ols(cbind(age[tr], sex[tr], Xp[tr, keep, drop=FALSE]), y[tr])
    pred[te, j] <- cbind(1, age[te], sex[te], Xp[te, keep, drop=FALSE]) %*% cf
    selected[[length(selected)+1]] <- data.table(fold=f, k=K[j], protein=prot[keep]) }
  cat(sprintf("fold %d done in %.1f min\n", f, as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console()
}
r2 <- function(p) 1 - sum((y-p)^2)/sum((y-mean(y))^2)
rmse <- function(p) sqrt(mean((y-p)^2))
full <- fread(file.path(SRC,"lasso_WC.csv")); setnames(full,1,"id"); full <- full[match(id, full$id)]
perf <- data.table(k=c(0, K, 1549), label=c("age and sex", paste(K,"proteins"), "full model (1,549 proteins)"),
                   R2=c(r2(pred_as), apply(pred,2,r2), r2(full$pred_WC)), RMSE=c(rmse(pred_as), apply(pred,2,rmse), rmse(full$pred_WC)))
print(perf); fwrite(perf, file.path(OUT,"topk_performance.csv"))
sel <- rbindlist(selected); stab <- sel[, .(folds_selected=.N), by=.(k, protein)]
fwrite(stab, file.path(OUT,"topk_selection_stability.csv"))
cat("k=20: proteins selected in all 10 folds:", stab[k==20 & folds_selected==10, .N], "; in >=8 folds:", stab[k==20 & folds_selected>=8, .N], "\n")
## deployable reduced scores: selection and OLS refit in the full data
fitall <- glmnet(X, y, alpha=1, penalty.factor=pf, nlambda=300, lambda.min.ratio=1e-4, dfmax=650)
for (kk in c(20, 50)) { keep <- pick(fitall, Xp, kk); cf <- ols(cbind(age, sex, Xp[, keep, drop=FALSE]), y)
  fwrite(data.table(term=c("(Intercept)","Age","Sex", prot[keep]), coefficient=as.numeric(cf)), file.path(OUT, sprintf("reduced_score_%d_coefficients.csv", kk))) }
oof <- data.table(id=id, WC=y, pWC_k20=pred[, "k20"], pWC_k50=pred[, "k50"], pWC_full=full$pred_WC, dlt_full=full$BioX_Delta)
for (v in c("k20","k50")) { a <- lm(oof[[paste0("pWC_",v)]] ~ oof$WC); oof[[paste0("dlt_",v)]] <- oof[[paste0("pWC_",v)]] - fitted(a) }
cat("cor(proWCdelta_k20, full) =", round(cor(oof$dlt_k20, oof$dlt_full),3), "; cor(k50, full) =", round(cor(oof$dlt_k50, oof$dlt_full),3), "\n")
fwrite(oof, file.path(OUT,"topk_oof_predictions.csv"))
