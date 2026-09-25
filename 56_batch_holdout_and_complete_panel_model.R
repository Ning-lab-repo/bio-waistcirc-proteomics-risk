## 56_batch_holdout_and_complete_panel_model.R
## (1) leave-one-Olink-batch-out validation: for each processing batch (1-7; batch 0 has 69 participants and is kept in
##     training), the LASSO model (2,920 proteins, age and sex; lambda.min by ten-fold CV within the training batches)
##     is fitted without that batch and applied to it; R2 and RMSE in the held-out batch, and agreement with the main
##     out-of-fold predictions.
## (2) a complete-panel version of the model for use in other cohorts: fitted to the 44,073 participants with at most 20%
##     of protein values missing (lambda.min by ten-fold CV), with its anchoring coefficients.
## Output: T78 (batch hold-out), model_proWC_complete_panel_coefficients.csv and model_proWC_complete_panel_anchoring.csv.
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(doParallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(123); registerDoParallel(cores = 24)
dat <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(dat, 1, "id"); prot <- setdiff(names(dat), c("id", "WC", "Age", "Sex"))
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
raw <- fread(bf, select = c(hdr[1], pcols, "Batch")); setnames(raw, 1, "id"); raw[, pna := rowMeans(is.na(as.matrix(.SD))), .SDcols = pcols]; raw <- raw[, .(id, Batch, pna)]
dat <- merge(dat, raw, by = "id"); lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id"); dat <- merge(dat, lasso[, .(id, pred_main = pred_WC)], by = "id")
X <- cbind(as.matrix(dat[, ..prot]), Age = num(dat$Age), Sex = as.integer(dat$Sex)); y <- num(dat$WC)
r2 <- function(p, yy) 1 - sum((yy - p)^2) / sum((yy - mean(yy))^2)
out <- list()
for (b in 1:7) { te <- which(dat$Batch == b); tr <- setdiff(seq_along(y), te); t0 <- Sys.time()
  cvm <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = 1, standardize = TRUE, nfolds = 10, parallel = TRUE)
  p <- as.numeric(predict(cvm, s = "lambda.min", newx = X[te, ]))
  out[[b]] <- data.table(held_out_batch = b, n = length(te), R2 = r2(p, y[te]), RMSE = sqrt(mean((y[te] - p)^2)), R2_main_oof_same_participants = r2(dat$pred_main[te], y[te]),
                         r_with_main_oof = cor(p, dat$pred_main[te]), nzero = cvm$nzero[which(cvm$lambda == cvm$lambda.min)])
  say("batch %d held out (n=%d): R2 %.3f (main OOF %.3f), r with main %.3f, %.1f min", b, length(te), out[[b]]$R2, out[[b]]$R2_main_oof_same_participants, out[[b]]$r_with_main_oof, as.numeric(difftime(Sys.time(), t0, units = "mins"))) }
t78 <- rbindlist(out); print(t78); fwrite(t78, file.path(TB, "T78_leave_one_batch_out.csv"))

## complete-panel model
cp <- which(dat$pna <= 0.2); say("complete-panel set n = %d", length(cp))
cvm <- cv.glmnet(X[cp, ], y[cp], family = "gaussian", alpha = 1, standardize = TRUE, nfolds = 10, parallel = TRUE)
cf <- as.matrix(coef(cvm, s = "lambda.min")); cfd <- data.table(term = rownames(cf), coefficient = cf[, 1])[coefficient != 0]
cfd[term == "HLA.DRA", term := "HLA-DRA"]; cfd[term == "ERVV.1", term := "ERVV-1"]   ## gene symbols as in the other tables
say("complete-panel model: lambda %.4f, %d non-zero terms (incl. intercept)", cvm$lambda.min, nrow(cfd))
p <- as.numeric(predict(cvm, s = "lambda.min", newx = X[cp, ])); a <- coef(lm(p ~ y[cp]))
fwrite(cfd, file.path(W, "output", "model_proWC_complete_panel_coefficients.csv"))
fwrite(data.table(parameter = c("lambda", "anchoring intercept b0", "anchoring slope b1", "n"), value = c(cvm$lambda.min, a[[1]], a[[2]], length(cp))), file.path(W, "output", "model_proWC_complete_panel_anchoring.csv"))
full <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/04_supplementary_tables/SupplementaryTable1_proWC_LASSO_coefficients.csv"); setnames(full, c("term", "coef_full"))
mm <- merge(cfd, full, by = "term", all = TRUE); mm[is.na(mm)] <- 0; say("coefficients: complete-panel vs deposited full model r = %.3f; terms in both %d", cor(mm$coefficient, mm$coef_full), sum(mm$coefficient != 0 & mm$coef_full != 0))
say("DONE")
