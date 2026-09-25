## 52_missingness_block_refit.R
## Protein missingness comes in a block: about half of the proteins are missing together in a subset of participants.
## (1) characterise the block (proteins, participants, Olink processing batch);
## (2) refit proWC with the same nested ten-fold cross-validation (LASSO over the 2,920 proteins with age and sex; inner
##     ten-fold CV for lambda.min) using only participants with at most 20% of protein values missing, i.e. with the
##     block measured, and compare the out-of-fold predictions, proWCdelta and its disease associations at fixed WC and
##     BMI with those of the main model in the same participants.
## Output: T67 (block summary), T68 (refit comparison).
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(doParallel); library(caret); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(123); registerDoParallel(cores = 24)

## (1) the missingness block
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
raw <- fread(bf, select = c(hdr[1], pcols, "Batch")); setnames(raw, 1, "id")
M <- is.na(as.matrix(raw[, ..pcols])); pm <- colMeans(M)
blk <- pcols[pm > 0.12]; say("proteins: %d; with >12%% missing (block): %d; missing fraction in block %.3f-%.3f; outside block max %.3f", length(pcols), length(blk), min(pm[blk]), max(pm[blk]), max(pm[!(pcols %in% blk)]))
say("block proteins by column position: first half %d, second half %d", sum(match(blk, pcols) <= length(pcols) / 2), sum(match(blk, pcols) > length(pcols) / 2))
fb <- rowMeans(M[, pcols %in% blk, drop = FALSE]); fo <- rowMeans(M[, !(pcols %in% blk), drop = FALSE]); pna <- rowMeans(M)
raw[, `:=`(block_missing = fb, other_missing = fo, pna = pna)]
say("participants: %d; block entirely or almost entirely missing (>90%%): %d; block complete (0%%): %d; >20%% of all proteins missing: %d", nrow(raw), sum(fb > 0.9), sum(fb == 0), sum(pna > 0.2))
bt <- raw[, .(n = .N, block_missing_gt90 = sum(block_missing > 0.9), pct = 100 * mean(block_missing > 0.9)), by = Batch][order(Batch)]; print(bt)
t67 <- rbind(data.table(item = c("proteins", "block proteins (>12% missing)", "block missing fraction min", "block missing fraction max", "max missing fraction outside block",
                                 "participants", "participants with >90% of block missing", "participants with block complete", "participants with >20% of all proteins missing"),
                        value = c(length(pcols), length(blk), min(pm[blk]), max(pm[blk]), max(pm[!(pcols %in% blk)]), nrow(raw), sum(fb > 0.9), sum(fb == 0), sum(pna > 0.2))),
             bt[, .(item = paste0("Olink batch ", Batch, ": n ", n, ", block missing (>90%)"), value = block_missing_gt90)])
fwrite(t67, file.path(TB, "T67_protein_missingness_block.csv"))

## (2) refit in participants with the block measured (<=20% of all proteins missing)
dat <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(dat, 1, "id"); prot <- setdiff(names(dat), c("id", "WC", "Age", "Sex"))
dat <- merge(dat, raw[, .(id, pna)], by = "id"); dc <- dat[pna <= 0.2]; say("refit set n = %d (of %d)", nrow(dc), nrow(dat))
X <- cbind(as.matrix(dc[, ..prot]), Age = num(dc$Age), Sex = as.integer(dc$Sex)); y <- num(dc$WC)
folds <- createFolds(y, k = 10); pred <- rep(NA_real_, length(y))
for (k in seq_along(folds)) { te <- folds[[k]]; tr <- setdiff(seq_along(y), te); t0 <- Sys.time()
  cvm <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = 1, standardize = TRUE, nfolds = 10, parallel = TRUE)
  pred[te] <- as.numeric(predict(cvm, s = "lambda.min", newx = X[te, ]))
  say("  fold %2d nzero %d %.1f min", k, cvm$nzero[which(cvm$lambda == cvm$lambda.min)], as.numeric(difftime(Sys.time(), t0, units = "mins"))) }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
e <- merge(data.table(id = dc$id, WC = y, pred_c = pred), lasso[, .(id, pred_m = pred_WC)], by = "id")
e[, `:=`(dlt_c = pred_c - fitted(lm(pred_c ~ WC)), dlt_m = pred_m - fitted(lm(pred_m ~ WC)))]
r2 <- function(p, y) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
say("same participants: R2 refit %.3f, main %.3f; r(pWC) %.4f; r(proWCdelta) %.4f", r2(e$pred_c, e$WC), r2(e$pred_m, e$WC), cor(e$pred_c, e$pred_m), cor(e$dlt_c, e$dlt_m))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "BMI", "yu_ten_need_diagnosis")); setnames(phen, c("id", "Sex", "Age", "BMI", "dx10"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
e <- merge(merge(e, phen, by = "id"), cov, by = "id"); e[, `:=`(Age = num(Age), BMI = num(BMI), Sex = as.integer(Sex))]
e <- e[complete.cases(e[, .(Age, Sex, BMI, tdi, smoking)])]; x <- as.character(e$dx10); x[is.na(x)] <- ""
dis <- data.table(disease = c("Type 2 diabetes", "Obesity diagnosis", "Dyslipidemia", "Hypertension", "Ischaemic heart disease", "Heart failure", "Liver disease", "Chronic kidney disease"), code = c("E11", "E66", "E78", "I10", "I25", "I50", "K76", "N18"))
for (i in seq_len(nrow(dis))) e[[dis$code[i]]] <- as.integer(grepl(paste0("(^|[|])", dis$code[i]), x))
e[, `:=`(zc = dlt_c / sd(dlt_c), zm = dlt_m / sd(dlt_m))]
out <- list(); for (i in seq_len(nrow(dis))) for (z in c("zm", "zc")) {
  s <- summary(glm(as.formula(paste(dis$code[i], "~", z, "+ Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)")), e, family = binomial()))$coefficients[z, ]
  out[[length(out) + 1]] <- data.table(disease = dis$disease[i], score = c(zm = "main model (all participants)", zc = "refitted without block-missing participants")[[z]], n = nrow(e), events = sum(e[[dis$code[i]]]),
                                       OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]])) }
t68 <- rbindlist(out); print(dcast(t68, disease ~ score, value.var = "OR"))
t68 <- rbind(t68, data.table(disease = c("R2, main model", "R2, refitted model", "r(pWC main, pWC refit)", "r(proWCdelta main, proWCdelta refit)"), score = "accuracy and agreement", n = nrow(e), events = NA,
                             OR = c(r2(e$pred_m, e$WC), r2(e$pred_c, e$WC), cor(e$pred_c, e$pred_m), cor(e$dlt_c, e$dlt_m)), lo = NA, hi = NA))
fwrite(t68, file.path(TB, "T68_refit_without_block_missing.csv")); say("DONE")
