## 66_fold_internal_imputation.R
## The main model replaces missing protein values by the median of each protein over all participants, once, before the
## cross-validation folds are formed. This script repeats the nested cross-validation with the medians computed inside
## each training fold and applied to the held-out fold, so that no information from the held-out participants enters
## the imputation, and compares the resulting proWCdelta and its disease associations with the main ones.
## The ten outer folds are run as separate processes and then combined:
##   for k in 1 2 3 4 5 6 7 8 9 10; do Rscript 66_fold_internal_imputation.R $k & done; wait
##   Rscript 66_fold_internal_imputation.R combine
## Output: T100 (accuracy and agreement), T101 (disease associations).
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(doParallel); library(caret); library(splines) })
arg <- commandArgs(trailingOnly = TRUE)[1]; stopifnot(arg %in% c(as.character(1:10), "combine"))
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables"); O <- file.path(W, "output", "fold_imputation")
dir.create(O, recursive = TRUE, showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
## the analysis set and the outcome, as in the main model
dat <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(dat, 1, "id"); prot <- setdiff(names(dat), c("id", "WC", "Age", "Sex"))
keep <- dat[, .(id, WC = num(WC), Age = num(Age), Sex = as.integer(Sex))]
## the raw NPX values, with missing values still missing
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
raw <- fread(bf, select = c(hdr[1], pcols)); setnames(raw, 1, "id")
stopifnot(all(make.names(pcols) %in% c(prot, make.names(prot))) || length(pcols) == length(prot))
setnames(raw, c("id", make.names(pcols)))
d <- merge(keep, raw, by = "id"); pv <- setdiff(names(d), c("id", "WC", "Age", "Sex")); say("n = %d; proteins = %d", nrow(d), length(pv))
set.seed(123); folds <- createFolds(d$WC, k = 10)

if (arg != "combine") {
  k <- as.integer(arg); registerDoParallel(cores = 5)
  te <- folds[[k]]; tr <- setdiff(seq_len(nrow(d)), te); t0 <- Sys.time()
  X <- as.matrix(d[, ..pv])
  med <- apply(X[tr, , drop = FALSE], 2, median, na.rm = TRUE)          ## medians from the training fold only
  for (j in seq_along(pv)) { m <- is.na(X[, j]); if (any(m)) X[m, j] <- med[j] }
  X <- cbind(X, Age = d$Age, Sex = d$Sex)
  set.seed(2000 + k)
  cvm <- cv.glmnet(X[tr, , drop = FALSE], d$WC[tr], alpha = 1, standardize = TRUE, nfolds = 10, parallel = TRUE)
  p <- as.numeric(predict(cvm, s = "lambda.min", newx = X[te, , drop = FALSE]))
  fwrite(data.table(id = d$id[te], WC = d$WC[te], pWC_fold = p), file.path(O, sprintf("pred_fold%02d.csv", k)))
  fwrite(data.table(fold = k, n_test = length(te), lambda = cvm$lambda.min, nzero = cvm$nzero[which(cvm$lambda == cvm$lambda.min)], minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))), file.path(O, sprintf("info_fold%02d.csv", k)))
  say("fold %d done in %.1f min", k, as.numeric(difftime(Sys.time(), t0, units = "mins"))); quit(save = "no")
}

## ---------------- combine and compare ----------------
fl <- file.path(O, sprintf("pred_fold%02d.csv", 1:10)); stopifnot(all(file.exists(fl)))
pr <- rbindlist(lapply(fl, fread)); info <- rbindlist(lapply(file.path(O, sprintf("info_fold%02d.csv", 1:10)), fread))
a <- lm(pWC_fold ~ WC, pr); pr[, dlt_fold := pWC_fold - fitted(a)]
r2 <- 1 - sum((pr$WC - pr$pWC_fold)^2) / sum((pr$WC - mean(pr$WC))^2)
main <- fread(file.path(F12, "lasso_WC.csv")); setnames(main, 1, "id"); main <- main[, .(id, dlt_main = BioX_Delta)]
m <- merge(pr, main, by = "id")
t100 <- data.table(quantity = c("out-of-fold R2 with fold-internal imputation", "out-of-fold R2 of the main model", "median number of non-zero terms", "median lambda", "correlation of the two proWCdelta values", "mean absolute difference (cm)"),
  value = c(r2, 0.791, median(info$nzero), median(info$lambda), cor(m$dlt_fold, m$dlt_main), mean(abs(m$dlt_fold - m$dlt_main))))
print(t100); fwrite(t100, file.path(TB, "T100_fold_internal_imputation.csv"))

## disease associations in participants free of each endpoint, as in the primary analysis
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
x <- Reduce(function(p, q) merge(p, q, by = "id", all.x = TRUE), list(merge(phen, m[, .(id, dlt_fold, dlt_main)], by = "id"), cov, fr))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) x[[v]] <- num(x[[v]]); x[, Sex := as.integer(Sex)]
x <- x[complete.cases(x[, .(Age, Sex, WC, BMI, tdi, smoking, dlt_fold, dlt_main)])]
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
x[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
x[, scr := cre / 88.4]; x[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
x[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(v, pat) grepl(pat, v, perl = TRUE)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(x$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(x$med, "insulin") | has(x$medname, gd) | (!is.na(x$hba1c) & x$hba1c >= 48),
  I50 = has(x$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(x$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(x$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(x$egfr) & x$egfr < 60))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease"), code = c("E11", "I50", "K76", "N18"))
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
x[, `:=`(zf = dlt_fold / sd(dlt_fold), zm = dlt_main / sd(dlt_main))]
o <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", cd); inc <- grepl(pt, x$dx10) | grepl(pt, x$dxmore)
  x[[cd]] <- as.integer(grepl(pt, x$dx10)); fr_ok <- !((grepl(pt, x$dxall) & !inc) | extra[[cd]]); y <- x[fr_ok]
  for (zz in c("zf", "zm")) { cf <- summary(glm(as.formula(paste(cd, "~", zz, "+", BASE)), y, family = binomial()))$coefficients[zz, ]
    o[[length(o) + 1]] <- data.table(disease = dis$disease[i], score = c(zf = "fold-internal imputation", zm = "main model")[[zz]], n = nrow(y), events = sum(y[[cd]]),
      OR = exp(cf[[1]]), lo = exp(cf[[1]] - 1.96 * cf[[2]]), hi = exp(cf[[1]] + 1.96 * cf[[2]])) } }
t101 <- rbindlist(o); print(t101[, .(disease, score, n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t101, file.path(TB, "T101_fold_internal_imputation_associations.csv")); say("DONE")
