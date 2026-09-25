## G6_sample_size_control.R
## The proteomic waist trained in the 485 Guangzhou participants with a waist in centimetres reaches an out-of-fold R2
## of 0.158, against 0.791 here. Two explanations have to be separated before anything is said about it: the Chinese
## serum measurements may carry less information about the waist than the plasma panel used here, or 485 participants
## and 361 proteins may simply be too few to train the model. This script answers that by handicapping the UK Biobank
## data in the same way and refitting from scratch: the same 485 participants drawn at random, the same out-of-fold
## procedure, first with the full panel and then with only the proteins that the two studies share. It also asks how
## well the flag-trained score of G4, which had the whole discovery cohort of 1,783 to learn from, tracks the measured
## waist in centimetres among the 485, which is the fairer estimate of what those serum proteins know about the waist.
## Output: G6_sample_size_control.csv.
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(readxl) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

## ---------------- which proteins the two studies share ----------------
ov <- fread(file.path(TB, "G_overlap.csv")); shared <- ov[in_ukb == TRUE]$gene
say("proteins measured in both studies: %d", length(shared))

## ---------------- UK Biobank, handicapped to the same size ----------------
ukb <- fread(file.path(F12, "analysis_data_WC.csv"))
idc <- names(ukb)[1]; prot <- setdiff(names(ukb), c(idc, "WC", "Age", "Sex"))
say("UK Biobank: %d participants, %d proteins; waist %.1f cm (SD %.1f)",
    nrow(ukb), length(prot), mean(ukb$WC, na.rm = TRUE), sd(ukb$WC, na.rm = TRUE))
ukb <- ukb[complete.cases(ukb[, .(WC, Age, Sex)])]

oof_r2 <- function(dat, cols, n, tag) {
  s <- dat[sample(.N, n)]
  X <- as.matrix(s[, ..cols]); for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); X[, j] <- v }
  X <- cbind(X, Age = as.numeric(s$Age), Sex = as.integer(s$Sex)); X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
  y <- as.numeric(s$WC); fold <- sample(rep(1:10, length.out = n)); p <- rep(NA_real_, n)
  for (k in 1:10) { tr <- fold != k
    cv <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = 1, nfolds = 10)
    p[!tr] <- as.numeric(predict(cv, s = "lambda.min", newx = X[!tr, , drop = FALSE])) }
  data.table(setting = tag, n = n, proteins = ncol(X) - 2, R2 = 1 - sum((p - y)^2) / sum((y - mean(y))^2),
             r = cor(p, y), MAE = mean(abs(p - y)), SD_WC = sd(y)) }

REP <- 10
out <- list()
for (i in seq_len(REP)) out[[length(out) + 1]] <- oof_r2(ukb, prot, 485, "UK Biobank, n = 485, full panel")
for (i in seq_len(REP)) out[[length(out) + 1]] <- oof_r2(ukb, intersect(prot, shared), 485, "UK Biobank, n = 485, shared proteins only")
for (i in seq_len(3))   out[[length(out) + 1]] <- oof_r2(ukb, prot, 1783, "UK Biobank, n = 1783, full panel")
r <- rbindlist(out)
sm <- r[, .(replicates = .N, R2 = sprintf("%.3f (%.3f-%.3f)", median(R2), min(R2), max(R2)),
            r = sprintf("%.3f", median(r)), MAE_cm = sprintf("%.2f", median(MAE)), SD_WC = sprintf("%.1f", median(SD_WC))), by = .(setting, n, proteins)]
print(sm, width = 200)
fwrite(r, file.path(TB, "G6_sample_size_control.csv"))

## ---------------- how well does the larger flag-trained score track the waist in centimetres? ----------------
E <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "E", na = c("", "NA", "NaN")))
sc <- fread(file.path(TB, "G_score_discovery.csv"))
cm <- merge(E[, .(Patient_ID = pat_ID, wc = as.numeric(wc), BMI = as.numeric(BMI), sex = as.integer(sex))],
            sc[, .(Patient_ID, pCO)], by = "Patient_ID")
say("the flag-trained score of G4 against the measured waist in the %d participants who have both:", nrow(cm))
say("  Pearson %.3f, Spearman %.3f; R2 of a linear fit %.3f; adding age and sex %.3f",
    cor(cm$pCO, cm$wc), cor(cm$pCO, cm$wc, method = "spearman"),
    summary(lm(wc ~ pCO, cm))$r.squared, summary(lm(wc ~ pCO + sex, cm))$r.squared)
say("  for comparison, BMI against the measured waist in the same people: Pearson %.3f, R2 %.3f",
    cor(cm$BMI, cm$wc), summary(lm(wc ~ BMI, cm))$r.squared)
say("  the waist in these 485 has SD %.1f cm, against %.1f cm in the UK Biobank cohort", sd(cm$wc), sd(ukb$WC, na.rm = TRUE))
say("DONE")
