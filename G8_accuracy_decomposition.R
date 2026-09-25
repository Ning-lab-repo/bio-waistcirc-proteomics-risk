## G8_accuracy_decomposition.R
## Why does the proteomic waist reach an out-of-fold R2 of 0.158 in the 485 Guangzhou participants against 0.791 here?
## Three candidate explanations are separated by handicapping the UK Biobank data one step at a time and refitting
## from scratch each time: the smaller sample, the small overlap between the two protein panels (which also removes
## leptin, not among the 132 shared proteins) and the much narrower spread of the waist in those 485 participants.
## What is left after all three is the platform difference itself.
## The comparison is also made in centimetres, because R2 is the share of the variance explained and is therefore not
## comparable between cohorts whose waists vary by different amounts; the mean absolute error is.
## The Supplementary Results report the three settings below (numbered 1, 2 and 5 as in earlier versions).
## Output: G8_accuracy_decomposition.csv.
suppressPackageStartupMessages({ library(data.table); library(readxl); library(glmnet) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

## ---------------- which of the informative proteins the Chinese serum panel carries at all ----------------
E <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "E", n_max = 2))
gset <- sub("^[^_]*_", "", names(E))
chk <- c("LEP", "FABP4", "GDF15", "ADM", "IGFBP1", "IGFBP2", "INHBC", "CRP", "IL6", "TNF", "ADIPOQ", "LPA", "MUC16", "SHBG")
say("informative proteins present in the Guangzhou panel: %s", paste(intersect(chk, gset), collapse = ", "))
say("absent: %s", paste(setdiff(chk, gset), collapse = ", "))

u <- fread(file.path(F12, "analysis_data_WC.csv")); prot <- setdiff(names(u), c(names(u)[1], "WC", "Age", "Sex"))
u <- u[complete.cases(u[, .(WC, Age, Sex)])]
ov <- fread(file.path(TB, "G_overlap.csv")); shared <- intersect(ov[in_ukb == TRUE]$gene, prot)
say("UK Biobank: %d participants, %d proteins, waist %.1f cm (SD %.1f); shared with Guangzhou: %d",
    nrow(u), length(prot), mean(u$WC), sd(u$WC), length(shared))
say("leptin alone with age and sex explains %.3f of the waist here; age and sex alone %.3f",
    summary(lm(WC ~ LEP + Age + factor(Sex), u))$r.squared, summary(lm(WC ~ Age + factor(Sex), u))$r.squared)

## ---------------- narrow the waist window until its spread matches those 485 participants ----------------
med <- median(u$WC); w <- 30
repeat { s <- u[abs(WC - med) <= w]; if (sd(s$WC) <= 7.6 || w < 2) break; w <- w - 0.5 }
say("restricted to a waist within %.1f cm of the median: %d participants, waist %.1f cm (SD %.1f)", w, nrow(s), mean(s$WC), sd(s$WC))

## ---------------- refit from scratch under each handicap; 71.1% women, as in the Guangzhou set ----------------
run <- function(dat, cols, tag, n = 485, fw = 0.711) {
  idx <- c(sample(which(dat$Sex == 0), round(n * fw)), sample(which(dat$Sex == 1), n - round(n * fw))); x <- dat[idx]
  X <- as.matrix(x[, ..cols]); for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); X[, j] <- v }
  X <- cbind(X, Age = as.numeric(x$Age), Sex = as.integer(x$Sex)); X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
  y <- as.numeric(x$WC); fo <- sample(rep(1:10, length.out = n)); p <- rep(NA_real_, n)
  for (k in 1:10) { tr <- fo != k; cv <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = 1, nfolds = 10)
    p[!tr] <- as.numeric(predict(cv, s = "lambda.min", newx = X[!tr, , drop = FALSE])) }
  data.table(setting = tag, n = n, proteins = ncol(X) - 2, SD_WC = sd(y),
             R2 = 1 - sum((p - y)^2) / sum((y - mean(y))^2), RMSE = sqrt(mean((p - y)^2)), MAE = mean(abs(p - y))) }

REP <- 6
r <- rbindlist(c(
  lapply(1:REP, function(i) run(u, prot,                   "1. sample size only")),
  lapply(1:REP, function(i) run(u, shared,                 "2. and only the shared proteins")),
  lapply(1:REP, function(i) run(s, shared,                 "3. narrow waist and only the shared proteins"))))
sm <- r[, .(replicates = .N, SD_WC = sprintf("%.1f", median(SD_WC)),
            R2 = sprintf("%.3f (%.3f-%.3f)", median(R2), min(R2), max(R2)),
            RMSE_cm = sprintf("%.2f", median(RMSE)), MAE_cm = sprintf("%.2f", median(MAE))), by = .(setting, n, proteins)]
print(sm, width = 220)
say("Guangzhou, for comparison: n = 485, 361 proteins, waist SD 7.6, R2 0.158, RMSE 6.99 cm, MAE 5.49 cm")
fwrite(r, file.path(TB, "G8_accuracy_decomposition.csv"))
say("DONE")
