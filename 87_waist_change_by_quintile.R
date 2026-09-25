## 87_waist_change_by_quintile.R
## Change in measured waist after baseline, by baseline proWCdelta (Supplementary Results, section 2):
##  (a) mean change from baseline to the imaging visit (instance 2) across quintiles of baseline proWCdelta, in the
##      participants re-measured there (T9; the quintile means quoted in the text);
##  (b) change from baseline to each later visit regressed on proWCdelta with baseline WC, age and sex (T8), without the
##      adjustment for the time between measurements that 29_waist_change_later_intervals.R adds; the text quotes the
##      time-adjusted slopes of 29.
## These two tables were first written by exploratory scripts that are not part of this package; this script
## reproduces them exactly. Output: T8_proWCdelta_predicts_future_WC.csv, T9_future_WC_change_by_proWCdelta_quintile.csv
suppressPackageStartupMessages(library(data.table))
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
BSZ <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/6-1Body_size_measures_participant.csv"
OUT <- if (nzchar(Sys.getenv("OUT"))) Sys.getenv("OUT") else file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, names(lasso)[1], "participant_id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC"))
setnames(phen, "Participant ID", "participant_id")
ht <- fread(BSZ, select = c("Participant ID", paste0("Waist circumference | Instance ", 1:3)))
setnames(ht, c("participant_id", "wc1", "wc2", "wc3"))
d <- merge(merge(phen, lasso[, .(participant_id, proWCd = BioX_Delta)], by = "participant_id"), ht, by = "participant_id", all.x = TRUE)
d[, `:=`(Sex = as.integer(Sex), Age = num(Age), WC = num(WC), proWCd = num(proWCd), wc1 = num(wc1), wc2 = num(wc2), wc3 = num(wc3))]

## (b) slopes on proWCdelta, baseline WC, age and sex
rep <- rbindlist(lapply(c("wc1", "wc2", "wc3"), function(v) {
  s <- d[!is.na(get(v)) & !is.na(WC) & !is.na(proWCd)]
  s[, dWC := get(v) - WC]
  f <- lm(dWC ~ proWCd + WC + Age + Sex, data = s)
  co <- summary(f)$coefficients; ci <- confint(f)
  data.table(instance = v, n = nrow(s), mean_dWC = mean(s$dWC),
             beta_per_cm_proWCd = co["proWCd", "Estimate"], lo = ci["proWCd", 1], hi = ci["proWCd", 2],
             p = co["proWCd", "Pr(>|t|)"], beta_per_SD = co["proWCd", "Estimate"] * sd(s$proWCd))
}))
print(rep); fwrite(rep, file.path(OUT, "T8_proWCdelta_predicts_future_WC.csv"))

## (a) quintiles of baseline proWCdelta among the participants re-measured at the imaging visit
q2 <- d[!is.na(wc2) & !is.na(WC) & !is.na(proWCd)]
q2[, dWC := wc2 - WC][, q := cut(proWCd, quantile(proWCd, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
qs <- q2[, .(n = .N, m = mean(dWC), se = sd(dWC) / sqrt(.N)), by = q][order(q)]
print(qs); fwrite(qs, file.path(OUT, "T9_future_WC_change_by_proWCdelta_quintile.csv"))
cat("DONE\n")
