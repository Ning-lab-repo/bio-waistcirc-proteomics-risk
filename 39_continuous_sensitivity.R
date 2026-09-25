## 39_continuous_sensitivity.R
## Sensitivity analyses for the continuous proWCdelta models (per SD, at fixed WC and BMI; covariates age, sex,
## deprivation, smoking, natural splines of WC and BMI) across the eight obesity-related endpoints:
##  main      logistic, 10-year first diagnosis, prevalent cases counted as non-events (as in the main analysis)
##  noprev    logistic, participants with a hospital-recorded diagnosis before baseline excluded
##  cox       Cox, time to first diagnosis, censored at death or 10 years, prevalent cases excluded
##  finegray  Fine-Gray, death as a competing risk, prevalent cases excluded
##  cox_2y    Cox as above, excluding events and censoring in the first 2 years (reverse causation)
##  sexspec   logistic as main, with proWCdelta from the sex-specific score
## Output: T38.
suppressPackageStartupMessages({ library(data.table); library(splines); library(survival) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
ss <- fread(file.path(W, "output", "oof_prowc_all_models.csv"), select = c("participant_id", "proWCd_C")); setnames(ss, c("id", "dlt_sex"))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "yu_ten_need_time", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "date_attending_assessment_centre", "Date_death_instance0"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "t10", "dxmore", "dxall", "d0", "ddeath"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- merge(merge(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), ss, by = "id"), cov, by = "id")
d[, `:=`(Sex = as.integer(Sex), Age = num(Age), WC = num(WC), BMI = num(BMI))]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt, dlt_sex)])]
d[, `:=`(z = dlt / sd(dlt), zs = dlt_sex / sd(dlt_sex), b0 = as.IDate(d0), dd = as.IDate(ddeath))]
s <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = s(dx10), t10 = s(t10), dxmore = s(dxmore), dxall = s(dxall))]
cat("n =", nrow(d), "\n")
dis <- data.table(disease = c("Type 2 diabetes", "Obesity", "Dyslipidemia", "Hypertension", "Ischaemic heart disease", "Heart failure", "Liver disease", "Chronic kidney disease"),
                  code = c("E11", "E66", "E78", "I10", "I25", "I50", "K76", "N18"))
first_time <- function(codes, times, code) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]; k <- which(startsWith(cs, code)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(ts[k])), na.rm = TRUE) }
adj <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
row <- function(an, dz, est, lo, hi, n, ev) data.table(analysis = an, disease = dz, estimate = est, lo = lo, hi = hi, n = n, events = ev)
out <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; nm <- dis$disease[i]; pat <- paste0("(^|[|])", cd)
  d[, y := as.integer(grepl(pat, dx10))]
  g <- summary(glm(as.formula(paste("y ~ z +", adj)), d, family = binomial()))$coefficients["z", ]
  out[[length(out) + 1]] <- row("main", nm, exp(g[1]), exp(g[1] - 1.96 * g[2]), exp(g[1] + 1.96 * g[2]), nrow(d), sum(d$y))
  g <- summary(glm(as.formula(paste("y ~ zs +", adj)), d, family = binomial()))$coefficients["zs", ]
  out[[length(out) + 1]] <- row("sexspec", nm, exp(g[1]), exp(g[1] - 1.96 * g[2]), exp(g[1] + 1.96 * g[2]), nrow(d), sum(d$y))
  incl <- grepl(pat, d$dx10) | grepl(pat, d$dxmore); prev <- grepl(pat, d$dxall) & !incl
  x <- d[!prev]
  g <- summary(glm(as.formula(paste("y ~ z +", adj)), x, family = binomial()))$coefficients["z", ]
  out[[length(out) + 1]] <- row("noprev", nm, exp(g[1]), exp(g[1] - 1.96 * g[2]), exp(g[1] + 1.96 * g[2]), nrow(x), sum(x$y))
  ft <- mapply(first_time, x$dx10, x$t10, MoreArgs = list(code = cd))
  x[, tev := (ft - as.numeric(b0)) / 365.25]; x[, tdth := (as.numeric(dd) - as.numeric(b0)) / 365.25]
  x[, time := pmin(fifelse(y == 1, tev, 10), fifelse(is.na(tdth), 10, tdth), 10, na.rm = TRUE)]
  x[, status := fifelse(y == 1 & tev <= time + 1e-9, 1L, fifelse(!is.na(tdth) & tdth <= 10 & tdth <= time + 1e-9, 2L, 0L))]
  x <- x[time > 0]
  cx <- summary(coxph(as.formula(paste("Surv(time, status == 1) ~ z +", adj)), x))$conf.int["z", ]
  out[[length(out) + 1]] <- row("cox", nm, cx[1], cx[3], cx[4], nrow(x), sum(x$status == 1))
  x2 <- x[time > 2]
  cx2 <- summary(coxph(as.formula(paste("Surv(time, status == 1) ~ z +", adj)), x2))$conf.int["z", ]
  out[[length(out) + 1]] <- row("cox_2y", nm, cx2[1], cx2[3], cx2[4], nrow(x2), sum(x2$status == 1))
  x[, st := factor(status, 0:2, c("censor", "event", "death"))]
  fg <- finegray(Surv(time, st) ~ ., data = x[, .(time, st, z, Age, Sex, tdi, smoking, WC, BMI)], etype = "event")
  fm <- summary(coxph(Surv(fgstart, fgstop, fgstatus) ~ z + Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3), data = fg, weights = fgwt))$conf.int["z", ]
  out[[length(out) + 1]] <- row("finegray", nm, fm[1], fm[3], fm[4], nrow(x), sum(x$status == 1))
  cat(sprintf("%-24s done\n", nm)); flush.console()
}
res <- rbindlist(out)
print(dcast(res[, .(disease, analysis, v = sprintf("%.2f (%.2f-%.2f)", estimate, lo, hi))], disease ~ analysis, value.var = "v"), width = 300)
fwrite(res, file.path(W, "tables", "T38_continuous_proWCdelta_sensitivity.csv"))
cat("DONE\n")
