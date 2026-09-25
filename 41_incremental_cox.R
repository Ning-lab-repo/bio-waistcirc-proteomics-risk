## 41_incremental_cox.R
## Incremental prediction by proWCdelta with time-to-event models (replaces the class-weighted logistic analysis).
## Endpoints: type 2 diabetes, chronic kidney disease, heart failure, composite cardiovascular disease (any I00-I99 or G45)
## and all-cause death, each within 10 years of baseline; participants with a hospital-recorded diagnosis before baseline
## are excluded for the disease endpoints. Models are fitted in a random 80% of participants and evaluated in the other 20%
## by Harrell's C; the difference in C is bootstrapped over the test set (1,000 replicates).
##  M0 age, sex, deprivation, smoking, splines of WC and BMI        M1 M0 + proWCdelta
##  M2 M0 + HbA1c, HDL and LDL cholesterol, log triglycerides, log CRP, systolic blood pressure, lipid-lowering,
##     antihypertensive and insulin treatment                        M3 M2 + proWCdelta
## All models use the same complete-case participants. Output: T40.
suppressPackageStartupMessages({ library(data.table); library(splines); library(survival) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
num <- function(x) suppressWarnings(as.numeric(x))
set.seed(20260925)
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "yu_ten_need_time", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "date_attending_assessment_centre", "Date_death_instance0", "HBA1C", "HDL", "LDLD", "TRIG", "CRP"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "t10", "dxmore", "dxall", "d0", "ddeath", "hba1c", "hdl", "ldl", "tg", "crp"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "SBP_auto_average", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "sbp", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- merge(merge(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, by = "id"), fr, by = "id", all.x = TRUE)
for (v in c("Age", "WC", "BMI", "hba1c", "hdl", "ldl", "tg", "crp", "sbp")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
s <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
d[, med := tolower(paste(s(med_m), s(med_f)))]
d[, `:=`(lipid_med = as.integer(grepl("cholesterol", med)), bp_med = as.integer(grepl("blood pressure", med)), insulin = as.integer(grepl("insulin", med)))]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt, hba1c, hdl, ldl, tg, crp, sbp)]) & tg > 0 & crp > 0]
d[, z := dlt / sd(dlt)]
d[, `:=`(dx10 = s(dx10), t10 = s(t10), dxmore = s(dxmore), dxall = s(dxall), b0 = as.IDate(d0), dd = as.IDate(ddeath))]
d[, tdth := (as.numeric(dd) - as.numeric(b0)) / 365.25]
d[, test := FALSE]; d[sample(.N, round(0.2 * .N)), test := TRUE]
cat("complete-case n =", nrow(d), "; test n =", sum(d$test), "\n")
first_time <- function(codes, times, pred) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]; k <- which(pred(cs)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(ts[k])), na.rm = TRUE) }
eps <- list(
  "Type 2 diabetes" = function(cs) startsWith(cs, "E11"),
  "Chronic kidney disease" = function(cs) startsWith(cs, "N18"),
  "Heart failure" = function(cs) startsWith(cs, "I50"),
  "Any circulatory diagnosis" = function(cs) startsWith(cs, "I") | startsWith(cs, "G45"))
m0 <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
m2 <- paste(m0, "+ ns(hba1c,3) + ns(hdl,3) + ns(ldl,3) + ns(log(tg),3) + ns(log(crp),3) + ns(sbp,3) + lipid_med + bp_med + insulin")
cidx <- function(fit, te) { lp <- predict(fit, newdata = te, type = "lp"); concordance(Surv(te$time, te$status) ~ lp, reverse = TRUE)$concordance }
res <- list()
run <- function(nm, x) {
  tr <- x[test == FALSE]; te <- x[test == TRUE]
  f <- lapply(list(M0 = m0, M1 = paste(m0, "+ z"), M2 = m2, M3 = paste(m2, "+ z")), function(r) coxph(as.formula(paste("Surv(time, status) ~", r)), tr))
  lp <- lapply(f, function(fit) predict(fit, newdata = te, type = "lp"))
  C <- function(idx, which) concordance(Surv(te$time[idx], te$status[idx]) ~ lp[[which]][idx], reverse = TRUE)$concordance
  all <- seq_len(nrow(te)); c0 <- C(all, "M0"); c1 <- C(all, "M1"); c2 <- C(all, "M2"); c3 <- C(all, "M3")
  bs <- replicate(1000, { b <- sample(all, replace = TRUE); c(C(b, "M1") - C(b, "M0"), C(b, "M3") - C(b, "M2")) })
  res[[length(res) + 1]] <<- data.table(endpoint = nm, n_train = nrow(tr), n_test = nrow(te), events_test = sum(te$status),
    C_M0 = c0, C_M1 = c1, dC_1v0 = c1 - c0, dC_1v0_lo = quantile(bs[1, ], 0.025), dC_1v0_hi = quantile(bs[1, ], 0.975),
    C_M2 = c2, C_M3 = c3, dC_3v2 = c3 - c2, dC_3v2_lo = quantile(bs[2, ], 0.025), dC_3v2_hi = quantile(bs[2, ], 0.975))
  cat(sprintf("%-32s test events %4d | C %.4f -> %.4f (dC %.4f, %.4f to %.4f) | with clinical %.4f -> %.4f (dC %.4f, %.4f to %.4f)\n", nm, sum(te$status), c0, c1, c1 - c0, quantile(bs[1, ], 0.025), quantile(bs[1, ], 0.975), c2, c3, c3 - c2, quantile(bs[2, ], 0.025), quantile(bs[2, ], 0.975))); flush.console()
}
for (nm in names(eps)) { pr <- eps[[nm]]
  inc <- vapply(strsplit(d$dx10, "[|]"), function(cs) any(pr(cs)), TRUE)
  incl <- inc | vapply(strsplit(d$dxmore, "[|]"), function(cs) any(pr(cs)), TRUE)
  prev <- vapply(strsplit(d$dxall, "[|]"), function(cs) any(pr(cs)), TRUE) & !incl
  x <- d[!prev]; x[, y := inc[!prev]]
  ft <- mapply(first_time, x$dx10, x$t10, MoreArgs = list(pred = pr))
  x[, tev := (ft - as.numeric(b0)) / 365.25]
  x[, time := pmin(fifelse(y, tev, 10), fifelse(is.na(tdth), 10, tdth), 10, na.rm = TRUE)]
  x[, status := as.integer(y & tev <= time + 1e-9)]; x <- x[time > 0]
  run(nm, x) }
x <- copy(d); x[, time := pmin(fifelse(is.na(tdth), 10, tdth), 10)]; x[, status := as.integer(!is.na(tdth) & tdth <= 10)]; x <- x[time > 0]
run("All-cause death", x)
fwrite(rbindlist(res), file.path(W, "tables", "T40_incremental_cox_cindex.csv"))
cat("DONE\n")
