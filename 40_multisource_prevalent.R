## 40_multisource_prevalent.R
## Stricter definition of baseline (prevalent) disease from several sources, then re-estimation of the main
## associations after excluding prevalent cases:
##   hospital record before baseline (all endpoints); self-reported illness at baseline; medication at baseline
##   (insulin or glucose-lowering drugs for type 2 diabetes, blood-pressure drugs for hypertension, lipid-lowering drugs for
##   dyslipidemia); HbA1c >= 48 mmol/mol (type 2 diabetes); eGFR < 60 (chronic kidney disease, CKD-EPI 2021 cr-cys).
## Estimates: continuous proWCdelta per SD (logistic and Cox, 10 years, covariates age, sex, deprivation, smoking,
## splines of WC and BMI) and the clinical-threshold comparison (normal WC + high proWC vs normal WC + normal proWC, M3).
## Output: T39.
suppressPackageStartupMessages({ library(data.table); library(splines); library(survival) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "yu_ten_need_time", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "date_attending_assessment_centre", "Date_death_instance0", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "t10", "dxmore", "dxall", "d0", "ddeath", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- merge(merge(merge(phen, lasso[, .(id, dlt = BioX_Delta, proWC = BioX_Adjusted)], by = "id"), cov, by = "id"), fr, by = "id", all.x = TRUE)
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]
d[, z := dlt / sd(dlt)]
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = s(dx10), t10 = s(t10), dxmore = s(dxmore), dxall = s(dxall), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname), b0 = as.IDate(d0), dd = as.IDate(ddeath))]
d[, dx10 := toupper(dx10)]; d[, dxmore := toupper(dxmore)]; d[, dxall := toupper(dxall)]; d[, t10 := t10]
## eGFR, CKD-EPI 2021 creatinine-cystatin C (Sex 0 = female)
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
glucose_drugs <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(
  E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, glucose_drugs) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I10 = has(d$selfrep, "(^|\\|)(hypertension|essential hypertension)(\\||$)") | has(d$med, "blood pressure medication"),
  E78 = has(d$selfrep, "(^|\\|)high cholesterol(\\||$)") | has(d$med, "cholesterol lowering medication"),
  I25 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"),
  K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60),
  E66 = rep(FALSE, nrow(d)))
dis <- data.table(disease = c("Type 2 diabetes", "Obesity", "Dyslipidemia", "Hypertension", "Ischaemic heart disease", "Heart failure", "Liver disease", "Chronic kidney disease"),
                  code = c("E11", "E66", "E78", "I10", "I25", "I50", "K76", "N18"))
first_time <- function(codes, times, code) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]; k <- which(startsWith(cs, code)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(ts[k])), na.rm = TRUE) }
adj <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
d[, thr := fifelse(Sex == 0, 88, 102)]; d[, hi := as.integer(proWC >= thr)]
out <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; nm <- dis$disease[i]; pat <- paste0("(^|[|])", cd)
  incl <- grepl(pat, d$dx10) | grepl(pat, d$dxmore); hosp_prev <- grepl(pat, d$dxall) & !incl
  prev <- hosp_prev | extra[[cd]]
  x <- d[!prev]; x[, y := as.integer(grepl(pat, dx10))]
  g <- summary(glm(as.formula(paste("y ~ z +", adj)), x, family = binomial()))$coefficients["z", ]
  ft <- mapply(first_time, x$dx10, x$t10, MoreArgs = list(code = cd))
  x[, tev := (ft - as.numeric(b0)) / 365.25]; x[, tdth := (as.numeric(dd) - as.numeric(b0)) / 365.25]
  x[, time := pmin(fifelse(y == 1, tev, 10), fifelse(is.na(tdth), 10, tdth), 10, na.rm = TRUE)]
  x[, status := fifelse(y == 1 & tev <= time + 1e-9, 1L, 0L)]; x <- x[time > 0]
  cx <- summary(coxph(as.formula(paste("Surv(time, status) ~ z +", adj)), x))$conf.int["z", ]
  th <- x[WC < thr]
  gt <- summary(glm(as.formula(paste("y ~ hi +", adj)), th, family = binomial()))$coefficients["hi", ]
  out[[i]] <- data.table(disease = nm, excluded_hospital = sum(hosp_prev), excluded_any_source = sum(prev), n = nrow(x), events = sum(x$y),
    OR_per_SD = exp(g[1]), OR_lo = exp(g[1] - 1.96 * g[2]), OR_hi = exp(g[1] + 1.96 * g[2]),
    HR_per_SD = cx[1], HR_lo = cx[3], HR_hi = cx[4],
    n_threshold = nrow(th), events_threshold = sum(th$y), OR_threshold = exp(gt[1]), ORt_lo = exp(gt[1] - 1.96 * gt[2]), ORt_hi = exp(gt[1] + 1.96 * gt[2]))
  cat(sprintf("%-24s excl %6d (hosp %5d) | n %d ev %d | OR/SD %.2f (%.2f-%.2f) | HR/SD %.2f (%.2f-%.2f) | threshold OR %.2f (%.2f-%.2f) [ev %d]\n", nm, sum(prev), sum(hosp_prev), nrow(x), sum(x$y), exp(g[1]), exp(g[1]-1.96*g[2]), exp(g[1]+1.96*g[2]), cx[1], cx[3], cx[4], exp(gt[1]), exp(gt[1]-1.96*gt[2]), exp(gt[1]+1.96*gt[2]), sum(th$y))); flush.console()
}
## the clinical-threshold columns (an earlier analysis, not reported) are printed above but not written
fwrite(rbindlist(out)[, !c("n_threshold", "events_threshold", "OR_threshold", "ORt_lo", "ORt_hi")], file.path(W, "tables", "T39_multisource_prevalent_exclusion.csv"))
cat("DONE\n")
