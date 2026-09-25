## 68b_five_year_landmark.R
## Five-year landmark for the primary endpoints, in participants free of each endpoint at baseline by all five sources
## (as in 68_biochem_nohba1c_splines_landmark.R): participants with a first diagnosis of the endpoint, or who died, in
## the first 5 years after baseline are excluded, and diagnoses in years 5 to 10 are modelled by logistic regression with
## the primary adjustment, with and without a natural spline of the estimated glomerular filtration rate.
## This replaces section (3) of script 68, which took the time to an event from survival_time, the time to death or
## censoring, instead of the date of the first diagnosis. The time to diagnosis is taken, as in
## 39_continuous_sensitivity.R, from the dates recorded with each ICD-10 code (yu_ten_need_time).
## Output: T105_five_year_landmark.csv (replaces the earlier version).
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "yu_ten_need_time", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS", "date_attending_assessment_centre", "Date_death_instance0"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "t10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys", "d0", "ddeath"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; d[, z := dlt / sd(dlt)]; say("modelling set n = %d", nrow(d))
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(t10 = as.character(t10), dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[is.na(t10), t10 := ""]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, gd) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease"), code = c("E11", "I50", "K76", "N18"))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", cd); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
ci <- function(m, term) { s <- summary(m)$coefficients[term, ]; c(OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]])) }
## date of the first diagnosis with the given code among the diagnoses recorded within 10 years (39_continuous_sensitivity.R)
first_time <- function(codes, times, code) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]
  k <- which(startsWith(cs, code)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(ts[k])), na.rm = TRUE) }
d[, `:=`(b0 = as.numeric(as.IDate(d0)), tdth = (as.numeric(as.IDate(ddeath)) - as.numeric(as.IDate(d0))) / 365.25)]

o <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]
  x <- d[get(paste0("free_", cd)) == TRUE]
  x[, tev := (mapply(first_time, x$dx10, x$t10, MoreArgs = list(code = cd)) - b0) / 365.25]
  nmiss <- sum(x[[cd]] == 1 & is.na(x$tev))
  early <- (x[[cd]] == 1 & !is.na(x$tev) & x$tev < 5) | (!is.na(x$tdth) & x$tdth < 5)
  say("%s: %d free of the endpoint, %d events, %d events in the first 5 years, %d deaths in the first 5 years without an earlier event, %d events without a date",
      dis$disease[i], nrow(x), sum(x[[cd]]), sum(x[[cd]] == 1 & !is.na(x$tev) & x$tev < 5),
      sum(!is.na(x$tdth) & x$tdth < 5 & !(x[[cd]] == 1 & !is.na(x$tev) & x$tev < 5)), nmiss)
  r <- ci(glm(as.formula(paste(cd, "~ z +", BASE)), x, family = binomial()), "z")
  o[[length(o) + 1]] <- data.table(disease = dis$disease[i], analysis = "all events within 10 years", n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
  y <- x[!early]
  r <- ci(glm(as.formula(paste(cd, "~ z +", BASE)), y, family = binomial()), "z")
  o[[length(o) + 1]] <- data.table(disease = dis$disease[i], analysis = "event or death in the first 5 years excluded", n = nrow(y), events = sum(y[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
  ye <- y[!is.na(egfr)]
  r <- ci(glm(as.formula(paste(cd, "~ z + ns(egfr,3) +", BASE)), ye, family = binomial()), "z")
  o[[length(o) + 1]] <- data.table(disease = dis$disease[i], analysis = "first 5 years excluded and adjusted for eGFR", n = nrow(ye), events = sum(ye[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t105 <- rbindlist(o); print(t105[, .(disease, analysis = substr(analysis, 1, 45), n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t105, file.path(TB, "T105_five_year_landmark.csv")); say("DONE")
