## 57_primary_incident_table2.R
## Primary analysis in participants free of each endpoint at baseline by all available sources (definitions as in
## 40_multisource_prevalent.R): per SD of proWCdelta at fixed WC and BMI, (a) per marginal SD, (b) per SD conditional on WC,
## BMI, age and sex, (c) regression-calibrated for error in single measurements of WC and BMI (first repeat assessment;
## bootstrap CI), and robustness of the primary estimate to (d) Olink plate, (e) excluding participants whose
## protein-predicted sex disagreed with recorded sex, (f) excluding participants with imputed Townsend index or smoking
## status; plus (g) standardised 10-year risks at the 10th and 90th conditional percentiles of proWCdelta.
## Output: T79 (Table 2 source), T80 (robustness of the primary estimate), T81 (conditional-percentile risks).
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
bt <- fread(bf, select = c(hdr[1], "PlateID", "TDI", "Smoking status")); setnames(bt, c("id", "plate", "tdi_raw", "smk_raw"))
bsz <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"), select = c("Participant ID", "Waist circumference | Instance 1", "Body mass index (BMI) | Instance 1")); setnames(bsz, c("id", "wc1", "bmi1"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id"), cov, fr, bt, bsz))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys", "wc1", "bmi1")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]; say("modelling set n = %d; plates %d", nrow(d), uniqueN(d$plate))
## conditional residual of proWCdelta given WC, BMI, age and sex
d[, e := resid(lm(dlt ~ ns(WC, 3) + ns(BMI, 3) + Age + Sex, d))]; sdc <- sd(d$e); d[, mu := dlt - e]; say("SD marginal %.2f, conditional %.2f cm", sdd, sdc)
## protein-predicted sex (as in 50_additional_analyses.R)
pp <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(pp, 1, "id"); prot <- setdiff(names(pp), c("id", "WC", "Age", "Sex")); ys <- as.integer(pp$Sex)
sdiff <- sapply(prot, function(v) { x <- pp[[v]]; abs(mean(x[ys == 1]) - mean(x[ys == 0])) / sd(x) }); top <- names(sort(sdiff, decreasing = TRUE))[1:30]
X <- as.data.frame(pp[, ..top]); X$ys <- ys; set.seed(7); fo <- sample(rep(1:5, length.out = nrow(pp))); ps <- numeric(nrow(pp))
for (k in 1:5) { f <- glm(ys ~ ., data = X[fo != k, ], family = binomial()); ps[fo == k] <- predict(f, X[fo == k, ], type = "response") }
disc <- pp$id[(ps >= 0.5) != (ys == 1)]; rm(pp, X); gc(); d[, sexdisc := id %in% disc]; say("sex-discordant in modelling set: %d", sum(d$sexdisc))
d[, imputed_cov := is.na(num(tdi_raw)) | is.na(num(smk_raw))]; say("imputed Townsend or smoking: %d", sum(d$imputed_cov))
## prevalent disease by all sources
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
glucose_drugs <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, glucose_drugs) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I10 = has(d$selfrep, "(^|\\|)(hypertension|essential hypertension)(\\||$)") | has(d$med, "blood pressure medication"),
  E78 = has(d$selfrep, "(^|\\|)high cholesterol(\\||$)") | has(d$med, "cholesterol lowering medication"),
  I25 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"), I2025 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60), E66 = rep(FALSE, nrow(d)))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease", "Ischaemic heart disease", "Ischaemic heart disease (I20-I25)", "Dyslipidemia", "Hypertension", "Obesity diagnosis"),
                  code = c("E11", "I50", "K76", "N18", "I25", "I2025", "E78", "I10", "E66"), pat = c("E11", "I50", "K76", "N18", "I25", "I2[0-5]", "E78", "I10", "E66"))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", dis$pat[i]); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
fz <- function(y, f, dd, term = "z") { s <- summary(glm(as.formula(paste(y, "~", term, "+", f)), dd, family = binomial()))$coefficients[term, ]; c(b = s[[1]], se = s[[2]]) }
## regression calibration model
cs <- d[!is.na(wc1) & !is.na(bmi1)]; cw <- lm(wc1 ~ WC + BMI + pWC + Age + Sex, cs); cb <- lm(bmi1 ~ WC + BMI + pWC + Age + Sex, cs)
d[, `:=`(WC_rc = predict(cw, d), BMI_rc = predict(cb, d))]; RCF <- "Age + Sex + tdi + smoking + ns(WC_rc,3) + ns(BMI_rc,3)"
out <- list(); rob <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[get(paste0("free_", cd)) == TRUE]
  all <- fz(cd, BASE, d); m <- fz(cd, BASE, x); rc <- fz(cd, RCF, x)
  rcb <- unlist(mclapply(1:200, function(b) { set.seed(b); ii <- sample.int(nrow(cs), replace = TRUE); c1 <- lm(wc1 ~ WC + BMI + pWC + Age + Sex, cs[ii]); c2 <- lm(bmi1 ~ WC + BMI + pWC + Age + Sex, cs[ii])
    xx <- x[sample.int(nrow(x), replace = TRUE)]; xx[, `:=`(WC_rc = predict(c1, xx), BMI_rc = predict(c2, xx))]; coef(glm(as.formula(paste(cd, "~ z +", RCF)), xx, family = binomial()))[["z"]] }, mc.cores = 20))
  kc <- sdc / sdd
  out[[i]] <- data.table(disease = dis$disease[i], n_all = nrow(d), events_all = sum(d[[cd]]), OR_all = exp(all[["b"]]), lo_all = exp(all[["b"]] - 1.96 * all[["se"]]), hi_all = exp(all[["b"]] + 1.96 * all[["se"]]),
    n_free = nrow(x), events_free = sum(x[[cd]]), OR_free = exp(m[["b"]]), lo_free = exp(m[["b"]] - 1.96 * m[["se"]]), hi_free = exp(m[["b"]] + 1.96 * m[["se"]]),
    OR_free_conditional_SD = exp(m[["b"]] * kc), lo_cond = exp((m[["b"]] - 1.96 * m[["se"]]) * kc), hi_cond = exp((m[["b"]] + 1.96 * m[["se"]]) * kc),
    OR_free_rc = exp(rc[["b"]]), lo_rc = exp(quantile(rcb, 0.025)), hi_rc = exp(quantile(rcb, 0.975)))
  say("%-34s all %.2f | free n=%d ev=%d OR %.2f (%.2f-%.2f) | cond %.2f | rc %.2f (%.2f-%.2f)", dis$disease[i], exp(all[["b"]]), nrow(x), sum(x[[cd]]), exp(m[["b"]]), exp(m[["b"]] - 1.96 * m[["se"]]), exp(m[["b"]] + 1.96 * m[["se"]]), exp(m[["b"]] * kc), exp(rc[["b"]]), exp(quantile(rcb, 0.025)), exp(quantile(rcb, 0.975)))
  if (cd %in% c("E11", "I50", "N18", "K76")) {
    pl <- fz(cd, paste(BASE, "+ factor(plate)"), x); sx <- fz(cd, BASE, x[sexdisc == FALSE]); ic <- fz(cd, BASE, x[imputed_cov == FALSE])
    rob[[cd]] <- data.table(disease = dis$disease[i], analysis = c("primary", "plus Olink plate (fixed effects)", "excluding protein-predicted sex discordance", "excluding imputed Townsend index or smoking"),
      n = c(nrow(x), nrow(x), sum(!x$sexdisc), sum(!x$imputed_cov)), OR = exp(c(m[["b"]], pl[["b"]], sx[["b"]], ic[["b"]])),
      lo = exp(c(m[["b"]], pl[["b"]], sx[["b"]], ic[["b"]]) - 1.96 * c(m[["se"]], pl[["se"]], sx[["se"]], ic[["se"]])), hi = exp(c(m[["b"]], pl[["b"]], sx[["b"]], ic[["b"]]) + 1.96 * c(m[["se"]], pl[["se"]], sx[["se"]], ic[["se"]])))
    print(rob[[cd]]) } }
t79 <- rbindlist(out); fwrite(t79, file.path(TB, "T79_table2_primary_incident.csv")); t80 <- rbindlist(rob); fwrite(t80, file.path(TB, "T80_primary_robustness_plate_sex_imputation.csv"))
## (g) standardised risks at conditional percentiles (type 2 diabetes, heart failure; free of disease)
qc <- quantile(d$e, c(0.1, 0.9)); o <- list()
for (cd in c("E11", "I50", "N18", "K76")) { x <- d[get(paste0("free_", cd)) == TRUE]; mm <- glm(as.formula(paste(cd, "~ z +", BASE)), x, family = binomial())
  r <- sapply(qc, function(q) { nd <- copy(x); nd[, z := (mu + q) / sdd]; mean(predict(mm, nd, type = "response")) })
  o[[cd]] <- data.table(disease = cd, conditional_p10_cm = qc[[1]], conditional_p90_cm = qc[[2]], risk_p10 = r[[1]], risk_p90 = r[[2]], n = nrow(x))
  say("%s: risk at conditional 10th/90th percentile (%.1f / %.1f cm from expected): %.2f%% / %.2f%%", cd, qc[[1]], qc[[2]], 100 * r[[1]], 100 * r[[2]]) }
fwrite(rbindlist(o), file.path(TB, "T81_conditional_percentile_risks.csv")); say("DONE")
