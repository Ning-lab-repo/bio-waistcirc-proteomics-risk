## 47_further_robustness.R
## Further robustness analyses:
##  (A) matched pairs: conditional logistic and cluster-robust estimates; the same after excluding prevalent disease defined
##      from five sources (hospital records, self-report, medication, HbA1c >= 48 mmol/mol, eGFR < 60), as in 40_multisource_prevalent.R
##  (B) incremental value by 5-fold cross-validation repeated twice (pooled out-of-fold linear predictors, bootstrap CI,
##      calibration slope), hospital-record definition and, for type 2 diabetes and CKD, the five-source definition
##  (C) type 2 diabetes in normoglycaemic participants (HbA1c < 42 mmol/mol, no diabetes by any source); fasting-time strata
##  (D) adjustment for prior contact with hospital care (number of distinct hospital diagnoses before baseline) and
##      number of self-reported non-cancer illnesses, including for injury and fracture
##  (E) sex-stratified per-SD estimates and interaction tests
##  (F) threshold comparison excluding participants within 2 cm below the WC threshold
##  (G) DXA measures adjusted for WC, BMI and height measured at the imaging visit
##  (H) bootstrap 95% CI for the out-of-fold and hold-out R2; R2 by protein missingness
## Output: T50-T57.
suppressPackageStartupMessages({ library(data.table); library(splines); library(survival); library(MatchIt); library(sandwich) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; F2 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure2"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
BSZ <- file.path(RAW, "6-1Body_size_measures_participant.csv")
num <- function(x) suppressWarnings(as.numeric(x))
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
s <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }

## ---------------- data ----------------
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "yu_ten_need_time",
  "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "data_Diagnoses_ICD10", "date_attending_assessment_centre", "Date_death_instance0",
  "Noncancer_illne_code_elfreported_Instance0", "Number of self-reported non-cancer illnesses", "HBA1C", "HDL", "LDLD", "TRIG", "CRP", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "t10", "dxmore", "dxall", "tall", "d0", "ddeath", "selfrep", "n_selfrep", "hba1c", "hdl", "ldl", "tg", "crp", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "SBP_auto_average", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "sbp", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
fh <- names(fread(file.path(F8, "wc_pro_Batch.csv"), nrows = 0)); fast <- fread(file.path(F8, "wc_pro_Batch.csv"), select = c(fh[1], "Fasting_time")); setnames(fast, c("id", "fasting"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, dlt = BioX_Delta, proWC = BioX_Adjusted)], by = "id"), cov, fr, fast))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "hdl", "ldl", "tg", "crp", "cre", "cys", "sbp", "fasting", "n_selfrep")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]
sdd <- sd(d$dlt); d[, z := dlt / sdd]
say("analysis n = %d", nrow(d))
d[, `:=`(dx10 = toupper(s(dx10)), t10 = s(t10), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), tall = s(tall), selfrep = tolower(s(selfrep)),
         med = tolower(paste(s(med_m), s(med_f))), medname = tolower(s(medname)), b0 = as.IDate(d0), dd = as.IDate(ddeath))]
d[, tdth := (as.numeric(dd) - as.numeric(b0)) / 365.25]
d[, `:=`(lipid_med = as.integer(grepl("cholesterol", med)), bp_med = as.integer(grepl("blood pressure", med)), insulin = as.integer(grepl("insulin", med)))]
## eGFR (CKD-EPI 2021, creatinine-cystatin C; Sex 0 = female)
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
## number of distinct hospital diagnoses (3-character codes) recorded before baseline
prior_n <- function(codes, dates, b0) { if (codes == "") return(0L); cs <- strsplit(codes, "[|]")[[1]]; ds <- as.IDate(strsplit(dates, "[|]")[[1]]); k <- !is.na(ds) & ds < b0; length(unique(substr(cs[k], 1, 3))) }
d[, n_prior := mapply(prior_n, dxall, tall, b0)]
say("prior hospital diagnoses: median %.0f, IQR %.0f-%.0f, zero %.1f%%", median(d$n_prior), quantile(d$n_prior, .25), quantile(d$n_prior, .75), 100 * mean(d$n_prior == 0))
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
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pat <- paste0("(^|[|])", cd)
  inc <- grepl(pat, d$dx10); incl <- inc | grepl(pat, d$dxmore); hosp_prev <- grepl(pat, d$dxall) & !incl
  d[[cd]] <- as.integer(inc); d[[paste0("prevH_", cd)]] <- hosp_prev; d[[paste0("prevM_", cd)]] <- hosp_prev | extra[[cd]] }
d[, INJ := as.integer(grepl("(^|[|])[ST]", dx10))]; d[, FRAC := as.integer(grepl("(^|[|])(S[0-9]2|T0[28]|T1[02])", dx10))]
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
fz <- function(y, f, dat, term = "z") { m <- glm(as.formula(paste(y, "~", term, "+", f)), dat, family = binomial()); co <- summary(m)$coefficients[term, ]
  c(OR = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) }

## ---------------- (A) matched pairs ----------------
d[, thr := fifelse(Sex == 0, 88, 102)]
d[, grp := fifelse(WC < thr & proWC < thr, "N/N", fifelse(WC < thr & proWC >= thr, "N/H", fifelse(WC >= thr & proWC < thr, "H/N", "H/H")))]
m <- d[grp %in% c("N/N", "N/H")]; m[, NH := as.integer(grp == "N/H")]
set.seed(123)
mt <- matchit(NH ~ Age + WC + BMI, data = as.data.frame(m), method = "nearest", exact = ~Sex, distance = "mahalanobis", caliper = c(WC = 0.15, BMI = 0.15), std.caliper = TRUE, ratio = 1)
md <- as.data.table(match.data(mt)); say("matched pairs %d", sum(md$NH == 1))
oA <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]
  for (def in c("hospital-record definition (main analysis)", "five-source definition of prevalent disease")) {
    x <- if (startsWith(def, "hospital")) copy(md) else md[get(paste0("prevM_", cd)) == FALSE]
    g <- glm(as.formula(paste(cd, "~ NH + Age + Sex + tdi + smoking")), x, family = binomial())
    b <- coef(g)[["NH"]]; se_n <- sqrt(vcov(g)["NH", "NH"]); se_c <- sqrt(vcovCL(g, cluster = ~subclass)["NH", "NH"])
    pairs <- x[, .N, by = subclass][N == 2, subclass]; xc <- x[subclass %in% pairs]
    cl <- tryCatch(summary(clogit(as.formula(paste(cd, "~ NH + strata(subclass)")), xc))$coefficients["NH", ], error = function(e) rep(NA, 5))
    oA[[length(oA) + 1]] <- data.table(disease = dis$disease[i], definition = def, n = nrow(x), events = sum(x[[cd]]), intact_pairs = length(pairs),
      OR = exp(b), lo_model = exp(b - 1.96 * se_n), hi_model = exp(b + 1.96 * se_n), lo_cluster = exp(b - 1.96 * se_c), hi_cluster = exp(b + 1.96 * se_c),
      OR_conditional = exp(cl[1]), lo_conditional = exp(cl[1] - 1.96 * cl[3]), hi_conditional = exp(cl[1] + 1.96 * cl[3])) } }
tA <- rbindlist(oA); print(tA[, .(disease, def = substr(definition, 1, 8), events, OR = round(OR, 2), clus = sprintf("%.2f-%.2f", lo_cluster, hi_cluster), cond = sprintf("%.2f (%.2f-%.2f)", OR_conditional, lo_conditional, hi_conditional))])
fwrite(tA, file.path(TB, "T50_matched_pairs_robust_and_multisource.csv"))

## ---------------- (C) normoglycaemic participants; fasting strata ----------------
ng <- d[!is.na(hba1c) & hba1c < 42 & prevM_E11 == FALSE]
oC <- list(); r <- fz("E11", BASE, ng)
oC[[1]] <- data.table(analysis = "HbA1c < 42 mmol/mol and no diabetes by any source", n = nrow(ng), events = sum(ng$E11), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
r <- fz("E11", paste(BASE, "+ ns(hba1c,3)"), ng)
oC[[2]] <- data.table(analysis = "same, additionally adjusted for HbA1c", n = nrow(ng), events = sum(ng$E11), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
for (fs in c("fasting < 4 h", "fasting >= 4 h")) { x <- if (fs == "fasting < 4 h") d[!is.na(fasting) & fasting < 4] else d[!is.na(fasting) & fasting >= 4]; r <- fz("E11", BASE, x)
  oC[[length(oC) + 1]] <- data.table(analysis = paste("all participants,", fs), n = nrow(x), events = sum(x$E11), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
xf <- d[!is.na(fasting)]; xf[, f4 := as.integer(fasting >= 4)]
pint <- summary(glm(as.formula(paste("E11 ~ z * f4 +", BASE)), xf, family = binomial()))$coefficients["z:f4", 4]
tC <- rbindlist(oC); tC[, p_interaction_fasting := pint]; print(tC)
fwrite(tC, file.path(TB, "T51_normoglycaemic_and_fasting.csv"))

## ---------------- (D) prior hospital contact and self-reported illness ----------------
x <- d[!is.na(n_selfrep)]; x[, `:=`(lprior = log1p(n_prior), nsr = factor(pmin(n_selfrep, 3)))]
oD <- list()
for (y in c(dis$code, "INJ", "FRAC")) { lab <- if (y == "INJ") "Injury and poisoning (S00-T98)" else if (y == "FRAC") "Fracture" else dis[code == y]$disease
  r0 <- fz(y, BASE, x); r1 <- fz(y, paste(BASE, "+ ns(lprior,3) + nsr"), x)
  oD[[length(oD) + 1]] <- data.table(outcome = lab, n = nrow(x), events = sum(x[[y]]), OR_base = r0[["OR"]], lo_base = r0[["lo"]], hi_base = r0[["hi"]], OR_adjusted = r1[["OR"]], lo_adjusted = r1[["lo"]], hi_adjusted = r1[["hi"]]) }
tD <- rbindlist(oD); print(tD[, .(outcome, base = sprintf("%.2f", OR_base), adj = sprintf("%.2f (%.2f-%.2f)", OR_adjusted, lo_adjusted, hi_adjusted))])
fwrite(tD, file.path(TB, "T52_prior_hospital_contact_adjustment.csv"))

## ---------------- (E) sex-stratified ----------------
oE <- list(); BASE_NS <- "Age + tdi + smoking + ns(WC,3) + ns(BMI,3)"
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]
  rm_ <- fz(cd, BASE_NS, d[Sex == 1]); rw <- fz(cd, BASE_NS, d[Sex == 0])
  pint <- summary(glm(as.formula(paste(cd, "~ z * Sex +", BASE)), d, family = binomial()))$coefficients["z:Sex", 4]
  oE[[i]] <- data.table(disease = dis$disease[i], OR_men = rm_[["OR"]], lo_men = rm_[["lo"]], hi_men = rm_[["hi"]], events_men = sum(d[Sex == 1][[cd]]),
                        OR_women = rw[["OR"]], lo_women = rw[["lo"]], hi_women = rw[["hi"]], events_women = sum(d[Sex == 0][[cd]]), p_interaction = pint) }
tE <- rbindlist(oE); print(tE[, .(disease, men = sprintf("%.2f (%.2f-%.2f)", OR_men, lo_men, hi_men), women = sprintf("%.2f (%.2f-%.2f)", OR_women, lo_women, hi_women), p_int = signif(p_interaction, 2))])
fwrite(tE, file.path(TB, "T53_sex_stratified_continuous.csv"))

## ---------------- (F) threshold comparison excluding WC within 2 cm below the threshold ----------------
oF <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; nw <- d[grp %in% c("N/N", "N/H")]; nw[, NH := as.integer(grp == "N/H")]
  for (lab in c("all normal-WC participants", "excluding WC within 2 cm below the threshold")) { x <- if (startsWith(lab, "all")) nw else nw[WC < thr - 2]
    r <- fz(cd, BASE, x, term = "NH")
    oF[[length(oF) + 1]] <- data.table(disease = dis$disease[i], sample = lab, n = nrow(x), n_NH = sum(x$NH), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } }
tF <- rbindlist(oF); print(dcast(tF[, .(disease, sample = substr(sample, 1, 12), txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease ~ sample, value.var = "txt"))
fwrite(tF, file.path(TB, "T54_threshold_excluding_near_threshold.csv"))

## ---------------- (G) DXA adjusted for anthropometry at the imaging visit ----------------
bs <- fread(BSZ, select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2", "Height | Instance 2", "Standing height | Instance 0"))
setnames(bs, c("id", "wc2", "bmi2", "ht2", "ht0"))
bc <- fread(file.path(RAW, "bodycomp_subset.csv"), select = c("Participant ID", "VAT (visceral adipose tissue) mass | Instance 2", "Total fat mass | Instance 2", "Trunk fat mass | Instance 2", "Android fat mass | Instance 2", "Total lean mass | Instance 2"))
setnames(bc, c("id", "vat", "tfm", "trunk", "android", "lean"))
g <- merge(merge(d[, .(id, Age, Sex, WC, BMI, dlt)], bs, by = "id"), bc, by = "id")
for (v in c("wc2", "bmi2", "ht2", "ht0", "vat", "tfm", "trunk", "android", "lean")) g[[v]] <- num(g[[v]])
oG <- list()
for (v in c("tfm", "trunk", "android", "vat", "lean")) { lab <- c(tfm = "Total fat mass", trunk = "Trunk fat mass", android = "Android fat mass", vat = "Visceral fat mass (DXA estimate)", lean = "Total lean mass")[[v]]
  for (adj in c("baseline WC, BMI and height", "imaging-visit WC, BMI and height")) {
    x <- g[!is.na(get(v)) & !is.na(wc2) & !is.na(bmi2) & !is.na(ht2) & !is.na(ht0)]; x[, yz := as.numeric(scale(get(v)))]
    f <- if (startsWith(adj, "baseline")) "ns(WC,3) + ns(BMI,3) + Age + factor(Sex) + ht0" else "ns(wc2,3) + ns(bmi2,3) + Age + factor(Sex) + ht2"
    fit <- lm(as.formula(paste("yz ~ dlt +", f)), x); co <- summary(fit)$coefficients["dlt", ]
    oG[[length(oG) + 1]] <- data.table(measure = lab, adjustment = adj, n = nrow(x), beta_per_SD = co[[1]] * sdd, lo = (co[[1]] - 1.96 * co[[2]]) * sdd, hi = (co[[1]] + 1.96 * co[[2]]) * sdd, p = co[[4]]) } }
tG <- rbindlist(oG); print(dcast(tG[, .(measure, adj = substr(adjustment, 1, 8), n, txt = sprintf("%.3f (%.3f to %.3f)", beta_per_SD, lo, hi))], measure + n ~ adj, value.var = "txt"))
fwrite(tG, file.path(TB, "T55_DXA_imaging_visit_adjustment.csv"))

## ---------------- (H) R2 with bootstrap CI; R2 by protein missingness ----------------
set.seed(2026)
r2c <- function(a, b) cor(a, b)^2
ho <- fread(file.path(F2, "lasso_test.csv")); setnames(ho, 1, "id")
bt_ci <- function(a, b, B = 1000) { n <- length(a); v <- replicate(B, { k <- sample.int(n, n, TRUE); r2c(a[k], b[k]) }); quantile(v, c(.025, .975)) }
full <- lasso[, .(a = num(Actual_WC), b = num(pred_WC))][complete.cases(a, b)]
oH <- list(data.table(set = "out-of-fold, full analysis set", n = nrow(full), R2 = r2c(full$a, full$b), lo = bt_ci(full$a, full$b)[[1]], hi = bt_ci(full$a, full$b)[[2]]),
           data.table(set = "geographic hold-out (England-trained model)", n = nrow(ho), R2 = r2c(ho$Actual_WC, ho$pred_WC), lo = bt_ci(ho$Actual_WC, ho$pred_WC)[[1]], hi = bt_ci(ho$Actual_WC, ho$pred_WC)[[2]]))
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1
pm <- fread(bf, select = c(hdr[1], hdr[2:pe])); setnames(pm, 1, "id"); pm <- data.table(id = pm$id, pna = rowMeans(is.na(as.matrix(pm[, -1]))))
full2 <- merge(lasso[, .(id, a = num(Actual_WC), b = num(pred_WC))], pm, by = "id")
for (lab in c("<= 20% of proteins missing", "> 20% of proteins missing")) { x <- if (startsWith(lab, "<=")) full2[pna <= .2] else full2[pna > .2]
  ci <- bt_ci(x$a, x$b); oH[[length(oH) + 1]] <- data.table(set = paste("out-of-fold,", lab), n = nrow(x), R2 = r2c(x$a, x$b), lo = ci[[1]], hi = ci[[2]]) }
tH <- rbindlist(oH); print(tH)
fwrite(tH, file.path(TB, "T56_R2_bootstrap_and_missingness.csv"))

## ---------------- (B) cross-validated incremental value ----------------
cc <- d[complete.cases(d[, .(hba1c, hdl, ldl, tg, crp, sbp)]) & tg > 0 & crp > 0]
say("incremental-value complete-case n = %d", nrow(cc))
first_time <- function(codes, times, pred) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]; k <- which(pred(cs)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(ts[k])), na.rm = TRUE) }
m0 <- BASE; m2 <- paste(m0, "+ ns(hba1c,3) + ns(hdl,3) + ns(ldl,3) + ns(log(tg),3) + ns(log(crp),3) + ns(sbp,3) + lipid_med + bp_med + insulin")
eps <- list("Type 2 diabetes" = function(cs) startsWith(cs, "E11"), "Chronic kidney disease" = function(cs) startsWith(cs, "N18"),
            "Heart failure" = function(cs) startsWith(cs, "I50"), "Any circulatory diagnosis" = function(cs) startsWith(cs, "I") | startsWith(cs, "G45"))
mk <- function(pr, drop_ms = NULL) { inc <- vapply(strsplit(cc$dx10, "[|]"), function(cs) any(pr(cs)), TRUE)
  incl <- inc | vapply(strsplit(cc$dxmore, "[|]"), function(cs) any(pr(cs)), TRUE); prev <- vapply(strsplit(cc$dxall, "[|]"), function(cs) any(pr(cs)), TRUE) & !incl
  if (!is.null(drop_ms)) prev <- prev | cc[[drop_ms]]
  x <- cc[!prev]; x[, y := inc[!prev]]; ft <- mapply(first_time, x$dx10, x$t10, MoreArgs = list(pred = pr)); x[, tev := (ft - as.numeric(b0)) / 365.25]
  x[, time := pmin(fifelse(y, tev, 10), fifelse(is.na(tdth), 10, tdth), 10, na.rm = TRUE)]; x[, status := as.integer(y & tev <= time + 1e-9)]; x[time > 0] }
cvrun <- function(nm, x, reps = 2, K = 5, B = 500) {
  lps <- list()
  for (r in seq_len(reps)) { set.seed(1000 + r); fold <- sample(rep_len(1:K, nrow(x)))
    lp <- matrix(NA_real_, nrow(x), 4, dimnames = list(NULL, c("M0", "M1", "M2", "M3")))
    for (k in 1:K) { tr <- x[fold != k]; te <- x[fold == k]
      for (mn in c("M0", "M1", "M2", "M3")) { f <- switch(mn, M0 = m0, M1 = paste(m0, "+ z"), M2 = m2, M3 = paste(m2, "+ z"))
        lp[fold == k, mn] <- predict(coxph(as.formula(paste("Surv(time, status) ~", f)), tr), newdata = te, type = "lp") } }
    lps[[r]] <- lp }
  C <- function(lp, idx) concordance(Surv(x$time[idx], x$status[idx]) ~ lp[idx], reverse = TRUE)$concordance
  all <- seq_len(nrow(x)); cm <- sapply(c("M0", "M1", "M2", "M3"), function(mn) mean(sapply(lps, function(lp) C(lp[, mn], all))))
  lp1 <- lps[[1]]; bsv <- replicate(B, { b <- sample(all, replace = TRUE); c(C(lp1[, "M1"], b) - C(lp1[, "M0"], b), C(lp1[, "M3"], b) - C(lp1[, "M2"], b)) })
  cal <- coef(coxph(Surv(x$time, x$status) ~ lp1[, "M1"]))[[1]]
  res <- data.table(endpoint = nm, n = nrow(x), events = sum(x$status), C_M0 = cm[["M0"]], C_M1 = cm[["M1"]], dC_1v0 = cm[["M1"]] - cm[["M0"]], dC_1v0_lo = quantile(bsv[1, ], .025), dC_1v0_hi = quantile(bsv[1, ], .975),
                    C_M2 = cm[["M2"]], C_M3 = cm[["M3"]], dC_3v2 = cm[["M3"]] - cm[["M2"]], dC_3v2_lo = quantile(bsv[2, ], .025), dC_3v2_hi = quantile(bsv[2, ], .975), calibration_slope_M1 = cal)
  say("%-45s ev %5d | C %.4f->%.4f dC %.4f (%.4f-%.4f) | clin %.4f->%.4f dC %.4f (%.4f-%.4f) | cal %.2f", nm, sum(x$status), cm[["M0"]], cm[["M1"]], res$dC_1v0, res$dC_1v0_lo, res$dC_1v0_hi, cm[["M2"]], cm[["M3"]], res$dC_3v2, res$dC_3v2_lo, res$dC_3v2_hi, cal)
  res }
oB <- list()
for (nm in names(eps)) oB[[length(oB) + 1]] <- cvrun(nm, mk(eps[[nm]]))
xd <- copy(cc); xd[, time := pmin(fifelse(is.na(tdth), 10, tdth), 10)]; xd[, status := as.integer(!is.na(tdth) & tdth <= 10)]; xd <- xd[time > 0]
oB[[length(oB) + 1]] <- cvrun("All-cause death", xd)
oB[[length(oB) + 1]] <- cvrun("Type 2 diabetes, five-source prevalent exclusion", mk(eps[["Type 2 diabetes"]], "prevM_E11"))
oB[[length(oB) + 1]] <- cvrun("Chronic kidney disease, five-source prevalent exclusion", mk(eps[["Chronic kidney disease"]], "prevM_N18"))
fwrite(rbindlist(oB), file.path(TB, "T57_incremental_cox_crossvalidated.csv"))
say("DONE")
