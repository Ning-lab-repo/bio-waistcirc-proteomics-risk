## 68_biochem_nohba1c_splines_landmark.R
## NOTE (23 September 2026): section (3) took the time to an event from survival_time, the time to death or censoring,
## instead of the date of the diagnosis, so it did not exclude early events; it is superseded by
## 68b_five_year_landmark.R, which writes T105. Sections (1) and (2) are unaffected.
## (1) the routine-biochemistry comparison repeated with a score built without HbA1c, which is also one of the criteria
##     used to define prevalent diabetes and the strongest single predictor of a later diabetes code;
## (2) sensitivity of the primary estimate to the form of the adjustment for body size: natural splines with 3, 5 and 7
##     degrees of freedom, a tensor-style interaction between the WC and BMI splines, and stratification on a grid of
##     WC and BMI;
## (3) landmark analyses excluding events in the first 5 years, for heart failure, chronic kidney disease and liver
##     disease, with and without adjustment for the estimated glomerular filtration rate.
## Output: T103 (biochemistry without HbA1c), T104 (adjustment form), T105 (5-year landmark).
suppressPackageStartupMessages({ library(data.table); library(splines); library(caret); library(survival) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS", "HDL", "TRIG", "CRP"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys", "hdl", "tg", "crp"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
mt <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure10/WC_death_time.csv", select = c("V1", "all_cause_death", "survival_time")); setnames(mt, 1, "id")
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr, mt))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys", "hdl", "tg", "crp", "survival_time")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; d[, z := dlt / sd(dlt)]; say("modelling set n = %d", nrow(d))
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
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

## ---------------- (1) the biochemical score with and without HbA1c ----------------
set.seed(20260923); bd <- d[complete.cases(d[, .(hba1c, hdl, tg, crp)]) & tg > 0 & crp > 0]
folds <- createFolds(bd$WC, k = 10); bd[, `:=`(bWC = NA_real_, bWC3 = NA_real_)]
for (f in folds) {
  bd[f, bWC := predict(lm(WC ~ Age + Sex + ns(hba1c,3) + ns(hdl,3) + ns(log(tg),3) + ns(log(crp),3), bd[-f]), bd[f])]
  bd[f, bWC3 := predict(lm(WC ~ Age + Sex + ns(hdl,3) + ns(log(tg),3) + ns(log(crp),3), bd[-f]), bd[f])] }
bd[, `:=`(zb = (bWC - fitted(lm(bWC ~ WC, bd))), zb3 = (bWC3 - fitted(lm(bWC3 ~ WC, bd))))]
bd[, `:=`(zb = zb / sd(zb), zb3 = zb3 / sd(zb3), zp = dlt / sd(dlt))]
say("biochemical scores: R2 for WC with HbA1c %.3f, without %.3f; r(with, without) %.3f", summary(lm(WC ~ bWC, bd))$r.squared, summary(lm(WC ~ bWC3, bd))$r.squared, cor(bd$zb, bd$zb3))
o1 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- bd[get(paste0("free_", cd)) == TRUE]
  for (sp in list("zp", "zb", "zb3", c("zp", "zb3"))) { m <- glm(as.formula(paste(cd, "~", paste(sp, collapse = " + "), "+", BASE)), x, family = binomial())
    for (tt in sp) { r <- ci(m, tt); o1[[length(o1) + 1]] <- data.table(disease = dis$disease[i], model = paste(sp, collapse = "+"),
      term = c(zp = "proWCdelta", zb = "biochemical score (four tests)", zb3 = "biochemical score without HbA1c")[[tt]], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } } }
t103 <- rbindlist(o1); print(t103[, .(disease, model, term = substr(term, 1, 26), n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t103, file.path(TB, "T103_biochemical_score_without_hba1c.csv"))

## ---------------- (2) the form of the adjustment for body size ----------------
o2 <- list(); x <- d[free_E11 == TRUE]
forms <- list(c("natural splines, 3 df (primary)", "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"),
              c("natural splines, 5 df", "Age + Sex + tdi + smoking + ns(WC,5) + ns(BMI,5)"),
              c("natural splines, 7 df", "Age + Sex + tdi + smoking + ns(WC,7) + ns(BMI,7)"),
              c("splines with a WC-by-BMI interaction", "Age + Sex + tdi + smoking + ns(WC,3) * ns(BMI,3)"),
              c("linear WC and BMI", "Age + Sex + tdi + smoking + WC + BMI"))
for (f in forms) { r <- ci(glm(as.formula(paste("E11 ~ z +", f[2])), x, family = binomial()), "z")
  o2[[length(o2) + 1]] <- data.table(analysis = f[1], n = nrow(x), events = sum(x$E11), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
## stratification on a grid of WC and BMI (2 cm by 1 kg/m2 cells), pooled by conditional logistic regression
x[, cell := paste(Sex, floor(WC / 2), floor(BMI), sep = "_")]
kk <- x[, .(n = .N, ev = sum(E11)), by = cell][n >= 10 & ev >= 1]
cells <- kk[["cell"]]; xs <- x[cell %in% cells]; cm <- clogit(E11 ~ z + Age + strata(cell), data = xs)
sc <- summary(cm)$coefficients["z", ]
o2[[length(o2) + 1]] <- data.table(analysis = "stratified on sex and cells of 2 cm WC by 1 kg/m2 BMI", n = nrow(xs), events = sum(xs$E11), OR = exp(sc[[1]]), lo = exp(sc[[1]] - 1.96 * sc[[3]]), hi = exp(sc[[1]] + 1.96 * sc[[3]]))
t104 <- rbindlist(o2); print(t104[, .(analysis, n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t104, file.path(TB, "T104_adjustment_form_sensitivity.csv"))

## ---------------- (3) five-year landmark ----------------
o3 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]
  x <- d[get(paste0("free_", cd)) == TRUE]
  for (g in list(c("all events within 10 years", 0), c("excluding events in the first 5 years", 5))) {
    y <- copy(x); if (as.numeric(g[2]) > 0) { y <- y[!(get(cd) == 1 & survival_time < 5)]; y <- y[!(all_cause_death == 1 & survival_time < 5)] }
    r <- ci(glm(as.formula(paste(cd, "~ z +", BASE)), y, family = binomial()), "z")
    o3[[length(o3) + 1]] <- data.table(disease = dis$disease[i], analysis = g[1], n = nrow(y), events = sum(y[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
  ye <- x[!is.na(egfr)]; ye <- ye[!(get(cd) == 1 & survival_time < 5)][!(all_cause_death == 1 & survival_time < 5)]
  r <- ci(glm(as.formula(paste(cd, "~ z + ns(egfr,3) +", BASE)), ye, family = binomial()), "z")
  o3[[length(o3) + 1]] <- data.table(disease = dis$disease[i], analysis = "first 5 years excluded and adjusted for eGFR", n = nrow(ye), events = sum(ye[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t105 <- rbindlist(o3); print(t105[, .(disease, analysis = substr(analysis, 1, 40), n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t105, file.path(TB, "T105_five_year_landmark_superseded.csv")); say("DONE")
