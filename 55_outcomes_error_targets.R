## 55_outcomes_error_targets.R
## (1) absolute risks in participants free of the endpoint at baseline by all available sources (hospital record before
##     baseline, self-report, medication, HbA1c >= 48 mmol/mol for type 2 diabetes, eGFR < 60 for chronic kidney disease;
##     definitions as in 40_multisource_prevalent.R): standardised 10-year risks by proWCdelta at fixed WC and BMI,
##     matched-pair risks in pairs whose two members are both free of the endpoint, and the Table rows for type 2 diabetes
## (2) why the odds ratio rises from M1 to M3: model with the WC spline only, and within-sex correlations of proWCdelta
##     with BMI, height and hip circumference
## (3) regression calibration for error in single baseline measurements of WC and BMI, using the first repeat assessment
## (4) training targets: bootstrap CIs for differences in log odds ratios between proWCdelta and the BMI, body-fat and
##     waist-to-hip measures, and odds ratios per SD conditional on WC and BMI
## Output: T73-T77.
suppressPackageStartupMessages({ library(data.table); library(splines); library(MatchIt); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
TB <- file.path(W, "tables"); O <- file.path(W, "output", "comparator_scores")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bsz <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"), select = c("Participant ID", "Waist circumference | Instance 1", "Body mass index (BMI) | Instance 1", "Standing height | Instance 0"))
setnames(bsz, c("id", "wc1", "bmi1", "height"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, proWC = BioX_Adjusted, dlt = BioX_Delta)], by = "id"), cov, fr, bsz))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "cre", "cys", "wc1", "bmi1", "height")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]; say("modelling set n = %d", nrow(d))
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
glucose_drugs <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, glucose_drugs) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Chronic kidney disease", "Liver disease"), code = c("E11", "I50", "N18", "K76"))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pat <- paste0("(^|[|])", cd); inc <- grepl(pat, d$dx10) | grepl(pat, d$dxmore)
  d[[cd]] <- as.integer(grepl(pat, d$dx10)); d[[paste0("prev_", cd)]] <- (grepl(pat, d$dxall) & !inc) | extra[[cd]] }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"

## (1a) standardised absolute risks in participants free of the endpoint
q <- quantile(d$z, c(0.1, 0.5, 0.9)); pts <- c(p10 = q[[1]], p50 = q[[2]], p90 = q[[3]], `+2 SD` = 2)
stdrisk <- function(dd, cd) { m <- glm(as.formula(paste(cd, "~ z +", BASE)), dd, family = binomial())
  r <- sapply(pts, function(v) { nd <- copy(dd); nd[, z := v]; mean(predict(m, nd, type = "response")) }); nd1 <- copy(dd); nd1[, z := z + 1]
  c(r, rd_per_SD = mean(predict(m, nd1, type = "response")) - mean(predict(m, dd, type = "response")), OR = exp(coef(m)[["z"]])) }
set.seed(2026); o1 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[!get(paste0("prev_", cd))]; est <- stdrisk(x, cd)
  bs <- do.call(rbind, mclapply(1:200, function(b) { set.seed(b); stdrisk(x[sample.int(nrow(x), replace = TRUE)], cd) }, mc.cores = 20))
  o1[[i]] <- data.table(disease = dis$disease[i], quantity = names(est), estimate = est, lo = apply(bs, 2, quantile, 0.025), hi = apply(bs, 2, quantile, 0.975), n = nrow(x), events = sum(x[[cd]]))
  say("%s free at baseline n=%d events=%d: OR %.2f; risk p10 %.1f%%, p50 %.1f%%, p90 %.1f%%; RD/SD %.2f", dis$disease[i], nrow(x), sum(x[[cd]]), est[["OR"]], 100 * est[["p10"]], 100 * est[["p50"]], 100 * est[["p90"]], 100 * est[["rd_per_SD"]]) }
t73 <- rbindlist(o1); fwrite(t73, file.path(TB, "T73_standardised_risk_free_of_prevalent_disease.csv"))

## (1b) matched pairs (same matching as 08/50) restricted to pairs with both members free of the endpoint
d[, thr := fifelse(Sex == 0, 88, 102)]
m <- d[WC < thr]; m[, NH := as.integer(proWC >= thr)]; set.seed(123)
mt <- matchit(NH ~ Age + WC + BMI, data = as.data.frame(m), method = "nearest", exact = ~Sex, distance = "mahalanobis", caliper = c(WC = 0.15, BMI = 0.15), std.caliper = TRUE, ratio = 1)
md <- as.data.table(match.data(mt)); say("pairs %d", sum(md$NH == 1))
o2 <- list(); for (cd in dis$code) { w <- dcast(md[, .(subclass, NH, y = get(cd), p = get(paste0("prev_", cd)))], subclass ~ NH, value.var = c("y", "p")); w <- w[!p_0 & !p_1]
  dif <- w$y_1 - w$y_0; se <- sd(dif) / sqrt(nrow(w))
  o2[[cd]] <- data.table(disease = dis[code == cd]$disease, pairs = nrow(w), risk_discordant = mean(w$y_1), risk_concordant = mean(w$y_0), difference = mean(dif), lo = mean(dif) - 1.96 * se, hi = mean(dif) + 1.96 * se,
                         events_discordant = sum(w$y_1), events_concordant = sum(w$y_0))
  say("matched %s: pairs %d, risk %.1f%% vs %.1f%%", cd, nrow(w), 100 * mean(w$y_1), 100 * mean(w$y_0)) }
t74 <- rbindlist(o2); fwrite(t74, file.path(TB, "T74_matched_pairs_risk_free_of_prevalent_disease.csv"))

## (1c) Table rows for type 2 diabetes among participants free of diabetes by all sources (all 52,879-equivalent modelling set)
x <- d[prev_E11 == FALSE]; x[, grp := factor(fifelse(WC < thr & proWC < thr, "N/N", fifelse(WC < thr, "N/H", fifelse(proWC < thr, "H/N", "H/H"))), levels = c("N/N", "N/H", "H/N", "H/H"))]
fq <- x[, .(n = .N, events = sum(E11), pct = 100 * mean(E11)), by = grp][order(grp)]; print(fq)
g1 <- summary(glm(E11 ~ grp + Age + Sex + tdi + smoking, x, family = binomial()))$coefficients
xn <- x[grp %in% c("N/N", "N/H")]; xn[, hi := as.integer(grp == "N/H")]; g3 <- summary(glm(as.formula(paste("E11 ~ hi +", BASE)), xn, family = binomial()))$coefficients["hi", ]
f <- function(b, se) sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
t75 <- data.table(row = c("event frequency", "OR age, sex, deprivation, smoking", "OR plus splines of WC and BMI"),
  `N/N` = c(sprintf("%.1f%%", fq$pct[1]), "Reference", "Reference"), `N/H` = c(sprintf("%.1f%%", fq$pct[2]), f(g1["grpN/H", 1], g1["grpN/H", 2]), f(g3[[1]], g3[[2]])),
  `H/N` = c(sprintf("%.1f%%", fq$pct[3]), f(g1["grpH/N", 1], g1["grpH/N", 2]), "-"), `H/H` = c(sprintf("%.1f%%", fq$pct[4]), f(g1["grpH/H", 1], g1["grpH/H", 2]), "-"), n = nrow(x))
print(t75); fwrite(rbind(t75, data.table(row = "n by group", `N/N` = fq$n[1], `N/H` = fq$n[2], `H/N` = fq$n[3], `H/H` = fq$n[4], n = nrow(x)), fill = TRUE), file.path(TB, "T75_table_T2D_free_of_prevalent_diabetes.csv"))

## (2) from M1 to M3
fz <- function(f, dd = d, y = "E11") { s <- summary(glm(as.formula(paste(y, "~ z +", f)), dd, family = binomial()))$coefficients["z", ]; exp(c(s[[1]], s[[1]] - 1.96 * s[[2]], s[[1]] + 1.96 * s[[2]])) }
m1 <- fz("Age + Sex + tdi + smoking"); m2 <- fz("Age + Sex + tdi + smoking + ns(WC,3)"); m2b <- fz("Age + Sex + tdi + smoking + ns(BMI,3)"); m3 <- fz(BASE)
cr <- rbindlist(lapply(0:1, function(sx) { e <- d[Sex == sx]; data.table(sex = c("women", "men")[sx + 1], r_BMI = cor(e$dlt, e$BMI), r_WC = cor(e$dlt, e$WC), r_height = cor(e$dlt, e$height, use = "complete.obs"), r_hip = cor(e$dlt, e$HC, use = "complete.obs"),
  partial_r_BMI_at_fixed_WC = cor(resid(lm(dlt ~ ns(WC, 3) + Age, e)), resid(lm(BMI ~ ns(WC, 3) + Age, e)))) })); print(cr)
say("T2D OR per SD: M1 %.2f; + WC spline %.2f; + BMI spline only %.2f; + both %.2f", m1[1], m2[1], m2b[1], m3[1])

## (3) regression calibration with the first repeat measurement
cs <- d[!is.na(wc1) & !is.na(bmi1)]; say("calibration sample (first repeat assessment) n = %d", nrow(cs))
cw <- lm(wc1 ~ WC + BMI + pWC + Age + Sex, cs); cb <- lm(bmi1 ~ WC + BMI + pWC + Age + Sex, cs); print(summary(cw)$coefficients); print(summary(cb)$coefficients)
d[, `:=`(WC_rc = predict(cw, d), BMI_rc = predict(cb, d))]
rc <- fz("Age + Sex + tdi + smoking + ns(WC_rc,3) + ns(BMI_rc,3)")
set.seed(7); rcb <- unlist(mclapply(1:200, function(b) { set.seed(b); i <- sample.int(nrow(cs), replace = TRUE); c1 <- lm(wc1 ~ WC + BMI + pWC + Age + Sex, cs[i]); c2 <- lm(bmi1 ~ WC + BMI + pWC + Age + Sex, cs[i])
  j <- sample.int(nrow(d), replace = TRUE); dd <- d[j]; dd[, `:=`(WC_rc = predict(c1, dd), BMI_rc = predict(c2, dd))]
  coef(glm(E11 ~ z + Age + Sex + tdi + smoking + ns(WC_rc, 3) + ns(BMI_rc, 3), dd, family = binomial()))[["z"]] }, mc.cores = 20))
say("regression-calibrated T2D OR per SD %.2f (bootstrap 95%% CI %.2f-%.2f); naive %.2f", rc[1], exp(quantile(rcb, 0.025)), exp(quantile(rcb, 0.975)), m3[1])
rc_hf <- fz("Age + Sex + tdi + smoking + ns(WC_rc,3) + ns(BMI_rc,3)", y = "I50"); rc_ckd <- fz("Age + Sex + tdi + smoking + ns(WC_rc,3) + ns(BMI_rc,3)", y = "N18")
t76 <- data.table(quantity = c("M1 age, sex, deprivation, smoking", "M2 + WC spline", "M2b + BMI spline only", "M3 + WC and BMI splines", "regression-calibrated WC and BMI (T2D)", "regression-calibrated WC and BMI (heart failure)", "regression-calibrated WC and BMI (chronic kidney disease)"),
  OR = c(m1[1], m2[1], m2b[1], m3[1], rc[1], rc_hf[1], rc_ckd[1]), lo = c(m1[2], m2[2], m2b[2], m3[2], exp(quantile(rcb, 0.025)), rc_hf[2], rc_ckd[2]), hi = c(m1[3], m2[3], m2b[3], m3[3], exp(quantile(rcb, 0.975)), rc_hf[3], rc_ckd[3]),
  note = c(rep("type 2 diabetes", 4), "CI from 200 bootstrap samples of calibration and main data", "Wald CI", "Wald CI"))
t76 <- rbind(t76, data.table(quantity = c("calibration sample n", "coefficient of pWC for repeat WC", "coefficient of measured WC for repeat WC", "coefficient of pWC for repeat BMI"),
  OR = c(nrow(cs), coef(cw)[["pWC"]], coef(cw)[["WC"]], coef(cb)[["pWC"]]), lo = NA, hi = NA, note = "calibration model"), fill = TRUE)
t76 <- rbind(t76, cr[, .(quantity = paste("within-sex correlations,", sex), OR = r_BMI, lo = r_height, hi = r_hip, note = sprintf("OR column = r with BMI; lo = r with height; hi = r with hip; partial r with BMI at fixed WC %.3f", partial_r_BMI_at_fixed_WC))])
fwrite(t76, file.path(TB, "T76_M1_to_M3_and_regression_calibration.csv"))

## (4) training targets: bootstrap differences in log OR and ORs per conditional SD
b <- fread(file.path(O, "oof_BMI.csv")); f2 <- fread(file.path(O, "oof_BFP.csv")); w <- fread(file.path(O, "oof_WHR.csv"))
bc <- fread(file.path(RAW, "bodycomp_subset.csv"), select = c("Participant ID", "Body fat percentage | Instance 0")); setnames(bc, c("id", "bfp"))
e <- Reduce(function(a, c) merge(a, c, by = "id"), list(d[, .(id, Age, Sex, WC, BMI, HC, tdi, smoking, dlt, E11, I50, N18)], b[, .(id, dlt_bmi = delta)], f2[, .(id, dlt_bfp = delta)], w[, .(id, dlt_whr = delta)], bc))
e[, bfp := num(bfp)]; e <- e[complete.cases(e[, .(HC, bfp, dlt_bmi, dlt_bfp, dlt_whr)])]; say("training-target sample n = %d", nrow(e))
ms <- c(WC = "dlt", BMI = "dlt_bmi", BFP = "dlt_bfp", WHR = "dlt_whr"); ADJ <- c(fixed_WC_BMI = BASE, fixed_WC_BMI_fat_hip = paste(BASE, "+ ns(bfp,3) + ns(HC,3)"))
for (k in names(ms)) { v <- ms[[k]]; e[[paste0("z_", k)]] <- e[[v]] / sd(e[[v]]); e[[paste0("c_", k)]] <- e[[v]] / sd(resid(lm(as.formula(paste(v, "~ ns(WC,3) + ns(BMI,3) + Age + Sex")), e))) }
lor <- function(dd) unlist(lapply(c("E11", "I50", "N18"), function(cd) unlist(lapply(names(ADJ), function(a) sapply(names(ms), function(k) coef(glm(as.formula(paste(cd, "~", paste0("z_", k), "+", ADJ[[a]])), dd, family = binomial()))[[2]])))))
est <- lor(e); nm <- as.vector(outer(names(ms), outer(names(ADJ), c("E11", "I50", "N18"), paste, sep = "|"), paste, sep = "|"))
bs <- do.call(rbind, mclapply(1:200, function(bb) { set.seed(1000 + bb); lor(e[sample.int(nrow(e), replace = TRUE)]) }, mc.cores = 20))
o4 <- list(); for (cd in c("E11", "I50", "N18")) for (a in names(ADJ)) { ix <- which(grepl(paste0("\\|", a, "\\|", cd, "$"), nm)); kk <- sub("\\|.*", "", nm[ix]); wc <- ix[kk == "WC"]
  for (j in ix[kk != "WC"]) { dlt <- est[wc] - est[j]; bd <- bs[, wc] - bs[, j]
    o4[[length(o4) + 1]] <- data.table(endpoint = cd, adjustment = a, comparison = paste("WC minus", sub("\\|.*", "", nm[j])), OR_WC = exp(est[wc]), OR_other = exp(est[j]), ratio_of_ORs = exp(dlt), lo = exp(quantile(bd, 0.025)), hi = exp(quantile(bd, 0.975))) } }
t77 <- rbindlist(o4); print(t77[, .(endpoint, adjustment, comparison, OR_WC = round(OR_WC, 2), OR_other = round(OR_other, 2), ratio = round(ratio_of_ORs, 2), lo = round(lo, 2), hi = round(hi, 2))])
cond <- rbindlist(lapply(c("E11", "I50", "N18"), function(cd) rbindlist(lapply(names(ms), function(k) { sm <- summary(glm(as.formula(paste(cd, "~", paste0("c_", k), "+", BASE)), e, family = binomial()))$coefficients[2, ]
  data.table(endpoint = cd, measure = k, OR_per_conditional_SD = exp(sm[[1]]), lo = exp(sm[[1]] - 1.96 * sm[[2]]), hi = exp(sm[[1]] + 1.96 * sm[[2]])) }))))
print(cond); fwrite(rbind(t77, cond, fill = TRUE), file.path(TB, "T77_training_target_differences_bootstrap.csv")); say("DONE")
