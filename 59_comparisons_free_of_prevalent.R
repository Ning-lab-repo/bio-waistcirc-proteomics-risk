## 59_comparisons_free_of_prevalent.R
## Re-estimation of the comparative analyses in participants free of each endpoint at baseline by all available sources
## (definitions as in 40_multisource_prevalent.R and 57_primary_incident_table2.R), in the formats of the original tables:
##   T4f progressive adjustment (two normal-WC groups), T6f matched pairs (pairs with both members free), T7f IDF thresholds,
##   T24f body fat and clinical markers (with separate "clinical markers" and "body fat + clinical markers" models and the
##   share of the log odds ratio removed), T32f routine-biochemistry comparison and top-tenth overlap, T65f training targets
##   (both adjustment sets, with likelihood-ratio statistics), T47f geographic hold-out.
suppressPackageStartupMessages({ library(data.table); library(splines); library(MatchIt); library(caret) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"; F2 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure2"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables"); O <- file.path(W, "output", "comparator_scores")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS", "HDL", "LDLD", "TRIG", "CRP"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys", "hdl", "ldl", "tg", "crp"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "SBP_auto_average", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "sbp", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bc <- fread(file.path(RAW, "bodycomp_subset.csv"), select = c("Participant ID", "Body fat percentage | Instance 0", "Trunk fat percentage | Instance 0")); setnames(bc, c("id", "bfp", "tfp"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, proWC = BioX_Adjusted, dlt = BioX_Delta)], by = "id"), cov, fr, bc))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "cre", "cys", "hdl", "ldl", "tg", "crp", "sbp", "bfp", "tfp")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; d[, z := dlt / sd(dlt)]; say("modelling set n = %d", nrow(d))
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, `:=`(lipid_med = as.integer(grepl("cholesterol", med) | grepl("(^|[| ])1([| ]|$)", med)), bp_med = as.integer(grepl("blood pressure", med) | grepl("(^|[| ])2([| ]|$)", med)), insulin = as.integer(grepl("insulin", med) | grepl("(^|[| ])3([| ]|$)", med)))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, gd) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I10 = has(d$selfrep, "(^|\\|)(hypertension|essential hypertension)(\\||$)") | has(d$med, "blood pressure medication"),
  E78 = has(d$selfrep, "(^|\\|)high cholesterol(\\||$)") | has(d$med, "cholesterol lowering medication"),
  I25 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"), I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"),
  K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"), N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60), E66 = rep(FALSE, nrow(d)))
dis <- data.table(disease = c("Type 2 diabetes", "Obesity", "Dyslipidemia", "Hypertension", "Ischaemic heart disease", "Heart failure", "Liver disease", "Chronic kidney disease"), code = c("E11", "E66", "E78", "I10", "I25", "I50", "K76", "N18"))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", cd); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
ci <- function(m, term) { s <- summary(m)$coefficients[term, ]; c(OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]]), p = s[[4]]) }
fr_of <- function(cd, x) x[get(paste0("free_", cd)) == TRUE]

## ---- A, C: threshold groups (clinical and IDF thresholds) ----
d[, thr := fifelse(Sex == 0, 88, 102)]; d[, thr2 := fifelse(Sex == 0, 80, 94)]
oA <- list(); oC <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- fr_of(cd, d)
  a <- x[WC < thr]; a[, NH := as.integer(proWC >= thr)]
  for (mm in list(c("M1 age,sex,TDI,smoking", "Age + Sex + tdi + smoking"), c("M2 + WC spline", "Age + Sex + tdi + smoking + ns(WC,3)"), c("M3 + WC spline + BMI spline", BASE))) {
    r <- ci(glm(as.formula(paste(cd, "~ NH +", mm[2])), a, family = binomial()), "NH"); oA[[length(oA) + 1]] <- data.table(disease = dis$disease[i], model = mm[1], n = nrow(a), events = sum(a[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
  b <- x[WC < thr2]; b[, NH := as.integer(proWC >= thr2)]
  for (mm in list(c("M1", "Age + Sex + tdi + smoking"), c("M3 + WC + BMI", BASE))) { r <- ci(glm(as.formula(paste(cd, "~ NH +", mm[2])), b, family = binomial()), "NH")
    oC[[length(oC) + 1]] <- data.table(disease = dis$disease[i], model = mm[1], n = nrow(b), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } }
fwrite(rbindlist(oA), file.path(TB, "T4f_progressive_adjustment_free.csv")); fwrite(rbindlist(oC), file.path(TB, "T7f_IDF_thresholds_free.csv"))
print(dcast(rbindlist(oA)[, .(disease, model, txt = sprintf("%.2f", OR))], disease ~ model, value.var = "txt"))

## ---- B: matched pairs (same matching as 05/08/50), pairs with both members free ----
m <- d[WC < thr]; m[, NH := as.integer(proWC >= thr)]; set.seed(123)
mt <- matchit(NH ~ Age + WC + BMI, data = as.data.frame(m), method = "nearest", exact = ~Sex, distance = "mahalanobis", caliper = c(WC = 0.15, BMI = 0.15), std.caliper = TRUE, ratio = 1)
md <- as.data.table(match.data(mt)); say("pairs %d", sum(md$NH == 1)); oB <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; fl <- md[, .(both = all(get(paste0("free_", cd)))), by = subclass][both == TRUE]$subclass; x <- md[subclass %in% fl]
  r <- ci(glm(as.formula(paste(cd, "~ NH + Age + Sex + tdi + smoking")), x, family = binomial()), "NH")
  oB[[i]] <- data.table(disease = dis$disease[i], pairs = length(fl), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
tB <- rbindlist(oB); print(tB); fwrite(tB, file.path(TB, "T6f_matched_pairs_free.csv"))

## ---- F: body fat and clinical markers ----
fat <- "ns(bfp,3) + ns(tfp,3)"; clin <- "ns(hba1c,3) + ns(hdl,3) + ns(ldl,3) + ns(log(tg),3) + ns(log(crp),3) + ns(sbp,3) + lipid_med + bp_med + insulin"
oF <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- fr_of(cd, d); x <- x[complete.cases(x[, .(bfp, tfp, hba1c, hdl, ldl, tg, crp, sbp)]) & tg > 0 & crp > 0]
  for (mm in list(c("WC and BMI", BASE), c("+ body fat", paste(BASE, "+", fat)), c("+ clinical markers", paste(BASE, "+", clin)), c("+ body fat + clinical markers", paste(BASE, "+", fat, "+", clin)))) {
    r <- ci(glm(as.formula(paste(cd, "~ z +", mm[2])), x, family = binomial()), "z"); oF[[length(oF) + 1]] <- data.table(disease = dis$disease[i], model = mm[1], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) } }
tF <- rbindlist(oF); fwrite(tF, file.path(TB, "T24f_fat_clinical_free.csv"))
sh <- dcast(tF, disease ~ model, value.var = "OR"); setnames(sh, c("disease", "fat", "fat_clin", "clin", "base")); sh[, `:=`(pct_clin = 100 * (1 - log(clin) / log(base)), pct_fat_clin = 100 * (1 - log(fat_clin) / log(base)))]; print(sh)
fwrite(sh, file.path(TB, "T71f_share_logOR_free.csv"))

## ---- G: routine-biochemistry comparison (as 36/50: same seed and folds) and top-tenth overlap ----
set.seed(20260923)
o <- fread(file.path(W, "output", "reduced_scores", "topk_oof_predictions.csv"))
bd <- merge(d, o[, .(id, dlt_k20)], by = "id"); bd <- bd[complete.cases(bd[, .(hba1c, hdl, tg, crp)]) & tg > 0 & crp > 0]
folds <- createFolds(bd$WC, k = 10); bd[, bWC := NA_real_]
for (f in folds) { mm <- lm(WC ~ Age + Sex + ns(hba1c,3) + ns(hdl,3) + ns(log(tg),3) + ns(log(crp),3), bd[-f]); bd[f, bWC := predict(mm, bd[f])] }
bd[, bdlt := bWC - fitted(lm(bWC ~ WC, bd))]; bd[, `:=`(z_pro = dlt / sd(dlt), z_20 = dlt_k20 / sd(dlt_k20), z_bio = bdlt / sd(bdlt))]
bd[, `:=`(topP = dlt >= quantile(dlt, 0.9), topB = bdlt >= quantile(bdlt, 0.9))]
bd[, grp := factor(fifelse(topP & topB, "both", fifelse(topP, "proteomic only", fifelse(topB, "biochemical only", "neither"))), levels = c("neither", "proteomic only", "biochemical only", "both"))]
oG <- list(); oT <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- fr_of(cd, bd)
  for (tt in list(c("Full proteomic score", "z_pro", "z_pro"), c("20-protein score", "z_20", "z_20"), c("Biochemical score", "z_bio", "z_bio"), c("Full proteomic score, adjusted for biochemical score", "z_pro + z_bio", "z_pro"), c("Biochemical score, adjusted for full proteomic score", "z_pro + z_bio", "z_bio"))) {
    r <- ci(glm(as.formula(paste(cd, "~", tt[2], "+", BASE)), x, family = binomial()), tt[3]); oG[[length(oG) + 1]] <- data.table(disease = dis$disease[i], series = tt[1], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
  mm <- summary(glm(as.formula(paste(cd, "~ grp +", BASE)), x, family = binomial()))$coefficients
  for (g in levels(x$grp)) { k <- paste0("grp", g); oT[[length(oT) + 1]] <- data.table(disease = dis$disease[i], group = g, n = sum(x$grp == g), events = sum(x[grp == g][[cd]]), risk = mean(x[grp == g][[cd]]),
    OR = if (g == "neither") 1 else exp(mm[k, 1]), lo = if (g == "neither") NA else exp(mm[k, 1] - 1.96 * mm[k, 2]), hi = if (g == "neither") NA else exp(mm[k, 1] + 1.96 * mm[k, 2])) } }
tG <- rbindlist(oG); fwrite(tG, file.path(TB, "T32f_biochemical_comparison_free.csv")); fwrite(rbindlist(oT), file.path(TB, "T60f_top_tenth_overlap_free.csv"))
print(dcast(tG[, .(disease, series = substr(series, 1, 22), txt = sprintf("%.2f", OR))], disease ~ series, value.var = "txt"), width = 250)
print(rbindlist(oT)[disease %in% c("Type 2 diabetes", "Heart failure", "Liver disease"), .(disease, group, risk = round(100 * risk, 1), OR = round(OR, 2))])

## ---- H: training targets, both adjustment sets ----
b <- fread(file.path(O, "oof_BMI.csv")); f2 <- fread(file.path(O, "oof_BFP.csv")); w <- fread(file.path(O, "oof_WHR.csv"))
e <- Reduce(function(a, c) merge(a, c, by = "id"), list(d, b[, .(id, dlt_bmi = delta)], f2[, .(id, dlt_bfp = delta)], w[, .(id, dlt_whr = delta)])); e <- e[complete.cases(e[, .(HC, bfp, dlt_bmi, dlt_bfp, dlt_whr)])]
for (v in c("dlt", "dlt_bmi", "dlt_bfp", "dlt_whr")) e[[paste0("z_", v)]] <- e[[v]] / sd(e[[v]])
SYM <- paste(BASE, "+ ns(bfp,3) + ns(HC,3)"); oH <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- fr_of(cd, e)
  for (adj in c("WC and BMI", "WC, BMI, body fat and hip")) { f0 <- if (adj == "WC and BMI") BASE else SYM; m0 <- glm(as.formula(paste(cd, "~", f0)), x, family = binomial())
    for (zz in c("z_dlt", "z_dlt_bmi", "z_dlt_bfp", "z_dlt_whr")) { mz <- glm(as.formula(paste(cd, "~", zz, "+", f0)), x, family = binomial()); r <- ci(mz, zz)
      oH[[length(oH) + 1]] <- data.table(disease = dis$disease[i], adjustment = adj, measure = c(z_dlt = "proWCdelta", z_dlt_bmi = "proBMIdelta", z_dlt_bfp = "proBFPdelta", z_dlt_whr = "proWHRdelta")[[zz]],
        n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], LR_chisq = 2 * (logLik(mz) - logLik(m0))[1]) } } }
tH <- rbindlist(oH); fwrite(tH, file.path(TB, "T65f_training_targets_free.csv"))
print(dcast(tH[, .(disease, adj = substr(adjustment, 1, 10), measure, txt = sprintf("%.2f", OR))], disease + adj ~ measure, value.var = "txt"), width = 200)

## ---- geographic hold-out ----
ho <- fread(file.path(F2, "lasso_test.csv")); setnames(ho, 1, "id"); h <- merge(d, ho[, .(id, dlt_ho = BioX_Delta)], by = "id"); h[, z_ho := dlt_ho / sd(dlt_ho)]; oO <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- fr_of(cd, h); if (sum(x[[cd]]) < 40) next
  r <- ci(glm(as.formula(paste(cd, "~ z_ho +", BASE)), x, family = binomial()), "z_ho"); oO[[length(oO) + 1]] <- data.table(disease = dis$disease[i], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
tO <- rbindlist(oO); print(tO); fwrite(tO, file.path(TB, "T47f_holdout_free.csv")); say("DONE")
