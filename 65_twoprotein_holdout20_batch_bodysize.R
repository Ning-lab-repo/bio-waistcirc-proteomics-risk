## 65_twoprotein_holdout20_batch_bodysize.R
## (1) a waist score built in the same way from age, sex, GDF15 and adrenomedullin only, and its discordance measure
##     compared with proWCdelta for heart failure;
## (2) the 20-protein score applied to the Scotland and Wales hold-out, both with the deposited coefficients (fitted in
##     the full analysis set) and with coefficients refitted in England only;
## (3) the processing batch that lacks the block of 1,461 proteins: characteristics and a batch-by-proWCdelta interaction;
## (4) primary estimates with height, weight and hip circumference also held fixed;
## (5) the contribution of each protein to the variance of pWC (coefficient times the protein's standard deviation).
## Output: T93 (two-protein score), T94 (20-protein hold-out), T95 (batch), T96 (body size), T97 (protein contributions).
suppressPackageStartupMessages({ library(data.table); library(splines); library(caret) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
F2 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure2"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; D <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
mf <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(mf, 1, "id"); prot <- setdiff(names(mf), c("id", "WC", "Age", "Sex"))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bsz <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"), select = c("Participant ID", "Standing height | Instance 0", "Weight | Instance 0")); setnames(bsz, c("id", "height", "weight"))
bf <- file.path(F8, "wc_pro_Batch.csv"); tech <- fread(bf, select = c(names(fread(bf, nrows = 0))[1], "Batch", "sample_age_days")); setnames(tech, c("id", "batch", "storage"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr, bsz, tech))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "cre", "cys", "height", "weight", "storage")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
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
ci <- function(m, term) { s <- summary(m)$coefficients[term, ]; c(OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]]), p = s[[4]]) }

## ---------------- (1) a two-protein (GDF15 and adrenomedullin) waist score ----------------
pp <- intersect(c("GDF15", "ADM"), prot); stopifnot(length(pp) == 2)
tp <- merge(d, mf[, c("id", ..pp)], by = "id"); set.seed(20260923)
folds <- createFolds(tp$WC, k = 10); tp[, pw2 := NA_real_]
for (f in folds) { m <- lm(as.formula(paste("WC ~ Age + Sex +", paste(pp, collapse = " + "))), tp[-f]); tp[f, pw2 := predict(m, tp[f])] }
tp[, d2 := pw2 - fitted(lm(pw2 ~ WC, tp))]; tp[, z2 := d2 / sd(d2)]
say("two-protein score: R2 for WC %.3f; r(discordance, proWCdelta) %.3f", summary(lm(WC ~ pw2, tp))$r.squared, cor(tp$d2, tp$dlt))
o1 <- list(); for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- tp[get(paste0("free_", cd)) == TRUE]
  for (sp in list("z", "z2", c("z", "z2"))) { m <- glm(as.formula(paste(cd, "~", paste(sp, collapse = " + "), "+", BASE)), x, family = binomial())
    for (tt in sp) { r <- ci(m, tt); o1[[length(o1) + 1]] <- data.table(disease = dis$disease[i], model = paste(sp, collapse = "+"), term = c(z = "proWCdelta", z2 = "GDF15 and adrenomedullin score")[[tt]], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } } }
t93 <- rbindlist(o1); print(t93[, .(disease, model, term = substr(term, 1, 24), n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t93, file.path(TB, "T93_two_protein_stress_score.csv"))

## ---------------- (2) the 20-protein score in the geographic hold-out ----------------
ho <- fread(file.path(F2, "lasso_test.csv")); setnames(ho, 1, "id"); hoid <- ho$id
c20 <- fread(file.path(D, "05_source_data", "model_reduced_score_20_proteins_coefficients.csv")); setnames(c20, 1:2, c("term", "coef"))
p20 <- setdiff(c20$term, c("(Intercept)", "Age", "Sex")); stopifnot(all(p20 %in% prot))
m20 <- merge(d[, .(id, Age, Sex, WC, BMI, tdi, smoking, E11, I50, K76, N18, free_E11, free_I50, free_K76, free_N18)], mf[, c("id", ..p20)], by = "id")
b <- setNames(c20$coef, c20$term)
m20[, pw20 := b[["(Intercept)"]] + b[["Age"]] * Age + b[["Sex"]] * Sex + as.matrix(.SD) %*% b[p20], .SDcols = p20]
tr <- m20[!(id %in% hoid)]; te <- m20[id %in% hoid]; say("hold-out n = %d; training n = %d", nrow(te), nrow(tr))
fit_en <- lm(as.formula(paste("WC ~ Age + Sex +", paste(p20, collapse = " + "))), tr)
te[, pw20_en := predict(fit_en, te)]
o2 <- list()
for (v in c("pw20", "pw20_en")) { te[, dd := get(v) - fitted(lm(get(v) ~ WC, te))]; te[, zz := dd / sd(dd)]
  r2 <- summary(lm(as.formula(paste("WC ~", v)), te))$r.squared
  for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- te[get(paste0("free_", cd)) == TRUE]
    r <- ci(glm(as.formula(paste(cd, "~ zz +", BASE)), x, family = binomial()), "zz")
    o2[[length(o2) + 1]] <- data.table(score = c(pw20 = "20-protein score, deposited coefficients (fitted in the full set)", pw20_en = "20-protein score, coefficients refitted in England only")[[v]],
      disease = dis$disease[i], n = nrow(x), events = sum(x[[cd]]), R2_for_WC = r2, OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } }
t94 <- rbindlist(o2); print(t94[, .(score = substr(score, 1, 40), disease, n, events, R2 = round(R2_for_WC, 3), txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t94, file.path(TB, "T94_reduced_score_holdout.csv"))

## ---------------- (3) the batch that lacks the protein block ----------------
d[, b7 := as.integer(batch == 7)]
t95a <- d[!is.na(b7), .(n = .N, age = mean(Age), men_pct = 100 * mean(Sex == 1), WC = mean(WC), BMI = mean(BMI), storage_days = mean(storage, na.rm = TRUE), proWCdelta = mean(dlt), T2D_pct = 100 * mean(E11)), by = .(batch_7 = b7)]
print(t95a)
x <- d[free_E11 == TRUE & !is.na(b7)]
m <- glm(as.formula(paste("E11 ~ z * b7 +", BASE)), x, family = binomial()); ii <- summary(m)$coefficients["z:b7", ]
say("batch-7 interaction with proWCdelta for type 2 diabetes: ratio of odds ratios %.2f (%.2f-%.2f), P = %.3f", exp(ii[[1]]), exp(ii[[1]] - 1.96 * ii[[2]]), exp(ii[[1]] + 1.96 * ii[[2]]), ii[[4]])
## the pooled interaction term above is printed only; the text reports the estimates within and outside batch 7
## (source_table_Results_missing_protein_block_strata.csv)
t95 <- t95a[, .(quantity = paste0("batch 7 = ", batch_7), n, age, men_pct, WC, BMI, storage_days, proWCdelta, T2D_pct)]
fwrite(t95, file.path(TB, "T95_batch7_characteristics.csv"))

## ---------------- (4) height, weight and hip circumference also held fixed ----------------
o4 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[get(paste0("free_", cd)) == TRUE & !is.na(height) & !is.na(weight) & !is.na(HC)]
  for (mm in list(c("WC and BMI", BASE), c("+ height, weight and hip", paste(BASE, "+ ns(height,3) + ns(weight,3) + ns(HC,3)")))) {
    r <- ci(glm(as.formula(paste(cd, "~ z +", mm[2])), x, family = binomial()), "z")
    o4[[length(o4) + 1]] <- data.table(disease = dis$disease[i], model = mm[1], n = nrow(x), events = sum(x[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) } }
t96 <- rbindlist(o4); print(t96[, .(disease, model, n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
fwrite(t96, file.path(TB, "T96_height_weight_hip.csv"))

## ---------------- (5) contribution of each protein to the variance of pWC ----------------
full <- fread(file.path(D, "05_source_data", "model_proWC_full_coefficients.csv")); setnames(full, 1:2, c("term", "coef"))
fp <- full[!term %in% c("(Intercept)", "Age", "Sex")]
sds <- mf[, lapply(.SD, function(v) sd(num(v), na.rm = TRUE)), .SDcols = fp$term]
fp[, sd_protein := as.numeric(sds[1])]; fp[, contribution := abs(coef) * sd_protein]
fp <- fp[order(-contribution)]; tot <- sum(fp$contribution)
fp[, share_percent := 100 * contribution / tot]
say("top proteins by contribution to the spread of pWC (coefficient x SD): %s", paste(sprintf("%s %.2f cm (%.1f%%)", fp$term[1:8], fp$contribution[1:8], fp$share_percent[1:8]), collapse = "; "))
say("the 20 largest contributions account for %.0f%% of the total", sum(fp$share_percent[1:20]))
fwrite(fp[, .(protein = term, coefficient = coef, sd_protein, contribution_cm = contribution, share_percent)], file.path(TB, "T97_protein_contributions.csv"))
say("DONE")
