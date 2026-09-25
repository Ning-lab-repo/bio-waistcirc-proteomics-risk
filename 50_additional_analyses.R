## 50_additional_analyses.R
## Additional analyses:
##  (1) standard ICD-10 code sets: ischaemic heart disease I20-I25, liver disease K70/K74-K76, type 2 diabetes E11 or E14
##  (2) proWCdelta additionally adjusted for PRAP1, RTN4R, IGSF3 and IGSF9, singly and together, and for all 20 proteins of
##      the reduced score together
##  (3) overlap with the biochemical discordance measure (rebuilt exactly as in 36_biochemical_score_head_to_head.R):
##      risk in the top decile of each measure alone and of both
##  (4) matched pairs: absolute 10-year risks and risk differences; paired differences in bioimpedance body composition
##  (5) record-linkage check: sex predicted from plasma proteins against recorded sex
##  (6) baseline characteristics overall, by clinical-threshold group and by region (Supplementary Table 4), with
##      medians (IQR) for C-reactive protein and triglycerides
##  (7) clinical-threshold groups with anchoring estimated separately in men and women
## Output: T58-T64.
suppressPackageStartupMessages({ library(data.table); library(splines); library(MatchIt); library(caret); library(glmnet) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); s <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "s_Diagnoses_ICD10", "data_Diagnoses_ICD10", "date_attending_assessment_centre", "HBA1C", "HDL", "LDLD", "TRIG", "CRP", "UK Biobank assessment centre | Instance 0"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "dxall", "tall", "d0", "hba1c", "hdl", "ldl", "tg", "crp", "centre"))
cv0 <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cv0, 1, "id")
cov <- cv0[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
fh <- names(fread(file.path(F8, "wc_pro_Batch.csv"), nrows = 0))
eth <- fread(file.path(F8, "wc_pro_Batch.csv"), select = c(fh[1], "ethnicity", "Smoking status")); setnames(eth, c("id", "white", "smk_raw"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "med_m", "med_f")); fr[, med := tolower(paste(s(med_m), s(med_f)))]
fr[, `:=`(lipid_med = as.integer(grepl("cholesterol", med) | grepl("(^|[| ])1([| ]|$)", med)), bp_med = as.integer(grepl("blood pressure", med) | grepl("(^|[| ])2([| ]|$)", med)), insulin = as.integer(grepl("insulin", med) | grepl("(^|[| ])3([| ]|$)", med)))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, proWC = BioX_Adjusted, dlt = BioX_Delta)], by = "id"), cov, eth, fr[, .(id, lipid_med, bp_med, insulin)]))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "hdl", "ldl", "tg", "crp")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d[, `:=`(dx10 = toupper(s(dx10)), dxall = toupper(s(dxall)), tall = s(tall), b0 = as.IDate(d0))]
all0 <- copy(d)                                   # 52,879 for descriptive tables
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]
say("modelling set n = %d", nrow(d))
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
fz <- function(y, f, dat, term = "z") { m <- glm(as.formula(paste(y, "~", term, "+", f)), dat, family = binomial()); co <- summary(m)$coefficients[term, ]
  c(OR = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]])) }
hasc <- function(x, pat) as.integer(grepl(pat, x, perl = TRUE))
d[, `:=`(E11 = hasc(dx10, "(^|\\|)E11"), I50 = hasc(dx10, "(^|\\|)I50"), N18 = hasc(dx10, "(^|\\|)N18"), I25 = hasc(dx10, "(^|\\|)I25"), K76 = hasc(dx10, "(^|\\|)K76"),
         E1114 = hasc(dx10, "(^|\\|)E1[14]"), I2025 = hasc(dx10, "(^|\\|)I2[0-5]"), LIV = hasc(dx10, "(^|\\|)K7[0456]"))]

## (1) standard code sets
o1 <- list(); for (pr in list(c("E11", "Type 2 diabetes, E11"), c("E1114", "Type 2 diabetes, E11 or E14"), c("I25", "Ischaemic heart disease, I25"), c("I2025", "Ischaemic heart disease, I20-I25"), c("K76", "Liver disease, K76"), c("LIV", "Liver disease, K70 or K74-K76"))) {
  r <- fz(pr[1], BASE, d); o1[[length(o1) + 1]] <- data.table(outcome = pr[2], n = nrow(d), events = sum(d[[pr[1]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t1 <- rbindlist(o1); print(t1); fwrite(t1, file.path(TB, "T58_standard_code_sets.csv"))

## (2) protein adjustment
pd <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(pd, 1, "id")
p20 <- c("LEP", "FABP4", "NCAN", "IL1RN", "WFIKKN2", "PON3", "SLITRK1", "IGFBP2", "IGSF9", "ADM", "IGFBP1", "SSC4D", "CHGB", "FURIN", "RTN4R", "CFH", "IGSF3", "OPTC", "SEZ6L", "PRAP1")
b <- merge(d, pd[, c("id", p20), with = FALSE], by = "id"); for (v in p20) b[[paste0("s_", v)]] <- as.numeric(scale(b[[v]]))
adj <- list("none" = "", "PRAP1" = "s_PRAP1", "RTN4R" = "s_RTN4R", "IGSF3" = "s_IGSF3", "IGSF9" = "s_IGSF9", "PRAP1, RTN4R, IGSF3 and IGSF9" = "s_PRAP1 + s_RTN4R + s_IGSF3 + s_IGSF9",
            "all 20 proteins of the reduced score" = paste(paste0("s_", p20), collapse = " + "))
o2 <- list(); for (cd in c("E11", "I50", "N18")) for (an in names(adj)) { f <- if (adj[[an]] == "") BASE else paste(BASE, "+", adj[[an]]); r <- fz(cd, f, b)
  o2[[length(o2) + 1]] <- data.table(disease = c(E11 = "Type 2 diabetes", I50 = "Heart failure", N18 = "Chronic kidney disease")[[cd]], adjustment = an, n = nrow(b), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t2 <- rbindlist(o2); print(t2[, .(disease, adjustment, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))]); fwrite(t2, file.path(TB, "T59_protein_adjustment_top_single_and_all20.csv"))
rm(pd); gc()

## (3) overlap with the biochemical discordance measure (same steps and seed as script 36)
set.seed(20260923)
o <- fread(file.path(W, "output", "reduced_scores", "topk_oof_predictions.csv"))
bd <- merge(d[, .(id, Age, Sex, WC, BMI, hba1c, hdl, tg, crp, tdi, smoking, dlt, E11, I50, N18, K76, I25)], o[, .(id, dlt_k20)], by = "id")
bd <- bd[complete.cases(bd[, .(Age, Sex, WC, BMI, hba1c, hdl, tg, crp, tdi, smoking)]) & tg > 0 & crp > 0]
say("biochemical comparison n = %d", nrow(bd))
folds <- createFolds(bd$WC, k = 10); bd[, bWC := NA_real_]
for (f in folds) { m <- lm(WC ~ Age + Sex + ns(hba1c,3) + ns(hdl,3) + ns(log(tg),3) + ns(log(crp),3), bd[-f]); bd[f, bWC := predict(m, bd[f])] }
bd[, bdlt := bWC - fitted(lm(bWC ~ WC, bd))]
say("r(proWCdelta, biochemical) = %.3f", cor(bd$dlt, bd$bdlt))
bd[, `:=`(topP = dlt >= quantile(dlt, 0.9), topB = bdlt >= quantile(bdlt, 0.9))]
bd[, grp := factor(fifelse(topP & topB, "both", fifelse(topP, "proteomic only", fifelse(topB, "biochemical only", "neither"))), levels = c("neither", "proteomic only", "biochemical only", "both"))]
o3 <- list(); for (cd in c("E11", "I50", "N18", "K76", "I25")) { m <- summary(glm(as.formula(paste(cd, "~ grp +", BASE)), bd, family = binomial()))$coefficients
  for (g in levels(bd$grp)) { k <- paste0("grp", g); o3[[length(o3) + 1]] <- data.table(disease = cd, group = g, n = sum(bd$grp == g), events = sum(bd[grp == g][[cd]]), risk = mean(bd[grp == g][[cd]]),
    OR = if (g == "neither") 1 else exp(m[k, 1]), lo = if (g == "neither") NA else exp(m[k, 1] - 1.96 * m[k, 2]), hi = if (g == "neither") NA else exp(m[k, 1] + 1.96 * m[k, 2])) } }
t3 <- rbindlist(o3); print(t3[, .(disease, group, n, risk = round(100 * risk, 1), txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))]); fwrite(t3, file.path(TB, "T60_overlap_proteomic_biochemical_top_decile.csv"))

## (4) matched pairs: absolute risks and paired body composition (same matching as 08_bodycomp_anchoring.R)
d[, thr := fifelse(Sex == 0, 88, 102)]
d[, grp := fifelse(WC < thr & proWC < thr, "N/N", fifelse(WC < thr & proWC >= thr, "N/H", fifelse(WC >= thr & proWC < thr, "H/N", "H/H")))]
m <- d[grp %in% c("N/N", "N/H")]; m[, NH := as.integer(grp == "N/H")]
say("matching set: normal WC %d, discordant %d", nrow(m), sum(m$NH))
set.seed(123)
mt <- matchit(NH ~ Age + WC + BMI, data = as.data.frame(m), method = "nearest", exact = ~Sex, distance = "mahalanobis", caliper = c(WC = 0.15, BMI = 0.15), std.caliper = TRUE, ratio = 1)
md <- as.data.table(match.data(mt)); say("pairs %d; unmatched discordant %d", sum(md$NH == 1), sum(m$NH) - sum(md$NH == 1))
o4 <- list(); for (cd in c("E11", "I50", "N18")) { w <- dcast(md[, .(subclass, NH, y = get(cd))], subclass ~ NH, value.var = "y"); setnames(w, c("subclass", "NN", "NH"))
  dif <- w$NH - w$NN; se <- sd(dif) / sqrt(nrow(w))
  o4[[length(o4) + 1]] <- data.table(outcome = cd, pairs = nrow(w), risk_NH = mean(w$NH), risk_NN = mean(w$NN), risk_difference = mean(dif), lo = mean(dif) - 1.96 * se, hi = mean(dif) + 1.96 * se) }
bc <- fread(file.path(RAW, "bodycomp_subset.csv"), select = c("Participant ID", "Body fat percentage | Instance 0", "Trunk fat percentage | Instance 0", "Whole body fat mass | Instance 0", "Trunk fat mass | Instance 0", "Whole body fat-free mass | Instance 0", "Trunk fat-free mass | Instance 0", "Whole body water mass | Instance 0"))
setnames(bc, c("id", "bfp", "tfp", "fm", "tfm", "ffm", "tffm", "water")); for (v in names(bc)[-1]) bc[[v]] <- num(bc[[v]])
mb <- merge(md[, .(id, NH, subclass)], bc, by = "id", all.x = TRUE)
for (v in c("bfp", "tfp", "fm", "tfm", "ffm", "tffm", "water")) { w <- dcast(mb[, .(subclass, NH, y = get(v))], subclass ~ NH, value.var = "y"); setnames(w, c("subclass", "NN", "NH")); w <- w[complete.cases(w)]
  tt <- t.test(w$NH, w$NN, paired = TRUE)
  o4[[length(o4) + 1]] <- data.table(outcome = v, pairs = nrow(w), risk_NH = mean(w$NH), risk_NN = mean(w$NN), risk_difference = unname(tt$estimate), lo = tt$conf.int[1], hi = tt$conf.int[2]) }
t4 <- rbindlist(o4); t4[, outcome := c(E11 = "Type 2 diabetes (10-year risk)", I50 = "Heart failure (10-year risk)", N18 = "Chronic kidney disease (10-year risk)", bfp = "Body fat percentage (points)", tfp = "Trunk fat percentage (points)",
  fm = "Whole-body fat mass (kg)", tfm = "Trunk fat mass (kg)", ffm = "Whole-body fat-free mass (kg)", tffm = "Trunk fat-free mass (kg)", water = "Whole-body water mass (kg)")[outcome]]
setnames(t4, c("risk_NH", "risk_NN", "risk_difference"), c("mean_discordant", "mean_concordant", "paired_difference")); print(t4); fwrite(t4, file.path(TB, "T61_matched_pairs_absolute_risk_and_paired_bodycomp.csv"))

## (5) record-linkage check: sex predicted from proteins (5-fold cross-validated LASSO logistic)
pp <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(pp, 1, "id"); prot <- setdiff(names(pp), c("id", "WC", "Age", "Sex"))
ys <- as.integer(pp$Sex); sdiff <- sapply(prot, function(v) { x <- pp[[v]]; abs(mean(x[ys == 1]) - mean(x[ys == 0])) / sd(x) }); top <- names(sort(sdiff, decreasing = TRUE))[1:30]
say("most sex-differentiated proteins: %s", paste(head(top, 8), collapse = ", "))
X <- as.data.frame(pp[, ..top]); X$ys <- ys; set.seed(7); fo <- sample(rep(1:5, length.out = nrow(pp))); ps <- numeric(nrow(pp))
for (k in 1:5) { f <- glm(ys ~ ., data = X[fo != k, ], family = binomial()); ps[fo == k] <- predict(f, X[fo == k, ], type = "response") }
mis <- mean((ps >= 0.5) != (ys == 1)); t5 <- data.table(n = nrow(pp), agreement = 1 - mis, discordant = sum((ps >= 0.5) != (ys == 1)), discordant_confident = sum((ps > 0.99 & ys == 0) | (ps < 0.01 & ys == 1)))
print(t5); fwrite(t5, file.path(TB, "T62_protein_predicted_sex_check.csv")); rm(pp, X); gc()

## (6) baseline characteristics (Supplementary Table 4)
a <- all0; a[, thr := fifelse(Sex == 0, 88, 102)]
a[, grp := fifelse(WC < thr & proWC < thr, "Normal WC + normal proWC", fifelse(WC < thr & proWC >= thr, "Normal WC + high proWC", fifelse(WC >= thr & proWC < thr, "High WC + normal proWC", "High WC + high proWC")))]
a[, region := fifelse(as.integer(centre) %in% c(11004, 11005, 11003, 11022, 11023), "Scotland and Wales", "England")]
a[, smk := as.character(smk_raw)]
msd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)); miq <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
pct <- function(x) sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))
prevI <- function(x, tl, b0) { if (x == "") return(FALSE); cs <- strsplit(x, "[|]")[[1]]; ds <- as.IDate(strsplit(tl, "[|]")[[1]]); any(grepl("^I", cs) & !is.na(ds) & ds < b0) }
a[, prior_cvd := mapply(prevI, dxall, tall, b0)]
desc <- function(x) data.table(characteristic = c("N", "Age, years, mean (SD)", "Women", "BMI, kg/m2, mean (SD)", "Waist circumference, cm, mean (SD)", "proWC, cm, mean (SD)", "Townsend deprivation index, mean (SD)",
    "Current smoking", "White ethnicity", "Lipid-lowering treatment", "Antihypertensive treatment", "Insulin treatment", "Hospital-recorded circulatory diagnosis before baseline", "HbA1c, mmol/mol, median (IQR)", "C-reactive protein, mg/L, median (IQR)", "Triglycerides, mmol/L, median (IQR)"),
  value = c(format(nrow(x), big.mark = ","), msd(x$Age), pct(x$Sex == 0), msd(x$BMI), msd(x$WC), msd(x$proWC), msd(x$tdi), pct(x$smk == "2"), pct(x$white == 1), pct(x$lipid_med == 1), pct(x$bp_med == 1), pct(x$insulin == 1), pct(x$prior_cvd), miq(x$hba1c), miq(x$crp), miq(x$tg)))
cols <- list("All participants" = a, "Normal WC + normal proWC" = a[grp == "Normal WC + normal proWC"], "Normal WC + high proWC" = a[grp == "Normal WC + high proWC"], "High WC + normal proWC" = a[grp == "High WC + normal proWC"],
             "High WC + high proWC" = a[grp == "High WC + high proWC"], "England (development)" = a[region == "England"], "Scotland and Wales (hold-out)" = a[region == "Scotland and Wales"])
t6 <- Reduce(function(x, y) merge(x, y, by = "characteristic", sort = FALSE), lapply(names(cols), function(nm) { z <- desc(cols[[nm]]); setnames(z, "value", nm); z }))
print(t6, width = 300); fwrite(t6, file.path(TB, "T63_baseline_characteristics_supplementary_table4.csv"))
say("smoking status codes: %s", paste(names(table(a$smk)), collapse = ","))

## (7) anchoring estimated within each sex
d[, pa := NA_real_]; for (sx in 0:1) { k <- d$Sex == sx; ft <- lm(pWC ~ WC, d[k]); d[k, pa := WC + pWC - predict(ft, d[k])] }
say("mean proWCdelta by sex with sex-specific anchoring: women %.2f, men %.2f", mean((d$pa - d$WC)[d$Sex == 0]), mean((d$pa - d$WC)[d$Sex == 1]))
nw <- d[WC < thr]; nw[, hi := as.integer(pa >= thr)]
o7 <- list(); for (cd in c("E11", "I50", "N18")) { r <- fz(cd, BASE, nw, term = "hi"); o7[[length(o7) + 1]] <- data.table(disease = cd, n_normalWC = nrow(nw), n_high = sum(nw$hi), male_pct_high = 100 * mean(nw[hi == 1]$Sex == 1), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t7 <- rbindlist(o7); print(t7); fwrite(t7, file.path(TB, "T64_sex_specific_anchoring_threshold.csv"))
say("DONE")
