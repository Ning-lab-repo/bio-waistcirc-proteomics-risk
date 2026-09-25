## 64_reliability_missingness_egfr_overlap.R
## (1) repeat-measurement reliability of the anchoring traits (baseline versus first repeat assessment): WC, hip
##     circumference, waist-to-hip ratio and BMI (body fat percentage was measured only at baseline);
## (2) primary estimates in participants with the complete protein panel and in those missing the block of 1,461
##     proteins (Olink processing batch 7), with a Wald test of the difference;
## (3) heart failure by kidney function (estimated glomerular filtration rate above or below 60 mL/min/1.73 m2) and
##     with adjustment for eGFR;
## (4) deaths within 10 years without the endpoint (competing events) and administrative censoring in the primary sets;
## (5) overlap between the proteins retained by the proWC model and by the diabetes-trained score (single fit at the
##     median lambda of the nested cross-validation).
## Output: T88 (reliability), T89 (missingness strata), T90 (heart failure and kidney function), T91 (competing events),
##         T92 (protein overlap).
suppressPackageStartupMessages({ library(data.table); library(splines); library(glmnet) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HC", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "HC", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bsz <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"), select = c("Participant ID", "Waist circumference | Instance 1", "Hip circumference | Instance 1", "Body mass index (BMI) | Instance 1")); setnames(bsz, c("id", "wc1", "hc1", "bmi1"))
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; bf <- file.path(F8, "wc_pro_Batch.csv")
bat <- fread(bf, select = c(names(fread(bf, nrows = 0))[1], "Batch")); setnames(bat, c("id", "batch"))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr, bsz, bat))
for (v in c("Age", "WC", "BMI", "HC", "hba1c", "cre", "cys", "wc1", "hc1", "bmi1")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; d[, z := dlt / sd(dlt)]; say("modelling set n = %d", nrow(d))

## ---------------- (1) repeat-measurement reliability ----------------
rel <- list(); add <- function(trait, x0, x1) { k <- !is.na(x0) & !is.na(x1)
  rel[[length(rel) + 1]] <<- data.table(trait = trait, n = sum(k), r = cor(x0[k], x1[k]), sd_baseline = sd(x0[k]), sd_within = sd(x1[k] - x0[k]) / sqrt(2), mean_change = mean(x1[k] - x0[k])) }
add("Waist circumference (cm)", d$WC, d$wc1); add("Hip circumference (cm)", d$HC, d$hc1)
add("Waist-to-hip ratio", d$WC / d$HC, d$wc1 / d$hc1); add("Body mass index (kg/m2)", d$BMI, d$bmi1)
t88 <- rbindlist(rel); print(t88); fwrite(t88, file.path(TB, "T88_repeat_measurement_reliability.csv"))

## ---------------- endpoint definitions (as 57/59) ----------------
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
ci <- function(m, term) { s <- summary(m)$coefficients[term, ]; c(OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]]), b = s[[1]], se = s[[2]]) }

## ---------------- (2) complete panel versus the missing block ----------------
mf <- fread(file.path(F12, "analysis_data_WC.csv"))
setnames(mf, 1, "id"); prot <- setdiff(names(mf), c("id", "WC", "Age", "Sex"))
say("batch 7 in modelling set: %d", sum(d$batch == 7, na.rm = TRUE))
o2 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[get(paste0("free_", cd)) == TRUE & !is.na(batch)]
  for (g in list(c("complete panel (batches 0-6)", quote(batch != 7)), c("missing block (batch 7)", quote(batch == 7)))) {
    y <- x[eval(g[[2]])]; r <- ci(glm(as.formula(paste(cd, "~ z +", BASE)), y, family = binomial()), "z")
    o2[[length(o2) + 1]] <- data.table(disease = dis$disease[i], group = g[[1]], n = nrow(y), events = sum(y[[cd]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], b = r[["b"]], se = r[["se"]]) } }
t89 <- rbindlist(o2); t89[, p_difference := NA_real_]
for (dd in unique(t89$disease)) { aa <- t89[disease == dd & group %like% "complete"]; bb <- t89[disease == dd & group %like% "missing"]
  pval <- 2 * pnorm(-abs(aa$b - bb$b) / sqrt(aa$se^2 + bb$se^2)); t89[disease == dd, p_difference := pval] }
print(t89[, .(disease, group, n, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi), p_difference = signif(p_difference, 2))])
fwrite(t89[, .(disease, group, n, events, OR, lo, hi, p_difference)], file.path(TB, "T89_missing_block_strata.csv"))

## ---------------- (3) heart failure and kidney function ----------------
o3 <- list(); x <- d[free_I50 == TRUE]
r <- ci(glm(as.formula(paste("I50 ~ z +", BASE)), x, family = binomial()), "z"); o3[[1]] <- data.table(analysis = "all participants free of heart failure", n = nrow(x), events = sum(x$I50), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
xe <- x[!is.na(egfr)]
r <- ci(glm(as.formula(paste("I50 ~ z + ns(egfr,3) +", BASE)), xe, family = binomial()), "z"); o3[[2]] <- data.table(analysis = "adjusted for eGFR", n = nrow(xe), events = sum(xe$I50), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
for (g in list(c("eGFR >= 60", quote(egfr >= 60)), c("eGFR < 60", quote(egfr < 60)))) { y <- xe[eval(g[[2]])]
  r <- ci(glm(as.formula(paste("I50 ~ z +", BASE)), y, family = binomial()), "z")
  o3[[length(o3) + 1]] <- data.table(analysis = g[[1]], n = nrow(y), events = sum(y$I50), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t90 <- rbindlist(o3); print(t90); fwrite(t90, file.path(TB, "T90_heart_failure_kidney_function.csv"))

## ---------------- (4) competing events within 10 years ----------------
mt <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure10/WC_death_time.csv", select = c("V1", "all_cause_death", "survival_time")); setnames(mt, 1, "id")
e <- merge(d, mt, by = "id", all.x = TRUE); o4 <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- e[get(paste0("free_", cd)) == TRUE]
  o4[[i]] <- data.table(disease = dis$disease[i], n = nrow(x), events = sum(x[[cd]]),
    deaths_within_10y_without_event = sum(x$all_cause_death == 1 & x$survival_time < 10 & x[[cd]] == 0, na.rm = TRUE),
    median_followup_y = median(x$survival_time, na.rm = TRUE)) }
t91 <- rbindlist(o4); print(t91); fwrite(t91, file.path(TB, "T91_competing_deaths.csv"))

## ---------------- (5) protein overlap of the two scores ----------------
ph2 <- d[, .(id, y = fifelse(free_E11, as.numeric(E11), NA_real_))]
m2 <- merge(mf, ph2, by = "id"); m2 <- m2[!is.na(y)]
X <- cbind(as.matrix(m2[, ..prot]), Age = num(m2$Age), Sex = as.integer(m2$Sex))
fit <- glmnet(X, m2$y, family = "binomial", alpha = 1, standardize = TRUE, lambda = 0.00117534586458836)
cf <- as.matrix(coef(fit)); sel <- rownames(cf)[cf[, 1] != 0]; sel <- setdiff(sel, c("(Intercept)", "Age", "Sex"))
wc <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data/model_proWC_full_coefficients.csv")
wcp <- setdiff(wc[[1]], c("(Intercept)", "Age", "Sex")); ov <- intersect(sel, wcp)
say("diabetes score: %d proteins; proWC model: %d; overlap %d (%.0f%% of the diabetes score)", length(sel), length(wcp), length(ov), 100 * length(ov) / length(sel))
t92 <- data.table(quantity = c("proteins in the diabetes-trained score", "proteins in the proWC model", "proteins in both", "share of the diabetes score also in proWC (%)"),
  value = c(length(sel), length(wcp), length(ov), 100 * length(ov) / length(sel)))
print(t92); fwrite(t92, file.path(TB, "T92_protein_overlap_scores.csv"))
fwrite(data.table(protein = sel, in_proWC_model = sel %in% wcp), file.path(TB, "T92b_diabetes_score_proteins.csv"))
say("DONE")
