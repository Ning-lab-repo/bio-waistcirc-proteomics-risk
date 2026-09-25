## 61_hf_plate_sparse_sex_mortality.R
## (1) heart failure in participants free of heart failure at baseline: adjustment for cardiac stress and injury proteins
##     (NPPB, NT-proBNP, TNNI3, GDF15), adrenomedullin and leptin, and exclusion of events in the first 2 years;
## (2) Olink plate: share of proWCdelta variance explained by plate and batch; random-intercept logistic model for the
##     primary type 2 diabetes estimate;
## (3) sparser model: nested out-of-fold predictions at lambda.1se (R2, number of proteins, T2D association);
## (4) sex interaction as a Wald test of the difference between the sex-stratified estimates (consistent with them);
## (5) mortality: per SD of proWCdelta at fixed WC and BMI, Cox models for all-cause death and underlying-cause
##     cardiovascular, cancer and diabetes death, after excluding baseline cardiovascular disease, cancer and diabetes and
##     the first 2 years of follow-up.
## Output: T83-T87.
suppressPackageStartupMessages({ library(data.table); library(splines); library(survival); library(lme4); library(glmnet); library(doParallel); library(caret) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "yu_ten_need_time", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "Number of self-reported cancers", "HBA1C", "date_attending_assessment_centre", "new_primary_death"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "t10", "dxmore", "dxall", "selfrep", "ncancer", "hba1c", "d0", "cod"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0")); setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bf <- file.path(F8, "wc_pro_Batch.csv"); bt <- fread(bf, select = c(names(fread(bf, nrows = 0))[1], "PlateID", "Batch")); setnames(bt, c("id", "plate", "batch"))
pr <- fread(file.path(F12, "analysis_data_WC.csv"), select = c("V1", "NPPB", "NTproBNP", "TNNI3", "GDF15", "ADM", "LEP")); setnames(pr, 1, "id")
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id"), cov, fr, bt, pr))
for (v in c("Age", "WC", "BMI", "hba1c", "ncancer")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; d[, z := dlt / sd(dlt)]; say("n = %d", nrow(d))
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname), b0 = as.IDate(gsub("/", "-", d0)))]
prior <- function(pat) grepl(pat, d$dxall, perl = TRUE) & !(grepl(pat, d$dx10, perl = TRUE) | grepl(pat, d$dxmore, perl = TRUE))
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
ci <- function(m, term) { s <- summary(m)$coefficients[term, ]; c(OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]])) }

## (1) heart failure
d[, I50 := as.integer(grepl("(^|[|])I50", dx10))]
d[, free_I50 := !(prior("(^|[|])I50") | grepl("(^|\\|)heart failure/pulmonary odema(\\||$)", selfrep, perl = TRUE))]
first_time <- function(codes, times, code) { if (codes == "") return(NA_real_); cs <- strsplit(codes, "[|]")[[1]]; ts <- strsplit(times, "[|]")[[1]]; k <- which(startsWith(cs, code)); if (!length(k)) return(NA_real_); min(as.numeric(as.IDate(gsub("/", "-", ts[k]))), na.rm = TRUE) }
h <- d[free_I50 == TRUE]; h[, tI50 := (mapply(first_time, dx10, s(t10), MoreArgs = list(code = "I50")) - as.numeric(b0)) / 365.25]
for (v in c("NPPB", "NTproBNP", "TNNI3", "GDF15", "ADM", "LEP")) h[[v]] <- num(h[[v]])
o1 <- list(); add <- list(none = "", NPPB = "+ NPPB", `NT-proBNP` = "+ NTproBNP", TNNI3 = "+ TNNI3", GDF15 = "+ GDF15", `NT-proBNP, TNNI3 and GDF15` = "+ NTproBNP + TNNI3 + GDF15",
  ADM = "+ ADM", LEP = "+ LEP", `ADM and LEP` = "+ ADM + LEP")
for (k in names(add)) { r <- ci(glm(as.formula(paste("I50 ~ z +", BASE, add[[k]])), h, family = binomial()), "z"); o1[[k]] <- data.table(analysis = paste("adjusted for", k), n = nrow(h), events = sum(h$I50), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
h2 <- h[!(I50 == 1 & !is.na(tI50) & tI50 < 2)]; r <- ci(glm(as.formula(paste("I50 ~ z +", BASE)), h2, family = binomial()), "z")
o1[["excl2y"]] <- data.table(analysis = "excluding events in the first 2 years", n = nrow(h2), events = sum(h2$I50), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]])
t83 <- rbindlist(o1); print(t83); fwrite(t83, file.path(TB, "T83_heart_failure_cardiac_proteins.csv"))
say("cor(proWCdelta, NT-proBNP) %.3f; cor(proWCdelta, ADM) %.3f; cor(ADM, LEP) %.3f", cor(h$dlt, h$NTproBNP, use = "complete.obs"), cor(h$dlt, h$ADM, use = "complete.obs"), cor(h$ADM, h$LEP, use = "complete.obs"))

## (2) plate and batch
say("share of proWCdelta variance explained: plate %.4f; batch %.4f", summary(lm(dlt ~ factor(plate), d))$r.squared, summary(lm(dlt ~ factor(batch), d))$r.squared)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
d[, E11 := as.integer(grepl("(^|[|])E11", dx10))]
d[, free_E11 := !(prior("(^|[|])E11") | grepl("(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)", selfrep, perl = TRUE) | grepl("insulin", med) | grepl(gd, medname) | (!is.na(hba1c) & hba1c >= 48))]
x <- d[free_E11 == TRUE]
X1 <- model.matrix(~ ns(WC, 3) + ns(BMI, 3), x)[, -1]; colnames(X1) <- paste0("s", seq_len(ncol(X1))); x <- cbind(x, X1)
gm <- glmer(as.formula(paste("E11 ~ z + Age + Sex + tdi + smoking +", paste(colnames(X1), collapse = " + "), "+ (1 | plate)")), x, family = binomial(), nAGQ = 0)
sg <- summary(gm)$coefficients["z", ]; say("T2D free, random intercept for plate: OR %.2f (%.2f-%.2f); plate SD %.3f", exp(sg[[1]]), exp(sg[[1]] - 1.96 * sg[[2]]), exp(sg[[1]] + 1.96 * sg[[2]]), attr(VarCorr(gm)$plate, "stddev"))
t84 <- data.table(quantity = c("R2 of proWCdelta on plate", "R2 of proWCdelta on batch", "T2D OR, random intercept for plate", "lo", "hi", "plate random-intercept SD"),
  value = c(summary(lm(dlt ~ factor(plate), d))$r.squared, summary(lm(dlt ~ factor(batch), d))$r.squared, exp(sg[[1]]), exp(sg[[1]] - 1.96 * sg[[2]]), exp(sg[[1]] + 1.96 * sg[[2]]), attr(VarCorr(gm)$plate, "stddev")))

## (4) sex interaction, Wald test of stratified estimates (all participants, main definition; as T53)
t53 <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data/source_table_Results_sex_stratified_continuous.csv")
t53[, `:=`(se_m = (log(hi_men) - log(lo_men)) / (2 * 1.96), se_w = (log(hi_women) - log(lo_women)) / (2 * 1.96))]
t53[, p_interaction_stratified := 2 * pnorm(-abs(log(OR_men) - log(OR_women)) / sqrt(se_m^2 + se_w^2))]
print(t53[, .(disease, OR_men = round(OR_men, 2), OR_women = round(OR_women, 2), p_product_term = signif(p_interaction, 2), p_stratified = signif(p_interaction_stratified, 2))])
fwrite(t53[, .(disease, OR_men, lo_men, hi_men, OR_women, lo_women, hi_women, p_interaction_product_term = p_interaction, p_interaction_stratified)], file.path(TB, "T85_sex_interaction_stratified.csv"))

## (5) mortality
pc_cvd <- prior("(^|[|])(I2[0-5]|I50|I6[0-9]|G45)"); pc_can <- prior("(^|[|])C")
mt <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure10/WC_death_time.csv", select = c("V1", "all_cause_death", "cvd_death", "cancer_death", "diabetes_death", "survival_time")); setnames(mt, 1, "id")
e <- merge(d, mt, by = "id"); e[, cod := toupper(s(cod))]
e[, `:=`(cvd_u = as.integer(all_cause_death == 1 & grepl("^I", cod)), cancer_u = as.integer(all_cause_death == 1 & grepl("^C|^D[0-4]", cod)), diab_u = as.integer(all_cause_death == 1 & grepl("^E1[0-4]", cod)))]
e[, prior_cvd := pc_cvd[match(e$id, d$id)] | grepl("(^|\\|)(angina|heart attack/myocardial infarction|stroke|heart failure/pulmonary odema|transient ischaemic attack \\(tia\\))(\\||$)", selfrep, perl = TRUE)]
e[, prior_cancer := pc_can[match(e$id, d$id)] | (!is.na(ncancer) & ncancer > 0)]
e[, excl := prior_cvd | prior_cancer | !free_E11]
COX <- paste(BASE); o5 <- list()
for (oc in list(c("All-cause death", "all_cause_death"), c("Cardiovascular death (anywhere)", "cvd_death"), c("Cardiovascular death (underlying)", "cvd_u"), c("Cancer death (underlying)", "cancer_u"), c("Diabetes-related death (anywhere)", "diabetes_death"), c("Diabetes death (underlying)", "diab_u"))) {
  for (setn in c("all", "restricted")) { x <- if (setn == "all") e else e[excl == FALSE & !(all_cause_death == 1 & survival_time < 2)]
    if (setn == "restricted") { x <- copy(x); x[, st := survival_time - 2]; x <- x[st > 0] } else { x <- copy(x); x[, st := survival_time] }
    cm <- coxph(as.formula(paste("Surv(st,", oc[2], ") ~ z +", COX)), x); sc <- summary(cm)$conf.int["z", ]
    o5[[length(o5) + 1]] <- data.table(outcome = oc[1], sample = if (setn == "all") "all participants" else "excluding baseline CVD, cancer and diabetes and the first 2 years", n = nrow(x), deaths = sum(x[[oc[2]]]), HR = sc[[1]], lo = sc[[3]], hi = sc[[4]]) } }
t86 <- rbindlist(o5); print(t86); fwrite(t86, file.path(TB, "T86_mortality_sensitivity.csv"))

## (3) sparser model at lambda.1se (nested out of fold)
registerDoParallel(cores = 24); set.seed(123)
dat <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(dat, 1, "id"); prot <- setdiff(names(dat), c("id", "WC", "Age", "Sex"))
X <- cbind(as.matrix(dat[, ..prot]), Age = num(dat$Age), Sex = as.integer(dat$Sex)); y <- num(dat$WC); folds <- createFolds(y, k = 10); p1 <- rep(NA_real_, length(y)); nz <- c()
for (k in seq_along(folds)) { te <- folds[[k]]; cvm <- cv.glmnet(X[-te, ], y[-te], alpha = 1, nfolds = 10, parallel = TRUE); p1[te] <- as.numeric(predict(cvm, s = "lambda.1se", newx = X[te, ])); nz <- c(nz, cvm$nzero[which(cvm$lambda == cvm$lambda.1se)]) }
r2 <- 1 - sum((y - p1)^2) / sum((y - mean(y))^2); a <- lm(p1 ~ y); dl <- data.table(id = dat$id, dlt1 = p1 - fitted(a))
x <- merge(d[free_E11 == TRUE], dl, by = "id"); x[, z1 := dlt1 / sd(dlt1)]; r <- ci(glm(as.formula(paste("E11 ~ z1 +", BASE)), x, family = binomial()), "z1")
say("lambda.1se: R2 %.3f, median proteins %.0f, T2D OR %.2f (%.2f-%.2f), r with main proWCdelta %.3f", r2, median(nz), r[["OR"]], r[["lo"]], r[["hi"]], cor(x$dlt, x$dlt1))
t84 <- rbind(t84, data.table(quantity = c("lambda.1se out-of-fold R2", "lambda.1se median number of non-zero terms", "lambda.1se T2D OR", "lo", "hi", "r with main proWCdelta"), value = c(r2, median(nz), r[["OR"]], r[["lo"]], r[["hi"]], cor(x$dlt, x$dlt1))))
fwrite(t84, file.path(TB, "T84_plate_and_sparse_model.csv")); say("DONE")
