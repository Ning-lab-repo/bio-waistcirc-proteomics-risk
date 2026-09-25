## 73_observed_covariates.R
## The covariate file used throughout carries values that were filled during data preparation, and the manuscript
## could only describe them rather than name the method. The observed values survive in the study's master file, so
## the primary models can simply be refitted in the participants who have them, and the extended adjustment in the
## participants who have education and physical activity as recorded.
## This replaces a description of an unknown procedure with a test of whether it mattered.
## Output: T109_observed_covariates.csv.
suppressPackageStartupMessages({ library(data.table); library(splines) })
M   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8/wc_pro_Batch.csv"
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

raw <- fread(M, select = c("Participant ID", "TDI", "Smoking status", "Education", "MET Physical activity"), showProgress = FALSE)
setnames(raw, c("id", "tdi_o", "smoke_o", "edu_o", "met_o"))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI",
  "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0",
  "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name",
  "Medication for cholesterol, blood pressure or diabetes | Instance 0",
  "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE),
            list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr, raw))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys", "tdi_o", "met_o")) d[[v]] <- num(d[[v]])
d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]
d[, z := dlt / sd(dlt)]

s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)),
         selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
prev <- function(cd) grepl(paste0("(^|[|])", cd), d$dxall) & !(grepl(paste0("(^|[|])", cd), d$dx10) | grepl(paste0("(^|[|])", cd), d$dxmore))
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) *
            0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
d[, free_E11 := !(prev("E11") | has(selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") |
                  has(med, "insulin") | has(medname, gd) | (!is.na(hba1c) & hba1c >= 48))]
d[, free_I50 := !(prev("I50") | has(selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"))]
d[, free_N18 := !(prev("N18") | has(selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") |
                  (!is.na(egfr) & egfr < 60))]
d[, `:=`(E11 = as.integer(grepl("(^|[|])E11", dx10)), I50 = as.integer(grepl("(^|[|])I50", dx10)),
         N18 = as.integer(grepl("(^|[|])N18", dx10)))]

blank <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA")
d[, ok_base := !blank(tdi_o) & !blank(smoke_o)]
d[, ok_ext  := ok_base & !blank(edu_o) & !blank(met_o)]
say("observed deprivation and smoking: %d of %d (%.1f%%)", sum(d$ok_base), nrow(d), 100 * mean(d$ok_base))
say("also observed education and physical activity: %d (%.1f%%)", sum(d$ok_ext), 100 * mean(d$ok_ext))

fit <- function(cd, rhs, keep, lab) {
  x <- d[which(d[[paste0("free_", cd)]] & keep)]
  m <- glm(as.formula(paste(cd, "~ z +", rhs)), x, family = binomial()); co <- summary(m)$coefficients
  data.table(endpoint = cd, model = lab, n = nrow(x), events = sum(x[[cd]]),
             OR = exp(co["z", 1]), lo = exp(co["z", 1] - 1.96 * co["z", 2]), hi = exp(co["z", 1] + 1.96 * co["z", 2])) }

AS_USED <- "Age + factor(Sex) + tdi + factor(smoking) + ns(WC,3) + ns(BMI,3)"           ## the values the manuscript used
OBSERVED <- "Age + factor(Sex) + tdi_o + factor(smoke_o) + ns(WC,3) + ns(BMI,3)"        ## the recorded values only
EXTENDED <- paste(OBSERVED, "+ factor(edu_o) + met_o")
res <- rbindlist(lapply(c("E11", "I50", "N18"), function(cd) rbindlist(list(
  fit(cd, AS_USED,  rep(TRUE, nrow(d)), "values as used in the manuscript"),
  fit(cd, OBSERVED, d$ok_base,          "recorded values only"),
  fit(cd, EXTENDED, d$ok_ext,           "recorded values, with education and physical activity")))))
res[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]
NM <- c(E11 = "Type 2 diabetes", I50 = "Heart failure", N18 = "Chronic kidney disease")
res[, endpoint := NM[endpoint]]
print(res[, .(endpoint, model, n, events, `OR per SD` = txt)], width = 200)
fwrite(res, file.path(TB, "T109_observed_covariates.csv"))
say("DONE")
