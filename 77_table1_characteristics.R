## 77_table1_characteristics.R
## Table 1. Baseline characteristics of the proteomic cohort and of the participants who later had abdominal MRI.
## Everything is measured at the baseline visit, when the blood for proteomics was drawn, except the imaging measures
## and the interval to the scan, which are given for the imaging subsample only.
## Output: T113_table1.csv (formatted cells, ready for the manuscript table)
suppressPackageStartupMessages({ library(data.table) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HBA1C",
  "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "yu_ten_need_diagnosis", "more_ten_new_diagnosis"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "hba1c", "dxall", "selfrep", "dx10", "dxmore"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name",
  "Medication for cholesterol, blood pressure or diabetes | Instance 0",
  "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
la  <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.character(get("Smoking status")),
               eth = as.character(get("Ethnic background")))]
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
rec <- fread(file.path(RAW, "2Recruitment.csv"), select = c("Participant ID", "Date of attending assessment centre | Instance 0",
             "Date of attending assessment centre | Instance 2")); setnames(rec, c("id", "date0", "date2"))
## ethnic background as recorded (study master file): 1 = White; 0 = any other group, including "prefer not to answer" (-3)
## and "do not know" (-1), which the covariate file had filled
MF <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8/wc_pro_Batch.csv"
eth <- fread(MF, select = c(names(fread(MF, nrows = 0))[1], "ethnicity", "Ethnic background | Instance 0")); setnames(eth, c("id", "white_rec", "eth_code"))

d <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(d, cov, mri, rec, fr, eth))
for (v in c("Age", "WC", "BMI", "hba1c", "vat", "asat")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, dlt)])]
lo <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
## diabetes before the blood sample, by the same sources as the primary analysis (57_primary_incident_table2.R):
## a hospital record dated before baseline, self-report, glucose-lowering medication or HbA1c of 48 mmol/mol or more
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
d[, `:=`(dx10 = toupper(lo(dx10)), dxmore = toupper(lo(dxmore)), dxall = toupper(lo(dxall)), med = paste(lo(med_m), lo(med_f)), medname = lo(medname))]
d[, diab := (grepl("(^|[|])E11", dxall) & !(grepl("(^|[|])E11", dx10) | grepl("(^|[|])E11", dxmore))) |
            grepl("(^|[|])(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)([|]|$)", lo(selfrep), perl = TRUE) |
            grepl("insulin", med) | grepl(gd, medname) | (!is.na(hba1c) & hba1c >= 48)]
d[, mri := !is.na(vat) & !is.na(asat) & asat > 0]
d[, years := as.numeric(as.IDate(date2) - as.IDate(date0)) / 365.25]
say("proteomic cohort %d; with abdominal MRI %d", nrow(d), sum(d$mri))

msd <- function(x, k = 1) sprintf(paste0("%.", k, "f (%.", k, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
mid <- function(x, k = 1) sprintf(paste0("%.", k, "f (%.", k, "f–%.", k, "f)"), median(x, na.rm = TRUE),
                                  quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
np  <- function(b) sprintf("%s (%.1f)", format(sum(b, na.rm = TRUE), big.mark = ","), 100 * mean(b, na.rm = TRUE))
col <- function(x, img) c(
  "Participants" = format(nrow(x), big.mark = ","),
  "Age, years" = msd(x$Age),
  "Women" = np(x$Sex == 0),
  "White ethnic background" = np(x$white_rec == 1),     ## as recorded
  "Ethnic background not stated" = np(x$eth_code %in% c(-1, -3) | is.na(x$eth_code)),
  "Townsend deprivation index" = msd(x$tdi, 2),
  "Current smoker" = np(x$smoking == "2"),              ## UK Biobank coding: 0 never, 1 previous, 2 current
  "Body-mass index, kg/m²" = msd(x$BMI),
  "Waist circumference, cm" = msd(x$WC),
  "proWCΔ, cm" = msd(x$dlt, 2),
  "HbA1c, mmol/mol" = mid(x$hba1c),
  "Diabetes at baseline" = np(x$diab),
  "Years from blood sample to MRI" = if (img) mid(x$years) else "–",
  "Visceral adipose tissue, L" = if (img) msd(x$vat, 2) else "–",
  "Abdominal subcutaneous adipose tissue, L" = if (img) msd(x$asat, 2) else "–")
t1 <- data.table(characteristic = names(col(d, FALSE)), proteomic_cohort = col(d, FALSE), mri_subsample = col(d[mri == TRUE], TRUE))
print(t1, width = 160)
fwrite(t1, file.path(TB, "T113_table1.csv"))
say("DONE")
