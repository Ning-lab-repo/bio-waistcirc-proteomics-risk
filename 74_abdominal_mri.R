## 74_abdominal_mri.R
## Where does the fat sit that the proteome registers and the tape measure does not? Abdominal magnetic resonance
## imaging measures visceral and abdominal subcutaneous adipose tissue directly.
##
## Timing. The proteins were measured in blood drawn at the baseline visit (instance 0), and proWCdelta is defined
## against the waist measured at that visit, so the design-consistent comparison holds the BASELINE waist and BMI
## fixed, as every other analysis in the paper does. The imaging was done at a later visit (instance 2). The time
## between the two is handled by repeating the analysis in participants whose waist and BMI barely changed between
## the visits: in them the later scan approximates the body at the time of the blood sample, and the baseline and
## imaging-visit anthropometry coincide. The imaging-visit adjustment is not used as the primary analysis, because it
## conditions a baseline exposure on a body size measured years later.
##
## The imaging block is checked against external quantities before use (an earlier extract was attached to the
## wrong participants and passed only internal checks).
## Output: T110_abdominal_mri.csv
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
BC  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ---------------- data ----------------
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2",
                        "VAT..visceral.adipose.tissue..mass...Instance.2"))
setnames(mri, c("id", "vat", "asat", "dxa_vat"))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI"))
la  <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
img <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"),
             select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2"))
setnames(img, c("id", "WC2", "BMI2"))
rec <- fread(file.path(RAW, "2Recruitment.csv"),
             select = c("Participant ID", "Date of attending assessment centre | Instance 0",
                        "Date of attending assessment centre | Instance 2"))
setnames(rec, c("id", "date0", "date2"))

d <- Reduce(function(a, b) merge(a, b, by = "id"), list(phen, la[, .(id, dlt = BioX_Delta)], mri))
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(d, img, rec))
for (v in c("Age", "WC", "BMI", "vat", "asat", "dxa_vat", "WC2", "BMI2")) d[[v]] <- num(d[[v]])
d[, Sex := as.integer(Sex)]

## ---------------- linkage check ----------------
r_w <- cor(d$vat, d$WC, use = "pair"); r_d <- cor(d$vat, d$dxa_vat, use = "pair")
say("LINKAGE CHECK: r(MRI visceral fat, baseline waist) %+.3f; r(MRI, DXA visceral fat) %+.3f", r_w, r_d)
stopifnot(r_w > 0.4, r_d > 0.5)

## per SD of proWCdelta in the whole proteomic cohort, the unit used for every disease association in the paper
coh <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI")) coh[[v]] <- num(coh[[v]])
sdz <- sd(coh[complete.cases(coh[, .(Age, Sex, WC, BMI, dlt)])]$dlt)
d <- d[!is.na(dlt) & !is.na(WC) & !is.na(BMI) & !is.na(Age) & !is.na(Sex) & !is.na(vat) & !is.na(asat) & asat > 0]
d[, z := dlt / sdz]
d[, ratio := vat / asat]
d[, years := as.numeric(as.IDate(date2) - as.IDate(date0)) / 365.25]
d[, `:=`(dWC = WC2 - WC, dBMI = BMI2 - BMI)]

## ---------------- timing ----------------
say("")
say("participants with a proteomic score and abdominal MRI: %d", nrow(d))
say("years from the blood sample to the scan: median %.1f (IQR %.1f to %.1f)",
    median(d$years, na.rm = TRUE), quantile(d$years, .25, na.rm = TRUE), quantile(d$years, .75, na.rm = TRUE))
say("change in waist over that interval: median %+.1f cm (IQR %+.1f to %+.1f); change in BMI median %+.2f",
    median(d$dWC, na.rm = TRUE), quantile(d$dWC, .25, na.rm = TRUE), quantile(d$dWC, .75, na.rm = TRUE), median(d$dBMI, na.rm = TRUE))
say("correlation of proWCdelta with the later change in waist: %+.3f", cor(d$z, d$dWC, use = "pair"))
## stable: waist within 5 cm (about one within-person SD of repeat measurement, 4.96 cm) and BMI within 1.5 kg/m2
d[, stable := !is.na(dWC) & !is.na(dBMI) & abs(dWC) < 5 & abs(dBMI) < 1.5]
say("participants whose waist and BMI were stable between the visits: %d (%.0f%%)", sum(d$stable), 100 * mean(d$stable))

## ---------------- associations ----------------
OUT <- c(vat = "visceral adipose tissue (L)", asat = "abdominal subcutaneous adipose tissue (L)",
         ratio = "visceral to subcutaneous ratio")
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
fit <- function(y, dat, lab) {
  m  <- lm(as.formula(paste(y, "~ z +", B0)), dat); co <- summary(m)$coefficients
  mw <- lm(as.formula(paste(y, "~ scale(WC) + Age + factor(Sex)")), dat); cw <- summary(mw)$coefficients
  data.table(outcome = OUT[[y]], population = lab, n = nrow(dat),
             beta = co["z", 1], lo = co["z", 1] - 1.96 * co["z", 2], hi = co["z", 1] + 1.96 * co["z", 2], p = co["z", 4],
             per_SD_measured_WC = cw[2, 1]) }
res <- rbindlist(lapply(names(OUT), function(y) rbindlist(list(
  fit(y, d,          "all, baseline waist and BMI held fixed"),
  fit(y, d[stable == TRUE], "waist and BMI stable between the visits")))))
res[, txt := sprintf("%+.3f (%+.3f to %+.3f)", beta, lo, hi)]
say("")
print(res[, .(outcome, population, n, `per SD of proWCdelta` = txt, p = signif(p, 2),
              `per SD of measured WC` = sprintf("%+.3f", per_SD_measured_WC))], width = 220)
fwrite(res, file.path(TB, "T110_abdominal_mri.csv"))
say("DONE")
