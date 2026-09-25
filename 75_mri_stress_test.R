## 75_mri_stress_test.R
## Before the paper claims that, at the same measured waist and BMI, a higher proteomic waist marks more visceral fat,
## every explanation other than a real difference in fat distribution is tested here.
##
## 1. Measurement error in the tape. proWCdelta is high partly where the baseline waist was measured too small, and
##    such people truly have larger waists, and hence more abdominal fat, than their measurement says. That would give
##    more fat of the composition a larger waist has at the same BMI; it would not shift the fat towards the viscera.
##    So the visceral share of the extra fat that goes with proWCdelta is compared with the visceral share of the extra
##    fat that goes with a larger measured waist at the same BMI. If the two are equal the finding is tape error; if the
##    proteomic share is clearly higher the finding is not.
## 2. The ten years between the blood sample and the scan: participants whose size was stable, with both visits'
##    waist and BMI held fixed (two measurements of a size that did not change halve the measurement error as well).
## 3. The form of the body-size adjustment.
## 4. Selection into the imaging visit.
## Output: T111_mri_stress_test.csv
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
BC  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI"))
la  <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
bs  <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"),
             select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2",
                        "Standing height | Instance 0"))
setnames(bs, c("id", "WC2", "BMI2", "H0"))
all <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI")) all[[v]] <- num(all[[v]]); all[, Sex := as.integer(Sex)]
all <- all[complete.cases(all[, .(Age, Sex, WC, BMI, dlt)])]
sdz <- sd(all$dlt)                                     ## SD of proWCdelta in the whole proteomic cohort
d <- merge(all, mri, by = "id"); d <- merge(d, bs, by = "id", all.x = TRUE)
for (v in c("vat", "asat", "WC2", "BMI2", "H0")) d[[v]] <- num(d[[v]])
d <- d[!is.na(vat) & !is.na(asat) & asat > 0]
d[, z := dlt / sdz]
d[, `:=`(dWC = WC2 - WC, dBMI = BMI2 - BMI)]
d[, stable := !is.na(dWC) & !is.na(dBMI) & abs(dWC) < 5 & abs(dBMI) < 1.5]
say("participants with a proteomic score and abdominal MRI: %d (stable size between visits: %d)", nrow(d), sum(d$stable))

B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
coefz <- function(y, rhs, x) { co <- summary(lm(as.formula(paste(y, "~ z +", rhs)), x))$coefficients
  c(b = co["z", 1], lo = co["z", 1] - 1.96 * co["z", 2], hi = co["z", 1] + 1.96 * co["z", 2]) }
rows <- list(); add <- function(test, spec, n, v) rows[[length(rows) + 1]] <<- data.table(test = test, specification = spec, n = n,
  estimate = v[["b"]], lo = v[["lo"]], hi = v[["hi"]])

## ---------------- 1. measurement error: the visceral share of the extra fat ----------------
## share = visceral increment / (visceral + subcutaneous increment), for proWCdelta at fixed WC and BMI, and for
## measured WC at fixed BMI (the composition that a larger true waist would bring)
share <- function(x) {
  bv <- coef(lm(as.formula(paste("vat ~ z +", B0)), x))[["z"]]
  bs <- coef(lm(as.formula(paste("asat ~ z +", B0)), x))[["z"]]
  wv <- coef(lm(vat  ~ scale(WC) + ns(BMI, 3) + Age + factor(Sex), x))[[2]]
  ws <- coef(lm(asat ~ scale(WC) + ns(BMI, 3) + Age + factor(Sex), x))[[2]]
  c(pro = bv / (bv + bs), wc = wv / (wv + ws), diff = bv / (bv + bs) - wv / (wv + ws)) }
for (pop in c("all", "stable")) {
  x <- if (pop == "all") d else d[stable == TRUE]
  s0 <- share(x); bt <- t(replicate(500, share(x[sample(.N, .N, replace = TRUE)])))
  ci <- function(k) quantile(bt[, k], c(.025, .975))
  say("")
  say("[1] visceral share of the extra abdominal fat (%s, n = %d)", pop, nrow(x))
  say("    with proWCdelta, at the same measured waist and BMI : %.1f%% (%.1f to %.1f)", 100 * s0[["pro"]], 100 * ci("pro")[1], 100 * ci("pro")[2])
  say("    with a larger measured waist, at the same BMI       : %.1f%% (%.1f to %.1f)", 100 * s0[["wc"]], 100 * ci("wc")[1], 100 * ci("wc")[2])
  say("    difference                                          : %+.1f points (%+.1f to %+.1f)", 100 * s0[["diff"]], 100 * ci("diff")[1], 100 * ci("diff")[2])
  add("1. tape measurement error", paste("visceral share, proWCdelta:", pop), nrow(x), c(b = s0[["pro"]], lo = ci("pro")[[1]], hi = ci("pro")[[2]]))
  add("1. tape measurement error", paste("visceral share, larger measured waist:", pop), nrow(x), c(b = s0[["wc"]], lo = ci("wc")[[1]], hi = ci("wc")[[2]]))
  add("1. tape measurement error", paste("difference in visceral share:", pop), nrow(x), c(b = s0[["diff"]], lo = ci("diff")[[1]], hi = ci("diff")[[2]]))
}

## ---------------- 2. timing and a second measurement of size ----------------
st <- d[stable == TRUE & !is.na(WC2) & !is.na(BMI2)]
say("")
say("[2] visceral fat per SD of proWCdelta")
v <- coefz("vat", B0, d);  say("    all, baseline size fixed                                  %+.3f L (%+.3f to %+.3f)", v[1], v[2], v[3]); add("2. timing", "all, baseline size fixed", nrow(d), v)
v <- coefz("vat", B0, st); say("    stable size, baseline size fixed                          %+.3f L (%+.3f to %+.3f)", v[1], v[2], v[3]); add("2. timing", "stable, baseline size fixed", nrow(st), v)
v <- coefz("vat", paste(B0, "+ ns(WC2,3) + ns(BMI2,3)"), st)
say("    stable size, size at BOTH visits fixed                    %+.3f L (%+.3f to %+.3f)", v[1], v[2], v[3]); add("2. timing", "stable, size at both visits fixed", nrow(st), v)
v <- coefz("log(vat/asat)", B0, st)
say("    stable size, log visceral:subcutaneous ratio              %+.3f (%+.3f to %+.3f)", v[1], v[2], v[3]); add("2. timing", "stable, log visceral:subcutaneous", nrow(st), v)

## ---------------- 3. form of the body-size adjustment ----------------
say("")
say("[3] visceral fat per SD of proWCdelta, other adjustments for body size")
for (sp in list(c("splines with 5 degrees of freedom", "Age + factor(Sex) + ns(WC,5) + ns(BMI,5)"),
                c("waist x BMI interaction",           "Age + factor(Sex) + ns(WC,3) * ns(BMI,3)"),
                c("height added",                      "Age + factor(Sex) + ns(WC,3) + ns(BMI,3) + ns(H0,3)"),
                c("within sex, sex-specific splines",  "Age + factor(Sex) * (ns(WC,3) + ns(BMI,3))"))) {
  x <- if (grepl("H0", sp[2])) d[!is.na(H0)] else d
  v <- coefz("vat", sp[2], x); say("    %-40s %+.3f L (%+.3f to %+.3f)", sp[1], v[1], v[2], v[3]); add("3. adjustment form", sp[1], nrow(x), v) }

## ---------------- 3b. confounding, and each sex separately ----------------
F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
cv <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cv, 1, "id")
cv <- cv[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoke = as.factor(get("Smoking status")),
             alcohol = as.factor(get("Alcohol intake")), eth = as.factor(get("Ethnic background")))]
dx <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "s_Diagnoses_ICD10",
            "Noncancer_illne_code_elfreported_Instance0", "HBA1C"))
setnames(dx, c("id", "dxall", "selfrep", "hba1c"))
dc <- merge(merge(d, cv, by = "id"), dx, by = "id")
lo_ <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
## diabetes at baseline (self-report, HbA1c of 48 mmol/mol or more) or a hospital diabetes code (E10-E14) at any time,
## before or after baseline: the exclusion removes everyone whose diabetes could have altered fat or proteins
dc[, diab := grepl("(^|[|])E1[0-4]", toupper(lo_(dxall))) | grepl("diabetes", lo_(selfrep)) | (!is.na(num(hba1c)) & num(hba1c) >= 48)]
dc <- dc[!is.na(tdi) & !is.na(smoke)]
BM <- paste(B0, "+ tdi + smoke")
say("")
say("[3b] visceral fat per SD of proWCdelta: confounding and sex")
for (sp in list(c("+ Townsend index and smoking", BM, "all"),
                c("+ ethnic background",          paste(BM, "+ eth"), "all"),
                c("+ alcohol intake",             paste(BM, "+ alcohol"), "all"),
                c("excluding diabetes at baseline or diagnosed later", BM, "nodiab"),
                c("women",                        "Age + ns(WC,3) + ns(BMI,3) + tdi + smoke", "women"),
                c("men",                          "Age + ns(WC,3) + ns(BMI,3) + tdi + smoke", "men"))) {
  x <- switch(sp[3], all = dc, nodiab = dc[diab == FALSE], women = dc[Sex == 0], men = dc[Sex == 1])
  v <- coefz("vat", sp[2], x); say("    %-40s %+.3f L (%+.3f to %+.3f)  n=%d", sp[1], v[1], v[2], v[3], nrow(x))
  add("3b. confounding and sex", sp[1], nrow(x), v) }

## ---------------- 4. selection into the imaging visit ----------------
say("")
say("[4] who came back for imaging")
all[, imaged := id %in% d$id]
say("    proWCdelta mean: imaged %+.2f cm, not imaged %+.2f cm", mean(all[imaged == TRUE]$dlt), mean(all[imaged == FALSE]$dlt))
m <- glm(imaged ~ scale(dlt) + Age + factor(Sex) + ns(WC,3) + ns(BMI,3), all, family = binomial()); co <- summary(m)$coefficients
say("    odds of attending imaging per SD of proWCdelta, same body size: %.2f (%.2f to %.2f)",
    exp(co[2, 1]), exp(co[2, 1] - 1.96 * co[2, 2]), exp(co[2, 1] + 1.96 * co[2, 2]))
say("    (fewer of those with high proWCdelta returning, because of illness, would bias the visceral association")
say("     towards the null, not away from it)")
add("4. selection", "odds of attending imaging per SD of proWCdelta", nrow(all),
    c(b = exp(co[2, 1]), lo = exp(co[2, 1] - 1.96 * co[2, 2]), hi = exp(co[2, 1] + 1.96 * co[2, 2])))

res <- rbindlist(rows); fwrite(res, file.path(TB, "T111_mri_stress_test.csv")); say("DONE")
