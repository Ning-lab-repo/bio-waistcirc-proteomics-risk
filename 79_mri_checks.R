## 79_mri_checks.R
## Further checks of the abdominal MRI findings:
##  1. linkage: MRI visceral fat against waist measured at the imaging visit, which comes from a separate data file
##     (the DXA estimate sits in the same file as the MRI volumes and is therefore not an independent check)
##  2. the WC and BMI of both visits held fixed in all imaging participants, and time to the scan added
##  3. the visceral share of the extra fat within each sex, and with the same (linear) form of the waist adjustment in
##     both arms of the comparison
##  4. partial R2 of proWCdelta for visceral fat at fixed WC and BMI
##  5. visceral fat with routine blood measurements held fixed (triglycerides, HDL cholesterol, HbA1c, C-reactive
##     protein, ALT, GGT): do the proteins add to what routine blood tests show?
##  6. inverse-probability weighting for attendance at imaging; robust standard errors
##  7. per-protein associations with visceral fat at fixed WC and BMI, and how they line up with proWCdelta
## Output: T114_mri_checks.csv, T115_protein_vat_associations.csv
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
SD  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260924)
have_sandwich <- requireNamespace("sandwich", quietly = TRUE)

## ---------------- data ----------------
prot <- fread(file.path(SD, "source_table_Methods_protein_missingness_by_protein.csv"))$protein
stopifnot(length(prot) == 2920)
hdr <- names(fread(file.path(F12, "pro53013_新诊_newdead.csv"), nrows = 0))
prot <- intersect(prot, hdr); say("protein columns found: %d", length(prot))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"),
              select = c("Participant ID", "Sex", "Age", "WC", "BMI", "HBA1C", "HDL", "TRIG", "CRP", "ALT", "GGT", prot))
setnames(phen, 1:11, c("id", "Sex", "Age", "WC", "BMI", "hba1c", "hdl", "tg", "crp", "alt", "ggt"))
la  <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
bs  <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"),
             select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2"))
setnames(bs, c("id", "WC2", "BMI2"))
rec <- fread(file.path(RAW, "2Recruitment.csv"), select = c("Participant ID", "Date of attending assessment centre | Instance 0",
             "Date of attending assessment centre | Instance 2")); setnames(rec, c("id", "date0", "date2"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]

all <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI", "hba1c", "hdl", "tg", "crp", "alt", "ggt")) all[[v]] <- num(all[[v]])
all[, Sex := as.integer(Sex)]
all <- all[complete.cases(all[, .(Age, Sex, WC, BMI, dlt)])]
all <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(all, mri, bs, rec, cov))
for (v in c("vat", "asat", "WC2", "BMI2")) all[[v]] <- num(all[[v]])
sdz <- sd(all$dlt); all[, z := dlt / sdz]
all[, imaged := !is.na(vat) & !is.na(asat) & asat > 0]
d <- all[imaged == TRUE]
d[, years := as.numeric(as.IDate(date2) - as.IDate(date0)) / 365.25]
say("analysis set %d; imaging participants %d", nrow(all), nrow(d))

B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
rows <- list(); add <- function(check, spec, n, est, lo, hi) rows[[length(rows) + 1]] <<- data.table(check, spec, n, est, lo, hi)
fitz <- function(f, x, w = NULL, robust = FALSE) {
  m <- if (is.null(w)) lm(as.formula(f), x) else lm(as.formula(f), x, weights = w)
  b <- coef(m)[["z"]]
  se <- if (robust && have_sandwich) sqrt(sandwich::vcovHC(m, type = "HC3")["z", "z"]) else summary(m)$coefficients["z", 2]
  c(b = b, lo = b - 1.96 * se, hi = b + 1.96 * se) }

## ---------------- 1. linkage against a separate file ----------------
r2 <- cor(d$vat, d$WC2, use = "pair"); r0 <- cor(d$vat, d$WC, use = "pair")
say("[1] r(VAT, waist at the imaging visit) = %.3f (n = %d); r(VAT, baseline waist) = %.3f", r2, sum(!is.na(d$WC2)), r0)
add("1. linkage", "r(VAT, waist measured at the imaging visit)", sum(!is.na(d$WC2)), r2, NA, NA)
add("1. linkage", "r(VAT, waist measured at baseline)", nrow(d), r0, NA, NA)

## ---------------- 2. both visits fixed, all participants; time to the scan ----------------
x2 <- d[!is.na(WC2) & !is.na(BMI2)]
v <- fitz(paste("vat ~ z +", B0), x2);                                   add("2. timing", "baseline WC and BMI fixed (same participants)", nrow(x2), v[1], v[2], v[3])
v <- fitz(paste("vat ~ z +", B0, "+ ns(WC2,3) + ns(BMI2,3)"), x2);       add("2. timing", "WC and BMI of both visits fixed, all participants", nrow(x2), v[1], v[2], v[3])
say("[2] both visits fixed, all participants: %.3f (%.3f to %.3f), n = %d", v[1], v[2], v[3], nrow(x2))
x3 <- d[!is.na(years)]
v <- fitz(paste("vat ~ z +", B0, "+ years"), x3);                        add("2. timing", "time from blood sample to scan added", nrow(x3), v[1], v[2], v[3])

## ---------------- 3. visceral share within each sex; same form in both arms ----------------
share <- function(x, sexadj = TRUE, linear = FALSE) {
  sx <- if (sexadj) " + factor(Sex)" else ""
  pro <- if (linear) paste0("~ z + scale(WC) + ns(BMI,3) + Age", sx) else paste0("~ z + ns(WC,3) + ns(BMI,3) + Age", sx)
  bv <- coef(lm(as.formula(paste("vat", pro)), x))[["z"]]; bs <- coef(lm(as.formula(paste("asat", pro)), x))[["z"]]
  wv <- coef(lm(as.formula(paste0("vat ~ scale(WC) + ns(BMI,3) + Age", sx)), x))[[2]]
  ws <- coef(lm(as.formula(paste0("asat ~ scale(WC) + ns(BMI,3) + Age", sx)), x))[[2]]
  c(pro = bv / (bv + bs), wc = wv / (wv + ws), diff = bv / (bv + bs) - wv / (wv + ws)) }
for (spec in list(list("women", d[Sex == 0], FALSE, FALSE), list("men", d[Sex == 1], FALSE, FALSE),
                  list("all, linear WC term in both arms", d, TRUE, TRUE))) {
  x <- spec[[2]]; s0 <- share(x, spec[[3]], spec[[4]])
  bt <- do.call(rbind, mclapply(1:500, function(b) { set.seed(20260924 + b); share(x[sample(.N, .N, replace = TRUE)], spec[[3]], spec[[4]]) }, mc.cores = 20))
  q <- function(k) quantile(bt[, k], c(.025, .975))
  for (k in c("pro", "wc", "diff")) add("3. visceral share", paste0(spec[[1]], ": ", c(pro = "with proWCdelta", wc = "with a larger measured WC", diff = "difference")[[k]]),
                                         nrow(x), s0[[k]], q(k)[[1]], q(k)[[2]])
  say("[3] %s: %.1f%% vs %.1f%%, difference %.1f (%.1f to %.1f)", spec[[1]], 100 * s0[["pro"]], 100 * s0[["wc"]], 100 * s0[["diff"]], 100 * q("diff")[1], 100 * q("diff")[2]) }

## ---------------- 4. partial R2 ----------------
for (y in c("vat", "asat")) {
  m0 <- lm(as.formula(paste(y, "~", B0)), d); m1 <- lm(as.formula(paste(y, "~ z +", B0)), d)
  pr2 <- (sum(resid(m0)^2) - sum(resid(m1)^2)) / sum(resid(m0)^2)
  add("4. partial R2", paste0(y, ": partial R2 of proWCdelta at fixed WC, BMI, age and sex"), nrow(d), pr2, NA, NA)
  add("4. partial R2", paste0(y, ": R2 of age, sex, WC and BMI"), nrow(d), summary(m0)$r.squared, NA, NA)
  add("4. partial R2", paste0(y, ": R2 with proWCdelta added"), nrow(d), summary(m1)$r.squared, NA, NA)
  say("[4] %s: R2 %.3f -> %.3f; partial R2 of proWCdelta %.3f", y, summary(m0)$r.squared, summary(m1)$r.squared, pr2) }

## ---------------- 5. routine blood measurements held fixed ----------------
x5 <- d[complete.cases(d[, .(hba1c, hdl, tg, crp, alt, ggt)]) & tg > 0 & crp > 0 & alt > 0 & ggt > 0]
RB <- "ns(log(tg),3) + ns(hdl,3) + ns(hba1c,3) + ns(log(crp),3) + ns(log(alt),3) + ns(log(ggt),3)"
v <- fitz(paste("vat ~ z +", B0), x5);             add("5. routine blood", "VAT, primary model (same participants)", nrow(x5), v[1], v[2], v[3])
v5 <- fitz(paste("vat ~ z +", B0, "+", RB), x5);   add("5. routine blood", "VAT, + triglycerides, HDL, HbA1c, CRP, ALT, GGT", nrow(x5), v5[1], v5[2], v5[3])
say("[5] VAT per SD: %.3f in the same participants; %.3f (%.3f to %.3f) with routine blood measurements fixed (n = %d)", v[1], v5[1], v5[2], v5[3], nrow(x5))
mb <- lm(as.formula(paste("vat ~", B0, "+", RB)), x5); mbz <- lm(as.formula(paste("vat ~ z +", B0, "+", RB)), x5)
add("5. routine blood", "VAT: partial R2 of proWCdelta beyond routine blood measurements", nrow(x5),
    (sum(resid(mb)^2) - sum(resid(mbz)^2)) / sum(resid(mb)^2), NA, NA)
m0 <- lm(as.formula(paste("vat ~", B0)), x5)
add("5. routine blood", "VAT: partial R2 of routine blood measurements at fixed WC, BMI, age and sex", nrow(x5),
    (sum(resid(m0)^2) - sum(resid(mb)^2)) / sum(resid(m0)^2), NA, NA)

## ---------------- 6. attendance weighting; robust standard errors ----------------
a6 <- all[!is.na(tdi) & !is.na(smoking)]
pm <- glm(imaged ~ z + Age + factor(Sex) + ns(WC,3) + ns(BMI,3) + tdi + smoking, a6, family = binomial())
a6[, p := fitted(pm)]; a6[, w := mean(imaged) / p]
x6 <- a6[imaged == TRUE]
v <- fitz(paste("vat ~ z +", B0), x6, w = x6$w, robust = TRUE); add("6. selection", "VAT, weighted for attendance at imaging (robust SE)", nrow(x6), v[1], v[2], v[3])
say("[6] weighted for attendance: %.3f (%.3f to %.3f)%s", v[1], v[2], v[3], if (have_sandwich) "" else " [model SE; sandwich not installed]")
v <- fitz(paste("vat ~ z +", B0), d, robust = TRUE);             add("6. selection", "VAT, primary model with robust SE", nrow(d), v[1], v[2], v[3])

res <- rbindlist(rows); setnames(res, c("check", "specification", "n", "estimate", "lo", "hi"))
fwrite(res, file.path(TB, "T114_mri_checks.csv")); print(res, digits = 3)

## ---------------- 7. proteins and visceral fat ----------------
for (p in prot) set(d, j = p, value = num(d[[p]]))
pv <- rbindlist(mclapply(prot, function(p) { x <- d[!is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
  co <- summary(lm(as.formula(paste("vat ~ pz +", B0)), x))$coefficients["pz", ]
  data.table(protein = p, n_vat = nrow(x), beta_vat = co[1], se_vat = co[2], p_vat = co[4]) }, mc.cores = 20))
for (p in prot) set(all, j = p, value = num(all[[p]]))
pd <- rbindlist(mclapply(prot, function(p) { x <- all[!is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
  co <- summary(lm(as.formula(paste("pz ~ z +", B0)), x))$coefficients["z", ]
  data.table(protein = p, beta_delta = co[1], p_delta = co[4]) }, mc.cores = 20))
co <- fread(file.path(SD, "model_proWC_full_coefficients.csv")); setnames(co, 1:2, c("protein", "coef"))
nm <- fread(file.path(SD, "source_table_Methods_protein_names.csv")); setnames(nm, 1:2, c("protein", "name"))
pv <- Reduce(function(a, b) merge(a, b, by = "protein", all.x = TRUE), list(pv, pd, co[, .(protein, coef)], nm[, .(protein, name)]))
pv[is.na(coef), coef := 0]; pv[, fdr_vat := p.adjust(p_vat, "BH")]
setorder(pv, p_vat)
say("[7] proteins associated with VAT at fixed WC and BMI (FDR < 0.05): %d of %d; correlation across proteins of the",
    sum(pv$fdr_vat < 0.05), nrow(pv))
say("    association with VAT and the association with proWCdelta (both at fixed WC, BMI, age and sex): r = %.2f",
    cor(pv$beta_vat, pv$beta_delta))
say("    of the 20 proteins most strongly associated with VAT, %d are in the proWC model", sum(pv$coef[1:20] != 0))
print(pv[1:20, .(protein, beta_vat = round(beta_vat, 3), p_vat = signif(p_vat, 2), beta_delta = round(beta_delta, 3), coef = signif(coef, 3))])
fwrite(pv, file.path(TB, "T115_protein_vat_associations.csv"))
say("DONE")
