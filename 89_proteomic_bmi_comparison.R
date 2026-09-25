## 89_proteomic_bmi_comparison.R
## Is the waist the more informative training label? proWCdelta is compared with the discordance of a proteomic BMI
## built in exactly the same way (38a_oof_comparator_score.R: LASSO on the same 2,920 proteins with age and sex, nested
## ten-fold cross-validation, out-of-fold predictions, residual anchoring on measured BMI), per SD at fixed WC, BMI, age
## and sex (the variation left at fixed body size), alone and together, for abdominal fat on MRI (as in 74 and 75) and
## for incident disease in participants free of each endpoint (definitions as in 57). The script first reproduces the
## published estimates for proWCdelta and stops if any differs.
## Output: T122_proteomic_BMI_comparison.csv; source_table_SupplementaryTable7_proteomic_BMI.csv (formatted)
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
SD <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
res <- list(); add <- function(...) res[[length(res) + 1]] <<- data.table(...)

## ---------------- the two discordances ----------------
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI")); setnames(phen, c("id", "Sex", "Age", "WC", "BMI"))
la <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
bm <- fread(file.path(W, "output", "comparator_scores", "oof_BMI.csv"))[, .(id, bdl = delta)]      ## 38a_oof_comparator_score.R BMI
pb <- fread(file.path(W, "output", "comparator_scores", "perf_BMI.csv"))
coh <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id"); for (v in c("Age", "WC", "BMI")) coh[[v]] <- num(coh[[v]]); coh[, Sex := as.integer(Sex)]
coh <- coh[complete.cases(coh[, .(Age, Sex, WC, BMI, dlt)])]; sdz <- sd(coh$dlt)
coh <- merge(coh, bm, by = "id"); stopifnot(nrow(coh) == 52742, abs(sdz - 5.4703) < 0.001)
rz <- function(v) resid(lm(as.formula(paste(v, "~ ns(WC,3) + ns(BMI,3) + Age + Sex")), coh))
S <- c(dlt = sd(rz("dlt")), bdl = sd(rz("bdl"))); r_fixed <- cor(rz("dlt"), rz("bdl"))
r2_wc <- with(la, 1 - sum((Actual_WC - pred_WC)^2) / sum((Actual_WC - mean(Actual_WC))^2))
say("out-of-fold R2: WC %.4f (n %d), BMI %.4f (n %d); SD at fixed size: proWCdelta %.3f cm, proBMIdelta %.3f kg/m2; r at fixed size %.3f",
    r2_wc, nrow(la), pb$r2_oof, pb$n, S[["dlt"]], S[["bdl"]], r_fixed)
add(section = "accuracy", measure = "out-of-fold R2 for the measured trait", model = "", term = "proWCdelta", n = nrow(la), events = NA_integer_, estimate = r2_wc, lo = NA_real_, hi = NA_real_, p = NA_real_)
add(section = "accuracy", measure = "out-of-fold R2 for the measured trait", model = "", term = "proBMIdelta", n = pb$n, events = NA_integer_, estimate = pb$r2_oof, lo = NA_real_, hi = NA_real_, p = NA_real_)
add(section = "accuracy", measure = "SD at fixed WC, BMI, age and sex", model = "", term = "proWCdelta", n = nrow(coh), events = NA_integer_, estimate = S[["dlt"]], lo = NA_real_, hi = NA_real_, p = NA_real_)
add(section = "accuracy", measure = "SD at fixed WC, BMI, age and sex", model = "", term = "proBMIdelta", n = nrow(coh), events = NA_integer_, estimate = S[["bdl"]], lo = NA_real_, hi = NA_real_, p = NA_real_)
add(section = "accuracy", measure = "correlation of the two discordances at fixed WC, BMI, age and sex", model = "", term = "both", n = nrow(coh), events = NA_integer_, estimate = r_fixed, lo = NA_real_, hi = NA_real_, p = NA_real_)
LAB <- c(c_dlt = "proWCdelta", c_bdl = "proBMIdelta")

## ---------------- abdominal MRI, as in 74_abdominal_mri.R ----------------
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2", "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
im <- Reduce(function(a, b) merge(a, b, by = "id"), list(phen, la[, .(id, dlt = BioX_Delta)], mri))
for (v in c("Age", "WC", "BMI", "vat", "asat")) im[[v]] <- num(im[[v]]); im[, Sex := as.integer(Sex)]
im <- im[!is.na(dlt) & !is.na(WC) & !is.na(BMI) & !is.na(Age) & !is.na(Sex) & !is.na(vat) & !is.na(asat) & asat > 0]
im <- merge(im, bm, by = "id"); stopifnot(nrow(im) == 6901, cor(im$vat, im$WC) > 0.7)
im[, z := dlt / sdz]; im[, `:=`(c_dlt = dlt / S[["dlt"]], c_bdl = bdl / S[["bdl"]])]
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
b_rep <- coef(lm(as.formula(paste("vat ~ z +", B0)), im))[["z"]]
say("REPRODUCTION: visceral fat per SD of proWCdelta %.6f (published 0.545079)", b_rep); stopifnot(abs(b_rep - 0.545078589469762) < 1e-8)
lmfit <- function(y, terms, dat) { co <- summary(lm(as.formula(paste(y, "~", paste(terms, collapse = " + "), "+", B0)), dat))$coefficients
  lapply(terms, function(t) c(b = co[t, 1], se = co[t, 2], p = co[t, 4])) }
for (y in c("vat", "asat")) for (mdl in list("c_dlt", "c_bdl", c("c_dlt", "c_bdl"))) { f <- lmfit(y, mdl, im)
  for (j in seq_along(mdl)) add(section = "abdominal MRI, per SD at fixed WC, BMI, age and sex", measure = ifelse(y == "vat", "visceral adipose tissue (L)", "abdominal subcutaneous adipose tissue (L)"),
    model = ifelse(length(mdl) == 1, "alone", "both in one model"), term = LAB[[mdl[j]]], n = nrow(im), events = NA_integer_,
    estimate = f[[j]][["b"]], lo = f[[j]][["b"]] - 1.96 * f[[j]][["se"]], hi = f[[j]][["b"]] + 1.96 * f[[j]][["se"]], p = f[[j]][["p"]]) }
r2b <- summary(lm(as.formula(paste("vat ~", B0)), im))$r.squared
for (t in names(LAB)) { r2 <- summary(lm(as.formula(paste("vat ~", t, "+", B0)), im))$r.squared
  add(section = "abdominal MRI, per SD at fixed WC, BMI, age and sex", measure = "partial R2 for visceral adipose tissue", model = "alone", term = LAB[[t]], n = nrow(im), events = NA_integer_, estimate = (r2 - r2b) / (1 - r2b), lo = NA_real_, hi = NA_real_, p = NA_real_) }
share <- function(t, dat) { bv <- coef(lm(as.formula(paste("vat ~", t, "+", B0)), dat))[[t]]; ba <- coef(lm(as.formula(paste("asat ~", t, "+", B0)), dat))[[t]]; bv / (bv + ba) }
set.seed(2026); BS <- replicate(500, { ii <- sample.int(nrow(im), replace = TRUE); sapply(names(LAB), function(t) share(t, im[ii])) })
for (t in names(LAB)) add(section = "abdominal MRI, per SD at fixed WC, BMI, age and sex", measure = "visceral share of the extra abdominal fat", model = "alone", term = LAB[[t]],
  n = nrow(im), events = NA_integer_, estimate = share(t, im), lo = quantile(BS[t, ], 0.025), hi = quantile(BS[t, ], 0.975), p = NA_real_)

## ---------------- incident disease, as in 57_primary_incident_table2.R ----------------
phen2 <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen2, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen2, la[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id"), cov, fr))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]; stopifnot(nrow(d) == 52742)
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
glucose_drugs <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, glucose_drugs) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I10 = has(d$selfrep, "(^|\\|)(hypertension|essential hypertension)(\\||$)") | has(d$med, "blood pressure medication"),
  E78 = has(d$selfrep, "(^|\\|)high cholesterol(\\||$)") | has(d$med, "cholesterol lowering medication"),
  I25 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60), E66 = rep(FALSE, nrow(d)))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease", "Ischaemic heart disease", "Dyslipidaemia", "Hypertension", "Obesity diagnosis"),
                  code = c("E11", "I50", "K76", "N18", "I25", "E78", "I10", "E66"),
                  published = c(1.82105304412415, 1.50033695791199, 1.43514603957416, 1.29353725419005, NA, NA, NA, NA))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", cd); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
d <- merge(d, bm, by = "id"); stopifnot(nrow(d) == 52742)
d[, `:=`(c_dlt = dlt / S[["dlt"]], c_bdl = bdl / S[["bdl"]])]
gfit <- function(y, terms, dat) { m <- glm(as.formula(paste(y, "~", paste(terms, collapse = " + "), "+", BASE)), dat, family = binomial()); co <- summary(m)$coefficients
  list(m = m, co = lapply(terms, function(t) c(b = co[t, 1], se = co[t, 2], p = co[t, 4]))) }
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[get(paste0("free_", cd)) == TRUE]; ev <- sum(x[[cd]])
  rep <- exp(gfit(cd, "z", x)$co[[1]][["b"]])
  if (!is.na(dis$published[i])) { say("REPRODUCTION: %s n = %d, events = %d, OR per SD %.6f (published %.6f)", dis$disease[i], nrow(x), ev, rep, dis$published[i]); stopifnot(abs(rep - dis$published[i]) < 1e-8) }
  fits <- list(a = gfit(cd, "c_dlt", x), b = gfit(cd, "c_bdl", x), ab = gfit(cd, c("c_dlt", "c_bdl"), x))
  for (k in c("a", "b", "ab")) { tt <- switch(k, a = "c_dlt", b = "c_bdl", ab = c("c_dlt", "c_bdl"))
    for (j in seq_along(tt)) { cc <- fits[[k]]$co[[j]]
      add(section = "incident disease, OR per SD at fixed WC, BMI, age and sex", measure = dis$disease[i], model = ifelse(k == "ab", "both in one model", "alone"), term = LAB[[tt[j]]],
          n = nrow(x), events = ev, estimate = exp(cc[["b"]]), lo = exp(cc[["b"]] - 1.96 * cc[["se"]]), hi = exp(cc[["b"]] + 1.96 * cc[["se"]]), p = cc[["p"]]) } }
  for (k in c("a", "b")) { base <- fits[[ifelse(k == "a", "b", "a")]]$m; chi <- as.numeric(2 * (logLik(fits$ab$m) - logLik(base)))
    add(section = "incident disease, likelihood-ratio chi-square (1 df) for adding one discordance to the other", measure = dis$disease[i], model = "both in one model",
        term = LAB[[ifelse(k == "a", "c_dlt", "c_bdl")]], n = nrow(x), events = ev, estimate = chi, lo = NA_real_, hi = NA_real_, p = pchisq(chi, 1, lower.tail = FALSE)) } }

tab <- rbindlist(res); fwrite(tab, file.path(TB, "T122_proteomic_BMI_comparison.csv"))
options(width = 250); print(tab[, .(section = substr(section, 1, 40), measure, model, term, n, events, est = signif(estimate, 4), lo = signif(lo, 4), hi = signif(hi, 4), p = signif(p, 2))], nrows = 200)

## ---------------- Supplementary Table 7 (formatted; section rows have empty columns) ----------------
mn <- function(x) gsub("-", "−", x)
ci <- function(e, l, h, f = "%.2f") mn(sprintf(paste0(f, " (", f, " to ", f, ")"), e, l, h))
g <- function(sec, mea, mod, trm) tab[section == sec & measure == mea & model == mod & term == trm]
cm <- function(x) formatC(x, format = "d", big.mark = ",")
rows <- list(c("Accuracy and spread", "", "", ""),
  c("Out-of-fold R² of the proteomic prediction for the measured trait", sprintf("%.3f", g("accuracy", "out-of-fold R2 for the measured trait", "", "proWCdelta")$estimate),
    sprintf("%.3f", g("accuracy", "out-of-fold R2 for the measured trait", "", "proBMIdelta")$estimate), paste0(cm(nrow(la)), "; ", cm(pb$n))),
  c("SD at fixed WC, BMI, age and sex", sprintf("%.2f cm", S[["dlt"]]), sprintf("%.2f kg/m²", S[["bdl"]]), cm(nrow(coh))))
sec_mri <- "abdominal MRI, per SD at fixed WC, BMI, age and sex"
rows <- c(rows, list(c("Abdominal MRI, per SD at fixed WC, BMI, age and sex", "", "", "")))
for (mea in c("visceral adipose tissue (L)", "abdominal subcutaneous adipose tissue (L)")) { a <- g(sec_mri, mea, "alone", "proWCdelta"); b <- g(sec_mri, mea, "alone", "proBMIdelta")
  rows <- c(rows, list(c(paste0(toupper(substr(mea, 1, 1)), substr(mea, 2, nchar(mea))), ci(a$estimate, a$lo, a$hi), ci(b$estimate, b$lo, b$hi), cm(a$n)))) }
a <- g(sec_mri, "visceral share of the extra abdominal fat", "alone", "proWCdelta"); b <- g(sec_mri, "visceral share of the extra abdominal fat", "alone", "proBMIdelta")
rows <- c(rows, list(c("Visceral share of the extra abdominal fat (%)", ci(100 * a$estimate, 100 * a$lo, 100 * a$hi, "%.1f"), ci(100 * b$estimate, 100 * b$lo, 100 * b$hi, "%.1f"), cm(a$n))))
a <- g(sec_mri, "partial R2 for visceral adipose tissue", "alone", "proWCdelta"); b <- g(sec_mri, "partial R2 for visceral adipose tissue", "alone", "proBMIdelta")
rows <- c(rows, list(c("Partial R² for visceral adipose tissue (%)", sprintf("%.1f", 100 * a$estimate), sprintf("%.1f", 100 * b$estimate), cm(a$n))))
a <- g(sec_mri, "visceral adipose tissue (L)", "both in one model", "proWCdelta"); b <- g(sec_mri, "visceral adipose tissue (L)", "both in one model", "proBMIdelta")
rows <- c(rows, list(c("Visceral adipose tissue (L), both in one model", ci(a$estimate, a$lo, a$hi), ci(b$estimate, b$lo, b$hi), cm(a$n))))
sec_d <- "incident disease, OR per SD at fixed WC, BMI, age and sex"; sec_lr <- "incident disease, likelihood-ratio chi-square (1 df) for adding one discordance to the other"
for (blk in list(list("Incident disease, odds ratio per SD at fixed WC, BMI, age and sex", "alone"), list("Incident disease, odds ratio per SD, both in one model", "both in one model"))) {
  rows <- c(rows, list(c(blk[[1]], "", "", "")))
  for (dz in dis$disease) { a <- g(sec_d, dz, blk[[2]], "proWCdelta"); b <- g(sec_d, dz, blk[[2]], "proBMIdelta")
    rows <- c(rows, list(c(dz, ci(a$estimate, a$lo, a$hi), ci(b$estimate, b$lo, b$hi), paste0(cm(a$n), " (", cm(a$events), ")")))) } }
rows <- c(rows, list(c("Incident disease, likelihood-ratio χ² (1 df) for adding each to the other", "", "", "")))
for (dz in dis$disease) { a <- g(sec_lr, dz, "both in one model", "proWCdelta"); b <- g(sec_lr, dz, "both in one model", "proBMIdelta")
  rows <- c(rows, list(c(dz, sprintf("%.1f", a$estimate), sprintf("%.1f", b$estimate), paste0(cm(a$n), " (", cm(a$events), ")")))) }
st7 <- as.data.table(do.call(rbind, rows)); setnames(st7, c("Measure", "proWCΔ", "Proteomic BMI discordance", "Participants (events)"))
fwrite(st7, file.path(SD, "source_table_SupplementaryTable7_proteomic_BMI.csv")); print(st7)
say("DONE")
