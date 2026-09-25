## 44_confounding_technical_holdout.R
## Further robustness of the continuous proWCdelta associations (per SD, at fixed WC and BMI):
##  (1) extended adjustment: ethnicity, education, alcohol, physical activity and television time; technical factors
##      (Olink processing batch, sample storage time, fasting time, season, assessment centre); medication
##      (lipid-lowering, antihypertensive, insulin, hormone replacement therapy)
##  (2) outcomes without an expected adiposity pathway (injury chapter S00-T98 and external causes V01-Y98)
##  (3) E-values for the base-model estimates
##  (4) contribution of individual proteins: proWCdelta adjusted for ADM, IGFBP1, IGFBP2, cystatin C (CST3) and leptin,
##      and each of the 20 proteins of the reduced score modelled alone at fixed WC and BMI
##  (5) protein missingness: participants with <=10% and <=20% missing proteins, and adjustment for the missing fraction
##  (6) geographic hold-out: proWCdelta from the model trained in England, associations in Scotland and Wales
##  (7) height, weight and impedance in the WC- and BMI-matched pairs
## Output: T43-T48.
suppressPackageStartupMessages({ library(data.table); library(splines); library(MatchIt) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"; F2 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure2"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x))
wins <- function(x) { q <- quantile(x, c(.01, .99), na.rm = TRUE); pmin(pmax(x, q[1]), q[2]) }

lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "UK Biobank assessment centre | Instance 0"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "centre"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
tech <- fread(bf, select = c(hdr[1], "Batch", "sample_age_days", "Fasting_time", "season_binary", "ethnicity", "Alcohol intake frequency", "Education", "MET Physical activity", "Time_TV"))
setnames(tech, c("id", "batch", "storage", "fasting", "season", "white", "alcohol", "education", "met", "tv"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "med_m", "med_f")); fr[, med := paste(med_m, med_f)]
has <- function(k) as.integer(grepl(paste0("(^|[| ])", k, "([| ]|$)"), fr$med))
fr[, `:=`(lipid_med = has(1), bp_med = has(2), insulin = has(3), hrt = has(4))]

d <- merge(phen, lasso[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id")
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(d, cov, tech, fr[, .(id, lipid_med, bp_med, insulin, hrt)]))
for (v in c("Age", "WC", "BMI", "storage", "fasting", "met", "tv")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
sdd <- sd(d$dlt); d[, z := dlt / sdd]
cat("cohort n", nrow(d), " SD(proWCdelta)", round(sdd, 3), "\n")
x <- as.character(d$dx10); x[is.na(x)] <- ""
dis <- data.table(disease = c("Type 2 diabetes", "Obesity", "Dyslipidemia", "Hypertension", "Ischaemic heart disease", "Heart failure", "Liver disease", "Chronic kidney disease"),
                  code = c("E11", "E66", "E78", "I10", "I25", "I50", "K76", "N18"))
for (i in seq_len(nrow(dis))) d[[dis$code[i]]] <- as.integer(grepl(paste0("(^|[|])", dis$code[i]), x))
d[, INJ := as.integer(grepl("(^|[|])[ST]", x))]; d[, EXT := as.integer(grepl("(^|[|])[VWXY]", x))]
neg <- data.table(disease = c("Injury, poisoning (S00-T98)", "External causes (V01-Y98)"), code = c("INJ", "EXT"))

## covariate coding; education and physical activity carry a missing-value category
d[, educ_c := factor(fifelse(is.na(education), "missing", as.character(education)))]
d[, met_c := factor(fifelse(is.na(met), "missing", as.character(cut(met, quantile(met, 0:3 / 3, na.rm = TRUE), include.lowest = TRUE, labels = c("low", "mid", "high")))))]
d[, `:=`(alcohol = factor(alcohol), white = factor(white), batch = factor(batch), season = factor(season), centre = factor(centre),
         storage_w = wins(storage), fasting_w = wins(fasting), tv_w = wins(tv))]
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
SOC  <- "white + educ_c + alcohol + met_c + ns(tv_w,3)"
TECH <- "batch + ns(storage_w,3) + ns(fasting_w,3) + season + centre"
MED  <- "lipid_med + bp_med + insulin + hrt"
need <- c("Age", "Sex", "tdi", "smoking", "WC", "BMI", "alcohol", "white", "tv_w", "batch", "storage_w", "fasting_w", "season", "centre", "lipid_med", "bp_med", "insulin", "hrt")
s <- d[complete.cases(d[, ..need])]
cat("extended-adjustment sample n", nrow(s), "\n")
models <- list("Base (age, sex, deprivation, smoking, WC and BMI splines)" = BASE,
               "+ ethnicity, education, alcohol, physical activity, television" = paste(BASE, "+", SOC),
               "+ batch, storage time, fasting time, season, centre" = paste(BASE, "+", TECH),
               "+ medication" = paste(BASE, "+", MED),
               "All of the above" = paste(BASE, "+", SOC, "+", TECH, "+", MED))
fitz <- function(y, f, dat) { m <- glm(as.formula(paste(y, "~ z +", f)), dat, family = binomial()); co <- summary(m)$coefficients["z", ]
  c(OR = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) }
out <- list()
for (i in seq_len(nrow(dis))) for (mn in names(models)) { r <- fitz(dis$code[i], models[[mn]], s)
  out[[length(out) + 1]] <- data.table(analysis = "extended adjustment", disease = dis$disease[i], model = mn, n = nrow(s), events = sum(s[[dis$code[i]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
cc <- s[educ_c != "missing" & met_c != "missing"]
for (i in seq_len(nrow(dis))) { r <- fitz(dis$code[i], models[["All of the above"]], cc)
  out[[length(out) + 1]] <- data.table(analysis = "extended adjustment", disease = dis$disease[i], model = "All of the above, complete cases for education and physical activity", n = nrow(cc), events = sum(cc[[dis$code[i]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
## outcomes without an expected adiposity pathway
for (i in seq_len(nrow(neg))) for (mn in names(models)[c(1, 5)]) { r <- fitz(neg$code[i], models[[mn]], s)
  out[[length(out) + 1]] <- data.table(analysis = "outcomes without an expected adiposity pathway", disease = neg$disease[i], model = mn, n = nrow(s), events = sum(s[[neg$code[i]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
t43 <- rbindlist(out)
## E-values (VanderWeele and Ding): OR taken as RR when the 10-year risk is below 15%, otherwise square-root transformed
ev <- function(rr) ifelse(rr <= 1, 1, rr + sqrt(rr * (rr - 1)))
t43[, risk := events / n]; t43[, `:=`(RR = fifelse(risk < .15, OR, sqrt(OR)), RR_lo = fifelse(risk < .15, lo, sqrt(lo)))]
t43[, `:=`(E_value = ev(RR), E_value_CI = ev(RR_lo))]
t43[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]
print(dcast(t43, disease ~ model, value.var = "txt"), width = 300)
print(t43[model == names(models)[1], .(disease, events, risk = round(risk, 3), E_value = round(E_value, 2), E_value_CI = round(E_value_CI, 2))])
fwrite(t43[, !c("txt")], file.path(TB, "T43_extended_adjustment_negative_outcomes_evalues.csv"))

## (4) protein components
pd <- fread(file.path(F12, "analysis_data_WC.csv"), select = c(1, which(names(fread(file.path(F12, "analysis_data_WC.csv"), nrows = 0)) %in% c("LEP", "FABP4", "NCAN", "IL1RN", "WFIKKN2", "PON3", "SLITRK1", "IGFBP2", "IGSF9", "ADM", "IGFBP1", "SSC4D", "CHGB", "FURIN", "RTN4R", "CFH", "IGSF3", "OPTC", "SEZ6L", "PRAP1", "CST3"))))
setnames(pd, 1, "id"); pv <- setdiff(names(pd), "id")
b <- merge(d, pd, by = "id"); b <- b[complete.cases(b[, .(Age, Sex, tdi, smoking, WC, BMI)])]
for (v in pv) b[[paste0("s_", v)]] <- as.numeric(scale(b[[v]]))
adj <- list("none" = "", "ADM" = "s_ADM", "IGFBP1 and IGFBP2" = "s_IGFBP1 + s_IGFBP2", "cystatin C (CST3)" = "s_CST3", "leptin" = "s_LEP",
            "ADM, IGFBP1, IGFBP2, CST3 and leptin" = "s_ADM + s_IGFBP1 + s_IGFBP2 + s_CST3 + s_LEP")
o4 <- list()
for (cd in c("E11", "I50", "N18")) for (an in names(adj)) { f <- if (adj[[an]] == "") BASE else paste(BASE, "+", adj[[an]]); r <- fitz(cd, f, b)
  o4[[length(o4) + 1]] <- data.table(analysis = "proWCdelta adjusted for single proteins", disease = dis[code == cd]$disease, term = paste("adjusted for", an), n = nrow(b), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) }
sp <- setdiff(pv, "CST3")
for (cd in c("E11", "I50", "N18")) for (v in sp) { m <- glm(as.formula(paste(cd, "~", paste0("s_", v), "+", BASE)), b, family = binomial()); co <- summary(m)$coefficients[paste0("s_", v), ]
  o4[[length(o4) + 1]] <- data.table(analysis = "reduced-score protein alone at fixed WC and BMI (per SD NPX)", disease = dis[code == cd]$disease, term = v, n = nrow(b), OR = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) }
t45 <- rbindlist(o4); t45[analysis != "proWCdelta adjusted for single proteins", p_BH := p.adjust(p, "BH")]
print(t45[analysis == "proWCdelta adjusted for single proteins", .(disease, term, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))])
print(dcast(t45[analysis != "proWCdelta adjusted for single proteins", .(disease, term, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], term ~ disease, value.var = "txt"))
fwrite(t45, file.path(TB, "T45_protein_components.csv"))

## (5) protein missingness, continuous models
pm <- fread(bf, select = c(hdr[1], pcols)); setnames(pm, 1, "id")
pm <- data.table(id = pm$id, pna = rowMeans(is.na(as.matrix(pm[, -1])))); rm(list = setdiff(ls(), c("pm", "d", "dis", "BASE", "fitz", "num", "sdd", "TB", "F2", "F12", "F7", "W", "wins"))); gc()
m5 <- merge(d, pm, by = "id"); m5 <- m5[complete.cases(m5[, .(Age, Sex, tdi, smoking, WC, BMI)])]
cat("missing fraction: median", round(median(m5$pna), 3), "; >10%:", sum(m5$pna > .1), "; >20%:", sum(m5$pna > .2), "\n")
o5 <- list()
sets <- list("all participants" = m5, "<=20% proteins missing" = m5[pna <= .2], "<=10% proteins missing" = m5[pna <= .1])
for (i in seq_len(nrow(dis))) { for (sn in names(sets)) { r <- fitz(dis$code[i], BASE, sets[[sn]])
    o5[[length(o5) + 1]] <- data.table(subset = sn, disease = dis$disease[i], n = nrow(sets[[sn]]), events = sum(sets[[sn]][[dis$code[i]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
  r <- fitz(dis$code[i], paste(BASE, "+ ns(pna,3)"), m5)
  o5[[length(o5) + 1]] <- data.table(subset = "all participants, adjusted for missing fraction", disease = dis$disease[i], n = nrow(m5), events = sum(m5[[dis$code[i]]]), OR = r[["OR"]], lo = r[["lo"]], hi = r[["hi"]]) }
t46 <- rbindlist(o5); print(dcast(t46[, .(disease, subset, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease ~ subset, value.var = "txt"), width = 250)
fwrite(t46, file.path(TB, "T46_protein_missingness_continuous.csv"))

## (6) geographic hold-out: model trained in England, applied to Scotland and Wales
ho <- fread(file.path(F2, "lasso_test.csv")); setnames(ho, 1, "id")
h <- merge(d[, !c("z", "pWC")], ho[, .(id, pWC_ho = pred_WC, dlt_ho = BioX_Delta)], by = "id")
h <- h[complete.cases(h[, .(Age, Sex, tdi, smoking, WC, BMI)])]
h[, `:=`(z_ho = dlt_ho / sd(dlt_ho), z_main = dlt / sd(dlt))]
cat("hold-out n", nrow(h), "; r(hold-out proWCdelta, cross-validated proWCdelta) =", round(cor(h$dlt_ho, h$dlt), 3), "; hold-out R2", round(cor(h$WC, h$pWC_ho)^2, 3), "\n")
o6 <- list()
for (i in seq_len(nrow(dis))) for (zz in c("z_ho", "z_main")) { m <- glm(as.formula(paste(dis$code[i], "~", zz, "+", BASE)), h, family = binomial()); co <- summary(m)$coefficients[zz, ]
  o6[[length(o6) + 1]] <- data.table(score = c(z_ho = "trained in England, applied to hold-out", z_main = "cross-validated score (main analysis), same participants")[[zz]], disease = dis$disease[i], n = nrow(h), events = sum(h[[dis$code[i]]]),
                                     OR = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) }
t47 <- rbindlist(o6); print(dcast(t47[, .(disease, events, score, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease + events ~ score, value.var = "txt"), width = 250)
fwrite(t47, file.path(TB, "T47_geographic_holdout_associations.csv"))

## (7) matched pairs (same procedure and seed as 08_bodycomp_anchoring.R): height, weight and impedance
BSZ <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/6-1Body_size_measures_participant.csv"
ht <- fread(BSZ, select = c("Participant ID", "Standing height | Instance 0")); setnames(ht, c("id", "height0"))
bc <- fread("/home/data/heamei/nmrLR/ukbnmr_met_rawdata/bodycomp_subset.csv", select = c("Participant ID", "Impedance of whole body | Instance 0", "Whole body fat mass | Instance 0", "Whole body fat-free mass | Instance 0"))
setnames(bc, c("id", "imp", "fm", "ffm"))
lz <- fread(file.path(F12, "lasso_WC.csv")); setnames(lz, 1, "id")
mm0 <- merge(merge(d[, .(participant_id = id, Sex, Age, WC, BMI, tdi, smoking)], lz[, .(participant_id = id, proWC = BioX_Adjusted)], by = "participant_id"), ht[, .(participant_id = id, height0 = num(height0))], by = "participant_id", all.x = TRUE)
mm0[, thr := fifelse(Sex == 0, 88, 102)]
mm0[, grp := fifelse(WC < thr & proWC < thr, "N/N", fifelse(WC < thr & proWC >= thr, "N/H", fifelse(WC >= thr & proWC < thr, "H/N", "H/H")))]
m <- mm0[grp %in% c("N/N", "N/H") & complete.cases(Age, Sex, WC, BMI, tdi, smoking)]; m[, NH := as.integer(grp == "N/H")]
set.seed(123)
mt <- matchit(NH ~ Age + WC + BMI, data = m, method = "nearest", exact = ~Sex, distance = "mahalanobis", caliper = c(WC = 0.15, BMI = 0.15), std.caliper = TRUE, ratio = 1)
md <- as.data.table(match.data(mt)); md <- merge(md, bc[, .(participant_id = id, imp = num(imp), fm = num(fm), ffm = num(ffm))], by = "participant_id", all.x = TRUE)
md[, `:=`(weight_bmi = BMI * (height0 / 100)^2, weight_bia = fm + ffm)]
ks <- readRDS(file.path(W, "output", "bodycomp_merge_keys.rds"))$matched_ids
cat("matched pairs", sum(md$NH == 1), "; identical to the published matching:", setequal(md$participant_id, ks), "\n")
o7 <- list()
for (v in c("Age", "WC", "BMI", "height0", "weight_bmi", "weight_bia", "fm", "ffm", "imp")) { dd <- md[!is.na(get(v))]
  pr <- dd[, .(y = get(v), NH, subclass)]; pw <- dcast(pr, subclass ~ NH, value.var = "y"); setnames(pw, c("subclass", "NN", "NH")); pw <- pw[complete.cases(pw)]
  tt <- t.test(pw$NH, pw$NN, paired = TRUE)
  o7[[length(o7) + 1]] <- data.table(variable = v, n_NH = sum(dd$NH == 1), n_NN = sum(dd$NH == 0), mean_NN = mean(dd[NH == 0][[v]]), mean_NH = mean(dd[NH == 1][[v]]),
                                     complete_pairs = nrow(pw), paired_diff = unname(tt$estimate), lo = tt$conf.int[1], hi = tt$conf.int[2], p = tt$p.value) }
t48 <- rbindlist(o7); print(t48)
fwrite(t48, file.path(TB, "T48_matched_pairs_height_weight_impedance.csv"))
cat("DONE\n")
