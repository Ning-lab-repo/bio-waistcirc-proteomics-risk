## 76_standardised_risk_ci.R
## The standardised 10-year risks at the 10th and 90th percentiles of proWCdelta, conditional on waist, BMI, age and
## sex, were reported by 57_primary_incident_table2.R (T81) without a measure of uncertainty. This script repeats that
## computation on the same population, with the same model, and adds bootstrap confidence intervals for each risk, for
## their difference and for their ratio. The residualisation that defines the conditional percentiles is refitted
## inside every resample, so the intervals cover that step too. The point estimates are checked against T81.
## Output: T112_standardised_risk_with_ci.csv
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

## ---------------- the population and endpoints of 57_primary_incident_table2.R ----------------
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis",
  "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name",
  "Medication for cholesterol, blood pressure or diabetes | Instance 0",
  "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, fr))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep),
         med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) *
            0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
gd <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") |
                    has(d$medname, gd) | (!is.na(d$hba1c) & d$hba1c >= 48),
              I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"),
              K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
              N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") |
                    (!is.na(d$egfr) & d$egfr < 60))
NAMES <- c(E11 = "Type 2 diabetes", I50 = "Heart failure", K76 = "Liver disease", N18 = "Chronic kidney disease")
for (cd in names(NAMES)) { pt <- paste0("(^|[|])", cd); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
say("modelling set n = %d", nrow(d))

## ---------------- standardised risks at the conditional 10th and 90th percentiles ----------------
std_risk <- function(dd) {
  dd <- copy(dd); dd[, e := resid(lm(dlt ~ ns(WC, 3) + ns(BMI, 3) + Age + Sex, dd))]; dd[, mu := dlt - e]
  qc <- quantile(dd$e, c(0.1, 0.9)); sz <- sd(dd$dlt)
  unlist(lapply(names(NAMES), function(cd) { x <- dd[get(paste0("free_", cd)) == TRUE]
    mm <- glm(as.formula(paste(cd, "~ z +", BASE)), x, family = binomial())
    r <- sapply(qc, function(q) { nd <- copy(x); nd[, z := (mu + q) / sdd]; mean(predict(mm, nd, type = "response")) })
    setNames(c(r[[1]], r[[2]]), paste0(cd, c("_p10", "_p90"))) })) }

obs <- std_risk(d)
t81 <- fread(file.path(TB, "T81_conditional_percentile_risks.csv"))
for (cd in names(NAMES)) { r <- t81[disease == cd]
  say("%s: %.4f%% / %.4f%% here, %.4f%% / %.4f%% in T81", cd, 100 * obs[[paste0(cd, "_p10")]], 100 * obs[[paste0(cd, "_p90")]],
      100 * r$risk_p10, 100 * r$risk_p90) }
stopifnot(all(abs(obs[paste0(t81$disease, "_p10")] - t81$risk_p10) < 1e-6), all(abs(obs[paste0(t81$disease, "_p90")] - t81$risk_p90) < 1e-6))
say("point estimates reproduce T81")

B <- 500
say("bootstrap, %d resamples", B)
bt <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(20260923 + b); std_risk(d[sample.int(nrow(d), replace = TRUE)]) },
                              mc.cores = 20))
q <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
res <- rbindlist(lapply(names(NAMES), function(cd) {
  a <- paste0(cd, "_p10"); b <- paste0(cd, "_p90")
  data.table(code = cd, disease = NAMES[[cd]], n = sum(d[[paste0("free_", cd)]]), events = sum(d[get(paste0("free_", cd)) == TRUE][[cd]]),
             risk_p10 = obs[[a]], lo_p10 = q(bt[, a])[[1]], hi_p10 = q(bt[, a])[[2]],
             risk_p90 = obs[[b]], lo_p90 = q(bt[, b])[[1]], hi_p90 = q(bt[, b])[[2]],
             difference = obs[[b]] - obs[[a]], lo_diff = q(bt[, b] - bt[, a])[[1]], hi_diff = q(bt[, b] - bt[, a])[[2]],
             ratio = obs[[b]] / obs[[a]], lo_ratio = q(bt[, b] / bt[, a])[[1]], hi_ratio = q(bt[, b] / bt[, a])[[2]]) }))
for (i in seq_len(nrow(res))) with(res[i], say("%-24s %.1f%% (%.1f-%.1f) vs %.1f%% (%.1f-%.1f); difference %.1f points (%.1f-%.1f); ratio %.2f (%.2f-%.2f)",
  disease, 100 * risk_p10, 100 * lo_p10, 100 * hi_p10, 100 * risk_p90, 100 * lo_p90, 100 * hi_p90,
  100 * difference, 100 * lo_diff, 100 * hi_diff, ratio, lo_ratio, hi_ratio))
fwrite(res, file.path(TB, "T112_standardised_risk_with_ci.csv"))
say("DONE")
