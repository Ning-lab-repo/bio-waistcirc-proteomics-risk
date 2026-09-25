## G5_gnhs_continuous_waist.R
## The one place in the published Guangzhou Nutrition and Health Study supplement where the waist is recorded in
## centimetres rather than as the sex-specific threshold flag: table S7, sheet E (Cai et al., Cell Reports Medicine
## 2023;4:101172). It holds 485 participants of the discovery cohort who had exactly one metabolic-syndrome component
## and no metabolic syndrome at baseline, with 438 serum proteins, age, sex, BMI and the measured waist, and table S7
## sheet A records which of them developed the syndrome by the second or third visit. That makes it possible to repeat
## the design of the present manuscript in full: train a proteomic waist on the measured waist, take its discordance
## from that measurement, and ask whether the discordance separates incident disease among people of the same measured
## size. It is small (88 events), so it tests whether the method reproduces, not how large the effect is.
## The same participants also carry the binary-flag score built in G4, so the two can be compared directly, which is
## what tells us whether the flag-based analysis in the larger cohorts distorts the answer.
## Output: G5_*.csv in derived/.
suppressPackageStartupMessages({ library(data.table); library(readxl); library(glmnet); library(splines) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

## ---------------- the analysis set ----------------
E <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "E", na = c("", "NA", "NaN")))
A <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "A", na = c("", "NA", "NaN")))
setnames(A, c("Patient_ID", "time", "Indicator", "MetS_BS", "MetS_F23"))
CLIN <- c("age","BMI","DBP","SBP","fbg","Glu","HDL","hdl","LDL","TG","tg","wc","waist","sex","sbpdbp","mets","metsscore",
          "TC","time","pat_ID","Sample_ID","DM_med","dyslipi_med","hyper_med","dm_cl","hyper_cl","tg_cl","tc_cl",
          "hdl_cl","ldl_cl","dys_cl","menopau_age")
pv <- setdiff(names(E), CLIN)
d <- merge(A[, .(Patient_ID, Indicator, ev = as.integer(MetS_F23 == 1))],
           E[, c(list(Patient_ID = pat_ID, wc = as.numeric(wc), BMI = as.numeric(BMI), age = as.numeric(age),
                      sex = as.integer(sex), waist = as.integer(waist)), lapply(.SD, as.numeric)), .SDcols = pv],
           by = "Patient_ID")
say("analysis set: %d participants, %d proteins, %d incident metabolic syndrome (%.1f%%)", nrow(d), length(pv), sum(d$ev), 100 * mean(d$ev))
say("measured waist %.1f cm (SD %.1f); BMI %.1f (SD %.1f); age %.1f (SD %.1f); women %.1f%%",
    mean(d$wc), sd(d$wc), mean(d$BMI), sd(d$BMI), mean(d$age), sd(d$age), 100 * mean(d$sex == 0))

keep <- pv[sapply(pv, function(p) mean(is.na(d[[p]])) <= 0.5)]
say("proteins with at most 50%% missing: %d", length(keep))
X <- as.matrix(d[, ..keep]); for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); X[, j] <- v }
X <- cbind(X, age = d$age, sex = d$sex); X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
say("predictors entering the model: %d (proteins, age and sex)", ncol(X))

## ---------------- proteomic waist, every prediction out of fold ----------------
## lambda is chosen by cross-validation inside each training fold, so the held-out tenth never informs the fit
y <- d$wc; fold <- sample(rep(1:10, length.out = nrow(d))); d[, pWC := NA_real_]; nz <- numeric(0)
for (k in 1:10) { tr <- fold != k
  cv <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = 1, nfolds = 10)
  d$pWC[!tr] <- as.numeric(predict(cv, s = "lambda.min", newx = X[!tr, , drop = FALSE]))
  nz <- c(nz, cv$nzero[which(cv$lambda == cv$lambda.min)]) }
r2 <- 1 - sum((d$pWC - y)^2) / sum((y - mean(y))^2)
say("out-of-fold R2 for the measured waist: %.3f (correlation %.3f); median %.0f terms retained; mean absolute error %.2f cm",
    r2, cor(d$pWC, y), median(nz), mean(abs(d$pWC - y)))

## ---------------- discordance, on the centimetre scale, exactly as in the manuscript ----------------
an <- lm(pWC ~ wc, d); say("anchoring regression: pWC = %.2f + %.3f x WC", coef(an)[1], coef(an)[2])
d[, dlt := pWC - (coef(an)[1] + coef(an)[2] * wc)]
d[, z := dlt / sd(dlt)]
say("proWCdelta: SD %.2f cm, range %.1f to %.1f; correlation with measured waist %.3f, with BMI %.3f",
    sd(d$dlt), min(d$dlt), max(d$dlt), cor(d$dlt, d$wc), cor(d$dlt, d$BMI))

## ---------------- does the discordance separate incident metabolic syndrome? ----------------
## Linear terms for waist and BMI are the primary specification: with 88 events, splines for both would leave
## fewer than ten events per parameter. The spline version is reported alongside as a sensitivity analysis.
fit <- function(rhs, lab, x = d) { m <- glm(as.formula(paste("ev ~", rhs)), x, family = binomial()); co <- summary(m)$coefficients
  if (!"z" %in% rownames(co)) return(NULL)
  data.table(model = lab, n = nrow(x), events = sum(x$ev), OR = exp(co["z", 1]),
             lo = exp(co["z", 1] - 1.96 * co["z", 2]), hi = exp(co["z", 1] + 1.96 * co["z", 2]), p = co["z", 4]) }
res <- rbindlist(list(
  fit("z", "proWCdelta alone"),
  fit("z + age + factor(sex)", "age and sex"),
  fit("z + age + factor(sex) + wc", "age, sex, measured waist"),
  fit("z + age + factor(sex) + wc + BMI", "age, sex, measured waist and BMI"),
  fit("z + age + factor(sex) + ns(wc,3) + ns(BMI,3)", "age, sex, splines of waist and BMI"),
  fit("z + age + factor(sex) + wc + BMI + factor(Indicator)", "age, sex, waist, BMI and the baseline component")))
res[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]; print(res[, .(model, n, events, txt, p = signif(p, 3))])
fwrite(res, file.path(TB, "G5_incident_mets.csv"))

## the same on the centimetre scale, so it can be compared with the manuscript directly
m5 <- glm(ev ~ I(dlt / 5) + age + factor(sex) + wc + BMI, d, family = binomial()); c5 <- summary(m5)$coefficients
say("per 5 cm of proWCdelta at the same measured waist and BMI: OR %.2f (%.2f-%.2f)",
    exp(c5[2, 1]), exp(c5[2, 1] - 1.96 * c5[2, 2]), exp(c5[2, 1] + 1.96 * c5[2, 2]))

## ---------------- the manuscript's headline contrast, among those below the waist threshold ----------------
nb <- d[waist == 0]; nb[, hi := as.integer(dlt > median(dlt))]
say("below the measured waist threshold: %d participants, %d events; higher proteomic waist %d of them",
    nrow(nb), sum(nb$ev), sum(nb$hi))
mh <- glm(ev ~ hi + age + factor(sex) + wc + BMI, nb, family = binomial()); ch <- summary(mh)$coefficients
say("normal measured waist, upper half of proWCdelta vs lower half: OR %.2f (%.2f-%.2f), p = %.3f",
    exp(ch["hi", 1]), exp(ch["hi", 1] - 1.96 * ch["hi", 2]), exp(ch["hi", 1] + 1.96 * ch["hi", 2]), ch["hi", 4])

## ---------------- does the binary flag used in the larger cohorts give the same thing? ----------------
## G4 built a score on the threshold flag in all 1,783 discovery participants; these 485 are a subset of them.
sc <- tryCatch(fread(file.path(TB, "G_score_discovery.csv")), error = function(e) NULL)
if (!is.null(sc)) { cm <- merge(d[, .(Patient_ID, dlt, z_cont = z, wc, BMI)], sc[, .(Patient_ID, z_flag = z)], by = "Patient_ID")
  say("participants with both scores: %d; correlation of the continuous-waist discordance with the flag-based one: Pearson %.3f, Spearman %.3f",
      nrow(cm), cor(cm$z_cont, cm$z_flag), cor(cm$z_cont, cm$z_flag, method = "spearman"))
  fwrite(cm, file.path(TB, "G5_score_comparison.csv")) } else say("G_score_discovery.csv not found; skipping the comparison")

fwrite(d[, .(Patient_ID, Indicator, ev, sex, age, wc, BMI, waist, pWC, dlt, z)], file.path(TB, "G5_scores.csv"))
fwrite(data.table(metric = c("n", "events", "oof_R2", "oof_r", "median_terms", "MAE_cm", "anchor_intercept", "anchor_slope", "SD_delta_cm"),
                  value = c(nrow(d), sum(d$ev), r2, cor(d$pWC, y), median(nz), mean(abs(d$pWC - y)), coef(an)[1], coef(an)[2], sd(d$dlt))),
       file.path(TB, "G5_model_summary.csv"))
say("DONE")
