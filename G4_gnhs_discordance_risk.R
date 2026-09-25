## G4_gnhs_discordance_risk.R
## External test of the method in the Guangzhou Nutrition and Health Study: build the proteomic estimate of central
## obesity in the same way as the proteomic waist circumference (penalised regression on the proteins with age and sex,
## all predictions out of fold), take the part of it that is not explained by the participant's measured body size, and
## ask whether that residual separates the metabolic-syndrome components that develop over the next two visits, among
## participants free of each component at baseline and at the same measured body size.
## Two differences from the UK Biobank analysis are unavoidable and limit what this can show: the waist is recorded only
## as the sex-specific threshold flag, so body size is held fixed by that flag and by BMI rather than by splines of a
## measured waist; and the platform measures a largely different, higher-abundance part of the proteome (128 of its
## proteins are on the Olink panel used here). Output: G_risk_*.csv.
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(splines) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

comp <- c(Glucose = "Glucose", TG = "TG", HDL = "HDL", BP = "BP", MetS = "MetS", Waist = "Waist")
run <- function(nm) {
  d <- fread(file.path(TB, sprintf("G_set_%s.csv", nm)))
  pv <- setdiff(names(d), c("Patient_ID", "Sample_ID", "Batch_ID", grep("^(Waist|Glucose|TG|HDL|BP|MetS|BMI|Age)", names(d), value = TRUE), "Sex"))
  keep <- pv[sapply(pv, function(p) mean(is.na(d[[p]])) <= 0.5)]
  say("[%s] %d participants, %d proteins with at most 50%% missing", nm, nrow(d), length(keep))

  ## median imputation, as in the main model; the outcome plays no part in it
  X <- as.matrix(d[, ..keep]); for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); X[, j] <- v }
  X <- cbind(X, Age = d$Age_BS, Sex = d$Sex); X <- X[, apply(X, 2, function(v) sd(v) > 0), drop = FALSE]
  y <- d$Waist_BS

  ## out-of-fold proteomic estimate of central obesity, lambda chosen inside each training fold
  fold <- sample(rep(1:10, length.out = nrow(d))); d[, pCO := NA_real_]; nz <- numeric(0)
  for (k in 1:10) { tr <- fold != k; cv <- cv.glmnet(X[tr, ], y[tr], family = "binomial", alpha = 1, nfolds = 10)
    d$pCO[!tr] <- as.numeric(predict(cv, s = "lambda.min", newx = X[!tr, , drop = FALSE], type = "link"))
    nz <- c(nz, cv$nzero[which(cv$lambda == cv$lambda.min)]) }
  auc <- function(p, o) { r <- rank(p); (sum(r[o == 1]) - sum(o) * (sum(o) + 1) / 2) / (sum(o) * sum(1 - o)) }
  say("[%s] out-of-fold AUC for central obesity %.3f; median %d terms retained", nm, auc(d$pCO, y), median(nz))

  ## the part of the proteomic estimate that measured body size does not explain
  d[, disc := resid(lm(pCO ~ factor(Waist_BS) + ns(BMI_BS, 3) + Age_BS + factor(Sex), d))]
  d[, z := disc / sd(disc)]
  say("[%s] discordance SD %.3f; correlation with BMI %.3f", nm, sd(d$disc), cor(d$z, d$BMI_BS))

  out <- list()
  for (cp in names(comp)) { b <- d[[paste0(comp[[cp]], "_BS")]]; f2 <- d[[paste0(comp[[cp]], "_F2")]]; f3 <- d[[paste0(comp[[cp]], "_F3")]]
    seen <- (!is.na(f2)) | (!is.na(f3)); inc <- as.integer((!is.na(f2) & f2 == 1) | (!is.na(f3) & f3 == 1))
    x <- d[!is.na(b) & b == 0 & seen]; x[, ev := inc[!is.na(b) & b == 0 & seen]]
    if (nrow(x) < 100 || sum(x$ev) < 20) { say("[%s] %s: too few (%d free, %d events)", nm, cp, nrow(x), sum(x$ev)); next }
    ## the waist flag drops out of the body-size adjustment when the endpoint is central obesity itself, because
    ## everyone in that analysis is below the threshold at baseline
    wt <- if (uniqueN(x$Waist_BS) > 1) "factor(Waist_BS) + " else ""
    for (sp in c("z", paste0("z + ", wt, "ns(BMI_BS,3)"))) {
      m <- glm(as.formula(paste("ev ~", sp, "+ Age_BS + factor(Sex)")), x, family = binomial()); co <- summary(m)$coefficients
      out[[length(out) + 1]] <- data.table(cohort = nm, component = cp, adjustment = if (sp == "z") "age and sex" else "age, sex, body size",
        n = nrow(x), events = sum(x$ev), OR = exp(co["z", 1]), lo = exp(co["z", 1] - 1.96 * co["z", 2]), hi = exp(co["z", 1] + 1.96 * co["z", 2]), p = co["z", 4]) } }
  r <- rbindlist(out); fwrite(d[, .(Patient_ID, Sex, Age_BS, BMI_BS, Waist_BS, pCO, disc, z)], file.path(TB, sprintf("G_score_%s.csv", nm))); r }

res <- rbindlist(lapply(c("discovery", "validation"), run))
res[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]
print(dcast(res, component + adjustment ~ cohort, value.var = c("n", "events", "txt")), width = 250)
fwrite(res, file.path(TB, "G_risk_by_cohort.csv")); say("DONE")
