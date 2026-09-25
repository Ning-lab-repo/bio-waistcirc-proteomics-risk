## G9_direct_transport.R
## Everything so far repeats the procedure inside the Guangzhou cohorts; nothing has carried coefficients across. This
## script does that, which is the stricter test and the one a reader will ask about. A waist model is trained here on
## the 132 proteins that both studies measure, with age and sex, on standardised values so that the weights do not
## depend on the units of either assay; those weights are then applied unchanged to the standardised Guangzhou
## measurements of the same proteins, and the resulting prediction is compared with the waist actually measured there.
## The discordance of that transported prediction is then tested against incident metabolic syndrome as before.
## A transported score cannot be expected to do well: the two assays measure different things about the same protein,
## and standardising within each cohort removes the mean and scale but not the difference in what is being captured.
## The point is to establish whether anything survives the crossing, and to report it either way.
## Output: G9_transport.csv.
suppressPackageStartupMessages({ library(data.table); library(readxl); library(glmnet) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260923)

## ---------------- the proteins both studies measure ----------------
E <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "E", na = c("", "NA", "NaN")))
CLIN <- c("age","BMI","DBP","SBP","fbg","Glu","HDL","hdl","LDL","TG","tg","wc","waist","sex","sbpdbp","mets","metsscore",
          "TC","time","pat_ID","Sample_ID","DM_med","dyslipi_med","hyper_med","dm_cl","hyper_cl","tg_cl","tc_cl",
          "hdl_cl","ldl_cl","dys_cl","menopau_age")
gcol <- setdiff(names(E), CLIN); names(gcol) <- sub("^[^_]*_", "", gcol)
u <- fread(file.path(F12, "analysis_data_WC.csv")); uprot <- setdiff(names(u), c(names(u)[1], "WC", "Age", "Sex"))
u <- u[complete.cases(u[, .(WC, Age, Sex)])]
sh <- intersect(names(gcol), uprot)
## keep only those the Guangzhou participants actually have measured often enough to use
gm <- sapply(sh, function(g) mean(is.na(as.numeric(E[[gcol[[g]]]]))))
sh <- sh[gm <= 0.5]
say("proteins on both platforms and measured in at least half the Guangzhou participants: %d", length(sh))

## ---------------- train here, on standardised values ----------------
Xu <- as.matrix(u[, ..sh]); for (j in seq_len(ncol(Xu))) { v <- Xu[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); Xu[, j] <- v }
Xu <- scale(Xu); Xu <- cbind(Xu, Age = as.numeric(scale(u$Age)), Sex = as.integer(u$Sex))
cv <- cv.glmnet(Xu, as.numeric(u$WC), family = "gaussian", alpha = 1, nfolds = 10)
b <- as.matrix(coef(cv, s = "lambda.min")); b <- b[b[, 1] != 0, , drop = FALSE]
say("UK Biobank model on those proteins: %d terms retained; in-cohort R2 %.3f",
    nrow(b) - 1, 1 - min(cv$cvm) / var(as.numeric(u$WC)))

## ---------------- apply the weights unchanged in Guangzhou ----------------
d <- merge(as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "A", na = c("", "NA", "NaN")))[
             , .(Patient_ID = .SD[[1]], ev = as.integer(.SD[[5]] == 1))],
           E[, c(list(Patient_ID = pat_ID, wc = as.numeric(wc), BMI = as.numeric(BMI), age = as.numeric(age), sex = as.integer(sex)),
                 lapply(.SD, as.numeric)), .SDcols = unname(gcol[sh])], by = "Patient_ID")
setnames(d, unname(gcol[sh]), sh)
Xg <- as.matrix(d[, ..sh]); for (j in seq_len(ncol(Xg))) { v <- Xg[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); Xg[, j] <- v }
Xg <- scale(Xg); Xg <- cbind(Xg, Age = as.numeric(scale(d$age)), Sex = d$sex)
keep <- intersect(rownames(b)[-1], colnames(Xg))
d[, pWC := b["(Intercept)", 1] + as.numeric(Xg[, keep, drop = FALSE] %*% b[keep, 1])]
say("weights carried across: %d of %d", length(keep), nrow(b) - 1)
say("transported prediction against the waist measured in Guangzhou: Pearson %.3f, Spearman %.3f, R2 of a linear fit %.3f",
    cor(d$pWC, d$wc), cor(d$pWC, d$wc, method = "spearman"), summary(lm(wc ~ pWC, d))$r.squared)
say("  mean predicted %.1f cm against %.1f cm measured; mean absolute error %.2f cm", mean(d$pWC), mean(d$wc), mean(abs(d$pWC - d$wc)))
say("  for reference, the score trained inside Guangzhou reached Pearson 0.398")

## ---------------- does the transported discordance carry risk? ----------------
an <- lm(pWC ~ wc, d); d[, dlt := pWC - (coef(an)[1] + coef(an)[2] * wc)]; d[, z := dlt / sd(dlt)]
out <- rbindlist(lapply(list(c("z", "age, sex"), c("z + wc + BMI", "age, sex, measured waist and BMI")), function(sp) {
  m <- glm(as.formula(paste("ev ~", sp[1], "+ age + factor(sex)")), d, family = binomial()); co <- summary(m)$coefficients
  data.table(model = sp[2], n = nrow(d), events = sum(d$ev), OR = exp(co["z", 1]),
             lo = exp(co["z", 1] - 1.96 * co["z", 2]), hi = exp(co["z", 1] + 1.96 * co["z", 2]), p = co["z", 4]) }))
out[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]; print(out[, .(model, n, events, txt, p = signif(p, 3))])
sc <- fread(file.path(TB, "G5_scores.csv"))
say("correlation of the transported discordance with the one trained inside Guangzhou: %.3f",
    cor(d$z, merge(d[, .(Patient_ID)], sc[, .(Patient_ID, zin = z)], by = "Patient_ID")$zin))
fwrite(out, file.path(TB, "G9_transport.csv")); say("DONE")
