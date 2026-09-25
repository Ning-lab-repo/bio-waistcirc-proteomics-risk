## 83_protein_axis_comparator.R
## The protein profile of proWCdelta lines up with that of visceral fat (r = 0.81) and not with that of subcutaneous
## fat (r = 0.12; 80_protein_axis_checks.R). Is that specific to proWCdelta, or would any adiposity signal at fixed BMI
## look the same? The comparator is a larger measured waist at the same BMI, which is also what error in the tape would
## produce. For each protein: its association with measured WC at fixed BMI (age, sex, spline of BMI; whole cohort),
## compared across proteins with its associations with visceral and subcutaneous fat (T115).
## Also the split-half reliability of the per-protein fat associations (random halves of the imaging participants).
## Output: T119_protein_axis_comparator.csv
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
set.seed(20260925)
pv <- fread(file.path(TB, "T115_protein_vat_associations.csv")); prot <- pv$protein
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", prot))
setnames(phen, 1:5, c("id", "Sex", "Age", "WC", "BMI"))
la <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
all <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI")) set(all, j = v, value = num(all[[v]])); all[, Sex := as.integer(Sex)]
all <- all[complete.cases(all[, .(Age, Sex, WC, BMI, dlt)])]
for (p in prot) set(all, j = p, value = num(all[[p]]))
all[, wz := as.numeric(scale(WC))]

## 1. association of each protein with measured WC at fixed BMI (the comparator axis)
pw <- rbindlist(mclapply(prot, function(p) { x <- all[!is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
  co <- summary(lm(pz ~ wz + ns(BMI, 3) + Age + factor(Sex), x))$coefficients["wz", ]
  data.table(protein = p, beta_wc = co[1]) }, mc.cores = 20))
pv <- merge(pv, pw, by = "protein")
r <- function(a, b) cor(pv[[a]], pv[[b]])
res <- data.table(
  comparison = c("proWCdelta vs visceral fat", "proWCdelta vs subcutaneous fat",
                 "larger measured WC at fixed BMI vs visceral fat", "larger measured WC at fixed BMI vs subcutaneous fat",
                 "proWCdelta vs larger measured WC at fixed BMI"),
  r = c(r("beta_delta", "beta_vat"), r("beta_delta", "beta_asat"), r("beta_wc", "beta_vat"), r("beta_wc", "beta_asat"),
        r("beta_delta", "beta_wc")))
say("correlations across %d proteins:", nrow(pv)); print(res)

## 2. split-half reliability of the per-protein associations with each fat depot (imaging participants)
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE, select = c("Participant.ID",
  "Visceral.adipose.tissue.volume..VAT....Instance.2", "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
d <- merge(all, mri, by = "id"); for (v in c("vat", "asat")) set(d, j = v, value = num(d[[v]])); d <- d[!is.na(vat) & !is.na(asat) & asat > 0]
d[, half := sample(rep(1:2, length.out = .N))]
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
sh <- rbindlist(mclapply(prot, function(p) { out <- list()
  for (h in 1:2) { x <- d[half == h & !is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
    out[[h]] <- c(coef(lm(as.formula(paste("vat ~ pz +", B0)), x))[["pz"]], coef(lm(as.formula(paste("asat ~ pz +", B0)), x))[["pz"]]) }
  data.table(protein = p, vat1 = out[[1]][1], vat2 = out[[2]][1], asat1 = out[[1]][2], asat2 = out[[2]][2]) }, mc.cores = 20))
rel <- data.table(comparison = c("split-half reliability, visceral profile", "split-half reliability, subcutaneous profile"),
                  r = c(cor(sh$vat1, sh$vat2), cor(sh$asat1, sh$asat2)))
say("split-half reliability of the per-protein associations:"); print(rel)
## correlations with proWCdelta corrected for the reliability of each depot profile (Spearman-Brown for the full sample)
sb <- function(r) 2 * r / (1 + r)
res2 <- data.table(comparison = c("proWCdelta vs visceral fat, corrected for reliability", "proWCdelta vs subcutaneous fat, corrected for reliability"),
                   r = c(r("beta_delta", "beta_vat") / sqrt(sb(rel$r[1])), r("beta_delta", "beta_asat") / sqrt(sb(rel$r[2]))))
print(res2)
fwrite(rbind(res, rel, res2), file.path(TB, "T119_protein_axis_comparator.csv"))
fwrite(pv, file.path(TB, "T115_protein_vat_associations.csv"))
say("DONE")
