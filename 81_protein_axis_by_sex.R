## 81_protein_axis_by_sex.R
## The alignment of the proWCdelta protein axis with visceral rather than subcutaneous fat (80_protein_axis_checks.R),
## repeated within women and within men. Output: T117_protein_axis_by_sex.csv
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
prot <- fread(file.path(TB, "T115_protein_vat_associations.csv"))$protein
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", prot))
setnames(phen, 1:5, c("id", "Sex", "Age", "WC", "BMI"))
la <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE, select = c("Participant.ID",
  "Visceral.adipose.tissue.volume..VAT....Instance.2", "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
all <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI")) set(all, j = v, value = num(all[[v]])); all[, Sex := as.integer(Sex)]
all <- all[complete.cases(all[, .(Age, Sex, WC, BMI, dlt)])]; sdz <- sd(all$dlt); all[, z := dlt / sdz]
for (p in prot) set(all, j = p, value = num(all[[p]]))
d <- merge(all, mri, by = "id"); for (v in c("vat", "asat")) set(d, j = v, value = num(d[[v]])); d <- d[!is.na(vat) & !is.na(asat) & asat > 0]
B <- "Age + ns(WC,3) + ns(BMI,3)"
rows <- list()
for (s in 0:1) { lab <- c("women", "men")[s + 1]; x0 <- d[Sex == s]; a0 <- all[Sex == s]
  r <- rbindlist(mclapply(prot, function(p) { x <- x0[!is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
    y <- a0[!is.na(get(p))]; y[, pz := as.numeric(scale(get(p)))]
    data.table(protein = p,
      bv = coef(lm(as.formula(paste("vat ~ pz +", B)), x))[["pz"]],
      ba = coef(lm(as.formula(paste("asat ~ pz +", B)), x))[["pz"]],
      bd = coef(lm(as.formula(paste("pz ~ z +", B)), y))[["z"]]) }, mc.cores = 20))
  bt <- replicate(2000, { i <- sample(nrow(r), replace = TRUE); cor(r$bv[i], r$bd[i]) - cor(r$ba[i], r$bd[i]) })
  say("%s (MRI n = %d): r(visceral, proWCdelta) %.2f; r(subcutaneous, proWCdelta) %.2f; difference %.2f (%.2f to %.2f)",
      lab, nrow(x0), cor(r$bv, r$bd), cor(r$ba, r$bd), cor(r$bv, r$bd) - cor(r$ba, r$bd), quantile(bt, .025), quantile(bt, .975))
  rows[[lab]] <- data.table(sex = lab, n_mri = nrow(x0), r_visceral = cor(r$bv, r$bd), r_subcutaneous = cor(r$ba, r$bd),
                            diff = cor(r$bv, r$bd) - cor(r$ba, r$bd), lo = quantile(bt, .025)[[1]], hi = quantile(bt, .975)[[1]]) }
fwrite(rbindlist(rows), file.path(TB, "T117_protein_axis_by_sex.csv")); say("DONE")
