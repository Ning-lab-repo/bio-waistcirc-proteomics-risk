## 80_protein_axis_checks.R
## Does the protein axis of proWCdelta line up with visceral fat specifically, and does the visceral finding hold in each
## sex? proWCdelta was trained on the tape measure only; the MRI volumes were never used to build it.
##  1. per-protein associations with abdominal subcutaneous fat at fixed WC and BMI, and with visceral fat at fixed WC,
##     BMI and subcutaneous fat, compared (across 2,920 proteins) with each protein's association with proWCdelta
##  2. partial R2 of proWCdelta for visceral and for subcutaneous fat, within each sex
##  3. visceral fat per SD of proWCdelta with the WC and BMI of both visits held fixed, within each sex
## Reads T115 (79_mri_checks.R). Output: T116_protein_axis_checks.csv, T115 extended with the new columns.
suppressPackageStartupMessages({ library(data.table); library(splines); library(parallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
SD  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

pv <- fread(file.path(TB, "T115_protein_vat_associations.csv")); prot <- pv$protein
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", prot))
setnames(phen, 1:5, c("id", "Sex", "Age", "WC", "BMI"))
la  <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
bs  <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"),
             select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2"))
setnames(bs, c("id", "WC2", "BMI2"))
all <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id")
for (v in c("Age", "WC", "BMI")) set(all, j = v, value = num(all[[v]])); all[, Sex := as.integer(Sex)]
all <- all[complete.cases(all[, .(Age, Sex, WC, BMI, dlt)])]
sdz <- sd(all$dlt)
d <- merge(merge(all, mri, by = "id"), bs, by = "id", all.x = TRUE)
for (v in c("vat", "asat", "WC2", "BMI2")) set(d, j = v, value = num(d[[v]]))
d <- d[!is.na(vat) & !is.na(asat) & asat > 0]; d[, z := dlt / sdz]
for (p in prot) set(d, j = p, value = num(d[[p]]))
say("imaging participants %d", nrow(d))
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"

## ---------------- 1. protein axes ----------------
pa <- rbindlist(mclapply(prot, function(p) { x <- d[!is.na(get(p))]; x[, pz := as.numeric(scale(get(p)))]
  a <- summary(lm(as.formula(paste("asat ~ pz +", B0)), x))$coefficients["pz", ]
  v <- summary(lm(as.formula(paste("vat ~ pz + ns(asat,3) +", B0)), x))$coefficients["pz", ]
  data.table(protein = p, beta_asat = a[1], p_asat = a[4], beta_vat_given_asat = v[1], p_vat_given_asat = v[4]) }, mc.cores = 20))
pv <- merge(pv[, setdiff(names(pv), setdiff(names(pa), "protein")), with = FALSE], pa, by = "protein")   ## rerun-safe
rows <- list(); add <- function(check, spec, n, est, lo = NA_real_, hi = NA_real_) rows[[length(rows) + 1]] <<- data.table(check, spec, n, est, lo, hi)
r1 <- cor(pv$beta_vat, pv$beta_delta); r2 <- cor(pv$beta_asat, pv$beta_delta); r3 <- cor(pv$beta_vat_given_asat, pv$beta_delta)
r4 <- cor(pv$beta_vat, pv$beta_asat)
## the comparison of r1 and r2 by bootstrap over proteins
bt <- replicate(2000, { i <- sample(nrow(pv), replace = TRUE); cor(pv$beta_vat[i], pv$beta_delta[i]) - cor(pv$beta_asat[i], pv$beta_delta[i]) })
say("[1] across %d proteins: r(visceral, proWCdelta) %.2f; r(subcutaneous, proWCdelta) %.2f (difference %.2f, %.2f to %.2f);",
    nrow(pv), r1, r2, r1 - r2, quantile(bt, .025), quantile(bt, .975))
say("    r(visceral given subcutaneous, proWCdelta) %.2f; r(visceral, subcutaneous) %.2f", r3, r4)
say("    proteins associated (FDR < 0.05): visceral %d; subcutaneous %d; visceral given subcutaneous %d",
    sum(p.adjust(pv$p_vat, "BH") < 0.05), sum(p.adjust(pv$p_asat, "BH") < 0.05), sum(p.adjust(pv$p_vat_given_asat, "BH") < 0.05))
add("1. protein axes", "r across proteins: visceral fat vs proWCdelta", nrow(pv), r1)
add("1. protein axes", "r across proteins: subcutaneous fat vs proWCdelta", nrow(pv), r2)
add("1. protein axes", "difference, lower 95% bootstrap limit", nrow(pv), quantile(bt, .025)[[1]])
add("1. protein axes", "difference, upper 95% bootstrap limit", nrow(pv), quantile(bt, .975)[[1]])
add("1. protein axes", "r across proteins: visceral fat at fixed subcutaneous fat vs proWCdelta", nrow(pv), r3)
add("1. protein axes", "r across proteins: visceral vs subcutaneous fat", nrow(pv), r4)
add("1. protein axes", "proteins associated with visceral fat (FDR < 0.05)", nrow(pv), sum(p.adjust(pv$p_vat, "BH") < 0.05))
add("1. protein axes", "proteins associated with subcutaneous fat (FDR < 0.05)", nrow(pv), sum(p.adjust(pv$p_asat, "BH") < 0.05))
add("1. protein axes", "proteins associated with visceral fat at fixed subcutaneous fat (FDR < 0.05)", nrow(pv), sum(p.adjust(pv$p_vat_given_asat, "BH") < 0.05))
setorder(pv, p_vat_given_asat)
say("    top proteins for visceral fat at fixed subcutaneous fat, WC and BMI:")
print(pv[1:15, .(protein, beta_vat_given_asat = round(beta_vat_given_asat, 3), p = signif(p_vat_given_asat, 2), beta_delta = round(beta_delta, 3))])

## ---------------- 2. partial R2 within each sex ----------------
pr2 <- function(y, x, f0) { m0 <- lm(as.formula(paste(y, "~", f0)), x); m1 <- lm(as.formula(paste(y, "~ z +", f0)), x)
  (sum(resid(m0)^2) - sum(resid(m1)^2)) / sum(resid(m0)^2) }
for (s in list(list("women", d[Sex == 0]), list("men", d[Sex == 1]))) {
  f0 <- "Age + ns(WC,3) + ns(BMI,3)"
  a <- pr2("vat", s[[2]], f0); b <- pr2("asat", s[[2]], f0)
  say("[2] %s: partial R2 of proWCdelta, visceral %.3f, subcutaneous %.3f", s[[1]], a, b)
  add("2. partial R2", paste0(s[[1]], ": visceral fat"), nrow(s[[2]]), a); add("2. partial R2", paste0(s[[1]], ": subcutaneous fat"), nrow(s[[2]]), b) }

## ---------------- 3. both visits fixed, within each sex ----------------
for (s in list(list("women", d[Sex == 0 & !is.na(WC2) & !is.na(BMI2)]), list("men", d[Sex == 1 & !is.na(WC2) & !is.na(BMI2)]))) {
  x <- s[[2]]
  for (spec in list(c("baseline WC and BMI fixed", "vat ~ z + Age + ns(WC,3) + ns(BMI,3)"),
                    c("WC and BMI of both visits fixed", "vat ~ z + Age + ns(WC,3) + ns(BMI,3) + ns(WC2,3) + ns(BMI2,3)"))) {
    co <- summary(lm(as.formula(spec[2]), x))$coefficients["z", ]
    say("[3] %s, %s: %.3f (%.3f to %.3f), n = %d", s[[1]], spec[1], co[1], co[1] - 1.96 * co[2], co[1] + 1.96 * co[2], nrow(x))
    add("3. both visits", paste0(s[[1]], ": ", spec[1]), nrow(x), co[1], co[1] - 1.96 * co[2], co[1] + 1.96 * co[2]) } }

fwrite(rbindlist(rows), file.path(TB, "T116_protein_axis_checks.csv"))
fwrite(pv, file.path(TB, "T115_protein_vat_associations.csv"))
say("DONE")
