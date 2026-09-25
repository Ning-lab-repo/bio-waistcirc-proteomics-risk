## G2_gnhs_analysis_set.R
## Builds the baseline analysis set of the Guangzhou Nutrition and Health Study serum proteomics data (Cai et al., Cell
## Reports Medicine 2023;4:101172, supplementary tables S1-S3) and compares the protein associations with central
## obesity there with the protein associations with waist circumference in UK Biobank. The two cohorts of that study
## (discovery, S2; validation, S3) were assayed on overlapping but different protein panels and are kept separate.
## Central obesity in that study is recorded as the metabolic-syndrome flag (1 if the waist exceeds the sex-specific
## threshold), not as a measurement in centimetres, so the waist cannot be modelled on a continuous scale here.
## Output: G_set_*.csv (analysis sets), G_protein_assoc.csv, G_overlap.csv.
suppressPackageStartupMessages({ library(data.table) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

cl <- fread(file.path(TB, "S1_clinical.csv"))
setnames(cl, c("SBP/DBP_BS", "SBP/DBP_F2", "SBP/DBP_F3", "BMI__F2", "BMI__F3", "Age__F2", "Age__F3"),
             c("BP_BS", "BP_F2", "BP_F3", "BMI_F2", "BMI_F3", "Age_F2", "Age_F3"))
say("999 codes in Waist_BS: %d", sum(cl$Waist_BS == 999, na.rm = TRUE)); cl[Waist_BS == 999, Waist_BS := NA]
for (v in grep("^(Waist|Glucose|TG|HDL|BP|MetS)_", names(cl), value = TRUE)) cl[[v]][cl[[v]] == 999] <- NA

## ---------------- one baseline proteomic sample per participant, per cohort ----------------
build <- function(n) { a <- fread(file.path(TB, sprintf("S%d_meta.csv", n))); b <- fread(file.path(TB, sprintf("S%d_matrix.csv", n)))
  setnames(a, names(a)[2], "time"); setnames(b, 1, "Sample_ID")
  a <- a[time == "baseline" & Sample_ID %in% b$Sample_ID]
  a <- a[!duplicated(Patient_ID)]                       ## one sample per participant; replicates dropped
  m <- merge(a[, .(Patient_ID, Sample_ID, Batch_ID)], b, by = "Sample_ID")
  say("S%d: %d baseline samples with proteins, %d participants", n, nrow(m), uniqueN(m$Patient_ID)); m }
s2 <- build(2); s3 <- build(3)
say("participants in both cohorts: %d", length(intersect(s2$Patient_ID, s3$Patient_ID)))

prot_of <- function(m) setdiff(names(m), c("Sample_ID", "Patient_ID", "Batch_ID"))
for (nm in c("s2", "s3")) { m <- get(nm); pv <- prot_of(m)
  miss <- sapply(pv, function(p) mean(is.na(m[[p]])))
  say("%s: %d proteins, median missing %.1f%%, proteins with more than 50%% missing: %d", nm, length(pv), 100 * median(miss), sum(miss > 0.5)) }

mk <- function(m, nm) { d <- merge(cl, m, by = "Patient_ID")
  d <- d[!is.na(Waist_BS) & !is.na(BMI_BS) & !is.na(Age_BS) & !is.na(Sex)]
  say("%s analysis set: %d participants, central obesity %.1f%%, women %.1f%%, age %.1f (SD %.1f), BMI %.1f (SD %.1f)",
      nm, nrow(d), 100 * mean(d$Waist_BS), 100 * mean(d$Sex == 0), mean(d$Age_BS), sd(d$Age_BS), mean(d$BMI_BS), sd(d$BMI_BS))
  fwrite(d, file.path(TB, sprintf("G_set_%s.csv", nm))); d }
d2 <- mk(s2, "discovery"); d3 <- mk(s3, "validation")

## ---------------- protein associations with central obesity, adjusted for age and sex ----------------
assoc <- function(d, nm) { pv <- prot_of(d)[prot_of(d) %in% names(d)]
  pv <- pv[sapply(pv, function(p) mean(is.na(d[[p]])) <= 0.5)]
  res <- rbindlist(lapply(pv, function(p) { x <- d[[p]]; ok <- !is.na(x); if (sum(ok) < 200) return(NULL)
    dd <- d[ok]; dd[, xx := scale(x[ok])]; m <- glm(Waist_BS ~ xx + Age_BS + factor(Sex), dd, family = binomial()); co <- summary(m)$coefficients
    data.table(cohort = nm, protein = p, n = sum(ok), beta = co[2, 1], se = co[2, 2], p = co[2, 4]) }))
  res[, fdr := p.adjust(p, "BH")]; say("%s: %d proteins tested, %d with FDR < 0.05", nm, nrow(res), sum(res$fdr < 0.05)); res }
a2 <- assoc(d2, "discovery"); a3 <- assoc(d3, "validation")
fwrite(rbind(a2, a3), file.path(TB, "G_protein_assoc.csv"))

## how well do the two cohorts of that study agree with each other
both <- merge(a2[, .(protein, b2 = beta, f2 = fdr)], a3[, .(protein, b3 = beta, f3 = fdr)], by = "protein")
say("proteins in both cohorts: %d; correlation of the estimates r = %.3f; same direction %.1f%%",
    nrow(both), cor(both$b2, both$b3), 100 * mean(sign(both$b2) == sign(both$b3)))
say("of the %d with FDR < 0.05 in discovery, %d replicate at FDR < 0.05 in validation with the same sign (%.1f%%)",
    sum(both$f2 < 0.05), sum(both$f2 < 0.05 & both$f3 < 0.05 & sign(both$b2) == sign(both$b3)),
    100 * sum(both$f2 < 0.05 & both$f3 < 0.05 & sign(both$b2) == sign(both$b3)) / sum(both$f2 < 0.05))

## ---------------- overlap with the UK Biobank waist model ----------------
gene <- function(p) sub("^[^_]*_", "", p)                       ## the columns are UniProt_GENE
ukb <- fread(file.path(F12, "analysis_data_WC.csv"), nrows = 0); uprot <- setdiff(names(ukb), c(names(ukb)[1], "WC", "Age", "Sex"))
gset <- unique(c(gene(a2$protein), gene(a3$protein)))
say("GNHS proteins mapped to a gene symbol: %d; also measured in UK Biobank: %d", length(gset), sum(gset %in% uprot))
fwrite(data.table(gene = gset, in_ukb = gset %in% uprot), file.path(TB, "G_overlap.csv"))
say("DONE")
