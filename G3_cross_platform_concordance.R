## G3_cross_platform_concordance.R
## Do the protein associations with central adiposity that we report in UK Biobank hold in an independently recruited
## cohort, on a different platform and in a different population? Compares the per-protein associations with measured
## waist circumference in UK Biobank (Olink, plasma, adjusted for age and sex) with the associations with the central
## obesity flag in the Guangzhou Nutrition and Health Study (mass spectrometry, serum, same adjustment), over the
## proteins measured in both. The two scales differ (cm per SD against log odds per SD), so only the direction and the
## rank of the associations can be compared. Output: G_cross_platform.csv.
suppressPackageStartupMessages({ library(data.table) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
D <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

u <- fread(file.path(D, "source_table_Figure2E_per_protein_WC_vs_proWC_coefficients.csv")); say("UK Biobank table: %s", paste(names(u), collapse = ", "))
g <- fread(file.path(TB, "G_protein_assoc.csv")); g[, gene := sub("^[^_]*_", "", protein)]
## drop the fits that did not converge (complete separation on a near-constant protein)
g <- g[is.finite(beta) & se < 10 & abs(beta) < 10]
g <- g[!duplicated(paste(cohort, gene))]        ## a few columns share a gene symbol; keep one per gene and cohort
say("GNHS estimates kept after dropping non-converged fits: discovery %d, validation %d", nrow(g[cohort == "discovery"]), nrow(g[cohort == "validation"]))

both <- merge(g[cohort == "discovery", .(gene, b2 = beta, f2 = fdr)], g[cohort == "validation", .(gene, b3 = beta, f3 = fdr)], by = "gene")
say("GNHS internal: %d proteins in both cohorts, Spearman rho = %.3f, same direction %.1f%%, and %.1f%% of those with FDR < 0.05 in the first replicate in the second",
    nrow(both), cor(both$b2, both$b3, method = "spearman"), 100 * mean(sign(both$b2) == sign(both$b3)),
    100 * mean(both[f2 < 0.05, f3 < 0.05 & sign(b2) == sign(b3)]))

uc <- "measured_beta"; pc <- "protein"     ## coefficient of each protein on measured WC, adjusted for age and sex
say("using UK Biobank columns: protein = %s, waist coefficient = %s", pc, uc)
u2 <- u[, .(gene = get(pc), bu = as.numeric(get(uc)))]

for (ch in c("discovery", "validation")) { m <- merge(u2, g[cohort == ch], by = "gene")
  say("%s: %d proteins measured in both studies; Spearman rho with the UK Biobank waist coefficient = %.3f; same direction %.1f%%",
      ch, nrow(m), cor(m$bu, m$beta, method = "spearman"), 100 * mean(sign(m$bu) == sign(m$beta)))
  sig <- m[fdr < 0.05]; say("  of the %d associated with central obesity in that cohort (FDR < 0.05), %.1f%% have the same direction in UK Biobank",
      nrow(sig), 100 * mean(sign(sig$bu) == sign(sig$beta)))
  fwrite(m, file.path(TB, sprintf("G_cross_platform_%s.csv", ch))) }

## the proteins that contribute most to the UK Biobank score, in the Chinese cohorts
con <- fread(file.path(D, "source_table_Results_protein_contributions.csv")); say("contribution table: %s", paste(names(con), collapse = ", "))
cn <- names(con)[grepl("protein|gene", tolower(names(con)))][1]; vn <- names(con)[grepl("contrib|abs|sd", tolower(names(con)))][1]
top <- head(con[order(-abs(as.numeric(get(vn))))][[cn]], 20)
say("of the 20 proteins contributing most to pWC in UK Biobank, %d are measured in the Chinese cohorts: %s",
    sum(top %in% g$gene), paste(intersect(top, g$gene), collapse = ", "))
pr <- merge(data.table(gene = top), g[cohort == "discovery", .(gene, beta, fdr)], by = "gene")
if (nrow(pr)) print(pr[order(-abs(beta))][, .(gene, OR_per_SD = sprintf("%.2f", exp(beta)), fdr = sprintf("%.3g", fdr))])
say("DONE")
