## 65c_figure2_panels_source_data.R: aggregate source data for the panels of Figure 2 that were previously drawn
## directly from individual-level data: the joint distribution of WC and proWC by quintile of proWCdelta (panel F),
## the distributions of WC and proWC by quintile (panels G and H) and the mean z-scores of the six displayed proteins
## by quintile (panel I). Output: T98 (quintile summaries), T99 (protein z-scores by quintile).
suppressPackageStartupMessages(library(data.table))
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
l <- fread(file.path(F12, "lasso_WC.csv")); setnames(l, 1, "id")
d <- l[, .(id, WC = num(Actual_WC), pWC = num(pred_WC), proWC = num(BioX_Adjusted), dlt = num(BioX_Delta))]
d <- d[complete.cases(d)]; d[, quintile := cut(dlt, quantile(dlt, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
q <- function(v) c(mean = mean(v), sd = sd(v), p25 = quantile(v, .25), median = median(v), p75 = quantile(v, .75), min = min(v), max = max(v))
t98 <- d[, as.list(c(n = .N, setNames(q(WC), paste0("WC_", names(q(WC)))), setNames(q(proWC), paste0("proWC_", names(q(proWC)))), setNames(q(dlt), paste0("proWCdelta_", names(q(dlt)))))), by = quintile][order(quintile)]
print(t98[, .(quintile, n, WC_mean = round(WC_mean, 2), proWC_mean = round(proWC_mean, 2), proWCdelta_mean = round(proWCdelta_mean, 2))])
fwrite(t98, file.path(W, "tables", "T98_figure2_quintile_distributions.csv"))
mf <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(mf, 1, "id")
pp <- c("LEP", "SSC4D", "IGFBP1", "OXT", "GH1", "FGF21"); pp <- intersect(make.names(pp), names(mf))
m <- merge(d[, .(id, quintile)], mf[, c("id", ..pp)], by = "id")
z <- melt(m, id.vars = c("id", "quintile"), variable.name = "protein", value.name = "npx")
z[, npx := num(npx)]; z[, zscore := (npx - mean(npx, na.rm = TRUE)) / sd(npx, na.rm = TRUE), by = protein]
t99 <- z[, .(n = sum(!is.na(zscore)), mean_z = mean(zscore, na.rm = TRUE), sd_z = sd(zscore, na.rm = TRUE),
             p25_z = quantile(zscore, .25, na.rm = TRUE), median_z = median(zscore, na.rm = TRUE), p75_z = quantile(zscore, .75, na.rm = TRUE)), by = .(protein, quintile)][order(protein, quintile)]
print(t99[, .(protein, quintile, n, mean_z = round(mean_z, 3))])
fwrite(t99, file.path(W, "tables", "T99_figure2I_protein_zscores_by_quintile.csv"))
say("DONE")
