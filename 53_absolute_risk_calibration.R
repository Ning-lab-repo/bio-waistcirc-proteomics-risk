## 53_absolute_risk_calibration.R
## (1) Absolute 10-year risks by continuous proWCdelta at fixed body size: from the primary logistic models (age, sex,
##     Townsend index, smoking, natural splines of WC and BMI), marginally standardised risks with proWCdelta set to its
##     10th, 50th and 90th percentiles and to -1, 0, +1 and +2 SD for every participant, and the average risk difference
##     per SD; 95% CIs from 200 bootstrap samples.
## (2) Calibration of the out-of-fold proWC model by tenth of predicted WC, and accuracy (R2, RMSE, calibration slope)
##     by sex, age, BMI category, ethnicity, region and protein missingness.
## (3) Share of the log odds ratio per SD of proWCdelta removed by adjustment for body fat and clinical markers
##     (from the Figure 5F estimates).
## Output: T69, T70, T71.
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"; F8 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
S <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis")); setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
bf <- file.path(F8, "wc_pro_Batch.csv"); hdr <- names(fread(bf, nrows = 0)); pe <- match("BA", hdr) - 1; pcols <- hdr[2:pe]
raw <- fread(bf, select = c(hdr[1], pcols, "ethnicity", "UK Biobank assessment centre | Instance 0")); setnames(raw, c("id", pcols, "white", "centre"))
raw[, pna := rowMeans(is.na(as.matrix(raw[, ..pcols])))]; raw <- raw[, .(id, pna, white, centre)]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen, lasso[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id"), cov, raw))
for (v in c("Age", "WC", "BMI")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
all0 <- copy(d)
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]
x <- as.character(d$dx10); x[is.na(x)] <- ""
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Chronic kidney disease", "Liver disease"), code = c("E11", "I50", "N18", "K76"))
for (i in seq_len(nrow(dis))) d[[dis$code[i]]] <- as.integer(grepl(paste0("(^|[|])", dis$code[i]), x))
say("modelling set n = %d; SD of proWCdelta %.2f cm", nrow(d), sdd)
q <- quantile(d$z, c(0.1, 0.5, 0.9)); say("proWCdelta percentiles (cm): 10th %.1f, 50th %.1f, 90th %.1f", q[1] * sdd, q[2] * sdd, q[3] * sdd)
pts <- c(p10 = q[[1]], p50 = q[[2]], p90 = q[[3]], `-1 SD` = -1, `0 SD` = 0, `+1 SD` = 1, `+2 SD` = 2)
F <- "~ z + Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
stdrisk <- function(dd, cd) { m <- glm(as.formula(paste(cd, F)), dd, family = binomial())
  r <- sapply(pts, function(v) { nd <- copy(dd); nd[, z := v]; mean(predict(m, nd, type = "response")) })
  nd1 <- copy(dd); nd1[, z := z + 1]; c(r, rd_per_SD = mean(predict(m, nd1, type = "response")) - mean(predict(m, dd, type = "response"))) }
set.seed(2026); B <- 200; out <- list()
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; est <- stdrisk(d, cd)
  bs <- t(sapply(seq_len(B), function(b) stdrisk(d[sample.int(nrow(d), replace = TRUE)], cd)))
  out[[i]] <- data.table(disease = dis$disease[i], quantity = names(est), proWCdelta_cm = c(pts * sdd, NA), estimate = est,
                         lo = apply(bs, 2, quantile, 0.025), hi = apply(bs, 2, quantile, 0.975), n = nrow(d), events = sum(d[[cd]]))
  say("%s: 10th %.1f%%, 50th %.1f%%, 90th %.1f%%; +2 SD %.1f%%; RD per SD %.2f points", dis$disease[i], 100 * est[["p10"]], 100 * est[["p50"]], 100 * est[["p90"]], 100 * est[["+2 SD"]], 100 * est[["rd_per_SD"]]) }
t69 <- rbindlist(out); fwrite(t69, file.path(TB, "T69_standardised_absolute_risk_by_proWCdelta.csv"))

## (2) calibration and subgroup accuracy of the out-of-fold model (all 52,879)
a <- all0[!is.na(WC) & !is.na(pWC)]
a[, dec := cut(pWC, quantile(pWC, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
cal <- a[, .(n = .N, mean_predicted = mean(pWC), mean_measured = mean(WC)), by = dec][order(dec)]; print(cal)
acc <- function(s, grp, lev) data.table(grouping = grp, subgroup = lev, n = nrow(s), R2 = 1 - sum((s$WC - s$pWC)^2) / sum((s$WC - mean(s$WC))^2), RMSE = sqrt(mean((s$WC - s$pWC)^2)),
                                       mean_error = mean(s$pWC - s$WC), calibration_slope = coef(lm(WC ~ pWC, s))[[2]])
a[, `:=`(sexg = fifelse(Sex == 1, "men", "women"), ageg = cut(Age, c(0, 50, 60, 100), right = FALSE, labels = c("<50", "50-59", ">=60")),
         bmig = cut(BMI, c(0, 25, 30, 100), right = FALSE, labels = c("<25", "25-29.9", ">=30")), eth = fifelse(white == 1, "White", "other"),
         region = fifelse(as.integer(centre) %in% c(11004, 11005, 11003, 11022, 11023), "Scotland and Wales", "England"), miss = fifelse(pna > 0.2, ">20% proteins missing", "<=20% proteins missing"))]
t70 <- rbindlist(c(list(acc(a, "all", "all")), lapply(c("sexg", "ageg", "bmig", "eth", "region", "miss"), function(g) rbindlist(lapply(sort(unique(na.omit(a[[g]]))), function(l) acc(a[get(g) == l], g, as.character(l)))))))
t70 <- rbind(t70, cal[, .(grouping = "tenth of predicted WC", subgroup = as.character(dec), n, R2 = NA, RMSE = NA, mean_error = mean_predicted - mean_measured, calibration_slope = NA)], fill = TRUE)
print(t70); fwrite(t70, file.path(TB, "T70_calibration_and_subgroup_accuracy.csv"))

## (3) share of the log odds ratio explained (Figure 5F estimates, n = 40,364)
g <- fread(file.path(S, "source_table_Figure5F_proWCdelta_beyond_fat_and_clinical_markers.csv"))
t71 <- dcast(g, disease ~ model, value.var = "OR"); setnames(t71, c("disease", "body_fat", "body_fat_clinical", "clinical", "base"))
t71[, `:=`(pct_logOR_removed_clinical = 100 * (1 - log(clinical) / log(base)), pct_logOR_removed_body_fat_clinical = 100 * (1 - log(body_fat_clinical) / log(base)))]
print(t71); fwrite(t71, file.path(TB, "T71_share_of_logOR_explained.csv")); say("DONE")
