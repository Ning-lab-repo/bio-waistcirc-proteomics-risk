suppressPackageStartupMessages({library(data.table); library(splines)})
f <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure8/wc_pro_Batch.csv"
hdr <- names(fread(f, nrows=0)); pe <- match("BA", hdr)-1; cols <- hdr[2:pe]
d0 <- fread(f, select=c(hdr[1], cols, "WC","Age","Sex")); setnames(d0,1,"id")
cat("wc_pro_Batch: NA WC", sum(is.na(d0$WC)), " NA Age", sum(is.na(d0$Age)), " NA Sex", sum(is.na(d0$Sex)), " rows", nrow(d0), "\n")
d0[, pna := rowMeans(is.na(as.matrix(d0[, ..cols])))]
m <- d0[, .(id, pna)]
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12,"lasso_WC.csv")); setnames(lasso,1,"id")
phen <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","Sex","Age","WC","BMI","yu_ten_need_diagnosis")); setnames(phen,"Participant ID","id")
cov <- fread(file.path(F7,"complete_data_imputed.csv")); setnames(cov,1,"id"); cov <- cov[, .(id, tdi=num(get("Townsend deprivation index at recruitment")), smoking=as.factor(get("Smoking status")))]
d <- merge(merge(merge(phen, lasso[, .(id, pWC=pred_WC, proWC=BioX_Adjusted, dlt=BioX_Delta)], by="id"), cov, by="id", all.x=TRUE), m, by="id")
d[, `:=`(Sex=as.integer(Sex), Age=num(Age), WC=num(WC), BMI=num(BMI))]
cat("analysis set n", nrow(d), "; participants with >20% missing proteins:", sum(d$pna>.2), " >50%:", sum(d$pna>.5), "\n")
r2 <- function(s) 1 - sum((s$WC-s$pWC)^2)/sum((s$WC-mean(s$WC))^2)
cat(sprintf("OOF R2 all %.3f; <=20%% missing %.3f; >20%% missing %.3f; >50%% missing %.3f\n", r2(d), r2(d[pna<=.2]), r2(d[pna>.2]), r2(d[pna>.5])))
## (an earlier comparison of groups defined by clinical waist thresholds, not reported in the paper, was removed)
fwrite(data.table(protein=cols, missing_fraction=sapply(cols, function(v) mean(is.na(d0[[v]])))), "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data/source_table_Methods_protein_missingness_by_protein.csv")
