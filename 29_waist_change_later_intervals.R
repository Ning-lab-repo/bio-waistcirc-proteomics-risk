## Is the association of baseline proWCdelta with later waist change explained by baseline measurement error?
## proWCdelta loads on the baseline (instance 0) measurement error. Change between two LATER measurements
## (instance 1 -> 2, 1 -> 3, 2 -> 3) does not contain the baseline error, so an association with it cannot be
## produced by baseline measurement error.
suppressPackageStartupMessages(library(data.table))
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
BSZ <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/6-1Body_size_measures_participant.csv"
REC <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/2Recruitment.csv"
num <- function(x) suppressWarnings(as.numeric(x))
l <- fread(file.path(F12,"lasso_WC.csv")); setnames(l,1,"id")
ph <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","WC","Age","Sex")); setnames(ph,"Participant ID","id")
ht <- fread(BSZ, select=c("Participant ID", paste0("Waist circumference | Instance ", 1:3))); setnames(ht, c("id","wc1","wc2","wc3"))
rec <- fread(REC, select=c("Participant ID", paste0("Date of attending assessment centre | Instance ", 0:3))); setnames(rec, c("id","d0","d1","d2","d3"))
d <- merge(merge(merge(ph, l[, .(id, dlt=BioX_Delta)], by="id"), ht, by="id"), rec, by="id")
d[, `:=`(wc0=num(WC), wc1=num(wc1), wc2=num(wc2), wc3=num(wc3), Age=num(Age))]
for (v in c("d0","d1","d2","d3")) d[[v]] <- as.IDate(d[[v]])
sdd <- sd(d$dlt)
res <- list()
fitit <- function(a, b, label) {
  s <- d[!is.na(get(paste0("wc",a))) & !is.na(get(paste0("wc",b)))]
  s[, dWC := get(paste0("wc",b)) - get(paste0("wc",a))]
  s[, yrs := as.numeric(get(paste0("d",b)) - get(paste0("d",a)))/365.25]
  s[, base := get(paste0("wc",a))]
  m <- lm(dWC ~ dlt + base + Age + Sex + yrs, s)
  ci <- confint(m)["dlt",]; p <- summary(m)$coefficients["dlt",4]
  res[[length(res)+1]] <<- data.table(interval=label, n=nrow(s), median_years=round(median(s$yrs),1), mean_change_cm=round(mean(s$dWC),2),
     beta_per_cm=coef(m)["dlt"], lo=ci[1], hi=ci[2], p=p, beta_per_SD=coef(m)["dlt"]*sdd,
     contains_baseline_error=(a==0))
}
fitit(0,1,"baseline -> first repeat (instance 1)")
fitit(0,2,"baseline -> imaging visit (instance 2)")
fitit(0,3,"baseline -> repeat imaging (instance 3)")
fitit(1,2,"first repeat -> imaging visit")
fitit(1,3,"first repeat -> repeat imaging")
fitit(2,3,"imaging visit -> repeat imaging")
r <- rbindlist(res); print(r[, .(interval, n, median_years, mean_change_cm, b=round(beta_per_cm,3), lo=round(lo,3), hi=round(hi,3), p=signif(p,2), perSD=round(beta_per_SD,2), contains_baseline_error)])
## reliability of a single waist measurement from instance 0 vs 1 residual correlation is not identifiable separately from true change; report correlation only
cat("cor(WC0, WC1):", round(d[!is.na(wc0)&!is.na(wc1), cor(wc0,wc1)],3), "  cor(WC2, WC3):", round(d[!is.na(wc2)&!is.na(wc3), cor(wc2,wc3)],3), "\n")
fwrite(r, "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data/source_table_Results_waist_change_later_intervals.csv")
