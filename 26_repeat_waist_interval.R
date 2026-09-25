## Interval between the baseline assessment and each repeat waist measurement (Results: median 10.1 years at the imaging visit).
suppressPackageStartupMessages(library(data.table))
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
BSZ <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/6-1Body_size_measures_participant.csv"
REC <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/2Recruitment.csv"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12,"lasso_WC.csv")); setnames(lasso, names(lasso)[1], "participant_id")
phen <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","WC","Age","Sex")); setnames(phen,"Participant ID","participant_id")
ht <- fread(BSZ, select=c("Participant ID","Waist circumference | Instance 1","Waist circumference | Instance 2","Waist circumference | Instance 3"))
setnames(ht, c("participant_id","wc1","wc2","wc3"))
rec <- fread(REC, select=c("Participant ID", paste0("Date of attending assessment centre | Instance ", 0:3)))
setnames(rec, c("participant_id","d0","d1","d2","d3"))
d <- merge(merge(phen, lasso[, .(participant_id)], by="participant_id"), ht, by="participant_id")
d <- merge(d, rec, by="participant_id", all.x=TRUE)
d[, `:=`(WC=num(WC), wc1=num(wc1), wc2=num(wc2), wc3=num(wc3))]
for (v in c("d0","d1","d2","d3")) d[[v]] <- as.IDate(d[[v]])
out <- rbindlist(lapply(c("wc1","wc2","wc3"), function(v) {
  dv <- sub("wc","d",v); s <- d[!is.na(WC) & !is.na(get(v)) & !is.na(Age)]
  yrs <- as.numeric(s[[dv]] - s$d0)/365.25
  data.table(instance=v, n=nrow(s), n_with_dates=sum(!is.na(yrs)),
             median_years=round(median(yrs, na.rm=TRUE),1), q1=round(quantile(yrs,.25,na.rm=TRUE),1), q3=round(quantile(yrs,.75,na.rm=TRUE),1),
             min=round(min(yrs,na.rm=TRUE),1), max=round(max(yrs,na.rm=TRUE),1),
             mean_years=round(mean(yrs,na.rm=TRUE),1))
}))
print(out)
fwrite(out, file.path("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data", "source_table_Results_repeat_waist_interval.csv"))
