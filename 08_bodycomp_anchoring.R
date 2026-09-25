## Does proWCdelta correspond to real body composition at the same waist and BMI?
suppressPackageStartupMessages({library(data.table); library(splines); library(MatchIt)})
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
F7  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
BC  <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/bodycomp_subset.csv"
BSZ <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata/6-1Body_size_measures_participant.csv"
num <- function(x) suppressWarnings(as.numeric(x))

lasso <- fread(file.path(F12,"lasso_WC.csv")); setnames(lasso, names(lasso)[1],"participant_id")
lasso <- lasso[, .(participant_id, proWC=BioX_Adjusted, proWCd=BioX_Delta)]
phen  <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","Sex","Age","WC","BMI","HC"))
setnames(phen,"Participant ID","participant_id")
cov <- fread(file.path(F7,"complete_data_imputed.csv")); setnames(cov, names(cov)[1],"participant_id")
cov <- cov[, .(participant_id, tdi=num(get("Townsend deprivation index at recruitment")), smoking=as.factor(get("Smoking status")))]
ht  <- fread(BSZ, select=c("Participant ID","Standing height | Instance 0")); setnames(ht, c("participant_id","height0"))
bc  <- fread(BC); setnames(bc, names(bc)[1], "participant_id")

d <- merge(merge(merge(phen, lasso, by="participant_id"), cov, by="participant_id", all.x=TRUE), ht, by="participant_id", all.x=TRUE)
d <- merge(d, bc, by="participant_id", all.x=TRUE)
d[, `:=`(Sex=as.integer(Sex), Age=num(Age), WC=num(WC), BMI=num(BMI), HC=num(HC),
         proWC=num(proWC), proWCd=num(proWCd), height0=num(height0))]
cat(sprintf("cohort merged: %d\n", nrow(d)))

OUT <- list(
 "Body fat percentage"            = "Body fat percentage | Instance 0",
 "Whole body fat mass (kg)"       = "Whole body fat mass | Instance 0",
 "Whole body fat-free mass (kg)"  = "Whole body fat-free mass | Instance 0",
 "Whole body water mass (kg)"     = "Whole body water mass | Instance 0",
 "Trunk fat percentage"           = "Trunk fat percentage | Instance 0",
 "Trunk fat mass (kg)"            = "Trunk fat mass | Instance 0",
 "Trunk fat-free mass (kg)"       = "Trunk fat-free mass | Instance 0",
 "Basal metabolic rate (kJ)"      = "Basal metabolic rate | Instance 0",
 "Impedance of whole body"        = "Impedance of whole body | Instance 0",
 "VAT mass, DXA (kg)"             = "VAT (visceral adipose tissue) mass | Instance 2",
 "VAT volume, DXA (L)"            = "VAT (visceral adipose tissue) volume | Instance 2",
 "Total fat mass, DXA (kg)"       = "Total fat mass | Instance 2",
 "Total lean mass, DXA (kg)"      = "Total lean mass | Instance 2",
 "Trunk fat mass, DXA (kg)"       = "Trunk fat mass | Instance 2",
 "Android fat mass, DXA (kg)"     = "Android fat mass | Instance 2"
# "VAT volume, MRI (L)"            = "Visceral adipose tissue volume (VAT) | Instance 2",   ## REMOVED: MRI block not correctly keyed to participant ID (see 07_audit/README_bodycomp_subset.txt)
# "ASAT volume, MRI (L)"           = "Abdominal subcutaneous adipose tissue volume (ASAT) | Instance 2",   ## REMOVED, same reason
# "Total adipose tissue vol, MRI"  = "Total adipose tissue volume | Instance 2",   ## REMOVED, same reason
# "Total lean tissue vol, MRI"     = "Total lean tissue volume | Instance 2",   ## REMOVED, same reason
# "Liver PDFF (%)"                 = "10P Liver PDFF (proton density fat fraction) | Instance 2",   ## REMOVED, same reason
# "Thigh fat-free muscle vol, MRI" = "Total thigh fat-free muscle volume | Instance 2"   ## REMOVED, same reason
)

## ---- (1) proWCdelta vs body composition, adjusted for WC, BMI, age, sex, height ----
res <- rbindlist(lapply(names(OUT), function(lab){
  v <- OUT[[lab]]; if (!v %in% names(d)) return(NULL)
  dd <- copy(d); dd[, y := num(get(v))]
  dd <- dd[!is.na(y) & !is.na(proWCd) & !is.na(WC) & !is.na(BMI) & !is.na(height0) & !is.na(Age) & !is.na(Sex)]
  if (nrow(dd) < 300) return(data.table(outcome=lab, n=nrow(dd), beta_per_SD=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_))
  dd[, z := as.numeric(scale(y))]
  f <- lm(z ~ proWCd + ns(WC,3) + ns(BMI,3) + Age + factor(Sex) + height0, data=dd)
  co <- summary(f)$coefficients; ci <- confint(f); s <- sd(dd$proWCd)
  data.table(outcome=lab, n=nrow(dd), beta_per_SD=co["proWCd","Estimate"]*s,
             lo=ci["proWCd",1]*s, hi=ci["proWCd",2]*s, p=co["proWCd","Pr(>|t|)"])
}))
res[, p_fdr := p.adjust(p,"BH")]
cat("\n=== (1) proWCdelta vs body composition, adjusted for WC spline, BMI spline, age, sex, height ===\n")
print(res[, .(outcome, n, beta_per_SD=round(beta_per_SD,3), CI=sprintf("%.3f to %.3f", lo, hi), p_fdr=signif(p_fdr,2))])
fwrite(res, file.path(W,"tables","T11_proWCdelta_vs_bodycomposition.csv"))

## ---- (2) matched pairs: same sex, age, WC, BMI ----
d[, thr := fifelse(Sex==0,88,102)]
d[, grp := fifelse(WC<thr & proWC<thr,"N/N", fifelse(WC<thr & proWC>=thr,"N/H",
            fifelse(WC>=thr & proWC<thr,"H/N","H/H")))]
m <- d[grp %in% c("N/N","N/H") & complete.cases(Age,Sex,WC,BMI,tdi,smoking)]
m[, NH := as.integer(grp=="N/H")]
set.seed(123)
mm <- matchit(NH ~ Age + WC + BMI, data=m, method="nearest", exact=~Sex,
              distance="mahalanobis", caliper=c(WC=0.15,BMI=0.15), std.caliper=TRUE, ratio=1)
md <- as.data.table(match.data(mm))
cat(sprintf("\nmatched pairs: %d\n", sum(md$NH==1)))
mres <- rbindlist(lapply(names(OUT), function(lab){
  v <- OUT[[lab]]; if (!v %in% names(md)) return(NULL)
  dd <- copy(md); dd[, y := num(get(v))]; dd <- dd[!is.na(y)]
  if (sum(dd$NH==1) < 50 || sum(dd$NH==0) < 50) return(NULL)
  f <- lm(y ~ NH + Age + factor(Sex), data=dd)
  co <- summary(f)$coefficients; ci <- confint(f)
  data.table(outcome=lab, n_NH=sum(dd$NH==1), n_NN=sum(dd$NH==0),
             mean_NN=mean(dd[NH==0]$y), mean_NH=mean(dd[NH==1]$y),
             diff=co["NH","Estimate"], lo=ci["NH",1], hi=ci["NH",2], p=co["NH","Pr(>|t|)"])
}))
mres[, p_fdr := p.adjust(p,"BH")]
cat("\n=== (2) body composition in WC- and BMI-matched pairs (N/H vs N/N) ===\n")
print(mres[, .(outcome, n_NH, mean_NN=round(mean_NN,2), mean_NH=round(mean_NH,2),
               diff=round(diff,3), CI=sprintf("%.3f to %.3f", lo, hi), p_fdr=signif(p_fdr,2))])
fwrite(mres, file.path(W,"tables","T12_matched_pairs_bodycomposition.csv"))
saveRDS(list(d=d[, .(participant_id, Sex, Age, WC, BMI, height0, proWC, proWCd, grp)], matched_ids=md$participant_id),
        file.path(W,"output","bodycomp_merge_keys.rds"))
cat("\nDONE\n")
