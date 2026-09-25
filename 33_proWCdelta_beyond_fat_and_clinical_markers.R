## Does proWCdelta carry risk information beyond measured body fat and routine clinical markers?
## Continuous proWCdelta (per SD of the full cohort), 10-year incident endpoints, logistic regression,
## all models fitted in the same complete-case sample.
suppressPackageStartupMessages({library(data.table); library(splines)})
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"
OUT <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/tables"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12,"lasso_WC.csv")); setnames(lasso,1,"id")
phen <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","Sex","Age","WC","BMI","yu_ten_need_diagnosis","HBA1C","HDL","LDLD","TRIG","CRP"))
setnames(phen, c("id","Sex","Age","WC","BMI","dx10","hba1c","hdl","ldl","tg","crp"))
cov <- fread(file.path(F7,"complete_data_imputed.csv")); setnames(cov,1,"id"); cov <- cov[, .(id, tdi=num(get("Townsend deprivation index at recruitment")), smoking=as.factor(get("Smoking status")))]
bc <- fread(file.path(RAW,"bodycomp_subset.csv"), select=c("Participant ID","Body fat percentage | Instance 0","Trunk fat percentage | Instance 0")); setnames(bc, c("id","bfp","tfp"))
fr <- fread(file.path(RAW,"framingham_inputDATA.csv"), select=c("Participant ID","SBP_auto_average","Medication for cholesterol, blood pressure or diabetes | Instance 0","Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id","sbp","med_m","med_f")); fr[, med := paste(med_m, med_f)]
fr[, `:=`(lipid_med=as.integer(grepl("(^|[| ])1([| ]|$)", med)), bp_med=as.integer(grepl("(^|[| ])2([| ]|$)", med)), insulin=as.integer(grepl("(^|[| ])3([| ]|$)", med)))]
d <- Reduce(function(a,b) merge(a,b,by="id",all.x=TRUE), list(merge(phen, lasso[, .(id, proWC=BioX_Adjusted, dlt=BioX_Delta)], by="id"), cov, bc, fr[, .(id, sbp, lipid_med, bp_med, insulin)]))
for (v in c("Age","WC","BMI","hba1c","hdl","ldl","tg","crp","bfp","tfp","sbp")) d[[v]] <- num(d[[v]])
d[, Sex := as.integer(Sex)]
sdd <- sd(d$dlt); d[, z := dlt/sdd]
cat("full cohort n", nrow(d), " SD(proWCdelta)", round(sdd,2), "\n")
x <- as.character(d$dx10); x[is.na(x)] <- ""
dis <- data.table(disease=c("Type 2 diabetes","Obesity","Dyslipidemia","Hypertension","Ischaemic heart disease","Heart failure","Liver disease","Chronic kidney disease"), code=c("E11","E66","E78","I10","I25","I50","K76","N18"))
for (i in seq_len(nrow(dis))) d[[dis$code[i]]] <- as.integer(grepl(paste0("(^|[|])",dis$code[i]), x))
need <- c("Age","Sex","tdi","smoking","WC","BMI","bfp","tfp","hba1c","hdl","ldl","tg","crp","sbp","lipid_med","bp_med","insulin")
cc <- d[complete.cases(d[, ..need]) & tg > 0 & crp > 0]
cat("complete-case sample n", nrow(cc), " (", round(100*nrow(cc)/nrow(d),1), "% of cohort)\n")
base <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
fat  <- "ns(bfp,3) + ns(tfp,3)"
clin <- "ns(hba1c,3) + ns(hdl,3) + ns(ldl,3) + ns(log(tg),3) + ns(log(crp),3) + ns(sbp,3) + lipid_med + bp_med + insulin"
models <- list("WC and BMI"=base, "+ body fat"=paste(base,"+",fat), "+ clinical markers"=paste(base,"+",clin), "+ body fat + clinical markers"=paste(base,"+",fat,"+",clin))
res <- list()
for (i in seq_len(nrow(dis))) for (mn in names(models)) {
  f <- as.formula(paste(dis$code[i], "~ z +", models[[mn]]))
  s <- summary(glm(f, cc, family=binomial()))$coefficients["z",]
  res[[length(res)+1]] <- data.table(disease=dis$disease[i], model=mn, n=nrow(cc), events=sum(cc[[dis$code[i]]]), OR=exp(s[1]), lo=exp(s[1]-1.96*s[2]), hi=exp(s[1]+1.96*s[2]), p=s[4]) }
res <- rbindlist(res); res[, model := factor(model, levels=names(models))]
print(dcast(res[, .(disease, model, txt=sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease ~ model, value.var="txt"))
fwrite(res, file.path(OUT,"T24_proWCdelta_beyond_fat_and_clinical_markers.csv"))
## how much of WC and of proWCdelta do routine measurements capture?
r2 <- function(f) summary(lm(as.formula(f), cc))$r.squared
desc <- data.table(quantity=c("WC ~ age + sex","WC ~ age + sex + clinical markers","proWCdelta ~ WC + BMI + age + sex","proWCdelta ~ WC + BMI + age + sex + body fat","proWCdelta ~ WC + BMI + age + sex + body fat + clinical markers"),
  R2=c(r2("WC ~ Age + Sex"), r2(paste("WC ~ Age + Sex +", clin)), r2("dlt ~ ns(WC,3) + ns(BMI,3) + Age + Sex"), r2(paste("dlt ~ ns(WC,3) + ns(BMI,3) + Age + Sex +", fat)), r2(paste("dlt ~ ns(WC,3) + ns(BMI,3) + Age + Sex +", fat, "+", clin))))
print(desc); fwrite(desc, file.path(OUT,"T28_variance_captured_by_routine_measures.csv"))
## the phenotype contrast (normal WC + high proWC vs normal on both) under the same adjustments
cc[, thr := fifelse(Sex==0,88,102)]; g <- cc[WC < thr]; g[, hi := as.integer(proWC >= thr)]
ph <- list()
for (i in seq_len(nrow(dis))) for (mn in names(models)) {
  s <- summary(glm(as.formula(paste(dis$code[i], "~ hi +", models[[mn]])), g, family=binomial()))$coefficients["hi",]
  ph[[length(ph)+1]] <- data.table(disease=dis$disease[i], model=mn, n=nrow(g), n_high=sum(g$hi), OR=exp(s[1]), lo=exp(s[1]-1.96*s[2]), hi=exp(s[1]+1.96*s[2])) }
ph <- rbindlist(ph)
print(dcast(ph[, .(disease, model, txt=sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease ~ factor(model, levels=names(models)), value.var="txt"))
fwrite(ph, file.path(OUT,"T27_phenotype_beyond_fat_and_clinical_markers.csv"))
