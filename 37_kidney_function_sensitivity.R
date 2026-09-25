## Kidney-function sensitivity: many plasma proteins are cleared by the kidney, so the continuous proWCdelta models are
## additionally adjusted for natural splines of log creatinine and log cystatin C, with and without participants whose
## eGFR (CKD-EPI 2021 creatinine-cystatin C equation) is below 60 mL/min/1.73 m2. Output: T35.
suppressPackageStartupMessages({library(data.table); library(splines)})
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12,"lasso_WC.csv")); setnames(lasso,1,"id")
phen <- fread(file.path(F12,"pro53013_新诊_newdead.csv"), select=c("Participant ID","Sex","Age","WC","BMI","yu_ten_need_diagnosis","CRE","CYS"))
setnames(phen, c("id","Sex","Age","WC","BMI","dx10","cre","cys"))
cov <- fread(file.path(F7,"complete_data_imputed.csv")); setnames(cov,1,"id"); cov <- cov[, .(id, tdi=num(get("Townsend deprivation index at recruitment")), smoking=as.factor(get("Smoking status")))]
d <- merge(merge(phen, lasso[, .(id, dlt=BioX_Delta)], by="id"), cov, by="id")
for (v in c("Age","WC","BMI","cre","cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, cre, cys)]) & cre > 0 & cys > 0]
## CKD-EPI 2021 creatinine-cystatin C equation (Sex 0 = female)
d[, scr := cre/88.4]; d[, `:=`(k=ifelse(Sex==0,0.7,0.9), a=ifelse(Sex==0,-0.219,-0.144))]
d[, egfr := 135 * pmin(scr/k,1)^a * pmax(scr/k,1)^(-0.544) * pmin(cys/0.8,1)^(-0.323) * pmax(cys/0.8,1)^(-0.778) * 0.9961^Age * ifelse(Sex==0,0.963,1)]
d[, z_pro := dlt/sd(dlt)]
cat(sprintf("n %d; median eGFR %.0f; eGFR < 60: %d (%.1f%%); cor(proWCdelta, eGFR) = %.3f; cor(proWCdelta, log cystatin C) = %.3f\n",
  nrow(d), median(d$egfr), d[egfr<60,.N], 100*d[,mean(egfr<60)], cor(d$z_pro, d$egfr), cor(d$z_pro, log(d$cys))))
base <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"; kid <- "ns(log(cre),3) + ns(log(cys),3)"
x <- as.character(d$dx10); x[is.na(x)] <- ""
f1 <- function(dd, cd, extra) { s <- summary(glm(as.formula(paste(cd, "~ z_pro +", base, extra)), dd, family=binomial()))$coefficients["z_pro",]
  sprintf("%.2f (%.2f-%.2f)", exp(s[1]), exp(s[1]-1.96*s[2]), exp(s[1]+1.96*s[2])) }
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
x <- as.character(d$dx10); x[is.na(x)] <- ""
fit1 <- function(dd, cd, extra, label) { s <- summary(glm(as.formula(paste(cd, "~ z_pro +", base, extra)), dd, family=binomial()))$coefficients["z_pro",]
  data.table(model=label, n=nrow(dd), events=sum(dd[[cd]]), OR=exp(s[1]), lo=exp(s[1]-1.96*s[2]), hi=exp(s[1]+1.96*s[2])) }
dis <- c("Type 2 diabetes"="E11","Obesity"="E66","Dyslipidemia"="E78","Hypertension"="I10","Ischaemic heart disease"="I25","Heart failure"="I50","Liver disease"="K76","Chronic kidney disease"="N18")
res <- rbindlist(lapply(names(dis), function(nm) { cd <- dis[[nm]]; d[[cd]] <<- as.integer(grepl(paste0("(^|[|])",cd), x))
  r <- rbind(fit1(d, cd, "", "WC and BMI"), fit1(d, cd, paste("+", kid), "+ creatinine and cystatin C"),
             fit1(d[egfr >= 60], cd, paste("+", kid), "eGFR >= 60 only, + creatinine and cystatin C"))
  r[, disease := nm]; r }))
print(dcast(res[, .(disease, model, v=sprintf("%.2f (%.2f-%.2f)", OR, lo, hi))], disease ~ model, value.var="v"), width=200)
fwrite(res[, .(disease, model, n, events, OR, lo, hi)], file.path(W,"tables","T35_kidney_function_sensitivity.csv"))
cat("DONE\n")
