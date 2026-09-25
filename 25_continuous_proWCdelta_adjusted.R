suppressPackageStartupMessages({library(data.table); library(splines)})
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
z <- readRDS(file.path(W,"tables",".t1cache.rds")); d <- z$d; x <- z$x
has <- function(code) as.integer(grepl(paste0("(^|[|])",code,"([.]|[|]|$)"), x))
dis <- data.table(disease=c("Type 2 diabetes","Obesity","Dyslipidemia","Hypertension",
                            "Ischaemic heart disease","Heart failure","Liver disease","Chronic kidney disease"),
                  code=c("E11","E66","E78","I10","I25","I50","K76","N18"))
for (i in seq_len(nrow(dis))) d[, (paste0("o_",dis$code[i])) := has(dis$code[i])]
d[, zd := scale(proWCd)[,1]]
dd <- d[!is.na(tdi)&!is.na(smoking)&!is.na(BMI)&!is.na(WC)&!is.na(Age)&!is.na(zd)]
cat("n =", nrow(dd), "\n\n")
run <- function(rhs, lab) rbindlist(lapply(seq_len(nrow(dis)), function(i){
  s <- summary(glm(as.formula(paste0("o_",dis$code[i]," ~ zd + ",rhs)), data=dd, family=binomial()))$coefficients
  data.table(model=lab, disease=dis$disease[i], OR=exp(s["zd",1]),
             lo=exp(s["zd",1]-1.96*s["zd",2]), hi=exp(s["zd",1]+1.96*s["zd",2]), p=s["zd",4]) }))
M1 <- "Age + factor(Sex) + tdi + smoking"
M3 <- "Age + factor(Sex) + tdi + smoking + ns(WC,3) + ns(BMI,3)"
r <- rbind(run(M1,"M1 age,sex,TDI,smoking"), run(M3,"M3 + WC spline + BMI spline"))
r[, p_fdr := p.adjust(p, method="BH"), by=model]
r[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]
cat("=== odds ratio per SD of continuous proWC-delta ===\n")
print(dcast(r, disease ~ model, value.var="txt"))
cat("\nall BH-FDR < 0.05 in the fully adjusted model: ", all(r[model!="M1 age,sex,TDI,smoking", p_fdr] < 0.05), "\n")
print(r[model!="M1 age,sex,TDI,smoking", .(disease, p_fdr=signif(p_fdr,2))])
fwrite(r, file.path(W,"tables","T15_continuous_proWCdelta_WC_BMI_adjusted.csv"))
fwrite(r, "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data/source_table_Results_continuous_proWCdelta_adjusted.csv")
