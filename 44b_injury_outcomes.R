## 44b_injury_outcomes.R
## Injury outcomes with complications of medical care removed (T80-T98 and Y40-Y98 are largely recorded during
## hospital care), per SD of proWCdelta at fixed WC and BMI; same sample and covariates as 44_confounding_technical_holdout.R.
suppressPackageStartupMessages({ library(data.table); library(splines) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"
num <- function(x) suppressWarnings(as.numeric(x))
lasso <- fread(file.path(F12, "lasso_WC.csv")); setnames(lasso, 1, "id")
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis")); setnames(phen, c("id", "Sex", "Age", "WC", "BMI", "dx10"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id")
cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- merge(merge(phen, lasso[, .(id, dlt = BioX_Delta)], by = "id"), cov, by = "id", all.x = TRUE)
for (v in c("Age", "WC", "BMI")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d[, z := dlt / sd(dlt)]; d <- d[complete.cases(d[, .(Age, Sex, tdi, smoking, WC, BMI)])]
x <- as.character(d$dx10); x[is.na(x)] <- ""
codes <- strsplit(x, "|", fixed = TRUE)
anyc <- function(f) vapply(codes, function(v) any(f(v)), logical(1))
d[, inj_all := as.integer(anyc(function(v) grepl("^[ST]", v)))]
d[, inj_acc := as.integer(anyc(function(v) grepl("^S", v) | grepl("^T([0-6][0-9]|7[0-9])", v)))]
d[, ext_acc := as.integer(anyc(function(v) grepl("^[VWX]", v)))]
d[, fracture := as.integer(anyc(function(v) grepl("^S[0-9]2", v) | grepl("^T(02|08|10|12)", v)))]
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
labs <- c(inj_all = "Injury and poisoning, S00-T98", inj_acc = "Injury and poisoning excluding complications of care, S00-T79",
          ext_acc = "Accidents and other external causes, V01-X59", fracture = "Fracture")
out <- rbindlist(lapply(names(labs), function(y) { s <- summary(glm(as.formula(paste(y, "~ z +", BASE)), d, family = binomial()))$coefficients["z", ]
  data.table(outcome = labs[[y]], n = nrow(d), events = sum(d[[y]]), OR = exp(s[[1]]), lo = exp(s[[1]] - 1.96 * s[[2]]), hi = exp(s[[1]] + 1.96 * s[[2]]), p = s[[4]]) }))
print(out[, .(outcome, events, txt = sprintf("%.2f (%.2f-%.2f)", OR, lo, hi), p = signif(p, 2))])
fwrite(out, file.path(W, "tables", "T44_injury_outcomes.csv"))
cat("DONE\n")
