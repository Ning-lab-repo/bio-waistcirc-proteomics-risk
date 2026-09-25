## 90_same_visit_mri_and_selection.R
## Analyses with the Olink data measured at the imaging visits (UKB-PPP participants of the COVID-19 repeat-imaging
## study: instance 2, first imaging visit; instance 3, repeat imaging visit) and the UKB-PPP sample-selection fields.
##  A. Selection. The published MRI and disease estimates for proWCdelta are reproduced (the script stops otherwise) and
##     then re-estimated without the participants selected by UKB-PPP consortium members ("UKB-PPP Consortium selected
##     participant | Instance 0"), without the repeat-imaging participants, and without both.
##  B. Same visit. A proWC model restricted to the proteins measured at the imaging visits is trained at baseline (LASSO
##     on those proteins with age and sex, lambda.min of ten-fold cross-validation, as for proWC) in participants
##     WITHOUT imaging-visit proteomics, then applied to the repeat-imaging participants' proteins from the first
##     imaging visit. Its discordance (anchored on the waist measured at that visit) is related to visceral and
##     subcutaneous fat on MRI at the same visit, with waist, BMI (same visit), age and sex held fixed, and compared with
##     (i) the same model applied to the same participants' baseline proteins (the decade-long design of the paper),
##     (ii) the out-of-fold proWCdelta of the paper in the same participants, and (iii) the composition of the extra
##     fat of a larger measured waist at the same visit.
##  C. Stability. In participants with proteins at both imaging visits, the discordance at each visit (each anchored on
##     that visit's own tape measurement) is correlated within person: random error in a single tape measurement does
##     not repeat across visits, so the correlation measures the part of the discordance that is a stable trait.
##  D. The imaging-visit and baseline discordances compared with the same covariates (age, sex, WC and BMI at both
##     visits), and the Olink plates of each participant's two imaging-visit samples (added after A-C, which it leaves
##     unchanged; its rows follow theirs in T123).
## Data: as in 74/75/79/89 (figure12 analysis files, bodycompositionInstance.2.csv, 6-1Body_size_measures_participant.csv,
## complete_data_imputed.csv, framingham_inputDATA.csv) plus yuxuan_prowc/protomics (instance 2 and 3 Olink NPX; the
## UKB-PPP selection, plate and protein-count fields).
## Output: tables/T123_same_visit_mri.csv, tables/T124_selection_sensitivity.csv; everything else is printed.
suppressPackageStartupMessages({ library(data.table); library(splines); library(glmnet); library(doParallel) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"; F7 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure7"
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; BC <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
PR  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/protomics"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
set.seed(2026); registerDoParallel(cores = 24)
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
res <- list(); add <- function(...) res[[length(res) + 1]] <<- data.table(...)
sel_res <- list(); add_sel <- function(...) sel_res[[length(sel_res) + 1]] <<- data.table(...)
hdr <- function(f) names(fread(f, nrows = 0))

## ================= 0. what the files hold =================
say("---- imaging-visit anthropometry and MRI fields available")
bs_file <- file.path(RAW, "6-1Body_size_measures_participant.csv"); bsh <- hdr(bs_file)
print(grep("Instance [23]", bsh, value = TRUE))
for (f in c(list.files(BC, full.names = TRUE, pattern = "\\.csv$"), list.files(RAW, full.names = TRUE, pattern = "\\.csv$"))) {
  h <- tryCatch(hdr(f), error = function(e) character(0))
  hit <- grep("(visceral|VAT|subcutaneous|ASAT).*Instance.?3|Instance.?3.*(visceral|VAT)", h, value = TRUE, ignore.case = TRUE)
  if (length(hit)) say("instance-3 MRI fields in %s: %s", basename(f), paste(hit, collapse = " | ")) }
say("---- search for instance-3 MRI fields finished")

## ================= 1. selection fields and imaging-visit proteins =================
sel <- fread(file.path(PR, "1Proteomics_biomarkers_participant.csv"), colClasses = "character"); setnames(sel, 1, "id")
sel[, id := as.integer(id)]
for (v in grep("Consortium selected", names(sel), value = TRUE)) { say("%s:", v); print(table(sel[[v]], useNA = "ifany")) }
for (v in grep("Number of proteins measured", names(sel), value = TRUE)) {
  x <- num(sel[[v]]); say("%s: %d participants; median %s proteins", v, sum(!is.na(x)), median(x, na.rm = TRUE)) }
v0 <- "UKB-PPP Consortium selected participant | Instance 0"; stopifnot(v0 %in% names(sel))
cons_ids <- sel[tolower(get(v0)) %in% c("yes", "1", "true"), id]
say("consortium-selected at baseline: %d participants", length(cons_ids)); stopifnot(length(cons_ids) > 100)

read_inst <- function(i) {
  fs <- sort(list.files(file.path(PR, paste0("Instance", i)), full.names = TRUE, pattern = "\\.csv$")); stopifnot(length(fs) == 8)
  parts <- lapply(fs, function(f) { d <- fread(f, colClasses = "character"); stopifnot(all(d[[1]] == d[[2]]))
    d <- d[, -2, with = FALSE]; setnames(d, 1, "id"); d })
  d <- Reduce(function(a, b) merge(a, b, by = "id"), parts)
  gene <- sub(";.*$", "", names(d)[-1]); stopifnot(!anyDuplicated(gene)); setnames(d, c("id", gene))
  for (g in gene) set(d, j = g, value = num(d[[g]]))
  d[, id := as.integer(id)]; stopifnot(!anyDuplicated(d$id)); d }
i2 <- read_inst(2); i3 <- read_inst(3)
g2 <- setdiff(names(i2), "id"); stopifnot(setequal(g2, setdiff(names(i3), "id")))
say("imaging-visit proteomics: instance 2 %d participants, instance 3 %d, both %d; %d proteins",
    nrow(i2), nrow(i3), length(intersect(i2$id, i3$id)), length(g2))
rep_ids <- union(i2$id, i3$id)     ## every repeat-imaging participant: left out of training and of the "without" subsets
miss2 <- rowMeans(is.na(as.matrix(i2[, ..g2]))); miss3 <- rowMeans(is.na(as.matrix(i3[, ..g2])))
say("share of proteins missing per participant: instance 2 median %.3f (max %.3f; >50%%: %d); instance 3 median %.3f (max %.3f; >50%%: %d)",
    median(miss2), max(miss2), sum(miss2 > 0.5), median(miss3), max(miss3), sum(miss3 > 0.5))
i2 <- i2[miss2 <= 0.5]; i3 <- i3[miss3 <= 0.5]
st1 <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/04_supplementary_tables/SupplementaryTable1_proWC_LASSO_coefficients.csv")
say("Supplementary Table 1: %d rows; columns %s", nrow(st1), paste(names(st1), collapse = " | "))
m1 <- norm(unlist(st1[, 1])) %in% norm(g2); say("proteins of the deposited 1,549-protein model also measured at the imaging visits: %d", sum(m1))

## ================= 2. baseline data and the published estimates (as in 89) =================
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI")); setnames(phen, c("id", "Sex", "Age", "WC", "BMI"))
la <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
coh <- merge(phen, la[, .(id, dlt = BioX_Delta)], by = "id"); for (v in c("Age", "WC", "BMI")) coh[[v]] <- num(coh[[v]]); coh[, Sex := as.integer(Sex)]
coh <- coh[complete.cases(coh[, .(Age, Sex, WC, BMI, dlt)])]; sdz <- sd(coh$dlt); stopifnot(nrow(coh) == 52742, abs(sdz - 5.4703) < 0.001)
say("repeat-imaging participants in the analysed cohort of 52,742: %d; consortium-selected: %d", sum(coh$id %in% rep_ids), sum(coh$id %in% cons_ids))

mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2", "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
im <- Reduce(function(a, b) merge(a, b, by = "id"), list(phen, la[, .(id, dlt = BioX_Delta)], mri))
for (v in c("Age", "WC", "BMI", "vat", "asat")) im[[v]] <- num(im[[v]]); im[, Sex := as.integer(Sex)]
im <- im[!is.na(dlt) & !is.na(WC) & !is.na(BMI) & !is.na(Age) & !is.na(Sex) & !is.na(vat) & !is.na(asat) & asat > 0]
stopifnot(nrow(im) == 6901); im[, z := dlt / sdz]
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
b_rep <- coef(lm(as.formula(paste("vat ~ z +", B0)), im))[["z"]]
say("REPRODUCTION: visceral fat per SD of proWCdelta %.6f (published 0.545079)", b_rep); stopifnot(abs(b_rep - 0.545078589469762) < 1e-8)

phen2 <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI", "yu_ten_need_diagnosis", "more_ten_new_diagnosis", "s_Diagnoses_ICD10", "Noncancer_illne_code_elfreported_Instance0", "HBA1C", "CRE", "CYS"))
setnames(phen2, c("id", "Sex", "Age", "WC", "BMI", "dx10", "dxmore", "dxall", "selfrep", "hba1c", "cre", "cys"))
fr <- fread(file.path(RAW, "framingham_inputDATA.csv"), select = c("Participant ID", "medication_name", "Medication for cholesterol, blood pressure or diabetes | Instance 0", "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance 0"))
setnames(fr, c("id", "medname", "med_m", "med_f"))
cov <- fread(file.path(F7, "complete_data_imputed.csv")); setnames(cov, 1, "id"); cov <- cov[, .(id, tdi = num(get("Townsend deprivation index at recruitment")), smoking = as.factor(get("Smoking status")))]
d <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(merge(phen2, la[, .(id, pWC = pred_WC, dlt = BioX_Delta)], by = "id"), cov, fr))
for (v in c("Age", "WC", "BMI", "hba1c", "cre", "cys")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, tdi, smoking, dlt)])]; sdd <- sd(d$dlt); d[, z := dlt / sdd]; stopifnot(nrow(d) == 52742)
s <- function(v) { v <- tolower(as.character(v)); v[is.na(v)] <- ""; v }
d[, `:=`(dx10 = toupper(s(dx10)), dxmore = toupper(s(dxmore)), dxall = toupper(s(dxall)), selfrep = s(selfrep), med = paste(s(med_m), s(med_f)), medname = s(medname))]
d[, scr := cre / 88.4]; d[, `:=`(kk = ifelse(Sex == 0, 0.7, 0.9), aa = ifelse(Sex == 0, -0.219, -0.144))]
d[, egfr := 135 * pmin(scr / kk, 1)^aa * pmax(scr / kk, 1)^(-0.544) * pmin(cys / 0.8, 1)^(-0.323) * pmax(cys / 0.8, 1)^(-0.778) * 0.9961^Age * ifelse(Sex == 0, 0.963, 1)]
has <- function(x, pat) grepl(pat, x, perl = TRUE)
glucose_drugs <- "metformin|gliclazide|glibenclamide|glipizide|glimepiride|pioglitazone|rosiglitazone|sitagliptin|saxagliptin|linagliptin|vildagliptin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|insulin"
extra <- list(E11 = has(d$selfrep, "(^|\\|)(diabetes|type 2 diabetes|type 1 diabetes|diabetic [a-z ]+)(\\||$)") | has(d$med, "insulin") | has(d$medname, glucose_drugs) | (!is.na(d$hba1c) & d$hba1c >= 48),
  I10 = has(d$selfrep, "(^|\\|)(hypertension|essential hypertension)(\\||$)") | has(d$med, "blood pressure medication"),
  E78 = has(d$selfrep, "(^|\\|)high cholesterol(\\||$)") | has(d$med, "cholesterol lowering medication"),
  I25 = has(d$selfrep, "(^|\\|)(angina|heart attack/myocardial infarction)(\\||$)"),
  I50 = has(d$selfrep, "(^|\\|)heart failure/pulmonary odema(\\||$)"), K76 = has(d$selfrep, "(^|\\|)liver failure/cirrhosis(\\||$)"),
  N18 = has(d$selfrep, "(^|\\|)(renal failure not requiring dialysis|renal failure requiring dialysis|renal/kidney failure)(\\||$)") | (!is.na(d$egfr) & d$egfr < 60), E66 = rep(FALSE, nrow(d)))
dis <- data.table(disease = c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease", "Ischaemic heart disease", "Dyslipidaemia", "Hypertension", "Obesity diagnosis"),
                  code = c("E11", "I50", "K76", "N18", "I25", "E78", "I10", "E66"),
                  published = c(1.82105304412415, 1.50033695791199, 1.43514603957416, 1.29353725419005, NA, NA, NA, NA))
for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; pt <- paste0("(^|[|])", cd); inc <- grepl(pt, d$dx10) | grepl(pt, d$dxmore)
  d[[cd]] <- as.integer(grepl(pt, d$dx10)); d[[paste0("free_", cd)]] <- !((grepl(pt, d$dxall) & !inc) | extra[[cd]]) }
BASE <- "Age + Sex + tdi + smoking + ns(WC,3) + ns(BMI,3)"
or1 <- function(cd, x) { co <- summary(glm(as.formula(paste(cd, "~ z +", BASE)), x, family = binomial()))$coefficients["z", ]; c(or = exp(co[[1]]), lo = exp(co[[1]] - 1.96 * co[[2]]), hi = exp(co[[1]] + 1.96 * co[[2]]), p = co[[4]]) }
for (i in 1:4) { x <- d[get(paste0("free_", dis$code[i])) == TRUE]; r <- or1(dis$code[i], x)[["or"]]
  say("REPRODUCTION: %s OR per SD %.6f (published %.6f)", dis$disease[i], r, dis$published[i]); stopifnot(abs(r - dis$published[i]) < 1e-8) }

## ================= A. selection =================
subsets <- list("all (published)" = function(ids) rep(TRUE, length(ids)),
                "without consortium-selected" = function(ids) !(ids %in% cons_ids),
                "without repeat-imaging participants" = function(ids) !(ids %in% rep_ids),
                "without both" = function(ids) !(ids %in% cons_ids) & !(ids %in% rep_ids))
for (nm in names(subsets)) {
  x <- im[subsets[[nm]](id)]
  for (y in c("vat", "asat")) { co <- summary(lm(as.formula(paste(y, "~ z +", B0)), x))$coefficients["z", ]
    add_sel(subset = nm, outcome = ifelse(y == "vat", "visceral adipose tissue (L)", "abdominal subcutaneous adipose tissue (L)"), n = nrow(x), events = NA_integer_,
            estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) }
  for (i in seq_len(nrow(dis))) { cd <- dis$code[i]; x <- d[get(paste0("free_", cd)) == TRUE][subsets[[nm]](id)]; r <- or1(cd, x)
    add_sel(subset = nm, outcome = dis$disease[i], n = nrow(x), events = sum(x[[cd]]), estimate = r[["or"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]]) } }
st <- rbindlist(sel_res); fwrite(st, file.path(TB, "T124_selection_sensitivity.csv"))
say("---- A. selection (visceral fat in L per SD; odds ratios per SD)")
print(st[, .(subset, outcome, n, events, est = round(estimate, 3), lo = round(lo, 3), hi = round(hi, 3))], nrows = 100)

## ================= B. same visit =================
dat <- fread(file.path(F12, "analysis_data_WC.csv")); setnames(dat, 1, "id")
prot0 <- setdiff(names(dat), c("id", "WC", "Age", "Sex")); stopifnot(length(prot0) == 2920, !anyDuplicated(norm(prot0)))
m <- match(norm(g2), norm(prot0)); say("imaging-visit proteins matched to baseline columns: %d of %d", sum(!is.na(m)), length(g2))
if (any(is.na(m))) say("not matched: %s", paste(g2[is.na(m)], collapse = ", "))
cm <- data.table(gene = g2[!is.na(m)], base = prot0[m[!is.na(m)]])
tr <- dat[!(id %in% rep_ids) & !is.na(WC)]
X <- as.matrix(tr[, cm$base, with = FALSE]); med0 <- apply(X, 2, median, na.rm = TRUE)
say("missing values in the baseline training matrix: %d", sum(is.na(X)))
for (j in which(colSums(is.na(X)) > 0)) X[is.na(X[, j]), j] <- med0[j]
X <- cbind(X, Age = num(tr$Age), Sex = num(tr$Sex))
grid <- 10^seq(0.2, -4, length.out = 100)
cvf <- cv.glmnet(X, tr$WC, alpha = 1, lambda = grid, nfolds = 10, parallel = TRUE)
bb <- as.matrix(coef(cvf, s = "lambda.min")); nz <- setdiff(rownames(bb)[bb[, 1] != 0], "(Intercept)")
r2cv <- 1 - cvf$cvm[cvf$lambda == cvf$lambda.min] / mean((tr$WC - mean(tr$WC))^2)
say("restricted model: trained in %d baseline participants without imaging-visit proteomics; lambda.min %.5f; %d proteins retained (age %s, sex %s); cross-validated R2 %.3f",
    nrow(tr), cvf$lambda.min, length(setdiff(nz, c("Age", "Sex"))), "Age" %in% nz, "Sex" %in% nz, r2cv)
add(part = "B", measure = "restricted model: cross-validated R2 at baseline", sample = "baseline, without repeat-imaging participants", n = nrow(tr), estimate = r2cv, lo = NA_real_, hi = NA_real_, p = NA_real_)
add(part = "B", measure = "restricted model: proteins retained", sample = "", n = nrow(tr), estimate = length(setdiff(nz, c("Age", "Sex"))), lo = NA_real_, hi = NA_real_, p = NA_real_)
predict_at <- function(P, age, sex) { M <- as.matrix(P[, cm$gene, with = FALSE]); md <- apply(M, 2, median, na.rm = TRUE)
  for (j in which(colSums(is.na(M)) > 0)) M[is.na(M[, j]), j] <- md[j]
  as.numeric(predict(cvf, newx = cbind(M, Age = age, Sex = sex), s = "lambda.min")) }

## ages and sizes at the imaging visits (searched in the raw extracts; the file used is printed)
cand <- c(file.path(RAW, "BioX_total_full_with_age_months.csv"), file.path(RAW, "1Baseline_characteristics_participant.csv"),
          file.path(RAW, "2Recruitment.csv"), bs_file, list.files(RAW, full.names = TRUE, pattern = "\\.csv$"))
cand <- unique(cand[file.exists(cand)]); H <- lapply(cand, function(f) tryCatch(hdr(f), error = function(e) character(0)))
find_col <- function(label) { for (k in seq_along(cand)) if (label %in% H[[k]]) {
    x <- fread(cand[k], select = c("Participant ID", label)); setnames(x, c("id", "v")); say("%s: from %s", label, basename(cand[k])); return(x) }
  NULL }
age_at <- function(i) {
  x <- find_col(paste0("Age when attended assessment centre | Instance ", i)); if (!is.null(x)) { x[, v := num(v)]; setnames(x, "v", "age"); return(x) }
  dt <- find_col(paste0("Date of attending assessment centre | Instance ", i)); yb <- find_col("Year of birth"); mb <- find_col("Month of birth")
  if (is.null(dt) || is.null(yb) || is.null(mb)) stop(sprintf("age at instance %d not found", i))
  x <- Reduce(function(a, b) merge(a, b, by = "id"), list(setnames(dt, "v", "date"), setnames(yb, "v", "yb"), setnames(mb, "v", "mb")))
  mnum <- ifelse(is.na(suppressWarnings(as.integer(x$mb))), match(x$mb, month.name), suppressWarnings(as.integer(x$mb)))
  x[, age := as.numeric(as.Date(date) - as.Date(sprintf("%s-%02d-15", yb, mnum))) / 365.25]
  say("age at instance %d computed from the visit date and the year and month of birth", i); x[, .(id, age)] }
a2 <- age_at(2); a3 <- tryCatch(age_at(3), error = function(e) { say("%s", conditionMessage(e)); NULL })
get_size <- function(i) { w <- find_col(paste0("Waist circumference | Instance ", i)); b <- find_col(paste0("Body mass index (BMI) | Instance ", i))
  if (is.null(w) || is.null(b)) return(NULL)
  x <- merge(setnames(w, "v", "wc"), setnames(b, "v", "bmi"), by = "id"); x[, `:=`(wc = num(wc), bmi = num(bmi))]
  setnames(x, c("wc", "bmi"), paste0(c("wc", "bmi"), i)); x }
bs <- get_size(2); stopifnot(!is.null(bs)); s3 <- get_size(3)
if (!is.null(s3)) bs <- merge(bs, s3, by = "id", all.x = TRUE) else say("waist or BMI at instance 3 not found: section C skipped")
say("waist and BMI recorded: instance 2 %d; instance 3 %d", bs[!is.na(wc2) & !is.na(bmi2), .N], if (is.null(s3)) 0L else bs[!is.na(wc3) & !is.na(bmi3), .N])

## the first imaging visit: proteins, MRI, size and age on the same day, and the same participants at baseline
base <- merge(dat[, .(id, WC0 = WC, Age0 = num(Age), Sex = num(Sex))], phen[, .(id, BMI0 = num(BMI))], by = "id")
sv <- Reduce(function(a, b) merge(a, b, by = "id"), list(i2[, .(id)], base, a2, bs[, .(id, wc2, bmi2)], mri, la[, .(id, dlt = BioX_Delta)]))
for (v in c("vat", "asat")) sv[[v]] <- num(sv[[v]])
sv <- sv[complete.cases(sv[, .(WC0, Age0, Sex, BMI0, age, wc2, bmi2, vat, asat, dlt)]) & asat > 0]
say("same-visit analysis set: %d participants (proteins, MRI, waist, BMI and age at the first imaging visit, and baseline data)", nrow(sv))
sv[, pwc2 := predict_at(i2[match(sv$id, i2$id)], sv$age, sv$Sex)]
X0 <- as.matrix(dat[match(sv$id, dat$id), cm$base, with = FALSE]); for (j in which(colSums(is.na(X0)) > 0)) X0[is.na(X0[, j]), j] <- med0[j]
sv[, pwc0 := as.numeric(predict(cvf, newx = cbind(X0, Age = sv$Age0, Sex = sv$Sex), s = "lambda.min"))]
sv[, `:=`(d2 = resid(lm(pwc2 ~ wc2)), d0 = resid(lm(pwc0 ~ WC0)))]
say("in these participants: r(pWC at the imaging visit, waist then) %.3f; r(pWC at baseline, baseline waist) %.3f; SD of the discordance %.2f cm (imaging visit) and %.2f cm (baseline)",
    cor(sv$pwc2, sv$wc2), cor(sv$pwc0, sv$WC0), sd(sv$d2), sd(sv$d0))
sv[, `:=`(z2 = d2 / sd(d2), z0 = d0 / sd(d0), zp = dlt / sd(dlt), zw = wc2 / sd(wc2))]
F_SV <- "age + factor(Sex) + ns(wc2,3) + ns(bmi2,3)"; F_B0 <- "Age0 + factor(Sex) + ns(WC0,3) + ns(BMI0,3)"
fits <- list(
  "same visit: imaging-visit proteins, size and age at that visit"                = list(t = "z2", f = F_SV),
  "decade apart: baseline proteins (same restricted model), baseline size"          = list(t = "z0", f = F_B0),
  "decade apart: baseline proteins, baseline and imaging-visit size"              = list(t = "z0", f = paste(F_B0, "+ ns(wc2,3) + ns(bmi2,3)")),
  "decade apart: published out-of-fold proWCdelta, baseline size"                 = list(t = "zp", f = F_B0),
  "comparator: larger measured waist at the imaging visit, BMI fixed"              = list(t = "zw", f = "age + factor(Sex) + ns(bmi2,3)"))
for (nm in names(fits)) { t <- fits[[nm]]$t; f <- fits[[nm]]$f
  for (y in c("vat", "asat")) { co <- summary(lm(as.formula(paste(y, "~", t, "+", f)), sv))$coefficients[t, ]
    add(part = "B", measure = paste(ifelse(y == "vat", "visceral adipose tissue (L) per SD", "abdominal subcutaneous adipose tissue (L) per SD"), "-", nm), sample = "first imaging visit", n = nrow(sv),
        estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) }
  r2b <- summary(lm(as.formula(paste("vat ~", f)), sv))$r.squared; r2 <- summary(lm(as.formula(paste("vat ~", t, "+", f)), sv))$r.squared
  add(part = "B", measure = paste("partial R2 for visceral fat -", nm), sample = "first imaging visit", n = nrow(sv), estimate = (r2 - r2b) / (1 - r2b), lo = NA_real_, hi = NA_real_, p = NA_real_) }
for (nm in names(fits)[c(1, 2, 5)]) { t <- fits[[nm]]$t; f <- sub("factor\\(Sex\\) \\+ ", "", fits[[nm]]$f)
  for (sx in list(list("women", 0), list("men", 1))) { x <- sv[Sex == sx[[2]]]
    for (y in c("vat", "asat")) { co <- summary(lm(as.formula(paste(y, "~", t, "+", f)), x))$coefficients[t, ]
      add(part = "B", measure = paste(ifelse(y == "vat", "visceral adipose tissue (L) per SD", "abdominal subcutaneous adipose tissue (L) per SD"), "-", nm), sample = sx[[1]], n = nrow(x),
          estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) } } }
share <- function(t, f, x) { bv <- coef(lm(as.formula(paste("vat ~", t, "+", f)), x))[[t]]; ba <- coef(lm(as.formula(paste("asat ~", t, "+", f)), x))[[t]]; bv / (bv + ba) }
for (nm in names(fits)[c(1, 2, 5)]) { t <- fits[[nm]]$t; f <- fits[[nm]]$f
  for (sx in list(list("all", c(0, 1)), list("women", 0), list("men", 1))) { x <- sv[Sex %in% sx[[2]]]; ff <- if (sx[[1]] == "all") f else sub("factor\\(Sex\\) \\+ ", "", f)
    bs_ <- replicate(500, { ii <- sample.int(nrow(x), replace = TRUE); share(t, ff, x[ii]) })
    add(part = "B", measure = paste("visceral share of the extra abdominal fat -", nm), sample = sx[[1]], n = nrow(x), estimate = share(t, ff, x),
        lo = quantile(bs_, 0.025), hi = quantile(bs_, 0.975), p = NA_real_) } }

## same visit against a decade apart, in the same participants: per cm of each discordance, and paired bootstrap
## differences (500 resamples) in the visceral-fat coefficient per cm and in the partial R2
pr2 <- function(t, f, x) { r2b <- summary(lm(as.formula(paste("vat ~", f)), x))$r.squared; r2 <- summary(lm(as.formula(paste("vat ~", t, "+", f)), x))$r.squared; (r2 - r2b) / (1 - r2b) }
bcm <- function(t, f, x) coef(lm(as.formula(paste("vat ~", t, "+", f)), x))[[t]]
for (t in c("d2", "d0")) { f <- if (t == "d2") F_SV else F_B0; co <- summary(lm(as.formula(paste("vat ~", t, "+", f)), sv))$coefficients[t, ]
  add(part = "B", measure = paste("visceral adipose tissue (L) per cm of the discordance -", ifelse(t == "d2", "same visit", "decade apart (same restricted model)")), sample = "first imaging visit",
      n = nrow(sv), estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) }
bt_ <- replicate(500, { ii <- sample.int(nrow(sv), replace = TRUE); x <- sv[ii]
  c(dpr2 = pr2("d2", F_SV, x) - pr2("d0", F_B0, x), dcm = bcm("d2", F_SV, x) - bcm("d0", F_B0, x)) })
add(part = "B", measure = "difference, same visit minus decade apart: partial R2 for visceral fat", sample = "first imaging visit", n = nrow(sv),
    estimate = pr2("d2", F_SV, sv) - pr2("d0", F_B0, sv), lo = quantile(bt_["dpr2", ], 0.025), hi = quantile(bt_["dpr2", ], 0.975), p = 2 * min(mean(bt_["dpr2", ] <= 0), mean(bt_["dpr2", ] >= 0)))
add(part = "B", measure = "difference, same visit minus decade apart: visceral adipose tissue (L) per cm", sample = "first imaging visit", n = nrow(sv),
    estimate = bcm("d2", F_SV, sv) - bcm("d0", F_B0, sv), lo = quantile(bt_["dcm", ], 0.025), hi = quantile(bt_["dcm", ], 0.975), p = 2 * min(mean(bt_["dcm", ] <= 0), mean(bt_["dcm", ] >= 0)))

## ================= C. stability across the two imaging visits =================
if (!is.null(s3) && !is.null(a3)) {
st2 <- Reduce(function(a, b) merge(a, b, by = "id"), list(i2[, .(id)], i3[, .(id)], base[, .(id, Sex)], a2, setnames(copy(a3), "age", "age3"), bs))
st2 <- st2[complete.cases(st2[, .(Sex, age, age3, wc2, bmi2, wc3, bmi3)])]
st2[, `:=`(pwc2 = predict_at(i2[match(st2$id, i2$id)], st2$age, st2$Sex), pwc3 = predict_at(i3[match(st2$id, i3$id)], st2$age3, st2$Sex))]
st2[, `:=`(d2 = resid(lm(pwc2 ~ wc2)), d3 = resid(lm(pwc3 ~ wc3)))]
st2[, `:=`(f2 = resid(lm(d2 ~ age + factor(Sex) + ns(wc2,3) + ns(bmi2,3))), f3 = resid(lm(d3 ~ age3 + factor(Sex) + ns(wc3,3) + ns(bmi3,3))))]
fz <- function(a, b) { r <- cor(a, b); se <- 1 / sqrt(length(a) - 3); c(r = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se)) }
say("stability set: %d participants with proteins, waist, BMI and age at both imaging visits; median interval %.1f years", nrow(st2), median(st2$age3 - st2$age))
for (pr in list(list("measured waist", "wc2", "wc3"), list("proteomic waist (pWC)", "pwc2", "pwc3"), list("discordance (each visit anchored on its own tape measurement)", "d2", "d3"),
                list("discordance at fixed waist, BMI, age and sex at each visit", "f2", "f3"))) {
  r <- fz(st2[[pr[[2]]]], st2[[pr[[3]]]])
  add(part = "C", measure = paste("correlation between the two imaging visits -", pr[[1]]), sample = "both imaging visits", n = nrow(st2), estimate = r[["r"]], lo = r[["lo"]], hi = r[["hi"]], p = NA_real_) }
## the stable difference between men and women inflates both correlations: within each sex, and adjusted for sex and age
st2[, `:=`(aw2 = resid(lm(wc2 ~ factor(Sex) + age)), aw3 = resid(lm(wc3 ~ factor(Sex) + age3)), ap2 = resid(lm(pwc2 ~ factor(Sex) + age)), ap3 = resid(lm(pwc3 ~ factor(Sex) + age3)),
           ab2 = resid(lm(bmi2 ~ factor(Sex) + age)), ab3 = resid(lm(bmi3 ~ factor(Sex) + age3)))]
for (pr in list(list("measured waist, adjusted for sex and age", "aw2", "aw3"), list("proteomic waist (pWC), adjusted for sex and age", "ap2", "ap3"),
                list("BMI, adjusted for sex and age", "ab2", "ab3"))) {
  r <- fz(st2[[pr[[2]]]], st2[[pr[[3]]]])
  add(part = "C", measure = paste("correlation between the two imaging visits -", pr[[1]]), sample = "both imaging visits", n = nrow(st2), estimate = r[["r"]], lo = r[["lo"]], hi = r[["hi"]], p = NA_real_) }
## how stable would the discordance be if it were nothing but tape error? If the proteins tracked the true waist and the
## tape added independent error at each visit, the discordance at fixed WC, BMI, age and sex would keep only the part
## of the true waist that the noisy tape fails to adjust for, and its correlation between the visits would be at most
## 1 - r, with r the between-visit correlation of the waist at fixed BMI, age and sex (91_tape_error_bound.R gives r
## in everyone measured at both imaging visits)
st2[, `:=`(cw2 = resid(lm(wc2 ~ ns(bmi2,3) + factor(Sex) + age)), cw3 = resid(lm(wc3 ~ ns(bmi3,3) + factor(Sex) + age3)))]
r <- fz(st2$cw2, st2$cw3)
add(part = "C", measure = "correlation between the two imaging visits - measured waist at fixed BMI, age and sex", sample = "both imaging visits", n = nrow(st2), estimate = r[["r"]], lo = r[["lo"]], hi = r[["hi"]], p = NA_real_)
add(part = "C", measure = "largest correlation of the size-adjusted discordance that tape error alone could produce (1 - r)", sample = "both imaging visits", n = nrow(st2), estimate = 1 - r[["r"]], lo = 1 - r[["hi"]], hi = 1 - r[["lo"]], p = NA_real_)
for (sx in list(list("women", 0), list("men", 1))) { x <- st2[Sex == sx[[2]]]
  for (pr in list(list("measured waist", "wc2", "wc3"), list("proteomic waist (pWC)", "pwc2", "pwc3"), list("discordance at fixed waist, BMI, age and sex at each visit", "f2", "f3"))) {
    r <- fz(x[[pr[[2]]]], x[[pr[[3]]]])
    add(part = "C", measure = paste("correlation between the two imaging visits -", pr[[1]]), sample = sx[[1]], n = nrow(x), estimate = r[["r"]], lo = r[["lo"]], hi = r[["hi"]], p = NA_real_) } }
}

## ================= D. the two protein samples on equal terms; Olink plates of the imaging-visit samples =================
## In B each discordance is adjusted for the size measured with its own proteins. Here both are related to VAT with the
## same covariates, age at both visits, sex, and WC and BMI at baseline and at the first imaging visit, so that neither
## the visit at which size was measured nor the anchoring differs between them; paired bootstrap (500 resamples of
## participants, with its own seed, after everything above, so that A-C are unchanged). Then the plates: if the two
## imaging-visit samples of a participant shared an Olink plate, plate effects could raise the correlation in C.
set.seed(2027)
F_BOTH <- "Age0 + age + factor(Sex) + ns(WC0,3) + ns(BMI0,3) + ns(wc2,3) + ns(bmi2,3)"
LB <- c(z2 = "imaging-visit proteins", z0 = "baseline proteins")
for (t in c("z2", "z0")) {
  for (y in c("vat", "asat")) { co <- summary(lm(as.formula(paste(y, "~", t, "+", F_BOTH)), sv))$coefficients[t, ]
    add(part = "D", measure = paste(ifelse(y == "vat", "visceral adipose tissue (L) per SD", "abdominal subcutaneous adipose tissue (L) per SD"), "- same covariates (age, sex, WC and BMI at both visits):", LB[[t]]),
        sample = "first imaging visit", n = nrow(sv), estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) }
  add(part = "D", measure = paste("partial R2 for visceral fat - same covariates (age, sex, WC and BMI at both visits):", LB[[t]]), sample = "first imaging visit", n = nrow(sv),
      estimate = pr2(t, F_BOTH, sv), lo = NA_real_, hi = NA_real_, p = NA_real_) }
for (t in c("d2", "d0")) { co <- summary(lm(as.formula(paste("vat ~", t, "+", F_BOTH)), sv))$coefficients[t, ]
  add(part = "D", measure = paste("visceral adipose tissue (L) per cm of the discordance - same covariates (age, sex, WC and BMI at both visits):", LB[[sub("d", "z", t)]]),
      sample = "first imaging visit", n = nrow(sv), estimate = co[[1]], lo = co[[1]] - 1.96 * co[[2]], hi = co[[1]] + 1.96 * co[[2]], p = co[[4]]) }
bt4 <- replicate(500, { ii <- sample.int(nrow(sv), replace = TRUE); x <- sv[ii]
  c(dpr2 = pr2("d2", F_BOTH, x) - pr2("d0", F_BOTH, x), dcm = bcm("d2", F_BOTH, x) - bcm("d0", F_BOTH, x)) })
for (k in list(list("dpr2", "partial R2 for visceral fat", pr2("d2", F_BOTH, sv) - pr2("d0", F_BOTH, sv)),
               list("dcm", "visceral adipose tissue (L) per cm", bcm("d2", F_BOTH, sv) - bcm("d0", F_BOTH, sv))))
  add(part = "D", measure = paste("difference, imaging-visit minus baseline proteins, same covariates (age, sex, WC and BMI at both visits):", k[[2]]), sample = "first imaging visit", n = nrow(sv),
      estimate = k[[3]], lo = quantile(bt4[k[[1]], ], 0.025), hi = quantile(bt4[k[[1]], ], 0.975), p = 2 * min(mean(bt4[k[[1]], ] <= 0), mean(bt4[k[[1]], ] >= 0)))
if (exists("st2")) {
  pv <- c("Plate used for sample run | Instance 2", "Plate used for sample run | Instance 3"); stopifnot(all(pv %in% names(sel)))
  pl <- sel[match(st2$id, sel$id), .(p2 = get(pv[1]), p3 = get(pv[2]))]
  ok <- !is.na(pl$p2) & !is.na(pl$p3) & pl$p2 != "" & pl$p3 != ""
  if (!all(ok)) say("plate not recorded for %d of the %d participants (left out of the plate comparison)", sum(!ok), nrow(st2))
  same <- ok & pl$p2 == pl$p3
  say("plates in the stability set (%d participants): first imaging visit %d plates, repeat imaging visit %d plates; both samples on one plate: %d",
      sum(ok), uniqueN(pl$p2[ok]), uniqueN(pl$p3[ok]), sum(same))
  for (k in list(list("participants whose two imaging-visit samples were measured on the same Olink plate", sum(same)),
                 list("Olink plates holding the first imaging-visit samples", uniqueN(pl$p2[ok])), list("Olink plates holding the repeat imaging-visit samples", uniqueN(pl$p3[ok]))))
    add(part = "D", measure = k[[1]], sample = "both imaging visits", n = sum(ok), estimate = k[[2]], lo = NA_real_, hi = NA_real_, p = NA_real_)
  x <- st2[ok & !same]
  x[, `:=`(g2 = resid(lm(d2 ~ age + factor(Sex) + ns(wc2,3) + ns(bmi2,3))), g3 = resid(lm(d3 ~ age3 + factor(Sex) + ns(wc3,3) + ns(bmi3,3))))]
  r <- fz(x$g2, x$g3)
  add(part = "D", measure = "correlation between the two imaging visits - discordance at fixed waist, BMI, age and sex at each visit, the two samples on different plates",
      sample = "both imaging visits", n = nrow(x), estimate = r[["r"]], lo = r[["lo"]], hi = r[["hi"]], p = NA_real_) }

tab <- rbindlist(res); fwrite(tab, file.path(TB, "T123_same_visit_mri.csv"))
say("---- B and C"); options(width = 250)
print(tab[, .(part, measure = substr(measure, 1, 120), sample, n, est = signif(estimate, 3), lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 2))], nrows = 200)
say("DONE")
