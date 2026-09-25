## 91_tape_error_bound.R
## How stable would the size-adjusted discordance be if it were nothing but error in the tape measurement?
## If the proteins tracked the true waist exactly and the tape added independent error at each visit, the discordance
## at fixed WC, BMI, age and sex would still keep the part of the true waist that the noisy tape fails to adjust for,
## and that part persists. Its correlation between two visits is then at most 1 - r, where r is the between-visit
## correlation of the waist at fixed BMI, age and sex (the reliability of the tape given BMI, times the stability of
## the true waist). This script estimates r between the first and the repeat imaging visit, in everyone with the
## measurements and in the participants of 90_same_visit_mri_and_selection.R (stability set), and compares 1 - r
## with the observed correlation of the discordance (0.695).
## Output: tables/T125_tape_error_bound.csv
suppressPackageStartupMessages({ library(data.table); library(splines) })
RAW <- "/home/data/heamei/nmrLR/ukbnmr_met_rawdata"; PR <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/protomics"
W <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
num <- function(x) suppressWarnings(as.numeric(x)); say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
bs <- fread(file.path(RAW, "6-1Body_size_measures_participant.csv"), select = c("Participant ID", "Waist circumference | Instance 2", "Body mass index (BMI) | Instance 2",
            "Waist circumference | Instance 3", "Body mass index (BMI) | Instance 3")); setnames(bs, c("id", "wc2", "bmi2", "wc3", "bmi3"))
ag <- fread(file.path(RAW, "2Recruitment.csv"), select = c("Participant ID", "Age when attended assessment centre | Instance 2", "Age when attended assessment centre | Instance 3")); setnames(ag, c("id", "age2", "age3"))
sx <- fread(file.path(RAW, "BioX_total_full_with_age_months.csv"), select = c("Participant ID", "Sex")); setnames(sx, c("id", "sex"))
ph <- fread("/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12/pro53013_新诊_newdead.csv", select = c("Participant ID", "Sex")); setnames(ph, c("id", "sex0"))
sx <- merge(sx, ph, by = "id", all = TRUE); sx[, sex := num(sex)]; sx[is.na(sex), sex := num(sex0)]; sx[, sex0 := NULL]
say("sex available for %d participants (baseline analysis file used where the first file lacks it)", sum(!is.na(sx$sex)))
d <- Reduce(function(a, b) merge(a, b, by = "id"), list(bs, ag, sx))
for (v in c("wc2", "bmi2", "wc3", "bmi3", "age2", "age3")) d[[v]] <- num(d[[v]])
d <- d[complete.cases(d[, .(wc2, bmi2, wc3, bmi3, age2, age3, sex)])]
say("sex coding: %s", paste(names(table(d$sex)), collapse = ", "))
ids23 <- intersect(fread(file.path(PR, "Instance2", "2Olink_Instance2_part1_olink_instance_2.csv"), select = 1)[[1]],
                   fread(file.path(PR, "Instance3", "2Olink_Instance3_part1_olink_instance_3.csv"), select = 1)[[1]])
res <- list()
est <- function(x, lab) {
  a2 <- resid(lm(wc2 ~ factor(sex) + age2, x)); a3 <- resid(lm(wc3 ~ factor(sex) + age3, x))
  c2 <- resid(lm(wc2 ~ ns(bmi2, 3) + factor(sex) + age2, x)); c3 <- resid(lm(wc3 ~ ns(bmi3, 3) + factor(sex) + age3, x))
  r_a <- cor(a2, a3); r_c <- cor(c2, c3); se <- 1 / sqrt(nrow(x) - 3)
  ci <- function(r) tanh(atanh(r) + c(-1.96, 1.96) * se)
  say("%s: n %d; waist between the imaging visits, adjusted for sex and age r = %.3f; at fixed BMI, sex and age r = %.3f (%.3f to %.3f); bound for pure tape error 1 - r = %.3f (%.3f to %.3f)",
      lab, nrow(x), r_a, r_c, ci(r_c)[1], ci(r_c)[2], 1 - r_c, 1 - ci(r_c)[2], 1 - ci(r_c)[1])
  res[[length(res) + 1]] <<- data.table(sample = lab, n = nrow(x), r_waist_sex_age = r_a, r_waist_fixed_bmi = r_c, lo = ci(r_c)[1], hi = ci(r_c)[2],
                                         bound_pure_error = 1 - r_c, bound_lo = 1 - ci(r_c)[2], bound_hi = 1 - ci(r_c)[1]) }
est(d, "everyone with waist, BMI and age at both imaging visits")
est(d[id %in% ids23], "participants with proteins at both imaging visits (before the missing-protein filter of 90)")
tab <- rbindlist(res); fwrite(tab, file.path(TB, "T125_tape_error_bound.csv")); print(tab)
say("observed correlation of the size-adjusted discordance between the imaging visits: 0.695 (0.661 to 0.727), T123")
say("DONE")
