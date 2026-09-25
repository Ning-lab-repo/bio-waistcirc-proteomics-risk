## 84_supplementary_table6_gnhs.R
## Supplementary Table 6: the external comparison in the Guangzhou Nutrition and Health Study, formatted from the
## source tables written by G5_gnhs_continuous_waist.R, G9_direct_transport.R, G4_gnhs_discordance_risk.R and
## G7_gnhs_metabolic_profile.R (source_table_SupplementaryResults6_*). No new analysis.
## Output: source_table_SupplementaryTable6_GNHS.csv (three columns, as Supplementary Table 4; rows with the second and
## third columns empty are section headings).
suppressPackageStartupMessages(library(data.table))
SD <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
rd <- function(f) fread(file.path(SD, paste0("source_table_SupplementaryResults6_", f, ".csv")))
mi <- function(x) gsub("-", "−", x, fixed = TRUE)                    ## typographic minus
ci <- function(e, lo, hi, d = 2) mi(sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), e, lo, hi))
nev <- function(n, ev) sprintf("%s (%s)", format(n, big.mark = ","), format(ev, big.mark = ","))
rows <- list(); add <- function(a, b = "", c = "") rows[[length(rows) + 1]] <<- data.table(a, b, c)

## A. waist recorded in centimetres (485 participants of the discovery cohort)
cw <- rd("GNHS_cm_waist_MetS"); tw <- rd("GNHS_transported_weights")
add("Incident metabolic syndrome, waist recorded in centimetres (odds ratio per SD of proWCΔ)")
pick <- function(d, m) { r <- d[model == m]; stopifnot(nrow(r) == 1); r }
for (p in list(c("age, sex, measured waist and BMI", "Age, sex, measured waist and BMI (primary analysis)"),
               c("age, sex, splines of waist and BMI", "Age, sex and natural splines of waist and BMI"),
               c("age, sex, waist, BMI and the baseline component", "Age, sex, measured waist, BMI and the baseline component of the syndrome"))) {
  r <- pick(cw, p[1]); add(p[2], nev(r$n, r$events), ci(r$OR, r$lo, r$hi)) }
r <- pick(tw, "age, sex, measured waist and BMI")
add("Weights fitted in UK Biobank on the shared proteins and applied unchanged; age, sex, measured waist and BMI", nev(r$n, r$events), ci(r$OR, r$lo, r$hi))

## B. waist recorded only as the sex-specific threshold flag (discovery and validation cohorts)
rc <- rd("GNHS_risk_by_cohort")[adjustment == "age, sex, body size"]
add("Incident events, waist recorded only as the threshold flag (odds ratio per SD; age, sex, the flag and a spline of BMI held fixed)")
lab <- c(MetS = "Metabolic syndrome", TG = "Raised triglycerides", Glucose = "Raised fasting glucose", HDL = "Low HDL cholesterol",
         BP = "Raised blood pressure", Waist = "Central obesity")
for (k in names(lab)) for (co in c("discovery", "validation")) {
  r <- rc[component == k & cohort == co]
  if (nrow(r) == 1) add(sprintf("%s, %s cohort", lab[[k]], co), nev(r$n, r$events), ci(r$OR, r$lo, r$hi))
  else add(sprintf("%s, %s cohort", lab[[k]], co), "–", "too few events") }

## C. baseline metabolic measures at the same age, sex, measured waist and BMI (485 participants)
mp <- rd("GNHS_metabolic_profile")
add("Baseline metabolic measures at the same age, sex, measured waist and BMI (difference per SD of proWCΔ)")
for (q in list(c("TG", "Triglycerides, log scale", 3), c("HDL", "HDL cholesterol, mmol/L", 3), c("SBP", "Systolic blood pressure, mmHg", 2),
               c("DBP", "Diastolic blood pressure, mmHg", 2), c("Glu", "Fasting glucose, mmol/L", 3), c("LDL", "LDL cholesterol, mmol/L", 3),
               c("TC", "Total cholesterol, mmol/L", 3))) {
  mk <- q[1]; r <- mp[marker == mk]; stopifnot(nrow(r) == 1)
  add(q[2], format(r$n, big.mark = ","), ci(r$beta, r$lo, r$hi, as.integer(q[3]))) }

out <- rbindlist(rows); setnames(out, c("Analysis", "Participants (events)", "Estimate (95% CI)"))
fwrite(out, file.path(SD, "source_table_SupplementaryTable6_GNHS.csv"))
print(out, nrows = 100)
