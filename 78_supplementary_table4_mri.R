## 78_supplementary_table4_mri.R
## Supplementary Table 4: every abdominal MRI estimate reported in the paper, in one table, formatted from the result
## tables of 74_abdominal_mri.R (T110), 75_mri_stress_test.R (T111), 79_mri_checks.R (T114) and
## 80_protein_axis_checks.R (T116), so that the table and the text share one source; the rows on proteins measured at the
## imaging visits and on selection into the proteomics project come from 90_same_visit_mri_and_selection.R (T123, T124).
## Supplementary Table 5: the 30 proteins most strongly associated with visceral fat at fixed WC and BMI (T115). Protein names from
##   86_protein_names_all.R (T119).
## Output: source_table_SupplementaryTable4_abdominal_MRI.csv (section rows have empty second and third columns),
##         source_table_SupplementaryTable5_proteins_visceral_fat.csv
suppressPackageStartupMessages(library(data.table))
W  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
SD <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
t10 <- fread(file.path(TB, "T110_abdominal_mri.csv")); t11 <- fread(file.path(TB, "T111_mri_stress_test.csv"))
t14 <- fread(file.path(TB, "T114_mri_checks.csv")); t16 <- fread(file.path(TB, "T116_protein_axis_checks.csv"))
mn <- function(s) gsub("-", "−", s, fixed = TRUE)
f  <- function(b, lo, hi, k = 3) if (is.na(lo)) mn(sprintf(paste0("%.", k, "f"), b)) else mn(sprintf(paste0("%.", k, "f (%.", k, "f to %.", k, "f)"), b, lo, hi))
one <- function(tab, keep) { r <- tab[which(keep)]; stopifnot(nrow(r) == 1); r }
rows <- list(); sec <- function(s) rows[[length(rows) + 1]] <<- data.table(s, "", "")
add <- function(lab, r, k = 3, scale = 1) rows[[length(rows) + 1]] <<- data.table(lab, format(r$n, big.mark = ","),
  f(scale * r$b, scale * r$lo, scale * r$hi, k))
e11 <- function(tst, spc) { k <- t11[["test"]] == tst & t11[["specification"]] == spc; one(t11, k)[, .(n, b = estimate, lo, hi)] }
e10 <- function(out, pop) { k <- grepl(out, t10[["outcome"]], fixed = TRUE) & grepl(pop, t10[["population"]], fixed = TRUE); one(t10, k)[, .(n, b = beta, lo, hi)] }
e14 <- function(spc) { k <- t14[["specification"]] == spc; one(t14, k)[, .(n, b = estimate, lo, hi)] }
e16 <- function(spc) { k <- t16[["spec"]] == spc; one(t16, k)[, .(n, b = est, lo, hi)] }

sec("Visceral adipose tissue, L per SD of proWCΔ")
add("Primary: age, sex and splines of baseline WC and BMI",       e11("2. timing", "all, baseline size fixed"))
add("Primary, robust standard errors",                            e14("VAT, primary model with robust SE"))
add("+ Townsend index and smoking",                               e11("3b. confounding and sex", "+ Townsend index and smoking"))
add("+ Townsend index, smoking and ethnic background",            e11("3b. confounding and sex", "+ ethnic background"))
add("+ Townsend index, smoking and alcohol intake",               e11("3b. confounding and sex", "+ alcohol intake"))
add("Primary, participants with routine blood measurements",      e14("VAT, primary model (same participants)"))
add("+ triglycerides, HDL cholesterol, HbA1c, CRP, ALT and GGT",   e14("VAT, + triglycerides, HDL, HbA1c, CRP, ALT, GGT"))
add("Excluding diabetes at baseline or diagnosed later",          e11("3b. confounding and sex", "excluding diabetes at baseline or diagnosed later"))
add("Weighted for attendance at imaging (robust standard errors)", e14("VAT, weighted for attendance at imaging (robust SE)"))
add("Splines with five degrees of freedom",                       e11("3. adjustment form", "splines with 5 degrees of freedom"))
add("Interaction between the WC and BMI splines",                 e11("3. adjustment form", "waist x BMI interaction"))
add("Spline of standing height added",                            e11("3. adjustment form", "height added"))
add("Sex-specific splines of WC and BMI",                         e11("3. adjustment form", "within sex, sex-specific splines"))
add("Time from blood sample to scan added",                       e14("time from blood sample to scan added"))
add("Baseline WC and BMI fixed, participants with imaging-visit WC and BMI", e14("baseline WC and BMI fixed (same participants)"))
add("WC and BMI of both visits fixed",                            e14("WC and BMI of both visits fixed, all participants"))
add("WC and BMI of both visits fixed, women",                     e16("women: WC and BMI of both visits fixed"))
add("WC and BMI of both visits fixed, men",                       e16("men: WC and BMI of both visits fixed"))
add("Stable size between the visits, baseline WC and BMI fixed",  e11("2. timing", "stable, baseline size fixed"))
add("Stable size, WC and BMI of both visits fixed",               e11("2. timing", "stable, size at both visits fixed"))
add("Women",                                                      e11("3b. confounding and sex", "women"))
add("Men",                                                        e11("3b. confounding and sex", "men"))
sec("Abdominal subcutaneous adipose tissue, L per SD of proWCΔ")
add("Primary",                                                    e10("abdominal subcutaneous", "all, baseline"))
add("Stable size between the visits",                             e10("abdominal subcutaneous", "stable"))
sec("Ratio of visceral to subcutaneous volume per SD of proWCΔ")
add("Primary",                                                    e10("visceral to subcutaneous ratio", "all, baseline"))
add("Stable size between the visits",                             e10("visceral to subcutaneous ratio", "stable"))
add("Stable size, log ratio",                                     e11("2. timing", "stable, log visceral:subcutaneous"))
sec("Partial R² of proWCΔ at fixed WC, BMI, age and sex")
add("Visceral fat",                                               e14("vat: partial R2 of proWCdelta at fixed WC, BMI, age and sex"))
add("Subcutaneous fat",                                           e14("asat: partial R2 of proWCdelta at fixed WC, BMI, age and sex"))
add("Visceral fat, women",                                        e16("women: visceral fat"))
add("Subcutaneous fat, women",                                    e16("women: subcutaneous fat"))
add("Visceral fat, men",                                          e16("men: visceral fat"))
add("Subcutaneous fat, men",                                      e16("men: subcutaneous fat"))
add("Visceral fat, beyond routine blood measurements",            e14("VAT: partial R2 of proWCdelta beyond routine blood measurements"))
sec("Visceral share of the extra abdominal fat, %")
add("With proWCΔ, at fixed WC and BMI",                      e11("1. tape measurement error", "visceral share, proWCdelta: all"), 1, 100)
add("With a larger measured WC, at fixed BMI",                    e11("1. tape measurement error", "visceral share, larger measured waist: all"), 1, 100)
add("Difference, percentage points",                              e11("1. tape measurement error", "difference in visceral share: all"), 1, 100)
add("Linear WC term in both arms: difference",                    e14("all, linear WC term in both arms: difference"), 1, 100)
add("Stable size: with proWCΔ",                              e11("1. tape measurement error", "visceral share, proWCdelta: stable"), 1, 100)
add("Stable size: with a larger measured WC",                     e11("1. tape measurement error", "visceral share, larger measured waist: stable"), 1, 100)
add("Stable size: difference, percentage points",                 e11("1. tape measurement error", "difference in visceral share: stable"), 1, 100)
add("Women: with proWCΔ",                                    e14("women: with proWCdelta"), 1, 100)
add("Women: with a larger measured WC",                           e14("women: with a larger measured WC"), 1, 100)
add("Women: difference, percentage points",                       e14("women: difference"), 1, 100)
add("Men: with proWCΔ",                                      e14("men: with proWCdelta"), 1, 100)
add("Men: with a larger measured WC",                             e14("men: with a larger measured WC"), 1, 100)
add("Men: difference, percentage points",                         e14("men: difference"), 1, 100)
sec("Selection into imaging")
add("Odds ratio for attending imaging per SD of proWCΔ, at the same age, sex, WC and BMI",
    e11("4. selection", "odds of attending imaging per SD of proWCdelta"), 2)
## proteins measured at the imaging visits, and selection into the proteomics project (90_same_visit_mri_and_selection.R)
t23 <- fread(file.path(TB, "T123_same_visit_mri.csv")); t24 <- fread(file.path(TB, "T124_selection_sensitivity.csv"))
e23 <- function(msr, smp) { k <- t23[["measure"]] == msr & t23[["sample"]] == smp; one(t23, k)[, .(n, b = estimate, lo, hi)] }
e24 <- function(sub, out) { k <- t24[["subset"]] == sub & t24[["outcome"]] == out; one(t24, k)[, .(n, b = estimate, lo, hi)] }
V1 <- "first imaging visit"; SV <- "same visit: imaging-visit proteins, size and age at that visit"; DA <- "decade apart: baseline proteins (same restricted model), baseline size"
sec("Proteins measured at the first imaging visit (model restricted to that panel), per SD of its discordance")
add("Visceral fat, L; WC, BMI and age at that visit fixed",        e23(paste("visceral adipose tissue (L) per SD -", SV), V1))
add("Subcutaneous fat, L; WC, BMI and age at that visit fixed",    e23(paste("abdominal subcutaneous adipose tissue (L) per SD -", SV), V1))
add("Visceral fat, L; same model with baseline proteins, baseline size fixed",     e23(paste("visceral adipose tissue (L) per SD -", DA), V1))
add("Subcutaneous fat, L; same model with baseline proteins, baseline size fixed", e23(paste("abdominal subcutaneous adipose tissue (L) per SD -", DA), V1))
add("Visceral share of the extra fat, %; imaging-visit proteins",  e23(paste("visceral share of the extra abdominal fat -", SV), "all"), 1, 100)
add("Visceral share, %; larger measured WC at that visit",         e23("visceral share of the extra abdominal fat - comparator: larger measured waist at the imaging visit, BMI fixed", "all"), 1, 100)
## the two protein samples on equal terms: the same covariates for both (section D of 90)
EQ <- "same covariates (age, sex, WC and BMI at both visits):"
sec("Imaging-visit against baseline proteins, with age, sex and the WC and BMI of both visits fixed for both")
add("Visceral fat, L per SD; imaging-visit proteins",              e23(paste("visceral adipose tissue (L) per SD -", EQ, "imaging-visit proteins"), V1))
add("Visceral fat, L per SD; baseline proteins",                   e23(paste("visceral adipose tissue (L) per SD -", EQ, "baseline proteins"), V1))
add("Visceral fat, L per cm; imaging-visit proteins",              e23(paste("visceral adipose tissue (L) per cm of the discordance -", EQ, "imaging-visit proteins"), V1))
add("Visceral fat, L per cm; baseline proteins",                   e23(paste("visceral adipose tissue (L) per cm of the discordance -", EQ, "baseline proteins"), V1))
add("Difference per cm, imaging visit minus baseline",             e23(paste("difference, imaging-visit minus baseline proteins,", EQ, "visceral adipose tissue (L) per cm"), V1))
add("Partial R² for visceral fat; imaging-visit proteins",         e23(paste("partial R2 for visceral fat -", EQ, "imaging-visit proteins"), V1))
add("Partial R² for visceral fat; baseline proteins",              e23(paste("partial R2 for visceral fat -", EQ, "baseline proteins"), V1))
add("Difference in partial R², imaging visit minus baseline",      e23(paste("difference, imaging-visit minus baseline proteins,", EQ, "partial R2 for visceral fat"), V1))
sec("Stability between the two imaging visits (median 3.0 years apart), correlation")
add("Discordance at fixed WC, BMI, age and sex",                   e23("correlation between the two imaging visits - discordance at fixed waist, BMI, age and sex at each visit", "both imaging visits"))
add("The same, both samples on different Olink plates",            e23("correlation between the two imaging visits - discordance at fixed waist, BMI, age and sex at each visit, the two samples on different plates", "both imaging visits"))
add("Largest correlation that tape error alone could produce",     e23("largest correlation of the size-adjusted discordance that tape error alone could produce (1 - r)", "both imaging visits"))
add("Proteomic waist, adjusted for sex and age",                   e23("correlation between the two imaging visits - proteomic waist (pWC), adjusted for sex and age", "both imaging visits"))
add("Measured waist, adjusted for sex and age",                    e23("correlation between the two imaging visits - measured waist, adjusted for sex and age", "both imaging visits"))
sec("Selection into the proteomics project: visceral fat, L per SD of proWCΔ")
add("Without the participants selected by consortium members",     e24("without consortium-selected", "visceral adipose tissue (L)"))
add("Without the repeat-imaging participants",                     e24("without repeat-imaging participants", "visceral adipose tissue (L)"))
add("Without both",                                                e24("without both", "visceral adipose tissue (L)"))
tab <- rbindlist(rows); setnames(tab, c("Specification", "Participants", "Estimate (95% CI)"))
print(tab)
fwrite(tab, file.path(SD, "source_table_SupplementaryTable4_abdominal_MRI.csv"))

## ---------------- Supplementary Table 5 ----------------
pv <- fread(file.path(TB, "T115_protein_vat_associations.csv"))
un <- fread(file.path(TB, "T119_protein_names_all.csv"))   ## 86_protein_names_all.R: UniProt names, curated for the proteins named in the paper
pv <- merge(pv, un[, .(protein = gene_symbol, short = protein_name)], by = "protein", all.x = TRUE); pv[!is.na(short) & short != "", name := short]
setorder(pv, p_vat); t5 <- pv[1:30]
stopifnot(all(!is.na(t5$name) & t5$name != ""))
fmt <- function(b, se) mn(sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se))
t5 <- t5[, .(Protein = protein, Name = name, `Visceral fat, L per SD (95% CI)` = fmt(beta_vat, se_vat),
             `P, visceral fat` = mn(sprintf("%.1e", p_vat)), `Subcutaneous fat, L per SD` = mn(sprintf("%.3f", beta_asat)),
             `Association with proWCΔ, SD per SD` = mn(sprintf("%.3f", beta_delta)),
             `proWC model coefficient` = mn(ifelse(coef == 0, "0", sprintf("%.3f", coef))))]
print(t5)
fwrite(t5, file.path(SD, "source_table_SupplementaryTable5_proteins_visceral_fat.csv"))
cat("DONE\n")
