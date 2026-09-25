## G7_gnhs_metabolic_profile.R
## Two questions that table S7 sheet E can answer and the flag-based tables cannot, because it records the metabolic
## measurements themselves rather than whether they crossed a threshold.
## First: among the 485 Guangzhou participants of the same age, sex, measured waist and BMI, is a higher proteomic
## waist accompanied by a worse metabolic state at that same visit? This is the counterpart of the finding here that
## people matched on the tape measure differ in fat and in metabolic markers.
## Second: is the association of the discordance with incident metabolic syndrome explained by those measurements? The
## manuscript reports that routine clinical markers explain part of it but not all; this repeats that test in an
## independent cohort on a different platform.
## Output: G7_profile.csv, G7_adjusted.csv.
suppressPackageStartupMessages({ library(data.table); library(readxl) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

sc <- fread(file.path(TB, "G5_scores.csv"))
E <- as.data.table(read_excel(file.path(G, "tables", "tableS7.xlsx"), sheet = "E", na = c("", "NA", "NaN")))
mk <- c("SBP", "DBP", "TG", "HDL", "LDL", "TC", "Glu")
d <- merge(sc, E[, c(list(Patient_ID = pat_ID), lapply(.SD, as.numeric)), .SDcols = mk], by = "Patient_ID")
say("participants with the score and the metabolic measurements: %d", nrow(d))

## ---------------- the metabolic state at the same measured body size ----------------
## TG is modelled on the log scale, as it is skewed; the others are modelled as measured.
prof <- rbindlist(lapply(mk, function(v) {
  y <- if (v == "TG") log(d[[v]]) else d[[v]]; x <- copy(d)[, yy := y][!is.na(yy)]
  m <- lm(yy ~ z + age + factor(sex) + wc + BMI, x); co <- summary(m)$coefficients
  data.table(marker = v, scale = if (v == "TG") "log" else "measured", n = nrow(x),
             beta = co["z", 1], lo = co["z", 1] - 1.96 * co["z", 2], hi = co["z", 1] + 1.96 * co["z", 2], p = co["z", 4]) }))
prof[, txt := sprintf("%+.3f (%+.3f to %+.3f)", beta, lo, hi)]
print(prof[, .(marker, scale, n, `per SD of proWCdelta` = txt, p = signif(p, 3))])
fwrite(prof, file.path(TB, "G7_profile.csv"))

## ---------------- is the incident association explained by those measurements? ----------------
fit <- function(rhs, lab) { m <- glm(as.formula(paste("ev ~", rhs)), d, family = binomial()); co <- summary(m)$coefficients
  data.table(model = lab, n = nrow(d), events = sum(d$ev), OR = exp(co["z", 1]),
             lo = exp(co["z", 1] - 1.96 * co["z", 2]), hi = exp(co["z", 1] + 1.96 * co["z", 2]), p = co["z", 4]) }
B <- "age + factor(sex) + wc + BMI"
adj <- rbindlist(list(
  fit(paste("z +", B), "age, sex, measured waist and BMI"),
  fit(paste("z +", B, "+ log(TG) + HDL"), "+ triglycerides and HDL cholesterol"),
  fit(paste("z +", B, "+ SBP + DBP"), "+ blood pressure"),
  fit(paste("z +", B, "+ Glu"), "+ fasting glucose"),
  fit(paste("z +", B, "+ log(TG) + HDL + SBP + DBP + Glu"), "+ all of the above"),
  fit(paste("z +", B, "+ log(TG) + HDL + SBP + DBP + Glu + LDL"), "+ all of the above and LDL cholesterol")))
adj[, txt := sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)]
print(adj[, .(model, events, `OR per SD` = txt, p = signif(p, 3))])
fwrite(adj, file.path(TB, "G7_adjusted.csv"))
say("DONE")
