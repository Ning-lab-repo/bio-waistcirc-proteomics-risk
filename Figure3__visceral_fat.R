## Figure3__visceral_fat.R
## Figure 3. At the same measured waist and BMI, a higher proteomic waist marks more visceral fat.
##   A. Visceral and abdominal subcutaneous fat (MRI) across quintiles of proWCdelta, relative to the lowest quintile, at the
##      same age, sex, measured waist and BMI.
##   B. The visceral share of the extra abdominal fat that goes with proWCdelta, against the visceral share of the extra
##      fat that goes with a larger measured waist at the same BMI (what measurement error in the tape would produce),
##      overall, in participants of stable size, and in women and men.
##   C. Protein profiles: for each of the 2,920 proteins, its association with proWCdelta against its association with
##      visceral and with subcutaneous fat, all at fixed waist and BMI. A larger measured waist at the same BMI gives
##      almost the same pattern (83_protein_axis_comparator.R).
##   D. Visceral fat per SD of proWCdelta under other adjustments, with the size at the scan fixed, in participants of
##      stable size, and by sex.
## Panels B-D are read from the tables written by 75_mri_stress_test.R, 79_mri_checks.R and
## 80_protein_axis_checks.R, so the figure and the text share one source. Output: Figure3.pdf and source tables.
suppressPackageStartupMessages({ library(data.table); library(splines); library(ggplot2); library(patchwork) })
F12 <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/figure12"
BC  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_prowc/bodycomposition"
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
OUT <- file.path(W, "figures"); dir.create(OUT, showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

COL <- c(vat = "#6A3D9A", asat = "#DDA0DD", pro = "#6A3D9A", wc = "#F4B88D")
theme_prowc <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
          axis.line = element_line(color = "black", linewidth = 0.8), axis.ticks = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(color = "black", size = 9), axis.title = element_text(color = "black", size = 10),
          plot.title = element_text(face = "bold", hjust = 0, size = 11), plot.title.position = "plot", legend.text = element_text(size = 9),
          plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA)) }

## ---------------- A. fifths of proWCdelta ----------------
mri <- fread(file.path(BC, "bodycompositionInstance.2.csv"), showProgress = FALSE,
             select = c("Participant.ID", "Visceral.adipose.tissue.volume..VAT....Instance.2",
                        "Abdominal.subcutaneous.adipose.tissue.volume..ASAT....Instance.2"))
setnames(mri, c("id", "vat", "asat"))
phen <- fread(file.path(F12, "pro53013_新诊_newdead.csv"), select = c("Participant ID", "Sex", "Age", "WC", "BMI"))
setnames(phen, c("id", "Sex", "Age", "WC", "BMI"))
la <- fread(file.path(F12, "lasso_WC.csv")); setnames(la, 1, "id")
d <- Reduce(function(a, b) merge(a, b, by = "id"), list(phen, la[, .(id, dlt = BioX_Delta)], mri))
for (v in c("Age", "WC", "BMI", "vat", "asat")) d[[v]] <- num(d[[v]]); d[, Sex := as.integer(Sex)]
d <- d[complete.cases(d[, .(Age, Sex, WC, BMI, dlt, vat, asat)]) & asat > 0]
d[, q := cut(dlt, quantile(dlt, 0:5 / 5), include.lowest = TRUE, labels = 1:5)]
med <- d[, .(median_proWCdelta = median(dlt)), by = q][order(q)]
B0 <- "Age + factor(Sex) + ns(WC,3) + ns(BMI,3)"
pa <- rbindlist(lapply(c("vat", "asat"), function(y) {
  co <- summary(lm(as.formula(paste(y, "~ q +", B0)), d))$coefficients
  rbind(data.table(tissue = y, q = 1, diff = 0, lo = 0, hi = 0),
        data.table(tissue = y, q = 2:5, diff = co[paste0("q", 2:5), 1],
                   lo = co[paste0("q", 2:5), 1] - 1.96 * co[paste0("q", 2:5), 2],
                   hi = co[paste0("q", 2:5), 1] + 1.96 * co[paste0("q", 2:5), 2])) }))
pa <- merge(pa, med[, .(q = as.integer(as.character(q)), median_proWCdelta)], by = "q")
pa[, tissue_label := factor(ifelse(tissue == "vat", "Visceral", "Abdominal subcutaneous"), c("Visceral", "Abdominal subcutaneous"))]
fwrite(pa, file.path(OUT, "source_table_Figure3A_fat_by_fifth_of_proWCdelta.csv"))
xl <- sprintf("%d\n(%s)", 1:5, sub("-", "−", sprintf("%+.1f", med$median_proWCdelta), fixed = TRUE))
pA <- ggplot(pa, aes(x = q, y = diff, colour = tissue_label, group = tissue_label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_line(position = position_dodge(0.25), linewidth = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, position = position_dodge(0.25), linewidth = 0.6) +
  geom_point(position = position_dodge(0.25), size = 2) +
  scale_colour_manual(values = c(Visceral = COL[["vat"]], `Abdominal subcutaneous` = COL[["asat"]]), name = NULL) +
  scale_x_continuous(breaks = 1:5, labels = xl) +
  labs(x = "Quintile of proWCΔ (median, cm)", y = "Difference from\nlowest quintile (L)",
       title = "A  Fat across quintiles of proWCΔ") +
  theme_prowc() + theme(legend.position = c(0.36, 0.88), legend.background = element_blank())

## ---------------- B. visceral share of the extra fat ----------------
st  <- fread(file.path(TB, "T111_mri_stress_test.csv"))
rc  <- fread(file.path(TB, "T114_mri_checks.csv"))
g <- function(tab, key) { r <- tab[specification == key]; stopifnot(nrow(r) == 1); r[, .(estimate, lo, hi)] }
sh <- rbindlist(list(
  data.table(pop = "All participants", who = "Higher proWCΔ",         g(st, "visceral share, proWCdelta: all")),
  data.table(pop = "All participants", who = "Larger measured waist", g(st, "visceral share, larger measured waist: all")),
  data.table(pop = "All participants", who = "difference",            g(st, "difference in visceral share: all")),
  data.table(pop = "Stable size",      who = "Higher proWCΔ",         g(st, "visceral share, proWCdelta: stable")),
  data.table(pop = "Stable size",      who = "Larger measured waist", g(st, "visceral share, larger measured waist: stable")),
  data.table(pop = "Stable size",      who = "difference",            g(st, "difference in visceral share: stable")),
  data.table(pop = "Women",            who = "Higher proWCΔ",         g(rc, "women: with proWCdelta")),
  data.table(pop = "Women",            who = "Larger measured waist", g(rc, "women: with a larger measured WC")),
  data.table(pop = "Women",            who = "difference",            g(rc, "women: difference")),
  data.table(pop = "Men",              who = "Higher proWCΔ",         g(rc, "men: with proWCdelta")),
  data.table(pop = "Men",              who = "Larger measured waist", g(rc, "men: with a larger measured WC")),
  data.table(pop = "Men",              who = "difference",            g(rc, "men: difference"))))
fwrite(sh, file.path(OUT, "source_table_Figure3B_visceral_share.csv"))
df <- sh[who == "difference"]
sgn <- function(v) sub("-", "−", sprintf("%+.1f", v), fixed = TRUE)
df[, lab := sprintf("%s\n%s (%s to %s)", pop,sgn(100 * estimate), sub("+", "", sgn(100 * lo), fixed = TRUE), sub("+", "", sgn(100 * hi), fixed = TRUE))]
sb <- merge(sh[who != "difference"], df[, .(pop, lab)], by = "pop")
sb[, lab := factor(lab, rev(df$lab))]
pB <- ggplot(sb, aes(x = 100 * estimate, y = lab, colour = who)) +
  geom_errorbar(aes(xmin = 100 * lo, xmax = 100 * hi), width = 0.18, orientation = "y", position = position_dodge(0.5), linewidth = 0.6) +
  geom_point(size = 2.2, position = position_dodge(0.5)) +
  scale_colour_manual(values = c(`Higher proWCΔ` = COL[["pro"]], `Larger measured waist` = COL[["wc"]]), name = NULL,
                      breaks = c("Higher proWCΔ", "Larger measured waist")) +
  coord_cartesian(xlim = c(40, 67)) +
  labs(x = "Visceral share of the extra abdominal fat (%)", y = NULL, title = "B  Composition of the extra abdominal fat") +
  theme_prowc() + theme(legend.position = "bottom", legend.key.height = unit(9, "pt"), axis.text.y = element_text(size = 8),
                        plot.margin = margin(2, 14, 2, 2)) +
  guides(colour = guide_legend(nrow = 1))

## ---------------- C. the protein axis ----------------
pv <- fread(file.path(TB, "T115_protein_vat_associations.csv"))
ax <- fread(file.path(TB, "T116_protein_axis_checks.csv"))
rv <- ax[spec == "r across proteins: visceral fat vs proWCdelta"]$est
rs <- ax[spec == "r across proteins: subcutaneous fat vs proWCdelta"]$est
pc <- rbind(pv[, .(protein, depot = "vat", x = beta_delta, y = beta_vat, p = p_vat)],
            pv[, .(protein, depot = "asat", x = beta_delta, y = beta_asat, p = p_asat)])
pc[, depot_lab := factor(ifelse(depot == "vat", sprintf("Visceral fat (r = %.2f)", rv), sprintf("Abdominal subcutaneous fat (r = %.2f)", rs)),
                         c(sprintf("Visceral fat (r = %.2f)", rv), sprintf("Abdominal subcutaneous fat (r = %.2f)", rs)))]
fwrite(pv[, .(protein, name, beta_delta, p_delta, beta_vat, se_vat, p_vat, fdr_vat, beta_asat, p_asat, beta_vat_given_asat, p_vat_given_asat, coef)],
       file.path(OUT, "source_table_Figure3C_protein_axis.csv"))
top <- pv[order(p_vat)][1:8]$protein
lbl <- pc[protein %in% top & depot == "vat"]
pC <- ggplot(pc, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(colour = depot), size = 0.7, alpha = 0.45, show.legend = FALSE) +
  scale_colour_manual(values = c(vat = COL[["vat"]], asat = "#B57EDC")) +
  facet_wrap(~ depot_lab, nrow = 1) +
  labs(x = "Association with proWCΔ (SD of protein per SD of proWCΔ)", y = "Association with the fat depot\n(L per SD of protein)",
       title = "C  Protein profiles of proWCΔ and of each fat depot") +
  theme_prowc() + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9.5))
minus_lab <- function(v) sub("^-", "−", format(v, trim = TRUE))                   ## true minus signs on the axes
pC <- pC + scale_x_continuous(labels = minus_lab) + scale_y_continuous(labels = minus_lab)
pC <- if (have_repel) pC + ggrepel::geom_label_repel(data = lbl, aes(label = protein), size = 2.7, family = "Arial", min.segment.length = 0,
                                                    fill = "white", label.size = NA, label.padding = unit(0.06, "lines"), label.r = unit(0, "lines"),
                                                    segment.size = 0.25, max.overlaps = 50, force = 2, seed = 1) else
  pC + geom_text(data = lbl, aes(label = protein), size = 2.7, family = "Arial", nudge_y = 0.03, check_overlap = TRUE)

## ---------------- D. robustness ----------------
keep <- list(
  list(st, "all, baseline size fixed",                                   "Primary",               "Age, sex, splines of waist and BMI"),
  list(st, "+ Townsend index and smoking",                               "Further adjustment",    "+ Townsend index and smoking"),
  list(st, "+ ethnic background",                                        "Further adjustment",    "  and ethnic background"),
  list(st, "+ alcohol intake",                                           "Further adjustment",    "  and alcohol intake"),
  list(rc, "VAT, + triglycerides, HDL, HbA1c, CRP, ALT, GGT",            "Further adjustment",    "+ Routine blood measurements"),
  list(st, "excluding diabetes at baseline or diagnosed later",          "Exclusion, weighting",  "Excluding diabetes, baseline or later"),
  list(rc, "VAT, weighted for attendance at imaging (robust SE)",        "Exclusion, weighting",  "Weighted for attendance at imaging"),
  list(st, "splines with 5 degrees of freedom",                          "Body-size adjustment",  "Splines with 5 df"),
  list(st, "waist x BMI interaction",                                    "Body-size adjustment",  "Waist × BMI interaction"),
  list(st, "height added",                                               "Body-size adjustment",  "Height added"),
  list(st, "within sex, sex-specific splines",                           "Body-size adjustment",  "Sex-specific splines"),
  list(rc, "WC and BMI of both visits fixed, all participants",          "Time to the scan",      "Size at the scan also fixed"),
  list(st, "stable, baseline size fixed",                                "Time to the scan",      "Size stable between visits"),
  list(st, "stable, size at both visits fixed",                          "Time to the scan",      "Stable, size at both visits fixed"),
  list(st, "women",                                                      "Sex",                   "Women"),
  list(st, "men",                                                        "Sex",                   "Men"))
pdd <- rbindlist(lapply(keep, function(k) { r <- k[[1]][specification == k[[2]]]; stopifnot(nrow(r) == 1)
  data.table(group = k[[3]], label = k[[4]], n = r$n, estimate = r$estimate, lo = r$lo, hi = r$hi) }))
pdd[, group := factor(group, c("Primary", "Further adjustment", "Exclusion, weighting", "Body-size adjustment", "Time to the scan", "Sex"))]
pdd[, label := factor(label, rev(label))]
fwrite(pdd, file.path(OUT, "source_table_Figure3D_robustness.csv"))
prim <- pdd[group == "Primary"]$estimate
pD <- ggplot(pdd, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = prim, linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", colour = COL[["vat"]], linewidth = 0.6) +
  geom_point(colour = COL[["vat"]], size = 2) +
  geom_text(aes(x = 0.86, label = sprintf("%.2f (%.2f–%.2f)", estimate, lo, hi)), hjust = 0, size = 2.7, family = "Arial") +
  geom_text(aes(x = 1.50, label = format(n, big.mark = ",")), hjust = 1, size = 2.7, family = "Arial") +
  ## column headers above the first row group
  geom_text(data = data.table(group = factor("Primary", levels(pdd$group)), x = c(0.86, 1.50), h = c(0, 1),
                              txt = c("Estimate (95% CI)", "n")),
            aes(x = x, y = Inf, label = txt, hjust = h), inherit.aes = FALSE, vjust = -0.6, size = 2.7, fontface = "bold", family = "Arial") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  coord_cartesian(xlim = c(0, 0.85), clip = "off") +
  labs(x = "Visceral fat per SD of proWCΔ (L)", y = NULL, title = "D  Visceral fat per SD of proWCΔ") +
  theme_prowc() +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 8.5),
        panel.spacing.y = unit(2, "pt"), plot.margin = margin(2, 130, 2, 2))

## rows are laid out independently: aligning them would force narrow panels onto the plotting area of the wide labels
row1 <- wrap_elements(full = (pA | pB) + plot_layout(widths = c(1.1, 1)))
fig <- row1 / wrap_elements(full = pC) / wrap_elements(full = pD) + plot_layout(heights = c(1.1, 0.95, 1.55))
ggsave(file.path(OUT, "Figure3.pdf"), fig, width = 180, height = 240, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure3.png"), fig, width = 180, height = 240, units = "mm", dpi = 300)
cat("Figure 3 written to", OUT, "\n")
