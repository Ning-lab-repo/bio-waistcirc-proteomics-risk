## Figure4__risk_at_same_size.R
## Figure 4. At the same measured waist and BMI, a higher proteomic waist marks a higher risk of cardiometabolic disease.
##   A. Odds of a first recorded diagnosis (hospital or death records) within 10 years per SD of proWCdelta, in participants free of each endpoint at
##      baseline, with age, sex, deprivation, smoking and splines of measured waist and BMI held fixed (Table 2).
##   B. Standardised 10-year risk at the 10th and 90th percentiles of proWCdelta conditional on waist, BMI, age and sex,
##      with bootstrap 95% confidence intervals.
## Reads T79 (57_primary_incident_table2.R) and T112 (76_standardised_risk_ci.R). Output: Figure4.pdf and source tables.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
W   <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
OUT <- file.path(W, "figures")
COL <- c(pro = "#6A3D9A", low = "#C9B6DD")
theme_prowc <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
          axis.line = element_line(color = "black", linewidth = 0.8), axis.ticks = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(color = "black", size = 9), axis.title = element_text(color = "black", size = 10),
          plot.title = element_text(face = "bold", hjust = 0, size = 11), plot.title.position = "plot",
          legend.text = element_text(size = 9),
          plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA)) }

## ---------------- A. eight endpoints ----------------
t2 <- fread(file.path(TB, "T79_table2_primary_incident.csv"))
t2 <- t2[disease != "Ischaemic heart disease (I20-I25)"]         ## alternative code set, reported in the supplement
t2[disease == "Dyslipidemia", disease := "Dyslipidaemia"]                ## British spelling, as in the text
ord <- c("Type 2 diabetes", "Heart failure", "Liver disease", "Chronic kidney disease",
         "Ischaemic heart disease", "Dyslipidaemia", "Hypertension", "Obesity diagnosis")
pa <- t2[, .(disease, n = n_free, events = events_free, OR = OR_free, lo = lo_free, hi = hi_free)]
pa[, disease := factor(disease, rev(ord))]
pa[, lab := sprintf("%.2f (%.2f–%.2f)   %s / %s", OR, lo, hi, format(events, big.mark = ","), format(n, big.mark = ","))]
fwrite(pa, file.path(OUT, "source_table_Figure4A_odds_ratios.csv"))
pA <- ggplot(pa, aes(x = OR, y = disease)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", colour = COL[["pro"]], linewidth = 0.6) +
  geom_point(colour = COL[["pro"]], size = 2.2) +
  geom_text(aes(x = 2.15, label = lab), hjust = 0, size = 2.8, family = "Arial") +
  annotate("text", x = 2.15, y = 8.75, label = "OR (95% CI)   events / participants", hjust = 0, size = 2.8,
           family = "Arial", fontface = "bold") +
  scale_x_log10(breaks = c(1, 1.25, 1.5, 2)) + coord_cartesian(xlim = c(0.95, 2.05), ylim = c(0.6, 8.9), clip = "off") +
  labs(x = "Odds ratio per SD of proWCΔ (log scale)", y = NULL,
       title = "A  First recorded diagnosis within 10 years") +
  theme_prowc() + theme(plot.margin = margin(4, 175, 4, 4))

## ---------------- B. standardised absolute risk ----------------
sr <- fread(file.path(TB, "T112_standardised_risk_with_ci.csv"))
pb <- rbind(sr[, .(disease, pct = "10th", risk = risk_p10, lo = lo_p10, hi = hi_p10)],
            sr[, .(disease, pct = "90th", risk = risk_p90, lo = lo_p90, hi = hi_p90)])
pb[, disease := factor(disease, ord[1:4])]
pb[, pct := factor(pct, c("10th", "90th"), labels = c("10th percentile", "90th percentile"))]
rt <- sr[, .(disease = factor(disease, ord[1:4]), lab = sprintf("×%.2f\n(%.2f–%.2f)", ratio, lo_ratio, hi_ratio),
             y = 100 * hi_p90 + 0.55)]
fwrite(pb, file.path(OUT, "source_table_Figure4B_standardised_risk.csv"))
pB <- ggplot(pb, aes(x = disease, y = 100 * risk, colour = pct)) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.18, position = position_dodge(0.5), linewidth = 0.6) +
  geom_point(size = 2.4, position = position_dodge(0.5)) +
  geom_text(data = rt, aes(x = disease, y = y, label = lab), inherit.aes = FALSE, size = 2.7, family = "Arial", lineheight = 0.9) +
  scale_colour_manual(values = c(`10th percentile` = COL[["low"]], `90th percentile` = COL[["pro"]]),
                      name = "proWCΔ, conditional on waist, BMI, age and sex") +
  scale_x_discrete(labels = function(x) vapply(strsplit(x, " "), function(w) paste(paste(head(w, -1), collapse = " "), tail(w, 1), sep = "\n"), "")) +   # break before the last word
  scale_y_continuous(limits = c(0, 6.1), breaks = 0:6) +
  labs(x = NULL, y = "Standardised 10-year risk (%)", title = "B  Absolute risk at the same waist and BMI") +
  theme_prowc() + theme(legend.position = "bottom", legend.title = element_text(size = 9)) +
  guides(colour = guide_legend(title.position = "top"))

fig <- wrap_elements(full = pA) / wrap_elements(full = pB) + plot_layout(heights = c(1, 1.05))
ggsave(file.path(OUT, "Figure4.pdf"), fig, width = 180, height = 170, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure4.png"), fig, width = 180, height = 170, units = "mm", dpi = 300)
cat("Figure 4 written to", OUT, "\n")
