suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(mgcv)
  library(png)
  library(grid)
})

trim_png_whitespace <- function(img, tol = 0.97, alpha_tol = 0.01, pad = 1L) {
  d <- dim(img)
  if (length(d) == 2) {
    mask <- img < tol
  } else {
    ch <- d[3]
    if (ch >= 3) {
      rgb <- img[, , 1:3, drop = FALSE]
      bright <- apply(rgb, c(1, 2), function(v) all(v > tol))
      if (ch >= 4) {
        alpha <- img[, , 4]
        mask <- (!bright) & (alpha > alpha_tol)
      } else {
        mask <- !bright
      }
    } else {
      mask <- img[, , 1] < tol
    }
  }

  idx <- which(mask, arr.ind = TRUE)
  if (nrow(idx) == 0) return(img)

  rmin <- max(min(idx[, 1]) - pad, 1)
  rmax <- min(max(idx[, 1]) + pad, nrow(mask))
  cmin <- max(min(idx[, 2]) - pad, 1)
  cmax <- min(max(idx[, 2]) + pad, ncol(mask))

  if (length(d) == 2) {
    img[rmin:rmax, cmin:cmax, drop = FALSE]
  } else {
    img[rmin:rmax, cmin:cmax, , drop = FALSE]
  }
}

read_fast_select <- function(file, select_cols, id_col = 1L) {
  id_name <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    header <- names(
      data.table::fread(
        file,
        nrows = 0,
        data.table = FALSE,
        showProgress = FALSE,
        check.names = FALSE
      )
    )
    id_name <- if (is.numeric(id_col)) header[id_col] else id_col
    keep <- intersect(unique(c(id_name, select_cols)), header)
    df <- data.table::fread(
      file,
      select = keep,
      data.table = FALSE,
      showProgress = FALSE,
      check.names = FALSE
    )
  } else {
    message("Package 'data.table' not found, falling back to read.csv (slower).")
    df <- read.csv(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    id_name <- if (is.numeric(id_col)) names(df)[id_col] else id_col
    keep <- intersect(unique(c(id_name, select_cols)), names(df))
    df <- df[, keep, drop = FALSE]
  }

  if (is.null(id_name)) {
    id_name <- if (is.numeric(id_col)) names(df)[1] else id_col
  }
  rownames(df) <- as.character(df[[id_name]])
  df[[id_name]] <- NULL
  df
}

fmt_p <- function(p) {
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "p=NA"
  ok <- !is.na(p)
  out[ok & p < 0.001] <- "p<0.001"
  out[ok & p >= 0.001] <- paste0("p=", formatC(p[ok & p >= 0.001], format = "f", digits = 3))
  out
}

p_to_star <- function(p, show_ns = FALSE) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01 & p < 0.05] <- "*"
  if (show_ns) out[!is.na(p) & out == ""] <- "ns"
  out
}

get_curve_stats <- function(df_one_sex) {
  dat <- df_one_sex %>%
    dplyr::select(Age_month, WC, proWC) %>%
    dplyr::filter(!is.na(Age_month), !is.na(WC), !is.na(proWC)) %>%
    dplyr::mutate(Diff = proWC - WC)

  # level difference (paired, robust)
  p_level <- suppressWarnings(
    wilcox.test(dat$proWC, dat$WC, paired = TRUE, exact = FALSE)$p.value
  )

  # shape difference: whether Diff varies with age
  k_use <- min(10, max(5, floor(length(unique(dat$Age_month)) / 4)))
  m_const <- gam(Diff ~ 1, data = dat, method = "REML")
  m_shape <- gam(Diff ~ s(Age_month, k = k_use, bs = "cs"), data = dat, method = "REML")

  cmp <- anova(m_const, m_shape, test = "Chisq")
  p_col <- grep("^Pr\\(", names(cmp), value = TRUE)[1]
  p_shape <- if (is.na(p_col) || !nzchar(p_col)) NA_real_ else cmp[[p_col]][2]
  if (is.na(p_shape)) {
    s_tab <- summary(m_shape)$s.table
    if (!is.null(s_tab) && nrow(s_tab) > 0) {
      p_shape <- s_tab[1, ncol(s_tab)]
    }
  }

  tibble(p_level = p_level, p_shape = p_shape)
}

pairwise_bmi_within_measure <- function(df_b) {
  all_res <- list()
  idx <- 1L
  for (sx in levels(df_b$Sex)) {
    for (ms in c("WC", "proWC")) {
      sub <- df_b %>%
        dplyr::filter(Sex == sx, !is.na(BMI_group))
      lv <- levels(sub$BMI_group)
      combs <- combn(lv, 2, simplify = FALSE)
      tmp <- lapply(combs, function(cp) {
        x <- sub[sub$BMI_group == cp[1], ms, drop = TRUE]
        y <- sub[sub$BMI_group == cp[2], ms, drop = TRUE]
        p <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
        data.frame(
          analysis = "BMI_pairwise_within_measure",
          Sex = sx,
          Measurement = ms,
          Group1 = cp[1],
          Group2 = cp[2],
          n1 = length(x),
          n2 = length(y),
          p_raw = p,
          stringsAsFactors = FALSE
        )
      })
      tb <- dplyr::bind_rows(tmp)
      tb$p_adj <- p.adjust(tb$p_raw, method = "BH")
      tb$star <- p_to_star(tb$p_adj, show_ns = TRUE)
      all_res[[idx]] <- tb
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(all_res)
}

pairwise_wc_vs_prowc_within_bmi <- function(df_b) {
  all_res <- list()
  idx <- 1L
  for (sx in levels(df_b$Sex)) {
    for (bg in levels(df_b$BMI_group)) {
      sub <- df_b %>%
        dplyr::filter(Sex == sx, BMI_group == bg) %>%
        dplyr::filter(!is.na(WC), !is.na(proWC))
      p <- suppressWarnings(wilcox.test(sub$proWC, sub$WC, paired = TRUE, exact = FALSE)$p.value)
      all_res[[idx]] <- data.frame(
        analysis = "WC_vs_proWC_within_BMI",
        Sex = sx,
        Measurement = "WC_vs_proWC",
        Group1 = bg,
        Group2 = bg,
        n1 = nrow(sub),
        n2 = nrow(sub),
        p_raw = p,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  res <- dplyr::bind_rows(all_res) %>%
    dplyr::group_by(Sex) %>%
    dplyr::mutate(p_adj = p.adjust(p_raw, method = "BH")) %>%
    dplyr::ungroup()
  res$star <- p_to_star(res$p_adj, show_ns = TRUE)
  res
}

# -----------------------------
# 1) Read input data
# -----------------------------
data_main <- read_fast_select(
  "pro53013_新诊_newdead.csv",
  select_cols = c("date_attending_assessment_centre", "BMI"),
  id_col = 1
)

data_birth <- read_fast_select(
  "1Baseline_characteristics_participant.csv",
  select_cols = c("Year of birth", "Month of birth"),
  id_col = 1
)

data_wc <- read_fast_select(
  "WC_death_time.csv",
  select_cols = c("Age", "Sex", "WC"),
  id_col = 1
)

data_lasso <- read_fast_select(
  "lasso_WC.csv",
  select_cols = c("BioX_Adjusted"),
  id_col = 1
)

# -----------------------------
# 2) Build Age(month)
# -----------------------------
birth_data <- data.frame(
  Birth_Date = paste(data_birth$`Year of birth`, data_birth$`Month of birth`, sep = "/"),
  row.names = rownames(data_birth)
)

aas_data <- data.frame(
  aas_Date = as.character(data_main$date_attending_assessment_centre),
  row.names = rownames(data_main)
)
aas_data$aas_Date <- format(as.Date(aas_data$aas_Date), "%Y/%m")

common_rows <- intersect(rownames(birth_data), rownames(aas_data))
merged_data <- cbind(
  birth_data[common_rows, , drop = FALSE],
  aas_data[common_rows, , drop = FALSE]
)

aas_year <- as.numeric(substr(merged_data$aas_Date, 1, 4))
aas_month <- as.numeric(substr(merged_data$aas_Date, 6, 7))
birth_year <- as.numeric(substr(merged_data$Birth_Date, 1, 4))
birth_month <- as.numeric(substr(merged_data$Birth_Date, 6, 7))
merged_data$Age_month <- (aas_year - birth_year) * 12 + (aas_month - birth_month)

# -----------------------------
# 3) Build plotting table
# -----------------------------
common_temp <- intersect(rownames(data_wc), rownames(data_lasso))
temp_df <- cbind(
  data_wc[common_temp, , drop = FALSE],
  data_lasso[common_temp, , drop = FALSE]
)
temp_df$BMI <- data_main[rownames(temp_df), "BMI"]
temp_df$Age_month <- merged_data[rownames(temp_df), "Age_month"]

stopifnot(all(abs(as.numeric(temp_df$Age) - round(as.numeric(temp_df$Age))) < 1e-9, na.rm = TRUE))   ## whole years, so (45, 50] is 46-50
plot_data <- temp_df %>%
  dplyr::transmute(
    ID = rownames(temp_df),
    Age = as.numeric(Age),
    Sex = factor(Sex, levels = c(1, 0), labels = c("Men", "Women")),
    WC = as.numeric(WC),
    proWC = as.numeric(BioX_Adjusted),
    Age_month = as.numeric(Age_month),
    BMI = as.numeric(BMI)
  ) %>%
  dplyr::filter(!is.na(Sex), !is.na(WC), !is.na(proWC), !is.na(Age_month))

measure_colors <- c(   ## not blue and orange, which mark men and women in panel C and in Figure 2
  "proWC" = "#8E6BBF",
  "WC" = "#A6A6A6"
)

base_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9),
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

# -----------------------------
# 4) Panel a: smooth curves + curve-comparison stats
#    Method:
#    - level p: paired Wilcoxon test (proWC vs WC)
#    - shape p: GAM on Diff=(proWC-WC), testing whether Diff changes with age
# -----------------------------
curve_stats <- dplyr::bind_rows(
  lapply(levels(plot_data$Sex), function(sx) {
    out <- get_curve_stats(dplyr::filter(plot_data, Sex == sx))
    out$Sex <- sx
    out
  })
)

curve_stats <- curve_stats %>%
  dplyr::mutate(
    Sex = factor(Sex, levels = levels(plot_data$Sex))
  )

ann_pos <- plot_data %>%
  dplyr::group_by(Sex) %>%
  dplyr::summarise(
    x = if (Sex[1] == "Men") max(Age_month, na.rm = TRUE) - 6 else min(Age_month, na.rm = TRUE) + 18,
    y = if (Sex[1] == "Men") 79.2 else 101.4,
    hj = if (Sex[1] == "Men") 1 else 0,
    vj = if (Sex[1] == "Men") 0 else 1,
    .groups = "drop"
  )

curve_stats <- dplyr::left_join(curve_stats, ann_pos, by = "Sex") %>%
  dplyr::mutate(
    label = paste0("Level: ", sub("^p", "P", fmt_p(p_level)), "\nShape: ", sub("^p", "P", fmt_p(p_shape)))
  )

p_a <- ggplot(plot_data, aes(x = Age_month)) +
  geom_smooth(
    aes(y = proWC, color = "proWC"),
    method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 0.9
  ) +
  geom_smooth(
    aes(y = WC, color = "WC"),
    method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 0.9
  ) +
  geom_text(
    data = curve_stats,
    aes(x = x, y = y, label = label, hjust = hj, vjust = vj),
    inherit.aes = FALSE,
    size = 4.8, lineheight = 0.95,
    color = "black"
  ) +
  facet_wrap(~Sex, ncol = 2) +
  scale_x_continuous(breaks = seq(480, 840, by = 120), labels = function(v) round(v / 12)) +
  coord_cartesian(ylim = c(78.5, 102), clip = "on") +
  scale_y_continuous(
    breaks = c(80, 85, 90, 95, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(
    values = measure_colors,
    breaks = c("proWC", "WC"),
    labels = c("proWC", "WC")
  ) +
  labs(x = "Age (years)", y = "Waist circumference (cm)") +
  base_theme +
  theme(
    axis.title.x = element_text(face = "bold", size = 18, family = "Arial", margin = margin(t = -18)),
    axis.title.y = element_text(face = "bold", size = 18, family = "Arial", margin = margin(r = 8)),
    axis.text = element_text(size = 13.5, family = "Arial"),
    strip.text = element_text(face = "bold", size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial")
  )

# -----------------------------
# 5) Panel b: boxplot by BMI + pairwise comparisons
# -----------------------------
plot_data_b <- plot_data %>%
  dplyr::mutate(
    BMI_group = cut(
      BMI,
      breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
      labels = c("Underweight", "NW", "OW", "OB"),
      right = TRUE
    ),
    BMI_group = factor(BMI_group, levels = c("Underweight", "NW", "OW", "OB"))
  ) %>%
  dplyr::filter(!is.na(BMI_group))

if (nrow(plot_data_b) == 0) {
  stop("plot_data_b is empty after BMI grouping. Check participant ID matching and BMI parsing.")
}

# Pairwise tests to CSV:
# 1) Within each sex and BMI group: paired WC vs proWC
pair_wc_vs_pro <- pairwise_wc_vs_prowc_within_bmi(plot_data_b)

# 2) Within each sex and measurement: pairwise BMI group comparisons
pair_bmi <- pairwise_bmi_within_measure(plot_data_b)

pairwise_results <- dplyr::bind_rows(pair_wc_vs_pro, pair_bmi)
write.csv(pairwise_results, "figure10_pairwise_results.csv", row.names = FALSE)

# Mark pairwise WC vs proWC significance in panel b (one star per BMI group per sex)
ann_b <- plot_data_b %>%
  dplyr::group_by(Sex, BMI_group) %>%
  dplyr::summarise(
    whisker_bottom = boxplot.stats(c(WC, proWC))$stats[1],
    whisker_top = boxplot.stats(c(WC, proWC))$stats[5],
    # Stars are placed just above the upper whisker.
    y_pos = whisker_top + 3.5,
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    pair_wc_vs_pro %>%
      dplyr::transmute(
        Sex,
        BMI_group = factor(Group1, levels = c("Underweight", "NW", "OW", "OB")),
        label = star
      ),
    by = c("Sex", "BMI_group")
  ) %>%
  dplyr::filter(!is.na(label), label != "", label != "ns")

plot_data_b$Sex <- factor(as.character(plot_data_b$Sex), levels = c("Men","Women"))
plot_data_long <- plot_data_b %>%
  tidyr::pivot_longer(
    cols = c(proWC, WC),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Measurement = factor(Measurement, levels = c("proWC", "WC"))
  )

b_y_min_raw <- min(ann_b$whisker_bottom, na.rm = TRUE)
if (!is.finite(b_y_min_raw)) b_y_min_raw <- min(plot_data_long$Value, na.rm = TRUE)
b_y_min <- floor((b_y_min_raw - 2) / 5) * 5
b_y_max_raw <- ceiling((max(ann_b$y_pos, na.rm = TRUE) + 2) / 5) * 5
b_y_max <- min(155, b_y_max_raw)
if (b_y_max <= max(ann_b$y_pos, na.rm = TRUE)) {
  b_y_max <- ceiling((max(ann_b$y_pos, na.rm = TRUE) + 1) / 5) * 5
}
if (!is.finite(b_y_min) || !is.finite(b_y_max) || b_y_min >= b_y_max) {
  b_y_min <- 0
  b_y_max <- 150
}
b_breaks <- pretty(c(b_y_min, b_y_max), n = 5)

plot_data_long$Sex <- factor(as.character(plot_data_long$Sex), levels = c("Men", "Women"))
ann_b$Sex <- factor(as.character(ann_b$Sex), levels = c("Men", "Women"))
p_b <- ggplot(plot_data_long, aes(x = BMI_group, y = Value, fill = Measurement)) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.45,
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8),
    width = 0.6,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_b,
    aes(x = BMI_group, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.2, fontface = "bold", color = "black"
  ) +
  facet_wrap(~Sex, ncol = 2) +
  coord_cartesian(ylim = c(b_y_min, b_y_max), clip = "on") +
  scale_y_continuous(
    breaks = b_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = measure_colors,
    breaks = c("proWC", "WC"),
    labels = c("proWC", "WC")
  ) +
  labs(x = "", y = "Waist circumference (cm)") +
  base_theme +
  theme(
    axis.title.x = element_text(face = "bold", size = 18, family = "Arial", margin = margin(t = -1)),
    axis.title.y = element_text(face = "bold", size = 18, family = "Arial", margin = margin(r = 8)),
    axis.text.y = element_text(size = 13.5, family = "Arial"),
    axis.text.x = element_text(size = 13, family = "Arial", angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold", size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial")
  )

# -----------------------------
# 6) Panel c: Age-stratified violin plot by sex for proWCΔ = proWC - WC
# -----------------------------
delta_data <- plot_data %>%
  dplyr::mutate(
    WCdelta = proWC - WC,
    Age_group = cut(
      Age,
      breaks = c(-Inf, 45, 50, 55, 60, Inf),
      labels = c("≤45", "46–50", "51–55", "56–60", ">60"),
      right = FALSE
    ),
    Age_group = factor(Age_group, levels = c("≤45", "46–50", "51–55", "56–60", ">60"))
  ) %>%
  dplyr::filter(!is.na(Age_group), !is.na(WCdelta), !is.na(Sex))

if (nrow(delta_data) == 0) {
  stop("delta_data is empty. Check Age/WC/proWC parsing and ID matching.")
}

c_lims <- as.numeric(quantile(delta_data$WCdelta, probs = c(0.02, 0.98), na.rm = TRUE))
if (any(!is.finite(c_lims)) || c_lims[2] <= c_lims[1]) {
  c_lims <- range(delta_data$WCdelta, na.rm = TRUE)
}
if (!all(is.finite(c_lims)) || c_lims[2] <= c_lims[1]) {
  c_lims <- c(-1, 1)
}
c_pad <- max(0.6, 0.12 * diff(c_lims))
c_breaks <- pretty(c(c_lims[1] - c_pad, c_lims[2] + c_pad), n = 5)
c_dodge <- 0.86

p_c <- ggplot(delta_data, aes(x = Age_group, y = WCdelta, fill = Sex)) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.9, color = "gray40") +
  geom_violin(
    position = position_dodge(width = c_dodge),
    width = 0.56,
    trim = TRUE,
    scale = "width",
    color = "gray30",
    linewidth = 0.25,
    alpha = 0.52
  ) +
  geom_boxplot(
    position = position_dodge(width = c_dodge),
    width = 0.13,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.40,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(color = Sex),
    position = position_dodge(width = c_dodge),
    shape = 21,
    size = 1.7,
    stroke = 0.35,
    fill = "white",
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = c(c_lims[1] - c_pad, c_lims[2] + c_pad), clip = "on") +
  scale_y_continuous(
    breaks = c_breaks,
    labels = function(v) sub("^-", "−", format(v, trim = TRUE)),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_fill_manual(
    values = c("Men" = "#2171B5", "Women" = "#F16913"),
    breaks = c("Men", "Women")
  ) +
  scale_color_manual(
    values = c("Men" = "#2171B5", "Women" = "#F16913"),
    breaks = c("Men", "Women")
  ) +
  labs(x = "Age strata (years)", y = "proWCΔ (cm)", fill = NULL) +
  base_theme +
  theme(
    legend.position = "top",
    axis.title.x = element_text(face = "bold", size = 18, family = "Arial", margin = margin(t = -18)),
    axis.title.y = element_text(face = "bold", size = 18, family = "Arial", margin = margin(r = 8)),
    axis.text = element_text(size = 13.5, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.key.size = unit(0.8, "lines"),
    axis.text.x = element_text(face = "bold", size = 13, family = "Arial")
  )

# -----------------------------

out <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/S8_abc"
top_row <- p_a + p_b + p_c + patchwork::plot_layout(ncol = 3, widths = c(1.12, 1.12, 0.96))
fig <- top_row + patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 20, face = "bold"))
ggsave(paste0(out, ".png"), fig, width = 16.5, height = 5.6, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
grDevices::cairo_pdf(paste0(out, ".pdf"), width = 16.5, height = 5.6, onefile = TRUE)
print(fig); grDevices::dev.off()
cat("DONE
")
