suppressPackageStartupMessages({library(glmnet); library(ggplot2); library(patchwork); library(data.table)})
setwd("/path/to/UKB_data")   # directory holding the LASSO result objects
L <- readRDS("lasso_results1.rds"); R <- readRDS("ridge_results1.rds")
lw <- fread("lasso_WC.csv"); setnames(lw, names(lw)[1], "id")
## Supplementary Figure 1. The model deposited as Supplementary Table 1 is L$final_model (lambda = 0.0146); its
## coefficients and Supplementary Table 2 are written by 02_export_lasso_coefficients.R, which refits it. Panels A-C
## show the regularisation path stored with the model object, evaluated in the full analysis set (in sample).

## headline out-of-fold metrics
r2_l_p  <- L$r2; r2_r_p  <- R$r2
r2_l_a  <- cor(lw$Actual_WC, lw$BioX_Adjusted)^2
r2_r_a  <- cor(R$biox_metrics$biox_adjusted, lw$Actual_WC)^2
nf_l    <- sum(as.matrix(L$final_model$beta)[,1][!(rownames(L$final_model$beta) %in% c("Age","Sex"))] != 0); nf_r <- sum(as.matrix(R$final_model$beta)[,1][!(rownames(R$final_model$beta) %in% c("Age","Sex"))] != 0)
lam_l   <- L$final_model$lambda; lam_r <- R$final_model$lambda
cat(sprintf("LASSO feats=%d lambda=%.5f R2pWC=%.4f R2proWC=%.4f\n", nf_l, lam_l, r2_l_p, r2_l_a))
cat(sprintf("RIDGE feats=%d lambda=%.5f R2pWC=%.4f R2proWC=%.4f\n", nf_r, lam_r, r2_r_p, r2_r_a))

th <- theme_classic(base_size = 8.5) +
  theme(legend.position = "top", legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 7), legend.key.height = unit(0.5,"lines"), legend.key.width = unit(1.1, "lines"),
        axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 8), plot.tag = element_text(size = 11, face = "bold"))
COL <- c(LASSO = "#1F4E9C", Ridge = "#C0272D")
mlab <- function(v) sub("^-", "−", format(v, trim = TRUE, big.mark = ","))   ## true minus signs on the axes

d1 <- rbind(data.table(x = log10(L$lambda_path), y = L$n_features, m = "LASSO"),
            data.table(x = log10(R$lambda_path), y = R$n_features, m = "Ridge"))
pA <- ggplot(d1, aes(x, y, colour = m)) + geom_line(linewidth = 0.7) +
  geom_vline(xintercept = log10(lam_l), colour = COL["LASSO"], linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = log10(lam_r), colour = COL["Ridge"],  linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = COL, labels = c(sprintf("LASSO: %s proteins", formatC(nf_l, big.mark = ",", format = "d")), sprintf("Ridge: %s proteins", formatC(nf_r, big.mark = ",", format = "d")))) +
  guides(colour = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = mlab) + scale_y_continuous(labels = mlab) +
  labs(x = expression(log[10](lambda)), y = "Number of proteins retained") + th

d2 <- rbind(data.table(x = log10(L$lambda_path), y = L$r2_mbmi_lambda, m = "LASSO", k = "proWC"),
            data.table(x = log10(L$lambda_path), y = L$r2_pbmi_lambda, m = "LASSO", k = "pWC"),
            data.table(x = log10(R$lambda_path), y = R$r2_mbmi_lambda, m = "Ridge", k = "proWC"),
            data.table(x = log10(R$lambda_path), y = R$r2_pbmi_lambda, m = "Ridge", k = "pWC"))
pB <- ggplot(d2, aes(x, y, colour = m, linetype = k)) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = COL, guide = guide_legend(order = 1, nrow = 2)) +
  scale_linetype_manual(values = c(proWC = "solid", pWC = "dashed"), guide = guide_legend(order = 2, nrow = 2)) +
  scale_x_continuous(labels = mlab) +
  labs(x = expression(log[10](lambda)), y = expression(R^2~"against measured WC")) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.04, vjust = -0.35, size = 2.2, lineheight = 0.95,
           label = sprintf("Out of fold:\nLASSO pWC %.3f, proWC %.3f\nRidge pWC %.3f, proWC %.3f",
                           r2_l_p, r2_l_a, r2_r_p, r2_r_a)) + th

d3 <- rbind(data.table(x = log10(L$lambda_path), y = L$rmse_lambda, m = "LASSO"),
            data.table(x = log10(R$lambda_path), y = R$rmse_lambda, m = "Ridge"))
pC <- ggplot(d3, aes(x, y, colour = m)) + geom_line(linewidth = 0.7) +
  geom_vline(xintercept = log10(lam_l), colour = COL["LASSO"], linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = log10(lam_r), colour = COL["Ridge"],  linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = COL, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = mlab) +
  labs(x = expression(log[10](lambda)), y = "RMSE (cm)") + th

## D, E: coefficients ranked; the proteins of largest magnitude are marked and listed with their coefficients, in
## the two empty corners (positive top left, negative bottom right), instead of being labelled with leader lines
beta_plot <- function(fm, col, thr, lab, out) {
  b <- as.matrix(fm$beta)[,1]; b <- b[!(names(b) %in% c("Age","Sex"))]
  dt <- data.table(protein = names(b), beta = as.numeric(b))
  dt[protein == "HLA.DRA", protein := "HLA-DRA"]; dt[protein == "ERVV.1", protein := "ERVV-1"]   ## symbols as in Supplementary Table 1
  dt[, big := abs(beta) > thr]
  dt[, yy := rank(beta, ties.method = "first")]
  fwrite(dt[order(yy), .(protein, coefficient = beta, rank = yy, labelled = big)], out)
  fmt <- function(v) sub("^-", "−", sprintf("%.2f", v))
  pos <- dt[big == TRUE & beta > 0][order(-beta)]; neg <- dt[big == TRUE & beta < 0][order(beta)]
  lst <- function(d, head) paste(c(head, sprintf("%s  %s", d$protein, fmt(d$beta))), collapse = "\n")
  ggplot(dt, aes(beta, yy)) +
    geom_point(data = dt[big == FALSE], colour = "grey70", size = 0.5) +
    geom_point(data = dt[big == TRUE],  colour = col, size = 1.4) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
    annotate("text", x = min(dt$beta), y = nrow(dt), hjust = 0, vjust = 1, size = 2.2, lineheight = 0.95,
             label = lst(pos, "Largest positive")) +
    annotate("text", x = max(dt$beta), y = 1, hjust = 1, vjust = 0, size = 2.2, lineheight = 0.95,
             label = lst(neg, "Largest negative")) +
    scale_x_continuous(labels = mlab) +
    labs(x = bquote(beta~"("*.(lab)*" model)"), y = "Proteins (ranked)") +
    th + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
pD <- beta_plot(R$final_model, COL[["Ridge"]], 1.40, "ridge", "/tmp/prowc_s1/source_table_SupplementaryFigure1D_ridge_coefficients.csv")
pE <- beta_plot(L$final_model, COL[["LASSO"]], 1.26, "LASSO", "/tmp/prowc_s1/source_table_SupplementaryFigure1E_lasso_coefficients.csv")

## F: accuracy of reduced scores (nested cross-validation)
tk <- fread("/tmp/prowc_s1/topk_performance.csv")
tp <- tk[k > 0]
## reference: waist score from age, sex and routine biochemistry (HbA1c, HDL, triglycerides, CRP), out of fold, in the
## 43,533 participants with these measures (36_biochemical_score_head_to_head.R, table T33)
bio <- fread("/tmp/prowc_s1/T33_biochemical_score_summary.csv")
r2_bio <- bio[grepl("biochemical score", quantity), value]
pF <- ggplot(tp, aes(k, R2)) +
  geom_hline(yintercept = tk[k == 0, R2], linetype = "dashed", linewidth = 0.5, colour = "grey45") +
  annotate("text", x = 5, y = tk[k == 0, R2], label = sprintf("age and sex only: %.3f", tk[k == 0, R2]),
           hjust = 0, vjust = -0.6, size = 2.3, colour = "grey30") +
  geom_hline(yintercept = r2_bio, linetype = "dotted", linewidth = 0.6, colour = "grey30") +
  annotate("text", x = 5, y = r2_bio, label = sprintf("age, sex and routine biochemistry (HbA1c, HDL, triglycerides, CRP): %.3f", r2_bio),
           hjust = 0, vjust = -0.6, size = 2.3, colour = "grey30") +
  geom_line(colour = COL[["LASSO"]], linewidth = 0.7) + geom_point(colour = COL[["LASSO"]], size = 1.8) +
  geom_text(aes(label = sprintf("%.3f", R2)), vjust = -0.9, size = 2.3) +
  scale_x_log10(breaks = tp$k, labels = formatC(tp$k, big.mark = ",", format = "d")) + scale_y_continuous(limits = c(0.15, 0.86)) +
  labs(x = "Number of proteins in the score (log scale)", y = expression("Out-of-fold "*R^2*" for measured WC")) + th
## plotted values, written as source tables
fwrite(d1[, .(penalty = m, log10_lambda = x, proteins_retained = y)], "/tmp/prowc_s1/source_table_SupplementaryFigure1A_proteins_retained.csv")
fwrite(d2[, .(penalty = m, prediction = k, log10_lambda = x, R2_in_sample = y)], "/tmp/prowc_s1/source_table_SupplementaryFigure1B_R2_path.csv")
fwrite(d3[, .(penalty = m, log10_lambda = x, RMSE_cm_in_sample = y)], "/tmp/prowc_s1/source_table_SupplementaryFigure1C_RMSE_path.csv")
fwrite(data.table(penalty = c("LASSO", "Ridge"), lambda = c(lam_l, lam_r), log10_lambda = log10(c(lam_l, lam_r)),
                  proteins_retained = c(nf_l, nf_r), R2_pWC_out_of_fold = c(r2_l_p, r2_r_p), R2_proWC_out_of_fold = c(r2_l_a, r2_r_a)),
       "/tmp/prowc_s1/source_table_SupplementaryFigure1_dashed_lines.csv")
fwrite(rbind(tk[, .(score = ifelse(k == 0, "age and sex only", paste0("first ", k, " proteins")), k, R2)],
             data.table(score = "age, sex and routine biochemistry", k = NA_integer_, R2 = r2_bio)),
       "/tmp/prowc_s1/source_table_SupplementaryFigure1F_reduced_scores.csv")
## three path panels side by side, the two coefficient panels, then the reduced scores; about 180 mm square
fig <- (pA | pB | pC) / (pD | pE) / pF + plot_layout(heights = c(1, 1.45, 0.85)) +
  plot_annotation(tag_levels = "A")
ggsave("/tmp/prowc_s1/S1_model_development.pdf", fig, width = 7.1, height = 7.1, units = "in", device = cairo_pdf)
ggsave("/tmp/prowc_s1/S1_model_development.png", fig, width = 7.1, height = 7.1, units = "in", dpi = 300, bg = "white")
cat("DONE\n")
