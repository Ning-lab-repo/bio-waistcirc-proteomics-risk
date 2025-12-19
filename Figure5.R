# 剂量反应曲线
library(survival)
library(ggplot2)
library(ggExtra)

# 读取数据
complete_data_WC <- read.csv("WC_death_time.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results <- read.csv("lasso_WC.csv", header = T, row.names = 1)
temp_df_WC <- cbind(complete_data_WC, results)

# 处理数据
df <- na.omit(temp_df_WC)

# 修改函数：添加Age和Sex作为协变量
generate_plot_data_pspline_adjusted <- function(exposure_var, outcome_var, time_var, data, df = 4, reference = "sample") {
  # 获取年龄和性别变量名
  age_var <- ifelse("Age" %in% names(data), "Age", "age")
  sex_var <- ifelse("Sex" %in% names(data), "Sex", "sex")
  
  # 构建模型公式 - 添加Age和Sex作为协变量
  formula <- as.formula(paste("Surv(", time_var, ", ", outcome_var, ") ~ pspline(", exposure_var, ", df = ", df, ") + ", age_var, " + ", sex_var))
  
  # 使用惩罚样条拟合Cox模型
  cox_model <- coxph(
    formula = formula,
    data = data,
    x = TRUE, 
    y = TRUE
  )
  
  # 创建预测范围 - 使用原始暴露变量的范围
  exposure_vals <- seq(min(data[[exposure_var]], na.rm = TRUE), 
                       max(data[[exposure_var]], na.rm = TRUE), 
                       length.out = 52879)
  
  # 创建预测数据框 - 固定Age和Sex
  # Age取平均值，Sex固定为男性（假设1为男性）
  age_mean <- mean(data[[age_var]], na.rm = TRUE)
  sex_male <- 1  # 固定为男性
  
  pred_df <- data.frame(
    exposure_vals = exposure_vals,
    age_value = age_mean,
    sex_value = sex_male
  )
  names(pred_df)[1] <- exposure_var
  names(pred_df)[2] <- age_var
  names(pred_df)[3] <- sex_var
  
  # 进行预测
  pred <- predict(cox_model, 
                  newdata = pred_df, 
                  type = "lp", 
                  se.fit = TRUE,
                  reference = reference)
  
  # 计算风险比和置信区间
  hr <- exp(pred$fit)
  lower_ci <- exp(pred$fit - 1.96 * pred$se.fit)
  upper_ci <- exp(pred$fit + 1.96 * pred$se.fit)
  
  # 返回整理好的数据
  return(data.frame(
    Exposure = exposure_vals,
    HR = hr,
    Lower = lower_ci,
    Upper = upper_ci
  ))
}

# 修改画图函数 - 保留网格线，使用Arial字体
create_hr_plot_with_margins <- function(plot_df, data_df, exposure_var, title, x_lab = NULL) {
  # 如果没有提供x_lab，使用暴露变量名
  if (is.null(x_lab)) {
    x_lab <- exposure_var
  }
  
  # 创建一个专门用于ggMarginal的数据框，包含x和y变量
  marginal_data <- data.frame(
    x_value = data_df[[exposure_var]],
    y_value = plot_df$HR
  )
  
  # 创建基础散点图（使用透明点）- 保留网格线
  p <- ggplot(marginal_data, aes(x = x_value, y = y_value)) +
    geom_point(alpha = 0) +  # 透明点
    labs(
      title = title,
      x = x_lab,
      y = "Hazard Ratio for mortality"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),  # 设置Arial字体
      # 保留网格线（默认的theme_minimal有网格线）
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(face = "bold",size = 7),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
  
  # 现在添加剂量反应曲线和置信区间
  p <- p + 
    geom_line(data = plot_df, aes(x = Exposure, y = HR), linewidth = 1.2, inherit.aes = FALSE) +
    geom_ribbon(data = plot_df, aes(x = Exposure, ymin = Lower, ymax = Upper), 
                alpha = 0.2, fill = "blue", inherit.aes = FALSE) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red")
  
  # 创建带边际的图形
  p_with_margins <- ggMarginal(
    p,
    type = "histogram",
    margins = "both",
    size = 4,
    fill = "lightblue",
    color = "black",
    bins = 20
  )
  
  return(p_with_margins)
}


# 生成调整Age和Sex后的数据（Sex固定为男性）

# WC的数据
wc_ac_data <- generate_plot_data_pspline_adjusted("WC", "all_cause_death", "survival_time", df, df = 4, reference = "sample")
wc_cvd_data <- generate_plot_data_pspline_adjusted("WC", "cvd_death", "survival_time", df, df = 4, reference = "sample")
wc_noncvd_data <- generate_plot_data_pspline_adjusted("WC", "non_cvd_death", "survival_time", df, df = 4, reference = "sample")
wc_cancer_data <- generate_plot_data_pspline_adjusted("WC", "cancer_death", "survival_time", df, df = 4, reference = "sample")
wc_diabetes_data <- generate_plot_data_pspline_adjusted("WC", "diabetes_death", "survival_time", df, df = 4, reference = "sample")

# proWC的数据
prowc_ac_data <- generate_plot_data_pspline_adjusted("BioX_Adjusted", "all_cause_death", "survival_time", df, df = 4, reference = "sample")
prowc_cvd_data <- generate_plot_data_pspline_adjusted("BioX_Adjusted", "cvd_death", "survival_time", df, df = 4, reference = "sample")
prowc_noncvd_data <- generate_plot_data_pspline_adjusted("BioX_Adjusted", "non_cvd_death", "survival_time", df, df = 4, reference = "sample")
prowc_cancer_data <- generate_plot_data_pspline_adjusted("BioX_Adjusted", "cancer_death", "survival_time", df, df = 4, reference = "sample")
prowc_diabetes_data <- generate_plot_data_pspline_adjusted("BioX_Adjusted", "diabetes_death", "survival_time", df, df = 4, reference = "sample")

# proWCΔ的数据
prowcΔ_ac_data <- generate_plot_data_pspline_adjusted("BioX_Delta", "all_cause_death", "survival_time", df, df = 4, reference = "sample")
prowcΔ_cvd_data <- generate_plot_data_pspline_adjusted("BioX_Delta", "cvd_death", "survival_time", df, df = 4, reference = "sample")
prowcΔ_noncvd_data <- generate_plot_data_pspline_adjusted("BioX_Delta", "non_cvd_death", "survival_time", df, df = 4, reference = "sample")
prowcΔ_cancer_data <- generate_plot_data_pspline_adjusted("BioX_Delta", "cancer_death", "survival_time", df, df = 4, reference = "sample")
prowcΔ_diabetes_data <- generate_plot_data_pspline_adjusted("BioX_Delta", "diabetes_death", "survival_time", df, df = 4, reference = "sample")

# 创建图形
p_wc_ac <- create_hr_plot_with_margins(wc_ac_data, df, "WC", "All-cause Mortality", "WC")
p_wc_cvd <- create_hr_plot_with_margins(wc_cvd_data, df, "WC", "CVD Mortality", "WC")
p_wc_noncvd <- create_hr_plot_with_margins(wc_noncvd_data, df, "WC", "Non-CVD Mortality", "WC")
p_wc_cancer <- create_hr_plot_with_margins(wc_cancer_data, df, "WC", "Cancer Mortality", "WC")
p_wc_diabetes <- create_hr_plot_with_margins(wc_diabetes_data, df, "WC", "Diabetes Mortality", "WC")

pro_wc_ac <- create_hr_plot_with_margins(prowc_ac_data, df, "BioX_Adjusted", "All-cause Mortality", "proWC")
pro_wc_cvd <- create_hr_plot_with_margins(prowc_cvd_data, df, "BioX_Adjusted", "CVD Mortality", "proWC")
pro_wc_noncvd <- create_hr_plot_with_margins(prowc_noncvd_data, df, "BioX_Adjusted", "Non-CVD Mortality", "proWC")
pro_wc_cancer <- create_hr_plot_with_margins(prowc_cancer_data, df, "BioX_Adjusted", "Cancer Mortality", "proWC")
pro_wc_diabetes <- create_hr_plot_with_margins(prowc_diabetes_data, df, "BioX_Adjusted", "Diabetes Mortality", "proWC")

proΔ_wc_ac <- create_hr_plot_with_margins(prowcΔ_ac_data, df, "BioX_Delta", "All-cause Mortality", "proWCΔ")
proΔ_wc_cvd <- create_hr_plot_with_margins(prowcΔ_cvd_data, df, "BioX_Delta", "CVD Mortality", "proWCΔ")
proΔ_wc_noncvd <- create_hr_plot_with_margins(prowcΔ_noncvd_data, df, "BioX_Delta", "Non-CVD Mortality", "proWCΔ")
proΔ_wc_cancer <- create_hr_plot_with_margins(prowcΔ_cancer_data, df, "BioX_Delta", "Cancer Mortality", "proWCΔ")
proΔ_wc_diabetes <- create_hr_plot_with_margins(prowcΔ_diabetes_data, df, "BioX_Delta", "Diabetes Mortality", "proWCΔ")





# 加载必要的包
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# 死亡
or_result <- read.csv("multi_death_or_results.csv", header = TRUE, fileEncoding = "GBK")

or_result$Disease_Code <- ifelse(or_result$Disease_Code == "All-Cause Death", 
                                 "All-cause Death", 
                                 or_result$Disease_Code)

# 分别处理每种指标的数据（包含p值）
process_metric_data <- function(data, metric) {
  metric_data <- data %>%
    select(Disease_Code, quintile, 
           estimate = paste0("Odds_ratio_", metric), 
           conf.low = paste0("conf.low_", metric), 
           conf.high = paste0("conf.high_", metric),
           p.adj = paste0("p.adj_", metric)) %>%  # 添加调整后的p值
    mutate(metric = case_when(
      metric == "proBMIΔ" ~ "proBMIΔ",
      metric == "BMI" ~ "BMI",
      metric == "proBMI" ~ "proBMI",
      metric == "proWCΔ" ~ "proWCΔ",
      metric == "WC" ~ "WC",
      metric == "proWC" ~ "proWC",
      TRUE ~ metric
    ))
  return(metric_data)
}

# 只处理指定的三个指标
metrics <- c("WC", "proWC", "proWCΔ")
all_metric_data <- map_dfr(metrics, ~ process_metric_data(or_result, .x))

# 设置死亡结局的因子顺序（移除Non-CVD Death）
death_levels <- c("All-cause Death", "CVD Death", "Cancer Death", "Diabetes Death")
all_metric_data$Disease_Code <- factor(all_metric_data$Disease_Code, levels = death_levels)

# 设置指标的因子顺序
metric_levels <- c("WC", "proWC", "proWCΔ")
all_metric_data$metric <- factor(all_metric_data$metric, levels = metric_levels)

# 设置五分位数的因子顺序并分配颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",  # 蓝色
  "Q2" = "#ffcc00",  # 黄色
  "Q3" = "#2ca02c",  # 绿色
  "Q4" = "#e377c2",  # 粉色
  "Q5" = "#d62728"   # 红色
)

all_metric_data$quintile <- factor(all_metric_data$quintile, 
                                   levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# 添加形状变量：参考组Q1保持圆形，显著的其他组为三角形
all_metric_data <- all_metric_data %>%
  mutate(
    # 判断是否为显著 (通常p.adj < 0.05)
    is_significant = ifelse(!is.na(p.adj) & p.adj < 0.05, TRUE, FALSE),
    # 确定点的形状：参考组Q1永远为圆形，其他显著组为三角形
    point_shape = case_when(
      quintile == "Q1" ~ 21,      # 圆形，参考组
      is_significant == TRUE ~ 24, # 三角形，显著
      TRUE ~ 21                    # 圆形，不显著
    ),
    # 确定点的大小（保持与之前一致）
    point_size = 2
  )

# 首先对数据进行log2转换
all_metric_data_log2 <- all_metric_data %>%
  mutate(
    estimate_log2 = log2(estimate),
    conf.low_log2 = log2(conf.low),
    conf.high_log2 = log2(conf.high)
  )

# 过滤掉缺失的死亡结局和Non-CVD Death
all_metric_data_log2 <- all_metric_data_log2 %>%
  filter(Disease_Code %in% death_levels) %>%
  filter(metric %in% metric_levels) %>%
  drop_na(estimate_log2)

# 使用按五分位数分配颜色的绘图代码
OR_death <- ggplot(all_metric_data_log2, aes(
  x = estimate_log2,  # 使用log2转换后的OR值
  y = quintile
)) +
  geom_vline(
    xintercept = 0,  # log2(1) = 0
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low_log2,  # 使用log2转换后的置信区间下限
      xmax = conf.high_log2, # 使用log2转换后的置信区间上限
      color = quintile  # 按五分位数分配颜色
    ),
    height = 0.5,
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    aes(
      fill = quintile,  # 按五分位数分配填充色
      color = quintile, # 按五分位数分配边框色
      shape = factor(point_shape)  # 使用形状变量
    ),
    size = all_metric_data_log2$point_size,  # 使用预定义的点大小
    stroke = 1,
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(values = quintile_colors) +  # 使用五分位数颜色
  scale_fill_manual(values = quintile_colors) +   # 使用五分位数颜色
  scale_shape_manual(
    values = c(
      "21" = 21,  # 圆形
      "24" = 24   # 三角形
    ),
    guide = "none"  # 隐藏形状图例
  ) +
  scale_x_continuous(
    labels = function(x) format(round(x, 1), nsmall = 1)  # 统一保留1位小数
  ) +
  facet_grid(metric ~ Disease_Code, scales = "free_x") +
  labs(
    x = expression(bold(log[bold("2")] * "(Odds Ratio)")),  # 在expression中添加bold()
    y = "Quintile",
    title = "",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 11),  # 统一：坐标轴标题字号
    axis.text = element_text(color = "black", size = 9),  # 统一：坐标轴刻度字号
    axis.text.y = element_text(size = 9),                 # 统一：Y轴文字大小
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 11),  # 统一：分面标签大小
    strip.text.x = element_text(size = 11, face = "bold"), # 统一：X分面标签
    strip.text.y = element_text(size = 11, face = "bold"), # 统一：Y分面标签
    legend.position = "none",
    plot.margin = margin(-18, 0, 1, 0),
    axis.title.y = element_text(margin = margin(r = 15)), 
    text = element_text(family = "Arial"),  # 统一字体家族
    plot.title = element_text(hjust = 0.5, size = 12, family = "Arial")  # 统一标题字号
  )

OR_death




