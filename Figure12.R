setwd("D:/UKB_data")
data <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("WC_death_time.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df <- cbind(complete_data_WC, results)

# 直接匹配行名添加BMI列
temp_df$BMI <- data[rownames(temp_df), "BMI"]

# 准备数据
male_data <- subset(temp_df, Sex == "1")
female_data <- subset(temp_df, Sex == "0")

# 设置datadist
dd_male <- datadist(male_data)
dd_female <- datadist(female_data)

# 创建生存对象函数
create_survival_objects <- function(data) {
  data$all_cause_surv <- with(data, Surv(survival_time, all_cause_death))
  data$cvd_surv <- with(data, Surv(survival_time, cvd_death))
  data$cancer_surv <- with(data, Surv(survival_time, cancer_death))
  data$diabetes_surv <- with(data, Surv(survival_time, diabetes_death))
  return(data)
}

male_data <- create_survival_objects(male_data)
female_data <- create_survival_objects(female_data)

# 性别分组RCS分析函数 - 不限制X轴范围
rcs_analysis_sex_stratified <- function(surv_obj, outcome_name, male_data, female_data, biomarker_type = "WC", use_log2 = FALSE) {
  
  # 定义性别分组名称和颜色
  sex_groups <- c("Male", "Female")
  sex_colors <- c(
    "Male" = "#1F77B4",    # 蓝色
    "Female" = "#FF7F0E"   # 橙色
  )
  
  # 存储模型和预测结果
  models <- list()
  predictions <- list()
  p_values <- list()
  
  # 为每个性别分组拟合模型
  for (sex_group in sex_groups) {
    data <- if(sex_group == "Male") male_data else female_data
    
    # 设置datadist
    options(datadist = if(sex_group == "Male") "dd_male" else "dd_female")
    
    # 构建公式
    if (biomarker_type == "WC") {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(WC, 4) + Age + BMI"))
    } else {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(BioX_Adjusted, 4) + Age + BMI"))
    }
    
    # 拟合模型
    models[[sex_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
    
    # 计算非线性P值
    anova_result <- anova(models[[sex_group]])
    p_values[[sex_group]] <- ifelse(biomarker_type == "WC", 
                                    anova_result["WC", "P"], 
                                    anova_result["BioX_Adjusted", "P"])
  }
  
  # 生成预测数据 - 使用每个性别组的实际范围
  pred_combined <- data.frame()
  
  for (sex_group in names(models)) {
    data <- if(sex_group == "Male") male_data else female_data
    options(datadist = if(sex_group == "Male") "dd_male" else "dd_female")
    
    # 使用该性别组的实际数据范围
    if (biomarker_type == "WC") {
      actual_min <- min(data$WC,  na.rm = TRUE)
      actual_max <- max(data$WC, na.rm = TRUE)
      biomarker_range <- seq(actual_min, actual_max, length = 100)
      
      pred <- Predict(models[[sex_group]], 
                      WC = biomarker_range, 
                      Age = median(data$Age, na.rm = TRUE),
                      BMI = median(data$BMI, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$WC, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        Sex_Group = sex_group
      )
    } else {
      # 对于proWC，使用实际数据范围
      actual_min <- min(data$BioX_Adjusted, na.rm = TRUE)
      actual_max <- max(data$BioX_Adjusted, na.rm = TRUE)
      biomarker_range <- seq(actual_min, actual_max, length = 100)
      
      pred <- Predict(models[[sex_group]], 
                      BioX_Adjusted = biomarker_range, 
                      Age = median(data$Age, na.rm = TRUE),
                      BMI = median(data$BMI, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$BioX_Adjusted, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        Sex_Group = sex_group
      )
    }
    
    pred_combined <- rbind(pred_combined, pred_df)
  }
  
  # 如果使用log2转换
  if (use_log2) {
    pred_combined$HR <- log2(pred_combined$HR)
    pred_combined$lower <- log2(pred_combined$lower)
    pred_combined$upper <- log2(pred_combined$upper)
  }
  
  # 计算合适的P值标注位置
  y_max <- max(pred_combined$upper, na.rm = TRUE)
  y_min <- min(pred_combined$lower, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 计算整体的X轴范围（覆盖所有性别组）
  if (biomarker_type == "WC") {
    x_label <- "WC"
    # 计算整体的数据范围（不限制）
    wc_min_all <- min(pred_combined$Biomarker, na.rm = TRUE)
    wc_max_all <- max(pred_combined$Biomarker, na.rm = TRUE)
  } else {
    x_label <- "proWC"
    # 计算整体的数据范围
    prowc_min_all <- min(pred_combined$Biomarker, na.rm = TRUE)
    prowc_max_all <- max(pred_combined$Biomarker, na.rm = TRUE)
  }
  
  # 绘制图表 - 不设置X轴限制
  p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = Sex_Group, fill = Sex_Group)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_hline(yintercept = ifelse(use_log2, 0, 1), 
               linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    labs(
      title = outcome_name,
      x = x_label,
      y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio"),
      color = "Sex Group",
      fill = "Sex Group"
    ) +
    #scale_x_continuous(limits = c(40, 160)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),  # 全局字体
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Arial"),
      axis.title = element_text(size = 11, face = "bold", family = "Arial"),  # 改为11
      axis.text = element_text(size = 9, family = "Arial"),  # 改为9
      axis.line = element_line(size = 0.8, color = "black"),
      axis.ticks = element_line(size = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold",family = "Arial"),  # 改为8
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  # 在图上添加P值
  y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
  i <- 1
  for (sex_group in names(p_values)) {
    p_value <- p_values[[sex_group]]
    # 动态计算X位置（使用当前图表的数据范围）
    plot_data <- p$data[p$data$Sex_Group == sex_group, ]
    x_min <- min(plot_data$Biomarker, na.rm = TRUE)
    x_max <- max(plot_data$Biomarker, na.rm = TRUE)
    x_pos <- x_min + 0.05 * (x_max - x_min)
    
    p <- p + annotate("text", 
                      x = x_pos, 
                      y = y_max - y_positions[i] * y_range,
                      label = paste("P-nonlinear ", 
                                    ifelse(p_value < 0.001, "≤ 0.001",
                                           sprintf("%.3f", p_value))),
                      size = 2.8, fontface = "bold", hjust = 0,
                      color = sex_colors[sex_group],
                      family = "Arial")  # 添加字体
    i <- i + 1
  }
  
  return(list(plot = p, 
              predictions = pred_combined,
              p_values = p_values,
              biomarker_range = if(biomarker_type == "WC") {
                list(WC_min = wc_min_all, WC_max = wc_max_all)
              } else {
                list(proWC_min = prowc_min_all, proWC_max = prowc_max_all)
              }))
}

# 运行所有性别分组分析
all_cause_wc_sex <- rcs_analysis_sex_stratified(all_cause_surv, "All-cause mortality", male_data, female_data, "WC", use_log2 = TRUE)
cvd_wc_sex <- rcs_analysis_sex_stratified(cvd_surv, "CVD mortality", male_data, female_data, "WC", use_log2 = TRUE)
cancer_wc_sex <- rcs_analysis_sex_stratified(cancer_surv, "Cancer mortality", male_data, female_data, "WC", use_log2 = TRUE)
diabetes_wc_sex <- rcs_analysis_sex_stratified(diabetes_surv, "Diabetes mortality", male_data, female_data, "WC", use_log2 = TRUE)

all_cause_prowc_sex <- rcs_analysis_sex_stratified(all_cause_surv, "All-cause mortality", male_data, female_data, "proWC", use_log2 = TRUE)
cvd_prowc_sex <- rcs_analysis_sex_stratified(cvd_surv, "CVD mortality", male_data, female_data, "proWC", use_log2 = TRUE)
cancer_prowc_sex <- rcs_analysis_sex_stratified(cancer_surv, "Cancer mortality", male_data, female_data, "proWC", use_log2 = TRUE)
diabetes_prowc_sex <- rcs_analysis_sex_stratified(diabetes_surv, "Diabetes mortality", male_data, female_data, "proWC", use_log2 = TRUE)


# 使用patchwork组合图表
library(patchwork)

# 移除图例
all_cause_prowc_sex_nolegend <- all_cause_prowc_sex$plot + theme(legend.position = "none")
cvd_prowc_sex_nolegend <- cvd_prowc_sex$plot + theme(legend.position = "none")
cancer_prowc_sex_nolegend <- cancer_prowc_sex$plot + theme(legend.position = "none")
diabetes_prowc_sex_nolegend <- diabetes_prowc_sex$plot + theme(legend.position = "none")

all_cause_wc_sex_nolegend <- all_cause_wc_sex$plot + theme(plot.title = element_blank(), legend.position = "none")
cvd_wc_sex_nolegend <- cvd_wc_sex$plot + theme(plot.title = element_blank(), legend.position = "none")
cancer_wc_sex_nolegend <- cancer_wc_sex$plot + theme(plot.title = element_blank(), legend.position = "none")
diabetes_wc_sex_nolegend <- diabetes_wc_sex$plot + theme(plot.title = element_blank(), legend.position = "none")

# 组合图表：proWC在上面，WC在下面
combined_plot_sex <- (
  (all_cause_prowc_sex_nolegend / all_cause_wc_sex_nolegend) |
    (cvd_prowc_sex_nolegend / cvd_wc_sex_nolegend) |
    (cancer_prowc_sex_nolegend / cancer_wc_sex_nolegend) |
    (diabetes_prowc_sex_nolegend / diabetes_wc_sex_nolegend)
) + 
  plot_layout(ncol = 4) &
  theme(
    legend.position = "none",
    # 确保组合后的图表也有透明背景
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# 添加共享图例
legend_sex <- cowplot::get_legend(all_cause_prowc_sex$plot + theme(legend.position = "top", legend.background = element_rect(fill = "transparent", color = NA)))

# 最终组合
final_plot_sex <- combined_plot_sex / legend_sex + 
  plot_layout(heights = c(0.9, 0.1))&
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot_sex


# 顶部
library(patchwork)

# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend_sex, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot_sex <- top_panel / combined_plot_sex + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot_sex



# Age
temp_df$Sex <- as.factor(temp_df$Sex)
age_median <- median(temp_df$Age)
temp_df$age_group <- ifelse(temp_df$Age >= age_median, "Higher Age", "Lower Age")

higher_age_df <- subset(temp_df, age_group == "Higher Age")
lower_age_df <- subset(temp_df, age_group == "Lower Age")

# 设置datadist
dd_higher_age <- datadist(higher_age_df)
dd_lower_age <- datadist(lower_age_df)

# 创建生存对象函数
create_survival_objects <- function(data) {
  data$all_cause_surv <- with(data, Surv(survival_time, all_cause_death))
  data$cvd_surv <- with(data, Surv(survival_time, cvd_death))
  data$cancer_surv <- with(data, Surv(survival_time, cancer_death))
  data$diabetes_surv <- with(data, Surv(survival_time, diabetes_death))
  return(data)
}

higher_age_df <- create_survival_objects(higher_age_df)
lower_age_df <- create_survival_objects(lower_age_df)

# 年龄分组RCS分析函数 - 按实际范围画图，不限制X轴
rcs_analysis_age_stratified <- function(surv_obj, outcome_name, age_data_list, biomarker_type = "WC", use_log2 = TRUE) {
  
  # 定义年龄分组名称和颜色
  age_groups <- c("Higher Age", "Lower Age")
  age_colors <- c(
    "Higher Age" = "#F28FB1",    # 粉色
    "Lower Age" = "#9FB4C7"      # 蓝色
  )
  
  # 存储模型和预测结果
  models <- list()
  predictions <- list()
  p_values <- list()
  actual_ranges <- list()
  
  # 为每个年龄分组拟合模型
  for (age_group in age_groups) {
    data <- age_data_list[[age_group]]
    
    # 检查是否有足够的数据
    if (nrow(data) < 10) {
      warning(paste("Age group", age_group, "has insufficient data (<10 observations) for", outcome_name))
      next
    }
    
    # 设置datadist
    options(datadist = get(paste0("dd_", tolower(gsub(" ", "_", age_group)))))
    
    # 构建公式
    if (biomarker_type == "WC") {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(WC, 4) + Sex + BMI"))
    } else {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(BioX_Adjusted, 4) + Sex + BMI"))
    }
    
    # 拟合模型
    models[[age_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
    
    # 计算非线性P值
    anova_result <- anova(models[[age_group]])
    p_values[[age_group]] <- ifelse(biomarker_type == "WC", 
                                    anova_result["WC", "P"], 
                                    anova_result["BioX_Adjusted", "P"])
    
    # 记录实际数据范围（最小值到最大值）
    if (biomarker_type == "WC") {
      actual_ranges[[age_group]] <- c(min(data$WC, na.rm = TRUE), max(data$WC, na.rm = TRUE))
    } else {
      actual_ranges[[age_group]] <- c(min(data$BioX_Adjusted, na.rm = TRUE), max(data$BioX_Adjusted, na.rm = TRUE))
    }
  }
  
  # 如果没有模型，返回空结果
  if (length(models) == 0) {
    warning(paste("No valid models for", outcome_name))
    return(NULL)
  }
  
  # 生成预测数据 - 使用每个年龄分组的实际范围（最小值到最大值）
  pred_combined <- data.frame()
  
  for (age_group in names(models)) {
    data <- age_data_list[[age_group]]
    options(datadist = get(paste0("dd_", tolower(gsub(" ", "_", age_group)))))
    
    # 使用该年龄分组的实际数据范围（最小值到最大值）
    if (biomarker_type == "WC") {
      actual_range <- actual_ranges[[age_group]]
      biomarker_range <- seq(actual_range[1], actual_range[2], length = 100)
      
      pred <- Predict(models[[age_group]], 
                      WC = biomarker_range, 
                      Sex = levels(data$Sex)[1],
                      BMI = median(data$BMI, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$WC, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        Age_Group = age_group
      )
    } else {
      actual_range <- actual_ranges[[age_group]]
      biomarker_range <- seq(actual_range[1], actual_range[2], length = 100)
      
      pred <- Predict(models[[age_group]], 
                      BioX_Adjusted = biomarker_range, 
                      Sex = levels(data$Sex)[1],
                      BMI = median(data$BMI, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$BioX_Adjusted, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        Age_Group = age_group
      )
    }
    
    pred_combined <- rbind(pred_combined, pred_df)
  }
  
  # 如果使用log2转换
  if (use_log2) {
    pred_combined$HR <- log2(pred_combined$HR)
    pred_combined$lower <- log2(pred_combined$lower)
    pred_combined$upper <- log2(pred_combined$upper)
  }
  
  # 计算整体的X轴范围（覆盖所有年龄分组）
  if (biomarker_type == "WC") {
    x_label <- "WC"
    # 计算整体的数据范围（使用实际预测范围）
    x_min <- min(pred_combined$Biomarker, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker, na.rm = TRUE)
  } else {
    x_label <- "proWC"
    # 计算整体的数据范围
    x_min <- min(pred_combined$Biomarker, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker, na.rm = TRUE)
  }
  
  # 计算合适的P值标注位置
  y_max <- max(pred_combined$upper, na.rm = TRUE)
  y_min <- min(pred_combined$lower, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 绘制图表 - 不设置X轴限制
  p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = Age_Group, fill = Age_Group)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_hline(yintercept = ifelse(use_log2, 0, 1), 
               linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = age_colors) +
    scale_fill_manual(values = age_colors) +
    labs(
      title = outcome_name,
      x = x_label,
      y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio (95% CI)"),
      color = "Age Group",
      fill = "Age Group"
    ) +
    #scale_x_continuous(limits = c(40, 160)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),  # 全局字体
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Arial"),
      axis.title = element_text(size = 11, face = "bold", family = "Arial"),  # 改为11
      axis.text = element_text(size = 9, family = "Arial"),  # 改为9
      axis.line = element_line(size = 0.8, color = "black"),
      axis.ticks = element_line(size = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold",family = "Arial"),  # 改为8
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  
  # 在图上添加P值 - 动态计算位置
  if (length(p_values) > 0) {
    y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
    i <- 1
    
    for (age_group in names(p_values)) {
      p_value <- p_values[[age_group]]
      
      if (!is.null(p_value) && !is.na(p_value)) {
        # 动态计算X位置（使用当前分组的数据范围）
        plot_data <- pred_combined[pred_combined$Age_Group == age_group, ]
        if (nrow(plot_data) > 0) {
          x_min_group <- min(plot_data$Biomarker, na.rm = TRUE)
          x_max_group <- max(plot_data$Biomarker, na.rm = TRUE)
          x_pos <- x_min_group + 0.05 * (x_max_group - x_min_group)
          
          # P值标注逻辑
          if (p_value < 0.001) {
            p_label <- "P-nonlinear ≤ 0.001"
          } else if (p_value < 0.01) {
            p_label <- paste0("P-nonlinear = ", sprintf("%.3f", p_value))
          } else if (p_value < 0.05) {
            p_label <- paste0("P-nonlinear = ", sprintf("%.3f", p_value))
          } else {
            p_label <- paste0("P-nonlinear = ", sprintf("%.3f", p_value))
          }
          
          # 添加标注
          p <- p + annotate("text", 
                            x = x_pos, 
                            y = y_max - y_positions[i] * y_range,
                            label = p_label,
                            size = 2.8, fontface = "bold", hjust = 0,
                            color = age_colors[age_group],
                            family = "Arial")
          
          i <- i + 1
        }
      }
    }
  }
  
  # 打印实际数据范围信息
  cat(paste("\n", outcome_name, "实际数据范围（", ifelse(biomarker_type == "WC", "WC", "proWC"), "）:\n"))
  for (age_group in names(actual_ranges)) {
    if (age_group %in% names(models)) {
      range <- actual_ranges[[age_group]]
      cat(paste(age_group, ":", sprintf("%.1f - %.1f", range[1], range[2]), 
                ifelse(biomarker_type == "WC", "cm", ""), 
                sprintf(" (n = %d)", nrow(age_data_list[[age_group]])), "\n"))
    }
  }
  
  # 打印整体数据范围
  cat(paste("整体范围:", sprintf("%.1f - %.1f", x_min, x_max), 
            ifelse(biomarker_type == "WC", "cm", ""), "\n"))
  
  return(list(plot = p, 
              predictions = pred_combined,
              p_values = p_values,
              actual_ranges = actual_ranges))
}


# 创建年龄数据列表
age_data_list <- list(
  "Higher Age" = higher_age_df,
  "Lower Age" = lower_age_df
)

# 运行所有分析（统一使用log2转换）
all_cause_wc_age <- rcs_analysis_age_stratified(all_cause_surv, "All-cause mortality", age_data_list, "WC", use_log2 = TRUE)
cvd_wc_age <- rcs_analysis_age_stratified(cvd_surv, "CVD mortality", age_data_list, "WC", use_log2 = TRUE)
cancer_wc_age <- rcs_analysis_age_stratified(cancer_surv, "Cancer mortality", age_data_list, "WC", use_log2 = TRUE)
diabetes_wc_age <- rcs_analysis_age_stratified(diabetes_surv, "Diabetes mortality", age_data_list, "WC", use_log2 = TRUE)

all_cause_prowc_age <- rcs_analysis_age_stratified(all_cause_surv, "All-cause mortality", age_data_list, "proWC", use_log2 = TRUE)
cvd_prowc_age <- rcs_analysis_age_stratified(cvd_surv, "CVD mortality", age_data_list, "proWC", use_log2 = TRUE)
cancer_prowc_age <- rcs_analysis_age_stratified(cancer_surv, "Cancer mortality", age_data_list, "proWC", use_log2 = TRUE)
diabetes_prowc_age <- rcs_analysis_age_stratified(diabetes_surv, "Diabetes mortality", age_data_list, "proWC", use_log2 = TRUE)

# 使用patchwork组合图表
library(patchwork)

# 移除图例
all_cause_prowc_age_nolegend <- all_cause_prowc_age$plot + theme(legend.position = "none")
cvd_prowc_age_nolegend <- cvd_prowc_age$plot + theme(legend.position = "none")
cancer_prowc_age_nolegend <- cancer_prowc_age$plot + theme(legend.position = "none")
diabetes_prowc_age_nolegend <- diabetes_prowc_age$plot + theme(legend.position = "none")

all_cause_wc_age_nolegend <- all_cause_wc_age$plot + theme(plot.title = element_blank(), legend.position = "none")
cvd_wc_age_nolegend <- cvd_wc_age$plot + theme(plot.title = element_blank(), legend.position = "none")
cancer_wc_age_nolegend <- cancer_wc_age$plot + theme(plot.title = element_blank(), legend.position = "none")
diabetes_wc_age_nolegend <- diabetes_wc_age$plot + theme(plot.title = element_blank(), legend.position = "none")

# 组合图表：proWC在上面，WC在下面
combined_plot_age <- (
  (all_cause_prowc_age_nolegend / all_cause_wc_age_nolegend) |  # proWC在上，WC在下
    (cvd_prowc_age_nolegend / cvd_wc_age_nolegend) |
    (cancer_prowc_age_nolegend / cancer_wc_age_nolegend) |
    (diabetes_prowc_age_nolegend / diabetes_wc_age_nolegend)
) + 
  plot_layout(ncol = 4) &
  theme(
    legend.position = "none",
    # 确保组合后的图表也有透明背景
    plot.background = element_rect(fill = "transparent", color = NA)
  )


# 添加共享图例
legend_age <- cowplot::get_legend(all_cause_prowc_age$plot + theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent", color = NA)))

# 最终组合
final_plot_age <- combined_plot_age / legend_age + 
  plot_layout(heights = c(0.9, 0.1)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

final_plot_age



# 顶部
library(patchwork)

# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend_age, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot_age <- top_panel / combined_plot_age + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot_age








# BMI
temp_df$BMI_group <- cut(temp_df$BMI,
                         breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
                         labels = c("Underweight", "NW", "OW", "OB"),
                         right = TRUE)

# 提取四个BMI分组的数据集
underweight_df <- temp_df[temp_df$BMI_group == "Underweight" & !is.na(temp_df$BMI_group), ]
normal_df <- temp_df[temp_df$BMI_group == "NW" & !is.na(temp_df$BMI_group), ] 
overweight_df <- temp_df[temp_df$BMI_group == "OW" & !is.na(temp_df$BMI_group), ]
obese_df <- temp_df[temp_df$BMI_group == "OB" & !is.na(temp_df$BMI_group), ]

# 创建生存对象函数
create_survival_objects <- function(data) {
  data$all_cause_surv <- with(data, Surv(survival_time, all_cause_death))
  data$cvd_surv <- with(data, Surv(survival_time, cvd_death))
  data$cancer_surv <- with(data, Surv(survival_time, cancer_death))
  data$diabetes_surv <- with(data, Surv(survival_time, diabetes_death))
  return(data)
}

# 为每个BMI分组创建生存对象
underweight_df <- create_survival_objects(underweight_df)
normal_df <- create_survival_objects(normal_df)
overweight_df <- create_survival_objects(overweight_df)
obese_df <- create_survival_objects(obese_df)

# 设置datadist
dd_underweight <- datadist(underweight_df)
dd_normal <- datadist(normal_df)
dd_overweight <- datadist(overweight_df)
dd_obese <- datadist(obese_df)

# BMI分组RCS分析函数 - 按实际范围画图，不限制X轴，排除Underweight
rcs_analysis_bmi_stratified <- function(surv_obj, outcome_name, bmi_data_list, biomarker_type = "WC", use_log2 = FALSE) {
  
  # 定义BMI分组名称和颜色 - 排除Underweight
  bmi_groups <- c("Normal", "Overweight", "Obese")  # 移除了Underweight
  bmi_colors <- c(
    "Normal" = "#9932CC",         # 深紫色  
    "Overweight" = "#20B2AA",     # 浅海绿色
    "Obese" = "#FF4500"           # 橙红色
  )
  
  # 不跳过任何分组（除了Underweight）
  skip_underweight <- TRUE  # Underweight不画
  
  # 存储模型和预测结果
  models <- list()
  predictions <- list()
  p_values <- list()
  actual_ranges <- list()
  
  # 为每个BMI分组拟合模型
  for (bmi_group in bmi_groups) {
    data <- bmi_data_list[[bmi_group]]
    
    # 检查是否有足够的数据
    if (nrow(data) < 10) {
      warning(paste("BMI group", bmi_group, "has insufficient data (<10 observations) for", outcome_name))
      next
    }
    
    # 检查是否有足够的死亡事件
    if (outcome_name == "Diabetes mortality") {
      diabetes_deaths <- sum(data$diabetes_death, na.rm = TRUE)
      if (diabetes_deaths < 5) {
        warning(paste("BMI group", bmi_group, "has insufficient diabetes deaths (<5) for", outcome_name))
        next
      }
    }
    
    # 设置datadist
    options(datadist = get(paste0("dd_", tolower(bmi_group))))
    
    # 构建公式
    if (biomarker_type == "WC") {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(WC, 4) + Age + Sex"))
    } else {
      formula <- as.formula(paste(deparse(substitute(surv_obj)), "~ rcs(BioX_Adjusted, 4) + Age + Sex"))
    }
    
    # 拟合模型
    models[[bmi_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
    
    # 计算非线性P值
    anova_result <- anova(models[[bmi_group]])
    p_values[[bmi_group]] <- ifelse(biomarker_type == "WC", 
                                    anova_result["WC", "P"], 
                                    anova_result["BioX_Adjusted", "P"])
    
    # 记录实际数据范围（最小值到最大值）
    if (biomarker_type == "WC") {
      actual_ranges[[bmi_group]] <- c(min(data$WC, na.rm = TRUE), max(data$WC, na.rm = TRUE))
    } else {
      actual_ranges[[bmi_group]] <- c(min(data$BioX_Adjusted, na.rm = TRUE), max(data$BioX_Adjusted, na.rm = TRUE))
    }
  }
  
  # 如果没有模型，返回空结果
  if (length(models) == 0) {
    warning(paste("No valid models for", outcome_name))
    return(NULL)
  }
  
  # 生成预测数据 - 使用每个BMI分组的实际范围（最小值到最大值）
  pred_combined <- data.frame()
  
  for (bmi_group in names(models)) {
    data <- bmi_data_list[[bmi_group]]
    options(datadist = get(paste0("dd_", tolower(bmi_group))))
    
    # 使用该BMI分组的实际数据范围（最小值到最大值）
    if (biomarker_type == "WC") {
      actual_range <- actual_ranges[[bmi_group]]
      biomarker_range <- seq(actual_range[1], actual_range[2], length = 100)
      
      pred <- Predict(models[[bmi_group]], 
                      WC = biomarker_range, 
                      Age = median(data$Age, na.rm = TRUE),
                      Sex = median(data$Sex, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$WC, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        BMI_Group = bmi_group
      )
    } else {
      actual_range <- actual_ranges[[bmi_group]]
      biomarker_range <- seq(actual_range[1], actual_range[2], length = 100)
      
      pred <- Predict(models[[bmi_group]], 
                      BioX_Adjusted = biomarker_range, 
                      Age = median(data$Age, na.rm = TRUE),
                      Sex = median(data$Sex, na.rm = TRUE),
                      ref.zero = TRUE,
                      fun = exp)
      
      pred_df <- data.frame(
        Biomarker = pred$BioX_Adjusted, 
        HR = pred$yhat, 
        lower = pred$lower, 
        upper = pred$upper,
        BMI_Group = bmi_group
      )
    }
    
    pred_combined <- rbind(pred_combined, pred_df)
  }
  
  # 如果使用log2转换
  if (use_log2) {
    pred_combined$HR <- log2(pred_combined$HR)
    pred_combined$lower <- log2(pred_combined$lower)
    pred_combined$upper <- log2(pred_combined$upper)
  }
  
  # 关键修改：将BMI_Group转换为因子，并指定正确的顺序
  pred_combined$BMI_Group <- factor(pred_combined$BMI_Group, 
                                    levels = c("Normal", "Overweight", "Obese"))
  
  # 计算整体的X轴范围（覆盖所有BMI分组）
  if (biomarker_type == "WC") {
    x_label <- "WC"
    # 计算整体的数据范围（使用实际预测范围）
    x_min <- min(pred_combined$Biomarker, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker, na.rm = TRUE)
  } else {
    x_label <- "proWC"
    # 计算整体的数据范围
    x_min <- min(pred_combined$Biomarker, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker, na.rm = TRUE)
  }
  
  # 计算合适的P值标注位置
  y_max <- max(pred_combined$upper, na.rm = TRUE)
  y_min <- min(pred_combined$lower, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 绘制图表 - 不设置X轴限制
  p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = BMI_Group, fill = BMI_Group)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_hline(yintercept = ifelse(use_log2, 0, 1), 
               linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = bmi_colors, 
                       breaks = c("Normal", "Overweight", "Obese")) +  # 确保颜色与因子顺序匹配
    scale_fill_manual(values = bmi_colors,
                      breaks = c("Normal", "Overweight", "Obese")) +   # 确保填充色与因子顺序匹配
    labs(
      title = outcome_name,
      x = x_label,
      y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio (95% CI)"),
      color = "BMI Group",
      fill = "BMI Group"
    ) +
    #scale_x_continuous(limits = c(50, 160)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),  # 全局字体
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Arial"),
      axis.title = element_text(size = 11, face = "bold", family = "Arial"),  # 改为11
      axis.text = element_text(size = 9, family = "Arial"),  # 改为9
      axis.line = element_line(size = 0.8, color = "black"),
      axis.ticks = element_line(size = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, face = "bold", family = "Arial"),  # 改为8
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  # 在图上添加P值
  if (length(p_values) > 0) {
    y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
    i <- 1
    for (bmi_group in names(p_values)) {
      p_value <- p_values[[bmi_group]]
      
      # 动态计算X位置（使用当前分组的数据范围）
      plot_data <- pred_combined[pred_combined$BMI_Group == bmi_group, ]
      if (nrow(plot_data) > 0) {
        x_min_group <- min(plot_data$Biomarker, na.rm = TRUE)
        x_max_group <- max(plot_data$Biomarker, na.rm = TRUE)
        x_pos <- x_min_group + 0.05 * (x_max_group - x_min_group)
        
        # 修复P值标注逻辑
        p_label <- if (p_value < 0.001) {
          "P-nonlinear ≤ 0.001"
        } else if (p_value < 0.01) {
          paste0("P-nonlinear = ", sprintf("%.3f", p_value))
        } else if (p_value < 0.05) {
          paste0("P-nonlinear = ", sprintf("%.3f", p_value))
        } else {
          paste0("P-nonlinear = ", sprintf("%.3f", p_value))
        }
        
        p <- p + annotate("text", 
                          x = x_pos, 
                          y = y_max - y_positions[i] * y_range,
                          label = p_label,
                          size = 2.8, fontface = "bold", hjust = 0,
                          color = bmi_colors[bmi_group],
                          family = "Arial")
        i <- i + 1
      }
    }
  }
  
  # 打印实际数据范围信息
  cat(paste("\n", outcome_name, "实际数据范围（", ifelse(biomarker_type == "WC", "WC", "proWC"), "）:\n"))
  for (bmi_group in names(actual_ranges)) {
    if (bmi_group %in% names(models)) {
      range <- actual_ranges[[bmi_group]]
      cat(paste(bmi_group, ":", sprintf("%.1f - %.1f", range[1], range[2]), 
                ifelse(biomarker_type == "WC", "cm", ""), 
                sprintf(" (n = %d)", nrow(bmi_data_list[[bmi_group]])), "\n"))
    }
  }
  
  # 打印整体数据范围
  cat(paste("整体范围:", sprintf("%.1f - %.1f", x_min, x_max), 
            ifelse(biomarker_type == "WC", "cm", ""), "\n"))
  
  return(list(plot = p, 
              predictions = pred_combined,
              p_values = p_values,
              actual_ranges = actual_ranges))
}

# 创建BMI数据列表
bmi_data_list <- list(
  Underweight = underweight_df,
  Normal = normal_df,
  Overweight = overweight_df,
  Obese = obese_df
)

# 运行所有分析（统一使用log2转换）
all_cause_wc_bmi <- rcs_analysis_bmi_stratified(all_cause_surv, "All-cause mortality", bmi_data_list, "WC", use_log2 = TRUE)
cvd_wc_bmi <- rcs_analysis_bmi_stratified(cvd_surv, "CVD mortality", bmi_data_list, "WC", use_log2 = TRUE)
cancer_wc_bmi <- rcs_analysis_bmi_stratified(cancer_surv, "Cancer mortality", bmi_data_list, "WC", use_log2 = TRUE)
diabetes_wc_bmi <- rcs_analysis_bmi_stratified(diabetes_surv, "Diabetes mortality", bmi_data_list, "WC", use_log2 = TRUE)

all_cause_prowc_bmi <- rcs_analysis_bmi_stratified(all_cause_surv, "All-cause mortality", bmi_data_list, "proWC", use_log2 = TRUE)
cvd_prowc_bmi <- rcs_analysis_bmi_stratified(cvd_surv, "CVD mortality", bmi_data_list, "proWC", use_log2 = TRUE)
cancer_prowc_bmi <- rcs_analysis_bmi_stratified(cancer_surv, "Cancer mortality", bmi_data_list, "proWC", use_log2 = TRUE)
diabetes_prowc_bmi <- rcs_analysis_bmi_stratified(diabetes_surv, "Diabetes mortality", bmi_data_list, "proWC", use_log2 = TRUE)

# 使用patchwork组合图表
library(patchwork)

# 移除图例
all_cause_prowc_bmi_nolegend <- all_cause_prowc_bmi$plot + theme(legend.position = "none")
cvd_prowc_bmi_nolegend <- cvd_prowc_bmi$plot + theme(legend.position = "none")
cancer_prowc_bmi_nolegend <- cancer_prowc_bmi$plot + theme(legend.position = "none")
diabetes_prowc_bmi_nolegend <- diabetes_prowc_bmi$plot + theme(legend.position = "none")

all_cause_wc_bmi_nolegend <- all_cause_wc_bmi$plot + theme(plot.title = element_blank(), legend.position = "none")
cvd_wc_bmi_nolegend <- cvd_wc_bmi$plot + theme(plot.title = element_blank(), legend.position = "none")
cancer_wc_bmi_nolegend <- cancer_wc_bmi$plot + theme(plot.title = element_blank(), legend.position = "none")
diabetes_wc_bmi_nolegend <- diabetes_wc_bmi$plot + theme(plot.title = element_blank(), legend.position = "none")

# 组合图表：proWC在上面，WC在下面
combined_plot_bmi <- (
  (all_cause_prowc_bmi_nolegend / all_cause_wc_bmi_nolegend) |  # proWC在上，WC在下
    (cvd_prowc_bmi_nolegend / cvd_wc_bmi_nolegend) |
    (cancer_prowc_bmi_nolegend / cancer_wc_bmi_nolegend) |
    (diabetes_prowc_bmi_nolegend / diabetes_wc_bmi_nolegend)
) + 
  plot_layout(ncol = 4) &
  theme(
    legend.position = "none",
    # 确保组合后的图表也有透明背景
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# 添加共享图例
legend_bmi <- cowplot::get_legend(all_cause_prowc_bmi$plot + theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent", color = NA)))

# 最终组合
final_plot_bmi <- combined_plot_bmi / legend_bmi + 
  plot_layout(heights = c(0.9, 0.1)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

final_plot_bmi


# 顶部
library(patchwork)

# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend_bmi, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot_bmi <- top_panel / combined_plot_bmi + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot_bmi

