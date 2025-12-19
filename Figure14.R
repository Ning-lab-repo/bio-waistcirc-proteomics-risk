setwd("D:/UKB_data")
data <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("WC_death_time.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df <- cbind(complete_data_WC, results)

# 直接匹配行名添加BMI列
temp_df$BMI <- data[rownames(temp_df), "BMI"]

wc_median <- median(temp_df$WC, na.rm = TRUE)
prowc_median <- median(temp_df$BioX_Adjusted, na.rm = TRUE)



# 安装必要的包
library(forestplot)
library(dplyr)
library(ggplot2)

temp_df <- temp_df %>%
  mutate(
    WC_proWC_group = case_when(
      WC <= wc_median & BioX_Adjusted <= prowc_median ~ "Lower WC + Lower proWC",
      WC <= wc_median & BioX_Adjusted > prowc_median ~ "Lower WC + Higher proWC", 
      WC > wc_median & BioX_Adjusted <= prowc_median ~ "Higher WC + Lower proWC",
      WC > wc_median & BioX_Adjusted > prowc_median ~ "Higher WC + Higher proWC",
      TRUE ~ NA_character_
    )
  )

# 将分组变量转换为因子，设置新的顺序
temp_df$WC_proWC_group <- factor(temp_df$WC_proWC_group,
                                 levels = c("Lower WC + Lower proWC",
                                            "Higher WC + Lower proWC",
                                            "Lower WC + Higher proWC", 
                                            "Higher WC + Higher proWC"))

# 更新数据清理
temp_df_clean <- temp_df[!is.na(temp_df$BMI), ]
male_data <- subset(temp_df, Sex == "1")
female_data <- subset(temp_df, Sex == "0")

# 修改分析函数来处理四种组合
calculate_group_hr <- function(data, outcome_var, time_var) {
  # 创建生存对象
  surv_obj <- with(data, Surv(get(time_var), get(outcome_var)))
  
  # 构建Cox回归公式，以"Lower WC + Lower proWC"为参考
  formula <- as.formula(paste("surv_obj ~ WC_proWC_group + Age + BMI"))
  
  # 拟合Cox模型
  model <- coxph(formula, data = data)
  
  # 获取模型汇总
  model_summary <- summary(model)
  
  # 提取所有组的HR和置信区间
  all_groups <- levels(data$WC_proWC_group)[-1] # 排除参考组
  
  results <- list()
  for (group in all_groups) {
    coef_name <- paste0("WC_proWC_group", group)
    if (coef_name %in% rownames(model_summary$coefficients)) {
      coef_info <- model_summary$coefficients[coef_name, ]
      conf_int <- model_summary$conf.int[coef_name, ]
      
      results[[group]] <- list(
        hr = conf_int["exp(coef)"],
        lower = conf_int["lower .95"],
        upper = conf_int["upper .95"],
        p_value = coef_info["Pr(>|z|)"]
      )
    }
  }
  
  return(list(model = model, results = results))
}

create_forest_data_combined <- function(male_data, female_data) {
  outcomes <- list(
    "All-cause mortality" = list(outcome = "all_cause_death", time = "survival_time"),
    "CVD mortality" = list(outcome = "cvd_death", time = "survival_time"), 
    "Cancer mortality" = list(outcome = "cancer_death", time = "survival_time"),
    "Diabetes mortality" = list(outcome = "diabetes_death", time = "survival_time")
  )
  
  sex_groups <- c("Male", "Female")
  # 使用新的分组顺序
  all_groups <- c("Lower WC + Lower proWC",  # 参考组
                  "Higher WC + Lower proWC",
                  "Lower WC + Higher proWC", 
                  "Higher WC + Higher proWC")
  
  all_results <- data.frame()
  
  for (outcome_name in names(outcomes)) {
    outcome_info <- outcomes[[outcome_name]]
    
    for (sex in sex_groups) {
      # 选择数据
      if (sex == "Male") {
        data_use <- male_data
      } else {
        data_use <- female_data
      }
      
      # 确保数据中有足够样本
      if (nrow(data_use) > 50) {
        tryCatch({
          cox_result <- calculate_group_hr(
            data_use, 
            outcome_info$outcome, 
            outcome_info$time
          )
          
          # 为每个组创建结果，包括参考组
          for (group_name in all_groups) {
            if (group_name == "Lower WC + Lower proWC") {
              # 参考组：HR=1，无置信区间
              result <- list(
                hr = 1.0,
                lower = NA,
                upper = NA,
                p_value = NA
              )
            } else if (group_name %in% names(cox_result$results)) {
              result <- cox_result$results[[group_name]]
            } else {
              next  # 跳过不存在的组
            }
            
            all_results <- rbind(all_results, data.frame(
              outcome = outcome_name,
              sex = sex,
              group = group_name,
              hr = result$hr,
              lower = result$lower,
              upper = result$upper,
              p_value = result$p_value,
              stringsAsFactors = FALSE
            ))
          }
        }, error = function(e) {
          message(paste("Error in", outcome_name, sex, ":", e$message))
        })
      }
    }
  }
  
  return(all_results)
}

# 生成数据
forest_data_complete <- create_forest_data_combined(male_data, female_data)

# 确保分组顺序正确
forest_data_complete$group <- factor(forest_data_complete$group,
                                     levels = c("Lower WC + Lower proWC",
                                                "Higher WC + Lower proWC",
                                                "Lower WC + Higher proWC", 
                                                "Higher WC + Higher proWC"))

# 数据处理和可视化
forest_data_complete <- forest_data_complete %>%
  mutate(
    Significant = case_when(
      is.na(p_value) ~ "Not significant",  # 参考组归为Not significant
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

forest_data_complete$Significant <- factor(forest_data_complete$Significant, 
                                           levels = c("Significant", "Not significant"))

# 绘图
p1 <- ggplot(forest_data_complete, aes(x = hr, y = group, color = sex, shape = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = sex),
    height = 0.3,
    linewidth = 1,
    position = position_dodge(width = 0.7)
  ) +
  geom_point(
    aes(fill = sex, shape = Significant),
    size = 2.5,
    stroke = 1,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  
  # 使用rev()函数翻转y轴，让"Lower WC + Lower proWC"在顶部
  scale_y_discrete(limits = rev) +
  
  scale_color_manual(
    values = c("Male" = "#1f77B4", "Female" = "#ff7f0E"),
    name = ""
  ) +
  scale_fill_manual(
    values = c("Male" = "#a6c8e3", "Female" = "#ffd1a4"),
    name = "Sex"
  ) +
  scale_shape_manual(
    values = c("Significant" = 24, "Not significant" = 21),
    name = ""
  ) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    # 设置透明背景
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 调整面板边框 - 移除或调整以不重合
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    # 移除axis.line，因为panel.border已经提供了边框
    axis.line = element_blank(),  # 移除坐标轴线，避免重合
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    axis.text = element_text(color = "black", size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    # 设置facet背景
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    # 设置图例背景为透明
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(face = "bold", size = 8, family = "Arial"),
    plot.margin = margin(0, 5, 1, 0)
  ) +
  guides(
    color = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        shape = 21,
        fill = c("#ffd1a4", "#a6c8e3"),
        color = c("#ff7f0E", "#1f77B4")
      )
    ),
    shape = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        fill = "gray80",
        color = "black"
      )
    ),
    fill = "none"
  )





temp_df$BMI_group <- cut(temp_df$BMI,
                         breaks = c(18.5, 24.9, 29.9, Inf),
                         labels = c("Normal", "Overweight", "Obese"),
                         right = TRUE, include.lowest = FALSE)

# 排除BMI < 18.5的数据
temp_df <- temp_df[!is.na(temp_df$BMI_group), ]

# 提取三个BMI分组的数据集（排除Underweight）
normal_df <- temp_df[temp_df$BMI_group == "Normal" & !is.na(temp_df$BMI_group), ] 
overweight_df <- temp_df[temp_df$BMI_group == "Overweight" & !is.na(temp_df$BMI_group), ]
obese_df <- temp_df[temp_df$BMI_group == "Obese" & !is.na(temp_df$BMI_group), ]

# 修改calculate_group_hr_bmi函数（确保有参考组）
calculate_group_hr_bmi <- function(data, outcome_var, time_var) {
  # 确保分组变量存在且有序
  if (!"WC_proWC_group" %in% names(data)) {
    stop("WC_proWC_group variable not found in data")
  }
  
  # 设置正确的因子水平顺序
  data$WC_proWC_group <- factor(data$WC_proWC_group,
                                levels = c("Lower WC + Lower proWC",
                                           "Higher WC + Lower proWC",
                                           "Lower WC + Higher proWC", 
                                           "Higher WC + Higher proWC"))
  
  # 创建生存对象
  surv_obj <- with(data, Surv(get(time_var), get(outcome_var)))
  
  # 构建Cox回归公式，以"Lower WC + Lower proWC"为参考
  formula <- as.formula(paste("surv_obj ~ WC_proWC_group + Age + Sex"))
  
  # 拟合Cox模型
  model <- coxph(formula, data = data)
  
  # 获取模型汇总
  model_summary <- summary(model)
  
  # 提取所有组的结果，包括参考组
  all_groups <- levels(data$WC_proWC_group)
  
  results <- list()
  for (group in all_groups) {
    if (group == "Lower WC + Lower proWC") {
      # 参考组：HR=1
      results[[group]] <- list(
        hr = 1.0,
        lower = NA,
        upper = NA,
        p_value = NA
      )
    } else {
      coef_name <- paste0("WC_proWC_group", group)
      if (coef_name %in% rownames(model_summary$coefficients)) {
        coef_info <- model_summary$coefficients[coef_name, ]
        conf_int <- model_summary$conf.int[coef_name, ]
        
        results[[group]] <- list(
          hr = conf_int["exp(coef)"],
          lower = conf_int["lower .95"],
          upper = conf_int["upper .95"],
          p_value = coef_info["Pr(>|z|)"]
        )
      } else {
        # 如果该组在模型中不存在（可能由于样本量不足）
        results[[group]] <- list(
          hr = NA,
          lower = NA,
          upper = NA,
          p_value = NA
        )
      }
    }
  }
  
  return(list(model = model, results = results))
}

create_forest_data_bmi_combined <- function(normal_df, overweight_df, obese_df) {
  outcomes <- list(
    "All-cause mortality" = list(outcome = "all_cause_death", time = "survival_time"),
    "CVD mortality" = list(outcome = "cvd_death", time = "survival_time"), 
    "Cancer mortality" = list(outcome = "cancer_death", time = "survival_time"),
    "Diabetes mortality" = list(outcome = "diabetes_death", time = "survival_time")
  )
  
  bmi_groups <- list(
    "Normal" = normal_df,
    "Overweight" = overweight_df,
    "Obese" = obese_df
  )
  
  all_results <- data.frame()
  
  for (outcome_name in names(outcomes)) {
    outcome_info <- outcomes[[outcome_name]]
    
    for (bmi_group_name in names(bmi_groups)) {
      data_use <- bmi_groups[[bmi_group_name]]
      
      # 移除样本量不足的判断条件，让所有组合都尝试计算
      # 确保数据中有一定样本
      if (nrow(data_use) > 20) {  # 降低阈值到20，以包含更多组合
        tryCatch({
          cox_result <- calculate_group_hr_bmi(
            data_use, 
            outcome_info$outcome, 
            outcome_info$time
          )
          
          # 提取每个组的结果（包括参考组）
          for (group_name in names(cox_result$results)) {
            result <- cox_result$results[[group_name]]
            
            all_results <- rbind(all_results, data.frame(
              outcome = outcome_name,
              bmi_group = bmi_group_name,
              group = group_name,
              hr = result$hr,
              lower = result$lower,
              upper = result$upper,
              p_value = result$p_value,
              stringsAsFactors = FALSE
            ))
          }
        }, error = function(e) {
          # 如果模型拟合失败，仍然添加参考组
          if (outcome_name == "Diabetes mortality" && bmi_group_name == "Normal") {
            # 为Normal组的Diabetes mortality添加参考组
            for (group_name in c("Lower WC + Lower proWC")) {
              all_results <- rbind(all_results, data.frame(
                outcome = outcome_name,
                bmi_group = bmi_group_name,
                group = group_name,
                hr = 1.0,
                lower = NA,
                upper = NA,
                p_value = NA,
                stringsAsFactors = FALSE
              ))
            }
          }
          message(paste("Error in", outcome_name, bmi_group_name, ":", e$message))
        })
      }
    }
  }
  
  return(all_results)
}

# 生成数据 - 只传入三个BMI分组
forest_data_bmi_combined <- create_forest_data_bmi_combined(normal_df, overweight_df, obese_df)

# 确保组的顺序正确
forest_data_bmi_combined$group <- factor(forest_data_bmi_combined$group,
                                         levels = c("Lower WC + Lower proWC",
                                                    "Higher WC + Lower proWC",
                                                    "Lower WC + Higher proWC", 
                                                    "Higher WC + Higher proWC"))

# 确保BMI分组的顺序
forest_data_bmi_combined$bmi_group <- factor(forest_data_bmi_combined$bmi_group,
                                             levels = c("Normal", "Overweight", "Obese"))

# 数据处理和可视化
forest_data_bmi_combined <- forest_data_bmi_combined %>%
  mutate(
    Significant = case_when(
      is.na(p_value) ~ "Not significant",  # 参考组归为Not significant
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

forest_data_bmi_combined$Significant <- factor(forest_data_bmi_combined$Significant, 
                                               levels = c("Significant", "Not significant"))

# 绘图 - 按BMI分组显示四种WC-proWC组合
p2 <- ggplot(forest_data_bmi_combined, aes(x = hr, y = group, color = bmi_group, shape = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = bmi_group),
    height = 0.3,
    linewidth = 1,
    position = position_dodge(width = 0.7)
  ) +
  geom_point(
    aes(fill = bmi_group, shape = Significant),
    size = 2.5,
    stroke = 1,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  
  # 使用rev()函数翻转y轴，让"Lower WC + Lower proWC"在顶部
  scale_y_discrete(limits = rev) +
  
  scale_color_manual(
    values = c("Normal" = "#9932CC",         # darkorchid (深紫色)
               "Overweight" = "#20B2AA",     # lightseagreen (浅海绿色)
               "Obese" = "#FF4500"),         # orangered (橙红色)
    name = ""
  ) +
  scale_fill_manual(
    values = c("Normal" = "#D8BFD8",         # thistle (浅紫色)
               "Overweight" = "#AFEEEE",     # paleturquoise (浅海绿色)
               "Obese" = "#FFB6C1"),         # lightpink (浅粉色)
    name = "BMI Group"
  ) +
  scale_shape_manual(
    values = c("Significant" = 24, "Not significant" = 21),
    name = ""
  ) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    # 设置透明背景
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 调整面板边框
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    # 移除axis.line
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    axis.text = element_text(color = "black", size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    # 设置facet背景
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    # 设置图例背景为透明
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(face = "bold", size = 8, family = "Arial"),
    plot.margin = margin(0, 5, 1, 0)
  ) +
  guides(
    color = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        shape = 21,
        fill = c("#D8BFD8", "#AFEEEE", "#FFB6C1"),  # 浅色填充
        color = c("#9932CC", "#20B2AA", "#FF4500")  # 深色边框
      )
    ),
    shape = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        fill = "gray80",
        color = "black"
      )
    ),
    fill = "none"
  )

# 显示图形
print(p2)




# Age
# 计算中位年龄
temp_df$Sex <- as.factor(temp_df$Sex)
age_median <- median(temp_df$Age)

# 创建分组变量
temp_df$age_group <- ifelse(temp_df$Age >= age_median, "Higher Age", "Lower Age")

higher_age_df <- temp_df[temp_df$age_group == "Higher Age", ]
lower_age_df <- temp_df[temp_df$age_group == "Lower Age", ]

# 修改分析函数来处理四种组合（按年龄分组）
calculate_group_hr_age <- function(data, outcome_var, time_var) {
  # 首先确保WC_proWC_group有正确的因子水平顺序
  data$WC_proWC_group <- factor(data$WC_proWC_group,
                                levels = c("Lower WC + Lower proWC",
                                           "Higher WC + Lower proWC",
                                           "Lower WC + Higher proWC", 
                                           "Higher WC + Higher proWC"))
  
  # 创建生存对象
  surv_obj <- with(data, Surv(get(time_var), get(outcome_var)))
  
  # 构建Cox回归公式，以"Lower WC + Lower proWC"为参考
  formula <- as.formula(paste("surv_obj ~ WC_proWC_group + Sex + BMI"))
  
  # 拟合Cox模型
  model <- coxph(formula, data = data)
  
  # 获取模型汇总
  model_summary <- summary(model)
  
  # 提取所有组的HR和置信区间（包括参考组）
  all_groups <- levels(data$WC_proWC_group)
  
  results <- list()
  for (group in all_groups) {
    if (group == "Lower WC + Lower proWC") {
      # 参考组：HR = 1, 置信区间为[1,1]
      results[[group]] <- list(
        hr = 1.0,
        lower = 1.0,
        upper = 1.0,
        p_value = NA
      )
    } else {
      coef_name <- paste0("WC_proWC_group", group)
      if (coef_name %in% rownames(model_summary$coefficients)) {
        coef_info <- model_summary$coefficients[coef_name, ]
        conf_int <- model_summary$conf.int[coef_name, ]
        
        results[[group]] <- list(
          hr = conf_int["exp(coef)"],
          lower = conf_int["lower .95"],
          upper = conf_int["upper .95"],
          p_value = coef_info["Pr(>|z|)"]
        )
      } else {
        # 如果组不存在于模型中（可能由于样本量不足）
        results[[group]] <- list(
          hr = NA,
          lower = NA,
          upper = NA,
          p_value = NA
        )
      }
    }
  }
  
  return(list(model = model, results = results))
}

# 数据生成函数（按年龄分组）
create_forest_data_age_combined <- function(higher_age_df, lower_age_df) {
  outcomes <- list(
    "All-cause mortality" = list(outcome = "all_cause_death", time = "survival_time"),
    "CVD mortality" = list(outcome = "cvd_death", time = "survival_time"), 
    "Cancer mortality" = list(outcome = "cancer_death", time = "survival_time"),
    "Diabetes mortality" = list(outcome = "diabetes_death", time = "survival_time")
  )
  
  age_groups <- list(
    "Lower Age" = lower_age_df,
    "Higher Age" = higher_age_df
  )
  
  all_results <- data.frame()
  
  for (outcome_name in names(outcomes)) {
    outcome_info <- outcomes[[outcome_name]]
    
    for (age_group_name in names(age_groups)) {
      data_use <- age_groups[[age_group_name]]
      
      # 确保数据中有足够样本
      if (nrow(data_use) > 50) {
        tryCatch({
          cox_result <- calculate_group_hr_age(
            data_use, 
            outcome_info$outcome, 
            outcome_info$time
          )
          
          # 提取每个组的结果（包括参考组）
          for (group_name in names(cox_result$results)) {
            result <- cox_result$results[[group_name]]
            
            all_results <- rbind(all_results, data.frame(
              outcome = outcome_name,
              age_group = age_group_name,
              group = group_name,
              hr = result$hr,
              lower = result$lower,
              upper = result$upper,
              p_value = result$p_value,
              stringsAsFactors = FALSE
            ))
          }
        }, error = function(e) {
          message(paste("Error in", outcome_name, age_group_name, ":", e$message))
        })
      }
    }
  }
  
  return(all_results)
}

# 生成数据
forest_data_age_combined <- create_forest_data_age_combined(higher_age_df, lower_age_df)

# 确保组的顺序正确（按照要求的顺序）
forest_data_age_combined$group <- factor(forest_data_age_combined$group,
                                         levels = c("Lower WC + Lower proWC",
                                                    "Higher WC + Lower proWC",
                                                    "Lower WC + Higher proWC", 
                                                    "Higher WC + Higher proWC"))

# 确保年龄分组的顺序
forest_data_age_combined$age_group <- factor(forest_data_age_combined$age_group,
                                             levels = c("Lower Age", "Higher Age"))

# 数据处理和可视化
forest_data_age_combined <- forest_data_age_combined %>%
  mutate(
    Significant = case_when(
      is.na(p_value) ~ "Not significant",  # 参考组归为Not significant
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

forest_data_age_combined$Significant <- factor(forest_data_age_combined$Significant, 
                                               levels = c("Significant", "Not significant"))

# 绘图 - 按年龄分组显示四种WC-proWC组合
p3 <- ggplot(forest_data_age_combined, aes(x = hr, y = group, color = age_group, shape = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = age_group),
    height = 0.3,
    linewidth = 1,
    position = position_dodge(width = 0.7)
  ) +
  geom_point(
    aes(fill = age_group, shape = Significant),
    size = 2.5,
    stroke = 1,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  
  # 使用rev()函数翻转y轴，让"Lower WC + Lower proWC"在顶部
  scale_y_discrete(limits = rev) +
  
  scale_color_manual(
    values = c("Higher Age" = "#F28FB1", 
               "Lower Age" = "#9FB4C7"),
    name = ""
  ) +
  scale_fill_manual(
    values = c("Higher Age" = "#F8D1DC", 
               "Lower Age" = "#CFD8E3"),
    name = "Age Group"
  ) +
  scale_shape_manual(
    values = c("Significant" = 24, "Not significant" = 21),
    name = ""
  ) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    # 设置透明背景
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 调整面板边框 - 移除或调整以不重合
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    # 移除axis.line，因为panel.border已经提供了边框
    axis.line = element_blank(),  # 移除坐标轴线，避免重合
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    axis.text = element_text(color = "black", size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    # 设置facet背景
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    # 设置图例背景为透明
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(face = "bold", size = 8, family = "Arial"),
    plot.margin = margin(0, 5, 1, 0)
  ) +
  guides(
    color = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        shape = 21,
        fill = c("#CFD8E3", "#F8D1DC"),  # 注意顺序：先Lower Age，后Higher Age
        color = c("#9FB4C7", "#F28FB1")
      )
    ),
    shape = guide_legend(
      title.hjust = 0.5,
      override.aes = list(
        size = 3,
        fill = "gray80",
        color = "black"
      )
    ),
    fill = "none"
  )

# 显示图形
print(p3)


