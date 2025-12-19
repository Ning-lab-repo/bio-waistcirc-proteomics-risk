disease_codes <- c("E78", "I10", "K57", "K21", "I25")

setwd("D:/UKB_data")
# 1. 读取基础数据
data <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

date <- read.csv("2Summary_diagnoses_part1_participant.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# 假设date是一个数据框，且Diagnoses...ICD10是其中的一个字符列
date$Diagnoses...ICD10 <- sapply(date$Diagnoses...ICD10, function(x) {
  if (is.na(x) || x == "") {
    return(x)
  }
  
  # 分割每个编码
  codes <- unlist(strsplit(x, "\\|"))
  
  # 处理每个编码
  processed_codes <- sapply(codes, function(code) {
    if (nchar(code) <= 3) {
      return(code)
    } else {
      # 在前两位数字后加上小数点
      paste0(substr(code, 1, 3), ".", substr(code, 4, nchar(code)))
    }
  })
  
  # 重新组合编码
  paste(processed_codes, collapse = "|")
})



setwd("D:/UKB_data/NC")
complete_data <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
complete_data$BMI <- data[rownames(complete_data), "BMI"]
results <- read.csv("lasso_WC.csv", header = TRUE, row.names = 1)

# 2. 获取无疾病对照组的ID
no_disease_data <- data[str_trim(data$Diagnoses...ICD10) == "", ]
normal_ids <- rownames(no_disease_data)


# 对每个疾病运行分析
disease_results <- list()

setwd("D:/UKB_data/ID")

analyze_disease_risk <- function(disease_code) {
  # 前面的数据准备部分保持不变
  disease_ids <- read.csv(paste0("participant_ids_", disease_code, ".csv"))
  rownames(disease_ids) <- disease_ids$Participant_ID
  match_data <- data[rownames(disease_ids), ]
  disease_date <- date[rownames(disease_ids), ]
  
  disease_results <- lapply(1:nrow(disease_date), function(i) {
    diagnoses <- disease_date$Diagnoses...ICD10[i]
    diag_list <- strsplit(diagnoses, "\\|")[[1]]
    
    if (length(diag_list) == 0 || all(diag_list == "")) return(NULL)
    
    matches <- sapply(diag_list, function(x) any(startsWith(x, disease_code), na.rm = TRUE))
    if (all(!matches)) return(NULL)
    
    positions <- which(matches)
    result <- lapply(positions, function(pos) {
      date_col <- paste0("Date.of.first.in.patient.diagnosis...ICD10...Array.", pos - 1)
      if (date_col %in% names(disease_date)) {
        date_val <- disease_date[i, date_col]
        if (!is.na(date_val) && date_val != "") as.character(date_val) else NA
      } else {
        NA
      }
    })
    
    data.frame(
      Participant_ID = rep(rownames(disease_date)[i], length(positions)),
      Disease_Date = unlist(result),
      stringsAsFactors = FALSE
    )
  })
  
  first_date <- do.call(rbind, disease_results)
  if (is.null(first_date)) return(NULL)
  first_date <- first_date %>%
    mutate(Disease_Date = as.Date(Disease_Date)) %>%
    group_by(Participant_ID) %>%
    filter(Disease_Date == min(Disease_Date)) %>%
    sample_n(1) %>%  
    ungroup()
  rownames(first_date) <- first_date$Participant_ID
  
  start_date <- data.frame(
    date_attending_assessment_centre = match_data[rownames(first_date), "date_attending_assessment_centre"],
    row.names = rownames(first_date)
  )
  merge_date <- cbind(first_date, start_date)
  
  date_filtered <- merge_date %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_disease = as.Date(Disease_Date, format = "%Y/%m/%d"),
      onset_time_days = as.numeric(difftime(date_disease, date_assessment, units = "days")),
      onset_time_years = as.numeric(difftime(date_disease, date_assessment, units = "days")) / 365.25
    ) %>%
    filter(date_disease > date_assessment) %>%
    select(-date_assessment, -date_disease)
  
  if (nrow(date_filtered) == 0) return(NULL)
  
  # 病例数据
  case_ids <- rownames(date_filtered)
  results_case <- results[case_ids, ]
  complete_case <- complete_data[case_ids, ]
  disease_df <- cbind(complete_case, results_case, onset_time_years = date_filtered[case_ids, "onset_time_years"])
  disease_df$status <- 1  # 事件发生
  
  # 对照数据
  results_control <- results[normal_ids, ]
  complete_control <- complete_data[normal_ids, ]
  
  control_assessment_dates <- data.frame(
    date_attending_assessment_centre = data[normal_ids, "date_attending_assessment_centre"]
  )
  control_df <- cbind(complete_control, results_control, control_assessment_dates)
  
  control_df <- control_df %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_death = ifelse(is.na(data$new2025516_dead_data[normal_ids]) | data$new2025516_dead_data[normal_ids] == "", 
                          NA, 
                          as.Date(data$new2025516_dead_data[normal_ids], format = "%Y/%m/%d")),
      date_censored = as.Date("2024-07-08"),
      date_end = case_when(
        !is.na(date_death) & date_death <= date_censored ~ date_death,
        TRUE ~ date_censored
      ),
      onset_time_years = as.numeric(difftime(date_end, date_assessment, units = "days")) / 365.25,
      status = 0
    ) %>%
    select(-date_attending_assessment_centre, -date_assessment, -date_death, -date_censored, -date_end)
  
  # 合并病例和对照数据
  analysis_df <- rbind(disease_df, control_df)
  
  # 按性别分组
  male_data <- subset(analysis_df, Sex == "1")
  female_data <- subset(analysis_df, Sex == "0")
  
  # 加载必要的包
  library(survival)
  library(rms)
  library(ggplot2)
  
  # 创建生存对象函数
  create_survival_object <- function(data) {
    data$disease_surv <- with(data, Surv(onset_time_years, status))
    return(data)
  }
  
  male_data <- create_survival_object(male_data)
  female_data <- create_survival_object(female_data)
  
  # 设置datadist - 按照您要求的方式修改
  sex_data_list <- list(
    "male" = male_data,
    "female" = female_data
  )
  
  # 为每个性别分组创建datadist对象
  dd_male <- datadist(male_data)
  dd_female <- datadist(female_data)
  
  # RCS分析函数 - 按性别分层，支持WC和proWC
  rcs_analysis_disease <- function(sex_data_list, disease_name, biomarker_type = "WC", use_log2 = FALSE) {
    
    # 定义性别分组名称和颜色
    sex_groups <- c("male", "female")
    sex_colors <- c(
      "male" = "#1F77B4",    # 蓝色
      "female" = "#FF7F0E"   # 橙色
    )
    
    # 存储模型和预测结果
    models <- list()
    predictions <- list()
    p_values <- list()
    
    # 为每个性别分组拟合模型
    for (sex_group in sex_groups) {
      data <- sex_data_list[[sex_group]]
      
      # 设置datadist - 按照您要求的方式
      options(datadist = get(paste0("dd_", sex_group)))
      
      # 构建公式 - 根据生物标志物类型选择
      if (biomarker_type == "WC") {
        formula <- as.formula("disease_surv ~ rcs(Actual_WC, 4) + Age + BMI")
      } else {
        formula <- as.formula("disease_surv ~ rcs(BioX_Adjusted, 4) + Age + BMI")
      }
      
      # 拟合模型
      models[[sex_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
      
      # 计算非线性P值
      anova_result <- anova(models[[sex_group]])
      p_values[[sex_group]] <- ifelse(biomarker_type == "WC", 
                                      anova_result["Actual_WC", "P"], 
                                      anova_result["BioX_Adjusted", "P"])
    }
    
    # 生成预测数据 - 每个组使用自己的范围
    pred_combined <- data.frame()
    
    for (sex_group in names(models)) {
      data <- sex_data_list[[sex_group]]
      options(datadist = get(paste0("dd_", sex_group)))
      
      # 为每个组使用自己的实际范围
      if (biomarker_type == "WC") {
        wc_min <- min(data$Actual_WC, na.rm = TRUE)
        wc_max <- max(data$Actual_WC, na.rm = TRUE)
        biomarker_range <- seq(wc_min, wc_max, length = 100)
        
        pred <- Predict(models[[sex_group]], 
                        Actual_WC = biomarker_range, 
                        Age = median(data$Age, na.rm = TRUE),
                        BMI = median(data$BMI, na.rm = TRUE),
                        ref.zero = TRUE,
                        fun = exp)
        
        pred_df <- data.frame(
          Biomarker = pred$Actual_WC, 
          HR = pred$yhat, 
          lower = pred$lower, 
          upper = pred$upper,
          Sex_Group = sex_group,
          Biomarker_Range_Min = wc_min,
          Biomarker_Range_Max = wc_max
        )
      } else {
        prowc_min <- min(data$BioX_Adjusted, na.rm = TRUE)
        prowc_max <- max(data$BioX_Adjusted, na.rm = TRUE)
        biomarker_range <- seq(prowc_min, prowc_max, length = 100)
        
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
          Sex_Group = sex_group,
          Biomarker_Range_Min = prowc_min,
          Biomarker_Range_Max = prowc_max
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
    
    # 绘制图表 - 字体增大版本
    p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = Sex_Group, fill = Sex_Group)) +
      geom_line(size = 1.2) +  # 线条适当加粗
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
      geom_hline(yintercept = ifelse(use_log2, 0, 1), 
                 linetype = "dashed", color = "black", size = 1) +  # 参考线保持原样
      scale_color_manual(values = sex_colors, labels = c("male" = "Male", "female" = "Female")) +
      scale_fill_manual(values = sex_colors, labels = c("male" = "Male", "female" = "Female")) +
      labs(
        title = disease_name,
        x = ifelse(biomarker_type == "WC", "WC", "proWC"),
        y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio"),
        color = "Sex Group",
        fill = "Sex Group"
      ) +
      theme_minimal() +
      theme(
        text = element_text(family = "Arial"),  # 全局字体
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Arial"),  # 标题字体增大到16
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),  # 坐标轴标题增大到14
        axis.title.x = element_text(size = 14, face = "bold", family = "Arial"),  # x轴标题
        axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),  # y轴标题
        axis.text = element_text(size = 12, family = "Arial"),  # 坐标轴刻度文本增大到12
        axis.text.x = element_text(size = 12, family = "Arial"),  # x轴刻度
        axis.text.y = element_text(size = 12, family = "Arial"),  # y轴刻度
        axis.line = element_line(size = 0.8, color = "black"),  # 坐标轴线保持原样
        axis.ticks = element_line(size = 0.8, color = "black"),  # 刻度线保持原样
        axis.ticks.length = unit(0.2, "cm"),  # 刻度长度保持原样
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Arial"),  # 图例文本增大到12
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
        # 不修改整体边距，保持原样
      ) +
      # 设置x轴范围，让每个组的范围显示完整
      coord_cartesian(xlim = c(
        min(pred_combined$Biomarker_Range_Min, na.rm = TRUE),
        max(pred_combined$Biomarker_Range_Max, na.rm = TRUE)
      ))
    
    # 在图上添加P值 - 保持原样
    y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
    i <- 1
    for (sex_group in names(p_values)) {
      p_value <- p_values[[sex_group]]
      
      # 获取该组的数据用于标注位置
      group_data <- pred_combined[pred_combined$Sex_Group == sex_group, ]
      group_x_min <- min(group_data$Biomarker, na.rm = TRUE)
      group_x_max <- max(group_data$Biomarker, na.rm = TRUE)
      
      # 修正P值显示逻辑
      p_value_label <- if (p_value < 0.001) {
        "≤ 0.001"
      } else if (p_value < 0.01) {
        sprintf("= %.3f", p_value)
      } else if (p_value < 0.05) {
        sprintf("= %.3f", p_value)
      } else {
        sprintf("= %.3f", p_value)
      }
      
      p <- p + annotate("text", 
                        x = group_x_min + 0.1 * (group_x_max - group_x_min), 
                        y = y_max - y_positions[i] * y_range,
                        label = paste("P-nonlinear", p_value_label),
                        size = 3.2,  # P值标注保持原样
                        fontface = "bold", 
                        hjust = 0,
                        color = sex_colors[sex_group],
                        family = "Arial")
      
      i <- i + 1
    }
    
    return(list(plot = p, 
                predictions = pred_combined,
                p_values = p_values,
                models = models,
                biomarker_type = biomarker_type,
                biomarker_range_min = min(pred_combined$Biomarker_Range_Min, na.rm = TRUE),
                biomarker_range_max = max(pred_combined$Biomarker_Range_Max, na.rm = TRUE)))
  }
  
  # 运行RCS分析 - 分别分析WC和proWC
  disease_name <- disease_code
  
  # 分析WC
  rcs_result_wc <- rcs_analysis_disease(sex_data_list, disease_name, "WC", use_log2 = TRUE)
  
  # 分析proWC
  rcs_result_prowc <- rcs_analysis_disease(sex_data_list, disease_name, "proWC", use_log2 = TRUE)
  
  # 返回所有结果
  return(list(
    wc_plot = rcs_result_wc$plot,
    wc_predictions = rcs_result_wc$predictions,
    wc_p_values = rcs_result_wc$p_values,
    wc_models = rcs_result_wc$models,
    
    proWC_plot = rcs_result_prowc$plot,
    proWC_predictions = rcs_result_prowc$predictions,
    proWC_p_values = rcs_result_prowc$p_values,
    proWC_models = rcs_result_prowc$models,
    
    disease_name = disease_name,
    sample_info = data.frame(
      n_cases = nrow(disease_df),
      n_controls = nrow(control_df),
      total_n = nrow(analysis_df),
      n_male = nrow(male_data),
      n_female = nrow(female_data)
    )
  ))
}


significant_diseases <- list()
failed_diseases <- list()

for (disease_code in disease_codes) {
  cat("Analyzing disease:", disease_code, "\n")
  
  tryCatch({
    result <- analyze_disease_risk(disease_code)
    
    if (!is.null(result)) {
      disease_results[[disease_code]] <- result
      
      # 检查是否有错误信息
      if (!is.null(result$error)) {
        cat("  -> Model error:", result$error, "\n")
        failed_diseases[[disease_code]] <- result$error
        next
      }
      
      wc_p_vals <- unlist(result$wc_p_values)
      proWC_p_vals <- unlist(result$proWC_p_values)
      
      # 只考虑成功拟合的组
      valid_wc_p <- wc_p_vals[!is.na(wc_p_vals)]
      valid_proWC_p <- proWC_p_vals[!is.na(proWC_p_vals)]
      
      if (length(valid_wc_p) > 0 && length(valid_proWC_p) > 0 &&
          all(valid_wc_p < 0.05) && all(valid_proWC_p < 0.05)) {
        significant_diseases[[disease_code]] <- result
        cat("  -> Significant disease found!\n")
      }
    }
  }, error = function(e) {
    cat("  -> Fatal error:", e$message, "\n")
    failed_diseases[[disease_code]] <- e$message
  })
}

# 保存失败信息
if (length(failed_diseases) > 0) {
  failed_df <- data.frame(
    Disease_Code = names(failed_diseases),
    Error_Message = unlist(failed_diseases)
  )
}



# 组合所有图表 - 创建两行（WC和proWC）的布局
library(patchwork)

# 提取WC图表（移除图例）
wc_plots_no_legend <- lapply(disease_results, function(x) {
  x$wc_plot + theme(legend.position = "none", plot.title = element_blank())
})

# 提取proWC图表（移除图例和标题）
proWC_plots_no_legend <- lapply(disease_results, function(x) {
  x$proWC_plot + theme(legend.position = "none")
})

# 组合图表：proWC在上面，WC在下面
combined_plot <- (
  (proWC_plots_no_legend[[1]] / wc_plots_no_legend[[1]]) |  # 第一列：疾病1
    (proWC_plots_no_legend[[2]] / wc_plots_no_legend[[2]]) |  # 第二列：疾病2
    (proWC_plots_no_legend[[3]] / wc_plots_no_legend[[3]]) |  # 第三列：疾病3
    (proWC_plots_no_legend[[4]] / wc_plots_no_legend[[4]]) |   # 第四列：疾病4
    (proWC_plots_no_legend[[5]] / wc_plots_no_legend[[5]])
) + 
  plot_layout(ncol = 5) &
  theme(legend.position = "none")


# 组合图表：proWC在上面，WC在下面，增加到10个
combined_plot <- (
  (proWC_plots_no_legend[[6]] / wc_plots_no_legend[[6]]) |  # 第六列：疾病6
    (proWC_plots_no_legend[[7]] / wc_plots_no_legend[[7]]) |  # 第七列：疾病7
    (proWC_plots_no_legend[[8]] / wc_plots_no_legend[[8]]) |  # 第八列：疾病8
    (proWC_plots_no_legend[[9]] / wc_plots_no_legend[[9]]) |  # 第九列：疾病9
    (proWC_plots_no_legend[[10]] / wc_plots_no_legend[[10]])  # 第十列：疾病10
) + 
  plot_layout(ncol = 5) &  # 设置每行最多5列
  theme(legend.position = "none")  # 移除图例



# 添加共享图例
legend <- cowplot::get_legend(disease_results[[1]]$wc_plot + theme(legend.position = "bottom"))

# 最终组合
final_plot <- combined_plot / legend + 
  plot_layout(heights = c(0.9, 0.1))

# 显示最终图表
final_plot




# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot <- top_panel / combined_plot + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot






# BMI
# 对每个疾病运行分析
disease_results <- list()
significant_diseases <- list()

setwd("D:/UKB_data/ID")

analyze_disease_risk_bmi <- function(disease_code) {
  # 前面的数据准备部分保持不变
  disease_ids <- read.csv(paste0("participant_ids_", disease_code, ".csv"))
  rownames(disease_ids) <- disease_ids$Participant_ID
  match_data <- data[rownames(disease_ids), ]
  disease_date <- date[rownames(disease_ids), ]
  
  disease_results <- lapply(1:nrow(disease_date), function(i) {
    diagnoses <- disease_date$Diagnoses...ICD10[i]
    diag_list <- strsplit(diagnoses, "\\|")[[1]]
    
    if (length(diag_list) == 0 || all(diag_list == "")) return(NULL)
    
    matches <- sapply(diag_list, function(x) any(startsWith(x, disease_code), na.rm = TRUE))
    if (all(!matches)) return(NULL)
    
    positions <- which(matches)
    result <- lapply(positions, function(pos) {
      date_col <- paste0("Date.of.first.in.patient.diagnosis...ICD10...Array.", pos - 1)
      if (date_col %in% names(disease_date)) {
        date_val <- disease_date[i, date_col]
        if (!is.na(date_val) && date_val != "") as.character(date_val) else NA
      } else {
        NA
      }
    })
    
    data.frame(
      Participant_ID = rep(rownames(disease_date)[i], length(positions)),
      Disease_Date = unlist(result),
      stringsAsFactors = FALSE
    )
  })
  
  first_date <- do.call(rbind, disease_results)
  if (is.null(first_date)) return(NULL)
  first_date <- first_date %>%
    mutate(Disease_Date = as.Date(Disease_Date)) %>%
    group_by(Participant_ID) %>%
    filter(Disease_Date == min(Disease_Date)) %>%
    sample_n(1) %>%  
    ungroup()
  rownames(first_date) <- first_date$Participant_ID
  
  start_date <- data.frame(
    date_attending_assessment_centre = match_data[rownames(first_date), "date_attending_assessment_centre"],
    row.names = rownames(first_date)
  )
  merge_date <- cbind(first_date, start_date)
  
  date_filtered <- merge_date %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_disease = as.Date(Disease_Date, format = "%Y/%m/%d"),
      onset_time_days = as.numeric(difftime(date_disease, date_assessment, units = "days")),
      onset_time_years = as.numeric(difftime(date_disease, date_assessment, units = "days")) / 365.25
    ) %>%
    filter(date_disease > date_assessment) %>%
    select(-date_assessment, -date_disease)
  
  if (nrow(date_filtered) == 0) return(NULL)
  
  # 病例数据
  case_ids <- rownames(date_filtered)
  results_case <- results[case_ids, ]
  complete_case <- complete_data[case_ids, ]
  disease_df <- cbind(complete_case, results_case, onset_time_years = date_filtered[case_ids, "onset_time_years"])
  disease_df$status <- 1  # 事件发生
  
  # 对照数据
  results_control <- results[normal_ids, ]
  complete_control <- complete_data[normal_ids, ]
  
  control_assessment_dates <- data.frame(
    date_attending_assessment_centre = data[normal_ids, "date_attending_assessment_centre"]
  )
  control_df <- cbind(complete_control, results_control, control_assessment_dates)
  
  control_df <- control_df %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_death = ifelse(is.na(data$new2025516_dead_data[normal_ids]) | data$new2025516_dead_data[normal_ids] == "", 
                          NA, 
                          as.Date(data$new2025516_dead_data[normal_ids], format = "%Y/%m/%d")),
      date_censored = as.Date("2024-07-08"),
      date_end = case_when(
        !is.na(date_death) & date_death <= date_censored ~ date_death,
        TRUE ~ date_censored
      ),
      onset_time_years = as.numeric(difftime(date_end, date_assessment, units = "days")) / 365.25,
      status = 0
    ) %>%
    select(-date_attending_assessment_centre, -date_assessment, -date_death, -date_censored, -date_end)
  
  # 合并病例和对照数据
  analysis_df <- rbind(disease_df, control_df)
  
  # 按BMI分组 - 排除Underweight组
  analysis_df$BMI_group <- cut(analysis_df$BMI,
                               breaks = c(18.5, 24.9, 29.9, Inf),
                               labels = c("Normal", "Overweight", "Obese"),
                               right = TRUE, include.lowest = FALSE)
  
  # 排除BMI < 18.5的数据
  analysis_df <- analysis_df[!is.na(analysis_df$BMI_group), ]
  
  # 提取三个BMI分组的数据集（排除Underweight）
  normal_df <- analysis_df[analysis_df$BMI_group == "Normal" & !is.na(analysis_df$BMI_group), ] 
  overweight_df <- analysis_df[analysis_df$BMI_group == "Overweight" & !is.na(analysis_df$BMI_group), ]
  obese_df <- analysis_df[analysis_df$BMI_group == "Obese" & !is.na(analysis_df$BMI_group), ]
  
  # 加载必要的包
  library(survival)
  library(rms)
  library(ggplot2)
  
  # 创建生存对象函数
  create_survival_object <- function(data) {
    data$disease_surv <- with(data, Surv(onset_time_years, status))
    return(data)
  }
  
  normal_df <- create_survival_object(normal_df)
  overweight_df <- create_survival_object(overweight_df)
  obese_df <- create_survival_object(obese_df)
  
  # 设置datadist - 为每个BMI分组创建datadist对象
  dd_normal <- datadist(normal_df)
  dd_overweight <- datadist(overweight_df)
  dd_obese <- datadist(obese_df)
  
  # 创建BMI数据列表（只包含三个分组）
  bmi_data_list <- list(
    Normal = normal_df,
    Overweight = overweight_df,
    Obese = obese_df
  )
  
  # RCS分析函数 - 按BMI分层，支持WC和proWC
  rcs_analysis_disease_bmi <- function(bmi_data_list, disease_name, biomarker_type = "WC", use_log2 = FALSE) {
    
    # 定义BMI分组名称和颜色（排除Underweight）
    bmi_groups <- c("Normal", "Overweight", "Obese")
    bmi_colors <- c(
      "Normal" = "#9932CC",         # 深紫色  
      "Overweight" = "#20B2AA",     # 浅海绿色
      "Obese" = "#FF4500"           # 橙红色
    )
    
    # 存储模型和预测结果
    models <- list()
    predictions <- list()
    p_values <- list()
    
    # 为每个BMI分组拟合模型
    for (bmi_group in bmi_groups) {
      data <- bmi_data_list[[bmi_group]]
      
      # 跳过样本量太小的分组
      if (nrow(data) < 50) next
      
      # 设置datadist
      options(datadist = get(paste0("dd_", tolower(bmi_group))))
      
      # 构建公式 - 根据生物标志物类型选择
      if (biomarker_type == "WC") {
        formula <- as.formula("disease_surv ~ rcs(Actual_WC, 4) + Age + Sex")
      } else {
        formula <- as.formula("disease_surv ~ rcs(BioX_Adjusted, 4) + Age + Sex")
      }
      
      # 拟合模型
      models[[bmi_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
      
      # 计算非线性P值
      anova_result <- anova(models[[bmi_group]])
      p_values[[bmi_group]] <- ifelse(biomarker_type == "WC", 
                                      anova_result["Actual_WC", "P"], 
                                      anova_result["BioX_Adjusted", "P"])
    }
    
    # 生成预测数据 - 每个组使用自己的实际范围
    pred_combined <- data.frame()
    
    for (bmi_group in names(models)) {
      data <- bmi_data_list[[bmi_group]]
      
      # 为每个组使用自己的实际范围
      if (biomarker_type == "WC") {
        # 获取该组WC的实际范围
        wc_min <- min(data$Actual_WC, na.rm = TRUE)
        wc_max <- max(data$Actual_WC, na.rm = TRUE)
        biomarker_range <- seq(wc_min, wc_max, length = 100)
        
        # 设置datadist
        options(datadist = get(paste0("dd_", tolower(bmi_group))))
        
        pred <- Predict(models[[bmi_group]], 
                        Actual_WC = biomarker_range, 
                        Age = median(data$Age, na.rm = TRUE),
                        Sex = median(data$Sex, na.rm = TRUE),
                        ref.zero = TRUE,
                        fun = exp)
        
        pred_df <- data.frame(
          Biomarker = pred$Actual_WC, 
          HR = pred$yhat, 
          lower = pred$lower, 
          upper = pred$upper,
          BMI_Group = bmi_group,
          Biomarker_Range_Min = wc_min,
          Biomarker_Range_Max = wc_max
        )
      } else {
        # 获取该组proWC的实际范围
        prowc_min <- min(data$BioX_Adjusted, na.rm = TRUE)
        prowc_max <- max(data$BioX_Adjusted, na.rm = TRUE)
        biomarker_range <- seq(prowc_min, prowc_max, length = 100)
        
        # 设置datadist
        options(datadist = get(paste0("dd_", tolower(bmi_group))))
        
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
          BMI_Group = bmi_group,
          Biomarker_Range_Min = prowc_min,
          Biomarker_Range_Max = prowc_max
        )
      }
      
      pred_combined <- rbind(pred_combined, pred_df)
    }
    
    # 关键修改1: 将BMI_Group转换为因子并指定正确的顺序
    pred_combined$BMI_Group <- factor(pred_combined$BMI_Group, 
                                      levels = c("Normal", "Overweight", "Obese"))
    
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
    
    # 计算所有组的整体x轴范围
    x_min <- min(pred_combined$Biomarker_Range_Min, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker_Range_Max, na.rm = TRUE)
    
    # 绘制图表 - 字体增大版本
    p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = BMI_Group, fill = BMI_Group)) +
      geom_line(size = 1.2) +  # 线条适当加粗
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
      geom_hline(yintercept = ifelse(use_log2, 0, 1), 
                 linetype = "dashed", color = "black", size = 1) +
      # 关键修改2: 在scale_color_manual和scale_fill_manual中指定breaks顺序
      scale_color_manual(values = bmi_colors,
                         breaks = c("Normal", "Overweight", "Obese")) +
      scale_fill_manual(values = bmi_colors,
                        breaks = c("Normal", "Overweight", "Obese")) +
      labs(
        title = disease_name,
        x = ifelse(biomarker_type == "WC", "WC", "proWC"),
        y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio"),
        color = "BMI Group",
        fill = "BMI Group"
      ) +
      theme_minimal() +
      theme(
        text = element_text(family = "Arial"),  # 全局字体
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Arial"),  # 标题字体增大到16
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),  # 坐标轴标题增大到14
        axis.title.x = element_text(size = 14, face = "bold", family = "Arial"),  # x轴标题
        axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),  # y轴标题
        axis.text = element_text(size = 12, family = "Arial"),  # 坐标轴刻度文本增大到12
        axis.text.x = element_text(size = 12, family = "Arial"),  # x轴刻度
        axis.text.y = element_text(size = 12, family = "Arial"),  # y轴刻度
        axis.line = element_line(size = 0.8, color = "black"),  # 坐标轴线保持原样
        axis.ticks = element_line(size = 0.8, color = "black"),  # 刻度线保持原样
        axis.ticks.length = unit(0.2, "cm"),  # 刻度长度保持原样
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Arial"),  # 图例文本增大到12
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
        # 不修改整体边距，保持原样
      ) +
      # 设置x轴范围，包含所有组的范围
      coord_cartesian(xlim = c(x_min, x_max))
    
    # 在图上添加P值 - 保持原样
    y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
    i <- 1
    for (bmi_group in names(p_values)) {
      p_value <- p_values[[bmi_group]]
      
      # 获取该组的数据用于标注位置
      group_data <- pred_combined[pred_combined$BMI_Group == bmi_group, ]
      group_x_min <- min(group_data$Biomarker, na.rm = TRUE)
      group_x_max <- max(group_data$Biomarker, na.rm = TRUE)
      
      # 修正P值显示逻辑
      p_value_label <- if (p_value < 0.001) {
        "≤ 0.001"
      } else if (p_value < 0.01) {
        sprintf("= %.3f", p_value)
      } else if (p_value < 0.05) {
        sprintf("= %.3f", p_value)
      } else {
        sprintf("= %.3f", p_value)
      }
      
      p <- p + annotate("text", 
                        x = group_x_min + 0.1 * (group_x_max - group_x_min), 
                        y = y_max - y_positions[i] * y_range,
                        label = paste("P-nonlinear", p_value_label),
                        size = 3.2, fontface = "bold", hjust = 0,  # P值标注保持原样
                        color = bmi_colors[bmi_group],
                        family = "Arial")
      
      i <- i + 1
    }
    
    return(list(plot = p, 
                predictions = pred_combined,
                p_values = p_values,
                models = models,
                biomarker_type = biomarker_type,
                biomarker_range_min = x_min,
                biomarker_range_max = x_max))
  }
  
  # 运行RCS分析 - 分别分析WC和proWC
  disease_name <- disease_code
  
  # 分析WC
  rcs_result_wc <- rcs_analysis_disease_bmi(bmi_data_list, disease_name, "WC", use_log2 = TRUE)
  
  # 分析proWC
  rcs_result_prowc <- rcs_analysis_disease_bmi(bmi_data_list, disease_name, "proWC", use_log2 = TRUE)
  
  # 返回所有结果
  return(list(
    wc_plot = rcs_result_wc$plot,
    wc_predictions = rcs_result_wc$predictions,
    wc_p_values = rcs_result_wc$p_values,
    wc_models = rcs_result_wc$models,
    
    proWC_plot = rcs_result_prowc$plot,
    proWC_predictions = rcs_result_prowc$predictions,
    proWC_p_values = rcs_result_prowc$p_values,
    proWC_models = rcs_result_prowc$models,
    
    disease_name = disease_name,
    sample_info = data.frame(
      n_cases = nrow(disease_df),
      n_controls = nrow(control_df),
      total_n = nrow(analysis_df),
      n_normal = nrow(normal_df),
      n_overweight = nrow(overweight_df),
      n_obese = nrow(obese_df)
    )
  ))
}


for (disease_code in disease_codes) {
  cat("Analyzing disease:", disease_code, "\n")
  
  tryCatch({
    result <- analyze_disease_risk_bmi(disease_code)
    
    if (!is.null(result)) {
      disease_results[[disease_code]] <- result
      
      # 检查是否有错误信息
      if (!is.null(result$error)) {
        cat("  -> Model error:", result$error, "\n")
        failed_diseases[[disease_code]] <- result$error
        next
      }
      
      wc_p_vals <- unlist(result$wc_p_values)
      proWC_p_vals <- unlist(result$proWC_p_values)
      
      # 只考虑成功拟合的组
      valid_wc_p <- wc_p_vals[!is.na(wc_p_vals)]
      valid_proWC_p <- proWC_p_vals[!is.na(proWC_p_vals)]
      
      if (length(valid_wc_p) > 0 && length(valid_proWC_p) > 0 &&
          all(valid_wc_p < 0.05) && all(valid_proWC_p < 0.05)) {
        significant_diseases[[disease_code]] <- result
        cat("  -> Significant disease found!\n")
      }
    }
  }, error = function(e) {
    cat("  -> Fatal error:", e$message, "\n")
    failed_diseases[[disease_code]] <- e$message
  })
}

# 保存失败信息
if (length(failed_diseases) > 0) {
  failed_df <- data.frame(
    Disease_Code = names(failed_diseases),
    Error_Message = unlist(failed_diseases)
  )
}


# 组合所有图表 - 创建两行（WC和proWC）的布局
library(patchwork)

# 提取WC图表（移除图例）
wc_plots_no_legend <- lapply(disease_results, function(x) {
  x$wc_plot + theme(legend.position = "none", plot.title = element_blank())
})

# 提取proWC图表（移除图例和标题）
proWC_plots_no_legend <- lapply(disease_results, function(x) {
  x$proWC_plot + theme(legend.position = "none")
})

# 组合图表：proWC在上面，WC在下面
combined_plot <- (
  (proWC_plots_no_legend[[1]] / wc_plots_no_legend[[1]]) |  # 第一列：疾病1
    (proWC_plots_no_legend[[2]] / wc_plots_no_legend[[2]]) |  # 第二列：疾病2
    (proWC_plots_no_legend[[3]] / wc_plots_no_legend[[3]]) |  # 第三列：疾病3
    (proWC_plots_no_legend[[4]] / wc_plots_no_legend[[4]]) |  # 第四列：疾病4
    (proWC_plots_no_legend[[5]] / wc_plots_no_legend[[5]])
) + 
  plot_layout(ncol = 5) &
  theme(legend.position = "none")

# 组合图表：proWC在上面，WC在下面，增加到10个
combined_plot <- (
  (proWC_plots_no_legend[[6]] / wc_plots_no_legend[[6]]) |  # 第六列：疾病6
    (proWC_plots_no_legend[[7]] / wc_plots_no_legend[[7]]) |  # 第七列：疾病7
    (proWC_plots_no_legend[[8]] / wc_plots_no_legend[[8]]) |  # 第八列：疾病8
    (proWC_plots_no_legend[[9]] / wc_plots_no_legend[[9]]) |  # 第九列：疾病9
    (proWC_plots_no_legend[[10]] / wc_plots_no_legend[[10]])  # 第十列：疾病10
) + 
  plot_layout(ncol = 5) &  # 设置每行最多5列
  theme(legend.position = "none")  # 移除图例

# 添加共享图例
legend <- cowplot::get_legend(disease_results[[1]]$wc_plot + theme(legend.position = "bottom"))

# 最终组合
final_plot <- combined_plot / legend + 
  plot_layout(heights = c(0.9, 0.1))

# 显示最终图表
final_plot




# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot2 <- top_panel / combined_plot + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot2






# Age
# 对每个疾病运行分析
disease_results <- list()

setwd("D:/UKB_data/ID")

analyze_disease_risk_age <- function(disease_code) {
  # 前面的数据准备部分保持不变
  disease_ids <- read.csv(paste0("participant_ids_", disease_code, ".csv"))
  rownames(disease_ids) <- disease_ids$Participant_ID
  match_data <- data[rownames(disease_ids), ]
  disease_date <- date[rownames(disease_ids), ]
  
  disease_results <- lapply(1:nrow(disease_date), function(i) {
    diagnoses <- disease_date$Diagnoses...ICD10[i]
    diag_list <- strsplit(diagnoses, "\\|")[[1]]
    
    if (length(diag_list) == 0 || all(diag_list == "")) return(NULL)
    
    matches <- sapply(diag_list, function(x) any(startsWith(x, disease_code), na.rm = TRUE))
    if (all(!matches)) return(NULL)
    
    positions <- which(matches)
    result <- lapply(positions, function(pos) {
      date_col <- paste0("Date.of.first.in.patient.diagnosis...ICD10...Array.", pos - 1)
      if (date_col %in% names(disease_date)) {
        date_val <- disease_date[i, date_col]
        if (!is.na(date_val) && date_val != "") as.character(date_val) else NA
      } else {
        NA
      }
    })
    
    data.frame(
      Participant_ID = rep(rownames(disease_date)[i], length(positions)),
      Disease_Date = unlist(result),
      stringsAsFactors = FALSE
    )
  })
  
  first_date <- do.call(rbind, disease_results)
  if (is.null(first_date)) return(NULL)
  first_date <- first_date %>%
    mutate(Disease_Date = as.Date(Disease_Date)) %>%
    group_by(Participant_ID) %>%
    filter(Disease_Date == min(Disease_Date)) %>%
    sample_n(1) %>%  
    ungroup()
  rownames(first_date) <- first_date$Participant_ID
  
  start_date <- data.frame(
    date_attending_assessment_centre = match_data[rownames(first_date), "date_attending_assessment_centre"],
    row.names = rownames(first_date)
  )
  merge_date <- cbind(first_date, start_date)
  
  date_filtered <- merge_date %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_disease = as.Date(Disease_Date, format = "%Y/%m/%d"),
      onset_time_days = as.numeric(difftime(date_disease, date_assessment, units = "days")),
      onset_time_years = as.numeric(difftime(date_disease, date_assessment, units = "days")) / 365.25
    ) %>%
    filter(date_disease > date_assessment) %>%
    select(-date_assessment, -date_disease)
  
  if (nrow(date_filtered) == 0) return(NULL)
  
  # 病例数据
  case_ids <- rownames(date_filtered)
  results_case <- results[case_ids, ]
  complete_case <- complete_data[case_ids, ]
  disease_df <- cbind(complete_case, results_case, onset_time_years = date_filtered[case_ids, "onset_time_years"])
  disease_df$status <- 1  # 事件发生
  
  # 对照数据
  results_control <- results[normal_ids, ]
  complete_control <- complete_data[normal_ids, ]
  
  control_assessment_dates <- data.frame(
    date_attending_assessment_centre = data[normal_ids, "date_attending_assessment_centre"]
  )
  control_df <- cbind(complete_control, results_control, control_assessment_dates)
  
  control_df <- control_df %>%
    mutate(
      date_assessment = as.Date(date_attending_assessment_centre, format = "%Y/%m/%d"),
      date_death = ifelse(is.na(data$new2025516_dead_data[normal_ids]) | data$new2025516_dead_data[normal_ids] == "", 
                          NA, 
                          as.Date(data$new2025516_dead_data[normal_ids], format = "%Y/%m/%d")),
      date_censored = as.Date("2024-07-08"),
      date_end = case_when(
        !is.na(date_death) & date_death <= date_censored ~ date_death,
        TRUE ~ date_censored
      ),
      onset_time_years = as.numeric(difftime(date_end, date_assessment, units = "days")) / 365.25,
      status = 0
    ) %>%
    select(-date_attending_assessment_centre, -date_assessment, -date_death, -date_censored, -date_end)
  
  # 合并病例和对照数据
  analysis_df <- rbind(disease_df, control_df)
  
  # 按年龄分组（使用中位数分割）
  age_median <- median(analysis_df$Age, na.rm = TRUE)
  analysis_df$age_group <- ifelse(analysis_df$Age >= age_median, "Higher Age", "Lower Age")
  
  higher_age_data <- subset(analysis_df, age_group == "Higher Age")
  lower_age_data <- subset(analysis_df, age_group == "Lower Age")
  
  # 加载必要的包
  library(survival)
  library(rms)
  library(ggplot2)
  
  # 创建生存对象函数
  create_survival_object <- function(data) {
    data$disease_surv <- with(data, Surv(onset_time_years, status))
    return(data)
  }
  
  higher_age_data <- create_survival_object(higher_age_data)
  lower_age_data <- create_survival_object(lower_age_data)
  
  # 设置datadist - 为每个年龄分组创建datadist对象
  dd_higher_age <- datadist(higher_age_data)
  dd_lower_age <- datadist(lower_age_data)
  
  # 创建年龄数据列表
  age_data_list <- list(
    "Higher Age" = higher_age_data,
    "Lower Age" = lower_age_data
  )
  
  # RCS分析函数 - 按年龄分层，支持WC和proWC
  rcs_analysis_disease_age <- function(age_data_list, disease_name, biomarker_type = "WC", use_log2 = FALSE) {
    
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
    
    # 为每个年龄分组拟合模型
    for (age_group in age_groups) {
      data <- age_data_list[[age_group]]
      
      # 设置datadist
      options(datadist = get(paste0("dd_", tolower(gsub(" ", "_", age_group)))))
      
      # 构建公式 - 根据生物标志物类型选择
      if (biomarker_type == "WC") {
        formula <- as.formula("disease_surv ~ rcs(Actual_WC, 4) + Sex + BMI")
      } else {
        formula <- as.formula("disease_surv ~ rcs(BioX_Adjusted, 4) + Sex + BMI")
      }
      
      # 拟合模型
      models[[age_group]] <- cph(formula, data = data, x = TRUE, y = TRUE)
      
      # 计算非线性P值
      anova_result <- anova(models[[age_group]])
      p_values[[age_group]] <- ifelse(biomarker_type == "WC", 
                                      anova_result["Actual_WC", "P"], 
                                      anova_result["BioX_Adjusted", "P"])
    }
    
    # 生成预测数据 - 每个组使用自己的实际范围
    pred_combined <- data.frame()
    
    for (age_group in age_groups) {
      data <- age_data_list[[age_group]]
      
      # 为每个组使用自己的实际范围
      if (biomarker_type == "WC") {
        # 获取该组WC的实际范围
        wc_min <- min(data$Actual_WC, na.rm = TRUE)
        wc_max <- max(data$Actual_WC, na.rm = TRUE)
        biomarker_range <- seq(wc_min, wc_max, length = 100)
        
        # 设置datadist
        options(datadist = get(paste0("dd_", tolower(gsub(" ", "_", age_group)))))
        
        pred <- Predict(models[[age_group]], 
                        Actual_WC = biomarker_range, 
                        Sex = median(data$Sex, na.rm = TRUE),  # 使用中位数性别
                        BMI = median(data$BMI, na.rm = TRUE),
                        ref.zero = TRUE,
                        fun = exp)
        
        pred_df <- data.frame(
          Biomarker = pred$Actual_WC, 
          HR = pred$yhat, 
          lower = pred$lower, 
          upper = pred$upper,
          Age_Group = age_group,
          Biomarker_Range_Min = wc_min,
          Biomarker_Range_Max = wc_max
        )
      } else {
        # 获取该组proWC的实际范围
        prowc_min <- min(data$BioX_Adjusted, na.rm = TRUE)
        prowc_max <- max(data$BioX_Adjusted, na.rm = TRUE)
        biomarker_range <- seq(prowc_min, prowc_max, length = 100)
        
        # 设置datadist
        options(datadist = get(paste0("dd_", tolower(gsub(" ", "_", age_group)))))
        
        pred <- Predict(models[[age_group]], 
                        BioX_Adjusted = biomarker_range, 
                        Sex = median(data$Sex, na.rm = TRUE),  # 使用中位数性别
                        BMI = median(data$BMI, na.rm = TRUE),
                        ref.zero = TRUE,
                        fun = exp)
        
        pred_df <- data.frame(
          Biomarker = pred$BioX_Adjusted, 
          HR = pred$yhat, 
          lower = pred$lower, 
          upper = pred$upper,
          Age_Group = age_group,
          Biomarker_Range_Min = prowc_min,
          Biomarker_Range_Max = prowc_max
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
    
    # 计算所有组的整体x轴范围
    x_min <- min(pred_combined$Biomarker_Range_Min, na.rm = TRUE)
    x_max <- max(pred_combined$Biomarker_Range_Max, na.rm = TRUE)
    
    # 绘制图表 - 字体增大版本
    p <- ggplot(pred_combined, aes(x = Biomarker, y = HR, color = Age_Group, fill = Age_Group)) +
      geom_line(size = 1.2) +  # 线条适当加粗
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
      geom_hline(yintercept = ifelse(use_log2, 0, 1), 
                 linetype = "dashed", color = "black", size = 1) +
      scale_color_manual(values = age_colors) +
      scale_fill_manual(values = age_colors) +
      labs(
        title = disease_name,
        x = ifelse(biomarker_type == "WC", "WC", "proWC"),
        y = ifelse(use_log2, expression(bold(log[bold("2")] * "(Hazard Ratio)")), "Hazard Ratio"),
        color = "Age Group",
        fill = "Age Group"
      ) +
      theme_minimal() +
      theme(
        text = element_text(family = "Arial"),  # 全局字体
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Arial"),  # 标题字体增大到16
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),  # 坐标轴标题增大到14
        axis.title.x = element_text(size = 14, face = "bold", family = "Arial"),  # x轴标题
        axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),  # y轴标题
        axis.text = element_text(size = 12, family = "Arial"),  # 坐标轴刻度文本增大到12
        axis.text.x = element_text(size = 12, family = "Arial"),  # x轴刻度
        axis.text.y = element_text(size = 12, family = "Arial"),  # y轴刻度
        axis.line = element_line(size = 0.8, color = "black"),  # 坐标轴线保持原样
        axis.ticks = element_line(size = 0.8, color = "black"),  # 刻度线保持原样
        axis.ticks.length = unit(0.2, "cm"),  # 刻度长度保持原样
        legend.title = element_blank(),  # 图例标题设为空（与其他函数一致）
        legend.text = element_text(size = 12, face = "bold", family = "Arial"),  # 图例文本增大到12
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
        # 不修改整体边距，保持原样
      ) +
      # 设置x轴范围，包含所有组的范围
      coord_cartesian(xlim = c(x_min, x_max))
    
    # 在图上添加P值 - 保持原样
    y_positions <- seq(0.04, by = 0.08, length.out = length(p_values))
    i <- 1
    for (age_group in age_groups) {
      p_value <- p_values[[age_group]]
      
      # 获取该组的数据用于标注位置
      group_data <- pred_combined[pred_combined$Age_Group == age_group, ]
      group_x_min <- min(group_data$Biomarker, na.rm = TRUE)
      group_x_max <- max(group_data$Biomarker, na.rm = TRUE)
      
      # 修正P值显示逻辑
      p_value_label <- if (p_value < 0.001) {
        "≤ 0.001"
      } else if (p_value < 0.01) {
        sprintf("= %.3f", p_value)
      } else if (p_value < 0.05) {
        sprintf("= %.3f", p_value)
      } else {
        sprintf("= %.3f", p_value)
      }
      
      p <- p + annotate("text", 
                        x = group_x_min + 0.1 * (group_x_max - group_x_min), 
                        y = y_max - y_positions[i] * y_range,
                        label = paste("P-nonlinear", p_value_label),
                        size = 3.2, fontface = "bold", hjust = 0,  # P值标注保持原样
                        color = age_colors[age_group],
                        family = "Arial")
      
      i <- i + 1
    }
    
    return(list(plot = p, 
                predictions = pred_combined,
                p_values = p_values,
                models = models,
                biomarker_type = biomarker_type,
                biomarker_range_min = x_min,
                biomarker_range_max = x_max))
  }
  
  # 运行RCS分析 - 分别分析WC和proWC
  disease_name <- disease_code
  
  # 分析WC
  rcs_result_wc <- rcs_analysis_disease_age(age_data_list, disease_name, "WC", use_log2 = TRUE)
  
  # 分析proWC
  rcs_result_prowc <- rcs_analysis_disease_age(age_data_list, disease_name, "proWC", use_log2 = TRUE)
  
  # 返回所有结果
  return(list(
    wc_plot = rcs_result_wc$plot,
    wc_predictions = rcs_result_wc$predictions,
    wc_p_values = rcs_result_wc$p_values,
    wc_models = rcs_result_wc$models,
    
    proWC_plot = rcs_result_prowc$plot,
    proWC_predictions = rcs_result_prowc$predictions,
    proWC_p_values = rcs_result_prowc$p_values,
    proWC_models = rcs_result_prowc$models,
    
    disease_name = disease_name,
    sample_info = data.frame(
      n_cases = nrow(disease_df),
      n_controls = nrow(control_df),
      total_n = nrow(analysis_df),
      n_higher_age = nrow(higher_age_data),
      n_lower_age = nrow(lower_age_data),
      age_median = age_median
    )
  ))
}


significant_diseases <- list()

for (disease_code in disease_codes) {
  cat("Analyzing disease:", disease_code, "\n")
  result <- analyze_disease_risk_age(disease_code)
  if (!is.null(result)) {
    disease_results[[disease_code]] <- result
    
    # 检查是否满足显著性条件
    wc_p_vals <- unlist(result$wc_p_values)
    proWC_p_vals <- unlist(result$proWC_p_values)
    
    # 条件：两个年龄组在WC和proWC上的非线性P值都 < 0.05
    if (all(wc_p_vals < 0.05) && all(proWC_p_vals < 0.05)) {
      significant_diseases[[disease_code]] <- result
      cat("  -> Significant disease found!\n")
    }
  }
}




# 组合所有图表 - 创建两行（WC和proWC）的布局
library(patchwork)

# 提取WC图表（移除图例）
wc_plots_no_legend <- lapply(disease_results, function(x) {
  x$wc_plot + theme(legend.position = "none", plot.title = element_blank())
})

# 提取proWC图表（移除图例和标题）
proWC_plots_no_legend <- lapply(disease_results, function(x) {
  x$proWC_plot + theme(legend.position = "none")
})

# 组合图表：proWC在上面，WC在下面
combined_plot <- (
  (proWC_plots_no_legend[[1]] / wc_plots_no_legend[[1]]) |  # 第一列：疾病1
    (proWC_plots_no_legend[[2]] / wc_plots_no_legend[[2]]) |  # 第二列：疾病2
    (proWC_plots_no_legend[[3]] / wc_plots_no_legend[[3]]) |  # 第三列：疾病3
    (proWC_plots_no_legend[[4]] / wc_plots_no_legend[[4]]) |  # 第四列：疾病4
    (proWC_plots_no_legend[[5]] / wc_plots_no_legend[[5]])
) + 
  plot_layout(ncol = 5) &
  theme(legend.position = "none")

# 组合图表：proWC在上面，WC在下面，增加到10个
combined_plot <- (
  (proWC_plots_no_legend[[6]] / wc_plots_no_legend[[6]]) |  # 第六列：疾病6
    (proWC_plots_no_legend[[7]] / wc_plots_no_legend[[7]]) |  # 第七列：疾病7
    (proWC_plots_no_legend[[8]] / wc_plots_no_legend[[8]]) |  # 第八列：疾病8
    (proWC_plots_no_legend[[9]] / wc_plots_no_legend[[9]]) |  # 第九列：疾病9
    (proWC_plots_no_legend[[10]] / wc_plots_no_legend[[10]])  # 第十列：疾病10
) + 
  plot_layout(ncol = 5) &  # 设置每行最多5列
  theme(legend.position = "none")  # 移除图例

# 添加共享图例
legend <- cowplot::get_legend(disease_results[[1]]$wc_plot + theme(legend.position = "bottom"))

# 最终组合
final_plot <- combined_plot / legend + 
  plot_layout(heights = c(0.9, 0.1))

# 显示最终图表
final_plot



# 创建一个包含图例的顶部区域
top_panel <- plot_spacer() + 
  inset_element(legend, 
                left = 0, right = 1, 
                top = 1, bottom = 0) +
  plot_layout(heights = c(0.1))

# 组合图表
final_plot3 <- top_panel / combined_plot + 
  plot_layout(heights = c(0.1, 0.9)) &
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# 显示图表
final_plot3


