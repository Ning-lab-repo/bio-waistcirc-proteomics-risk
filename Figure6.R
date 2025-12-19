# 组合
complete_data_WC <- read.csv("WC_death_time.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df_WC <- cbind(complete_data_WC, results)

# 加载必要的包
library(survival)
library(tidyverse)
library(broom)

# 准备数据
analysis_data <- temp_df_WC[, c("WC", "BioX_Adjusted", "Age", "Sex", 
                                "all_cause_death", "cvd_death", "cancer_death", "diabetes_death",
                                "survival_time")] %>%
  rename(proWC = BioX_Adjusted) %>%
  na.omit()

# 计算中位数
wc_median <- median(analysis_data$WC, na.rm = TRUE)
prowc_median <- median(analysis_data$proWC, na.rm = TRUE)

cat("WC中位数:", wc_median, "\n")
cat("proWC中位数:", prowc_median, "\n")

# 创建分组变量
analysis_data <- analysis_data %>%
  mutate(
    WC_group = ifelse(WC <= wc_median, "Lower WC", "Higher WC"),
    proWC_group = ifelse(proWC <= prowc_median, "Lower proWC", "Higher proWC"),
    # 创建联合分组
    joint_group = case_when(
      WC_group == "Lower WC" & proWC_group == "Lower proWC" ~ "Lower WC + Lower proWC",
      WC_group == "Lower WC" & proWC_group == "Higher proWC" ~ "Lower WC + Higher proWC",
      WC_group == "Higher WC" & proWC_group == "Lower proWC" ~ "Higher WC + Lower proWC",
      WC_group == "Higher WC" & proWC_group == "Higher proWC" ~ "Higher WC + Higher proWC"
    ),
    # 将联合分组设置为因子，以Lower WC + Lower proWC为参考组（修改这里）
    joint_group = factor(joint_group, 
                         levels = c("Lower WC + Lower proWC",  # 参考组放在第一位
                                    "Lower WC + Higher proWC",
                                    "Higher WC + Lower proWC", 
                                    "Higher WC + Higher proWC"))
  )

# 检查各组人数
cat("\n各组人数分布:\n")
print(table(analysis_data$joint_group))

# 定义要分析的结局变量
outcomes <- c("all_cause_death", "cvd_death", "cancer_death", "diabetes_death")

# 存储结果的列表
results_list <- list()


for (outcome in outcomes) {
  
  # 创建生存对象
  if (outcome == "all_cause_death") {
    survival_obj <- with(analysis_data, Surv(survival_time, all_cause_death))
    outcome_name <- "All-cause mortality"
    death_var <- "all_cause_death"
  } else if (outcome == "cvd_death") {
    survival_obj <- with(analysis_data, Surv(survival_time, cvd_death))
    outcome_name <- "CVD mortality"
    death_var <- "cvd_death"
  } else if (outcome == "cancer_death") {
    survival_obj <- with(analysis_data, Surv(survival_time, cancer_death))
    outcome_name <- "Cancer mortality"
    death_var <- "cancer_death"
  } else if (outcome == "diabetes_death") {
    survival_obj <- with(analysis_data, Surv(survival_time, diabetes_death))
    outcome_name <- "Diabetes mortality"
    death_var <- "diabetes_death"
  }
  
  # 运行Cox模型（调整年龄和性别）
  cox_model <- coxph(survival_obj ~ joint_group + Age + Sex, data = analysis_data)
  
  # 提取模型结果
  model_summary <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
  
  # 计算各组的死亡人数
  death_counts <- analysis_data %>%
    group_by(joint_group) %>%
    summarise(
      n_deaths = sum(.data[[death_var]], na.rm = TRUE)
    )
  
  # 提取分组变量的结果
  group_results <- model_summary %>%
    filter(str_detect(term, "joint_group")) %>%
    mutate(
      term = str_remove(term, "joint_group"),
      hr = estimate,
      ci_lower = conf.low,
      ci_upper = conf.high,
      hr_ci = sprintf("%.2f (%.2f–%.2f)", hr, ci_lower, ci_upper),
      n = as.numeric(table(analysis_data$joint_group)[str_remove(term, "joint_group")]),
      outcome = outcome_name,
      group = term
    ) %>%
    dplyr::select(outcome, group, n, hr_ci, hr, ci_lower, ci_upper, p.value)
  
  # 为每个结果行添加死亡人数
  group_results <- group_results %>%
    left_join(death_counts, by = c("group" = "joint_group"))
  
  # 添加参考组（Lower WC + Lower proWC）
  reference_deaths <- death_counts %>%
    filter(joint_group == "Lower WC + Lower proWC") %>%
    pull(n_deaths)
  
  reference_group <- data.frame(
    outcome = outcome_name,
    group = "Lower WC + Lower proWC",
    n = as.numeric(table(analysis_data$joint_group)["Lower WC + Lower proWC"]),
    hr_ci = "1.00 (reference)",
    hr = 1,
    ci_lower = NA,
    ci_upper = NA,
    p.value = NA,
    n_deaths = reference_deaths
  )
  
  # 合并结果
  outcome_results <- bind_rows(reference_group, group_results)
  results_list[[outcome]] <- outcome_results
}

# 合并所有结果
final_results <- bind_rows(results_list)

# 准备森林图数据
forest_data <- final_results %>%
  mutate(
    # 设置分组顺序，参考组在底部
    group = factor(group, 
                   levels = c("Higher WC + Higher proWC",
                              "Lower WC + Higher proWC",
                              "Higher WC + Lower proWC",
                              "Lower WC + Lower proWC"),
                   labels = c("Higher WC + Higher proWC",
                              "Lower WC + Higher proWC",
                              "Higher WC + Lower proWC",
                              "Lower WC + Lower proWC")),
    outcome = factor(outcome, levels = c("All-cause mortality", "CVD mortality", 
                                         "Cancer mortality", "Diabetes mortality")),
    # 添加显著性标记（参考组为空）
    significance = case_when(
      group == "Lower WC + Lower proWC" ~ "",  # 参考组不显示显著性
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # 为参考组设置特殊的HR显示
    hr_display = ifelse(group == "Lower WC + Lower proWC", "1.00 (reference)", hr_ci),
    # 添加形状变量：显著的使用三角形(24)，不显著的使用圆形(21)，参考组使用圆形(21)
    point_shape = case_when(
      group == "Lower WC + Lower proWC" ~ 21,  # 参考组用圆形
      p.value < 0.05 ~ 24,  # 显著的用三角形
      TRUE ~ 21  # 不显著的用圆形
    )
  )

# 为四个组定义不同颜色
group_colors <- c(
  "Higher WC + Higher proWC" = "#E41A1C",    # 红色
  "Lower WC + Higher proWC" = "#4DAF4A",     # 绿色
  "Higher WC + Lower proWC" = "#377EB8",     # 蓝色
  "Lower WC + Lower proWC" = "#FF7F00"       # 橙色
)

# 创建显著性图例数据
significance_legend_data <- data.frame(
  shape = c(24, 21),
  label = c("P < 0.05", "P ≥ 0.05")
)

# 创建森林图（使用log2转换的HR值）
forest_plot <- ggplot(forest_data, aes(x = log2(hr), y = group)) +  # 修改这里：使用log2(hr)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +  # 修改这里：xintercept = 0
  
  # 为参考组以外的组添加误差线
  geom_errorbarh(aes(xmin = log2(ci_lower), xmax = log2(ci_upper)),  # 修改这里：对置信区间也进行log2转换
                 height = 0.2, color = "black", linewidth = 0.8,
                 data = filter(forest_data, group != "Lower WC + Lower proWC")) +
  
  # 使用形状21（圆形）和24（三角形），这些形状可以显示填充色
  geom_point(aes(fill = group, shape = factor(point_shape)), 
             size = 3, stroke = 0.8, color = "black") +
  
  # 分面显示不同结局
  facet_grid(outcome ~ ., scales = "free_y", space = "free_y") +
  
  # 设置x轴
  scale_x_continuous(
    name = expression(bold(log[bold("2")] * "(Hazard Ratio)")),  # 使用数学表达式
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  scale_y_discrete(name = NULL) +
  # 设置填充颜色
  scale_fill_manual(values = group_colors) +
  # 设置形状映射：21=圆形，24=三角形
  scale_shape_manual(
    name = NULL,
    values = c("21" = 21, "24" = 24),
    labels = c("Not Significant", "Significant")
  ) +
  
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.title.x = element_text(face = "bold", size = 11, family = "Arial", 
                                margin = margin(t = 5)),  # 增加x轴标题大小
    axis.title.y = element_text(face = "bold", size = 10, family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial"),
    axis.text.x = element_text(size = 8, family = "Arial"),
    axis.text.y = element_text(size = 8, family = "Arial"),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 9, family = "Arial"),
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0.2, "cm"),
    legend.text = element_text(size = 8, family = "Arial"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(0, 1, 0.3, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  ) +
  
  # 设置图例样式
  guides(
    fill = "none",
    shape = guide_legend(
      title = NULL,
      override.aes = list(
        fill = "black",
        color = "black",
        size = 2.5,
        stroke = 0.6
      ),
      nrow = 1,
      direction = "horizontal"
    )
  )

# 显示图形
print(forest_plot)


# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(flextable)
setwd("D:/UKB_data/NC")

# 读取数据
or_data <- read.csv("all_diseases_or_results_combine.csv", 
                    header = TRUE,
                    stringsAsFactors = FALSE)

class_data <- read.csv("class.csv", header = TRUE, stringsAsFactors = FALSE)

# 创建疾病代码到类别的映射
disease_to_category <- list()

# 构建映射关系
for(i in 1:nrow(class_data)){
  code_range <- class_data[i, "ICD.10.code"]
  category <- class_data[i, "Disease"]
  
  # 分割代码范围（例如"A01-B89"）
  code_parts <- strsplit(code_range, "-")[[1]]
  
  if(length(code_parts) == 2){
    start_code <- code_parts[1]
    end_code <- code_parts[2]
    
    # 提取字母和数字部分
    start_letter <- gsub("[0-9].*", "", start_code)
    start_num <- as.numeric(gsub("[^0-9]", "", start_code))
    end_letter <- gsub("[0-9].*", "", end_code)
    end_num <- as.numeric(gsub("[^0-9]", "", end_code))
    
    # 生成连续的代码
    if(start_letter == end_letter){
      # 相同字母的情况（如A01-A99）
      for(num in start_num:end_num){
        code_prefix <- paste0(start_letter, sprintf("%02d", num))
        disease_to_category[[code_prefix]] <- category
      }
    } else {
      # 不同字母的情况（如A01-B89）
      # 生成从起始字母到结束字母的所有组合
      all_letters <- LETTERS[which(LETTERS == start_letter):which(LETTERS == end_letter)]
      
      for(letter in all_letters){
        if(letter == start_letter){
          # 起始字母：从start_num到99
          for(num in start_num:99){
            code_prefix <- paste0(letter, sprintf("%02d", num))
            disease_to_category[[code_prefix]] <- category
          }
        } else if(letter == end_letter){
          # 结束字母：从00到end_num
          for(num in 0:end_num){
            code_prefix <- paste0(letter, sprintf("%02d", num))
            disease_to_category[[code_prefix]] <- category
          }
        } else {
          # 中间字母：所有00-99
          for(num in 0:99){
            code_prefix <- paste0(letter, sprintf("%02d", num))
            disease_to_category[[code_prefix]] <- category
          }
        }
      }
    }
  } else {
    # 单个代码
    disease_to_category[[code_range]] <- category
  }
}

all_code <- read.csv("D:\\zip\\赖煜萱8.25\\患病人数统计表_基线后首次患病.csv", fileEncoding = "GBK")
high_codes <- all_code[all_code$Number.of.participants > 100, ]

# 将high_codes中的疾病代码分类
categorized_diseases <- list()

for(disease_code in high_codes$Disease_Code){
  # 提取疾病代码的前三位（去掉小数点后的部分）
  base_code <- gsub("\\..*", "", disease_code)  # 移除小数点及后面的内容
  
  # 找到对应的疾病类别
  category <- NA
  for(code_pattern in names(disease_to_category)){
    if(grepl(paste0("^", code_pattern), base_code)){
      category <- disease_to_category[[code_pattern]]
      break
    }
  }
  
  if(!is.na(category)){
    if(is.null(categorized_diseases[[category]])){
      categorized_diseases[[category]] <- character()
    }
    categorized_diseases[[category]] <- c(categorized_diseases[[category]], disease_code)
  } else {
    cat("未找到疾病代码", disease_code, "的类别\n")
  }
}

# 假设已经有了 categorized_diseases 列表
# 将疾病代码分类到对应的类别中
for(category_name in names(categorized_diseases)) {
  or_data$category[or_data$Disease_Code %in% categorized_diseases[[category_name]]] <- category_name
}

# 移除未分类的数据
or_data_classified <- or_data %>% filter(!is.na(category))

major_data <- or_data_classified %>%
  filter(!grepl("\\.", Disease_Code) & Disease_Code != "I85")

# 准备绘图数据
plot_data <- major_data %>%
  # 定义组的顺序和颜色（不换行）
  mutate(Group = factor(Group, 
                        levels = c("Lower WC + Lower proWC", 
                                   "Lower WC + Higher proWC",
                                   "Higher WC + Lower proWC", 
                                   "Higher WC + Higher proWC"),
                        labels = c("Lower WC + Lower proWC", 
                                   "Lower WC + Higher proWC",
                                   "Higher WC + Lower proWC", 
                                   "Higher WC + Higher proWC"),
                        ordered = TRUE)) %>%
  # 计算log2 OR用于绘图
  mutate(log2_or = log2(estimate)) %>%
  # 根据p.adj判断显著性，如果NA则为不显著
  mutate(significance = ifelse(is.na(p.adj), "Not significant",
                               ifelse(p.adj < 0.05, "Significant", "Not significant")),
         point_shape = ifelse(is.na(p.adj), 16,
                              ifelse(p.adj < 0.05, 17, 16)))  # 17=三角形, 16=圆形

# 找出所有数据中log2 OR最高的前5个疾病（正值）- 不分组
top_diseases <- plot_data %>%
  filter(log2_or > 0) %>%  # 只考虑正值的OR
  arrange(desc(log2_or)) %>%
  slice_head(n = 20) %>%
  # 只使用疾病代码作为标签
  mutate(label = Disease_Code)

# 定义颜色方案（不换行标签）
group_colors <- c(
  "Higher WC + Higher proWC" = "#E41A1C",    # 红色
  "Lower WC + Higher proWC" = "#4DAF4A",     # 绿色
  "Higher WC + Lower proWC" = "#377EB8",     # 蓝色
  "Lower WC + Lower proWC" = "#FF7F00"       # 橙色
)

# 创建可视化图表（使用log2 OR）
p1 <- ggplot(plot_data, aes(x = category, y = log2_or, 
                            color = Group, shape = significance,
                            group = interaction(Group, Disease_Code))) +  # 关键行
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  
  # 为所有数据中前5个疾病添加标签（只标注疾病代码）
  geom_text(data = top_diseases,
            aes(x = category, y = log2_or, label = label),
            size = 2.5, color = "black", nudge_x = 0.25, nudge_y = 0.02,
            hjust = 0, show.legend = FALSE, fontface = "bold") +
  
  # 设置颜色和形状
  scale_color_manual(values = group_colors, name = "") +
  scale_shape_manual(values = c("Not significant" = 16, "Significant" = 17), 
                     name = "") +
  
  # 标签和主题
  labs(
    title = "",
    x = "",
    y = expression(bold(log[bold("2")] * "(Odds Ratio)"))
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 20, -10, 20)
  )

print(p1)



minor_data <- or_data_classified %>%
  filter(grepl("\\.", Disease_Code),  # 包含小数点的为小类
         !Disease_Code %in% c("I85.9", "K63.9", "I25.5", "N08.3", "E11.2","E11.5"))


# 准备绘图数据
plot_data <- minor_data %>%
  # 定义组的顺序和颜色（不换行）
  mutate(Group = factor(Group, 
                        levels = c("Lower WC + Lower proWC", 
                                   "Lower WC + Higher proWC",
                                   "Higher WC + Lower proWC", 
                                   "Higher WC + Higher proWC"),
                        labels = c("Lower WC + Lower proWC", 
                                   "Lower WC + Higher proWC",
                                   "Higher WC + Lower proWC", 
                                   "Higher WC + Higher proWC"),
                        ordered = TRUE)) %>%
  # 计算log2 OR用于绘图
  mutate(log2_or = log2(estimate)) %>%
  # 根据p.adj判断显著性，如果NA则为不显著
  mutate(significance = ifelse(is.na(p.adj), "Not significant",
                               ifelse(p.adj < 0.05, "Significant", "Not significant")),
         point_shape = ifelse(is.na(p.adj), 16,
                              ifelse(p.adj < 0.05, 17, 16)))  # 17=三角形, 16=圆形

# 找出所有数据中log2 OR最高的前5个疾病（正值）- 不分组
top_diseases <- plot_data %>%
  filter(log2_or > 0) %>%  # 只考虑正值的OR
  arrange(desc(log2_or)) %>%
  slice_head(n = 20) %>%
  # 只使用疾病代码作为标签
  mutate(label = Disease_Code)

# 定义颜色方案（不换行标签）
group_colors <- c(
  "Higher WC + Higher proWC" = "#E41A1C",    # 红色
  "Lower WC + Higher proWC" = "#4DAF4A",     # 绿色
  "Higher WC + Lower proWC" = "#377EB8",     # 蓝色
  "Lower WC + Lower proWC" = "#FF7F00"       # 橙色
)

# 创建可视化图表（使用log2 OR）
# 安装并加载ggrepel包
# install.packages("ggrepel")
library(ggrepel)

p2 <- ggplot(plot_data, aes(x = category, y = log2_or, 
                            color = Group, shape = significance,
                            group = interaction(Group, Disease_Code))) +  # 关键行
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  
  # 使用ggrepel避免标签重叠，不显示连接线，标签更靠近点
  geom_text_repel(data = top_diseases,
                  aes(x = category, y = log2_or, label = label),
                  size = 2.5, 
                  color = "black",
                  show.legend = FALSE, 
                  fontface = "bold",
                  max.overlaps = Inf,  # 允许无限重叠，自动调整
                  min.segment.length = 0,  # 总是显示连接线
                  box.padding = 0.2,  # 减小标签周围的填充（原为0.5）
                  point.padding = 0.2,  # 减小点到标签的填充（原为0.5）
                  force = 0.5,  # 减小排斥力（原为1）
                  direction = "both",  # 可以向任何方向移动
                  segment.color = NA,  # 隐藏连接线颜色
                  segment.size = 0) +  # 连接线大小为0
  
  # 设置颜色和形状
  scale_color_manual(values = group_colors, name = "") +
  scale_shape_manual(values = c("Not significant" = 16, "Significant" = 17), 
                     name = "") +
  
  # 标签和主题
  labs(
    title = "",
    x = "",
    y = expression(bold(log[bold("2")] * "(Odds Ratio)"))
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 20, -10, 20)
  )

print(p2)
