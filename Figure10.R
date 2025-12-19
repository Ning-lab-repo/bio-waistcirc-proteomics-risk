# 挑选蛋白
setwd("D:/UKB_data/NC/wc_m_v")

# 读取WC的OLS结果
data_ols <- read.csv("TableA_WC_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

# 1. 读取所有TableB_WC_vs_蛋白_Cox文件并合并
cox_file_list <- list.files(pattern = "^TableB_WC_vs_[A-Z][0-9]{2}_Cox\\.csv$")

# 检查Cox文件
if (length(cox_file_list) == 0) {
  warning("未找到匹配的Cox文件")
  data_cox_combined <- NULL
} else {
  cat("找到", length(cox_file_list), "个Cox文件:\n")
  print(cox_file_list)
  
  # 读取并合并所有Cox文件
  cox_data_list <- lapply(cox_file_list, function(file) {
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    data$source_file <- gsub("\\.csv$", "", file)  # 添加源文件名
    
    # 从文件名提取疾病代码：TableB_WC_vs_I08_Cox → I08
    disease_code <- gsub("^TableB_WC_vs_([A-Z][0-9]{2})_Cox\\.csv$", "\\1", file)
    data$Disease_Code <- disease_code
    
    return(data)
  })
  
  data_cox_combined <- do.call(rbind, cox_data_list)
}

setwd("D:/UKB_data/NC")
h1 <- read.csv("all_hr_results_with_WC.csv", header = TRUE, stringsAsFactors = FALSE)
h2 <- read.csv("all_hr_results_with_proWC.csv", header = TRUE, stringsAsFactors = FALSE)
h3 <- read.csv("all_hr_results_with_proWCΔ.csv", header = TRUE, stringsAsFactors = FALSE)

h1_s <- h1 %>% filter(p_BH < 0.05)
h2_s <- h2 %>% filter(p_BH < 0.05)
h3_s <- h3 %>% filter(p_BH < 0.05)

setwd("D:/UKB_data/NC")
# 读入疾病分类
class_data <- read.csv("class.csv", header = TRUE, stringsAsFactors = FALSE)

# 定义函数检查ICD-10编码是否在某个范围内
is_in_icd_range <- function(icd_code, range_str) {
  ranges <- strsplit(range_str, "-")[[1]]
  if (length(ranges) == 2) {
    return(icd_code >= ranges[1] & icd_code <= ranges[2])
  } else {
    return(icd_code == ranges[1])
  }
}

setwd("D:/UKB_data/NC/wc_m_v")

# 初始化一个空的数据框来存储所有结果
all_mediation_data <- data.frame()

# 读取所有TableC_WC_Mediation文件并合并
file_list <- list.files(pattern = "^TableC_WC_Mediation_[A-Z][0-9]{2}\\.csv$")

for (file in file_list) {
  # 从文件名提取ICD-10编码
  icd_code <- gsub("^TableC_WC_Mediation_([A-Z][0-9]{2})\\.csv$", "\\1", file)
  
  # 读取数据
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  data$ACME_p_BH <- p.adjust(data$ACME_p, method = "BH")
  data$ADE_p_BH <- p.adjust(data$ADE_p, method = "BH")
  data$PM_p_BH <- p.adjust(data$PM_p, method = "BH")
  data$ACME_p_bonferroni <- p.adjust(data$ACME_p, method = "bonferroni")
  data$ADE_p_bonferroni <- p.adjust(data$ADE_p, method = "bonferroni")
  data$PM_p_bonferroni <- p.adjust(data$PM_p, method = "bonferroni")
  
  # 检查必要的列是否存在
  if (!all(c("ACME", "ACME_p") %in% colnames(data))) {
    warning(paste("文件", file, "中未找到ACME或ACME_p列，跳过此文件"))
    next
  }
  
  # 添加ICD编码
  data$ICD_10_code <- icd_code
  
  # 确定蛋白质名称列 - 检查常见的蛋白质标识列
  protein_col <- NULL
  possible_protein_cols <- c("Protein", "protein", "SeqId", "seqId", "SOMAMER_ID", "Aptamer", "AptamerName")
  
  for (col in possible_protein_cols) {
    if (col %in% colnames(data)) {
      protein_col <- col
      break
    }
  }
  
  # 如果找到蛋白质列，重命名为统一的"protein"
  if (!is.null(protein_col)) {
    data$protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$protein <- as.character(1:nrow(data))
    warning(paste("文件", file, "中未找到蛋白质标识列，使用行号作为蛋白质名称"))
  }
  
  # 确定疾病分类
  data$Disease <- NA
  for (i in 1:nrow(class_data)) {
    if (is_in_icd_range(icd_code, class_data$ICD.10.code[i])) {
      data$Disease <- class_data$Disease[i]
      break
    }
  }
  
  # 如果未找到匹配的疾病分类，使用"Unknown"
  if (is.na(data$Disease[1])) {
    data$Disease <- "Unknown"
  }
  
  # 选择需要的列并合并
  required_cols <- c("protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH", "ACME_p_bonferroni", 
                     "ADE", "ADE_str", "ADE_p", "ADE_p_BH", "ADE_p_bonferroni", "PM", "PM_str", "PM_p", "PM_p_BH", "PM_p_bonferroni")
  
  # 只保留在required_cols中存在的列
  existing_cols <- required_cols[required_cols %in% colnames(data)]
  all_mediation_data <- rbind(all_mediation_data, data[, existing_cols])
}

# 筛选与h1_s相关的数据
all_mediation_data <- all_mediation_data %>%
  semi_join(h1_s, by = c("ICD_10_code" = "Disease_Code"))

# 筛选满足条件的蛋白：ACME_p_bonferroni和PM_p_bonferroni < 1e-5
significant_proteins <- all_mediation_data %>%
  filter(ACME_p_bonferroni < 1e-5) %>%
  select(protein, ICD_10_code) %>%
  distinct()



# 筛选I10、F17和N18这三个疾病
target_diseases <- c("I10", "F17", "N18")
significant_proteins_target <- significant_proteins %>%
  filter(ICD_10_code %in% target_diseases)



# 1. 从data_ols中筛选显著蛋白
if (nrow(significant_proteins_target) > 0) {
  filtered_ols <- data_ols[data_ols$protein %in% significant_proteins_target$protein, ]
  cat("从OLS数据中筛选到", nrow(filtered_ols), "个蛋白\n")
} else {
  filtered_ols <- data.frame()
  warning("没有显著蛋白，无法进行后续分析")
}

# 2. 从data_cox_combined中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0 && !is.null(data_cox_combined)) {
  filtered_cox_list <- list()
  
  for(i in 1:nrow(significant_proteins_target)) {
    protein_name <- significant_proteins_target$protein[i]
    disease_code <- significant_proteins_target$ICD_10_code[i]
    
    subset_data <- data_cox_combined[data_cox_combined$Disease_Code == disease_code & 
                                       data_cox_combined$protein == protein_name, ]
    
    if(nrow(subset_data) > 0) {
      filtered_cox_list[[paste(disease_code, protein_name, sep="_")]] <- subset_data
    }
  }
  
  if(length(filtered_cox_list) > 0) {
    filtered_cox_combined <- do.call(rbind, filtered_cox_list)
    cat("从Cox数据中筛选到", nrow(filtered_cox_combined), "个蛋白-疾病组合\n")
  } else {
    filtered_cox_combined <- data.frame()
    cat("Cox数据中没有匹配的蛋白-疾病组合\n")
  }
} else {
  filtered_cox_combined <- data.frame()
}

# 3. 从all_mediation_data中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0) {
  # 使用dplyr进行高效过滤
  library(dplyr)
  
  filtered_mediation_combined <- all_mediation_data %>%
    inner_join(significant_proteins_target, 
               by = c("ICD_10_code", "protein"))
  
  cat("从Mediation数据中筛选到", nrow(filtered_mediation_combined), "个蛋白-疾病组合\n")
  
} else {
  filtered_mediation_combined <- data.frame()
  cat("Mediation数据中没有匹配的蛋白-疾病组合\n")
}


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)  # 添加这个包来防止标签重叠

# 准备数据：合并三个数据源
# 1. 从OLS数据获取 Beta (βa)
ols_data <- filtered_ols %>%
  select(protein, Beta) %>%
  rename(beta_a = Beta)

# 2. 从Cox数据计算 Beta (βb) - 注意：HR需要转换为对数尺度
cox_data <- filtered_cox_combined %>%
  select(protein, Disease_Code, HR) %>%
  # 将HR转换为βb系数：βb = log(HR)
  mutate(beta_b = log(HR)) %>%
  select(protein, Disease_Code, beta_b)

# 3. 从Mediation数据获取中介比例
mediation_data <- filtered_mediation_combined %>%
  select(protein, ICD_10_code, PM) %>%
  rename(Disease_Code = ICD_10_code, proportion_mediated = PM)

# 合并所有数据
plot_data <- ols_data %>%
  inner_join(cox_data, by = "protein") %>%
  inner_join(mediation_data, by = c("protein", "Disease_Code")) %>%
  # 计算绝对值
  mutate(
    abs_beta_a = abs(beta_a),
    abs_beta_b = abs(beta_b),
    # 确定效应方向：根据beta_a的符号
    effect_direction = ifelse(beta_a < 0, "Negative", "Positive")
  )



# 确保只包含目标疾病
plot_data <- plot_data %>%
  filter(Disease_Code %in% target_diseases)

# 创建疾病标签
disease_labels <- c(
  "F17" = "F17",
  "I10" = "I10", 
  "N18" = "N18"
)

# 创建图形
p1 <- ggplot(plot_data, aes(x = abs_beta_a, y = abs_beta_b, 
                           color = effect_direction, size = proportion_mediated)) +
  geom_point(alpha = 0.7) +
  # 使用ggrepel防止标签重叠，设置黑色标签和Arial字体
  geom_text_repel(
    aes(label = protein), 
    color = "black",           # 设置标签为黑色
    size = 3, 
    family = "Arial",          # 设置蛋白标签字体为Arial
    fontface = "bold",         # 设置蛋白标签加粗
    max.overlaps = 50,         # 增加最大重叠容忍度
    min.segment.length = 0.1,  # 减少最小线段长度
    box.padding = 0.5,         # 标签周围的填充
    point.padding = 0.3,       # 点周围的填充
    segment.color = "grey50",  # 连接线颜色
    segment.alpha = 0.6,       # 连接线透明度
    force = 1,                 # 排斥力
    show.legend = FALSE
  ) +
  # 分面显示不同疾病
  facet_wrap(~ Disease_Code, labeller = labeller(Disease_Code = disease_labels)) +
  # 坐标轴标签 - 使用expression设置下标并加粗
  labs(
    x = expression(bold("Absolute value of " * bold(beta[a]) * " coefficient")),
    y = expression(bold("Absolute value of " * bold(beta[b]) * " coefficient")), 
    color = "Effect",
    size = "Proportion",
    title = NULL
  ) +
  # 颜色设置：蓝色表示负相关，红色表示正相关
  scale_color_manual(
    values = c("Negative" = "blue", "Positive" = "red"),
    labels = c("Negative" = "Negative", "Positive" = "Positive")
  ) +
  # 调整点的大小范围
  scale_size_continuous(
    range = c(0, 5),
    breaks = seq(0, max(plot_data$proportion_mediated, na.rm = TRUE), 
                 length.out = 5),
    name = "Proportion"
  ) +
  # 主题设置
  theme_minimal(base_family = "Arial") +  # 设置基础字体为Arial
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 主面板边框四边加粗（包含坐标轴线）
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    # 分面标题框四边加粗
    strip.background = element_rect(color = "black", fill = "grey95", linewidth = 1.2),
    strip.text = element_text(face = "bold", family = "Arial", size = 10),  # 分面标签加粗
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 12),
    
    # 坐标轴标签加粗
    axis.title = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.x = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
    
    # 坐标轴刻度加粗且为黑色
    axis.text = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.x = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.y = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    # 移除单独的坐标轴线设置，避免重复
    axis.ticks = element_line(color = "black", linewidth = 0.8),  # 只保留刻度线
    
    # 图例文字加粗
    legend.title = element_text(face = "bold", family = "Arial", size = 10),
    legend.text = element_text(face = "bold", family = "Arial", size = 9),
    
    # 数学符号字体设置（确保β符号正确显示）
    text = element_text(family = "Arial")
  ) +
  # 使用guides()函数调整图例顺序
  guides(
    color = guide_legend(order = 1,  # Effect图例在上面，顺序为1
                         override.aes = list(size = 2)),  # 调整图例中点的大小
    size = guide_legend(order = 2)   # Proportion图例在下面，顺序为2
  )

# 显示图形
print(p1)




# 挑选蛋白
setwd("D:/UKB_data/NC/wc_m_v")

# 读取BioX_Adjusted的OLS结果
data_ols <- read.csv("TableA_BioX_Adjusted_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

# 1. 读取所有TableB_BioX_Adjusted_vs_蛋白_Cox文件并合并
cox_file_list <- list.files(pattern = "^TableB_BioX_Adjusted_vs_[A-Z][0-9]{2}_Cox\\.csv$")

# 检查Cox文件
if (length(cox_file_list) == 0) {
  warning("未找到匹配的Cox文件")
  data_cox_combined <- NULL
} else {
  cat("找到", length(cox_file_list), "个Cox文件:\n")
  print(cox_file_list)
  
  # 读取并合并所有Cox文件
  cox_data_list <- lapply(cox_file_list, function(file) {
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    data$source_file <- gsub("\\.csv$", "", file)  # 添加源文件名
    
    # 从文件名提取疾病代码：TableB_BioX_Adjusted_vs_I08_Cox → I08
    disease_code <- gsub("^TableB_BioX_Adjusted_vs_([A-Z][0-9]{2})_Cox\\.csv$", "\\1", file)
    data$Disease_Code <- disease_code
    
    return(data)
  })
  
  data_cox_combined <- do.call(rbind, cox_data_list)
}

setwd("D:/UKB_data/NC")
h1 <- read.csv("all_hr_results_with_WC.csv", header = TRUE, stringsAsFactors = FALSE)
h2 <- read.csv("all_hr_results_with_proWC.csv", header = TRUE, stringsAsFactors = FALSE)
h3 <- read.csv("all_hr_results_with_proWCΔ.csv", header = TRUE, stringsAsFactors = FALSE)

h1_s <- h1 %>% filter(p_BH < 0.05)
h2_s <- h2 %>% filter(p_BH < 0.05)
h3_s <- h3 %>% filter(p_BH < 0.05)

setwd("D:/UKB_data/NC")
# 读入疾病分类
class_data <- read.csv("class.csv", header = TRUE, stringsAsFactors = FALSE)

# 定义函数检查ICD-10编码是否在某个范围内
is_in_icd_range <- function(icd_code, range_str) {
  ranges <- strsplit(range_str, "-")[[1]]
  if (length(ranges) == 2) {
    return(icd_code >= ranges[1] & icd_code <= ranges[2])
  } else {
    return(icd_code == ranges[1])
  }
}

setwd("D:/UKB_data/NC/wc_m_v")

# 初始化一个空的数据框来存储所有结果
all_mediation_data <- data.frame()

# 读取所有TableC_BioX_Adjusted_Mediation文件并合并
file_list <- list.files(pattern = "^TableC_BioX_Adjusted_Mediation_[A-Z][0-9]{2}\\.csv$")

for (file in file_list) {
  # 从文件名提取ICD-10编码
  icd_code <- gsub("^TableC_BioX_Adjusted_Mediation_([A-Z][0-9]{2})\\.csv$", "\\1", file)
  
  # 读取数据
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  data$ACME_p_BH <- p.adjust(data$ACME_p, method = "BH")
  data$ADE_p_BH <- p.adjust(data$ADE_p, method = "BH")
  data$PM_p_BH <- p.adjust(data$PM_p, method = "BH")
  data$ACME_p_bonferroni <- p.adjust(data$ACME_p, method = "bonferroni")
  data$ADE_p_bonferroni <- p.adjust(data$ADE_p, method = "bonferroni")
  data$PM_p_bonferroni <- p.adjust(data$PM_p, method = "bonferroni")
  
  # 检查必要的列是否存在
  if (!all(c("ACME", "ACME_p") %in% colnames(data))) {
    warning(paste("文件", file, "中未找到ACME或ACME_p列，跳过此文件"))
    next
  }
  
  # 添加ICD编码
  data$ICD_10_code <- icd_code
  
  # 确定蛋白质名称列 - 检查常见的蛋白质标识列
  protein_col <- NULL
  possible_protein_cols <- c("Protein", "protein", "SeqId", "seqId", "SOMAMER_ID", "Aptamer", "AptamerName")
  
  for (col in possible_protein_cols) {
    if (col %in% colnames(data)) {
      protein_col <- col
      break
    }
  }
  
  # 如果找到蛋白质列，重命名为统一的"protein"
  if (!is.null(protein_col)) {
    data$protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$protein <- as.character(1:nrow(data))
    warning(paste("文件", file, "中未找到蛋白质标识列，使用行号作为蛋白质名称"))
  }
  
  # 确定疾病分类
  data$Disease <- NA
  for (i in 1:nrow(class_data)) {
    if (is_in_icd_range(icd_code, class_data$ICD.10.code[i])) {
      data$Disease <- class_data$Disease[i]
      break
    }
  }
  
  # 如果未找到匹配的疾病分类，使用"Unknown"
  if (is.na(data$Disease[1])) {
    data$Disease <- "Unknown"
  }
  
  # 选择需要的列并合并
  required_cols <- c("protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH", "ACME_p_bonferroni", 
                     "ADE", "ADE_str", "ADE_p", "ADE_p_BH", "ADE_p_bonferroni", "PM", "PM_str", "PM_p", "PM_p_BH", "PM_p_bonferroni")
  
  # 只保留在required_cols中存在的列
  existing_cols <- required_cols[required_cols %in% colnames(data)]
  all_mediation_data <- rbind(all_mediation_data, data[, existing_cols])
}

# 筛选与h2_s相关的数据 (修改为h2_s)
all_mediation_data <- all_mediation_data %>%
  semi_join(h2_s, by = c("ICD_10_code" = "Disease_Code"))


# 筛选满足条件的蛋白：ACME_p_bonferroni和PM_p_bonferroni < 1e-5
significant_proteins <- all_mediation_data %>%
  filter(ACME_p_bonferroni < 1e-5) %>%
  select(protein, ICD_10_code) %>%
  distinct()


# 筛选I10、F17和N18这三个疾病
target_diseases <- c("I10", "F17", "N18")
significant_proteins_target <- significant_proteins %>%
  filter(ICD_10_code %in% target_diseases)


# 1. 从data_ols中筛选显著蛋白
if (nrow(significant_proteins_target) > 0) {
  filtered_ols <- data_ols[data_ols$protein %in% significant_proteins_target$protein, ]
  cat("从OLS数据中筛选到", nrow(filtered_ols), "个蛋白\n")
} else {
  filtered_ols <- data.frame()
  warning("没有显著蛋白，无法进行后续分析")
}

# 2. 从data_cox_combined中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0 && !is.null(data_cox_combined)) {
  filtered_cox_list <- list()
  
  for(i in 1:nrow(significant_proteins_target)) {
    protein_name <- significant_proteins_target$protein[i]
    disease_code <- significant_proteins_target$ICD_10_code[i]
    
    subset_data <- data_cox_combined[data_cox_combined$Disease_Code == disease_code & 
                                       data_cox_combined$protein == protein_name, ]
    
    if(nrow(subset_data) > 0) {
      filtered_cox_list[[paste(disease_code, protein_name, sep="_")]] <- subset_data
    }
  }
  
  if(length(filtered_cox_list) > 0) {
    filtered_cox_combined <- do.call(rbind, filtered_cox_list)
    cat("从Cox数据中筛选到", nrow(filtered_cox_combined), "个蛋白-疾病组合\n")
  } else {
    filtered_cox_combined <- data.frame()
    cat("Cox数据中没有匹配的蛋白-疾病组合\n")
  }
} else {
  filtered_cox_combined <- data.frame()
}

# 3. 从all_mediation_data中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0) {
  # 使用dplyr进行高效过滤
  library(dplyr)
  
  filtered_mediation_combined <- all_mediation_data %>%
    inner_join(significant_proteins_target, 
               by = c("ICD_10_code", "protein"))
  
  cat("从Mediation数据中筛选到", nrow(filtered_mediation_combined), "个蛋白-疾病组合\n")
  
} else {
  filtered_mediation_combined <- data.frame()
  cat("Mediation数据中没有匹配的蛋白-疾病组合\n")
}


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)  # 添加这个包来防止标签重叠

# 准备数据：合并三个数据源
# 1. 从OLS数据获取 Beta (βa)
ols_data <- filtered_ols %>%
  select(protein, Beta) %>%
  rename(beta_a = Beta)

# 2. 从Cox数据计算 Beta (βb) - 注意：HR需要转换为对数尺度
cox_data <- filtered_cox_combined %>%
  select(protein, Disease_Code, HR) %>%
  # 将HR转换为βb系数：βb = log(HR)
  mutate(beta_b = log(HR)) %>%
  select(protein, Disease_Code, beta_b)

# 3. 从Mediation数据获取中介比例
mediation_data <- filtered_mediation_combined %>%
  select(protein, ICD_10_code, PM) %>%
  rename(Disease_Code = ICD_10_code, proportion_mediated = PM)

# 合并所有数据
plot_data <- ols_data %>%
  inner_join(cox_data, by = "protein") %>%
  inner_join(mediation_data, by = c("protein", "Disease_Code")) %>%
  # 计算绝对值
  mutate(
    abs_beta_a = abs(beta_a),
    abs_beta_b = abs(beta_b),
    # 确定效应方向：根据beta_a的符号
    effect_direction = ifelse(beta_a < 0, "Negative", "Positive")
  )


# 确保只包含目标疾病
plot_data <- plot_data %>%
  filter(Disease_Code %in% target_diseases)

# 创建疾病标签
disease_labels <- c(
  "F17" = "F17",
  "I10" = "I10", 
  "N18" = "N18"
)

# 创建图形
p2 <- ggplot(plot_data, aes(x = abs_beta_a, y = abs_beta_b, 
                            color = effect_direction, size = proportion_mediated)) +
  geom_point(alpha = 0.7) +
  # 使用ggrepel防止标签重叠，设置黑色标签和Arial字体
  geom_text_repel(
    aes(label = protein), 
    color = "black",           # 设置标签为黑色
    size = 3, 
    family = "Arial",          # 设置蛋白标签字体为Arial
    fontface = "bold",         # 设置蛋白标签加粗
    max.overlaps = 50,         # 增加最大重叠容忍度
    min.segment.length = 0.1,  # 减少最小线段长度
    box.padding = 0.5,         # 标签周围的填充
    point.padding = 0.3,       # 点周围的填充
    segment.color = "grey50",  # 连接线颜色
    segment.alpha = 0.6,       # 连接线透明度
    force = 1,                 # 排斥力
    show.legend = FALSE
  ) +
  # 分面显示不同疾病
  facet_wrap(~ Disease_Code, labeller = labeller(Disease_Code = disease_labels)) +
  # 坐标轴标签 - 使用expression设置下标并加粗
  labs(
    x = expression(bold("Absolute value of " * bold(beta[a]) * " coefficient")),
    y = expression(bold("Absolute value of " * bold(beta[b]) * " coefficient")), 
    color = "Effect",
    size = "Proportion",
    title = NULL
  ) +
  # 颜色设置：蓝色表示负相关，红色表示正相关
  scale_color_manual(
    values = c("Negative" = "blue", "Positive" = "red"),
    labels = c("Negative" = "Negative", "Positive" = "Positive")
  ) +
  # 调整点的大小范围
  scale_size_continuous(
    range = c(0, 5),
    breaks = seq(0, max(plot_data$proportion_mediated, na.rm = TRUE), 
                 length.out = 5),
    name = "Proportion"
  ) +
  # 主题设置
  theme_minimal(base_family = "Arial") +  # 设置基础字体为Arial
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 主面板边框四边加粗（包含坐标轴线）
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    # 分面标题框四边加粗
    strip.background = element_rect(color = "black", fill = "grey95", linewidth = 1.2),
    strip.text = element_text(face = "bold", family = "Arial", size = 10),  # 分面标签加粗
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 12),
    
    # 坐标轴标签加粗
    axis.title = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.x = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
    
    # 坐标轴刻度加粗且为黑色
    axis.text = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.x = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.y = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    # 移除单独的坐标轴线设置，避免重复
    axis.ticks = element_line(color = "black", linewidth = 0.8),  # 只保留刻度线
    
    # 图例文字加粗
    legend.title = element_text(face = "bold", family = "Arial", size = 10),
    legend.text = element_text(face = "bold", family = "Arial", size = 9),
    
    # 数学符号字体设置（确保β符号正确显示）
    text = element_text(family = "Arial")
  )  +
  # 使用guides()函数调整图例顺序
  guides(
    color = guide_legend(order = 1,  # Effect图例在上面，顺序为1
                         override.aes = list(size = 2)),  # 调整图例中点的大小
    size = guide_legend(order = 2)   # Proportion图例在下面，顺序为2
  )

# 显示图形
print(p2)






# 挑选蛋白
setwd("D:/UKB_data/NC/wc_m_v")

# 读取BioX_Delta的OLS结果
data_ols <- read.csv("TableA_BioX_Delta_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

# 1. 读取所有TableB_BioX_Delta_vs_蛋白_Cox文件并合并
cox_file_list <- list.files(pattern = "^TableB_BioX_Delta_vs_[A-Z][0-9]{2}_Cox\\.csv$")

# 检查Cox文件
if (length(cox_file_list) == 0) {
  warning("未找到匹配的Cox文件")
  data_cox_combined <- NULL
} else {
  cat("找到", length(cox_file_list), "个Cox文件:\n")
  print(cox_file_list)
  
  # 读取并合并所有Cox文件
  cox_data_list <- lapply(cox_file_list, function(file) {
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    data$source_file <- gsub("\\.csv$", "", file)  # 添加源文件名
    
    # 从文件名提取疾病代码：TableB_BioX_Delta_vs_I08_Cox → I08
    disease_code <- gsub("^TableB_BioX_Delta_vs_([A-Z][0-9]{2})_Cox\\.csv$", "\\1", file)
    data$Disease_Code <- disease_code
    
    return(data)
  })
  
  data_cox_combined <- do.call(rbind, cox_data_list)
}

setwd("D:/UKB_data/NC")
h1 <- read.csv("all_hr_results_with_WC.csv", header = TRUE, stringsAsFactors = FALSE)
h2 <- read.csv("all_hr_results_with_proWC.csv", header = TRUE, stringsAsFactors = FALSE)
h3 <- read.csv("all_hr_results_with_proWCΔ.csv", header = TRUE, stringsAsFactors = FALSE)

h1_s <- h1 %>% filter(p_BH < 0.05)
h2_s <- h2 %>% filter(p_BH < 0.05)
h3_s <- h3 %>% filter(p_BH < 0.05)

setwd("D:/UKB_data/NC")
# 读入疾病分类
class_data <- read.csv("class.csv", header = TRUE, stringsAsFactors = FALSE)

# 定义函数检查ICD-10编码是否在某个范围内
is_in_icd_range <- function(icd_code, range_str) {
  ranges <- strsplit(range_str, "-")[[1]]
  if (length(ranges) == 2) {
    return(icd_code >= ranges[1] & icd_code <= ranges[2])
  } else {
    return(icd_code == ranges[1])
  }
}

setwd("D:/UKB_data/NC/wc_m_v")

# 初始化一个空的数据框来存储所有结果
all_mediation_data <- data.frame()

# 读取所有TableC_BioX_Delta_Mediation文件并合并
file_list <- list.files(pattern = "^TableC_BioX_Delta_Mediation_[A-Z][0-9]{2}\\.csv$")

for (file in file_list) {
  # 从文件名提取ICD-10编码
  icd_code <- gsub("^TableC_BioX_Delta_Mediation_([A-Z][0-9]{2})\\.csv$", "\\1", file)
  
  # 读取数据
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  data$ACME_p_BH <- p.adjust(data$ACME_p, method = "BH")
  data$ADE_p_BH <- p.adjust(data$ADE_p, method = "BH")
  data$PM_p_BH <- p.adjust(data$PM_p, method = "BH")
  data$ACME_p_bonferroni <- p.adjust(data$ACME_p, method = "bonferroni")
  data$ADE_p_bonferroni <- p.adjust(data$ADE_p, method = "bonferroni")
  data$PM_p_bonferroni <- p.adjust(data$PM_p, method = "bonferroni")
  
  # 检查必要的列是否存在
  if (!all(c("ACME", "ACME_p") %in% colnames(data))) {
    warning(paste("文件", file, "中未找到ACME或ACME_p列，跳过此文件"))
    next
  }
  
  # 添加ICD编码
  data$ICD_10_code <- icd_code
  
  # 确定蛋白质名称列 - 检查常见的蛋白质标识列
  protein_col <- NULL
  possible_protein_cols <- c("Protein", "protein", "SeqId", "seqId", "SOMAMER_ID", "Aptamer", "AptamerName")
  
  for (col in possible_protein_cols) {
    if (col %in% colnames(data)) {
      protein_col <- col
      break
    }
  }
  
  # 如果找到蛋白质列，重命名为统一的"protein"
  if (!is.null(protein_col)) {
    data$protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$protein <- as.character(1:nrow(data))
    warning(paste("文件", file, "中未找到蛋白质标识列，使用行号作为蛋白质名称"))
  }
  
  # 确定疾病分类
  data$Disease <- NA
  for (i in 1:nrow(class_data)) {
    if (is_in_icd_range(icd_code, class_data$ICD.10.code[i])) {
      data$Disease <- class_data$Disease[i]
      break
    }
  }
  
  # 如果未找到匹配的疾病分类，使用"Unknown"
  if (is.na(data$Disease[1])) {
    data$Disease <- "Unknown"
  }
  
  # 选择需要的列并合并
  required_cols <- c("protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH", "ACME_p_bonferroni", 
                     "ADE", "ADE_str", "ADE_p", "ADE_p_BH", "ADE_p_bonferroni", "PM", "PM_str", "PM_p", "PM_p_BH", "PM_p_bonferroni")
  
  # 只保留在required_cols中存在的列
  existing_cols <- required_cols[required_cols %in% colnames(data)]
  all_mediation_data <- rbind(all_mediation_data, data[, existing_cols])
}

# 筛选与h3_s相关的数据 (修改为h3_s)
all_mediation_data <- all_mediation_data %>%
  semi_join(h3_s, by = c("ICD_10_code" = "Disease_Code"))


# 筛选满足条件的蛋白：ACME_p_bonferroni和PM_p_bonferroni < 1e-5
significant_proteins <- all_mediation_data %>%
  filter(ACME_p_bonferroni < 1e-5) %>%
  select(protein, ICD_10_code) %>%
  distinct()


# 筛选I10、F17和N18这三个疾病
target_diseases <- c("I10", "F17", "N18")
significant_proteins_target <- significant_proteins %>%
  filter(ICD_10_code %in% target_diseases)


# 1. 从data_ols中筛选显著蛋白
if (nrow(significant_proteins_target) > 0) {
  filtered_ols <- data_ols[data_ols$protein %in% significant_proteins_target$protein, ]
  cat("从OLS数据中筛选到", nrow(filtered_ols), "个蛋白\n")
} else {
  filtered_ols <- data.frame()
  warning("没有显著蛋白，无法进行后续分析")
}

# 2. 从data_cox_combined中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0 && !is.null(data_cox_combined)) {
  filtered_cox_list <- list()
  
  for(i in 1:nrow(significant_proteins_target)) {
    protein_name <- significant_proteins_target$protein[i]
    disease_code <- significant_proteins_target$ICD_10_code[i]
    
    subset_data <- data_cox_combined[data_cox_combined$Disease_Code == disease_code & 
                                       data_cox_combined$protein == protein_name, ]
    
    if(nrow(subset_data) > 0) {
      filtered_cox_list[[paste(disease_code, protein_name, sep="_")]] <- subset_data
    }
  }
  
  if(length(filtered_cox_list) > 0) {
    filtered_cox_combined <- do.call(rbind, filtered_cox_list)
    cat("从Cox数据中筛选到", nrow(filtered_cox_combined), "个蛋白-疾病组合\n")
  } else {
    filtered_cox_combined <- data.frame()
    cat("Cox数据中没有匹配的蛋白-疾病组合\n")
  }
} else {
  filtered_cox_combined <- data.frame()
}

# 3. 从all_mediation_data中筛选显著蛋白（仅目标疾病）
if (nrow(significant_proteins_target) > 0) {
  # 使用dplyr进行高效过滤
  library(dplyr)
  
  filtered_mediation_combined <- all_mediation_data %>%
    inner_join(significant_proteins_target, 
               by = c("ICD_10_code", "protein"))
  
  cat("从Mediation数据中筛选到", nrow(filtered_mediation_combined), "个蛋白-疾病组合\n")
  
} else {
  filtered_mediation_combined <- data.frame()
  cat("Mediation数据中没有匹配的蛋白-疾病组合\n")
}



# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)  # 添加这个包来防止标签重叠

# 准备数据：合并三个数据源
# 1. 从OLS数据获取 Beta (βa)
ols_data <- filtered_ols %>%
  select(protein, Beta) %>%
  rename(beta_a = Beta)

# 2. 从Cox数据计算 Beta (βb) - 注意：HR需要转换为对数尺度
cox_data <- filtered_cox_combined %>%
  select(protein, Disease_Code, HR) %>%
  # 将HR转换为βb系数：βb = log(HR)
  mutate(beta_b = log(HR)) %>%
  select(protein, Disease_Code, beta_b)

# 3. 从Mediation数据获取中介比例
mediation_data <- filtered_mediation_combined %>%
  select(protein, ICD_10_code, PM) %>%
  rename(Disease_Code = ICD_10_code, proportion_mediated = PM)

# 合并所有数据
plot_data <- ols_data %>%
  inner_join(cox_data, by = "protein") %>%
  inner_join(mediation_data, by = c("protein", "Disease_Code")) %>%
  # 计算绝对值
  mutate(
    abs_beta_a = abs(beta_a),
    abs_beta_b = abs(beta_b),
    # 确定效应方向：根据beta_a的符号
    effect_direction = ifelse(beta_a < 0, "Negative", "Positive")
  )


# 确保只包含目标疾病
plot_data <- plot_data %>%
  filter(Disease_Code %in% target_diseases)

# 创建疾病标签
disease_labels <- c(
  "F17" = "F17",
  "I10" = "I10", 
  "N18" = "N18"
)

# 创建图形
p3 <- ggplot(plot_data, aes(x = abs_beta_a, y = abs_beta_b, 
                            color = effect_direction, size = proportion_mediated)) +
  geom_point(alpha = 0.7) +
  # 使用ggrepel防止标签重叠，设置黑色标签和Arial字体
  geom_text_repel(
    aes(label = protein), 
    color = "black",           # 设置标签为黑色
    size = 3, 
    family = "Arial",          # 设置蛋白标签字体为Arial
    fontface = "bold",         # 设置蛋白标签加粗
    max.overlaps = 50,         # 增加最大重叠容忍度
    min.segment.length = 0.1,  # 减少最小线段长度
    box.padding = 0.5,         # 标签周围的填充
    point.padding = 0.3,       # 点周围的填充
    segment.color = "grey50",  # 连接线颜色
    segment.alpha = 0.6,       # 连接线透明度
    force = 1,                 # 排斥力
    show.legend = FALSE
  ) +
  # 分面显示不同疾病
  facet_wrap(~ Disease_Code, labeller = labeller(Disease_Code = disease_labels)) +
  # 坐标轴标签 - 使用expression设置下标并加粗
  labs(
    x = expression(bold("Absolute value of " * bold(beta[a]) * " coefficient")),
    y = expression(bold("Absolute value of " * bold(beta[b]) * " coefficient")), 
    color = "Effect",
    size = "Proportion",
    title = NULL
  ) +
  # 颜色设置：蓝色表示负相关，红色表示正相关
  scale_color_manual(
    values = c("Negative" = "blue", "Positive" = "red"),
    labels = c("Negative" = "Negative", "Positive" = "Positive")
  ) +
  # 调整点的大小范围
  scale_size_continuous(
    range = c(0, 5),
    breaks = seq(0, max(plot_data$proportion_mediated, na.rm = TRUE), 
                 length.out = 5),
    name = "Proportion"
  ) +
  # 主题设置
  theme_minimal(base_family = "Arial") +  # 设置基础字体为Arial
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 主面板边框四边加粗（包含坐标轴线）
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    # 分面标题框四边加粗
    strip.background = element_rect(color = "black", fill = "grey95", linewidth = 1.2),
    strip.text = element_text(face = "bold", family = "Arial", size = 10),  # 分面标签加粗
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 12),
    
    # 坐标轴标签加粗
    axis.title = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.x = element_text(face = "bold", family = "Arial", size = 12),
    axis.title.y = element_text(face = "bold", family = "Arial", size = 12),
    
    # 坐标轴刻度加粗且为黑色
    axis.text = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.x = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    axis.text.y = element_text(face = "bold", family = "Arial", color = "black", size = 10),
    # 移除单独的坐标轴线设置，避免重复
    axis.ticks = element_line(color = "black", linewidth = 0.8),  # 只保留刻度线
    
    # 图例文字加粗
    legend.title = element_text(face = "bold", family = "Arial", size = 10),
    legend.text = element_text(face = "bold", family = "Arial", size = 9),
    
    # 数学符号字体设置（确保β符号正确显示）
    text = element_text(family = "Arial")
  ) +
  # 使用guides()函数调整图例顺序
  guides(
    color = guide_legend(order = 1,  # Effect图例在上面，顺序为1
                         override.aes = list(size = 2)),  # 调整图例中点的大小
    size = guide_legend(order = 2)   # Proportion图例在下面，顺序为2
  )

# 显示图形
print(p3)


library(cowplot)

# 提取图例
legend_only <- get_legend(p3)

# 将图例保存为图片
ggsave("legend_only.png", plot = legend_only, 
       width = 4, height = 3, dpi = 300, bg = "white")

cat("图例已保存为 legend_only.png\n")


















