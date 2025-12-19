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
    warning(paste("文件", file, "中未找到PM或PM_p列，跳过此文件"))
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
  
  # 如果找到蛋白质列，重命名为统一的"Protein"
  if (!is.null(protein_col)) {
    data$Protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$Protein <- as.character(1:nrow(data))
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
  required_cols <- c("Protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH","ACME_p_bonferroni", "ADE", "ADE_str", "ADE_p", "ADE_p_BH","ADE_p_bonferroni","PM","PM_str","PM_p","PM_p_BH","PM_p_bonferroni")
  all_mediation_data <- rbind(all_mediation_data, data[, required_cols])
}

all_mediation_data <- all_mediation_data %>%
  semi_join(h1_s, by = c("ICD_10_code" = "Disease_Code"))


all_mediation_data <- all_mediation_data %>%
  mutate(
    Significant = case_when(
      ACME_p_bonferroni < 0.05 ~ "Significant",
      is.na(ACME_p_bonferroni) ~ "Not significant",  # 处理NA值
      TRUE ~ "Not significant"
    )
  )

all_mediation_data$Significant <- factor(all_mediation_data$Significant, 
                                         levels = c("Significant", "Not significant"))


# 创建曼哈顿图 - 保留原来NC风格配色
pm1 <- ggplot(all_mediation_data, aes(x = Disease, y = ACME, 
                                      shape = Significant)) +
  geom_point(position = position_jitter(width = 0.3, height = 0), size = 2, 
             aes(color = Significant, alpha = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 0.8) +
  
  # 保留原来的NC风格配色
  scale_color_manual(values = c("Significant" = "#F4B88D",  # 橙黄色
                                "Not significant" = "#A9A9A9"),  # 浅灰色
                     na.value = "#A9A9A9",
                     name = "Significance") +
  
  # 设置形状：显著用实心三角形，不显著用圆形
  scale_shape_manual(values = c("Significant" = 17, "Not significant" = 16),
                     na.value = 16,
                     name = "Significance") +
  
  # 设置透明度：显著点更突出
  scale_alpha_manual(values = c("Significant" = 0.9, "Not significant" = 0.4),
                     na.value = 0.4,
                     name = "Significance") +
  
  labs(
    title = "WC",
    x = "",
    y = "ACME"
  ) +
  
  # 优化主题 - 统一字体和字号
  theme_minimal(base_family = "Arial") +
  theme(
    # 统一字体设置
    text = element_text(family = "Arial", size = 10),
    
    # 网格线设置：只保留y轴主要网格线
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # 坐标轴线条
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    # 坐标轴标题 - 统一字号和粗细
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 15)),
    
    # 坐标轴文本 - 统一字号
    axis.text.x = element_text(
      angle = 60, 
      hjust = 1, 
      size = 9,
      face = "bold",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 9,
      family = "Arial"
    ),
    
    # 标题设置（如需添加标题）
    plot.title = element_text(
      hjust = 0.5, 
      face = "bold", 
      size = 12,
      family = "Arial",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "black",
      family = "Arial",
      margin = margin(b = 15)
    ),
    
    # 图例
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 10),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 9),
    
    # 边距
    plot.margin = margin(20, 20, 20, 20),
    
    # 背景
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  
  # 图例设置 - 合并图例
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      title = NULL,
      override.aes = list(
        shape = c(17, 16), 
        size = 3,
        alpha = c(0.9, 0.4)
      )
    ),
    shape = "none",
    alpha = "none"
  )

# 显示图形
print(pm1)



# 初始化一个空的数据框来存储所有结果
all_mediation_data <- data.frame()

# 读取所有TableC_WC_Mediation文件并合并
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
    warning(paste("文件", file, "中未找到PM或PM_p列，跳过此文件"))
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
  
  # 如果找到蛋白质列，重命名为统一的"Protein"
  if (!is.null(protein_col)) {
    data$Protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$Protein <- as.character(1:nrow(data))
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
  required_cols <- c("Protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH","ACME_p_bonferroni", "ADE", "ADE_str", "ADE_p", "ADE_p_BH","ADE_p_bonferroni","PM","PM_str","PM_p","PM_p_BH","PM_p_bonferroni")
  all_mediation_data <- rbind(all_mediation_data, data[, required_cols])
}


all_mediation_data <- all_mediation_data %>%
  semi_join(h2_s, by = c("ICD_10_code" = "Disease_Code"))

all_mediation_data <- all_mediation_data %>%
  mutate(
    Significant = case_when(
      ACME_p_bonferroni < 0.05 ~ "Significant",
      is.na(ACME_p_bonferroni) ~ "Not significant",  # 处理NA值
      TRUE ~ "Not significant"
    )
  )

all_mediation_data$Significant <- factor(all_mediation_data$Significant, 
                                         levels = c("Significant", "Not significant"))

# 创建曼哈顿图 - NC风格配色
pm2 <- ggplot(all_mediation_data, aes(x = Disease, y = ACME, 
                                      shape = Significant)) +
  geom_point(position = position_jitter(width = 0.3, height = 0), size = 2, 
             aes(color = Significant, alpha = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 0.8) +
  
  # NC风格配色：显著用深蓝色，不显著用浅灰色
  scale_color_manual(values = c("Significant" = "#87CEEB",  # 深蓝色 - NC常用色
                                "Not significant" = "#A9A9A9"),  # 浅灰色
                     na.value = "#A9A9A9",
                     name = "Significance") +
  
  # 设置形状：显著用实心三角形，不显著用圆形
  scale_shape_manual(values = c("Significant" = 17, "Not significant" = 16),
                     na.value = 16,
                     name = "Significance") +
  
  # 设置透明度：显著点更突出
  scale_alpha_manual(values = c("Significant" = 0.9, "Not significant" = 0.4),
                     na.value = 0.4,
                     name = "Significance") +
  
  labs(
    title = "proWC",
    x = "",
    y = "ACME"
  ) +
  
  # 优化主题 - 统一字体和字号
  theme_minimal(base_family = "Arial") +
  theme(
    # 统一字体设置
    text = element_text(family = "Arial", size = 10),
    
    # 网格线设置：只保留y轴主要网格线
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # 坐标轴线条
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    # 坐标轴标题 - 统一字号和粗细
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 15)),
    
    # 坐标轴文本 - 统一字号
    axis.text.x = element_text(
      angle = 60, 
      hjust = 1, 
      size = 9,
      face = "bold",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 9,
      family = "Arial"
    ),
    
    # 标题设置
    plot.title = element_text(
      hjust = 0.5, 
      face = "bold", 
      size = 12,
      family = "Arial",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "black",
      family = "Arial",
      margin = margin(b = 15)
    ),
    
    # 图例 - 顶部水平排列
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 10),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 9),
    
    # 边距
    plot.margin = margin(20, 20, 20, 20),
    
    # 背景
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  
  # 图例设置 - 合并图例
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      title = NULL,
      override.aes = list(
        shape = c(17, 16), 
        size = 3,
        alpha = c(0.9, 0.4)
      )
    ),
    shape = "none",
    alpha = "none"
  )

# 显示图形
print(pm2)




# 初始化一个空的数据框来存储所有结果
all_mediation_data <- data.frame()

# 读取所有TableC_WC_Mediation文件并合并
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
    warning(paste("文件", file, "中未找到PM或PM_p列，跳过此文件"))
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
  
  # 如果找到蛋白质列，重命名为统一的"Protein"
  if (!is.null(protein_col)) {
    data$Protein <- data[[protein_col]]
  } else {
    # 如果没有找到蛋白质列，使用第一列作为蛋白质标识（通常是行名或标识符）
    data$Protein <- as.character(1:nrow(data))
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
  required_cols <- c("Protein", "ICD_10_code", "Disease", "ACME", "ACME_str", "ACME_p", "ACME_p_BH","ACME_p_bonferroni", "ADE", "ADE_str", "ADE_p", "ADE_p_BH","ADE_p_bonferroni","PM","PM_str","PM_p","PM_p_BH","PM_p_bonferroni")
  all_mediation_data <- rbind(all_mediation_data, data[, required_cols])
}


all_mediation_data <- all_mediation_data %>%
  semi_join(h3_s, by = c("ICD_10_code" = "Disease_Code"))


all_mediation_data <- all_mediation_data %>%
  mutate(
    Significant = case_when(
      ACME_p_bonferroni < 0.05 ~ "Significant",
      is.na(ACME_p_bonferroni) ~ "Not significant",  # 处理NA值
      TRUE ~ "Not significant"
    )
  )

all_mediation_data$Significant <- factor(all_mediation_data$Significant, 
                                         levels = c("Significant", "Not significant"))

# 创建曼哈顿图 - NC风格配色
pm3 <- ggplot(all_mediation_data, aes(x = Disease, y = ACME, 
                                      shape = Significant)) +
  geom_point(position = position_jitter(width = 0.3, height = 0), size = 2, 
             aes(color = Significant, alpha = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7, linewidth = 0.8) +
  
  # NC风格配色：显著用紫红色，不显著用浅灰色
  scale_color_manual(values = c("Significant" = "#DDA0DD",  # 紫红色
                                "Not significant" = "#A9A9A9"),  # 浅灰色
                     na.value = "#A9A9A9",
                     name = "Significance") +
  
  # 设置形状：显著用实心三角形，不显著用圆形
  scale_shape_manual(values = c("Significant" = 17, "Not significant" = 16),
                     na.value = 16,
                     name = "Significance") +
  
  # 设置透明度：显著点更突出
  scale_alpha_manual(values = c("Significant" = 0.9, "Not significant" = 0.4),
                     na.value = 0.4,
                     name = "Significance") +
  
  labs(
    title = "proWCΔ",
    x = "",
    y = "ACME"
  ) +
  
  # 优化主题 - 统一字体和字号
  theme_minimal(base_family = "Arial") +
  theme(
    # 统一字体设置
    text = element_text(family = "Arial", size = 10),
    
    # 网格线设置：只保留y轴主要网格线
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # 坐标轴线条
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    # 坐标轴标题 - 统一字号和粗细
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 15)),
    
    # 坐标轴文本 - 统一字号
    axis.text.x = element_text(
      angle = 60, 
      hjust = 1, 
      size = 9,
      face = "bold",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 9,
      family = "Arial"
    ),
    
    # 标题设置
    plot.title = element_text(
      hjust = 0.5, 
      face = "bold", 
      size = 12,
      family = "Arial",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "black",
      family = "Arial",
      margin = margin(b = 15)
    ),
    
    # 图例 - 顶部水平排列
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 10),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 9),
    
    # 边距
    plot.margin = margin(20, 20, 20, 20),
    
    # 背景
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  ) +
  
  # 图例设置 - 合并图例
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      title = NULL,
      override.aes = list(
        shape = c(17, 16), 
        size = 3,
        alpha = c(0.9, 0.4)
      )
    ),
    shape = "none",
    alpha = "none"
  )

# 显示图形
print(pm3)






