setwd("D:/UKB_data/NC")
diet_f <- read.csv("9Diet_participant.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 筛选列名中包含 "Instance.0" 的列
diet_f0 <- diet_f %>%
  select(
    # 选择包含 "Instance 0" 的列
    contains("Instance 0"),
    # 选择不包含任何 "Instance" 的列
    !contains("Instance")
  )

# 将列名中的 " | Instance 0" 替换为空字符串
colnames(diet_f0) <- gsub(" \\| Instance 0", "", colnames(diet_f0))



# 提取指定的饮食相关列
selected_columns <- c(
  "Cooked vegetable intake",
  "Salad / raw vegetable intake", 
  "Fresh fruit intake",
  "Dried fruit intake",
  "Oily fish intake",
  "Non-oily fish intake",
  "Processed meat intake",
  "Poultry intake",
  "Beef intake",
  "Lamb/mutton intake",
  "Pork intake",
  "Cheese intake",
  "Bread intake",
  "Cereal intake",
  "Tea intake",
  "Coffee intake",
  "Water intake",
  "Salt added to food"
)

# 提取这些列
diet_data <- diet_f0 %>%
  select(any_of(selected_columns))

setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
data <- read.csv("complete_data_imputed.csv", header = T, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

names(data)[names(data) == "Time spent watching television (TV)"] <- "Time spent watching television"
common_rows <- Reduce(intersect, list(rownames(complete_data_WC), rownames(results), rownames(data)))

# 按共同的行名合并
complete_data_WC <- complete_data_WC[common_rows, ]
results <- results[common_rows, ]
data <- data[common_rows, ]

# 使用cbind合并
temp_df <- cbind(complete_data_WC, results, data)

setwd("D:/UKB_data")

data_1 <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

temp_df$`UK Biobank assessment centre` <- data_1[rownames(temp_df), "UK Biobank assessment centre | Instance 0"]

temp_df$`UK Biobank assessment centre` <- as.factor(temp_df$`UK Biobank assessment centre`)

selected_df <- subset(temp_df, select = c("WC", "Age", "Sex", "BioX_Adjusted", "BioX_Delta",
                                          "Townsend deprivation index at recruitment",
                                          "Smoking status", "Alcohol intake",
                                          "Ethnic background", "Education","Time spent watching television", "MET minutes per week for Physical activity", "UK Biobank assessment centre", "Energy"))


selected_df_clean <- selected_df
colnames(selected_df_clean) <- gsub(" ", "_", colnames(selected_df_clean))

# 更新协变量列表
covariates_clean <- c("Age", "Sex", "Townsend_deprivation_index_at_recruitment",
                      "Smoking_status", "Alcohol_intake", 
                      "Ethnic_background", "Education","Time_spent_watching_television", "MET_minutes_per_week_for_Physical_activity", "UK_Biobank_assessment_centre", "Energy")


# 匹配
a <- rownames(complete_data_WC)
b <- rownames(diet_data)
c <- intersect(a, b)
WC_data <- complete_data_WC[c,]
dietary_lifestyle_factors <- diet_data[c,]

# 查看每一列的唯一值
#lapply(dietary_lifestyle_factors, unique)


# 对整个数据框进行处理
dietary_lifestyle_factors[] <- lapply(dietary_lifestyle_factors, function(x) {
  # 将-10替换为0.5
  x[x == -10] <- 0.5
  # 将-1和-3替换为0
  x[x %in% c(-1, -3)] <- NA
  return(x)
})

missing_ratio <- colMeans(is.na(dietary_lifestyle_factors))

# 显示每列的缺失值比例
cat("各列缺失值比例:\n")
print(round(missing_ratio, 4))

# 找出缺失值比例大于50%的列
high_missing_cols <- names(dietary_lifestyle_factors)[missing_ratio > 0.5]

# 无




# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 确保数据准备正确
# 假设WC_data包含腰围数据，dietary_lifestyle_factors包含饮食变量
# merged_df已经包含了所有需要的数据

# 定义要分析的饮食变量
diet_vars <- c(
  "Cooked vegetable intake",
  "Salad / raw vegetable intake", 
  "Fresh fruit intake",
  "Dried fruit intake",
  "Oily fish intake",
  "Non-oily fish intake",
  "Processed meat intake",
  "Poultry intake",
  "Beef intake",
  "Lamb/mutton intake",
  "Pork intake",
  "Cheese intake",
  "Bread intake",
  "Cereal intake",
  "Tea intake",
  "Coffee intake",
  "Water intake",
  "Salt added to food"
)

# 初始化结果数据框
results_df <- data.frame()

# 对每个具体食物进行线性回归分析
for (food_item in diet_vars) {
  
  # 创建临时数据框，合并协变量和当前食物摄入量
  temp_data <- data.frame(
    WC = selected_df_clean$WC,
    proWC = selected_df_clean$BioX_Adjusted,
    proWCΔ = selected_df_clean$BioX_Delta,
    food_intake = dietary_lifestyle_factors[[food_item]],
    selected_df_clean[, covariates_clean]
  )
  
  
  # 对三个结局变量分别进行回归
  outcomes <- c("WC", "proWC", "proWCΔ")
  
  for (outcome in outcomes) {
    # 构建回归公式
    formula <- as.formula(paste(outcome, "~ food_intake +", 
                                paste(covariates_clean, collapse = " + ")))
    
    tryCatch({
      # 运行线性回归
      lm_model <- lm(formula, data = temp_data)
      model_summary <- summary(lm_model)
      coef_table <- coef(model_summary)
      
      # 提取食物摄入量的系数
      if ("food_intake" %in% rownames(coef_table)) {
        food_coef <- coef_table["food_intake", ]
        
        # 计算95%置信区间
        conf_int <- confint(lm_model)["food_intake", ]
        
        # 存储结果
        result_row <- data.frame(
          Food_Item = food_item,
          Outcome = outcome,
          Beta = food_coef["Estimate"],
          SE = food_coef["Std. Error"],
          CI_lower = conf_int[1],
          CI_upper = conf_int[2],
          P_value = food_coef["Pr(>|t|)"],
          N = nrow(temp_data),
          stringsAsFactors = FALSE
        )
        
        results_df <- rbind(results_df, result_row)
      }
    }, error = function(e) {
      message(paste("Error in regression for", food_item, "and outcome", outcome, ":", e$message))
    })
  }
}

# 计算Bonferroni校正的P值
# 对每个结局变量分别进行校正
results_df$P_BH <- NA

for (outcome in unique(results_df$Outcome)) {
  outcome_indices <- which(results_df$Outcome == outcome)
  n_tests <- length(outcome_indices)
  results_df$P_BH[outcome_indices] <- p.adjust(results_df$P_value[outcome_indices], method = "BH")
}

results_df1 <- results_df %>%
  filter(!Food_Item %in% c("Tea intake", "Coffee intake", "Water intake"))

library(ggplot2)
library(dplyr)
library(tidyr)

# 准备数据
forest_data <- results_df1 %>%
  mutate(
    # 创建显著性变量
    Significant = ifelse(P_BH < 0.05, "Significant", "Not significant"),
    # 确保结局变量为因子并按正确顺序排列
    Outcome = factor(Outcome, levels = c("WC", "proWC", "proWCΔ")),
    # 从食物名称中移除" intake"，保留"Salt added to food"
    Food_Item_Display = ifelse(
      Food_Item == "Salt added to food", 
      "Salt added to food", 
      gsub(" intake", "", Food_Item)
    ),
    # 创建食物和结局的组合标签
    Food_Outcome = paste(Food_Item_Display, Outcome, sep = " - ")
  ) %>%
  # 按食物和结局排序
  arrange(Food_Item, Outcome)

# 创建浅色版本的颜色调色板
light_color_palette <- c("#F4B88D", "#87CEFA", "#DDA0DD")  # 淡紫色、淡蓝色、淡橙色
color_palette <- c("#FF8C00",  "#1E90FF","#9370DB")        # 原色：紫色、蓝色、橙色

library(ggh4x)

p <- ggplot(forest_data, aes(
  x = Beta,
  y = Food_Item_Display  # 使用修改后的显示名称
)) +
  # 添加零参考线
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  # 添加置信区间
  geom_errorbarh(
    aes(
      xmin = CI_lower,
      xmax = CI_upper,
      color = Outcome
    ),
    height = 0.3,  # 误差线高度
    linewidth = 1,  # 线条粗细
    position = position_dodge(width = 0.3)
  ) +
  # 添加点估计值
  geom_point(
    aes(
      fill = Outcome,
      color = Outcome,
      shape = Significant
    ),
    size = 2.5,       # 圆圈大小
    stroke = 1,     # 边框粗细
    position = position_dodge(width = 0.3)
  ) +
  # 设置颜色和填充
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = light_color_palette) +
  # 设置形状（显著性与非显著性）
  scale_shape_manual(
    values = c("Significant" = 24, "Not significant" = 21),  # 使用有填充的形状
    name = "Significance"  # 图例标题
  ) +
  # 分面显示，每个分面使用不同的x轴刻度
  facet_grid(. ~ Outcome, scales = "free_x") +
  # 为不同分面设置不同的x轴刻度
  facetted_pos_scales(
    x = list(
      # WC分面的x轴刻度
      Outcome == "WC" ~ scale_x_continuous(
        breaks = c(-0.5, 0.0, 0.5, 1.0, 1.5),
        labels = function(x) format(round(x, 1), nsmall = 1)
      ),
      # proWC分面的x轴刻度（与WC相同）
      Outcome == "proWC" ~ scale_x_continuous(
        breaks = c(-0.5, 0.0, 0.5, 1.0, 1.5),
        labels = function(x) format(round(x, 1), nsmall = 1)
      ),
      # proWCΔ分面的x轴刻度
      Outcome == "proWCΔ" ~ scale_x_continuous(
        breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2, 0.3),
        labels = function(x) format(round(x, 1), nsmall = 1)
      )
    )
  ) +
  # 坐标轴和标签
  labs(
    x = "Beta Coefficient",
    y = "",
    title = "",
    color = "Outcome",
    fill = "Outcome"
  ) +
  # 主题设置 - 统一字体和大小
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    axis.text = element_text(size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(size = 9, family = "Arial"),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    strip.text.x = element_text(size = 10, face = "bold", family = "Arial"),
    strip.text.y = element_text(size = 10, face = "bold", family = "Arial"),
    legend.position = "top",  # 图例在顶部
    legend.box = "horizontal",  # 水平排列图例
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(size = 8, family = "Arial"),
    plot.margin = margin(0, 8, 1, 10),  # 调整上边距为图例留空间
  ) +
  # 调整图例顺序和外观
  guides(
    shape = guide_legend(
      title = NULL,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(
        size = 3, 
        fill = "gray80", 
        color = "black"
      )
    ),
    color = "none",  # 隐藏颜色图例（因为分面已经显示）
    fill = "none"    # 隐藏填充图例
  )

# 显示图形
print(p)



# 添加BMI
setwd("D:/UKB_data/NC")
diet_f <- read.csv("9Diet_participant.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 筛选列名中包含 "Instance.0" 的列
diet_f0 <- diet_f %>%
  select(
    # 选择包含 "Instance 0" 的列
    contains("Instance 0"),
    # 选择不包含任何 "Instance" 的列
    !contains("Instance")
  )

# 将列名中的 " | Instance 0" 替换为空字符串
colnames(diet_f0) <- gsub(" \\| Instance 0", "", colnames(diet_f0))



# 提取指定的饮食相关列
selected_columns <- c(
  "Cooked vegetable intake",
  "Salad / raw vegetable intake", 
  "Fresh fruit intake",
  "Dried fruit intake",
  "Oily fish intake",
  "Non-oily fish intake",
  "Processed meat intake",
  "Poultry intake",
  "Beef intake",
  "Lamb/mutton intake",
  "Pork intake",
  "Cheese intake",
  "Bread intake",
  "Cereal intake",
  "Tea intake",
  "Coffee intake",
  "Water intake",
  "Salt added to food"
)

# 提取这些列
diet_data <- diet_f0 %>%
  select(any_of(selected_columns))

setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
data <- read.csv("complete_data_imputed.csv", header = T, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
names(data)[names(data) == "Time spent watching television (TV)"] <- "Time spent watching television"
common_rows <- Reduce(intersect, list(rownames(complete_data_WC), rownames(results), rownames(data)))

# 按共同的行名合并
complete_data_WC <- complete_data_WC[common_rows, ]
results <- results[common_rows, ]
data <- data[common_rows, ]

# 使用cbind合并
temp_df <- cbind(complete_data_WC, results, data)

setwd("D:/UKB_data")

data_1 <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

temp_df$`UK Biobank assessment centre` <- data_1[rownames(temp_df), "UK Biobank assessment centre | Instance 0"]

temp_df$`UK Biobank assessment centre` <- as.factor(temp_df$`UK Biobank assessment centre`)

temp_df$BMI <- data_1[rownames(temp_df), "BMI"]

selected_df <- subset(temp_df, select = c("WC", "Age", "Sex", "BioX_Adjusted", "BioX_Delta",
                                          "Townsend deprivation index at recruitment",
                                          "Smoking status", "Alcohol intake",
                                          "Ethnic background", "Education","Time spent watching television", "MET minutes per week for Physical activity", "UK Biobank assessment centre", "Energy", "BMI"))


selected_df_clean <- selected_df
colnames(selected_df_clean) <- gsub(" ", "_", colnames(selected_df_clean))

# 更新协变量列表
covariates_clean <- c("Age", "Sex", "Townsend_deprivation_index_at_recruitment",
                      "Smoking_status", "Alcohol_intake", 
                      "Ethnic_background", "Education","Time_spent_watching_television", "MET_minutes_per_week_for_Physical_activity", "UK_Biobank_assessment_centre", "Energy", "BMI")


# 匹配
a <- rownames(complete_data_WC)
b <- rownames(diet_data)
c <- intersect(a, b)
WC_data <- complete_data_WC[c,]
dietary_lifestyle_factors <- diet_data[c,]

# 查看每一列的唯一值
#lapply(dietary_lifestyle_factors, unique)


# 对整个数据框进行处理
dietary_lifestyle_factors[] <- lapply(dietary_lifestyle_factors, function(x) {
  # 将-10替换为0.5
  x[x == -10] <- 0.5
  # 将-1和-3替换为0
  x[x %in% c(-1, -3)] <- NA
  return(x)
})

missing_ratio <- colMeans(is.na(dietary_lifestyle_factors))

# 显示每列的缺失值比例
cat("各列缺失值比例:\n")
print(round(missing_ratio, 4))

# 找出缺失值比例大于50%的列
high_missing_cols <- names(dietary_lifestyle_factors)[missing_ratio > 0.5]

# 无



# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 确保数据准备正确
# 假设WC_data包含腰围数据，dietary_lifestyle_factors包含饮食变量
# merged_df已经包含了所有需要的数据

# 定义要分析的饮食变量
diet_vars <- c(
  "Cooked vegetable intake",
  "Salad / raw vegetable intake", 
  "Fresh fruit intake",
  "Dried fruit intake",
  "Oily fish intake",
  "Non-oily fish intake",
  "Processed meat intake",
  "Poultry intake",
  "Beef intake",
  "Lamb/mutton intake",
  "Pork intake",
  "Cheese intake",
  "Bread intake",
  "Cereal intake",
  "Tea intake",
  "Coffee intake",
  "Water intake",
  "Salt added to food"
)


# 初始化结果数据框
results_df <- data.frame()

# 对每个具体食物进行线性回归分析
for (food_item in diet_vars) {
  
  # 创建临时数据框，合并协变量和当前食物摄入量
  temp_data <- data.frame(
    WC = selected_df_clean$WC,
    proWC = selected_df_clean$BioX_Adjusted,
    proWCΔ = selected_df_clean$BioX_Delta,
    food_intake = dietary_lifestyle_factors[[food_item]],
    selected_df_clean[, covariates_clean]
  )
  
  
  # 对三个结局变量分别进行回归
  outcomes <- c("WC", "proWC", "proWCΔ")
  
  for (outcome in outcomes) {
    # 构建回归公式
    formula <- as.formula(paste(outcome, "~ food_intake +", 
                                paste(covariates_clean, collapse = " + ")))
    
    tryCatch({
      # 运行线性回归
      lm_model <- lm(formula, data = temp_data)
      model_summary <- summary(lm_model)
      coef_table <- coef(model_summary)
      
      # 提取食物摄入量的系数
      if ("food_intake" %in% rownames(coef_table)) {
        food_coef <- coef_table["food_intake", ]
        
        # 计算95%置信区间
        conf_int <- confint(lm_model)["food_intake", ]
        
        # 存储结果
        result_row <- data.frame(
          Food_Item = food_item,
          Outcome = outcome,
          Beta = food_coef["Estimate"],
          SE = food_coef["Std. Error"],
          CI_lower = conf_int[1],
          CI_upper = conf_int[2],
          P_value = food_coef["Pr(>|t|)"],
          N = nrow(temp_data),
          stringsAsFactors = FALSE
        )
        
        results_df <- rbind(results_df, result_row)
      }
    }, error = function(e) {
      message(paste("Error in regression for", food_item, "and outcome", outcome, ":", e$message))
    })
  }
}

# 计算Bonferroni校正的P值
# 对每个结局变量分别进行校正
results_df$P_BH <- NA

for (outcome in unique(results_df$Outcome)) {
  outcome_indices <- which(results_df$Outcome == outcome)
  n_tests <- length(outcome_indices)
  results_df$P_BH[outcome_indices] <- p.adjust(results_df$P_value[outcome_indices], method = "BH")
}

results_df1 <- results_df %>%
  filter(!Food_Item %in% c("Tea intake", "Coffee intake", "Water intake"))

library(ggplot2)
library(dplyr)
library(tidyr)

# 准备数据
forest_data <- results_df1 %>%
  mutate(
    # 创建显著性变量
    Significant = ifelse(P_BH < 0.05, "Significant", "Not significant"),
    # 确保结局变量为因子并按正确顺序排列
    Outcome = factor(Outcome, levels = c("WC", "proWC", "proWCΔ")),
    # 从食物名称中移除" intake"，保留"Salt added to food"
    Food_Item_Display = ifelse(
      Food_Item == "Salt added to food", 
      "Salt added to food", 
      gsub(" intake", "", Food_Item)
    ),
    # 创建食物和结局的组合标签
    Food_Outcome = paste(Food_Item_Display, Outcome, sep = " - ")
  ) %>%
  # 按食物和结局排序
  arrange(Food_Item, Outcome)

# 创建浅色版本的颜色调色板
light_color_palette <- c("#F4B88D", "#87CEFA", "#DDA0DD")  # 淡紫色、淡蓝色、淡橙色
color_palette <- c("#FF8C00",  "#1E90FF","#9370DB")        # 原色：紫色、蓝色、橙色

library(ggh4x)

p2 <- ggplot(forest_data, aes(
  x = Beta,
  y = Food_Item_Display
)) +
  # 添加零参考线
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  # 添加置信区间
  geom_errorbarh(
    aes(
      xmin = CI_lower,
      xmax = CI_upper,
      color = Outcome
    ),
    height = 0.3,  # 误差线高度
    linewidth = 1,  # 线条粗细
    position = position_dodge(width = 0.3)
  ) +
  # 添加点估计值
  geom_point(
    aes(
      fill = Outcome,
      color = Outcome,
      shape = Significant
    ),
    size = 2.5,       # 圆圈大小
    stroke = 1,     # 边框粗细
    position = position_dodge(width = 0.3)
  ) +
  # 设置颜色和填充
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = light_color_palette) +
  # 设置形状（显著性与非显著性）
  scale_shape_manual(
    values = c("Significant" = 24, "Not significant" = 21),  # 使用有填充的形状
    name = "Significance"  # 图例标题
  ) +
  # 分面显示，所有分面使用相同的x轴刻度
  facet_grid(. ~ Outcome, scales = "free_x") +
  # 坐标轴和标签
  labs(
    x = "Beta Coefficient",
    y = "",
    title = "",
    color = "Outcome",
    fill = "Outcome"
  ) +
  # 主题设置 - 统一字体和大小
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    axis.text = element_text(size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(size = 9, family = "Arial"),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    strip.text.x = element_text(size = 10, face = "bold", family = "Arial"),
    strip.text.y = element_text(size = 10, face = "bold", family = "Arial"),
    legend.position = "top",  # 图例在顶部
    legend.box = "horizontal",  # 水平排列图例
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(size = 8, family = "Arial"),
    plot.margin = margin(0, 8, 1, 10),  # 调整上边距为图例留空间
  ) +
  # 调整图例顺序和外观
  guides(
    shape = guide_legend(
      title = NULL,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(
        size = 3, 
        fill = "gray80", 
        color = "black"
      )
    ),
    color = "none",  # 隐藏颜色图例（因为分面已经显示）
    fill = "none"    # 隐藏填充图例
  )

# 显示图形
print(p2)







