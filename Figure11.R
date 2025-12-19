setwd("D:/UKB_data")
data <- read.csv("pro53013_新诊_newdead.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
setwd("D:/UKB_data/NC")
ft <- read.csv("1Baseline_characteristics_participant.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

birth_data <- data.frame(
  Birth_Date = paste(ft$`Year of birth`, ft$`Month of birth`, sep = "/"),
  row.names = rownames(ft)
)


aas_data <- new_data <- data.frame(
  aas_Date = paste(data$date_attending_assessment_centre),
  row.names = rownames(data)
)

aas_data$aas_Date <- format(as.Date(aas_data$aas_Date), "%Y/%m")

# 确保两个数据集的行的顺序和名称一致
common_rows <- intersect(rownames(birth_data), rownames(aas_data))
birth_data_common <- birth_data[common_rows, , drop = FALSE]
aas_data_common <- aas_data[common_rows, , drop = FALSE]

# 按列合并
merged_data <- cbind(birth_data_common, aas_data_common)

# 提取年份和月份
aas_year <- as.numeric(substr(merged_data$aas_Date, 1, 4))
aas_month <- as.numeric(substr(merged_data$aas_Date, 6, 7))

birth_year <- as.numeric(substr(merged_data$Birth_Date, 1, 4))
birth_month <- as.numeric(substr(merged_data$Birth_Date, 6, 7))

# 计算月份差
merged_data$Age_month <- (aas_year - birth_year) * 12 + (aas_month - birth_month)




setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("WC_death_time.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df <- cbind(complete_data_WC, results)

# 直接匹配行名添加BMI列
temp_df$BMI <- data[rownames(temp_df), "BMI"]
temp_df$Age_month <- merged_data[rownames(temp_df), "Age_month"]


# 加载必要的包
library(ggplot2)
library(dplyr)

# 准备数据
plot_data <- temp_df[, c("Age", "Sex", "WC", "BioX_Adjusted", "Age_month")]
colnames(plot_data)[4] <- "proWC"

# 将性别转换为因子变量
plot_data$Sex <- factor(plot_data$Sex, levels = c(1, 0), labels = c("Male", "Female"))



# 进行检验
wilcox_test_WC <- wilcox.test(WC ~ Sex, data = plot_data)
wilcox_test_proWC <- wilcox.test(proWC ~ Sex, data = plot_data)



# 如果您想要两个标注稍微错开一点（避免重叠），可以这样：
signif_data <- data.frame(
  Sex = factor(c("Female", "Female"), levels = c("Male", "Female")),
  Measurement = c("WC", "proWC"),
  x_pos = c(775, 775),  # x轴位置都是775
  y_pos = c(82, 81),    # y轴位置稍微错开，WC在82，proWC在78
  
  label = c(
    ifelse(wilcox_test_WC$p.value < 0.001, "***",
           ifelse(wilcox_test_WC$p.value < 0.01, "**",
                  ifelse(wilcox_test_WC$p.value < 0.05, "*", "ns"))),
    ifelse(wilcox_test_proWC$p.value < 0.001, "***",
           ifelse(wilcox_test_proWC$p.value < 0.01, "**",
                  ifelse(wilcox_test_proWC$p.value < 0.05, "*", "ns")))
  ),
  color_group = c("WC", "proWC")
)


p_smooth <- ggplot(plot_data, aes(x = Age_month)) +
  # 只保留平滑曲线层
  geom_smooth(aes(y = WC, color = "WC"), method = "loess", se = TRUE) +
  geom_smooth(aes(y = proWC, color = "proWC"), method = "loess", se = TRUE) +
  
  # 添加显著性星号标注
  geom_text(data = signif_data, 
            aes(x = x_pos, y = y_pos, label = label, color = color_group),
            size = 5, fontface = "bold", family = "Arial", show.legend = FALSE) +
  
  facet_wrap(~ Sex, ncol = 2) +
  labs(x = "Age", y = "Value", 
       color = "") +
  
  # 简化颜色标度 - 只需要两种颜色
  scale_color_manual(
    name = "",
    values = c(
      "WC" = "#F4B88D",    
      "proWC" = "#87CEEB"  
    ),
    labels = c(
      "WC" = "WC",
      "proWC" = "proWC"
    )
  ) +
  
  # 主题设置 - 统一字体和大小
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),  # 改为11
    axis.text = element_text(size = 9, family = "Arial"),  # 改为9
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.background = element_blank(),
    legend.text = element_text(size = 8, family = "Arial"),  # 改为8
    strip.background = element_rect(fill = "lightgray", linewidth = 0.8),
    strip.text = element_text(face = "bold", family = "Arial", size = 10),  # 改为10
    legend.key = element_rect(fill = NA, color = NA),
    legend.margin = margin(l = 0, r = 0, t = 0, b = 0),
    plot.margin = margin(10, 10, 10, 10)
  )

# 显示平滑曲线图
p_smooth

# 箱线图
# 准备数据
# 准备数据
plot_data <- temp_df[, c("Age", "Sex", "WC", "BioX_Adjusted", "Age_month", "BMI")]
colnames(plot_data)[4] <- "proWC"

# 将性别转换为因子变量
plot_data$Sex <- factor(plot_data$Sex, levels = c(1, 0), labels = c("Male", "Female"))

# 创建BMI分组
plot_data$BMI_group <- cut(plot_data$BMI,
                           breaks = c(-Inf, 18.5, 24.9, 29.9, Inf),
                           labels = c("Underweight", "NW", "OW", "OB"),
                           right = TRUE)

# 移除NA值
plot_data <- plot_data[complete.cases(plot_data$BMI_group), ]

# 分别对男性和女性进行Kruskal-Wallis检验
kruskal_test_WC_male <- kruskal.test(WC ~ BMI_group, data = plot_data[plot_data$Sex == "Male", ])
kruskal_test_proWC_male <- kruskal.test(proWC ~ BMI_group, data = plot_data[plot_data$Sex == "Male", ])
kruskal_test_WC_female <- kruskal.test(WC ~ BMI_group, data = plot_data[plot_data$Sex == "Female", ])
kruskal_test_proWC_female <- kruskal.test(proWC ~ BMI_group, data = plot_data[plot_data$Sex == "Female", ])

# 如果您想稍微错开显示两个标注（避免重叠），可以这样：
signif_data <- data.frame(
  Sex = factor(rep(c("Male", "Female"), each = 2), levels = c("Male", "Female")),
  Measurement = rep(c("WC", "proWC"), 2),
  x_pos = rep("OB", 4),  # x轴位置为OB
  
  # y轴位置稍微错开
  y_pos = c(
    42, 50,  # Male分面：WC在42，proWC在44
    42, 50   # Female分面：WC在42，proWC在44
  ),
  
  label = c(
    # Male分面的显著性
    ifelse(kruskal_test_WC_male$p.value < 0.001, "***",
           ifelse(kruskal_test_WC_male$p.value < 0.01, "**",
                  ifelse(kruskal_test_WC_male$p.value < 0.05, "*", ""))),
    ifelse(kruskal_test_proWC_male$p.value < 0.001, "***",
           ifelse(kruskal_test_proWC_male$p.value < 0.01, "**",
                  ifelse(kruskal_test_proWC_male$p.value < 0.05, "*", ""))),
    # Female分面的显著性
    ifelse(kruskal_test_WC_female$p.value < 0.001, "***",
           ifelse(kruskal_test_WC_female$p.value < 0.01, "**",
                  ifelse(kruskal_test_WC_female$p.value < 0.05, "*", ""))),
    ifelse(kruskal_test_proWC_female$p.value < 0.001, "***",
           ifelse(kruskal_test_proWC_female$p.value < 0.01, "**",
                  ifelse(kruskal_test_proWC_female$p.value < 0.05, "*", "")))
  )
)


# 将数据转换为长格式以便于绘图
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = c(WC, proWC),
    names_to = "Measurement",
    values_to = "Value"
  )

# 创建传统的box and whisker plot（箱体上下线表示最大值最小值）
p3 <- ggplot(plot_data_long, aes(x = BMI_group, y = Value, fill = Measurement)) +
  # 先绘制须线（黑色）
  stat_boxplot(
    geom = "errorbar",  # 绘制须线
    width = 0.5,
    position = position_dodge(0.8),
    color = "black"  # 须线改为黑色
  ) +
  # 再绘制箱体（覆盖须线）
  geom_boxplot(
    # 箱体设置：上下边缘表示最大值和最小值
    coef = 0,  # 设置为0表示不使用IQR计算须线，直接使用最大值最小值
    outlier.shape = NA,  # 不显示异常点，因为上下边缘就是最大值最小值
    alpha = 0.9,  # 增加不透明度以更好地覆盖须线
    position = position_dodge(0.8),
    width = 0.6,
    color = "black",  # 箱体边框为黑色
    size = 0.6,  # 稍微增加边框粗细
    fatten = 1.5  # 控制箱体中位线的粗细
  ) +
  
  # 添加显著性标注 - 在Male和Female分面中都显示，使用相应颜色
  geom_text(
    data = signif_data,
    aes(x = x_pos, y = y_pos, label = label, color = Measurement),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    family = "Arial",
    show.legend = FALSE
  ) +
  
  facet_wrap(~ Sex, ncol = 2) +
  
  # 添加y轴限制：从0到160
  coord_cartesian(ylim = c(0, 160)) +
  
  labs(
    title = "",
    x = "",
    y = "Value",
    fill = ""
  ) +
  scale_fill_manual(
    values = c(
      "WC" = "#F4B88D",    # 橙色
      "proWC" = "#87CEEB"  # 蓝色
    ),
    labels = c("WC", "proWC")
  ) +
  scale_color_manual(
    values = c(
      "WC" = "#F4B88D",    # 橙色
      "proWC" = "#87CEEB"  # 蓝色
    )
  ) +
  
  # 主题设置 - 统一字体和大小
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),  # 全局字体
    
    # 去除网格线
    panel.grid = element_blank(),
    
    # 添加坐标轴线
    axis.line = element_line(color = "black", linewidth = 0.8),
    
    # 添加刻度线
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    # 坐标轴标题样式
    axis.title = element_text(face = "bold", size = 11, family = "Arial"),
    
    # 坐标轴刻度文字样式
    axis.text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(size = 9, family = "Arial"),
    axis.text.y = element_text(size = 9, family = "Arial"),
    
    # 图例样式
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9, family = "Arial"),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    legend.margin = margin(l = 0, r = 0, t = 0, b = 5),
    
    # 分面标题样式
    strip.background = element_rect(fill = "lightgray", linewidth = 0.8),
    strip.text = element_text(face = "bold", family = "Arial", size = 10),
    
    # 图表外边距
    plot.margin = margin(10, 10, 10, 10)
  )







