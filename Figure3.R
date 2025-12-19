library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(grid)

setwd("D:/UKB_data/NC")
or_result <- read.csv("hbh_or_results.csv", header = T, fileEncoding = "GBK")
# 读入疾病分类
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
#diseases <- do.call(rbind, categorized_diseases)

# 将疾病代码分类到对应的类别中
for(category_name in names(categorized_diseases)) {
  or_result$category[or_result$Disease_Code %in% categorized_diseases[[category_name]]] <- category_name
}

# 移除未分类的数据
or_result_classified <- or_result %>% filter(!is.na(category))

# 定义所有OR值类型和对应的标签
or_types <- c("Odds_ratio_proBMIΔ", "Odds_ratio_BMI", "Odds_ratio_proBMI", 
              "Odds_ratio_proWCΔ", "Odds_ratio_WC", "Odds_ratio_proWC")

or_labels <- c(
  "Odds_ratio_proBMIΔ" = "proBMIΔ",
  "Odds_ratio_BMI" = "BMI", 
  "Odds_ratio_proBMI" = "proBMI",
  "Odds_ratio_proWCΔ" = "proWCΔ",
  "Odds_ratio_WC" = "WC",
  "Odds_ratio_proWC" = "proWC"
)

# 重塑OR数据（包含人数）
or_result_long <- or_result_classified %>%
  select(Disease_Code, category, quintile, 
         Odds_ratio_proBMIΔ, Odds_ratio_BMI, Odds_ratio_proBMI,
         Odds_ratio_proWCΔ, Odds_ratio_WC, Odds_ratio_proWC,
         starts_with("Number_of_participants")) %>%
  pivot_longer(
    cols = c(starts_with("Odds_ratio"), starts_with("Number_of_participants")),
    names_to = "variable_type",
    values_to = "value"
  ) %>%
  mutate(
    metric_type = ifelse(grepl("^Odds_ratio", variable_type), "or", "participants"),
    or_type = case_when(
      grepl("proBMIΔ", variable_type) ~ "proBMIΔ",
      grepl("proBMI", variable_type) ~ "proBMI",
      grepl("BMI", variable_type) & !grepl("pro", variable_type) ~ "BMI",
      grepl("proWCΔ", variable_type) ~ "proWCΔ",
      grepl("proWC", variable_type) ~ "proWC",
      grepl("WC", variable_type) & !grepl("pro", variable_type) ~ "WC",
      TRUE ~ "other"
    )
  ) %>%
  select(-variable_type) %>%
  pivot_wider(
    names_from = metric_type,
    values_from = value
  ) %>%
  rename(or_value = or, participants = participants)

# 重塑p值数据
p_value_long <- or_result_classified %>%
  select(Disease_Code, category, quintile, 
         p.adj_proBMIΔ, p.adj_BMI, p.adj_proBMI,
         p.adj_proWCΔ, p.adj_WC, p.adj_proWC) %>%
  pivot_longer(
    cols = starts_with("p.adj"),
    names_to = "p_type",
    values_to = "p_value"
  ) %>%
  mutate(p_type = gsub("^p\\.adj_", "", p_type))

# 合并OR值、人数和p值
final_result_long <- or_result_long %>%
  left_join(p_value_long, 
            by = c("Disease_Code", "category", "quintile", 
                   "or_type" = "p_type"))

# 方法2：根据是否包含小数点来区分，并排除N47
major_category_data <- final_result_long %>%
  filter(!grepl("\\.", Disease_Code) & Disease_Code != "N47")


# 定义类别顺序：将类别放在最后
category_order <- c(unique(or_result_classified$category))
major_category_data$category <- factor(major_category_data$category, 
                                       levels = category_order)

# 添加显著性标记列
major_category_data <- major_category_data %>%
  mutate(
    significant = case_when(
      p_value < 0.05 ~ "Significant",
      is.na(p_value) ~ "Not significant",  # 处理NA值
      TRUE ~ "Not significant"
    ),
    shape = case_when(
      p_value < 0.05 ~ "triangle",
      is.na(p_value) ~ "circle",  # 处理NA值
      TRUE ~ "circle"
    )
  )

# 为每个OR指标创建图形，展示不同疾病类别
plot_list_by_or_type <- list()

# 有图例
for (or_type in unique(major_category_data$or_type)) {
  # 筛选当前OR类型的数据
  or_type_data <- major_category_data %>% 
    filter(or_type == !!or_type) %>%
    arrange(category, Disease_Code, quintile) %>%
    mutate(log2_or = log2(or_value))
  
  # 找出OR为正值时最大的前5个疾病
  top_positive_diseases <- or_type_data %>%
    filter(log2_or > 0) %>%
    group_by(Disease_Code) %>%
    slice_max(order_by = log2_or, n = 1) %>%
    ungroup() %>%
    arrange(desc(log2_or)) %>%
    distinct(Disease_Code, .keep_all = TRUE) %>%
    top_n(20, log2_or)
  
  # 筛选出需要标注的数据点
  labels_data <- or_type_data %>%
    semi_join(top_positive_diseases, by = "Disease_Code") %>%
    group_by(Disease_Code) %>%
    filter(log2_or == max(log2_or)) %>%
    slice(1)
  
  # 创建图形
  p <- ggplot(or_type_data, aes(x = category, y = log2_or, 
                                color = quintile, 
                                shape = significant,
                                group = interaction(quintile, Disease_Code))) +
    geom_point(size = 1.5, alpha = 0.5, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    
    # 为前5个疾病添加标签（标注在左边，黑色字体）
    geom_text(
      data = labels_data,
      aes(label = Disease_Code),
      size = 2.5,
      nudge_x = -0.3,
      hjust = 1,
      show.legend = FALSE,
      check_overlap = TRUE,
      family = "Arial",
      color = "black"  # 添加黑色字体
    ) +
    
    scale_shape_manual(values = c("Significant" = 17, "Not significant" = 16)) +
    scale_color_manual(
      values = c(
        "Q1" = "#1f77b4",
        "Q2" = "#ffcc00", 
        "Q3" = "#2ca02c",
        "Q4" = "#e377c2",
        "Q5" = "#d62728"
      )
    ) +
    labs(
      title = or_type,
      x = "",
      y = expression(bold(log[bold("2")] * "(Odds Ratio)")),
      color = "Quintile :",
      shape = ""
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),
      axis.text = element_text(size = 9, family = "Arial"),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 9, family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 12),
      axis.title.y = element_text(margin = margin(r = 15)), 
      legend.position = "none",
      plot.margin = margin(10, 10, 2, 10)
    )
  
  plot_list_by_or_type[[or_type]] = p
}


setwd("D:/UKB_data/caogao")
# 保存所有图形到PNG文件
for (or_type in names(plot_list_by_or_type)) {
  safe_filename <- gsub("[^a-zA-Z0-9]", "_", or_type)
  png_file <- paste0("or_by_category_", safe_filename, ".png")
  
  png(png_file, width = 2200, height = 900, res = 150)  # 增加宽度以容纳图例
  print(plot_list_by_or_type[[or_type]])
  dev.off()
  
  cat("已保存:", png_file, "\n")
}

plot_list_by_or_type[["WC"]]


# 疾病数量
major_data <- or_result_classified %>%
  filter(!grepl("\\.", Disease_Code) & Disease_Code != "N47") %>%  # 不包含小数点的为大类
  filter(quintile != "Q1")  # 移除Q1

# 为所有OR类型创建图表
or_types <- list(
  c("Odds_ratio_proBMIΔ", "p.adj_proBMIΔ", "Number_of_participants_proBMIΔ", "proBMIΔ"),
  c("Odds_ratio_BMI", "p.adj_BMI", "Number_of_participants_BMI", "BMI"),
  c("Odds_ratio_proBMI", "p.adj_proBMI", "Number_of_participants_proBMI", "proBMI"),
  c("Odds_ratio_proWCΔ", "p.adj_proWCΔ", "Number_of_participants_proWCΔ", "proWCΔ"),
  c("Odds_ratio_WC", "p.adj_WC", "Number_of_participants_WC", "WC"),
  c("Odds_ratio_proWC", "p.adj_proWC", "Number_of_participants_proWC", "proWC")
)

# 修正函数，只计算显著疾病
create_stacked_plot_significant <- function(or_column, p_adj_column, participants_column, title_suffix) {
  # 预处理数据，对长标签进行换行
  stacked_data <- major_data %>%
    mutate(
      category = str_wrap(category, width = 20)
    ) %>%
    group_by(category, quintile) %>%
    summarise(
      sig_or_gt1 = sum(.data[[or_column]] > 1 & .data[[p_adj_column]] < 0.05, na.rm = TRUE),  # 显著且OR>1
      sig_or_lt1 = sum(.data[[or_column]] < 1 & .data[[p_adj_column]] < 0.05, na.rm = TRUE),  # 显著且OR<1
      .groups = 'drop'
    ) %>%
    pivot_longer(
      cols = c(sig_or_lt1, sig_or_gt1),
      names_to = "or_direction",
      values_to = "count"
    ) %>%
    mutate(or_direction = factor(or_direction, 
                                 levels = c("sig_or_lt1", "sig_or_gt1"),
                                 labels = c("OR < 1", "OR > 1")))
  
  # 蓝色和橙色配色方案
  nc_colors <- c("OR < 1" = "#FF7F0E",   # 蓝色
                 "OR > 1" = "#1F77B4")   # 橙色
  
  ggplot(stacked_data, aes(x = quintile, y = count, fill = or_direction)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    labs(
      title = "",
      x = paste("Quintile of", title_suffix),
      y = "Number of diseases",
      fill = ""
    ) +
    scale_fill_manual(values = nc_colors) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black", family = "Arial"),
      axis.text.y = element_text(size = 8, color = "black", family = "Arial"),
      axis.title.x = element_text(face = "bold", size = 10, color = "black", family = "Arial", 
                                  margin = margin(t = 5)),
      axis.title.y = element_text(face = "bold", size = 10, color = "black", family = "Arial"),
      strip.text = element_text(size = 6.6, face = "bold", margin = margin(0, 0, 0, 0), family = "Arial"),
      strip.background = element_blank(),
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 9, family = "Arial"),
      legend.text = element_text(size = 9, family = "Arial"),
      legend.key.size = unit(0.5, "cm"),
      legend.key = element_rect(fill = "white", colour = "grey50"),
      legend.margin = margin(t = -5, b = 0),
      legend.box.margin = margin(t = -10),
      panel.spacing = unit(0.3, "lines"),
      plot.margin = margin(5, 10, 5, 10)
    )
}

create_stacked_plot_significant <- function(or_column, p_adj_column, participants_column, title_suffix) {
  # 预处理数据，对长标签进行换行
  stacked_data <- major_data %>%
    mutate(
      category = str_wrap(category, width = 20)
    ) %>%
    group_by(category, quintile) %>%
    summarise(
      sig_or_gt1 = sum(.data[[or_column]] > 1 & .data[[p_adj_column]] < 0.05, na.rm = TRUE),  # 显著且OR>1
      sig_or_lt1 = sum(.data[[or_column]] < 1 & .data[[p_adj_column]] < 0.05, na.rm = TRUE),  # 显著且OR<1
      .groups = 'drop'
    ) %>%
    pivot_longer(
      cols = c(sig_or_lt1, sig_or_gt1),
      names_to = "or_direction",
      values_to = "count"
    ) %>%
    mutate(or_direction = factor(or_direction, 
                                 levels = c("sig_or_lt1", "sig_or_gt1"),
                                 labels = c("OR < 1", "OR > 1")))
  
  # 蓝色和橙色配色方案
  nc_colors <- c("OR < 1" = "#FF7F0E",   # 橙色
                 "OR > 1" = "#1F77B4")   # 蓝色
  
  ggplot(stacked_data, aes(x = quintile, y = count, fill = or_direction)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    labs(
      title = "",
      x = paste("Quintile of", title_suffix),
      y = "Number of diseases",
      fill = ""
    ) +
    scale_fill_manual(values = nc_colors) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 移除背景颜色设置
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black", family = "Arial"),
      axis.text.y = element_text(size = 8, color = "black", family = "Arial"),
      axis.title.x = element_text(face = "bold", size = 10, color = "black", family = "Arial", 
                                  margin = margin(t = 5)),
      axis.title.y = element_text(face = "bold", size = 10, color = "black", family = "Arial"),
      strip.text = element_text(size = 6.6, face = "bold", margin = margin(0, 0, 0, 0), family = "Arial"),
      strip.background = element_blank(),
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 9, family = "Arial"),
      legend.text = element_text(size = 9, family = "Arial"),
      legend.key.size = unit(0.5, "cm"),
      legend.key = element_rect(fill = "white", colour = "grey50"),
      legend.margin = margin(t = -5, b = 0),
      legend.box.margin = margin(t = -10),
      panel.spacing = unit(0.3, "lines"),
      plot.margin = margin(5, 10, 5, 10)
    )
}
# 创建图表列表
plots <- list()
for (i in seq_along(or_types)) {
  or_type <- or_types[[i]]
  plots[[i]] <- create_stacked_plot_significant(or_type[1], or_type[2], or_type[3], or_type[4])
}

# 查看第六个图（proWC）
plots[[5]]



# 韦恩图
major_data <- or_result_classified %>%
  filter(!grepl("\\.", Disease_Code) & Disease_Code != "N47") %>%  # 不包含小数点的为大类
  filter(quintile != "Q1")  # 移除Q1

# 极小的字体设置
venn_theme_tiny <- list(
  cat.cex = 0.6,           # 极小的类别名称字体
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  cex = 0.4,               # 极小的数字字体
  fontfamily = "Arial",
  lwd = 1.2,
  lty = "solid",
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  formatted = TRUE,
  # 紧凑的边距设置
  margin = 0.03,
  cat.dist = c(0.07, 0.07, 0.07),
  cat.pos = c(-30, 30, 180),
  rotation.degree = 0,
  # 更紧凑的布局
  cat.col = c("black", "black", "black")
)

# 使用极小字体的函数 - 修改为返回图形对象而不直接绘制
create_venn_tiny <- function(quintile_data, quintile_name) {
  # 数据筛选部分
  wc_sig <- quintile_data %>%
    filter(p.adj_WC < 0.05 & Odds_ratio_WC > 1) %>%
    mutate(Indicator = "WC", OR = Odds_ratio_WC) %>%
    select(Disease_Code, Indicator, OR)
  
  prowc_sig <- quintile_data %>%
    filter(p.adj_proWC < 0.05 & Odds_ratio_proWC > 1) %>%
    mutate(Indicator = "proWC", OR = Odds_ratio_proWC) %>%
    select(Disease_Code, Indicator, OR)
  
  prowcΔ_sig <- quintile_data %>%
    filter(p.adj_proWCΔ < 0.05 & Odds_ratio_proWCΔ > 1) %>%
    mutate(Indicator = "proWCΔ", OR = Odds_ratio_proWCΔ) %>%
    select(Disease_Code, Indicator, OR)
  
  venn_data <- list(
    WC = wc_sig$Disease_Code,
    proWC = prowc_sig$Disease_Code,
    proWCΔ = prowcΔ_sig$Disease_Code
  )
  
  # 创建图形对象但不立即绘制
  venn_plot <- venn.diagram(
    x = venn_data,
    category.names = c("WC", "proWC", "proWCΔ"),
    filename = NULL,
    fill = c("#F4B88D", "#87CEEB", "#DDA0DD"),
    alpha = 0.7,
    cat.cex = venn_theme_tiny$cat.cex,
    cat.fontface = venn_theme_tiny$cat.fontface,
    cat.fontfamily = venn_theme_tiny$cat.fontfamily,
    cex = venn_theme_tiny$cex,
    fontfamily = venn_theme_tiny$fontfamily,
    lwd = venn_theme_tiny$lwd,
    lty = venn_theme_tiny$lty,
    print.mode = venn_theme_tiny$print.mode,
    sigdigs = venn_theme_tiny$sigdigs,
    formatted = venn_theme_tiny$formatted,
    margin = venn_theme_tiny$margin,
    cat.dist = venn_theme_tiny$cat.dist,
    cat.pos = venn_theme_tiny$cat.pos,
    rotation.degree = venn_theme_tiny$rotation.degree,
    cat.col = venn_theme_tiny$cat.col
  )
  
  return(venn_plot)
}

# 分别创建并赋值给g1, g3, g5, g7
# Q5分位数 - 赋予g1
wc_data_Q5 <- major_data %>%
  filter(quintile == "Q5") %>%
  select(Disease_Code, 
         Odds_ratio_WC, p.adj_WC,
         Odds_ratio_proWC, p.adj_proWC,
         Odds_ratio_proWCΔ, p.adj_proWCΔ)

g1 <- create_venn_tiny(wc_data_Q5, "Q5")

# Q2分位数 - 赋予g3
wc_data_Q2 <- major_data %>%
  filter(quintile == "Q2") %>%
  select(Disease_Code, 
         Odds_ratio_WC, p.adj_WC,
         Odds_ratio_proWC, p.adj_proWC,
         Odds_ratio_proWCΔ, p.adj_proWCΔ)

g3 <- create_venn_tiny(wc_data_Q2, "Q2")

# Q3分位数 - 赋予g5
wc_data_Q3 <- major_data %>%
  filter(quintile == "Q3") %>%
  select(Disease_Code, 
         Odds_ratio_WC, p.adj_WC,
         Odds_ratio_proWC, p.adj_proWC,
         Odds_ratio_proWCΔ, p.adj_proWCΔ)

g5 <- create_venn_tiny(wc_data_Q3, "Q3")

# Q4分位数 - 赋予g7
wc_data_Q4 <- major_data %>%
  filter(quintile == "Q4") %>%
  select(Disease_Code, 
         Odds_ratio_WC, p.adj_WC,
         Odds_ratio_proWC, p.adj_proWC,
         Odds_ratio_proWCΔ, p.adj_proWCΔ)

g7 <- create_venn_tiny(wc_data_Q4, "Q4")

# 查看g1 (Q5分位数)
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8, x = 0.5, y = 0.5))
grid.draw(g1)
popViewport()

# 查看g3 (Q2分位数)
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8, x = 0.5, y = 0.5))
grid.draw(g3)
popViewport()

# 查看g5 (Q3分位数)
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8, x = 0.5, y = 0.5))
grid.draw(g5)
popViewport()

# 查看g7 (Q4分位数)
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8, x = 0.5, y = 0.5))
grid.draw(g7)
popViewport()