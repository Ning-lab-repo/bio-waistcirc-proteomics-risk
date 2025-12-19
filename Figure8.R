setwd("C:\\Users\\Administrator\\Desktop\\中介分析\\dead_wc")

# WC
data_ols <- read.csv("TableA_WC_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

data_ac <- read.csv("TableC_WC_Mediation_all_cause.csv", header = TRUE, stringsAsFactors = FALSE)
data_ac$ACME_p_bonf <- NULL
apply_bonferroni_correction <- function(df) {
  # 需要校正的p值列
  p_cols <- c("ADE_p", "ACME_p", "PM_p")
  
  # 计算总检验次数
  total_tests <- sum(sapply(df[p_cols], function(x) sum(!is.na(x))))
  cat("总检验次数:", total_tests, "\n")
  
  # 创建新的数据框，保持原始列顺序
  new_df <- df
  
  # 对每个p值列单独应用Bonferroni校正，并在其后添加新列
  for (col in p_cols) {
    n_tests <- sum(!is.na(df[[col]]))
    bonferroni_col <- paste0(col, "_BH")
    
    # 计算校正后的p值
    #bonferroni_p <- p.adjust(df[[col]], method = "bonferroni")
    bonferroni_p <- p.adjust(df[[col]], method = "BH")
    
    # 找到原始列的位置
    col_index <- which(names(new_df) == col)
    
    # 在原始列后插入校正列
    if (col_index < ncol(new_df)) {
      # 如果原始列不在最后，插入到其后
      new_df <- cbind(new_df[, 1:col_index, drop = FALSE], 
                      bonferroni_p,
                      new_df[, (col_index + 1):ncol(new_df), drop = FALSE])
      names(new_df)[col_index + 1] <- bonferroni_col
    } else {
      # 如果原始列在最后，直接添加
      new_df[[bonferroni_col]] <- bonferroni_p
    }
    
    # 统计信息
    orig_sig <- sum(df[[col]] < 0.05, na.rm = TRUE)
    adj_sig <- sum(bonferroni_p < 0.05, na.rm = TRUE)
    
    cat(sprintf("%s: 原始显著 %d/%d -> Bonferroni校正后 %d/%d\n", 
                col, orig_sig, n_tests, adj_sig, n_tests))
  }
  
  return(new_df)
}

# 执行校正
bonferroni_ac <- apply_bonferroni_correction(data_ac)

significant_data1 <- bonferroni_ac[bonferroni_ac$PM_p_BH < 0.05, ]


data_cvd <- read.csv("TableC_WC_Mediation_cvd.csv", header = TRUE, stringsAsFactors = FALSE)
data_cvd$ACME_p_bonf <- NULL
# 执行校正
bonferroni_cvd <- apply_bonferroni_correction(data_cvd)

significant_data2 <- bonferroni_cvd[bonferroni_cvd$PM_p_BH < 0.05, ]


data_cancer <- read.csv("TableC_WC_Mediation_cancer.csv", header = TRUE, stringsAsFactors = FALSE)
data_cancer$ACME_p_bonf <- NULL
# 执行校正
bonferroni_cancer <- apply_bonferroni_correction(data_cancer)

significant_data3 <- bonferroni_cancer[bonferroni_cancer$PM_p_BH < 0.05, ]


data_diabetes <- read.csv("TableC_WC_Mediation_diabetes.csv", header = TRUE, stringsAsFactors = FALSE)
data_diabetes$ACME_p_bonf <- NULL
# 执行校正
bonferroni_diabetes <- apply_bonferroni_correction(data_diabetes)

significant_data4 <- bonferroni_diabetes[bonferroni_diabetes$PM_p_BH < 0.05, ]

# 获取四个数据集的共有蛋白
common_proteins <- Reduce(intersect, list(
  significant_data1$protein,
  significant_data2$protein, 
  significant_data3$protein,
  significant_data4$protein
))

# 创建合并数据集
combined_data <- data.frame(protein = common_proteins)

# 添加每个数据集中对应蛋白的信息
combined_data <- merge(combined_data, significant_data1, by = "protein", all.x = TRUE, suffixes = c("", "_acd"))
combined_data <- merge(combined_data, significant_data2, by = "protein", all.x = TRUE, suffixes = c("", "_cvd"))
combined_data <- merge(combined_data, significant_data3, by = "protein", all.x = TRUE, suffixes = c("", "_cancer"))
combined_data <- merge(combined_data, significant_data4, by = "protein", all.x = TRUE, suffixes = c("", "_diabetes"))
write.csv(combined_data, "WC_common_mediator_proteins.csv")


# 合并四个数据集并计算综合显著性（使用最小的p值）
combined_data <- bind_rows(
  significant_data1 %>% mutate(dataset = "data1"),
  significant_data2 %>% mutate(dataset = "data2"), 
  significant_data3 %>% mutate(dataset = "data3"),
  significant_data4 %>% mutate(dataset = "data4")
) %>%
  filter(protein %in% common_proteins) %>%
  # 首先筛选显著的结果
  filter(ACME_p_BH < 0.05) %>%
  group_by(protein) %>%
  summarise(
    min_ACME_p = min(ACME_p_BH, na.rm = TRUE),
    min_ADE_p = min(ADE_p_BH, na.rm = TRUE),
    # 计算平均PM或最大PM
    mean_PM = mean(PM, na.rm = TRUE),
    max_PM = max(PM, na.rm = TRUE)
  ) %>%
  # 按PM从大到小排序，取前5个
  arrange(desc(mean_PM)) %>%
  head(5)

# 获取最显著的前五个蛋白
top_proteins <- combined_data$protein

graph_list <- list()

# 为每个蛋白创建图形
for(protein in top_proteins) {
  
  # 预先计算所有标签
  labels <- list()
  
  for(i in 1:4) {
    # 获取对应数据集的数据
    data_name <- paste0("significant_data", i)
    protein_data <- get(data_name) %>% filter(protein == !!protein)
    
    # 获取对应的OLS数据
    ols_data <- data_ols %>% filter(protein == !!protein)
    
    # 获取对应的结局名称
    outcomes <- c("all_cause_death", "cvd_death", "cancer_death", "diabetes_death")
    outcome_name <- outcomes[i]
    
    if(nrow(protein_data) > 0 && nrow(ols_data) > 0) {
      # 添加显著性标记函数 - 使用 3×10⁻⁴ 格式
      add_significance_scientific <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        
        # 将beta转换为 3×10⁻⁴ 格式
        if(beta == 0) return("0")
        
        exponent <- floor(log10(abs(beta)))
        coefficient <- beta / (10^exponent)
        
        # 格式化系数，保留1位小数
        coefficient_formatted <- sprintf("%.1f", round(coefficient, 1))
        
        # 处理上标数字
        superscript <- function(x) {
          x <- as.character(x)
          x <- gsub("-", "⁻", x)
          x <- gsub("0", "⁰", x)
          x <- gsub("1", "¹", x)
          x <- gsub("2", "²", x)
          x <- gsub("3", "³", x)
          x <- gsub("4", "⁴", x)
          x <- gsub("5", "⁵", x)
          x <- gsub("6", "⁶", x)
          x <- gsub("7", "⁷", x)
          x <- gsub("8", "⁸", x)
          x <- gsub("9", "⁹", x)
          return(x)
        }
        
        beta_formatted <- paste0(coefficient_formatted, "×10", superscript(exponent))
        
        if(p_value < 0.001) return(paste0(beta_formatted, "***"))
        if(p_value < 0.01) return(paste0(beta_formatted, "**"))
        if(p_value < 0.05) return(paste0(beta_formatted, "*"))
        return(beta_formatted)
      }
      
      # 普通数值的显著性标记函数（用于WC到蛋白质的路径）
      add_significance <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        if(p_value < 0.001) return(sprintf("%.3f***", beta))
        if(p_value < 0.01) return(sprintf("%.3f**", beta))
        if(p_value < 0.05) return(sprintf("%.3f*", beta))
        return(sprintf("%.3f", beta))
      }
      
      # 计算中介比例显示
      prop_mediated <- ifelse(!is.na(protein_data$PM), 
                              ifelse(protein_data$PM * 100 >= 0.1,
                                     sprintf("%.1f%%", protein_data$PM * 100),
                                     sprintf("%.2f%%", protein_data$PM * 100)),
                              "N/A")
      
      labels[[outcome_name]] <- list(
        color = ifelse(protein_data$ACME_p_BH < 0.05, "lightcoral", "lightgray"),
        # 蛋白质到结局的ACME（实线）- 使用 3×10⁻⁴ 格式
        acme_beta = add_significance_scientific(protein_data$ACME_p_BH, protein_data$ACME),
        # 直接效应ADE（虚线）- 使用 3×10⁻⁴ 格式
        ade_beta = add_significance_scientific(protein_data$ADE_p_BH, protein_data$ADE),
        # WC到蛋白质的线性回归系数 - 保持普通数值格式
        wc_to_protein = add_significance(ols_data$P_bonferroni, ols_data$Beta),
        # 中介比例
        prop_mediated = prop_mediated
      )
    } else {
      labels[[outcome_name]] <- list(
        color = "lightgray", 
        acme_beta = "N/A", 
        ade_beta = "N/A", 
        wc_to_protein = "N/A",
        prop_mediated = "N/A"
      )
    }
  }
  
  # 创建图形 - 调整字体、大小和图形尺寸
  graph_code <- paste0("
digraph mediation_", protein, " {
  graph [layout = dot, rankdir = TB, ranksep = 1.0]
  node [shape = ellipse, style = filled, fontname = Arial, fontsize = 14, height = 0.8, width = 1.2, fixedsize = true]
  edge [fontname = Arial, fontsize = 12]
  
  // 主节点
  WC [label = 'WC', fillcolor = '#F4B88D']
  Protein [label = '", protein, "', fillcolor = 'white']
  
  // 结局节点
  AllCause [label = 'All-cause\\nDeath', fillcolor = '", labels[["all_cause_death"]]$color, "']
  CVD [label = 'CVD\\nDeath', fillcolor = '", labels[["cvd_death"]]$color, "']
  Cancer [label = 'Cancer\\nDeath', fillcolor = '", labels[["cancer_death"]]$color, "']
  Diabetes [label = 'Diabetes\\nDeath', fillcolor = '", labels[["diabetes_death"]]$color, "']
  
  // 主要路径：WC → 蛋白质
  WC -> Protein [label = '", labels[["all_cause_death"]]$wc_to_protein, "']
  
  // 蛋白质到结局的实线路径（ACME效应），在中介比例前加<br/>换行
  Protein -> AllCause [label = <", labels[["all_cause_death"]]$acme_beta, "<BR/><B>(", labels[["all_cause_death"]]$prop_mediated, ")</B>>]
  Protein -> CVD [label = <", labels[["cvd_death"]]$acme_beta, "<BR/><B>(", labels[["cvd_death"]]$prop_mediated, ")</B>>]
  Protein -> Cancer [label = <", labels[["cancer_death"]]$acme_beta, "<BR/><B>(", labels[["cancer_death"]]$prop_mediated, ")</B>>]
  Protein -> Diabetes [label = <", labels[["diabetes_death"]]$acme_beta, "<BR/><B>(", labels[["diabetes_death"]]$prop_mediated, ")</B>>]
  
  // 直接效应路径：WC→结局（ADE效应，虚线）
  WC -> AllCause [label = '", labels[["all_cause_death"]]$ade_beta, "', style = dashed, color = gray50]
  WC -> CVD [label = '", labels[["cvd_death"]]$ade_beta, "', style = dashed, color = gray50]
  WC -> Cancer [label = '", labels[["cancer_death"]]$ade_beta, "', style = dashed, color = gray50]
  WC -> Diabetes [label = '", labels[["diabetes_death"]]$ade_beta, "', style = dashed, color = gray50]
}
")
  
  # 将图形保存到列表，并指定更大的图形尺寸
  graph_list[[protein]] <- grViz(graph_code) %>% 
    htmlwidgets::onRender("
      function(el) {
        el.style.width = '100%';
        el.style.height = '600px';
      }
    ")
}

if(length(graph_list) == 5) {
  library(webshot)
  library(magick)
  
  # 创建临时文件
  temp_files <- paste0("temp_", 1:5, ".png")
  
  # 保存每个图形为PNG
  for(i in 1:5) {
    # 保存为HTML临时文件
    html_file <- paste0("temp_", i, ".html")
    htmltools::save_html(htmltools::as.tags(graph_list[[i]]), html_file)
    
    # 使用webshot转换为PNG - 添加背景透明参数
    webshot(html_file, temp_files[i], vwidth = 800, vheight = 600)
    
    # 删除临时HTML文件
    file.remove(html_file)
  }
  
  # 使用magick加载并处理图片（将白色背景转为透明）
  images <- lapply(temp_files, function(file) {
    img <- image_read(file)
    
    # 先转换为RGBA格式
    img <- image_convert(img, "RGBA")
    
    # 将白色背景转为透明
    # 使用较小的fuzz值确保只移除纯白色背景
    img_transparent <- image_transparent(img, color = "white", fuzz = 5)
    
    # 确保背景完全透明
    img_transparent <- image_background(img_transparent, "none")
    
    return(img_transparent)
  })
  
  # 获取图片尺寸
  img_width <- image_info(images[[1]])$width
  img_height <- image_info(images[[1]])$height
  
  # 方法1：使用正确的空白图片尺寸计算
  # 每行总宽度 = 3 * 图片宽度
  total_row_width <- 3 * img_width
  
  # 第二行需要空白图片的宽度 = (总宽度 - 2*图片宽度) / 2
  blank_width <- (total_row_width - 2 * img_width) / 2
  
  # 创建完全透明的空白图片
  blank_image <- image_blank(blank_width, img_height, 
                             color = "none")
  
  # 第一行：3个图片
  row1 <- image_append(c(images[[1]], images[[2]], images[[3]]))
  
  # 第二行：空白 + 图片4 + 图片5 + 空白
  row2 <- image_append(c(blank_image, images[[4]], images[[5]], blank_image))
  
  # 组合行
  combined_image <- image_append(c(row1, row2), stack = TRUE)
  
  # 确保整个组合图片的背景透明
  combined_image <- image_background(combined_image, "none")
  
  # 保存最终图片
  image_write(combined_image, "combined_mediation_plots.png", 
              format = "png", 
              quality = 100)
  
  # 清理临时文件
  file.remove(temp_files)
  
  cat("组合图形已保存为 'combined_mediation_plots.png'\n")
  cat("第二行的两个图片已居中显示，背景已设为透明\n")
  
} else {
  warning("图形数量不是5个")
}








# 中介分析
data_ols <- read.csv("TableA_BioX_Adjusted_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

data_ac <- read.csv("TableC_BioX_Adjusted_Mediation_all_cause.csv", header = TRUE, stringsAsFactors = FALSE)
data_ac$ACME_p_bonf <- NULL

apply_bonferroni_correction <- function(df) {
  # 需要校正的p值列
  p_cols <- c("ADE_p", "ACME_p", "PM_p")
  
  # 计算总检验次数
  total_tests <- sum(sapply(df[p_cols], function(x) sum(!is.na(x))))
  cat("总检验次数:", total_tests, "\n")
  
  # 创建新的数据框，保持原始列顺序
  new_df <- df
  
  # 对每个p值列单独应用Bonferroni校正，并在其后添加新列
  for (col in p_cols) {
    n_tests <- sum(!is.na(df[[col]]))
    bonferroni_col <- paste0(col, "_BH")
    
    # 计算校正后的p值
    bonferroni_p <- p.adjust(df[[col]], method = "BH")
    
    # 找到原始列的位置
    col_index <- which(names(new_df) == col)
    
    # 在原始列后插入校正列
    if (col_index < ncol(new_df)) {
      # 如果原始列不在最后，插入到其后
      new_df <- cbind(new_df[, 1:col_index, drop = FALSE], 
                      bonferroni_p,
                      new_df[, (col_index + 1):ncol(new_df), drop = FALSE])
      names(new_df)[col_index + 1] <- bonferroni_col
    } else {
      # 如果原始列在最后，直接添加
      new_df[[bonferroni_col]] <- bonferroni_p
    }
    
    # 统计信息
    orig_sig <- sum(df[[col]] < 0.05, na.rm = TRUE)
    adj_sig <- sum(bonferroni_p < 0.05, na.rm = TRUE)
    
    cat(sprintf("%s: 原始显著 %d/%d -> Bonferroni校正后 %d/%d\n", 
                col, orig_sig, n_tests, adj_sig, n_tests))
  }
  
  return(new_df)
}

# 执行校正
bonferroni_ac <- apply_bonferroni_correction(data_ac)

significant_data1 <- bonferroni_ac[bonferroni_ac$PM_p_BH < 0.05, ]



data_cvd <- read.csv("TableC_BioX_Adjusted_Mediation_cvd.csv", header = TRUE, stringsAsFactors = FALSE)
data_cvd$ACME_p_bonf <- NULL
bonferroni_cvd <- apply_bonferroni_correction(data_cvd)

significant_data2 <- bonferroni_cvd[bonferroni_cvd$PM_p_BH < 0.05, ]


data_cancer <- read.csv("TableC_BioX_Adjusted_Mediation_cancer.csv", header = TRUE, stringsAsFactors = FALSE)
data_cancer$ACME_p_bonf <- NULL
bonferroni_cancer <- apply_bonferroni_correction(data_cancer)

significant_data3 <- bonferroni_cancer[bonferroni_cancer$PM_p_BH < 0.05, ]


data_diabetes <- read.csv("TableC_BioX_Adjusted_Mediation_diabetes.csv", header = TRUE, stringsAsFactors = FALSE)
data_diabetes$ACME_p_bonf <- NULL
bonferroni_diabetes <- apply_bonferroni_correction(data_diabetes)

significant_data4 <- bonferroni_diabetes[bonferroni_diabetes$PM_p_BH < 0.05, ]



# 获取四个数据集的共有蛋白
common_proteins <- Reduce(intersect, list(
  significant_data1$protein,
  significant_data2$protein, 
  significant_data3$protein,
  significant_data4$protein
))


# 创建合并数据集
combined_data <- data.frame(protein = common_proteins)

# 添加每个数据集中对应蛋白的信息
combined_data <- merge(combined_data, significant_data1, by = "protein", all.x = TRUE, suffixes = c("", "_acd"))
combined_data <- merge(combined_data, significant_data2, by = "protein", all.x = TRUE, suffixes = c("", "_cvd"))
combined_data <- merge(combined_data, significant_data3, by = "protein", all.x = TRUE, suffixes = c("", "_cancer"))
combined_data <- merge(combined_data, significant_data4, by = "protein", all.x = TRUE, suffixes = c("", "_diabetes"))
write.csv(combined_data, "proWC_common_mediator_proteins.csv")

# 合并四个数据集并计算综合显著性（使用最小的p值）
combined_data <- bind_rows(
  significant_data1 %>% mutate(dataset = "data1"),
  significant_data2 %>% mutate(dataset = "data2"), 
  significant_data3 %>% mutate(dataset = "data3"),
  significant_data4 %>% mutate(dataset = "data4")
) %>%
  filter(protein %in% common_proteins) %>%
  # 首先筛选显著的结果（ACME_p < 0.05）
  filter(ACME_p_BH < 0.05) %>%
  group_by(protein) %>%
  summarise(
    min_ACME_p = min(ACME_p_BH, na.rm = TRUE),
    min_ADE_p = min(ADE_p_BH, na.rm = TRUE),
    # 计算平均PM或最大PM
    mean_PM = mean(PM, na.rm = TRUE),
    max_PM = max(PM, na.rm = TRUE)
  ) %>%
  # 按PM从大到小排序，取前5个
  arrange(desc(mean_PM)) %>%
  head(5)

# 获取最显著的前五个蛋白
top_proteins <- combined_data$protein

graph_list <- list()

# 为每个蛋白创建图形
for(protein in top_proteins) {
  
  # 预先计算所有标签
  labels <- list()
  
  for(i in 1:4) {
    # 获取对应数据集的数据
    data_name <- paste0("significant_data", i)
    protein_data <- get(data_name) %>% filter(protein == !!protein)
    
    # 获取对应的OLS数据
    ols_data <- data_ols %>% filter(protein == !!protein)
    
    # 获取对应的结局名称
    outcomes <- c("all_cause_death", "cvd_death", "cancer_death", "diabetes_death")
    outcome_name <- outcomes[i]
    
    if(nrow(protein_data) > 0 && nrow(ols_data) > 0) {
      # 添加显著性标记函数 - 使用 3×10⁻⁴ 格式
      add_significance_scientific <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        
        # 将beta转换为 3×10⁻⁴ 格式
        if(beta == 0) return("0")
        
        exponent <- floor(log10(abs(beta)))
        coefficient <- beta / (10^exponent)
        
        # 格式化系数，保留1位小数
        coefficient_formatted <- sprintf("%.1f", round(coefficient, 1))
        
        # 处理上标数字
        superscript <- function(x) {
          x <- as.character(x)
          x <- gsub("-", "⁻", x)
          x <- gsub("0", "⁰", x)
          x <- gsub("1", "¹", x)
          x <- gsub("2", "²", x)
          x <- gsub("3", "³", x)
          x <- gsub("4", "⁴", x)
          x <- gsub("5", "⁵", x)
          x <- gsub("6", "⁶", x)
          x <- gsub("7", "⁷", x)
          x <- gsub("8", "⁸", x)
          x <- gsub("9", "⁹", x)
          return(x)
        }
        
        beta_formatted <- paste0(coefficient_formatted, "×10", superscript(exponent))
        
        if(p_value < 0.001) return(paste0(beta_formatted, "***"))
        if(p_value < 0.01) return(paste0(beta_formatted, "**"))
        if(p_value < 0.05) return(paste0(beta_formatted, "*"))
        return(beta_formatted)
      }
      
      # 普通数值的显著性标记函数（用于WC到蛋白质的路径）
      add_significance <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        if(p_value < 0.001) return(sprintf("%.3f***", beta))
        if(p_value < 0.01) return(sprintf("%.3f**", beta))
        if(p_value < 0.05) return(sprintf("%.3f*", beta))
        return(sprintf("%.3f", beta))
      }
      
      # 计算中介比例显示
      prop_mediated <- ifelse(!is.na(protein_data$PM), 
                              ifelse(protein_data$PM * 100 >= 0.1,
                                     sprintf("%.1f%%", protein_data$PM * 100),
                                     sprintf("%.2f%%", protein_data$PM * 100)),
                              "N/A")
      
      labels[[outcome_name]] <- list(
        color = ifelse(protein_data$ACME_p_BH < 0.05, "lightcoral", "lightgray"),
        # 蛋白质到结局的ACME（实线）- 使用 3×10⁻⁴ 格式
        acme_beta = add_significance_scientific(protein_data$ACME_p_BH, protein_data$ACME),
        # 直接效应ADE（虚线）- 使用 3×10⁻⁴ 格式
        ade_beta = add_significance_scientific(protein_data$ADE_p_BH, protein_data$ADE),
        # WC到蛋白质的线性回归系数 - 保持普通数值格式
        wc_to_protein = add_significance(ols_data$P_bonferroni, ols_data$Beta),
        # 中介比例
        prop_mediated = prop_mediated
      )
    } else {
      labels[[outcome_name]] <- list(
        color = "lightgray", 
        acme_beta = "N/A", 
        ade_beta = "N/A", 
        wc_to_protein = "N/A",
        prop_mediated = "N/A"
      )
    }
  }
  
  # 创建图形 - 调整字体、大小和图形尺寸
  graph_code <- paste0("
digraph mediation_", protein, " {
  graph [layout = dot, rankdir = TB, ranksep = 1.0]
  node [shape = ellipse, style = filled, fontname = Arial, fontsize = 14, height = 0.8, width = 1.2, fixedsize = true]
  edge [fontname = Arial, fontsize = 12]
  
  // 主节点
  proWC [label = 'proWC', fillcolor = '#87CEEB']
  Protein [label = '", protein, "', fillcolor = 'white']
  
  // 结局节点
  AllCause [label = 'All-cause\\nDeath', fillcolor = '", labels[["all_cause_death"]]$color, "']
  CVD [label = 'CVD\\nDeath', fillcolor = '", labels[["cvd_death"]]$color, "']
  Cancer [label = 'Cancer\\nDeath', fillcolor = '", labels[["cancer_death"]]$color, "']
  Diabetes [label = 'Diabetes\\nDeath', fillcolor = '", labels[["diabetes_death"]]$color, "']
  
  // 主要路径：proWC → 蛋白质
  proWC -> Protein [label = '", labels[["all_cause_death"]]$wc_to_protein, "']
  
  // 蛋白质到结局的实线路径（ACME效应），在中介比例前加<br/>换行
  Protein -> AllCause [label = <", labels[["all_cause_death"]]$acme_beta, "<BR/><B>(", labels[["all_cause_death"]]$prop_mediated, ")</B>>]
  Protein -> CVD [label = <", labels[["cvd_death"]]$acme_beta, "<BR/><B>(", labels[["cvd_death"]]$prop_mediated, ")</B>>]
  Protein -> Cancer [label = <", labels[["cancer_death"]]$acme_beta, "<BR/><B>(", labels[["cancer_death"]]$prop_mediated, ")</B>>]
  Protein -> Diabetes [label = <", labels[["diabetes_death"]]$acme_beta, "<BR/><B>(", labels[["diabetes_death"]]$prop_mediated, ")</B>>]
  
  // 直接效应路径：proWC→结局（ADE效应，虚线）
  proWC -> AllCause [label = '", labels[["all_cause_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWC -> CVD [label = '", labels[["cvd_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWC -> Cancer [label = '", labels[["cancer_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWC -> Diabetes [label = '", labels[["diabetes_death"]]$ade_beta, "', style = dashed, color = gray50]
}
")
  
  # 将图形保存到列表，并指定更大的图形尺寸
  graph_list[[protein]] <- grViz(graph_code) %>% 
    htmlwidgets::onRender("
      function(el) {
        el.style.width = '100%';
        el.style.height = '600px';
      }
    ")
}


if(length(graph_list) == 5) {
  library(webshot)
  library(magick)
  
  # 创建临时文件
  temp_files <- paste0("temp_", 1:5, ".png")
  
  # 保存每个图形为PNG
  for(i in 1:5) {
    # 保存为HTML临时文件
    html_file <- paste0("temp_", i, ".html")
    htmltools::save_html(htmltools::as.tags(graph_list[[i]]), html_file)
    
    # 使用webshot转换为PNG - 添加背景透明参数
    webshot(html_file, temp_files[i], vwidth = 800, vheight = 600)
    
    # 删除临时HTML文件
    file.remove(html_file)
  }
  
  # 使用magick加载并处理图片（将白色背景转为透明）
  images <- lapply(temp_files, function(file) {
    img <- image_read(file)
    
    # 先转换为RGBA格式
    img <- image_convert(img, "RGBA")
    
    # 将白色背景转为透明
    # 使用较小的fuzz值确保只移除纯白色背景
    img_transparent <- image_transparent(img, color = "white", fuzz = 5)
    
    # 确保背景完全透明
    img_transparent <- image_background(img_transparent, "none")
    
    return(img_transparent)
  })
  
  # 获取图片尺寸
  img_width <- image_info(images[[1]])$width
  img_height <- image_info(images[[1]])$height
  
  # 方法1：使用正确的空白图片尺寸计算
  # 每行总宽度 = 3 * 图片宽度
  total_row_width <- 3 * img_width
  
  # 第二行需要空白图片的宽度 = (总宽度 - 2*图片宽度) / 2
  blank_width <- (total_row_width - 2 * img_width) / 2
  
  # 创建完全透明的空白图片
  blank_image <- image_blank(blank_width, img_height, 
                             color = "none")
  
  # 第一行：3个图片
  row1 <- image_append(c(images[[1]], images[[2]], images[[3]]))
  
  # 第二行：空白 + 图片4 + 图片5 + 空白
  row2 <- image_append(c(blank_image, images[[4]], images[[5]], blank_image))
  
  # 组合行
  combined_image <- image_append(c(row1, row2), stack = TRUE)
  
  # 确保整个组合图片的背景透明
  combined_image <- image_background(combined_image, "none")
  
  # 保存最终图片
  image_write(combined_image, "combined_mediation_plots2.png", 
              format = "png", 
              quality = 100)
  
  # 清理临时文件
  file.remove(temp_files)
  
  cat("组合图形已保存为 'combined_mediation_plots.png'\n")
  cat("第二行的两个图片已居中显示，背景已设为透明\n")
  
} else {
  warning("图形数量不是5个")
}




# 中介分析proWCΔ
data_ols <- read.csv("TableA_BioX_Delta_vs_Proteins_OLS.csv", header = TRUE, stringsAsFactors = FALSE)

data_ac <- read.csv("TableC_BioX_Delta_Mediation_all_cause.csv", header = TRUE, stringsAsFactors = FALSE)
data_ac$ACME_p_bonf <- NULL
# 执行校正
bonferroni_ac <- apply_bonferroni_correction(data_ac)

significant_data1 <- bonferroni_ac[bonferroni_ac$PM_p_BH < 0.05, ]



data_cvd <- read.csv("TableC_BioX_Delta_Mediation_cvd.csv", header = TRUE, stringsAsFactors = FALSE)
data_cvd$ACME_p_bonf <- NULL
bonferroni_cvd <- apply_bonferroni_correction(data_cvd)

significant_data2 <- bonferroni_cvd[bonferroni_cvd$PM_p_BH < 0.05, ]

data_cancer <- read.csv("TableC_BioX_Delta_Mediation_cancer.csv", header = TRUE, stringsAsFactors = FALSE)
data_cancer$ACME_p_bonf <- NULL
bonferroni_cancer <- apply_bonferroni_correction(data_cancer)

significant_data3 <- bonferroni_cancer[bonferroni_cancer$PM_p_BH < 0.05, ]


data_diabetes <- read.csv("TableC_BioX_Delta_Mediation_diabetes.csv", header = TRUE, stringsAsFactors = FALSE)
data_diabetes$ACME_p_bonf <- NULL
bonferroni_diabetes <- apply_bonferroni_correction(data_diabetes)

significant_data4 <- bonferroni_diabetes[bonferroni_diabetes$PM_p_BH < 0.05, ]


# 获取四个数据集的共有蛋白
common_proteins <- Reduce(intersect, list(
  significant_data1$protein,
  significant_data2$protein, 
  significant_data3$protein,
  significant_data4$protein
))


# 创建合并数据集
combined_data <- data.frame(protein = common_proteins)

# 添加每个数据集中对应蛋白的信息
combined_data <- merge(combined_data, significant_data1, by = "protein", all.x = TRUE, suffixes = c("", "_acd"))
combined_data <- merge(combined_data, significant_data2, by = "protein", all.x = TRUE, suffixes = c("", "_cvd"))
combined_data <- merge(combined_data, significant_data3, by = "protein", all.x = TRUE, suffixes = c("", "_cancer"))
combined_data <- merge(combined_data, significant_data4, by = "protein", all.x = TRUE, suffixes = c("", "_diabetes"))
write.csv(combined_data, "proWCΔ_common_mediator_proteins.csv")

# 合并四个数据集并计算综合显著性（使用最小的p值）
combined_data <- bind_rows(
  significant_data1 %>% mutate(dataset = "data1"),
  significant_data2 %>% mutate(dataset = "data2"), 
  significant_data3 %>% mutate(dataset = "data3"),
  significant_data4 %>% mutate(dataset = "data4")
) %>%
  filter(protein %in% common_proteins) %>%
  # 首先筛选显著的结果（ACME_p < 0.05）
  filter(ACME_p_BH < 0.05) %>%
  group_by(protein) %>%
  summarise(
    min_ACME_p = min(ACME_p_BH, na.rm = TRUE),
    min_ADE_p = min(ADE_p_BH, na.rm = TRUE),
    # 计算平均PM或最大PM
    mean_PM = mean(PM, na.rm = TRUE),
    max_PM = max(PM, na.rm = TRUE)
  ) %>%
  # 按PM从大到小排序，取前5个
  arrange(desc(mean_PM)) %>%
  head(5)

# 获取最显著的前五个蛋白
top_proteins <- combined_data$protein

graph_list <- list()

# 为每个蛋白创建图形
for(protein in top_proteins) {
  
  # 预先计算所有标签
  labels <- list()
  
  for(i in 1:4) {
    # 获取对应数据集的数据
    data_name <- paste0("significant_data", i)
    protein_data <- get(data_name) %>% filter(protein == !!protein)
    
    # 获取对应的OLS数据
    ols_data <- data_ols %>% filter(protein == !!protein)
    
    # 获取对应的结局名称
    outcomes <- c("all_cause_death", "cvd_death", "cancer_death", "diabetes_death")
    outcome_name <- outcomes[i]
    
    if(nrow(protein_data) > 0 && nrow(ols_data) > 0) {
      # 添加显著性标记函数 - 使用 3×10⁻⁴ 格式
      add_significance_scientific <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        
        # 将beta转换为 3×10⁻⁴ 格式
        if(beta == 0) return("0")
        
        exponent <- floor(log10(abs(beta)))
        coefficient <- beta / (10^exponent)
        
        # 格式化系数，保留1位小数
        coefficient_formatted <- sprintf("%.1f", round(coefficient, 1))
        
        # 处理上标数字
        superscript <- function(x) {
          x <- as.character(x)
          x <- gsub("-", "⁻", x)
          x <- gsub("0", "⁰", x)
          x <- gsub("1", "¹", x)
          x <- gsub("2", "²", x)
          x <- gsub("3", "³", x)
          x <- gsub("4", "⁴", x)
          x <- gsub("5", "⁵", x)
          x <- gsub("6", "⁶", x)
          x <- gsub("7", "⁷", x)
          x <- gsub("8", "⁸", x)
          x <- gsub("9", "⁹", x)
          return(x)
        }
        
        beta_formatted <- paste0(coefficient_formatted, "×10", superscript(exponent))
        
        if(p_value < 0.001) return(paste0(beta_formatted, "***"))
        if(p_value < 0.01) return(paste0(beta_formatted, "**"))
        if(p_value < 0.05) return(paste0(beta_formatted, "*"))
        return(beta_formatted)
      }
      
      # 普通数值的显著性标记函数（用于WC到蛋白质的路径）
      add_significance <- function(p_value, beta) {
        if(is.na(p_value) || is.na(beta)) return("N/A")
        if(p_value < 0.001) return(sprintf("%.3f***", beta))
        if(p_value < 0.01) return(sprintf("%.3f**", beta))
        if(p_value < 0.05) return(sprintf("%.3f*", beta))
        return(sprintf("%.3f", beta))
      }
      
      # 计算中介比例显示
      prop_mediated <- ifelse(!is.na(protein_data$PM), 
                              ifelse(protein_data$PM * 100 >= 0.1,
                                     sprintf("%.1f%%", protein_data$PM * 100),
                                     sprintf("%.2f%%", protein_data$PM * 100)),
                              "N/A")
      
      labels[[outcome_name]] <- list(
        color = ifelse(protein_data$ACME_p_BH < 0.05, "lightcoral", "lightgray"),
        # 蛋白质到结局的ACME（实线）- 使用 3×10⁻⁴ 格式
        acme_beta = add_significance_scientific(protein_data$ACME_p_BH, protein_data$ACME),
        # 直接效应ADE（虚线）- 使用 3×10⁻⁴ 格式
        ade_beta = add_significance_scientific(protein_data$ADE_p_BH, protein_data$ADE),
        # proWC到蛋白质的线性回归系数 - 保持普通数值格式
        wc_to_protein = add_significance(ols_data$P_bonferroni, ols_data$Beta),
        # 中介比例
        prop_mediated = prop_mediated
      )
    } else {
      labels[[outcome_name]] <- list(
        color = "lightgray", 
        acme_beta = "N/A", 
        ade_beta = "N/A", 
        wc_to_protein = "N/A",
        prop_mediated = "N/A"
      )
    }
  }
  
  # 创建图形 - 恢复圆圈长宽设置
  graph_code <- paste0("
digraph mediation_", protein, " {
  graph [layout = dot, rankdir = TB, ranksep = 1.0]
  node [shape = ellipse, style = filled, fontname = Arial, fontsize = 14, height = 0.8, width = 1.2, fixedsize = true]
  edge [fontname = Arial, fontsize = 12]
  
  // 主节点
  proWCΔ [label = 'proWCΔ', fillcolor = '#DDA0DD']
  Protein [label = '", protein, "', fillcolor = 'white']
  
  // 结局节点
  AllCause [label = 'All-cause\\nDeath', fillcolor = '", labels[["all_cause_death"]]$color, "']
  CVD [label = 'CVD\\nDeath', fillcolor = '", labels[["cvd_death"]]$color, "']
  Cancer [label = 'Cancer\\nDeath', fillcolor = '", labels[["cancer_death"]]$color, "']
  Diabetes [label = 'Diabetes\\nDeath', fillcolor = '", labels[["diabetes_death"]]$color, "']
  
  // 主要路径：proWCΔ → 蛋白质
  proWCΔ -> Protein [label = '", labels[["all_cause_death"]]$wc_to_protein, "']
  
  // 蛋白质到结局的实线路径（ACME效应），在中介比例前加<br/>换行
  Protein -> AllCause [label = <", labels[["all_cause_death"]]$acme_beta, "<BR/><B>(", labels[["all_cause_death"]]$prop_mediated, ")</B>>]
  Protein -> CVD [label = <", labels[["cvd_death"]]$acme_beta, "<BR/><B>(", labels[["cvd_death"]]$prop_mediated, ")</B>>]
  Protein -> Cancer [label = <", labels[["cancer_death"]]$acme_beta, "<BR/><B>(", labels[["cancer_death"]]$prop_mediated, ")</B>>]
  Protein -> Diabetes [label = <", labels[["diabetes_death"]]$acme_beta, "<BR/><B>(", labels[["diabetes_death"]]$prop_mediated, ")</B>>]
  
  // 直接效应路径：proWCΔ→结局（ADE效应，虚线）
  proWCΔ -> AllCause [label = '", labels[["all_cause_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWCΔ -> CVD [label = '", labels[["cvd_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWCΔ -> Cancer [label = '", labels[["cancer_death"]]$ade_beta, "', style = dashed, color = gray50]
  proWCΔ -> Diabetes [label = '", labels[["diabetes_death"]]$ade_beta, "', style = dashed, color = gray50]
}
")
  
  # 将图形保存到列表，并指定图形尺寸
  graph_list[[protein]] <- grViz(graph_code) %>% 
    htmlwidgets::onRender("
      function(el) {
        el.style.width = '100%';
        el.style.height = '600px';
      }
    ")
}

if(length(graph_list) == 1) {
  library(webshot)
  library(magick)
  
  # 创建临时文件
  temp_file <- "temp.png"
  html_file <- "temp.html"
  
  # 保存为HTML
  htmltools::save_html(htmltools::as.tags(graph_list[[1]]), html_file)
  
  # 转换为PNG
  webshot(html_file, temp_file, vwidth = 800, vheight = 600)
  file.remove(html_file)
  
  # 处理背景透明
  img <- image_read(temp_file) %>%
    image_convert("RGBA") %>%
    image_transparent("white", fuzz = 5) %>%
    image_background("none")
  
  # 保存最终图片
  image_write(img, "mediation_plot.png", format = "png", quality = 100)
  
  # 清理
  file.remove(temp_file)
  
  cat("图形已保存为 'mediation_plot.png'（背景透明）\n")
  
} else {
  warning(paste("图形数量是", length(graph_list), "个，不是1个"))
}


