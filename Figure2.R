# 加载extrafont包
library(extrafont)

# 导入系统字体（只需要运行一次）
font_import()

# 查看可用的字体
fonts()

# 加载字体到R
loadfonts(device = "win")

# 然后运行您的函数


setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

protein_data <- complete_data_WC[ ,1:2920]
X <- as.matrix(protein_data)  # 蛋白质特征矩阵
y <- complete_data_WC$WC  # WC

# 设置随机种子保证可重复性
set.seed(123)
registerDoParallel(cores = 4)  # 使用4个CPU核心

# 自定义函数计算BioX指标
calculate_biox_metrics <- function(actual, predicted) {
  # 线性拟合实际值与预测值
  fit <- lm(predicted ~ actual)
  line_of_best_fit <- predict(fit, newdata = data.frame(actual = actual))
  
  # 计算BioX指标
  biox_adjusted <- actual + (predicted - line_of_best_fit)
  biox_delta <- predicted - line_of_best_fit
  
  return(list(
    biox = predicted,
    biox_adjusted = biox_adjusted,
    biox_delta = biox_delta,
    slope = coef(fit)[2],
    intercept = coef(fit)[1]
  ))
}

# 10折交叉验证函数 ridge
cv_ridge <- function(X, y, covariates, nfolds = 10) {
  
  # 确保协变量矩阵与X有相同的行数
  if (nrow(covariates) != nrow(X)) {
    stop("协变量矩阵的行数必须与特征矩阵X相同")
  }
  
  # 将协变量矩阵转换为数值矩阵
  covariates_matrix <- as.matrix(covariates)
  
  # 合并特征和协变量
  X_with_covariates <- cbind(X, covariates_matrix)
  
  fold_ids <- createFolds(y, k = nfolds)
  predictions <- rep(NA, length(y))
  best_models <- list()  # 存储每个fold的最佳模型
  
  for (fold in 1:nfolds) {
    # 划分训练集和验证集
    train_idx <- unlist(fold_ids[-fold])
    test_idx <- fold_ids[[fold]]
    
    # 训练Ridge模型（包含协变量）
    cv_model <- cv.glmnet(X_with_covariates[train_idx, ], y[train_idx], alpha = 0, standardize = TRUE)
    best_lambda <- cv_model$lambda.min
    
    # 预测（使用相同的协变量结构）
    predictions[test_idx] <- predict(cv_model, s = best_lambda, newx = X_with_covariates[test_idx, ])
  }
  
  # 计算性能指标
  r2 <- cor(y, predictions)^2
  rmse <- sqrt(mean((y - predictions)^2))
  
  # 计算BioX指标
  biox_metrics <- calculate_biox_metrics(y, predictions)
  
  final_model <- glmnet(X_with_covariates, y, alpha = 0, lambda = best_lambda, standardize = TRUE) #cv.glmnet
  
  return(list(
    predictions = predictions,
    r2 = r2,
    rmse = rmse,
    biox_metrics = biox_metrics,
    final_model = final_model  # 返回最终模型
  ))
}

# 提取年龄和性别协变量
age_sex_covariates <- complete_data_WC[, c("Age", "Sex")]  # 根据实际列名调整

# 运行交叉验证
ridge_results <- cv_ridge(
  X = X,  # 你的特征矩阵
  y = y,  # 你的目标变量
  covariates = age_sex_covariates,
  nfolds = 10
)

# 10折交叉验证函数 lasso
cv_lasso <- function(X, y, covariates, nfolds = 10) {
  
  # 确保协变量矩阵与X有相同的行数
  if (nrow(covariates) != nrow(X)) {
    stop("协变量矩阵的行数必须与特征矩阵X相同")
  }
  
  # 将协变量矩阵转换为数值矩阵
  covariates_matrix <- as.matrix(covariates)
  
  # 合并特征和协变量
  X_with_covariates <- cbind(X, covariates_matrix)
  
  fold_ids <- createFolds(y, k = nfolds)
  predictions <- rep(NA, length(y))
  best_models <- list()  # 存储每个fold的最佳模型
  
  for (fold in 1:nfolds) {
    # 划分训练集和验证集
    train_idx <- unlist(fold_ids[-fold])
    test_idx <- fold_ids[[fold]]
    
    # 训练Ridge模型（包含协变量）
    cv_model <- cv.glmnet(X_with_covariates[train_idx, ], y[train_idx], alpha = 1, standardize = TRUE)
    best_lambda <- cv_model$lambda.min
    
    # 预测（使用相同的协变量结构）
    predictions[test_idx] <- predict(cv_model, s = best_lambda, newx = X_with_covariates[test_idx, ])
  }
  
  # 计算性能指标
  r2 <- cor(y, predictions)^2
  rmse <- sqrt(mean((y - predictions)^2))
  
  # 计算BioX指标
  biox_metrics <- calculate_biox_metrics(y, predictions)
  
  final_model <- glmnet(X_with_covariates, y, alpha = 1, lambda = best_lambda, standardize = TRUE)
  
  return(list(
    predictions = predictions,
    r2 = r2,
    rmse = rmse,
    biox_metrics = biox_metrics,
    final_model = final_model # 返回最终模型
  ))
}

# 运行交叉验证
lasso_results <- cv_lasso(
  X = X,  # 你的特征矩阵
  y = y,  # 你的目标变量
  covariates = age_sex_covariates,
  nfolds = 10
)

# 结果整理
results_t <- data.frame(
  row.names = rownames(protein_data),
  Actual_WC = y,  # WC
  pred_WC = lasso_results$predictions,  # pWC
  BioX_Adjusted = lasso_results$biox_metrics$biox_adjusted,  # proWC
  BioX_Delta = lasso_results$biox_metrics$biox_delta  # proWCΔ
)

#write.csv(results_t,"lasso_WC.csv")


# 散点图
library(ggplot2)
library(ggExtra)

make_marginal_hist <- function(x, y, xlabel, ylabel, title, fit_line = TRUE) {
  # 创建数据框
  df <- data.frame(x = x, y = y)
  
  # 创建基础散点图
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, color = "orange", size = 1.5) +
    labs(x = xlabel, y = ylabel, title = NULL) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),
      axis.text = element_text(size = 9, family = "Arial"),
      plot.title = element_blank()
    )
  
  # 添加拟合线
  if (fit_line) {
    model <- lm(y ~ x, data = df)
    slope <- coef(model)[2]
    intercept <- coef(model)[1]
    r_squared <- summary(model)$r.squared
    
    p <- p +
      geom_smooth(method = "lm", se = FALSE, color = "brown", 
                  linewidth = 0.8, linetype = "dashed") +
      annotate("text", 
               x = min(x, na.rm = TRUE) + 0.05 * (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
               y = max(y, na.rm = TRUE) - 0.05 * (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)),  # 上移标签
               label = sprintf("y = %.2fx + %.2f\nR² = %.3f",
                               slope, intercept, r_squared),
               hjust = 0, vjust = 1, size = 3.5, color = "black",
               fontface = "bold", family = "Arial")
  }
  
  # 添加边缘直方图
  ggMarginal(
    p,
    type = "histogram",
    fill = "lightblue",
    color = "black",
    size = 4,
    margins = "both",
    bins = 20
  )
}

results_t <- read.csv("lasso_WC.csv", header = T,row.names = 1)

# pWC
scatter_hist1 <- make_marginal_hist(
  x = results_t$Actual_WC,
  y = results_t$pred_WC,
  xlabel = "WC",
  ylabel = "pWC",
  title = ""
)
scatter_hist1


# proWC
scatter_hist2 <- make_marginal_hist(
  x = results_t$Actual_WC,
  y = results_t$BioX_Adjusted,
  xlabel = "WC",
  ylabel = "proWC",
  title = ""
)
scatter_hist2


# 加载必要的库
library(glmnet)
library(caret)

# 读取数据
England_data <- read.csv("England_data.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
Wales_data <- read.csv("Wales_data.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
Scotland_data <- read.csv("Scotland_data.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
test_data <- rbind(Wales_data, Scotland_data)


# 准备test测试数据
protein_data_test <- test_data[, 1:2920]
X_test <- as.matrix(protein_data_test)
y_test <- test_data$WC
age_sex_covariates_test <- test_data[, c("Age", "Sex")]


# 修改函数以返回训练好的模型
cv_lasso_train <- function(X, y, covariates, nfolds = 10) {
  
  # 确保协变量矩阵与X有相同的行数
  if (nrow(covariates) != nrow(X)) {
    stop("协变量矩阵的行数必须与特征矩阵X相同")
  }
  
  # 将协变量矩阵转换为数值矩阵
  covariates_matrix <- as.matrix(covariates)
  
  # 合并特征和协变量
  X_with_covariates <- cbind(X, covariates_matrix)
  
  fold_ids <- createFolds(y, k = nfolds)
  predictions <- rep(NA, length(y))
  
  # 使用全部数据训练最终模型
  cv_model <- cv.glmnet(X_with_covariates, y, alpha = 1, standardize = TRUE)
  best_lambda <- cv_model$lambda.min
  final_model <- glmnet(X_with_covariates, y, alpha = 1, lambda = best_lambda, standardize = TRUE)
  
  # 交叉验证预测
  for (fold in 1:nfolds) {
    train_idx <- unlist(fold_ids[-fold])
    test_idx <- fold_ids[[fold]]
    
    cv_model_fold <- cv.glmnet(X_with_covariates[train_idx, ], y[train_idx], alpha = 1, standardize = TRUE)
    best_lambda_fold <- cv_model_fold$lambda.min
    
    predictions[test_idx] <- predict(cv_model_fold, s = best_lambda_fold, newx = X_with_covariates[test_idx, ])
  }
  
  # 计算性能指标
  r2 <- cor(y, predictions)^2
  rmse <- sqrt(mean((y - predictions)^2))
  
  # 计算BioX指标
  biox_metrics <- calculate_biox_metrics(y, predictions)
  
  return(list(
    predictions = predictions,
    r2 = r2,
    rmse = rmse,
    biox_metrics = biox_metrics,
    final_model = final_model,
    cv_model = cv_model,
    feature_names = colnames(X_with_covariates)
  ))
}

# 在England数据上训练模型
lasso_results <- cv_lasso_train(
  X = X_train,
  y = y_train,
  covariates = age_sex_covariates_train,
  nfolds = 10
)

# 测试函数
test_lasso_model <- function(model, X_test, y_test, covariates_test, dataset_name) {
  
  # 准备测试数据（确保与训练数据相同的特征顺序）
  covariates_matrix_test <- as.matrix(covariates_test)
  X_test_with_covariates <- cbind(X_test, covariates_matrix_test)
  
  # 确保测试数据特征与训练数据匹配
  train_features <- model$feature_names
  test_features <- colnames(X_test_with_covariates)
  
  # 检查特征是否匹配，如果不匹配则进行调整
  if (!identical(train_features, test_features)) {
    warning(paste("特征不匹配，正在调整", dataset_name, "数据"))
    # 只保留训练集中存在的特征
    common_features <- intersect(train_features, test_features)
    X_test_with_covariates <- X_test_with_covariates[, common_features, drop = FALSE]
    
    # 添加缺失的特征（设为0）
    missing_features <- setdiff(train_features, test_features)
    if (length(missing_features) > 0) {
      missing_matrix <- matrix(0, nrow = nrow(X_test_with_covariates), ncol = length(missing_features))
      colnames(missing_matrix) <- missing_features
      X_test_with_covariates <- cbind(X_test_with_covariates, missing_matrix)
    }
    
    # 按训练特征顺序重新排列
    X_test_with_covariates <- X_test_with_covariates[, train_features, drop = FALSE]
  }
  
  # 预测
  predictions <- predict(model$final_model, newx = X_test_with_covariates, s = model$final_model$lambda)
  predictions <- as.numeric(predictions)
  
  # 计算性能指标
  r2 <- cor(y_test, predictions)^2
  rmse <- sqrt(mean((y_test - predictions)^2))
  mae <- mean(abs(y_test - predictions))
  
  # 计算BioX指标
  biox_metrics <- calculate_biox_metrics(y_test, predictions)
  
  cat("=== ", dataset_name, "测试结果 ===\n")
  cat("R²:", round(r2, 4), "\n")
  cat("RMSE:", round(rmse, 4), "\n")
  cat("MAE:", round(mae, 4), "\n")
  cat("斜率:", round(biox_metrics$slope, 4), "\n")
  cat("截距:", round(biox_metrics$intercept, 4), "\n\n")
  
  return(list(
    predictions = predictions,
    r2 = r2,
    rmse = rmse,
    mae = mae,
    biox_metrics = biox_metrics,
    actual = y_test
  ))
}

# 在Wales数据上测试
test_results <- test_lasso_model(
  model = lasso_results,
  X_test = X_test,
  y_test = y_test,
  covariates_test = age_sex_covariates_test,
  dataset_name = "test"
)

results_t <- data.frame(
  row.names = rownames(protein_data_test),
  Actual_WC = y_test,  # BMI
  pred_WC = test_results$predictions,  # pBMI
  BioX_Adjusted = test_results$biox_metrics$biox_adjusted,  # mBMI
  BioX_Delta = test_results$biox_metrics$biox_delta  # mBMIΔ
)

#write.csv(results_t,"lasso_test.csv")


results_t <- read.csv("lasso_test.csv", header = T,row.names = 1)
scatter_hist3 <- make_marginal_hist(
  x = results_t$Actual_WC,
  y = results_t$pred_WC,
  xlabel = "WC",
  ylabel = "pWC",
  title = ""
)
scatter_hist3

scatter_hist4 <- make_marginal_hist(
  x = results_t$Actual_WC,
  y = results_t$BioX_Adjusted,
  xlabel = "WC",
  ylabel = "proWC",
  title = ""
)
scatter_hist4


# ------------------------- 蛋白质信号关联分析 -------------------------
complete_data_WC <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df <- cbind(complete_data_WC, results)

analyze_biox_delta_signals <- function(merged_data, target_col = "Actual_WC", protein_cols, covariate_cols = c("Age", "Sex")) {
  # 提取BioX_Delta列名
  biox_delta_col <- paste0("BioX_Adjusted")
  
  # 检查数据是否存在
  if (!biox_delta_col %in% colnames(merged_data)) {
    stop(paste("Column", biox_delta_col, "not found in data"))
  }
  
  # 创建分析用数据框
  analysis_df <- merged_data %>%
    select(all_of(c(target_col, biox_delta_col, protein_cols, covariate_cols))) %>%
    na.omit()
  
  if (nrow(analysis_df) == 0) {
    stop("No complete cases after removing NA values")
  }
  
  # 初始化结果存储
  result <- list(
    protein_names = protein_cols,
    measured_betas = numeric(length(protein_cols)),
    measured_pvalues = numeric(length(protein_cols)),
    delta_betas = numeric(length(protein_cols)),
    delta_pvalues = numeric(length(protein_cols))
  )
  
  # 对每个蛋白质进行关联分析
  for (i in seq_along(protein_cols)) {
    protein <- protein_cols[i]
    
    # 1. 测量值与蛋白质的关联 (Measured_X ~ Protein)
    fit_measured <- lm(as.formula(paste(protein, "~", target_col, "+", paste(covariate_cols, collapse = " + "))), data = analysis_df)
    result$measured_betas[i] <- tidy(fit_measured)$estimate[2]
    result$measured_pvalues[i] <- tidy(fit_measured)$p.value[2]
    
    # 2. BioX_Delta与蛋白质的关联 (BioX_Delta ~ Protein)
    fit_delta <- lm(as.formula(paste(protein, "~", biox_delta_col, "+", paste(covariate_cols, collapse = " + "))), data = analysis_df)
    result$delta_betas[i] <- tidy(fit_delta)$estimate[2]
    result$delta_pvalues[i] <- tidy(fit_delta)$p.value[2]
  }
  
  # 应用多重检验校正 (BH校正)
  measured_p_adj <- p.adjust(result$measured_pvalues, method = "BH")
  delta_p_adj <- p.adjust(result$delta_pvalues, method = "BH")
  
  # 转换为数据框便于分析
  result_df <- data.frame(
    protein = protein_cols,
    measured_beta = result$measured_betas,
    measured_p = result$measured_pvalues,
    measured_p_adj_BH = measured_p_adj,
    delta_beta = result$delta_betas,
    delta_p = result$delta_pvalues,
    delta_p_adj_BH = delta_p_adj
  )
  
  # 添加显著性标记（基于BH校正后的p值）
  result_df <- result_df %>%
    mutate(
      measured_sig_BH = ifelse(measured_p_adj_BH < 0.05, "Significant", "Non-significant"),
      delta_sig_BH = ifelse(delta_p_adj_BH < 0.05, "Significant", "Non-significant"),
      effect_type_BH = case_when(
        measured_p_adj_BH < 0.05 & delta_p_adj_BH < 0.05 ~ "Both",
        measured_p_adj_BH < 0.05 ~ "Measured only",
        delta_p_adj_BH < 0.05 ~ "Delta only",
        TRUE ~ "Non-significant"
      )
    )
  
  return(results = result_df)
}

# 执行分析
analysis_output <- analyze_biox_delta_signals(
  merged_data = temp_df,
  target_col = "Actual_WC",
  protein_cols = colnames(X),  # 替换为实际蛋白质列名
  covariate_cols = c("Age", "Sex")
)
write.csv(analysis_output,"analysis_output.csv")

#Protein Effect Sizes: WC vs proWCΔ
create_biox_delta_plot <- function(result_df, target_name = "WC") {
  # 计算回归统计量
  fit <- lm(delta_beta ~ measured_beta, data = result_df)
  slope <- round(coef(fit)[2], 3)
  intercept <- round(coef(fit)[1], 3)
  r_squared <- round(summary(fit)$r.squared, 3)
  
  ggplot(result_df, aes(x = measured_beta, y = delta_beta)) +
    geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
    
    # 添加统计信息文本 - 统一字号
    annotate("text",
             x = min(result_df$measured_beta, na.rm = TRUE),
             y = max(result_df$delta_beta, na.rm = TRUE),
             label = sprintf("y = %.3fx + %.3f\nR² = %.3f", slope, intercept, r_squared),
             hjust = 0, vjust = 1,
             size = 3.5,  # 统一文本字号，与第一个函数一致
             fontface = "bold",
             color = "black",
             family = "Arial") +  # 添加字体家族
    
    labs(
      x = paste("% difference per unit of", target_name),
      y = paste("% difference per unit of proWC")
    ) +
    scale_x_continuous(
      labels = function(x) round(x * 100),
      expand = expansion(mult = c(0.05, 0.05))) + 
    scale_y_continuous(
      labels = function(x) round(x * 100),
      expand = expansion(mult = c(0.05, 0.05))) + 
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),  # 统一字体家族
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11),  # 统一坐标轴标题字号
      axis.text = element_text(size = 9),  # 统一坐标轴刻度字号
      axis.title.y = element_text(margin = margin(r = 10)), 
      legend.position = "bottom",
      plot.margin = margin(0.7, 0.5, 0.2, 0.5, "cm")
    )
}

# 使用示例
p <- create_biox_delta_plot(analysis_output, "WC")
p


setwd("D:/UKB_data/NC")
complete_data_WC <- read.csv("analysis_data_WC.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
results<- read.csv("lasso_WC.csv", header = T,row.names = 1)
temp_df <- cbind(complete_data_WC, results)

# 修改 plot_mBMI_BMI_correlation 函数，进一步缩小字体
plot_proWC_WC_correlation <- function(data, 
                                      biox_delta_col = "BioX_Delta",
                                      biox_adj_col = "BioX_Adjusted",
                                      bmi_col = "Actual_WC") {
  
  # 检查必要列是否存在
  required_cols <- c(biox_delta_col, biox_adj_col, bmi_col)
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, names(data)), collapse = ", "))
  }
  
  # 数据预处理
  plot_df <- data %>%
    dplyr::select(all_of(required_cols)) %>%
    na.omit() %>%
    mutate(
      biox_quintile = cut(
        .data[[biox_delta_col]],
        breaks = quantile(.data[[biox_delta_col]], probs = seq(0, 1, 0.2), na.rm = TRUE),
        labels = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
        include.lowest = TRUE
      )
    ) %>%
    filter(!is.na(biox_quintile))
  
  # 定义颜色方案
  quintile_colors <- c(
    "Q1" = "#1f77b4",  
    "Q2" = "#ffcc00",  
    "Q3" = "#2ca02c",  
    "Q4" = "#e377c2",  
    "Q5" = "#d62728"   
  )
  
  # 获取Q1和Q5的数据用于绘制边框
  q1_data <- plot_df %>% filter(biox_quintile == "Q1")
  q5_data <- plot_df %>% filter(biox_quintile == "Q5")
  
  # 计算Q1和Q5的凸包(convex hull)用于绘制斜矩形
  find_hull <- function(df) {
    df[chull(df[[bmi_col]], df[[biox_adj_col]]), ]
  }
  q1_hull <- find_hull(q1_data)
  q5_hull <- find_hull(q5_data)
  
  # 绘制散点图
  ggplot(plot_df, aes(x = .data[[bmi_col]], y = .data[[biox_adj_col]])) +
    # 所有散点
    geom_point(
      aes(color = biox_quintile),
      alpha = 0.5,
      size = 1.2
    ) +
    # Q1组的斜虚线框(凸包)
    geom_polygon(
      data = q1_hull,
      aes(x = .data[[bmi_col]], y = .data[[biox_adj_col]]),
      fill = NA,
      color = quintile_colors["Q1"],
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Q5组的斜虚线框(凸包)
    geom_polygon(
      data = q5_hull,
      aes(x = .data[[bmi_col]], y = .data[[biox_adj_col]]),
      fill = NA,
      color = quintile_colors["Q5"],
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Q5标注 - 放在上面
    geom_text(
      data = data.frame(
        x = min(q5_hull[[bmi_col]]), 
        y = max(q5_hull[[biox_adj_col]]),
        label = "Quintile 5\n(proWC > WC)"
      ),
      aes(x = x, y = y, label = label),
      color = quintile_colors["Q5"],
      hjust = -0.1,  # 稍微向左偏移
      vjust = 1,     # 顶部对齐
      size = 3.5,
      fontface = "bold",
      family = "Arial"
    ) +
    # Q1标注 - 放在下面
    geom_text(
      data = data.frame(
        x = max(q1_hull[[bmi_col]]), 
        y = min(q1_hull[[biox_adj_col]]),
        label = "Quintile 1\n(proWC < WC)"
      ),
      aes(x = x, y = y, label = label),
      color = quintile_colors["Q1"],
      hjust = 1.1,   # 稍微向右偏移
      vjust = 0,     # 底部对齐
      size = 3.5,
      fontface = "bold",
      family = "Arial"
    ) +
    # 颜色映射
    scale_color_manual(
      name = "Quintile of proWCΔ",
      values = quintile_colors,
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    ) +
    # 标签和标题
    labs(
      x = "WC",
      y = "proWC",
    ) +
    
    # 主题美化
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),
      axis.text = element_text(size = 9, family = "Arial"),
      axis.title.y = element_text(margin = margin(r = 20)), 
      legend.position = "none"
    )
}
p_s <- plot_proWC_WC_correlation(temp_df)

# 箱线图+小提琴
plot_wc_distribution_violin <- function(data,
                                        biox_delta_col = "BioX_Delta",
                                        bmi_col = "Actual_WC") {
  
  # 检查必要列是否存在
  required_cols <- c(biox_delta_col, bmi_col)
  if (!all(required_cols %in% names(data))) {
    stop("缺少必要的列: ", paste(setdiff(required_cols, names(data)), collapse = ", "))
  }
  
  # 数据预处理
  plot_df <- data %>%
    dplyr::select(all_of(required_cols)) %>%
    na.omit() %>%
    dplyr::mutate(
      biox_quintile = cut(
        .data[[biox_delta_col]],
        breaks = quantile(.data[[biox_delta_col]], probs = seq(0, 1, 0.2), na.rm = TRUE),
        labels = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
        include.lowest = TRUE
      )
    ) %>%
    dplyr::filter(!is.na(biox_quintile))
  
  # 定义颜色方案
  quintile_colors <- c(
    "Q1" = "#1f77b4",
    "Q2" = "#ffcc00", 
    "Q3" = "#2ca02c",
    "Q4" = "#e377c2",
    "Q5" = "#d62728"
  )
  
  # 绘制小提琴图 + 箱形图
  p_wc <- ggplot2::ggplot(plot_df, aes(x = biox_quintile, y = .data[[bmi_col]], fill = biox_quintile)) +
    ggplot2::geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    ggplot2::scale_fill_manual(values = quintile_colors) +
    ggplot2::labs(
      title = "",
      x = "Quintile of proWCΔ",
      y = "WC"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = element_text(family = "Arial"),  # 统一字体家族
      panel.grid = element_blank(),  # 去掉所有网格线
      axis.line = element_line(color = "black", linewidth = 0.8),  # 统一线宽
      axis.ticks = element_line(color = "black", linewidth = 0.8),  # 统一线宽
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),  # 统一坐标轴标题字号
      axis.text = element_text(size = 9, family = "Arial"),  # 统一坐标轴刻度字号
      axis.title.y = element_text(margin = margin(r = 15)),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = "Arial")  # 统一标题字号
    )
  
  return(p_wc)
}

# 运行第一个图的小提琴图版本
p_wc_violin <- plot_wc_distribution_violin(temp_df)
print(p_wc_violin)

plot_adj_wc_distribution_violin <- function(data,
                                            biox_delta_col = "BioX_Delta", 
                                            biox_adj_col = "BioX_Adjusted") {
  
  # 检查必要列是否存在
  required_cols <- c(biox_delta_col, biox_adj_col)
  if (!all(required_cols %in% names(data))) {
    stop("缺少必要的列: ", paste(setdiff(required_cols, names(data)), collapse = ", "))
  }
  
  # 数据预处理
  plot_df <- data %>%
    dplyr::select(all_of(required_cols)) %>%
    na.omit() %>%
    dplyr::mutate(
      biox_quintile = cut(
        .data[[biox_delta_col]],
        breaks = quantile(.data[[biox_delta_col]], probs = seq(0, 1, 0.2), na.rm = TRUE),
        labels = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
        include.lowest = TRUE
      )
    ) %>%
    dplyr::filter(!is.na(biox_quintile))
  
  # 定义颜色方案
  quintile_colors <- c(
    "Q1" = "#1f77b4",
    "Q2" = "#ffcc00",
    "Q3" = "#2ca02c", 
    "Q4" = "#e377c2",
    "Q5" = "#d62728"
  )
  
  # 绘制Adjusted proteomic WC的小提琴图
  p_adj_wc <- ggplot2::ggplot(plot_df, aes(x = biox_quintile, y = .data[[biox_adj_col]], fill = biox_quintile)) +
    ggplot2::geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    ggplot2::scale_fill_manual(values = quintile_colors) +
    ggplot2::labs(
      title = "",
      x = "Quintile of proWCΔ",
      y = "proWC"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = element_text(family = "Arial"),  # 统一字体家族
      panel.grid = element_blank(),  # 去掉所有网格线
      axis.line = element_line(color = "black", linewidth = 0.8),  # 统一线宽
      axis.ticks = element_line(color = "black", linewidth = 0.8),  # 统一线宽
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),  # 统一坐标轴标题字号
      axis.text = element_text(size = 9, family = "Arial"),  # 统一坐标轴刻度字号
      axis.title.y = element_text(margin = margin(r = 15)), 
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = "Arial")  # 统一标题字号
    )
  
  return(p_adj_wc)
}

# 运行第二个图的小提琴图版本
p_adj_wc_violin <- plot_adj_wc_distribution_violin(temp_df)
print(p_adj_wc_violin)

# 修改蛋白箱线图函数，进一步缩小字体
plot_protein_boxplots_free_y <- function(data, biox_delta_col = "BioX_Delta", 
                                         covariates = c("Age", "Sex", "Actual_WC")) {
  
  # 定义代谢特征
  protein_data <- c("LEP", "SSC4D", "IGFBP1", "OXT", "GH1", "FGF21")
  
  # 确保数据是数据框格式
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  # 数据预处理 - 更安全的方式
  required_cols <- c(biox_delta_col, covariates, protein_data)
  
  # 检查所有必需的列是否存在
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # 选择并处理数据
  temp_df <- data[, required_cols, drop = FALSE]
  temp_df <- na.omit(temp_df)
  
  if (nrow(temp_df) < 500) stop("Insufficient data (need ≥500 observations)")
  
  # 创建五分位数分组
  biox_values <- temp_df[[biox_delta_col]]
  quintile_breaks <- quantile(biox_values, probs = seq(0, 1, 0.2), na.rm = TRUE)
  
  temp_df$biox_quintile <- cut(
    biox_values,
    breaks = quintile_breaks,
    labels = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
    include.lowest = TRUE
  )
  
  temp_df <- temp_df[!is.na(temp_df$biox_quintile), ]
  
  # 计算Z分数
  for (trait in protein_data) {
    temp_df[[paste0(trait, "_zscore")]] <- scale(temp_df[[trait]])[, 1]
  }
  
  # 准备绘图数据（长格式）
  zscore_cols <- paste0(protein_data, "_zscore")
  plot_data <- temp_df %>%
    dplyr::select(biox_quintile, dplyr::all_of(zscore_cols)) %>%
    tidyr::pivot_longer(-biox_quintile, 
                        names_to = "trait", 
                        values_to = "zscore") %>%
    dplyr::mutate(trait = gsub("_zscore", "", trait)) %>%
    dplyr::mutate(
      trait = factor(trait, levels = protein_data),
      biox_quintile = factor(biox_quintile, levels = c('Q1','Q2','Q3','Q4','Q5'))
    )
  
  # 计算合适的y轴范围
  y_range <- range(plot_data$zscore, na.rm = TRUE)
  y_limit <- max(abs(y_range)) * 1.1  # 扩大10%的范围
  
  quintile_colors <- c(
    "Q1" = "#1f77b4",
    "Q2" = "#ffcc00",
    "Q3" = "#2ca02c",
    "Q4" = "#e377c2",
    "Q5" = "#d62728"
  )
  
  # 绘制整合图形 - 使用须线样式
  ggplot2::ggplot(plot_data, aes(x = trait, y = zscore)) +
    # 使用stat_boxplot来绘制须线，并设置颜色
    ggplot2::stat_boxplot(
      aes(color = biox_quintile),
      geom = "errorbar",
      width = 0.8,
      position = position_dodge(0.7),
      linewidth = 0.8
    ) +
    # 箱体部分
    ggplot2::geom_boxplot(
      aes(fill = biox_quintile, color = biox_quintile),
      width = 0.6,
      position = position_dodge(0.7),
      outlier.shape = NA,
      linewidth = 0.4,
      alpha = 0.7
    ) +
    ggplot2::stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.2,
      stroke = 0.8,
      aes(fill = biox_quintile),
      color = "black",
      position = position_dodge(0.7)
    ) +
    # 使用动态计算的y轴范围
    ggplot2::scale_y_continuous(limits = c(-y_limit, y_limit)) +
    ggplot2::scale_fill_manual(
      values = quintile_colors,
      name = "Quintile of proWCΔ"
    ) +
    ggplot2::scale_color_manual(
      values = quintile_colors,
      name = "Quintile of proWCΔ"
    ) +
    ggplot2::labs(
      x = "",
      y = "Z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = element_text(family = "Arial"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),
      axis.text = element_text(size = 9, family = "Arial"),
      axis.text.x = element_text(size = 8, family = "Arial"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, family = "Arial"),
      axis.title.y = element_text(margin = margin(r = 10)), 
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.margin = margin(0.7, 0.2, -0.15, 0, "cm")
    ) +
    ggplot2::guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1))
}                                                                         

# 重新生成图形
p_f <- plot_protein_boxplots_free_y(temp_df)
p_f


