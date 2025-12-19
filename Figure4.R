setwd("D:/UKB_data/桑葚图14类疾病")
# -*- coding: utf-8 -*-

# Genuine Sankey diagram with horizontal link ends, fixed-size comparable flows.

# Only draw "significant" connections: p.adj < 0.05 AND OR > 1.

# Do NOT display diseases or quintile nodes with no significant connections.

# Horizontal gaps between columns expanded to 3x.

#

# Outputs: diseaseX.html and diseaseX.png in output_dir.

required_pkgs <- c("htmltools", "webshot")

for (p in required_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")

library(htmltools); library(webshot)

# phantomjs for webshot
if (!webshot::is_phantomjs_installed()) {
  message("phantomjs 未检测到，尝试自动安装（可能需要联网）...")
  tryCatch(webshot::install_phantomjs(), error = function(e) warning("phantomjs 安装失败，请手动安装。"))
}

# -----------------------------
# 配置（按需修改）
# -----------------------------
input_dir <- "D:/UKB_data/桑葚图14类疾病"            # 包含 disease1.csv ... 的目录
output_dir <- "genuine_sankey_v2" # 输出目录
dir.create(output_dir, showWarnings = FALSE)

# OR -> cm 基准，用于计算像素映射
cm_per_OR_base <- 1.0
px_per_cm <- 96 / 2.54   # CSS px per cm at 96dpi
px_per_OR_base <- cm_per_OR_base * px_per_cm

# 缩放因子
height_scale_factor <- 0.5
width_scale_factor  <- 2.0

# 视觉参数
horizontal_margin <- 40
vertical_margin <- 18
title_height <- 40
node_widths <- list(left = 120, center = 120, right = 120)
gap_between_columns_base <- 60
# 列间距是原来的 3 倍
gap_between_columns <- gap_between_columns_base * 21  # 调整为 *3 以匹配请求（原为*18可能过大）

node_gap_min <- 10   # 列内部节点间最小空隙（像素）
node_gap_max <- 50  # 新增：列内部节点间最大空隙（像素）

col_quintile <- c(Q2="#ffcc00", Q3="#2ca02c", Q4="#e377c2", Q5="#d62728")

# 中间节点颜色：Paul Tol的muted调色板（颜色盲友好，适合学术审稿）
center_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", 
                   "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9",
                   "#F8C471", "#82E0AA", "#F1948A", "#85C1E9", "#D7BDE2",
                   "#F9E79F", "#A9DFBF", "#F5B7B1", "#AED6F1", "#E8DAEF")

save_svg_file <- TRUE

# -----------------------------
# 贝塞尔 & ribbon 参数（调整为标准Sankey曲线）
# -----------------------------
# 移除原BEZIER_CTRL_MIN_PX等，使用标准中点控制，确保水平端
# 无需FRACTION，固定为 dx/2

RIBBON_SAMPLES <- 40  # 增加采样点以更平滑

# -----------------------------
# 辅助函数
# -----------------------------
safe_num <- function(x) {
  x2 <- gsub("[^0-9\\.eE\\-]", "", as.character(x))
  as.numeric(ifelse(x2 == "" | is.na(x2), NA, x2))
}

# cubic bezier point and derivative (scalar)
bezier_pt <- function(t, p0, p1, p2, p3) {
  (1-t)^3 * p0 + 3*(1-t)^2 * t * p1 + 3*(1-t) * t^2 * p2 + t^3 * p3
}
bezier_deriv <- function(t, p0, p1, p2, p3) {
  3*(1-t)^2 * (p1 - p0) + 6*(1-t)*t * (p2 - p1) + 3*t^2 * (p3 - p2)
}

# compute a ribbon polygon path string for a cubic bezier center curve with width w_px
make_ribbon_path <- function(start_x, start_y, end_x, end_y, c1x, c1y, c2x, c2y, w_px, samples = RIBBON_SAMPLES) {
  half_w <- w_px / 2
  ts <- seq(0, 1, length.out = samples)
  top_pts <- matrix(NA, nrow = samples, ncol = 2)
  bot_pts <- matrix(NA, nrow = samples, ncol = 2)
  for (k in seq_along(ts)) {
    t <- ts[k]
    x <- bezier_pt(t, start_x, c1x, c2x, end_x)
    y <- bezier_pt(t, start_y, c1y, c2y, end_y)
    dx <- bezier_deriv(t, start_x, c1x, c2x, end_x)
    dy <- bezier_deriv(t, start_y, c1y, c2y, end_y)
    nx <- -dy
    ny <- dx
    norm <- sqrt(nx*nx + ny*ny)
    if (norm < 1e-6) {
      nx <- -(end_y - start_y)
      ny <- (end_x - start_x)
      norm <- sqrt(nx*nx + ny*ny)
      if (norm < 1e-6) { nx <- 0; ny <- 1; norm <- 1 }
    }
    nx <- nx / norm
    ny <- ny / norm
    top_pts[k, ] <- c(x + nx * half_w, y + ny * half_w)
    bot_pts[k, ] <- c(x - nx * half_w, y - ny * half_w)
  }
  fmtpt <- function(pt) sprintf("%f %f", pt[1], pt[2])
  top_str <- paste(sapply(seq_len(samples), function(i) fmtpt(top_pts[i, ])), collapse = " L ")
  bot_str <- paste(sapply(seq(samples,1), function(i) fmtpt(bot_pts[i, ])), collapse = " L ")
  path_d <- paste0("M ", top_str, " L ", bot_str, " Z")
  path_d
}

# 获取文件列表
files <- list.files(input_dir, pattern="^disease\\d+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("input_dir 中未找到 diseaseX.csv 文件，请检查路径。")

# -----------------------------
# 1) 预扫描所有文件，计算"原始所需高度"
# -----------------------------
compute_sums_for_file <- function(fpath) {
  df <- tryCatch(read.csv(fpath, check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) { warning("读取失败: ", fpath, " -> ", e$message); return(NULL) })
  if (is.null(df)) return(NULL)
  for (col in c("Disease_Code","quintile","p.adj_proWC","Odds_ratio_proWC","p.adj_WC","Odds_ratio_WC")) {
    if (! (col %in% colnames(df))) df[[col]] <- NA
  }
  df$Odds_ratio_proWC <- safe_num(df$Odds_ratio_proWC)
  df$Odds_ratio_WC <- safe_num(df$Odds_ratio_WC)
  df$p.adj_proWC <- safe_num(df$p.adj_proWC)
  df$p.adj_WC <- safe_num(df$p.adj_WC)
  df$quintile <- as.character(df$quintile)
  df$Disease_Code <- as.character(df$Disease_Code)
  pro_sig <- df[!is.na(df$p.adj_proWC) & df$p.adj_proWC < 0.05 & !is.na(df$Odds_ratio_proWC) & df$Odds_ratio_proWC > 1, , drop = FALSE]
  wc_sig  <- df[!is.na(df$p.adj_WC)     & df$p.adj_WC < 0.05     & !is.na(df$Odds_ratio_WC)     & df$Odds_ratio_WC > 1, , drop = FALSE]
  quintiles <- c("Q2","Q3","Q4","Q5")
  left_sums <- setNames(numeric(length(quintiles)), quintiles)
  right_sums <- setNames(numeric(length(quintiles)), quintiles)
  for (q in quintiles) {
    left_sums[q] <- sum(pro_sig$Odds_ratio_proWC[pro_sig$quintile == q], na.rm=TRUE)
    right_sums[q] <- sum(wc_sig$Odds_ratio_WC[wc_sig$quintile == q], na.rm=TRUE)
  }
  disease_nodes <- sort(unique(c(pro_sig$Disease_Code, wc_sig$Disease_Code)))
  disease_sums <- setNames(numeric(length(disease_nodes)), disease_nodes)
  if (nrow(pro_sig) > 0) {
    for (i in seq_len(nrow(pro_sig))) {
      r <- pro_sig[i,]
      if (!is.na(r$Disease_Code) && !is.na(r$Odds_ratio_proWC)) disease_sums[r$Disease_Code] <- disease_sums[r$Disease_Code] + r$Odds_ratio_proWC
    }
  }
  
  list(left_sums = left_sums, right_sums = right_sums, disease_sums = disease_sums)
}

all_info <- list()
global_max_OR_sum <- 0
for (f in files) {
  info <- compute_sums_for_file(f)
  if (is.null(info)) next
  all_info[[basename(f)]] <- info
  max_col_sum_OR <- max(sum(info$left_sums), sum(info$disease_sums), sum(info$right_sums), na.rm=TRUE)
  if (is.na(max_col_sum_OR)) max_col_sum_OR <- 0
  global_max_OR_sum <- max(global_max_OR_sum, max_col_sum_OR, na.rm = TRUE)
}
if (global_max_OR_sum <= 0) stop("未在任何文件中检测到有效显著连接（p.adj < 0.05 且 OR > 1）。请检查 CSV 内容与列名。")

orig_required_px <- global_max_OR_sum * px_per_OR_base
px_per_OR_new <- px_per_OR_base * height_scale_factor

base_width <- horizontal_margin*2 + node_widths$left + node_widths$center + node_widths$right + gap_between_columns*2
fixed_width <- ceiling(base_width * width_scale_factor)
fixed_inner_height <- ceiling(orig_required_px * height_scale_factor)
fixed_height <- fixed_inner_height * 1.3 + vertical_margin*2 + title_height

message("全局 OR 最大总和 = ", global_max_OR_sum)
message("原始 px_per_OR = ", round(px_per_OR_base,2), " px")
message("新的 px_per_OR = ", round(px_per_OR_new,2), " px")
message("固定输出宽度 = ", fixed_width, " px; 固定输出高度 = ", fixed_height, " px (内部高度 = ", fixed_inner_height, ")")
message("gap_between_columns 已设置为原来的 3 倍：", gap_between_columns, " px")

# -----------------------------
# 2) 为每个文件布局并绘制 ribbon
# -----------------------------
for (f in files) {
  fname <- basename(f)
  cat("处理:", fname, "...\n")
  info <- all_info[[fname]]
  if (is.null(info)) { cat("  -> 预扫描无信息，跳过。\n"); next }
  df <- tryCatch(read.csv(file.path(input_dir, fname), check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) { warning("读取失败:", fname, e$message); NULL })
  if (is.null(df)) next
  
  # 获取分类标题
  chart_title <- if (!is.null(df$category) && length(df$category) > 0 && !is.na(df$category[1])) {
    as.character(df$category[1])
  } else {
    fname
  }
  
  for (col in c("Disease_Code","quintile","p.adj_proWC","Odds_ratio_proWC","p.adj_WC","Odds_ratio_WC")) {
    if (! (col %in% colnames(df))) df[[col]] <- NA
  }
  df$Odds_ratio_proWC <- safe_num(df$Odds_ratio_proWC)
  df$Odds_ratio_WC <- safe_num(df$Odds_ratio_WC)
  df$p.adj_proWC <- safe_num(df$p.adj_proWC)
  df$p.adj_WC <- safe_num(df$p.adj_WC)
  df$quintile <- as.character(df$quintile)
  df$Disease_Code <- as.character(df$Disease_Code)
  pro_sig <- df[!is.na(df$p.adj_proWC) & df$p.adj_proWC < 0.05 & !is.na(df$Odds_ratio_proWC) & df$Odds_ratio_proWC > 1, , drop = FALSE]
  wc_sig  <- df[!is.na(df$p.adj_WC)     & df$p.adj_WC < 0.05     & !is.na(df$Odds_ratio_WC)     & df$Odds_ratio_WC > 1, , drop = FALSE]
  if (nrow(pro_sig) == 0 && nrow(wc_sig) == 0) { cat("  -> 无显著连接，跳过。\n"); next }
  
  left_sums <- info$left_sums
  right_sums <- info$right_sums
  disease_sums <- info$disease_sums
  quintiles <- c("Q2","Q3","Q4","Q5")
  active_left_q <- quintiles[left_sums > 0]
  active_right_q <- quintiles[right_sums > 0]
  disease_nodes <- names(disease_sums)[disease_sums > 0]
  if (length(disease_nodes) == 0) { cat("  -> 没有显著疾病节点，跳过。\n"); next }
  
  left_h <- pmax(round(left_sums[active_left_q] * px_per_OR_new), 2)
  right_h <- pmax(round(right_sums[active_right_q] * px_per_OR_new), 2)
  disease_order <- names(sort(disease_sums[disease_nodes], decreasing = TRUE))
  center_h <- pmax(round(as.numeric(disease_sums[disease_order]) * px_per_OR_new), 2)
  
  compute_gap_for_column <- function(heights) {
    n <- length(heights)
    if (n == 0) return(list(gap = 0, top = 0))
    total_h <- sum(heights)
    n_gaps <- n + 1
    required_gap_total <- max(0, fixed_inner_height - total_h)
    gap_size <- required_gap_total / n_gaps
    
    if (gap_size < node_gap_min) {
      gap_size <- node_gap_min
      top_offset <- max(0, (fixed_inner_height - (total_h + (n-1)*gap_size)) / 2)
      return(list(gap = gap_size, top = top_offset))
    } else if (gap_size > node_gap_max) {
      gap_size <- node_gap_max
      top_offset <- max(0, (fixed_inner_height - (total_h + (n-1)*gap_size)) / 2)
      return(list(gap = gap_size, top = top_offset))
    } else {
      top_offset <- gap_size
      return(list(gap = gap_size, top = top_offset))
    }
  }
  
  left_gap_info <- compute_gap_for_column(left_h)
  center_gap_info <- compute_gap_for_column(center_h)
  right_gap_info <- compute_gap_for_column(right_h)
  
  center_node_width <- node_widths$left
  
  base_left_x <- horizontal_margin
  base_center_x <- base_left_x + node_widths$left + gap_between_columns
  base_right_x <- base_center_x + node_widths$center + gap_between_columns
  current_base_width <- horizontal_margin*2 + node_widths$left + node_widths$center + node_widths$right + gap_between_columns*2
  extra_width <- fixed_width - current_base_width
  left_x_shift <- floor(extra_width / 2)
  base_left_x <- base_left_x + left_x_shift
  base_center_x <- base_center_x + left_x_shift
  base_right_x <- base_right_x + left_x_shift
  
  # 修改：调整节点区域起始位置，为各层标题留出更多空间
  node_area_start_y <- vertical_margin + 120  # 增加顶部空间给各层标题
  
  left_nodes <- data.frame(name = paste0("proWC_", active_left_q), label = active_left_q, x = base_left_x, y = NA, h = as.numeric(left_h), stringsAsFactors = FALSE)
  cur_y <- node_area_start_y + left_gap_info$top
  for (i in seq_len(nrow(left_nodes))) {
    left_nodes$y[i] <- cur_y
    cur_y <- cur_y + left_nodes$h[i] + left_gap_info$gap
  }
  center_nodes <- data.frame(name = disease_order, label = disease_order, x = base_center_x, y = NA, h = as.numeric(center_h), stringsAsFactors = FALSE)
  cur_y <- node_area_start_y + center_gap_info$top
  for (i in seq_len(nrow(center_nodes))) {
    center_nodes$y[i] <- cur_y
    cur_y <- cur_y + center_nodes$h[i] + center_gap_info$gap
  }
  right_nodes <- data.frame(name = paste0("WC_", active_right_q), label = active_right_q, x = base_right_x, y = NA, h = as.numeric(right_h), stringsAsFactors = FALSE)
  cur_y <- node_area_start_y + right_gap_info$top
  for (i in seq_len(nrow(right_nodes))) {
    right_nodes$y[i] <- cur_y
    cur_y <- cur_y + right_nodes$h[i] + right_gap_info$gap
  }
  
  links <- list()
  if (nrow(pro_sig) > 0) {
    for (i in seq_len(nrow(pro_sig))) {
      r <- pro_sig[i,]
      if (is.na(r$quintile) || is.na(r$Odds_ratio_proWC) || is.na(r$Disease_Code)) next
      src <- paste0("proWC_", r$quintile)
      if (! (src %in% left_nodes$name)) next
      tgt <- r$Disease_Code
      if (! (tgt %in% center_nodes$name)) next
      w <- as.numeric(r$Odds_ratio_proWC)
      links[[length(links)+1]] <- list(source = src, target = tgt, weight = w, quintile = r$quintile, type = "pro")
    }
  }
  if (nrow(wc_sig) > 0) {
    for (i in seq_len(nrow(wc_sig))) {
      r <- wc_sig[i,]
      if (is.na(r$quintile) || is.na(r$Odds_ratio_WC) || is.na(r$Disease_Code)) next
      src <- r$Disease_Code
      tgt <- paste0("WC_", r$quintile)
      if (! (src %in% center_nodes$name)) next
      if (! (tgt %in% right_nodes$name)) next
      w <- as.numeric(r$Odds_ratio_WC)
      links[[length(links)+1]] <- list(source = src, target = tgt, weight = w, quintile = r$quintile, type = "wc")
    }
  }
  if (length(links) == 0) { cat("  -> 无显著连接（过滤后），跳过。\n"); next }
  
  node_info_lookup <- function(name) {
    if (startsWith(name, "proWC_")) {
      row <- left_nodes[left_nodes$name == name,]; list(xleft = row$x, xright = row$x + node_widths$left, ytop = row$y, h = row$h)
    } else if (startsWith(name, "WC_")) {
      row <- right_nodes[right_nodes$name == name,]; list(xleft = row$x, xright = row$x + node_widths$right, ytop = row$y, h = row$h)
    } else {
      row <- center_nodes[center_nodes$name == name,]; list(xleft = row$x, xright = row$x + node_widths$center, ytop = row$y, h = row$h)
    }
  }
  # 为proWC和WC分别创建flow_offset
  pro_flow_offset <- list()
  wc_flow_offset <- list()
  
  # 初始化flow_offset
  for (nm in c(left_nodes$name, center_nodes$name, right_nodes$name)) {
    pro_flow_offset[[nm]] <- 0
    wc_flow_offset[[nm]] <- 0
  }
  
  for (i in seq_along(links)) {
    L <- links[[i]]
    w_px <- max(round(L$weight * px_per_OR_new), 1)
    Sinfo <- node_info_lookup(L$source)
    Tinfo <- node_info_lookup(L$target)
    
    if (L$type == "pro") {
      # 使用pro_flow_offset计算proWC连接的位置
      start_y <- Sinfo$ytop + pro_flow_offset[[L$source]] + w_px/2
      end_y <- Tinfo$ytop + pro_flow_offset[[L$target]] + w_px/2
      
      # 只更新pro_flow_offset
      pro_flow_offset[[L$source]] <- pro_flow_offset[[L$source]] + w_px
      pro_flow_offset[[L$target]] <- pro_flow_offset[[L$target]] + w_px
    } else if (L$type == "wc") {
      # 使用wc_flow_offset计算WC连接的位置
      start_y <- Sinfo$ytop + wc_flow_offset[[L$source]] + w_px/2
      end_y <- Tinfo$ytop + wc_flow_offset[[L$target]] + w_px/2
      
      # 只更新wc_flow_offset
      wc_flow_offset[[L$source]] <- wc_flow_offset[[L$source]] + w_px
      wc_flow_offset[[L$target]] <- wc_flow_offset[[L$target]] + w_px
    }
    
    links[[i]]$w_px <- w_px
    links[[i]]$start_x <- Sinfo$xright
    links[[i]]$start_y <- start_y
    links[[i]]$end_x <- Tinfo$xleft
    links[[i]]$end_y <- end_y
  }
  
  svg_elems <- list()
  svg_width <- fixed_width
  svg_height <- fixed_height
  
  # 完全移除背景矩形，确保透明背景
  # svg_elems[[length(svg_elems)+1]] <- tags$rect(x=0, y=0, width=svg_width, height=svg_height, style="fill:none; stroke:none;")
  
  # 修改：主标题放在顶部中央
  # 修改：将主标题放在第一个节点上方30像素
  center_title_x <- base_center_x + node_widths$center / 2
  
  # 找到第一个节点的Y坐标（最小的Y值）
  first_node_y <- min(
    if(nrow(left_nodes) > 0) min(left_nodes$y) else Inf,
    if(nrow(center_nodes) > 0) min(center_nodes$y) else Inf,
    if(nrow(right_nodes) > 0) min(right_nodes$y) else Inf
  )
  
  # 如果找到了有效的节点位置，则在第一个节点上方30像素放置标题
  if (is.finite(first_node_y)) {
    title_y <- first_node_y - 30
  } else {
    # 备用方案：使用节点区域起始位置
    title_y <- node_area_start_y - 30
  }
  
  svg_elems[[length(svg_elems)+1]] <- tags$text(x = center_title_x, y = title_y, chart_title, 
                                                style = "font-weight:bold; font-size:90px; text-anchor:middle;")
  
  max_node_bottom <- max(
    if(nrow(left_nodes) > 0) max(left_nodes$y + left_nodes$h) else 0,
    if(nrow(center_nodes) > 0) max(center_nodes$y + center_nodes$h) else 0,
    if(nrow(right_nodes) > 0) max(right_nodes$y + right_nodes$h) else 0
  )
  
  title_y <- max_node_bottom + 160  # 在最后一个节点下方25像素
  
  # 左列标题 - proWC
  left_title_x <- base_left_x + node_widths$left / 2
  svg_elems[[length(svg_elems)+1]] <- tags$text(x = left_title_x, y = title_y, "proWC", 
                                                style = "font-weight:bold; font-size:70px; text-anchor:middle; fill:#333;")
  
  # 中间列标题 - Disease
  center_title_x <- base_center_x + center_node_width / 2
  svg_elems[[length(svg_elems)+1]] <- tags$text(x = center_title_x, y = title_y, "Disease", 
                                                style = "font-weight:bold; font-size:70px; text-anchor:middle; fill:#333;")
  
  # 右列标题 - WC
  right_title_x <- base_right_x + node_widths$right / 2
  svg_elems[[length(svg_elems)+1]] <- tags$text(x = right_title_x, y = title_y, "WC", 
                                                style = "font-weight:bold; font-size:70px; text-anchor:middle; fill:#333;")
  
  quintile_order <- c("Q2" = 1, "Q3" = 2, "Q4" = 3, "Q5" = 4)
  link_order_idx <- order(sapply(links, function(x) {
    quintile_order[[x$quintile]]
  }), decreasing = TRUE)
  
  for (idx in link_order_idx) {
    L <- links[[idx]]
    dx <- as.numeric(L$end_x - L$start_x)
    ctrl_off <- dx / 2
    c1x <- L$start_x + ctrl_off
    c1y <- L$start_y
    c2x <- L$start_x + ctrl_off
    c2y <- L$end_y
    path_d <- make_ribbon_path(L$start_x, L$start_y, L$end_x, L$end_y,
                               c1x, c1y, c2x, c2y, L$w_px, samples = RIBBON_SAMPLES)
    col <- col_quintile[[L$quintile]]
    if (is.null(col) || is.na(col)) col <- "#777777"
    svg_elems[[length(svg_elems)+1]] <- tags$path(d = path_d, style = paste0("fill:", col, "; fill-opacity:0.75; stroke:none;"))
  }
  
  for (i in seq_len(nrow(left_nodes))) {
    nd <- left_nodes[i,]
    svg_elems[[length(svg_elems)+1]] <- tags$rect(x = nd$x, y = nd$y, width = node_widths$left, height = nd$h, style = paste0("fill:", col_quintile[[nd$label]], "; stroke:#222; stroke-width:0.6;"))
    svg_elems[[length(svg_elems)+1]] <- tags$text(x = nd$x + node_widths$left + 6, y = nd$y + nd$h/2, nd$label, style="font-size:58px; fill:#000; dominant-baseline:middle;")
  }
  for (i in seq_len(nrow(center_nodes))) {
    nd <- center_nodes[i,]
    col_idx <- (i - 1) %% length(center_colors) + 1
    center_col <- center_colors[col_idx]
    svg_elems[[length(svg_elems)+1]] <- tags$rect(x = nd$x, y = nd$y, width = node_widths$center, height = nd$h, style = paste0("fill:", center_col, "; stroke:#333; stroke-width:0.4;"))
    svg_elems[[length(svg_elems)+1]] <- tags$text(x = nd$x + node_widths$center + 6, y = nd$y + nd$h/2, nd$label, style="font-size:58px; fill:#000; dominant-baseline:middle;")
  }
  for (i in seq_len(nrow(right_nodes))) {
    nd <- right_nodes[i,]
    svg_elems[[length(svg_elems)+1]] <- tags$rect(x = nd$x, y = nd$y, width = node_widths$right, height = nd$h, style = paste0("fill:", col_quintile[[nd$label]], "; stroke:#222; stroke-width:0.6;"))
    svg_elems[[length(svg_elems)+1]] <- tags$text(x = nd$x + node_widths$right + 6, y = nd$y + nd$h/2, nd$label, style="font-size:58px; fill:#000; dominant-baseline:middle;")
  }
  
  # 修改：确保SVG背景完全透明
  svg_tag <- tags$svg(xmlns="http://www.w3.org/2000/svg", width = svg_width, height = svg_height, style="background:transparent;", svg_elems)
  
  # 修改：HTML页面也设置为透明背景
  page <- tags$html(
    tags$head(
      tags$meta(charset="utf-8"),
      tags$style("body { background: transparent !important; margin: 0; padding: 0; }")
    ), 
    tags$body(svg_tag, style = "background: transparent; margin: 0; padding: 0;")
  )
  
  html_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(fname), ".html"))
  htmltools::save_html(page, file = html_file)
  cat("  -> 已保存 HTML:", html_file, "\n")
  if (save_svg_file) {
    svg_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(fname), ".svg"))
    writeLines(as.character(svg_tag), con = svg_file, useBytes = TRUE)
    cat("  -> 已保存 SVG:", svg_file, "\n")
  }
  png_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(fname), ".png"))
  tryCatch({
    # 修改：webshot使用透明背景
    webshot::webshot(
      url = html_file, 
      file = png_file, 
      vwidth = svg_width, 
      vheight = svg_height, 
      delay = 0.2,
      # 添加透明背景参数
      selector = "body",
      zoom = 1
    )
    cat("  -> 已保存 PNG:", png_file, "\n")
  }, error = function(e) {
    warning("  -> webshot 保存 PNG 失败: ", e$message, "（HTML 已保存：", html_file, "）")
  })
  cat("完成:", fname, "\n\n")
}

cat("全部完成。输出目录：", normalizePath(output_dir), "\n")


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(forestplot)
library(purrr)

# 各个和各类疾病的下的OR图  
or_result <- read.csv("hbh_or_results.csv", header = T, fileEncoding = "GBK")

top_diseases <- c("G44", "G83")


# 筛选这些疾病的所有Q1-Q5数据
top_diseases_data <- or_result %>%
  filter(Disease_Code %in% top_diseases)

# 分别处理每种指标的数据
# 分别处理每种指标的数据
process_metric_data <- function(data, metric) {
  # 检查所需的列是否存在
  required_cols <- c(paste0("Odds_ratio_", metric), 
                     paste0("conf.low_", metric), 
                     paste0("conf.high_", metric))
  
  # 如果缺少任何必需的列，返回空数据框
  if (!all(required_cols %in% names(data))) {
    warning(paste("Columns for metric", metric, "not found in data"))
    return(data.frame())
  }
  
  metric_data <- data %>%
    dplyr::select(Disease_Code, quintile, 
                  estimate = paste0("Odds_ratio_", metric), 
                  conf.low = paste0("conf.low_", metric), 
                  conf.high = paste0("conf.high_", metric)) %>%
    dplyr::mutate(metric = case_when(
      metric == "proBMIΔ" ~ "proBMIΔ",
      metric == "BMI" ~ "BMI",
      metric == "proBMI" ~ "proBMI",
      metric == "proWCΔ" ~ "proWCΔ",  # 修正这里
      metric == "WC" ~ "WC",
      metric == "proWC" ~ "proWC",
      TRUE ~ metric
    ))
  return(metric_data)
}

# 处理所有指标
# 筛选只保留 WC、proWC、proWCΔ 三个指标
selected_metrics <- c("WC", "proWC", "proWCΔ")

# 处理选定的指标
selected_metric_data <- purrr::map_dfr(selected_metrics, ~ process_metric_data(top_diseases_data, .x))

# 设置因子顺序 - 按 WC、proWC、proWCΔ 顺序
selected_metric_data$metric <- factor(selected_metric_data$metric, 
                                      levels = c("WC", "proWC", "proWCΔ"))  # 修改顺序
selected_metric_data$quintile <- factor(selected_metric_data$quintile, 
                                        levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# 为 Q1-Q5 分配颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",  # 蓝色
  "Q2" = "#ffcc00",  # 黄色
  "Q3" = "#2ca02c",  # 绿色
  "Q4" = "#e377c2",  # 粉色
  "Q5" = "#d62728"   # 红色
)

# 对数据进行log2转换
selected_metric_data_log2 <- selected_metric_data %>%
  mutate(
    estimate_log2 = log2(estimate),
    conf.low_log2 = log2(conf.low),
    conf.high_log2 = log2(conf.high)
  )

# 绘制图形
OR_plot1 <- ggplot(selected_metric_data_log2, aes(
  x = estimate_log2,  # 使用log2转换后的OR值
  y = quintile
)) +
  geom_vline(
    xintercept = 0,  # log2(1) = 0
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low_log2,  # 使用log2转换后的置信区间下限
      xmax = conf.high_log2, # 使用log2转换后的置信区间上限
      color = quintile
    ),
    height = 0.3,  # 与OR_death一致：误差线高度
    linewidth = 1,  # 与OR_death一致：线条粗细
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    aes(
      fill = quintile,
      color = quintile
    ),
    size = 4,       # 与OR_death一致：圆圈大小
    shape = 21,     # 与OR_death一致：带边框的圆形
    stroke = 1,     # 与OR_death一致：边框粗细
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(values = quintile_colors) +
  scale_fill_manual(values = quintile_colors) +
  scale_x_continuous(
    labels = function(x) format(round(x, 1), nsmall = 1)  # 与OR_death一致：统一保留1位小数
  ) +
  facet_grid(metric ~ Disease_Code, scales = "free_x") +
  labs(
    x = "Log2(Odds ratio)",  # 与OR_death一致：坐标轴标签
    y = "Quintile",
    title = "",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_minimal() +  # 移除base_size，与OR_death一致
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 9),  # 与OR_death一致：字体大小
    axis.text = element_text(color = "black", size = 7), # 与OR_death一致：字体大小
    axis.text.y = element_text(size = 8),  # 与OR_death一致：Y轴文字稍大
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 9),  # 与OR_death一致：分面标签大小
    strip.text.x = element_text(size = 9, face = "bold"), # 与OR_death一致
    strip.text.y = element_text(size = 9, face = "bold"), # 与OR_death一致
    legend.position = "none",  # 与OR_death一致：隐藏图例
    plot.margin = margin(-20, 0, 1, 0),  # 与OR_death一致：边距设置
  )

# 显示图形
print(OR_plot1)


# 筛选前5种疾病（基于Q5时mWCΔ的OR值，筛选OR值在4左右的疾病）
top_diseases <- c("I07", "I05","I89","I87","I87","I60")


# 筛选这些疾病的所有Q1-Q5数据
top_diseases_data <- or_result %>%
  filter(Disease_Code %in% top_diseases)

# 分别处理每种指标的数据
# 分别处理每种指标的数据
process_metric_data <- function(data, metric) {
  # 检查所需的列是否存在
  required_cols <- c(paste0("Odds_ratio_", metric), 
                     paste0("conf.low_", metric), 
                     paste0("conf.high_", metric))
  
  # 如果缺少任何必需的列，返回空数据框
  if (!all(required_cols %in% names(data))) {
    warning(paste("Columns for metric", metric, "not found in data"))
    return(data.frame())
  }
  
  metric_data <- data %>%
    dplyr::select(Disease_Code, quintile, 
                  estimate = paste0("Odds_ratio_", metric), 
                  conf.low = paste0("conf.low_", metric), 
                  conf.high = paste0("conf.high_", metric)) %>%
    dplyr::mutate(metric = case_when(
      metric == "proBMIΔ" ~ "proBMIΔ",
      metric == "BMI" ~ "BMI",
      metric == "proBMI" ~ "proBMI",
      metric == "proWCΔ" ~ "proWCΔ",  # 修正这里
      metric == "WC" ~ "WC",
      metric == "proWC" ~ "proWC",
      TRUE ~ metric
    ))
  return(metric_data)
}

# 处理所有指标
# 筛选只保留 WC、proWC、proWCΔ 三个指标
selected_metrics <- c("WC", "proWC", "proWCΔ")

# 处理选定的指标
selected_metric_data <- purrr::map_dfr(selected_metrics, ~ process_metric_data(top_diseases_data, .x))

# 设置因子顺序 - 按 WC、proWC、proWCΔ 顺序
selected_metric_data$metric <- factor(selected_metric_data$metric, 
                                      levels = c("WC", "proWC", "proWCΔ"))  # 修改顺序
selected_metric_data$quintile <- factor(selected_metric_data$quintile, 
                                        levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# 为 Q1-Q5 分配颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",  # 蓝色
  "Q2" = "#ffcc00",  # 黄色
  "Q3" = "#2ca02c",  # 绿色
  "Q4" = "#e377c2",  # 粉色
  "Q5" = "#d62728"   # 红色
)

# 对数据进行log2转换
selected_metric_data_log2 <- selected_metric_data %>%
  mutate(
    estimate_log2 = log2(estimate),
    conf.low_log2 = log2(conf.low),
    conf.high_log2 = log2(conf.high)
  )

# 绘制图形
OR_plot2 <- ggplot(selected_metric_data_log2, aes(
  x = estimate_log2,  # 使用log2转换后的OR值
  y = quintile
)) +
  geom_vline(
    xintercept = 0,  # log2(1) = 0
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low_log2,  # 使用log2转换后的置信区间下限
      xmax = conf.high_log2, # 使用log2转换后的置信区间上限
      color = quintile
    ),
    height = 0.3,  # 与OR_death一致：误差线高度
    linewidth = 1,  # 与OR_death一致：线条粗细
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    aes(
      fill = quintile,
      color = quintile
    ),
    size = 4,       # 与OR_death一致：圆圈大小
    shape = 21,     # 与OR_death一致：带边框的圆形
    stroke = 1,     # 与OR_death一致：边框粗细
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(values = quintile_colors) +
  scale_fill_manual(values = quintile_colors) +
  scale_x_continuous(
    labels = function(x) format(round(x, 1), nsmall = 1)  # 与OR_death一致：统一保留1位小数
  ) +
  facet_grid(metric ~ Disease_Code, scales = "free_x") +
  labs(
    x = "Log2(Odds ratio)",  # 与OR_death一致：坐标轴标签
    y = "Quintile",
    title = "",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_minimal() +  # 移除base_size，与OR_death一致
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 9),  # 与OR_death一致：字体大小
    axis.text = element_text(color = "black", size = 7), # 与OR_death一致：字体大小
    axis.text.y = element_text(size = 8),  # 与OR_death一致：Y轴文字稍大
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 9),  # 与OR_death一致：分面标签大小
    strip.text.x = element_text(size = 9, face = "bold"), # 与OR_death一致
    strip.text.y = element_text(size = 9, face = "bold"), # 与OR_death一致
    legend.position = "none",  # 与OR_death一致：隐藏图例
    plot.margin = margin(-20, 0, 1, 0),  # 与OR_death一致：边距设置
  )

# 显示图形
print(OR_plot2)

top_diseases <- c("J39")


# 筛选这些疾病的所有Q1-Q5数据
top_diseases_data <- or_result %>%
  filter(Disease_Code %in% top_diseases)

# 分别处理每种指标的数据
# 分别处理每种指标的数据
process_metric_data <- function(data, metric) {
  # 检查所需的列是否存在
  required_cols <- c(paste0("Odds_ratio_", metric), 
                     paste0("conf.low_", metric), 
                     paste0("conf.high_", metric))
  
  # 如果缺少任何必需的列，返回空数据框
  if (!all(required_cols %in% names(data))) {
    warning(paste("Columns for metric", metric, "not found in data"))
    return(data.frame())
  }
  
  metric_data <- data %>%
    dplyr::select(Disease_Code, quintile, 
                  estimate = paste0("Odds_ratio_", metric), 
                  conf.low = paste0("conf.low_", metric), 
                  conf.high = paste0("conf.high_", metric)) %>%
    dplyr::mutate(metric = case_when(
      metric == "proBMIΔ" ~ "proBMIΔ",
      metric == "BMI" ~ "BMI",
      metric == "proBMI" ~ "proBMI",
      metric == "proWCΔ" ~ "proWCΔ",  # 修正这里
      metric == "WC" ~ "WC",
      metric == "proWC" ~ "proWC",
      TRUE ~ metric
    ))
  return(metric_data)
}

# 处理所有指标
# 筛选只保留 WC、proWC、proWCΔ 三个指标
selected_metrics <- c("WC", "proWC", "proWCΔ")

# 处理选定的指标
selected_metric_data <- purrr::map_dfr(selected_metrics, ~ process_metric_data(top_diseases_data, .x))

# 设置因子顺序 - 按 WC、proWC、proWCΔ 顺序
selected_metric_data$metric <- factor(selected_metric_data$metric, 
                                      levels = c("WC", "proWC", "proWCΔ"))  # 修改顺序
selected_metric_data$quintile <- factor(selected_metric_data$quintile, 
                                        levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# 为 Q1-Q5 分配颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",  # 蓝色
  "Q2" = "#ffcc00",  # 黄色
  "Q3" = "#2ca02c",  # 绿色
  "Q4" = "#e377c2",  # 粉色
  "Q5" = "#d62728"   # 红色
)

# 对数据进行log2转换
selected_metric_data_log2 <- selected_metric_data %>%
  mutate(
    estimate_log2 = log2(estimate),
    conf.low_log2 = log2(conf.low),
    conf.high_log2 = log2(conf.high)
  )

# 绘制图形
OR_plot3 <- ggplot(selected_metric_data_log2, aes(
  x = estimate_log2,  # 使用log2转换后的OR值
  y = quintile
)) +
  geom_vline(
    xintercept = 0,  # log2(1) = 0
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low_log2,  # 使用log2转换后的置信区间下限
      xmax = conf.high_log2, # 使用log2转换后的置信区间上限
      color = quintile
    ),
    height = 0.3,  # 与OR_death一致：误差线高度
    linewidth = 1,  # 与OR_death一致：线条粗细
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    aes(
      fill = quintile,
      color = quintile
    ),
    size = 4,       # 与OR_death一致：圆圈大小
    shape = 21,     # 与OR_death一致：带边框的圆形
    stroke = 1,     # 与OR_death一致：边框粗细
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(values = quintile_colors) +
  scale_fill_manual(values = quintile_colors) +
  scale_x_continuous(
    labels = function(x) format(round(x, 1), nsmall = 1)  # 与OR_death一致：统一保留1位小数
  ) +
  facet_grid(metric ~ Disease_Code, scales = "free_x") +
  labs(
    x = "Log2(Odds ratio)",  # 与OR_death一致：坐标轴标签
    y = "Quintile",
    title = "",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_minimal() +  # 移除base_size，与OR_death一致
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 9),  # 与OR_death一致：字体大小
    axis.text = element_text(color = "black", size = 7), # 与OR_death一致：字体大小
    axis.text.y = element_text(size = 8),  # 与OR_death一致：Y轴文字稍大
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 9),  # 与OR_death一致：分面标签大小
    strip.text.x = element_text(size = 9, face = "bold"), # 与OR_death一致
    strip.text.y = element_text(size = 9, face = "bold"), # 与OR_death一致
    legend.position = "none",  # 与OR_death一致：隐藏图例
    plot.margin = margin(-20, 0, 1, 0),  # 与OR_death一致：边距设置
  )

# 显示图形
print(OR_plot3)


top_diseases <- c("K82", "K61","K01","K91")


# 筛选这些疾病的所有Q1-Q5数据
top_diseases_data <- or_result %>%
  filter(Disease_Code %in% top_diseases)

# 分别处理每种指标的数据
# 分别处理每种指标的数据
process_metric_data <- function(data, metric) {
  # 检查所需的列是否存在
  required_cols <- c(paste0("Odds_ratio_", metric), 
                     paste0("conf.low_", metric), 
                     paste0("conf.high_", metric))
  
  # 如果缺少任何必需的列，返回空数据框
  if (!all(required_cols %in% names(data))) {
    warning(paste("Columns for metric", metric, "not found in data"))
    return(data.frame())
  }
  
  metric_data <- data %>%
    dplyr::select(Disease_Code, quintile, 
                  estimate = paste0("Odds_ratio_", metric), 
                  conf.low = paste0("conf.low_", metric), 
                  conf.high = paste0("conf.high_", metric)) %>%
    dplyr::mutate(metric = case_when(
      metric == "proBMIΔ" ~ "proBMIΔ",
      metric == "BMI" ~ "BMI",
      metric == "proBMI" ~ "proBMI",
      metric == "proWCΔ" ~ "proWCΔ",  # 修正这里
      metric == "WC" ~ "WC",
      metric == "proWC" ~ "proWC",
      TRUE ~ metric
    ))
  return(metric_data)
}

# 处理所有指标
# 筛选只保留 WC、proWC、proWCΔ 三个指标
selected_metrics <- c("WC", "proWC", "proWCΔ")

# 处理选定的指标
selected_metric_data <- purrr::map_dfr(selected_metrics, ~ process_metric_data(top_diseases_data, .x))

# 设置因子顺序 - 按 WC、proWC、proWCΔ 顺序
selected_metric_data$metric <- factor(selected_metric_data$metric, 
                                      levels = c("WC", "proWC", "proWCΔ"))  # 修改顺序
selected_metric_data$quintile <- factor(selected_metric_data$quintile, 
                                        levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# 为 Q1-Q5 分配颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",  # 蓝色
  "Q2" = "#ffcc00",  # 黄色
  "Q3" = "#2ca02c",  # 绿色
  "Q4" = "#e377c2",  # 粉色
  "Q5" = "#d62728"   # 红色
)

# 对数据进行log2转换
selected_metric_data_log2 <- selected_metric_data %>%
  mutate(
    estimate_log2 = log2(estimate),
    conf.low_log2 = log2(conf.low),
    conf.high_log2 = log2(conf.high)
  )

# 绘制图形
OR_plot4 <- ggplot(selected_metric_data_log2, aes(
  x = estimate_log2,  # 使用log2转换后的OR值
  y = quintile
)) +
  geom_vline(
    xintercept = 0,  # log2(1) = 0
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low_log2,  # 使用log2转换后的置信区间下限
      xmax = conf.high_log2, # 使用log2转换后的置信区间上限
      color = quintile
    ),
    height = 0.3,  # 与OR_death一致：误差线高度
    linewidth = 1,  # 与OR_death一致：线条粗细
    position = position_dodge(width = 0.3)
  ) +
  geom_point(
    aes(
      fill = quintile,
      color = quintile
    ),
    size = 4,       # 与OR_death一致：圆圈大小
    shape = 21,     # 与OR_death一致：带边框的圆形
    stroke = 1,     # 与OR_death一致：边框粗细
    position = position_dodge(width = 0.3)
  ) +
  scale_color_manual(values = quintile_colors) +
  scale_fill_manual(values = quintile_colors) +
  scale_x_continuous(
    labels = function(x) format(round(x, 1), nsmall = 1)  # 与OR_death一致：统一保留1位小数
  ) +
  facet_grid(metric ~ Disease_Code, scales = "free_x") +
  labs(
    x = "Log2(Odds ratio)",  # 与OR_death一致：坐标轴标签
    y = "Quintile",
    title = "",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_minimal() +  # 移除base_size，与OR_death一致
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 9),  # 与OR_death一致：字体大小
    axis.text = element_text(color = "black", size = 7), # 与OR_death一致：字体大小
    axis.text.y = element_text(size = 8),  # 与OR_death一致：Y轴文字稍大
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    strip.text = element_text(face = "bold", size = 9),  # 与OR_death一致：分面标签大小
    strip.text.x = element_text(size = 9, face = "bold"), # 与OR_death一致
    strip.text.y = element_text(size = 9, face = "bold"), # 与OR_death一致
    legend.position = "none",  # 与OR_death一致：隐藏图例
    plot.margin = margin(-20, 0, 1, 0),  # 与OR_death一致：边距设置
  )

# 显示图形
print(OR_plot4)


library(tidyverse)
library(ggplot2)

or_result <- read.csv("hbh_or_results.csv", header = TRUE, fileEncoding = "GBK")

# 定义疾病分组
disease_groups <- list(
  group1 = c("G44", "G83"),
  group2 = c("I07", "I05", "I89", "I87", "I60"),
  group3 = c("J39"),
  group4 = c("K82", "K61", "K01", "K91")
)

# 定义颜色
quintile_colors <- c(
  "Q1" = "#1f77b4",
  "Q2" = "#ffcc00", 
  "Q3" = "#2ca02c",
  "Q4" = "#e377c2",
  "Q5" = "#d62728"
)

# 处理指标数据的函数 - 添加p值处理
process_metric_data <- function(data, metric) {
  required_cols <- c(
    paste0("Odds_ratio_", metric), 
    paste0("conf.low_", metric), 
    paste0("conf.high_", metric),
    paste0("p.adj_", metric)  # 添加p值列
  )
  
  if (!all(required_cols %in% names(data))) {
    return(data.frame())
  }
  
  metric_data <- data %>%
    dplyr::select(
      Disease_Code, quintile, 
      estimate = paste0("Odds_ratio_", metric), 
      conf.low = paste0("conf.low_", metric), 
      conf.high = paste0("conf.high_", metric),
      p.adj = paste0("p.adj_", metric)  # 选择p值
    ) %>%
    dplyr::mutate(
      metric = case_when(
        metric == "proBMIΔ" ~ "proBMIΔ",
        metric == "BMI" ~ "BMI",
        metric == "proBMI" ~ "proBMI",
        metric == "proWCΔ" ~ "proWCΔ",
        metric == "WC" ~ "WC",
        metric == "proWC" ~ "proWC",
        TRUE ~ metric
      )
    )
  return(metric_data)
}

# 修改后的OR图绘制函数
create_or_plot <- function(disease_codes, group_name = "") {
  # 筛选数据
  top_diseases_data <- or_result %>%
    filter(Disease_Code %in% disease_codes)
  
  # 处理选定的指标
  selected_metrics <- c("WC", "proWC", "proWCΔ")
  selected_metric_data <- purrr::map_dfr(selected_metrics, ~ process_metric_data(top_diseases_data, .x))
  
  # 设置因子顺序
  selected_metric_data$metric <- factor(selected_metric_data$metric, 
                                        levels = c("WC", "proWC", "proWCΔ"))
  selected_metric_data$quintile <- factor(selected_metric_data$quintile, 
                                          levels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
  
  # 对数据进行log2转换，并添加形状标识
  selected_metric_data_log2 <- selected_metric_data %>%
    mutate(
      estimate_log2 = log2(estimate),
      conf.low_log2 = log2(conf.low),
      conf.high_log2 = log2(conf.high),
      # 确定点形状：Q1保持圆形(16)，其他组p<0.05的用三角形(17)
      point_shape = case_when(
        quintile == "Q1" ~ 16,  # 圆形，参考组
        p.adj < 0.05 ~ 17,      # 三角形，显著
        TRUE ~ 16               # 圆形，不显著
      )
    )
  
  # 绘制图形
  ggplot(selected_metric_data_log2, aes(
    x = estimate_log2,
    y = quintile
  )) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    geom_errorbarh(
      aes(
        xmin = conf.low_log2,
        xmax = conf.high_log2,
        color = quintile
      ),
      height = 0.5,
      linewidth = 1,
      position = position_dodge(width = 0.3)
    ) +
    # 使用geom_point并根据形状参数绘制不同的点
    geom_point(
      aes(
        fill = quintile,
        color = quintile,
        shape = factor(point_shape)  # 将形状映射为因子
      ),
      size = 2,
      stroke = 1,
      position = position_dodge(width = 0.3)
    ) +
    # 手动设置形状：16=圆形，17=三角形
    scale_shape_manual(
      values = c("16" = 16, "17" = 17),
      guide = "none"  # 不显示图例
    ) +
    scale_color_manual(values = quintile_colors) +
    scale_fill_manual(values = quintile_colors) +
    scale_x_continuous(
      labels = function(x) format(round(x, 1), nsmall = 1)
    ) +
    facet_grid(metric ~ Disease_Code, scales = "free_x") +
    labs(
      x = expression(bold(log[bold("2")] * "(Odds Ratio)")),
      y = "Quintile",
      title = group_name,
      color = "Quintile",
      fill = "Quintile"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.8),
      axis.ticks.x = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 11, family = "Arial"),
      axis.text = element_text(color = "black", size = 9, family = "Arial"),
      axis.text.y = element_text(size = 9, family = "Arial"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
      strip.text = element_text(face = "bold", size = 11, family = "Arial"),
      strip.text.x = element_text(size = 11, face = "bold", family = "Arial"),
      strip.text.y = element_text(size = 11, face = "bold", family = "Arial"),
      legend.position = "none",
      plot.margin = margin(1, 0, 1, 0),
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, size = 12, family = "Arial", face = "bold"),
      legend.text = element_text(family = "Arial"),
      legend.title = element_text(family = "Arial")
    )
}

# 创建OR图
cat("正在生成OR图...\n")
OR_plot1 <- create_or_plot(disease_groups$group1, "Nervous system disorders")
OR_plot2 <- create_or_plot(disease_groups$group2, "Circulatory system disorders") 
OR_plot3 <- create_or_plot(disease_groups$group3, "Respiratory system disorders")
OR_plot4 <- create_or_plot(disease_groups$group4, "Digestive system disorders")

# 显示图形
print(OR_plot1)
print(OR_plot2)
print(OR_plot3)
print(OR_plot4)












