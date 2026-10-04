#!/usr/bin/env Rscript
# ============================================================================
# ggbiomed.R — R ggplot2 生物医学出版级作图 helper
# 配套 skill: r-ggplot2-plotting
# 依赖: ggplot2, patchwork, scales, ggrepel, RColorBrewer, viridis,
#       survminer, pROC, ggpubr, ggforce, ggridges, gridExtra
# 用户库: ~/Library/R/<major.minor>/library（helper 按 R.version 动态定位）
# ============================================================================

# 库路径处理：
#  - conda R（R.home 含 conda/miniconda，或 CONDA_PREFIX 非空）：直接走 conda 自带库，无需注入
#  - 系统 R：动态注入用户级库 ~/Library/R/<major.minor>/library，保证 source() 时能定位已装包
.is_conda <- (nzchar(Sys.getenv("CONDA_PREFIX")) ||
              grepl("conda|miniconda|mamba", R.home(), ignore.case = TRUE))
if (!.is_conda) {
  .rv <- paste0(R.version$major, ".", strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1])
  .user_lib <- file.path(Sys.getenv("HOME"), "Library/R", .rv, "library")
  if (dir.exists(.user_lib)) .libPaths(c(.user_lib, .libPaths()))
  rm(.rv, .user_lib)
}
rm(.is_conda)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})
# ggpubr 可选：依赖 car/lmtest（需 Fortran）在源码安装环境可能缺失；
# 缺失时仅损失 stat_compare_means 等辅助函数，核心作图不受影响
if (requireNamespace("ggpubr", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggpubr))
} else {
  cat("（提示）ggpubr 未安装，stat_* 辅助函数不可用，核心作图函数正常。\n")
}

# ---------------------------------------------------------------------------
# 0. 默认离散配色 = Cell 期刊五套配色（社区扒自人类巨噬细胞图谱论文）
#     总调度约定：DM/SCI 相关作图自动调用，无需每次手动指定。
# ---------------------------------------------------------------------------
.skill_dir <- "/Users/zuoxianbo/.workbuddy/skills/r-ggplot2-plotting"
if (file.exists(file.path(.skill_dir, "cell_palettes.R"))) {
  source(file.path(.skill_dir, "cell_palettes.R"))
  # 任何带离散 colour/fill 美学的 ggplot，若无显式 scale，默认走 Cell 配色
  options(ggplot2.discrete.colour = scale_color_cell,
          ggplot2.discrete.fill   = scale_fill_cell)
  # 便捷别名（SCI = Cell journal palettes）
  scale_color_sci <- scale_color_cell
  scale_fill_sci  <- scale_fill_cell
  # 切换默认套数/变体（如 >6 组改用 10/15 色）
  use_cell_palette <- function(n = 6, variant = "a") {
    options(ggplot2.discrete.colour = function(...) scale_color_cell(n, variant),
            ggplot2.discrete.fill   = function(...) scale_fill_cell(n, variant))
    invisible(NULL)
  }
  cat("（配色）已加载 Cell 期刊五套配色，离散 colour/fill 默认自动调用。\n")
} else {
  cat("（提示）cell_palettes.R 未找到，跳过 Cell 配色默认。\n")
}

# ---------------------------------------------------------------------------
# 1. 出版级主题 (Nature / Wiley / Cell 投稿风格)
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 11, font = "Arial") {
  theme_classic(base_size = base_size, base_family = font) +
    theme(
      plot.title    = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 1, color = "black"),
      axis.line     = element_line(linewidth = 0.6, color = "black"),
      axis.ticks    = element_line(linewidth = 0.6, color = "black"),
      legend.title  = element_text(size = base_size - 1, face = "bold"),
      legend.text   = element_text(size = base_size - 1),
      legend.key    = element_blank(),
      legend.background = element_blank(),
      panel.border  = element_blank(),
      strip.background = element_blank(),
      strip.text    = element_text(face = "bold"),
      plot.margin   = margin(0.4, 0.4, 0.4, 0.4, "cm")
    )
}

# 左上角加 a/b/c/d 面板标签
panel <- function(p, letter, x = 0.02, y = 0.98, size = 14) {
  p + annotate("text", x = x, y = y, label = letter,
               size = size / .pt, fontface = "bold", hjust = 0, vjust = 1,
               inherit.aes = FALSE, parse = FALSE)
}

# ---------------------------------------------------------------------------
# 2. 火山图 (volcano)
# ---------------------------------------------------------------------------
gg_volcano <- function(df, log2fc = "log2FC", pval = "pval",
                       label_col = NULL, fdr = 0.05, fc = 1,
                       top_n = 0, group = NULL,
                       xlab = "log2 fold change",
                       ylab = expression(-log[10](p))) {
  df <- as.data.frame(df)
  lfc <- df[[log2fc]]; pv <- df[[pval]]
  df$.neglog10p <- -log10(pv)
  if (!is.null(group)) {
    # 多组别火山图：按 group 列着色（来自 PDF 案例「多组别火山图」）
    .ng <- length(unique(na.omit(df[[group]])))
    p <- ggplot(df, aes(x = .data[[log2fc]], y = .neglog10p,
                        color = .data[[group]])) +
      geom_point(size = 1.6, alpha = 0.75) +
      scale_color_cell(n = if (.ng <= 3) 3 else if (.ng <= 6) 6 else if (.ng <= 10) 10 else 15) +
      geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey40") +
      theme_pub() + xlab(xlab) + ylab(ylab) + guides(color = guide_legend(title = ""))
  } else {
    df$.group <- ifelse(pv < fdr & lfc >  fc, "Up",
                 ifelse(pv < fdr & lfc < -fc, "Down", "NS"))
    p <- ggplot(df, aes(x = .data[[log2fc]], y = .neglog10p, color = .group)) +
      geom_point(size = 1.6, alpha = 0.75) +
      scale_color_manual(values = c(Up = "#E64B35", Down = "#3182BD", NS = "#BBBBBB")) +
      geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey40") +
      theme_pub() + xlab(xlab) + ylab(ylab) + guides(color = "none")
  }
  if (!is.null(label_col) && top_n > 0) {
    lab <- df[order(df$.neglog10p, decreasing = TRUE)[1:min(top_n, nrow(df))], ]
    p <- p + geom_text_repel(data = lab, aes(label = .data[[label_col]]),
                             size = 3, max.overlaps = 20)
  }
  p
}

# ---------------------------------------------------------------------------
# 3. 曼哈顿图 (GWAS / PGx)
# ---------------------------------------------------------------------------
gg_manhattan <- function(df, chr = "CHR", pos = "POS", pval = "P",
                         label_col = NULL, sig = 5e-8) {
  df <- as.data.frame(df)
  df$.logp <- -log10(df[[pval]])
  # 染色体累加偏移
  chr_order <- sort(unique(df[[chr]]))
  offsets <- setNames(cumsum(c(0, tapply(df[[pos]], df[[chr]], max, na.rm = TRUE)[-length(chr_order)])),
                      chr_order)
  df$.x <- df[[pos]] + offsets[as.character(df[[chr]])]
  cols <- rep(c("#1B4F72", "#E67E22"), length.out = length(chr_order))
  p <- ggplot(df, aes(x = .x, y = .logp, color = factor(.data[[chr]]))) +
    geom_point(size = 1.2) +
    scale_color_manual(values = cols) +
    geom_hline(yintercept = -log10(sig), linetype = "dashed", color = "red") +
    theme_pub() + guides(color = "none") +
    xlab("Chromosome") + ylab(expression(-log[10](p)))
  p
}

# ---------------------------------------------------------------------------
# 4. Kaplan-Meier 生存曲线（纯 R 实现，无 survminer/survival 依赖）
# ---------------------------------------------------------------------------
gg_km <- function(formula, data, risk.table = TRUE, pval = TRUE,
                  palette = c("#E64B35", "#3182BD", "#00A087")) {
  # 优先用 survminer（若已装）
  if (requireNamespace("survminer", quietly = TRUE) &&
      requireNamespace("survival", quietly = TRUE)) {
    suppressPackageStartupMessages({ library(survminer); library(survival) })
    p <- ggsurvplot(surv_fit(formula, data = data), data = data,
                    risk.table = risk.table, pval = pval, palette = palette,
                    ggtheme = theme_classic(), legend.title = "",
                    xlab = "Time", ylab = "Survival probability", fontsize = 4)
    return(p)
  }
  # 解析 Surv(time,status) ~ group（直接按列名取值，不依赖 survival 包）
  terms <- all.vars(formula)   # c("time","status",[group])
  time   <- data[[terms[1]]]
  status <- data[[terms[2]]]
  if (length(terms) >= 3 && terms[3] %in% names(data)) {
    group <- as.factor(data[[terms[3]]])
  } else {
    group <- factor(rep("All", nrow(data)))
  }

  # 单组 KM 计算
  km_one <- function(t, s) {
    ord <- order(t); t <- t[ord]; s <- s[ord]
    n <- length(t); at_risk <- n; surv <- 1
    rows <- data.frame(time = 0, surv = 1, at_risk = n)
    for (tt in unique(t)) {
      d <- sum(s[t == tt] == 1)
      if (d > 0) {
        surv <- surv * (at_risk - d) / at_risk
        rows <- rbind(rows, data.frame(time = tt, surv = surv, at_risk = at_risk - d))
      }
      at_risk <- at_risk - sum(t == tt)
    }
    rows
  }
  # log-rank 检验（两组）
  lr_p <- NA
  if (pval && length(levels(group)) == 2) {
    ut <- sort(unique(time[status == 1]))
    g1 <- levels(group)[1]
    O1 <- E1 <- 0; V <- 0
    for (tt in ut) {
      n_tot <- sum(time >= tt); d_tot <- sum(time == tt & status == 1)
      if (d_tot == 0 || n_tot < 2) next
      n1 <- sum(group == g1 & time >= tt)
      O1 <- O1 + sum(group == g1 & time == tt & status == 1)
      E1 <- E1 + d_tot * n1 / n_tot
      V  <- V  + d_tot * (n_tot - d_tot) / (n_tot * (n_tot - 1)) * n1 * (n_tot - n1)
    }
    if (V > 0) {
      chi2 <- (O1 - E1)^2 / V
      lr_p <- pchisq(chi2, df = 1, lower.tail = FALSE)
    }
  }
  # 绘图
  df_list <- list()
  for (g in levels(group)) {
    k <- km_one(time[group == g], status[group == g])
    k$group <- g
    df_list[[g]] <- k
  }
  kmdf <- do.call(rbind, df_list)
  p <- ggplot(kmdf, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 1) +
    scale_color_manual(values = palette) +
    theme_pub() + xlab("Time") + ylab("Survival probability") +
    ylim(0, 1) + guides(color = guide_legend(title = ""))
  if (!is.na(lr_p)) p <- p + ggtitle(sprintf("Log-rank p = %.3g", lr_p))
  if (risk.table) {
    rt <- kmdf[kmdf$time == 0 | kmdf$surv < 1, ]
    rt$txt <- paste0(rt$at_risk)
    p_rt <- ggplot(kmdf, aes(x = time, y = group, label = at_risk)) +
      geom_text(size = 3) + theme_void() + xlab("Number at risk")
    p <- p / p_rt
  }
  p
}

# ---------------------------------------------------------------------------
# 5. ROC 曲线
# ---------------------------------------------------------------------------
gg_roc <- function(response, score, title = "ROC") {
  suppressPackageStartupMessages(library(pROC))
  r <- as.numeric(as.factor(response)) - 1
  roc_obj <- roc(r, score)
  ci <- ci.auc(roc_obj)
  df <- data.frame(
    FPR = 1 - rev(roc_obj$specificities),
    TPR = rev(roc_obj$sensitivities)
  )
  ggplot(df, aes(x = FPR, y = TPR)) +
    geom_line(linewidth = 1, color = "#E64B35") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    theme_pub() +
    xlab("1 - Specificity (False Positive Rate)") + ylab("Sensitivity (True Positive Rate)") +
    ggtitle(sprintf("%s  AUC = %.3f (%.3f–%.3f)", title, auc(roc_obj), ci[1], ci[3]))
}

# ---------------------------------------------------------------------------
# 6. 森林图 (OR / HR)
# ---------------------------------------------------------------------------
gg_forest <- function(df, est = "est", lo = "lo", hi = "hi",
                      label = "label", xlab = "Odds Ratio (log scale)") {
  df <- as.data.frame(df)
  df$.y <- seq_len(nrow(df))
  ggplot(df, aes(x = .data[[est]], y = .y)) +
    geom_point(size = 2.5, color = "#1B4F72") +
    geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]), height = 0.2,
                  color = "#1B4F72", orientation = "y") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() + theme_pub() +
    ylab("") + xlab(xlab) +
    scale_y_continuous(breaks = df$.y, labels = df[[label]], trans = "reverse") +
    theme(axis.text.y = element_text(hjust = 0))
}

# ---------------------------------------------------------------------------
# 7. 热图 (纯 ggplot2 geom_tile, 支持行列聚类; 与 ggsave 兼容)
# ---------------------------------------------------------------------------
gg_heatmap <- function(mat, scale = "row", cluster_rows = TRUE,
                       cluster_cols = TRUE, annotation_col = NULL,
                       show_rownames = FALSE, color_low = "#3182BD",
                       color_mid = "#FFFFEE", color_high = "#E64B35") {
  mat <- as.matrix(mat)
  if (scale == "row")    mat <- t(scale(t(mat)))
  else if (scale == "column") mat <- scale(mat)
  if (cluster_rows && nrow(mat) > 1) mat <- mat[hclust(dist(mat))$order, , drop = FALSE]
  if (cluster_cols && ncol(mat) > 1) mat <- mat[, hclust(dist(t(mat)))$order, drop = FALSE]
  df <- as.data.frame(as.table(mat))
  names(df) <- c("row", "col", "value")
  df$row <- factor(df$row, levels = rev(unique(df$row)))  # 顶部=第一行
  p <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = color_low, mid = color_mid, high = color_high,
                         midpoint = if (scale == "none") mean(df$value) else 0) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if (!show_rownames) p <- p + theme(axis.text.y = element_blank())
  p
}

# ---------------------------------------------------------------------------
# 8. PCA 散点
# ---------------------------------------------------------------------------
gg_pca <- function(mat, group, center = TRUE, scale. = TRUE,
                  ellipse = TRUE, palette = c("#E64B35", "#3182BD", "#00A087",
                                             "#F0A202", "#7E57C2", "#26A69A")) {
  pca <- prcomp(mat, center = center, scale. = scale.)
  sd <- pca$sdev^2 / sum(pca$sdev^2) * 100
  coord <- as.data.frame(pca$x[, 1:2])
  coord$group <- group
  colnames(coord)[1:2] <- c("PC1", "PC2")
  p <- ggplot(coord, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = palette) +
    theme_pub() +
    xlab(sprintf("PC1 (%.1f%%)", sd[1])) +
    ylab(sprintf("PC2 (%.1f%%)", sd[2]))
  if (ellipse) p <- p + stat_ellipse(level = 0.95, alpha = 0.4)
  p
}

# ---------------------------------------------------------------------------
# 9. 剂量-反应 / IC50 (四参数 log-logistic)
# ---------------------------------------------------------------------------
gg_dose_response <- function(df, dose = "dose", resp = "resp",
                            group = NULL, logdose = TRUE) {
  df <- as.data.frame(df)
  p <- ggplot(df, aes(x = .data[[dose]], y = .data[[resp]]))
  if (!is.null(group)) p <- p + aes(color = .data[[group]])
  p <- p + geom_point(size = 2) +
    scale_y_continuous() + theme_pub() +
    xlab(if (logdose) "Dose (log10)" else "Dose") +
    ylab("Response")
  if (logdose) p <- p + scale_x_log10()
  p
}

# ---------------------------------------------------------------------------
# 10. 带误差棒柱状图
# ---------------------------------------------------------------------------
gg_bar_err <- function(df, x, y, group = NULL, err = NULL,
                       dodge = TRUE, ylab = "Value") {
  df <- as.data.frame(df)
  pos <- if (is.null(err)) position_identity() else position_dodge(0.9)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(group)) p <- p + aes(fill = .data[[group]])
  p <- p + geom_col(position = pos, width = 0.7, color = "black") +
    theme_pub() + ylab(ylab) + xlab("")
  if (!is.null(err)) {
    p <- p + geom_errorbar(aes(ymin = .data[[y]] - .data[[err]],
                               ymax = .data[[y]] + .data[[err]]),
                           position = pos, width = 0.25)
  }
  p
}

# ---------------------------------------------------------------------------
# 11. 箱线 + 小提琴 + 抖点
# ---------------------------------------------------------------------------
gg_box_violin <- function(df, x, y, group = NULL, jit = TRUE) {
  df <- as.data.frame(df)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(group)) p <- p + aes(fill = .data[[group]])
  p <- p + geom_violin(alpha = 0.3, trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    theme_pub() + xlab("") + ylab(y)
  if (jit) p <- p + geom_jitter(width = 0.1, size = 0.8, alpha = 0.5)
  p
}

# ---------------------------------------------------------------------------
# 12. 多面板拼图 (patchwork)
# ---------------------------------------------------------------------------
multi_panel <- function(..., ncol = NULL, nrow = NULL, widths = NULL,
                        heights = NULL, guides = "collect") {
  plots <- list(...)
  args <- c(plots, list(ncol = ncol, nrow = nrow, widths = widths,
                       heights = heights, guides = guides))
  args <- args[!sapply(args, is.null)]
  do.call(wrap_plots, args)
}

# ---------------------------------------------------------------------------
# 13. 哑铃图 (dumbbell) — PDF 案例「个性化哑铃图」
#     连接两端数值，让组间差异一眼可见
# ---------------------------------------------------------------------------
gg_dumbbell <- function(df, id = "id", x1 = "x1", x2 = "x2", group = NULL,
                        xlab = "Value", ylab = "", sort = TRUE) {
  df <- as.data.frame(df)
  if (sort) {
    ord <- order(df[[x2]] - df[[x1]], decreasing = TRUE)
    df <- df[ord, ]
    df[[id]] <- factor(df[[id]], levels = rev(unique(df[[id]])))
  }
  p <- ggplot(df, aes(y = .data[[id]])) +
    geom_segment(aes(x = .data[[x1]], xend = .data[[x2]]),
                 color = "grey50", linewidth = 0.8) +
    geom_point(aes(x = .data[[x1]]), color = "#7F7F7F", size = 3) +
    geom_point(aes(x = .data[[x2]]), color = "#E64B35", size = 3) +
    theme_pub() + xlab(xlab) + ylab(ylab)
  if (!is.null(group)) {
    p <- p + aes(color = .data[[group]]) +
      scale_color_brewer(palette = "Set2") +
      guides(color = guide_legend(title = ""))
  }
  p
}

# ---------------------------------------------------------------------------
# 14. 云雨图 (raincloud) — PDF 案例「云雨图进阶」
#     小提琴(云)+箱线(雨)+抖点，把分布/摘要/个体数据放在一起
# ---------------------------------------------------------------------------
gg_raincloud <- function(df, x, y, group = NULL, ylab = "Value") {
  df <- as.data.frame(df)
  p <- ggplot(df, aes(x = .data[[y]], y = .data[[x]]))
  if (!is.null(group)) p <- p + aes(fill = .data[[group]])
  p <- p +
    geom_violin(alpha = 0.3, trim = FALSE, scale = "width", color = "black", width = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", alpha = 0.8) +
    geom_jitter(height = 0, width = 0.08, size = 0.7, alpha = 0.5) +
    theme_pub() + xlab("") + ylab(ylab) + coord_flip()
  p
}

# ---------------------------------------------------------------------------
# 15. 游泳图 (swimmer) — PDF 案例「游泳图」
#     用时间轴组织个体过程与关键事件（肿瘤治疗时长/进展/死亡）
# ---------------------------------------------------------------------------
gg_swimmer <- function(df, id = "id", start = "start", end = "end",
                       event = NULL, event_time = NULL,
                       xlab = "Time", ylab = "") {
  df <- as.data.frame(df)
  p <- ggplot(df, aes(x = .data[[start]], xend = .data[[end]], y = .data[[id]])) +
    geom_segment(aes(x = .data[[start]], xend = .data[[end]], y = .data[[id]],
                     yend = .data[[id]]), color = "#3182BD", linewidth = 3) +
    theme_pub() + xlab(xlab) + ylab(ylab)
  if (!is.null(event) && !is.null(event_time)) {
    ev <- data.frame(id2 = df[[id]], et = df[[event_time]], type = df[[event]])
    p <- p + geom_point(data = ev, aes(x = et, y = id2, color = type),
                        size = 3, inherit.aes = FALSE) +
      guides(color = guide_legend(title = "Event"))
  }
  p
}

# ---------------------------------------------------------------------------
# 16. 面积图 (area) — PDF 案例「个性化面积图」
#     多组变化曲线与排序信息组合展示
# ---------------------------------------------------------------------------
gg_area <- function(df, x, y, group = NULL, stacked = FALSE,
                   xlab = "", ylab = "Value") {
  df <- as.data.frame(df)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(group)) p <- p + aes(fill = .data[[group]], color = .data[[group]])
  p <- p + geom_area(alpha = 0.6,
                     position = if (stacked) "stack" else "identity") +
    theme_pub() + xlab(xlab) + ylab(ylab)
  p
}

# ---------------------------------------------------------------------------
# 17. 相关性热图 (correlation heatmap) — PDF 案例「相关性热图」等
#     计算相关矩阵 + 可选显著性星标
# ---------------------------------------------------------------------------
gg_corr <- function(mat, method = "pearson", sig = TRUE,
                    color_low = "#3182BD", color_mid = "white", color_high = "#E64B35",
                    label = FALSE) {
  mat <- as.matrix(mat)
  cmat <- cor(mat, method = method, use = "pairwise.complete.obs")
  pmat <- matrix(NA, nrow(cmat), ncol(cmat))
  if (sig) {
    for (i in 1:nrow(cmat)) for (j in 1:ncol(cmat)) {
      if (i != j) {
        t <- try(cor.test(mat[, i], mat[, j], method = method), silent = TRUE)
        pmat[i, j] <- if (inherits(t, "try-error")) NA else t$p.value
      }
    }
  }
  df <- as.data.frame(as.table(cmat))
  names(df) <- c("Var1", "Var2", "value")
  if (sig) {
    pv <- as.data.frame(as.table(pmat)); names(pv) <- c("Var1", "Var2", "p")
    df <- merge(df, pv, by = c("Var1", "Var2"))
    df$stars <- ifelse(df$p < 0.001, "***",
                ifelse(df$p < 0.01, "**",
                ifelse(df$p < 0.05, "*", "")))
  }
  p <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = color_low, mid = color_mid, high = color_high,
                         midpoint = 0, limit = c(-1, 1)) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if (label) {
    lab_df <- df
    lab_df$lab <- if (sig) paste0(sprintf("%.2f", df$value), df$stars)
                  else sprintf("%.2f", df$value)
    p <- p + geom_text(data = lab_df, aes(label = lab), color = "black", size = 3)
  }
  p
}

# ---------------------------------------------------------------------------
# 导出提示
# ---------------------------------------------------------------------------
#' 统一保存: ggsave("out.png", p, width=6, height=5, dpi=300)
#' 矢量: ggsave("out.pdf", p, width=6, height=5)
#' 调用示例见 SKILL.md
cat("ggbiomed.R loaded. Functions: gg_volcano, gg_manhattan, gg_km, gg_roc,\n",
    "gg_forest, gg_heatmap, gg_pca, gg_dose_response, gg_bar_err,\n",
    "gg_box_violin, gg_dumbbell, gg_raincloud, gg_swimmer, gg_area, gg_corr,\n",
    "multi_panel, theme_pub, panel\n")
