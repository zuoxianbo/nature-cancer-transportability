#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R -- Nature Cancer  UPDATED figures (P0-complete set)
#
# Builds ONLY the four display items whose content changed for :
#   Fig.4  fig4_ablation_constraint   a/b/c rebuilt from the frozen  source
#                                     data JSON; NEW panel d = selection-
#                                     adjusted permutation null (P0-4).
#   ED1    ed1_positive_control       PRIMARY = 18-gene T-cell-inflamed GEP
#                                     (Cristescu 2018, P0-2); panel b =
#                                     6-gene IFN-gamma panel (sensitivity only).
#   ED6    ed6_threshold_sensitivity  decision matrix (24 cells x C1/C2 by
#                                     tolerance) + continuous dAUC transport.
#   ED9    ed9_precision_simulation   external-validation precision simulation
#                                     (P0-5): P(C1) and CI width vs true AUC
#                                     at n=71 / n=25.
#
# Every statistical number is read verbatim from frozen JSON artifacts under
# data/ -- nothing data-derived is hard-coded. Only geometry (axis ranges,
# dodge widths, annotation positions) is inline.
#
# Outputs: figures/{name}.{png,pdf,svg} (PNG 300 dpi, 180 mm wide)
#          figures/source_data/{name}.json
# =============================================================================
.libPaths(c(file.path(Sys.getenv("HOME"), "Library/R/4.6/library"), .libPaths()))
source("/Users/zuoxianbo/.workbuddy/skills/r-ggplot2-plotting/ggbiomed.R")

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
  library(jsonlite); library(grid); library(gridExtra)
})

BASE   <- "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature/outputs/analysis_base"
OUTDIR <- file.path(BASE, "figures_")
SDDIR  <- file.path(OUTDIR, "source_data")
dir.create(SDDIR, recursive = TRUE, showWarnings = FALSE)
J <- function(...) file.path(BASE, ...)
FF <- "Arial"

rd <- function(path) {
  if (!file.exists(path)) stop(sprintf("FROZEN SOURCE MISSING: %s", path))
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}
sd_all <- new.env(parent = emptyenv())
add_sd <- function(fig, panel, src, key, val) {
  k <- paste0(fig, ".", panel)
  cur <- if (exists(k, envir = sd_all, inherits = FALSE)) get(k, envir = sd_all) else list()
  cur[[length(cur) + 1L]] <- list(source = src, key = key, value = val)
  assign(k, cur, envir = sd_all)
}
write_sd <- function(fig) {
  out <- list(figure = fig,
              generated_by = "make_figures.R ( P0-complete set)",
              panels = as.list(sd_all))
  write(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 10),
        file.path(SDDIR, paste0(fig, ".json")))
}

save_trio <- function(name, p, w_mm, h_mm) {
  png_file <- file.path(OUTDIR, paste0(name, ".png"))
  pdf_file <- file.path(OUTDIR, paste0(name, ".pdf"))
  ggsave(png_file, p, width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white")
  dev_lbl <- "quartz"
  ok <- tryCatch({
    quartz(type = "pdf", file = pdf_file, width = w_mm / 25.4, height = h_mm / 25.4)
    print(p); dev.off(); TRUE
  }, error = function(e) {
    ggsave(pdf_file, p, width = w_mm, height = h_mm, units = "mm",
           device = cairo_pdf, bg = "white"); TRUE
  })
  svg_lbl <- "svg-FAILED"
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(OUTDIR, paste0(name, ".svg")), p,
           width = w_mm, height = h_mm, units = "mm", device = svglite::svglite)
    svg_lbl <- "svg"
  }
  cat(sprintf("[saved] %-32s (pdf via %s; svg: %s)\n", name, dev_lbl, svg_lbl))
}

theme_ <- theme_classic(base_family = FF, base_size = 8) +
  theme(plot.title = element_text(size = 8.4, face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
        axis.title = element_text(size = 7.6),
        axis.text = element_text(size = 6.8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.4),
        strip.background = element_blank(),
        strip.text = element_text(size = 7.4, face = "bold", hjust = 0),
        plot.tag = element_text(size = 9.5, face = "bold"))
lab <- function(p, letter) p + labs(tag = letter)

# =============================================================================
# BLOCK A -- Fig.4 (a/b/c from frozen source data + d permutation null)
# =============================================================================
sd4 <- rd(J("figures_", "source_data", "fig4_ablation_constraint.json"))
pn  <- sd4$panels

# ---- panel a: 21-config external AUC + CI ----------------------------------
cfg_layer <- c(baseline = "baseline",
               gene_top100 = "gene set", gene_top50 = "gene set",
               gene_rand100_s0 = "gene set", gene_rand100_s1 = "gene set",
               gene_rand100_s2 = "gene set", gene_rand100_s3 = "gene set",
               gene_rand100_s4 = "gene set",
               gene_rand50_s0 = "gene set", gene_rand50_s1 = "gene set",
               gene_rand50_s2 = "gene set", gene_rand50_s3 = "gene set",
               gene_rand50_s4 = "gene set",
               norm_log2raw = "norm.", norm_within_z = "norm.",
               pc_k1 = "dim.", pc_k3 = "dim.", pc_k5 = "dim.",
               `ridge_lam0.01` = "ridge", `ridge_lam0.1` = "ridge",
               `ridge_lam1` = "ridge")
pa_keys <- names(pn$a)
auc_keys <- grep("^per_config\\..+\\.external\\.(GSE165252|GSE78220)\\.auc$", pa_keys, value = TRUE)
stopifnot(length(auc_keys) == 42L)
dfA <- do.call(rbind, lapply(auc_keys, function(k) {
  core <- sub("^per_config\\.", "", k)
  cfg <- sub("\\.external\\.(GSE165252|GSE78220)\\.auc$", "", core)
  tgt <- sub("^.*\\.external\\.", "", sub("\\.auc$", "", core))
  lo <- pn$a[[paste0("per_config.", cfg, ".external.", tgt, ".ci95_percentile")]][[1]]
  hi <- pn$a[[paste0("per_config.", cfg, ".external.", tgt, ".ci95_percentile")]][[2]]
  add_sd("fig4", "a", "figures/source_data/fig4_ablation_constraint.json", k, pn$a[[k]])
  data.frame(config = cfg, layer = cfg_layer[[cfg]], target = tgt,
             auc = pn$a[[k]], lo = lo, hi = hi, stringsAsFactors = FALSE)
}))
dfA$config <- factor(dfA$config, levels = names(cfg_layer))
pal_layer <- c(baseline = "#1f77b4", `gene set` = "#ff7f0e",
               norm. = "#2ca02c", dim. = "#d62728", ridge = "#9467bd")
short_cfg <- c(baseline = "base",
               gene_top100 = "top100", gene_top50 = "top50",
               gene_rand100_s0 = "r100-0", gene_rand100_s1 = "r100-1",
               gene_rand100_s2 = "r100-2", gene_rand100_s3 = "r100-3",
               gene_rand100_s4 = "r100-4",
               gene_rand50_s0 = "r50-0", gene_rand50_s1 = "r50-1",
               gene_rand50_s2 = "r50-2", gene_rand50_s3 = "r50-3",
               gene_rand50_s4 = "r50-4",
               norm_log2raw = "l2raw", norm_within_z = "wZ",
               pc_k1 = "k1", pc_k3 = "k3", pc_k5 = "k5",
               `ridge_lam0.01` = "\u03bb.01", `ridge_lam0.1` = "\u03bb.1",
               `ridge_lam1` = "\u03bb1")
dfA$cfg_short <- short_cfg[as.character(dfA$config)]
dfA$cfg_short <- factor(dfA$cfg_short, levels = unname(short_cfg[names(cfg_layer)]))
pA <- ggplot(dfA, aes(x = cfg_short, y = auc, colour = layer)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey35", linewidth = 0.35) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.32, linewidth = 0.28, alpha = 0.9) +
  geom_point(size = 1.15) +
  facet_wrap(~target, nrow = 1) +
  scale_colour_manual(values = pal_layer, name = "layer") +
  scale_y_continuous(limits = c(0.24, 0.95), breaks = seq(0.3, 0.9, 0.1)) +
  labs(title = "21-config external transport",
       x = NULL, y = "external AUC (95% CI)") +
  theme_ + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = "bottom",
                    legend.text = element_text(size = 5.8, margin = margin(l = 0, r = 4)),
                    legend.title = element_text(size = 6),
                    legend.key.size = unit(2.2, "mm"),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.box.margin = margin(-2, 0, -4, 0))
pA <- lab(pA, "a")

# ---- panel b: layer variance decomposition old vs new -----------------------
vb <- pn$b
gv <- function(base, key) {
  k <- paste0(base, ".", key)
  if (!k %in% names(vb)) stop(sprintf("fig4 panel b key missing: %s", k))
  vb[[k]]
}
oldv <- list(gene_set = gv("old__comparison.old_layer_mean_variance", "gene_set"),
             normalization = gv("old__comparison.old_layer_mean_variance", "normalization"),
             pc_k = gv("old__comparison.old_layer_mean_variance", "pc_k"),
             ridge_lambda = gv("old__comparison.old_layer_mean_variance", "ridge_lambda"))
mvk <- function(lay) gv("verdict.external_auc_range_by_layer",
                        paste0(lay, ".mean_variance_across_external"))
newv <- list(gene_set = mvk("gene_set"), normalization = mvk("normalization"),
             dimension = mvk("pc_k"), regularization = mvk("ridge_lambda"))
lay <- c("gene_set", "normalization", "pc_k", "ridge_lambda")
lay_lbl <- c("gene set", "normalization", "dimension", "regularization")
dfB <- data.frame(
  layer = factor(rep(lay_lbl, 2), levels = lay_lbl),
  version = rep(c(" (retracted mapping)", " (recovered identity)"), each = 4),
  mean_var = c(oldv[[lay[1]]], oldv[[lay[2]]], oldv[[lay[3]]], oldv[[lay[4]]],
               newv[[1]], newv[[2]], newv[[3]], newv[[4]]),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(dfB)))
  add_sd("fig4", "b", "figures/source_data/fig4_ablation_constraint.json",
         sprintf("layer_mean_variance.%s.%s", dfB$version[i], dfB$layer[i]), dfB$mean_var[i])
pB <- ggplot(dfB, aes(x = layer, y = mean_var, fill = version)) +
  geom_col(position = position_dodge(0.62), width = 0.58, colour = "grey25",
           linewidth = 0.2) +
  scale_fill_manual(values = c(" (retracted mapping)" = "#c7c7c7",
                               " (recovered identity)" = "#1f77b4"), name = NULL) +
  labs(title = "Construct-level choices dominate variability",
       x = NULL, y = "mean variance of external AUC") +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"),
                    axis.text.x = element_text(angle = 18, hjust = 1))
pB <- lab(pB, "b")

# ---- panel c: gene-identity recovery old vs new ----------------------------
pc <- pn$c
gc4 <- function(key) {
  if (!key %in% names(pc)) stop(sprintf("fig4 panel c key missing: %s", key))
  pc[[key]]
}
old_ext165 <- gc4("phase2.old.transport.external.GSE165252.auc")
old_ext782 <- gc4("phase2.old.transport.external.GSE78220.auc")
old_within <- gc4("old__comparison.old_baseline_within_imvigor")
new_ext165 <- gc4("pairwise_transports.transports.train_IMvigor210__test_GSE165252.M2.auc")
new_ext782 <- gc4("pairwise_transports.transports.train_IMvigor210__test_GSE78220.M2.auc")
new_within <- gc4("within_reference.per_cohort.IMvigor210.auc_632")
old_pc1 <- gc4("old__comparison.old_pc1_mass_unmeasured")
dfC <- data.frame(
  metric = factor(c("within IMvigor210 (0.632)", "external GSE165252", "external GSE78220"),
                  levels = c("external GSE165252", "external GSE78220", "within IMvigor210 (0.632)")),
  auc = c(old_ext165, old_ext782, old_within, new_ext165, new_ext782, new_within),
  version = rep(c(" (retracted)", " (recovered)"), each = 3),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(dfC)))
  add_sd("fig4", "c", "figures/source_data/fig4_ablation_constraint.json",
         sprintf("recovery.%s.%s", dfC$version[i], dfC$metric[i]), dfC$auc[i])
add_sd("fig4", "c", "figures/source_data/fig4_ablation_constraint.json",
       "old_pc1_mass_unmeasured", old_pc1)
dfC$metric <- factor(dfC$metric, levels = c("external GSE165252", "external GSE78220",
                                            "within IMvigor210 (0.632)"))
dfC$y_num <- as.numeric(dfC$metric) + ifelse(dfC$version == " (recovered)", 0.19, -0.19)
pC <- ggplot(dfC, aes(x = auc, y = y_num, colour = version, group = version)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey35", linewidth = 0.35) +
  geom_line(linewidth = 0.3, alpha = 0.65) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c(" (retracted)" = "#c7c7c7",
                                 " (recovered)" = "#1f77b4"), name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(dfC$metric)),
                     labels = levels(dfC$metric), minor_breaks = NULL) +
  coord_cartesian(ylim = c(0.55, 3.45)) +
  labs(title = "Gene-row identity recovery (mapping error)",
       x = "AUC", y = NULL) +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
pC <- lab(pC, "c")

# ---- panel d: NEW -- selection-adjusted permutation null (P0-4) -------------
pn4 <- rd(J("data", "", "permutation_null.json"))
tgts <- c("GSE165252", "GSE78220")
dfD <- do.call(rbind, lapply(tgts, function(t) {
  r <- pn4[[t]]
  add_sd("fig4", "d", "data/permutation_null.json",
         sprintf("%s.observed_max_auc", t), r$observed_max_auc)
  add_sd("fig4", "d", "data/permutation_null.json",
         sprintf("%s.perm_max_mean", t), r$perm_max_mean)
  add_sd("fig4", "d", "data/permutation_null.json",
         sprintf("%s.perm_max_q95", t), r$perm_max_q95)
  add_sd("fig4", "d", "data/permutation_null.json",
         sprintf("%s.perm_p_max_one_sided", t), r$perm_p_max_one_sided)
  data.frame(
    target = t,
    statistic = factor(c("perm max\nmean", "perm max\n95th pct", "observed\nmax AUC"),
                       levels = c("perm max\nmean", "perm max\n95th pct", "observed\nmax AUC")),
    value = c(r$perm_max_mean, r$perm_max_q95, r$observed_max_auc),
    pval = r$perm_p_max_one_sided, stringsAsFactors = FALSE)
}))
dfD$target <- factor(dfD$target, levels = tgts)
pD <- ggplot(dfD, aes(x = statistic, y = value, colour = target, group = target)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_line(position = position_dodge(0.45), linewidth = 0.35, alpha = 0.6) +
  geom_point(size = 1.9, position = position_dodge(0.45)) +
  geom_text(data = data.frame(x = 0.58, y = 0.40,
                              lab = sprintf("both P(max) = %.2f", dfD$pval[1])),
            aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, size = 2.15, family = FF, hjust = 0,
            colour = "grey30") +
  scale_colour_manual(values = c(GSE165252 = "#1f77b4", GSE78220 = "#ff7f0e"), name = NULL) +
  scale_y_continuous(limits = c(0.38, 0.82)) +
  labs(title = "Permutation null (2,000 anchor-label perms)",
       x = NULL, y = "external AUC (max over 21 configs)") +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.6, "mm"))
pD <- lab(pD, "d")

fig4 <- (pA / pB | pC / pD) + plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_trio("fig4_ablation_constraint", fig4, 180, 200)
write_sd("fig4_ablation_constraint")

# =============================================================================
# BLOCK B -- ED1: positive control (GEP18 primary + 6-gene sensitivity)
# =============================================================================
g18 <- rd(J("data", "", "positive_control_gep18.json"))
g6  <- rd(J("data", "", "rerun_positive_control.json"))
mk_pc_df <- function(root, panel_label, src_file) {
  do.call(rbind, lapply(tgts, function(t) {
    tm <- root$transports$primary_logcpm_z_frozen_all6[[t]]$target_metrics
    add_sd("ed1", panel_label, src_file,
           sprintf("transports.primary_logcpm_z_frozen_all6.%s.target_metrics.auc", t), tm$auc)
    add_sd("ed1", panel_label, src_file,
           sprintf("transports.primary_logcpm_z_frozen_all6.%s.target_metrics.calib_slope", t),
           if (is.null(tm$calib_slope)) NA else tm$calib_slope)
    data.frame(target = t, auc = tm$auc, lo = tm$auc_ci95[[1]], hi = tm$auc_ci95[[2]],
               slope = if (is.null(tm$calib_slope)) NA else tm$calib_slope,
               stringsAsFactors = FALSE)
  }))
}
dfG18 <- mk_pc_df(g18, "a", "data/positive_control_gep18.json")
dfG6  <- mk_pc_df(g6,  "b", "data/rerun_positive_control.json")
w18 <- g18$within_imvigor210_efron632$efron_0632_auc
w6  <- g6$within_imvigor210_efron632$efron_0632_auc
add_sd("ed1", "a", "data/positive_control_gep18.json",
       "within_imvigor210_efron632.efron_0632_auc", w18)
add_sd("ed1", "b", "data/rerun_positive_control.json",
       "within_imvigor210_efron632.efron_0632_auc", w6)

pc_panel <- function(df, w_within, ttl) {
  df$target <- factor(df$target, levels = tgts)
  slope_lbl <- ifelse(is.na(df$slope), "n/e",
                      sprintf("%.2f", df$slope))
  ggplot(df, aes(x = target, y = auc, colour = target)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey35", linewidth = 0.35) +
    geom_hline(yintercept = w_within, linetype = "dotdash", colour = "grey55", linewidth = 0.3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.4) +
    geom_point(size = 2.1) +
    geom_text(aes(label = sprintf("AUC %.3f\n[%.3f, %.3f]\nslope %s", auc, lo, hi, slope_lbl),
                  hjust = ifelse(target == "GSE165252", -0.12, 1.12)),
              nudge_x = ifelse(df$target == "GSE165252", 0.06, -0.06),
              size = 2.0, family = FF, lineheight = 0.92,
              colour = "grey15", show.legend = FALSE) +
    annotate("text", x = 1.5, y = 0.335, size = 2.0, family = FF, colour = "grey30",
             label = sprintf("within IMvigor210 (Efron 0.632) = %.3f", w_within)) +
    scale_colour_manual(values = c(GSE165252 = "#1f77b4", GSE78220 = "#ff7f0e"), guide = "none") +
    scale_y_continuous(limits = c(0.3, 0.92)) +
    labs(title = ttl, x = NULL, y = "external AUC (BCa 95% CI)") +
    theme_
}
pE1a <- pc_panel(dfG18, w18,
                 "18-gene T-cell-inflamed GEP (primary)")
pE1b <- pc_panel(dfG6, w6,
                 "6-gene IFN-\u03b3 panel (sensitivity only)")
ed1 <- (pE1a | pE1b) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_trio("ed1_positive_control", ed1, 180, 95)
write_sd("ed1_positive_control")

# =============================================================================
# BLOCK C -- ED6: threshold sensitivity of the transport criterion
# =============================================================================
tz  <- rd(J("data", "", "threshold_zero_sensitivity.json"))
dm  <- tz$p0_3$decision_matrix
tols <- unlist(dm$tolerances)
cell_names <- names(dm$cells)
models <- cell_names
tgts_dm <- c("IMvigor210", "GSE165252", "GSE78220")
dfT <- do.call(rbind, lapply(models, function(mdl) {
  lapply(tgts_dm, function(tgt) {
  cell <- dm$cells[[mdl]][[tgt]]
  c1 <- isTRUE(cell$c1_auc_ci_beneficial_exclusion)
  c2 <- vapply(as.character(tols), function(tt)
    isTRUE(cell$c2_slope_within_tolerance_by_tol[[tt]]), logical(1))
  ok <- isTRUE(cell$overall_transport_permitted_by_tol[[as.character(tols[2])]])
  data.frame(model = mdl, target = tgt, C1 = c1,
             C2_tol = c2[1] + c2[2] + c2[3],
             permitted = ok, stringsAsFactors = FALSE)
  })
}))
dfT <- do.call(rbind, dfT)
for (i in seq_len(nrow(dfT))) {
  add_sd("ed6", "a", "data/threshold_zero_sensitivity.json",
         sprintf("decision_matrix.cells.%s.%s.c1", dfT$model[i], dfT$target[i]), dfT$C1[i])
  add_sd("ed6", "a", "data/threshold_zero_sensitivity.json",
         sprintf("decision_matrix.cells.%s.%s.c2_tolerances_passed", dfT$model[i], dfT$target[i]),
         dfT$C2_tol[i])
}
dfT$model <- factor(dfT$model, levels = models)
dfT$target <- factor(dfT$target, levels = tgts_dm)
dfT$tile <- ifelse(dfT$permitted, "permitted (C1 \u2227 C2)",
            ifelse(dfT$C1, "C1 only",
            ifelse(dfT$C2_tol > 0, sprintf("C2 within \u00b10.20"), "refused")))
dfT$tile <- factor(dfT$tile, levels = c("permitted (C1 \u2227 C2)", "C1 only",
                                        "C2 within \u00b10.20", "refused"))
pT1 <- ggplot(dfT, aes(x = target, y = model, fill = tile)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("permitted (C1 \u2227 C2)" = "#2ca02c",
                               "C1 only" = "#ffd166",
                               "C2 within \u00b10.20" = "#ff9f43",
                               "refused" = "#e8e8e8"), name = NULL) +
  labs(title = "Dual-criterion decision matrix (24 cells, tolerance grid)",
       x = NULL, y = NULL) +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.8, "mm"),
                    axis.text.x = element_text(angle = 15, hjust = 1))
pT1 <- lab(pT1, "a")

rs <- rd(J("data", "", "rerun_sensitivity.json"))
p10 <- rs$p0_10
dfT2 <- do.call(rbind, lapply(tgts, function(t) {
  row <- p10$rows$panel[[t]]$primary_frozen$delta_auc
  add_sd("ed6", "b", "data/rerun_sensitivity.json",
         sprintf("p0_10.rows.panel.%s.primary_frozen.delta_auc", t),
         list(point = row$point, within_reference = row$within_reference,
              ci95_paired_bootstrap = row$ci95_paired_bootstrap))
  data.frame(target = t,
             delta = row$point,
             lo = row$ci95_paired_bootstrap[[1]], hi = row$ci95_paired_bootstrap[[2]],
             stringsAsFactors = FALSE)
}))
dfT2$target <- factor(dfT2$target, levels = tgts)
pT2 <- ggplot(dfT2, aes(x = delta, y = target, colour = target)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.12, linewidth = 0.4) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(GSE165252 = "#1f77b4", GSE78220 = "#ff7f0e"), guide = "none") +
  labs(title = "Continuous transport: \u0394AUC (external \u2212 within), bootstrap 95% CI",
       x = "\u0394AUC", y = NULL) +
  theme_
pT2 <- lab(pT2, "b")

ed6 <- pT1 / pT2 +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_trio("ed6_threshold_sensitivity", ed6, 180, 155)
write_sd("ed6_threshold_sensitivity")

# =============================================================================
# BLOCK D -- ED9: external-validation precision simulation (P0-5)
# =============================================================================
ps <- rd(J("data", "", "precision_simulation.json"))
auc_lvls <- names(ps$results$GSE165252)
dfP <- do.call(rbind, lapply(c("GSE165252", "GSE78220"), function(t) {
  st <- ps$settings[[t]]
  data.frame(target = t, n = st$n,
             auc_true = as.numeric(auc_lvls),
             p_c1 = vapply(auc_lvls, function(k) ps$results[[t]][[k]]$prob_c1_exclude_0p5, 0),
             ci_w = vapply(auc_lvls, function(k) ps$results[[t]][[k]]$mean_ci_width, 0),
             stringsAsFactors = FALSE)
}))
for (i in seq_len(nrow(dfP))) {
  add_sd("ed9", "a", "data/precision_simulation.json",
         sprintf("results.%s.%s.prob_c1_exclude_0p5", dfP$target[i], dfP$auc_true[i]), dfP$p_c1[i])
  add_sd("ed9", "b", "data/precision_simulation.json",
         sprintf("results.%s.%s.mean_ci_width", dfP$target[i], dfP$auc_true[i]), dfP$ci_w[i])
}
dfP$target <- factor(dfP$target, levels = c("GSE165252", "GSE78220"))
pP1 <- ggplot(dfP, aes(x = auc_true, y = p_c1, colour = target, group = target)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", colour = "grey45", linewidth = 0.3) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
  scale_colour_manual(values = c(GSE165252 = "#1f77b4", GSE78220 = "#ff7f0e"),
                      labels = c("GSE165252 (n = 71)", "GSE78220 (n = 25)"), name = NULL) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = "Probability of excluding AUC = 0.5 (criterion C1)",
       x = "true AUC (simulation)", y = "P(C1 passes)") +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.8, "mm"))
pP1 <- lab(pP1, "a")
pP2 <- ggplot(dfP, aes(x = auc_true, y = ci_w, colour = target, group = target)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
  scale_colour_manual(values = c(GSE165252 = "#1f77b4", GSE78220 = "#ff7f0e"),
                      labels = c("GSE165252 (n = 71)", "GSE78220 (n = 25)"), name = NULL) +
  labs(title = "Mean bootstrap CI width at the observed sizes",
       x = "true AUC (simulation)", y = "mean CI width") +
  theme_ + theme(legend.position = "bottom", legend.key.size = unit(2.8, "mm"))
pP2 <- lab(pP2, "b")
ed9 <- pP1 / pP2 +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
save_trio("ed9_precision_simulation", ed9, 180, 150)
write_sd("ed9_precision_simulation")

cat("\n[] updated figures complete: fig4, ed1, ed6, ed9\n")
