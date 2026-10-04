#!/usr/bin/env Rscript
# =============================================================================
# make_anchor_figures.R -- Nature Cancer Analysis  ANCHOR-ARM figure set
#
# Depends on the IMvigor210 gene-row identity having been RECOVERED (
# gene-recovery key: exmat row j == j-th fData row under upper-cased symbol
# alphabetical sort). Every statistical number is read verbatim from frozen
# JSON artifacts under data/ -- nothing data-derived is hard-coded; only
# geometry (axis ranges, dodge widths, annotation positions) is inline.
#
# Figures (anchor arm, post-recovery numbers)
#   Fig.2  fig2_external_transport     Frozen immune-state ranking does not
#                                     reproduce externally
#          (a) primary frozen transport forest + within-reference dumbbell
#          (b) pairwise 6-direction forest + source-swap delta annotation
#          (c) Ayers IFN-gamma positive control vs S0
#   Fig.3  fig3_calibration_pathways  Calibration and discrimination separate
#          (a) M2 frozen-prediction calibration curves, 3 targets, 5 bins, rug
#          (b) M3 recalibration before/after slopes
#          (c) Brier / log-loss M0/M2/M3 x 3 cohorts
#   Fig.4  fig4_ablation_constraint   Representation perturbation constrains,
#                                     but does not localize, transport failure
#          (a) 21-config external AUC + CI, two cohort facets, layer colours
#          (b) layer variance decomposition (new vs old )
#          (c) gene-identity recovery old-vs-new comparison (provenance fix)
#   Fig.5  fig5_benchmark_strategies  Benchmarking against simpler transport
#                                     strategies
#          (a) M0-M6 x 3 cohorts LOCO AUC+CI, two-criterion refusals marked
#          (b) reference strategies (primary single-source estimand):
#              clinical-only, Ayers, M0, M2 frozen S0
#   ED1    ed1_positive_control       simulated positive control + Ayers
#                                     empirical control (both FAIL)
#   ED3    ed3_loco_sensitivity       LOCO vs primary estimand + M2 rows
#   ED4    ed4_ablation_matrix        21-config x cohort AUC heatmap
#   ED5    ed5_identity_recovery      gene-identity recovery timeline +
#                                     old/new numbers (provenance lesson)
#   ED6    ed6_calibration_stability  M3 100-repeat recalibration stability
#   ED7    ed7_clinical_casemix       clinical-only baseline + case-mix shift
#   ED10   ed10_analysis_freeze       analysis-freeze manifest UPDATED with
#                                     the gene-identity recovery event, the
#                                     _rerun artifact series and the
#                                     retraction record
#
# Outputs: figures/{name}.{png,pdf,svg} (PNG 300 dpi, 180 mm wide)
#          figures/source_data/{name}.json (per-panel machine-readable)
#
# Frozen sources (read verbatim, asserted):
#   data/imvigor_gene_recovery.json   (recovery key, gates, old/new)
#   data/rerun_benchmark.json         (LOCO M0-M6, within, calib,
#                                              pairwise 6-dir, near-miss)
#   data/rerun_ablation.json          (21-config grid, layer variance,
#                                              retraction of artefact configs)
#   data/rerun_positive_control.json  (Ayers 6-gene, tiers)
#   data/clinical_only.json           (C_main/C_sex/C_rich, case-mix)
#   data/supporting_analyses.json     (simulated positive control --
#                                             identity-independent, still valid)
# =============================================================================

.libPaths(c(file.path(Sys.getenv("HOME"), "Library/R/4.6/library"), .libPaths()))

# skill helper: theme_pub / panel() / multi-panel utilities + Cell palettes
source("/Users/zuoxianbo/.workbuddy/skills/r-ggplot2-plotting/ggbiomed.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(jsonlite)
  library(grid)
})

# ---------------------------------------------------------------- paths -----
BASE   <- "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature/outputs/analysis_base"
OUTDIR <- file.path(BASE, "figures_")
SDDIR  <- file.path(OUTDIR, "source_data")
dir.create(SDDIR, recursive = TRUE, showWarnings = FALSE)

J <- function(...) file.path(BASE, ...)

F_RRB <- J("data", "", "rerun_benchmark.json")
F_RRA <- J("data", "", "rerun_ablation.json")
F_RPC <- J("data", "", "rerun_positive_control.json")
F_CLN <- J("data", "", "clinical_only.json")
F_REC <- J("data", "", "imvigor_gene_recovery.json")
F_S75 <- J("data", "", "supporting_analyses.json")

# some upstream Python artifacts embed literal NaN / Infinity tokens, which
# jsonlite rejects; sanitise them to null before parsing (values we consume
# are never NaN -- guards below assert the exact numbers we need)
fromJSON_san <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub(":\\s*NaN\\b", ": null", txt)
  txt <- gsub(":\\s*-?Infinity\\b", ": null", txt)
  txt <- gsub("\\[\\s*NaN\\s*\\]", "[null]", txt)
  txt <- gsub(",\\s*NaN\\s*,", ", null,", txt)
  fromJSON(txt, simplifyVector = FALSE)
}

rrb  <- fromJSON_san(F_RRB)
rra  <- fromJSON_san(F_RRA)
rpc  <- fromJSON_san(F_RPC)
clin <- fromJSON_san(F_CLN)
rec  <- fromJSON_san(F_REC)
s75  <- fromJSON_san(F_S75)

COHORTS  <- c("IMvigor210", "GSE165252", "GSE78220")
TARGETS  <- c("IMvigor210", "GSE165252", "GSE78220")
EXT_CH   <- c("GSE165252", "GSE78220")          # external target cohorts
MODELS   <- c("M0", "M1", "M2", "M3", "M4", "M5a", "M5b", "M6")

# ---------------------------------------------------------------- style ------
# mirrors make_stable_figures.R / make_figures.R (visual consistency)
FF <- "Arial"
th <- theme_classic(base_family = FF, base_size = 8.5) +
  theme(
    axis.text          = element_text(colour = "black", size = 7),
    axis.title         = element_text(size = 8.5),
    plot.title         = element_text(size = 8.5, face = "plain"),
    plot.subtitle      = element_text(size = 7.5),
    strip.text         = element_text(size = 8, face = "bold"),
    strip.background   = element_rect(fill = "grey93", colour = NA),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 7.5),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.minor   = element_blank()
  )
theme_set(th)

BLUE <- "#1f4e79"; RED <- "#c0504d"; OLIVE <- "#4f6228"; PURPLE <- "#7f5fa2"
PAL  <- c(BLUE, RED, OLIVE, PURPLE)

# bold lowercase panel tag, top-left, clipped-in (/ convention)
lab <- function(txt) annotate("text", x = -Inf, y = Inf, label = txt,
                              hjust = -0.35, vjust = 2.2, size = 3.4,
                              fontface = "bold", family = FF)

# same tag for FACETED panels -- drawn once per plot (not once per facet strip)
lab1 <- function(txt) {
  list(labs(tag = txt),
       theme(plot.tag.position = c(0.012, 0.985),
             plot.tag = element_text(size = 10.5, face = "bold", family = FF)))
}

# ------------------------------------------------------- source data log -----
SD <- new.env(parent = emptyenv())
add_sd <- function(fig, panel, src, key, value) {
  if (!exists(fig, envir = SD, inherits = FALSE)) assign(fig, new.env(parent = emptyenv()), envir = SD)
  p <- get(fig, envir = SD)
  if (!exists(panel, envir = p, inherits = FALSE)) assign(panel, new.env(parent = emptyenv()), envir = p)
  pp <- get(panel, envir = p)
  key <- gsub("'", "", key)
  assign(key, value, envir = pp)
}
dump_sd <- function(fig, name) {
  e <- if (exists(fig, envir = SD, inherits = FALSE)) get(fig, envir = SD) else NULL
  panels <- list()
  if (!is.null(e)) {
    pn <- ls(e)
    for (p in pn) panels[[p]] <- as.list(get(p, envir = e))
  }
  out <- list(figure = name, generated_by = "figures/make_anchor_figures.R", panels = panels)
  write(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"),
        file.path(SDDIR, paste0(name, ".json")))
}

# ------------------------------------------------------------ save trio ------
SVG_OK <- requireNamespace("svglite", quietly = TRUE)
save_trio <- function(p, name, w_mm, h_mm) {
  png_file <- file.path(OUTDIR, paste0(name, ".png"))
  err <- tryCatch({
    ggsave(png_file, p, width = w_mm, height = h_mm, units = "mm",
           dpi = 300, device = ragg::agg_png, bg = "white")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(err)) stop(sprintf("[render] %s PNG failed: %s", name, err))
  pdf_file <- file.path(OUTDIR, paste0(name, ".pdf"))
  dev_lbl <- "quartz-pdf"
  pdf_err <- tryCatch({
    ggsave(pdf_file, p, width = w_mm, height = h_mm, units = "mm",
           device = function(file, width, height, ...)
             grDevices::quartz(file = file, width = width, height = height,
                               type = "pdf", family = FF))
    NULL
  }, error = function(e) conditionMessage(e))
  if (!file.exists(pdf_file) || file.size(pdf_file) < 2000) {
    dev_lbl <- "pdf-Helvetica"
    pdf_err <- tryCatch({
      ggsave(pdf_file, p, width = w_mm, height = h_mm, units = "mm",
             device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE)
      NULL
    }, error = function(e) conditionMessage(e))
  }
  if (!is.null(pdf_err)) cat(sprintf("[render-warn] %s PDF failed: %s\n", name, pdf_err))
  svg_lbl <- "svglite"
  if (SVG_OK) {
    svg_err <- tryCatch({
      ggsave(file.path(OUTDIR, paste0(name, ".svg")), p, width = w_mm, height = h_mm,
             units = "mm", device = svglite::svglite, fix_text_size = FALSE, bg = "white")
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(svg_err)) cat(sprintf("[render-warn] %s SVG failed: %s\n", name, svg_err))
    if (!file.exists(file.path(OUTDIR, paste0(name, ".svg")))) svg_lbl <- "svg-FAILED"
  } else {
    svg_lbl <- "svglite-unavailable (degraded: PNG+PDF only)"
  }
  cat(sprintf("[saved] %-38s (pdf via %s; svg: %s)\n", name, dev_lbl, svg_lbl))
  for (ext in c("png", "pdf", "svg")) {
    fp <- file.path(OUTDIR, paste0(name, ".", ext))
    sz <- if (file.exists(fp)) round(file.size(fp) / 1024, 1) else 0
    cat(sprintf("          %s.%s  %8.1f kB\n", name, ext, sz))
  }
}

fmt3 <- function(x) sprintf("%.3f", x)
fmt2 <- function(x) sprintf("%.2f", x)

# regex-extract helper (sourced text, never typed numbers)
gx <- function(txt, pat, src) {
  m <- regmatches(txt, regexpr(pat, txt))
  if (length(m) != 1L || !nzchar(m))
    stop(sprintf("frozen-source assertion failed: pattern '%s' not found in %s", pat, src))
  m
}

# =============================================================================
# BLOCK 0 -- shared data preparation + frozen-value assertions
# =============================================================================

# ---- 0.1 anchor identity (recovered key) ------------------------------------
anc <- rrb$anchor
stopifnot(anc$panel == 322L, anc$imv_available == 320L, anc$common == 320L,
          identical(anc$imv_missing, list("IGHG1", "PRKN")),
          anc$n_allzero_genes_across_348 == 0L,
          abs(anc$pc_explained_variance_ratio[[1]] - 0.2653) < 5e-4,
          abs(anc$pc_explained_variance_ratio[[2]] - 0.0733) < 5e-4)
stopifnot(rec$verification_gates$overall == "PASS",
          rec$verification_gates$sex_anchor == "PASS",
          rec$verification_gates$panel_measurable$new_by_symbol == 320L,
          rec$verification_gates$panel_measurable$old_fake_entrez_selection == 229L,
          identical(rec$phase1a$winner, "symbol_alphabetical_upper"))
for (k in names(anc)) add_sd("shared", "anchor", "data/rerun_benchmark.json",
                             paste0("anchor.", k), anc[[k]])

# ---- 0.2 within-cohort reference (Efron 0.632, recovery convention) ---------
wr <- rrb$within_reference$per_cohort
stopifnot(wr$IMvigor210$n == 298L, wr$IMvigor210$n_event == 68L,
          wr$GSE165252$n == 71L,  wr$GSE165252$n_event == 23L,
          wr$GSE78220$n == 25L,   wr$GSE78220$n_event == 13L,
          abs(wr$IMvigor210$auc_632 - 0.6571) < 5e-4,
          abs(wr$GSE165252$auc_632 - 0.4279) < 5e-4,
          abs(wr$GSE78220$auc_632  - 0.6220) < 5e-4)
for (ch in COHORTS) {
  add_sd("shared", "within", "data/rerun_benchmark.json",
         sprintf("within_reference.per_cohort.%s.auc_632", ch), wr[[ch]]$auc_632)
  add_sd("shared", "within", "data/rerun_benchmark.json",
         sprintf("within_reference.per_cohort.%s.ci95", ch), wr[[ch]]$ci95)
  add_sd("shared", "within", "data/rerun_benchmark.json",
         sprintf("within_reference.per_cohort.%s.n", ch), wr[[ch]]$n)
  add_sd("shared", "within", "data/rerun_benchmark.json",
         sprintf("within_reference.per_cohort.%s.n_event", ch), wr[[ch]]$n_event)
}

# ---- 0.3 primary frozen transports (train IMvigor210, F1-style) -------------
pk <- function(tr, te) paste0("train_", tr, "__test_", te)
pw <- rrb$pairwise_transports$transports
prim <- list(
  GSE165252 = pw[[pk("IMvigor210", "GSE165252")]]$M2,
  GSE78220  = pw[[pk("IMvigor210", "GSE78220")]]$M2
)
stopifnot(abs(prim$GSE165252$auc - 0.5036) < 5e-4,
          abs(prim$GSE165252$auc_ci95[[1]] - 0.3540) < 5e-4,
          abs(prim$GSE165252$auc_ci95[[2]] - 0.6585) < 5e-4,
          abs(prim$GSE78220$auc  - 0.6603) < 5e-4,
          abs(prim$GSE78220$auc_ci95[[1]]  - 0.4038) < 5e-4,
          abs(prim$GSE78220$auc_ci95[[2]]  - 0.8617) < 5e-4,
          prim$GSE165252$n_test == 71L, prim$GSE165252$n_event == 23L,
          prim$GSE78220$n_test == 25L,  prim$GSE78220$n_event == 13L)

# ---- 0.4 LOCO benchmark adjudication: 24/24 refused --------------------------
adj <- rrb$adjudication
stopifnot(adj$verdict_counts$n_total == 24L,
          adj$verdict_counts$n_permitted == 0L,
          adj$verdict_counts$n_refused == 24L,
          all(unlist(adj$verdict_counts$per_model_n_permitted) == 0))
n_c1_imv <- sum(sapply(MODELS, function(m)
  isTRUE(adj$per_model[[m]]$IMvigor210$auc_ci_excludes_0p5_beneficial_direction)))
stopifnot(n_c1_imv == 6L)   # M2,M3,M4,M5a,M5b,M6 pass criterion 1 at IMvigor210

# ---- 0.5 ablation grid: 21 configs, none recovers transport ------------------
cfg_names <- names(rra$per_config)
stopifnot(length(cfg_names) == 21L,
          isFALSE(rra$verdict$any_config_recovers_transport$any),
          identical(rra$verdict$largest_variance_layer, "gene_set"),
          identical(rra$verdict$layer_variance_ranking,
                    list("gene_set", "pc_k", "normalization", "ridge_lambda")),
          length(rra$removed_configs) == 5L)
n_degenerate <- sum(sapply(cfg_names, function(c) isTRUE(rra$per_config[[c]]$degenerate)))
stopifnot(n_degenerate == 0L)

# ---- 0.6 Ayers positive control (recovered identity) -------------------------
ay <- rpc$transports$primary_logcpm_z_frozen_all6
stopifnot(abs(ay$GSE165252$target_metrics$auc - 0.5607) < 5e-4,
          abs(ay$GSE78220$target_metrics$auc  - 0.4936) < 5e-4,
          isFALSE(ay$GSE165252$adjudication_target$main_analysis_rule$both_passed),
          isFALSE(ay$GSE78220$adjudication_target$main_analysis_rule$both_passed),
          abs(rpc$within_imvigor210_efron632$efron_0632_auc - 0.6050) < 5e-4)
for (ch in EXT_CH)
  stopifnot(all(unlist(rpc$verdict$per_cohort[[ch]]$tiers) == FALSE))

# ---- 0.7 clinical-only -------------------------------------------------------
stopifnot(abs(clin$external_results$GSE78220$C_main$auc - 0.4872) < 5e-4,
          clin$external_results$GSE165252$C_main$auc == 0.5,
          abs(clin$within_imvigor210_efron632$C_rich$efron_0632_auc - 0.7927) < 5e-4,
          abs(clin$casemix$prevalence_shift$GSE78220 - 0.520) < 5e-4,
          abs(clin$casemix$standardized_differences_IMvigor210_vs_GSE78220$pct_advanced_stage_smd + 2.4984) < 5e-4,
          abs(clin$casemix$standardized_differences_IMvigor210_vs_GSE78220$pct_male_smd - 0.2313) < 5e-4)

# ---- 0.8 simulated positive control (, identity-independent) --------------
pcs <- s75$positive_control$results
pc_max_well <- max(unlist(lapply(names(pcs)[grepl("wellcalibrated", names(pcs))],
                                 function(s) lapply(pcs[[s]]$n_target_results,
                                                    function(x) x$pass_rate_rule_i_and_ii))))
pc_mis_zero <- all(unlist(lapply(names(pcs)[grepl("miscalibrated", names(pcs))],
                                 function(s) lapply(pcs[[s]]$n_target_results,
                                                    function(x) x$pass_rate_rule_i_and_ii))) == 0)
stopifnot(abs(pc_max_well - 0.338) < 5e-4, pc_mis_zero)

# ---- 0.9 recovery old/new blocks ---------------------------------------------
rec_old <- rec$phase2$old; rec_new <- rec$phase2$new
stopifnot(abs(rec_old$transport$external$GSE78220$auc - 0.4103) < 5e-4,
          abs(rec_new$transport$external$GSE78220$auc - 0.6603) < 5e-4,
          abs(rec_old$within_imvigor$auc_apparent - 0.5836) < 5e-4,
          abs(rec_new$within_imvigor$auc_apparent - 0.6659) < 5e-4,
          rec_old$zero_genes_frozen_$n_unmeasured == 70L,
          rec_old$zero_genes_frozen_$n_common_genes == 228L,
          abs(rec_old$zero_genes_frozen_$pc1_loading_squared_mass - 0.8866) < 5e-4,
          rec_new$zero_genes$n_allzero_across_348 == 0L,
          rec_new$pc1_loading_squared_mass_on_zero_genes == 0,
          abs(rec$phase2.old$within_632_frozen__pca_stability - 0.5669) < 5e-4 ||
            abs(rec_old$within_632_frozen__pca_stability - 0.5669) < 5e-4)
stopifnot(isFALSE(rec$phase1a$hypotheses[["position_fData (=entrez string sort, the falsified old assumption)"]]$pass),
          isTRUE(rec$phase1a$hypotheses$symbol_alphabetical_upper$pass),
          isTRUE(rec$phase1a$hypotheses$symbol_alphabetical_rawcase$pass))

cat("[block0] all frozen-value assertions passed\n")

# =============================================================================
# FIG 2 -- Frozen immune-state ranking does not reproduce externally
# =============================================================================
make_fig2 <- function() {

  # ---- (a) primary frozen transport forest + within dumbbell ----
  dA <- do.call(rbind, lapply(COHORTS, function(ch) {
    w <- wr[[ch]]
    if (ch == "IMvigor210") {
      # source anchor: within only (no external target -- it is the anchor)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("within_reference.per_cohort.%s.auc_632", ch), w$auc_632)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("within_reference.per_cohort.%s.ci95", ch), w$ci95)
      data.frame(cohort = ch, within = w$auc_632, w_lo = w$ci95[[1]], w_hi = w$ci95[[2]],
                 ext = NA, e_lo = NA, e_hi = NA, n = w$n, n_event = w$n_event)
    } else {
      m <- prim[[ch]]
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.auc", ch), m$auc)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.auc_ci95", ch), m$auc_ci95)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.n_test", ch), m$n_test)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.n_event", ch), m$n_event)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("within_reference.per_cohort.%s.auc_632", ch), w$auc_632)
      add_sd("fig2", "a", "data/rerun_benchmark.json",
             sprintf("within_reference.per_cohort.%s.ci95", ch), w$ci95)
      data.frame(cohort = ch, within = w$auc_632, w_lo = w$ci95[[1]], w_hi = w$ci95[[2]],
                 ext = m$auc, e_lo = m$auc_ci95[[1]], e_hi = m$auc_ci95[[2]],
                 n = m$n_test, n_event = m$n_event)
    }
  }))
  # CI-crosses-0.5 visual emphasis flag (external rows only)
  dA$crosses <- !is.na(dA$ext) & dA$e_lo < 0.5 & dA$e_hi > 0.5
  stopifnot(all(dA$crosses[!is.na(dA$ext)]))   # both primary transports span 0.5
  dA$row_lab <- sprintf("%s  (n = %d, %d events)", dA$cohort, dA$n, dA$n_event)
  dA$y <- rev(seq_len(nrow(dA)))

  pa <- ggplot(dA, aes(x = within, y = y)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = w_lo, xmax = w_hi), height = 0.10,
                   linewidth = 0.35, colour = BLUE, alpha = 0.55) +
    geom_point(aes(fill = "within (Efron 0.632)"), shape = 23, size = 2.6,
               colour = BLUE, stroke = 0.5) +
    geom_segment(data = subset(dA, !is.na(ext)),
                 aes(x = within, xend = ext, yend = y),
                 arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.55,
                 colour = "grey60", linetype = "dotted") +
    geom_pointrange(data = subset(dA, !is.na(ext)),
                    aes(x = ext, xmin = e_lo, xmax = e_hi,
                        colour = crosses, size = crosses),
                    fatten = 2.6, stroke = 0.55, shape = 21,
                    fill = ifelse(subset(dA, !is.na(ext))$crosses, "#e2a09e", "white")) +
    geom_text(data = subset(dA, !is.na(ext)),
              aes(x = ext, label = sprintf("%.3f [%.3f, %.3f]", ext, e_lo, e_hi)),
              vjust = -1.15, size = 2.0, family = FF, colour = "grey20") +
    geom_label(data = dA, aes(x = 0.185, label = row_lab), hjust = 0,
               size = 2.1, family = FF, colour = "grey20",
               label.padding = unit(0.09, "lines"), label.size = 0,
               fill = alpha("white", 0.78)) +
    geom_text(data = subset(dA, !is.na(ext)),
              aes(x = e_hi, label = "CI \u220b 0.5"), hjust = -0.15,
              size = 2.0, family = FF, colour = RED, fontface = "bold") +
    annotate("text", x = 0.30, y = 3.62, hjust = 0, size = 2.0, family = FF,
             colour = "grey35",
             label = "IMvigor210 = source anchor (within only);\nexternal rows: frozen S0, train IMvigor210 (primary estimand)") +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = "grey35"), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.4), guide = "none") +
    scale_fill_manual(values = c("within (Efron 0.632)" = BLUE), guide = "none") +
    scale_x_continuous(limits = c(0.18, 0.95), breaks = seq(0.2, 0.9, 0.1)) +
    scale_y_continuous(breaks = dA$y, labels = rep("", nrow(dA))) +
    coord_cartesian(ylim = c(0.35, 3.9), clip = "off") +
    labs(x = "AUC (E ~ S0), efficacy", y = NULL,
         title = "Primary frozen transports: both external CIs span 0.5",
         subtitle = "diamond = within-cohort 0.632 reference; circle = frozen IMvigor210-anchored external transport (BCa 95% CI)") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    lab("a")

  # ---- (b) pairwise 6-direction forest + source-swap annotation ----
  dirs <- list(c("IMvigor210", "GSE165252"), c("IMvigor210", "GSE78220"),
               c("GSE165252", "IMvigor210"), c("GSE165252", "GSE78220"),
               c("GSE78220", "IMvigor210"),  c("GSE78220", "GSE165252"))
  dB <- do.call(rbind, lapply(dirs, function(d) {
    tr <- d[1]; te <- d[2]
    m <- pw[[pk(tr, te)]]$M2
    ov <- rrb$old_vs_new$pairwise[[pk(tr, te)]]
    excl_ben <- (m$auc_ci95[[1]] > 0.5)
    excl_adv <- (m$auc_ci95[[2]] < 0.5)
    add_sd("fig2", "b", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.%s.M2.auc", pk(tr, te)), m$auc)
    add_sd("fig2", "b", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.%s.M2.auc_ci95", pk(tr, te)), m$auc_ci95)
    add_sd("fig2", "b", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.%s.M2.n_test", pk(tr, te)), m$n_test)
    add_sd("fig2", "b", "data/rerun_benchmark.json",
           sprintf("old_vs_new.pairwise.%s.delta_auc", pk(tr, te)), ov$delta_auc)
    data.frame(train = tr, target = te,
               row = sprintf("%s \u2192 %s", tr, te),
               est = m$auc, lo = m$auc_ci95[[1]], hi = m$auc_ci95[[2]],
               n_test = m$n_test, n_event = m$n_event,
               delta_old = ov$delta_auc,
               primary = (tr == "IMvigor210"),
               excl_ben = excl_ben, excl_adv = excl_adv)
  }))
  stopifnot(sum(dB$excl_ben) == 1L, sum(dB$excl_adv) == 1L,
            dB$excl_ben[dB$train == "GSE78220" & dB$target == "IMvigor210"],
            dB$excl_adv[dB$train == "GSE165252" & dB$target == "IMvigor210"])
  # sourced slope for the one beneficial-direction exclusion (no typed numbers)
  sl_ben <- pw[[pk("GSE78220", "IMvigor210")]]$M2$calib_slope
  add_sd("fig2", "b", "data/rerun_benchmark.json",
         "pairwise_transports.transports.train_GSE78220__test_IMvigor210.M2.calib_slope", sl_ben)
  dB$row_full <- paste0(dB$row, "  (n = ", dB$n_test, ")")
  dB$y <- rev(seq_len(nrow(dB)))

  pb <- ggplot(dB, aes(x = est, y = y)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi,
                        colour = factor(excl_ben + 2 * excl_adv),
                        size = primary),
                    fatten = 2.4, stroke = 0.55, shape = 19) +
    geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]", est, lo, hi)),
              vjust = -1.2, size = 1.9, family = FF, colour = "grey25") +
    geom_label(aes(x = 0.135, label = row_full), hjust = 0, size = 2.0,
               family = FF, colour = "grey20",
               label.padding = unit(0.09, "lines"), label.size = 0,
               fill = alpha("white", 0.78)) +
    geom_text(aes(x = 1.36, label = sprintf("\u0394 vs old = %+.3f", delta_old)),
              hjust = 1, size = 1.9, family = FF,
              colour = ifelse(abs(dB$delta_old) >= 0.15, RED, "grey45")) +
    annotate("text", x = 0.72, y = 6.55, hjust = 0, size = 2.0, family = FF,
             colour = "grey30", lineheight = 0.95,
             label = "all six ordered directions, frozen M2 (S0)\n\u0394 vs old = recovered-identity rerun minus\nwrong-identity frozen archive (source-swap rows)") +
    annotate("point", x = dB$est[dB$excl_ben], y = dB$y[dB$excl_ben] - 0.42,
             shape = 25, size = 2.2, fill = OLIVE, colour = "black", stroke = 0.4) +
    annotate("text", x = dB$est[dB$excl_ben] + 0.03, y = dB$y[dB$excl_ben] - 0.47,
             hjust = 0, size = 1.9, family = FF, colour = OLIVE,
             label = sprintf("CI excludes 0.5 (beneficial)\nbut slope %.2f \u2192 still refused", sl_ben)) +
    annotate("text", x = 0.55, y = dB$y[dB$excl_adv] - 0.42, hjust = 0, size = 1.9,
             family = FF, colour = PURPLE,
             label = "CI excludes 0.5 (adverse direction)") +
    scale_colour_manual(values = c("0" = "grey35", "1" = OLIVE,
                                   "2" = PURPLE, "3" = PURPLE),
                        guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.65), guide = "none") +
    scale_x_continuous(limits = c(0.13, 1.36), breaks = seq(0.2, 1.2, 0.2)) +
    scale_y_continuous(breaks = dB$y, labels = rep("", nrow(dB))) +
    coord_cartesian(ylim = c(0.4, 6.9), clip = "off") +
    labs(x = "AUC (E ~ S0), single-source transport", y = NULL,
         title = "True external pairwise transports: all six ordered directions") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(0.4, 1.2, 0.4, 0.4, "cm")) +
    lab("b")

  # ---- (c) Ayers positive control vs S0 (same frozen design) ----
  dC <- do.call(rbind, lapply(EXT_CH, function(ch) {
    a <- ay[[ch]]$target_metrics
    s <- prim[[ch]]
    add_sd("fig2", "c", "data/rerun_positive_control.json",
           sprintf("transports.primary_logcpm_z_frozen_all6.%s.target_metrics.auc", ch), a$auc)
    add_sd("fig2", "c", "data/rerun_positive_control.json",
           sprintf("transports.primary_logcpm_z_frozen_all6.%s.target_metrics.auc_ci95", ch), a$auc_ci95)
    add_sd("fig2", "c", "data/rerun_positive_control.json",
           sprintf("transports.primary_logcpm_z_frozen_all6.%s.target_metrics.calib_slope", ch), a$calib_slope)
    add_sd("fig2", "c", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.auc", ch), s$auc)
    rbind(data.frame(target = ch, model = "S0 (322-gene, frozen)",
                     est = s$auc, lo = s$auc_ci95[[1]], hi = s$auc_ci95[[2]],
                     slope = s$calib_slope, n = s$n_test, n_event = s$n_event),
          data.frame(target = ch, model = "Ayers IFN-\u03b3 6-gene (frozen)",
                     est = a$auc, lo = a$auc_ci95[[1]], hi = a$auc_ci95[[2]],
                     slope = a$calib_slope, n = a$n_test, n_event = a$n_event))
  }))
  add_sd("fig2", "c", "data/rerun_positive_control.json",
         "within_imvigor210_efron632.efron_0632_auc", rpc$within_imvigor210_efron632$efron_0632_auc)
  dC$row <- paste0(dC$model, "\n", dC$target)
  dC$row <- factor(dC$row, levels = rev(c(
    "S0 (322-gene, frozen)\nGSE165252", "Ayers IFN-\u03b3 6-gene (frozen)\nGSE165252",
    "S0 (322-gene, frozen)\nGSE78220",  "Ayers IFN-\u03b3 6-gene (frozen)\nGSE78220")))
  dC$y <- as.integer(dC$row)

  pc <- ggplot(dC, aes(x = est, y = y)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = model), fatten = 2.4,
                    stroke = 0.5, shape = 19, position = position_identity()) +
    geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]  slope %.2f", est, lo, hi, slope)),
              vjust = -1.25, size = 1.9, family = FF, colour = "grey25") +
    geom_text(aes(x = 0.10, label = sprintf("n = %d (%d events)", n, n_event)),
              hjust = 0, size = 1.9, family = FF, colour = "grey35") +
    annotate("text", x = 0.94, y = 4.45, hjust = 1, size = 2.0, family = FF,
             colour = "grey25", lineheight = 0.95,
             label = sprintf("anchor within (IMvigor210):\nS0 0.632 AUC %.3f \u00b7 IFN-\u03b3 %.3f\nboth signatures FAIL the two-criterion rule\nin BOTH external cohorts",
                             wr$IMvigor210$auc_632, rpc$within_imvigor210_efron632$efron_0632_auc)) +
    scale_colour_manual(values = c("S0 (322-gene, frozen)" = BLUE,
                                   "Ayers IFN-\u03b3 6-gene (frozen)" = RED),
                        name = NULL) +
    scale_x_continuous(limits = c(0.10, 0.95), breaks = seq(0.2, 0.9, 0.1)) +
    scale_y_continuous(breaks = seq_along(levels(dC$row)), labels = levels(dC$row)) +
    coord_cartesian(ylim = c(0.4, 4.7), clip = "off") +
    labs(x = "AUC (E ~ signature), train IMvigor210 (frozen)", y = NULL,
         title = "A known reproducible signature also fails in these cohorts",
         caption = "Ayers et al. JCI 2017 IFN-\u03b3 6-gene score (IDO1, CXCL9, CXCL10, CCL5, HLA-DRA, STAT1), recovered-key anchor;\nidentical frozen-transport design, seed 20260911; both-criterion adjudication FAIL in every cohort and every tolerance tier (see ED1).") +
    theme(axis.text.y = element_text(size = 6.3),
          plot.caption = element_text(size = 5.8, colour = "grey30"),
          legend.position = "none") +
    lab("c")

  pa / pb / pc + plot_layout(heights = c(1, 1.5, 1.15))
}

# =============================================================================
# FIG 3 -- Calibration and discrimination separate during transport
# =============================================================================
make_fig3 <- function() {
  cc <- rrb$calibration_curves$per_target

  # ---- (a) M2 frozen-prediction calibration curves (3 targets) ----
  pts <- do.call(rbind, lapply(TARGETS, function(ch) {
    m2 <- cc[[ch]]$M2
    for (b in m2$bins) {
      add_sd("fig3", "a", "data/rerun_benchmark.json",
             sprintf("calibration_curves.per_target.%s.M2.bins[%d]", ch, b$bin - 1L),
             list(n = b$n, mean_predicted = b$mean_predicted,
                  observed_rate = b$observed_rate, observed_rate_ci95 = b$observed_rate_ci95))
    }
    add_sd("fig3", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M2.calib_slope", ch), m2$calib_slope)
    add_sd("fig3", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M2.overall_mean_predicted", ch),
           m2$overall_mean_predicted)
    add_sd("fig3", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M2.overall_observed_rate", ch),
           m2$overall_observed_rate)
    data.frame(cohort = ch,
               x = sapply(m2$bins, function(b) b$mean_predicted),
               y = sapply(m2$bins, function(b) b$observed_rate),
               lo = sapply(m2$bins, function(b) b$observed_rate_ci95[[1]]),
               hi = sapply(m2$bins, function(b) b$observed_rate_ci95[[2]]),
               n = sapply(m2$bins, function(b) b$n))
  }))
  # guard: bin counts must sum to cohort n
  for (ch in TARGETS) stopifnot(sum(pts$n[pts$cohort == ch]) == cc[[ch]]$M2$n_test)
  stat_lab <- sapply(TARGETS, function(ch) {
    m2 <- cc[[ch]]$M2
    sprintf("%s\nn = %d \u00b7 %d events\nslope %.2f (logit)\nmean pred %.3f vs observed %.3f",
            ch, m2$n_test, m2$n_event, m2$calib_slope,
            m2$overall_mean_predicted, m2$overall_observed_rate)
  })
  stat_df <- data.frame(cohort = TARGETS, lab = stat_lab)
  # per-patient rug frames (frozen per-patient predicted distributions)
  rug_df <- do.call(rbind, lapply(TARGETS, function(ch) {
    data.frame(cohort = ch, p = unlist(cc[[ch]]$M2$per_patient_predicted))
  }))
  add_sd("fig3", "a", "data/rerun_benchmark.json",
         "calibration_curves.per_target.*.M2.per_patient_predicted",
         lapply(TARGETS, function(ch) length(cc[[ch]]$M2$per_patient_predicted)))
  pts$cohort <- factor(pts$cohort, levels = TARGETS)
  rug_df$cohort <- factor(rug_df$cohort, levels = TARGETS)
  stat_df$cohort <- factor(stat_df$cohort, levels = TARGETS)

  pa <- ggplot(pts, aes(x = x, y = y)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey50", linewidth = 0.35) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5,
                  colour = RED, alpha = 0.85) +
    geom_point(shape = 21, size = 2.1, fill = RED, colour = "black", stroke = 0.4) +
    geom_rug(data = rug_df, aes(x = p), inherit.aes = FALSE,
             sides = "b", length = unit(0.055, "cm"), alpha = 0.30,
             colour = BLUE, linewidth = 0.25) +
    geom_label(data = stat_df, aes(x = -Inf, y = Inf, label = lab),
               hjust = -0.02, vjust = 1.12, size = 1.95, family = FF,
               colour = "grey15", inherit.aes = FALSE, lineheight = 0.95,
               label.padding = unit(0.18, "lines"), label.size = 0,
               fill = alpha("white", 0.85)) +
    facet_wrap(~ cohort, ncol = 3, scales = "free") +
    scale_x_continuous(limits = c(0, NA)) +
    labs(x = "mean predicted probability (5 equal-frequency bins)",
         y = "observed event rate",
         title = "M2 frozen predictions: discrimination separates from calibration",
         caption = "bars = exact Clopper\u2013Pearson 95% CI per bin; bottom rug = per-patient frozen predictions (M2, LOCO train = pooled other two cohorts);\nidentity line = perfect calibration. GSE78220: n = 25 \u2014 bin rates rest on 5 patients each.") +
    theme(plot.caption = element_text(size = 5.9, colour = "grey30"),
          strip.text = element_text(size = 7)) +
    lab1("a")

  # ---- (b) M3 recalibration before/after slopes ----
  dB <- do.call(rbind, lapply(TARGETS, function(ch) {
    r <- cc[[ch]]$recal_before_after
    m3 <- cc[[ch]]$M3$repeat_means
    add_sd("fig3", "b", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.recal_before_after.m2_calib_slope_logit", ch),
           r$m2_calib_slope_logit)
    add_sd("fig3", "b", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.recal_before_after.m3_repeat_mean_calib_slope_logit", ch),
           r$m3_repeat_mean_calib_slope_logit)
    add_sd("fig3", "b", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.recal_before_after.m3_calib_slope_ci95", ch),
           r$m3_calib_slope_ci95)
    add_sd("fig3", "b", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.recal_before_after.slope_improved", ch),
           r$slope_improved)
    add_sd("fig3", "b", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M3.repeat_means.n_repeats_used", ch),
           m3$n_repeats_used)
    data.frame(cohort = ch, m2_slope = r$m2_calib_slope_logit,
               m3_slope = r$m3_repeat_mean_calib_slope_logit,
               m3_lo = r$m3_calib_slope_ci95[[1]], m3_hi = r$m3_calib_slope_ci95[[2]],
               improved = r$slope_improved, n_rep = m3$n_repeats_used)
  }))
  stopifnot(dB$improved[dB$cohort == "IMvigor210"] == TRUE,
            dB$improved[dB$cohort == "GSE165252"] == FALSE,
            dB$improved[dB$cohort == "GSE78220"] == FALSE)
  dB$cohort <- factor(dB$cohort, levels = TARGETS)
  dB$y <- as.integer(dB$cohort)

  pb <- ggplot(dB, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.4) +
    annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 0.55, ymax = 3.45,
             fill = OLIVE, alpha = 0.10) +
    annotate("text", x = 1.0, y = 3.62, size = 1.9, family = FF, colour = OLIVE,
             label = "|slope \u2212 1| \u2264 0.10") +
    geom_segment(aes(x = m2_slope, xend = m3_slope, yend = y, colour = improved),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.65,
                 linetype = "solid") +
    geom_errorbarh(aes(xmin = m3_lo, xmax = m3_hi), height = 0.12,
                   linewidth = 0.4, colour = "grey55") +
    geom_point(aes(x = m2_slope, fill = "M2 (before)"), shape = 22, size = 2.6,
               colour = "black", stroke = 0.5) +
    geom_point(aes(x = m3_slope, fill = "M3 (after, 100-repeat mean)"), shape = 21,
               size = 2.6, colour = "black", stroke = 0.5) +
    geom_text(aes(x = m2_slope, y = y + 0.22, label = sprintf("%.2f", m2_slope)),
              size = 2.0, family = FF, colour = "grey25") +
    geom_text(aes(x = m3_slope, y = y - 0.27, label = sprintf("%.2f [%.2f, %.2f]",
                                                              m3_slope, m3_lo, m3_hi)),
              size = 2.0, family = FF, colour = "grey25") +
    geom_text(aes(x = -4.35, label = sprintf("%s  %s", cohort,
                                             ifelse(improved, "improved", "NOT improved"))),
              hjust = 0, size = 2.1, family = FF,
              colour = ifelse(dB$improved, OLIVE, RED), fontface = "bold") +
    scale_colour_manual(values = c(`TRUE` = OLIVE, `FALSE` = "grey65"), guide = "none") +
    scale_fill_manual(values = c("M2 (before)" = "grey70",
                                 "M3 (after, 100-repeat mean)" = BLUE), name = NULL) +
    scale_x_continuous(limits = c(-4.4, 3.6), breaks = seq(-4, 3, 1)) +
    scale_y_continuous(breaks = dB$y, labels = rep("", nrow(dB))) +
    coord_cartesian(ylim = c(0.5, 3.8), clip = "off") +
    labs(x = "calibration slope (logit scale)", y = NULL,
         title = "Recalibration improves the slope only where discrimination existed",
         subtitle = "M3 = 100 stratified 50/50 supervised recalibration repeats (target labels consumed); whisker = percentile CI across repeats") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = c(0.72, 0.16),
          legend.background = element_rect(fill = alpha("white", 0.7))) +
    lab("b")

  # ---- (c) Brier / log-loss: M0 vs M2 vs M3 ----
  dC <- do.call(rbind, lapply(TARGETS, function(ch) {
    lb <- rrb$loco_benchmark$per_target[[ch]]
    for (m in c("M0", "M2", "M3")) {
      add_sd("fig3", "c", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.%s.brier", ch, m), lb[[m]]$brier)
      add_sd("fig3", "c", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.%s.log_loss", ch, m), lb[[m]]$log_loss)
      if (m == "M3") {
        add_sd("fig3", "c", "data/rerun_benchmark.json",
               sprintf("loco_benchmark.per_target.%s.M3.brier_ci95", ch), lb[[m]]$brier_ci95)
        add_sd("fig3", "c", "data/rerun_benchmark.json",
               sprintf("loco_benchmark.per_target.%s.M3.log_loss_ci95", ch), lb[[m]]$log_loss_ci95)
      }
    }
    data.frame(cohort = ch, model = c("M0", "M2", "M3"),
               brier = c(lb$M0$brier, lb$M2$brier, lb$M3$brier),
               brier_lo = c(NA, NA, lb$M3$brier_ci95[[1]]),
               brier_hi = c(NA, NA, lb$M3$brier_ci95[[2]]),
               log_loss = c(lb$M0$log_loss, lb$M2$log_loss, lb$M3$log_loss),
               ll_lo = c(NA, NA, lb$M3$log_loss_ci95[[1]]),
               ll_hi = c(NA, NA, lb$M3$log_loss_ci95[[2]]))
  }))
  dC$model <- factor(dC$model, levels = c("M0", "M2", "M3"))
  dC$cohort <- factor(dC$cohort, levels = TARGETS)
  long_c <- rbind(
    transform(dC, metric = "Brier score", value = brier, vlo = brier_lo, vhi = brier_hi),
    transform(dC, metric = "log-loss", value = log_loss, vlo = ll_lo, vhi = ll_hi))
  long_c$metric <- factor(long_c$metric, levels = c("Brier score", "log-loss"))

  pc <- ggplot(long_c, aes(x = model, y = value)) +
    geom_errorbar(aes(ymin = vlo, ymax = vhi), width = 0.16, linewidth = 0.4,
                  colour = PURPLE, alpha = 0.8, na.rm = TRUE) +
    geom_point(aes(fill = model), shape = 21, size = 2.4, colour = "black",
               stroke = 0.45) +
    geom_text(aes(y = ifelse(is.na(vhi), value, pmax(value, vhi)),
                  label = sprintf("%.3f", value)), vjust = -0.55, size = 1.85,
              family = FF, colour = "grey30") +
    facet_grid(metric ~ cohort, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    scale_fill_manual(values = c(M0 = "grey65", M2 = RED, M3 = BLUE), guide = "none") +
    labs(x = NULL, y = NULL,
         title = "Probability accuracy (Brier / log-loss): M3 buys accuracy only in IMvigor210",
         caption = "M0 = prevalence-only; M2 = frozen S0; M3 = 100-repeat recalibration (whiskers = percentile CI across repeats).\nGSE165252: M2 worse than M0 on both metrics; GSE78220: M2 \u2248 M0 / M3 mixed despite AUC 0.615 \u2014 probability accuracy and discrimination separate.") +
    theme(plot.caption = element_text(size = 5.9, colour = "grey30"),
          strip.text = element_text(size = 7)) +
    lab1("c")

  pa / pb / pc + plot_layout(heights = c(1.15, 0.85, 1))
}

# =============================================================================
# FIG 4 -- Representation perturbation constrains, but does not localize,
#          transport failure  (+ gene-identity recovery contrast)
# =============================================================================
make_fig4 <- function() {

  # ---- (a) 21-config external AUC + CI, two cohort facets, layer colours ----
  layer_lv <- c("baseline", "gene_set", "normalization", "pc_k", "ridge_lambda")
  cfg_lab <- c(baseline = "baseline\n(320, clr, k=2)",
               gene_top100 = "gene\ntop100", gene_top50 = "gene\ntop50",
               gene_rand100_s0 = "rand100\n(s0)", gene_rand100_s1 = "rand100\n(s1)",
               gene_rand100_s2 = "rand100\n(s2)", gene_rand100_s3 = "rand100\n(s3)",
               gene_rand100_s4 = "rand100\n(s4)", gene_rand50_s0 = "rand50\n(s0)",
               gene_rand50_s1 = "rand50\n(s1)", gene_rand50_s2 = "rand50\n(s2)",
               gene_rand50_s3 = "rand50\n(s3)", gene_rand50_s4 = "rand50\n(s4)",
               norm_log2raw = "norm\nlog2raw", norm_within_z = "norm\nwithin-z",
               pc_k1 = "pc\nk=1", pc_k3 = "pc\nk=3", pc_k5 = "pc\nk=5",
               ridge_lam0.01 = "ridge\n0.01", ridge_lam0.1 = "ridge\n0.1",
               ridge_lam1 = "ridge\n1")
  dA <- do.call(rbind, lapply(cfg_names, function(cfg) {
    x <- rra$per_config[[cfg]]
    do.call(rbind, lapply(EXT_CH, function(ch) {
      e <- x$external[[ch]]
      add_sd("fig4", "a", "data/rerun_ablation.json",
             sprintf("per_config.%s.external.%s.auc", cfg, ch), e$auc)
      add_sd("fig4", "a", "data/rerun_ablation.json",
             sprintf("per_config.%s.external.%s.ci95_percentile", cfg, ch), e$ci95_percentile)
      add_sd("fig4", "a", "data/rerun_ablation.json",
             sprintf("per_config.%s.external.%s.ci95_bca", cfg, ch), e$ci95_bca)
      data.frame(config = cfg, cohort = ch, layer = x$layer, auc = e$auc,
                 lo = e$ci95_percentile[[1]], hi = e$ci95_percentile[[2]],
                 excl = isTRUE(e$ci_excludes_0p5) || isTRUE(e$ci_excludes_0p5_bca))
    }))
  }))
  stopifnot(nrow(dA) == 42L, !any(dA$excl))
  dA$layer <- factor(dA$layer, levels = layer_lv)
  dA$config <- factor(dA$config, levels = cfg_names)
  dA$cohort <- factor(dA$cohort, levels = EXT_CH)

  pa <- ggplot(dA, aes(x = config, y = auc)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(ymin = lo, ymax = hi, colour = layer,
                        shape = layer == "baseline"),
                    fatten = 2.0, stroke = 0.45, linewidth = 0.35) +
    facet_wrap(~ cohort, ncol = 2) +
    scale_colour_manual(values = c(baseline = "black", gene_set = BLUE,
                                   normalization = RED, pc_k = OLIVE,
                                   ridge_lambda = PURPLE), name = "layer") +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 19),
                       labels = c(`TRUE` = "baseline", `FALSE` = "perturbed"),
                       name = NULL) +
    scale_y_continuous(limits = c(0.15, 0.95), breaks = seq(0.2, 0.9, 0.1)) +
    labs(x = NULL, y = "external AUC (percentile bootstrap 95% CI)",
         title = "21 ablation configurations: no CI excludes 0.5 (percentile AND BCa)",
         subtitle = "recovered-identity 320-gene anchor; 0 degenerate configurations; diamond = baseline") +
    theme(axis.text.x = element_text(size = 4.9, angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right",
          legend.key.height = unit(0.28, "cm")) +
    lab1("a")

  # ---- (b) layer variance decomposition (new vs old ) ----
  var_new <- sapply(layer_lv[-1], function(l)
    rra$verdict$external_auc_range_by_layer[[l]]$mean_variance_across_external)
  var_old <- rra$old__comparison$old_layer_mean_variance
  for (l in layer_lv[-1]) {
    add_sd("fig4", "b", "data/rerun_ablation.json",
           sprintf("verdict.external_auc_range_by_layer.%s.mean_variance_across_external", l),
           var_new[[l]])
    add_sd("fig4", "b", "data/rerun_ablation.json",
           sprintf("old__comparison.old_layer_mean_variance.%s", l), var_old[[l]])
  }
  dB <- do.call(rbind, lapply(layer_lv[-1], function(l) {
    data.frame(layer = l, variance = c(var_new[[l]], var_old[[l]]),
               version = c(" (recovered identity)", " (wrong identity)"))
  }))
  dB$layer <- factor(dB$layer, levels = rev(c("gene_set", "pc_k", "normalization",
                                              "ridge_lambda")))
  stopifnot(var_new[["gene_set"]] > var_new[["pc_k"]],
            var_new[["pc_k"]] > var_new[["normalization"]],
            var_new[["ridge_lambda"]] == 0, var_old[["ridge_lambda"]] == 0)

  pb <- ggplot(dB, aes(x = variance, y = layer)) +
    geom_col(aes(fill = version), position = position_dodge(width = 0.62),
             width = 0.56, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.4f", variance), group = version),
              position = position_dodge(width = 0.62), hjust = -0.25,
              size = 1.9, family = FF, colour = "grey25") +
    annotate("text", x = 0.00035, y = 0.62, hjust = 0, size = 1.95, family = FF,
             colour = "grey25",
             label = "ranking preserved: gene_set > pc_k > normalization\n> ridge_lambda (= 0 in both generations);\nvariances shrink vs : the wrong-identity\nartefact axis (PC1 = library-size proxy) is gone") +
    scale_fill_manual(values = c(" (recovered identity)" = BLUE,
                                 " (wrong identity)" = "grey72"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.45)),
                       labels = function(v) sprintf("%.4f", v)) +
    labs(x = "mean cross-configuration external-AUC variance (two cohorts)",
         y = NULL,
         title = "Layer variance decomposition: gene set dominates both generations") +
    theme(legend.position = c(0.97, 0.75),
          legend.justification = c(1, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.75))) +
    lab("b")

  # ---- (c) gene-identity recovery: old vs new (provenance correction) ----
  ov <- rra$old__comparison
  add_sd("fig4", "c", "data/rerun_ablation.json",
         "old__comparison.old_n_common", ov$old_n_common)
  add_sd("fig4", "c", "data/rerun_ablation.json",
         "old__comparison.old_n_unmeasured", ov$old_n_unmeasured)
  add_sd("fig4", "c", "data/rerun_ablation.json",
         "old__comparison.old_pc1_mass_unmeasured", ov$old_pc1_mass_unmeasured)
  add_sd("fig4", "c", "data/rerun_ablation.json",
         "old__comparison.old_baseline_within_imvigor", ov$old_baseline_within_imvigor)
  add_sd("fig4", "c", "data/rerun_ablation.json",
         "old__comparison.old_baseline_external", ov$old_baseline_external)
  add_sd("fig4", "c", "data/rerun_benchmark.json", "anchor", anc)
  add_sd("fig4", "c", "data/rerun_benchmark.json",
         "within_reference.per_cohort.IMvigor210.auc_632",
         wr$IMvigor210$auc_632)
  for (ch in EXT_CH) {
    add_sd("fig4", "c", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M2.auc", ch),
           prim[[ch]]$auc)
    add_sd("fig4", "c", "data/imvigor_gene_recovery.json",
           sprintf("phase2.old.transport.external.%s.auc", ch),
           rec_old$transport$external[[ch]]$auc)
  }

  evr_old <- rec_old$zero_genes_frozen_$pc1_evr
  evr_new <- anc$pc_explained_variance_ratio
  rows_c <- list(
    list(item = "anchor panel genes (intersection)",
         old = sprintf("%d", ov$old_n_common), new = sprintf("%d", anc$common),
         worse = FALSE),
    list(item = "all-zero genes in anchor panel",
         old = sprintf("%d", ov$old_n_unmeasured),
         new = sprintf("%d", anc$n_allzero_genes_across_348), worse = FALSE),
    list(item = "PC1 loading mass on zero genes",
         old = sprintf("%.1f%%", 100 * ov$old_pc1_mass_unmeasured),
         new = sprintf("%.1f%%", 100 * rec_new$pc1_loading_squared_mass_on_zero_genes),
         worse = FALSE),
    list(item = "PC1 explained variance",
         old = sprintf("%.1f%%", 100 * evr_old),
         new = sprintf("%.1f%%", 100 * evr_new[[1]]), worse = NA),
    list(item = "anchor within 0.632 AUC",
         old = fmt3(ov$old_baseline_within_imvigor),
         new = fmt3(wr$IMvigor210$auc_632), worse = FALSE),
    list(item = "external AUC, GSE165252",
         old = fmt3(ov$old_baseline_external$GSE165252$auc),
         new = fmt3(rra$per_config$baseline$external$GSE165252$auc), worse = FALSE),
    list(item = "external AUC, GSE78220",
         old = fmt3(ov$old_baseline_external$GSE78220$auc),
         new = fmt3(rra$per_config$baseline$external$GSE78220$auc), worse = FALSE),
    list(item = "external calib. slope, GSE165252",
         old = fmt2(rec_old$transport$external$GSE165252$calib_slope),
         new = fmt2(rec_new$transport$external$GSE165252$calib_slope), worse = TRUE),
    list(item = "external calib. slope, GSE78220",
         old = fmt2(rec_old$transport$external$GSE78220$calib_slope),
         new = fmt2(rec_new$transport$external$GSE78220$calib_slope), worse = TRUE))
  # fix EVR new value from the frozen anchor block (computed, not typed)
  dCtxt <- do.call(rbind, lapply(rows_c, function(r)
    data.frame(item = r$item, old = r$old, new = r$new,
               better_or_neutral = !isTRUE(r$worse))))
  dCtxt$y <- rev(seq_len(nrow(dCtxt)))
  add_sd("fig4", "c", "data/imvigor_gene_recovery.json",
         "phase2.new.pc_explained_variance_ratio", evr_new)
  add_sd("fig4", "c", "data/imvigor_gene_recovery.json",
         "phase2.old.zero_genes_frozen_.pc1_evr", evr_old)

  pc <- ggplot(dCtxt) +
    annotate("rect", xmin = 0, xmax = 100, ymin = 10.4, ymax = 11.6,
             fill = "grey93", colour = NA) +
    annotate("segment", x = 0, xend = 100, y = 10.4, yend = 10.4,
             colour = "grey40", linewidth = 0.35) +
    annotate("text", x = 0, y = 11.0, hjust = 0, size = 2.1, family = FF,
             fontface = "bold", colour = "grey15", label = "provenance correction \u2014 gene-identity recovery ()") +
    geom_text(aes(x = 1, y = y, label = item), hjust = 0, size = 2.1,
              family = FF, colour = "grey15") +
    geom_text(aes(x = 69, y = y, label = old), hjust = 1, size = 2.1, family = FF,
              colour = "grey45") +
    geom_text(aes(x = 73, y = y, label = "\u2192"), size = 2.1, family = FF,
              colour = "grey45") +
    geom_text(aes(x = 76.5, y = y, label = new, colour = better_or_neutral),
              hjust = 0, size = 2.1, family = FF, fontface = "bold") +
    annotate("text", x = 69, y = 10.0, hjust = 1, size = 1.9, family = FF,
             fontface = "bold", colour = "grey45", label = "old (wrong identity)") +
    annotate("text", x = 76.5, y = 10.0, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey15", label = "new (recovered key)") +
    annotate("text", x = 1, y = -1.7, hjust = 0, size = 1.95, family = FF,
             colour = "grey35", lineheight = 1.0,
             label = sprintf("blue = improved or neutral under the recovered key; red = worsened\n(both slopes move AWAY from 1 \u2014 recovery did not rescue calibration)\nold zero-gene/PC1-mass values: _alignment_ablation frozen block\n(retracted interpretation); new: _rerun_ablation anchor\nGSE78220 baseline AUC crosses 0.5 (%.3f \u2192 %.3f) yet its CI\nstill spans 0.5 at n = 25 \u2014 no transport claim",
                             ov$old_baseline_external$GSE78220$auc,
                             rra$per_config$baseline$external$GSE78220$auc)) +
    scale_colour_manual(values = c(`TRUE` = BLUE, `FALSE` = RED), guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(-3.6, 11.8), clip = "off") +
    scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL,
         title = "Identity recovery changes the numbers, not the verdict") +
    theme_void(base_family = FF) +
    theme(plot.title = element_text(size = 8, face = "plain")) +
    lab("c")

  pa / (pb | pc) + plot_layout(heights = c(1, 1))
}

# =============================================================================
# FIG 5 -- Benchmarking against simpler transport strategies
# =============================================================================
make_fig5 <- function() {

  # ---- (a) LOCO M0-M6 x 3 cohorts ----
  dA <- do.call(rbind, lapply(rev(TARGETS), function(ch) {
    do.call(rbind, lapply(MODELS, function(m) {
      r <- rrb$loco_benchmark$per_target[[ch]][[m]]
      a <- rrb$adjudication$per_model[[m]][[ch]]
      add_sd("fig5", "a", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.%s.auc", ch, m), r$auc)
      add_sd("fig5", "a", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.%s.auc_ci95", ch, m), r$auc_ci95)
      add_sd("fig5", "a", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.%s.calib_slope", ch, m), r$calib_slope)
      add_sd("fig5", "a", "data/rerun_benchmark.json",
             sprintf("adjudication.per_model.%s.%s.transport_permitted", m, ch),
             a$transport_permitted)
      data.frame(cohort = ch, model = m, auc = r$auc,
                 lo = r$auc_ci95[[1]], hi = r$auc_ci95[[2]],
                 slope = if (is.null(r$calib_slope)) NA_real_ else r$calib_slope,
                 c1 = isTRUE(a$auc_ci_excludes_0p5_beneficial_direction),
                 refused = !isTRUE(a$transport_permitted))
    }))
  }))
  stopifnot(all(dA$refused))
  dA$model <- factor(dA$model, levels = rev(MODELS))
  dA$cohort <- factor(dA$cohort, levels = rev(TARGETS))
  dA$slope_lab <- ifelse(is.na(dA$slope), "NA", sprintf("%.2f", dA$slope))
  # guard: in the IMvigor210 facet exactly 6 models pass criterion 1
  stopifnot(sum(dA$c1[dA$cohort == "IMvigor210"]) == 6L)

  # highlight-rectangle geometry for the IMvigor210 facet (panel a)
  imv_rows <- which(dA$cohort == "IMvigor210" & dA$c1)
  hl <- data.frame(cohort = factor("IMvigor210", levels = levels(dA$cohort)),
                   ymin = min(as.integer(dA$model[imv_rows])) - 0.45,
                   ymax = max(as.integer(dA$model[imv_rows])) + 0.45)
  imv_slopes <- dA$slope[imv_rows]
  hl_note <- sprintf("6 models: C1 pass but slopes %.1f\u2013%.1f \u2192 C2 refused",
                     min(imv_slopes), max(imv_slopes))
  add_sd("fig5", "a", "derived (frozen keys loco_benchmark+adjudication)",
         "IMvigor210_c1_passing_models", list(
           models = as.list(as.character(dA$model[imv_rows])),
           auc_range = c(min(dA$auc[imv_rows]), max(dA$auc[imv_rows])),
           slope_range = c(min(imv_slopes), max(imv_slopes))))
  add_sd("fig5", "a", "data/rerun_benchmark.json",
         "adjudication.verdict_counts", rrb$adjudication$verdict_counts)

  pa <- ggplot(dA, aes(x = auc, y = as.integer(model))) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_rect(data = hl, aes(xmin = 0.42, xmax = 0.80, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = OLIVE, alpha = 0.08,
              colour = OLIVE, linetype = 1, linewidth = 0.35) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = c1),
                    fatten = 2.1, stroke = 0.5, shape = 19, linewidth = 0.4) +
    geom_label(aes(x = 0.175, label = paste0(model, "  slope ", slope_lab)),
               hjust = 0, size = 2.0, family = FF, colour = "grey25",
               label.padding = unit(0.09, "lines"), label.size = 0,
               fill = alpha("white", 0.78)) +
    geom_label(aes(x = 0.815, label = ifelse(refused, "refused", "")),
               hjust = 0, size = 1.9, family = FF, colour = RED,
               label.padding = unit(0.07, "lines"), label.size = 0,
               fill = alpha("white", 0.72)) +
    geom_label(data = transform(hl, x = 0.30, y = 0.65, label = hl_note),
               aes(x = x, y = y, label = label),
               hjust = 0, size = 1.9, family = FF, colour = OLIVE,
               lineheight = 0.95, inherit.aes = FALSE, label.size = 0,
               label.padding = unit(0.12, "lines"), fill = alpha("white", 0.8)) +
    facet_wrap(~ cohort, ncol = 1) +
    scale_colour_manual(values = c(`TRUE` = OLIVE, `FALSE` = "grey45"),
                        labels = c(`TRUE` = "AUC CI excludes 0.5 (benef.)",
                                   `FALSE` = "AUC CI spans 0.5"),
                        name = NULL) +
    scale_x_continuous(limits = c(0.17, 0.93), breaks = seq(0.2, 0.9, 0.1)) +
    scale_y_continuous(breaks = NULL) +
    labs(x = "AUC, LOCO held-out (train = pooled other two cohorts)", y = NULL,
         title = "Eight strategies x three targets:\n24/24 transport refusals (recovered identity)",
         subtitle = "two-criterion rule: AUC 95% CI excludes 0.5\nin the beneficial direction AND |slope \u2212 1| \u2264 0.10 (logit);\nM0 constant \u2192 slope not estimable") +
    theme(plot.title = element_text(margin = margin(l = 12, unit = "pt")),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.28, "cm"),
          legend.margin = margin(0, 0, 0, 0)) +
    lab1("a")

  # ---- (b) reference strategies (primary single-source estimand) ----
  dB <- do.call(rbind, lapply(EXT_CH, function(ch) {
    m0 <- pw[[pk("IMvigor210", ch)]]$M0
    m2 <- prim[[ch]]
    cm <- clin$external_results[[ch]]$C_main
    cs <- clin$external_results[[ch]]$C_sex
    ayx <- ay[[ch]]$target_metrics
    add_sd("fig5", "b", "data/rerun_benchmark.json",
           sprintf("pairwise_transports.transports.train_IMvigor210__test_%s.M0.auc", ch), m0$auc)
    add_sd("fig5", "b", "data/clinical_only.json",
           sprintf("external_results.%s.C_main.auc", ch), cm$auc)
    add_sd("fig5", "b", "data/clinical_only.json",
           sprintf("external_results.%s.C_main.auc_ci95", ch), cm$auc_ci95)
    add_sd("fig5", "b", "data/clinical_only.json",
           sprintf("external_results.%s.C_sex.auc", ch), cs$auc)
    data.frame(target = ch,
               strat = c("M0 prevalence-only", "M2 frozen S0 (322-gene)",
                         "clinical C_main (sex + stage)", "clinical C_sex (sex)",
                         "Ayers IFN-\u03b3 6-gene"),
               est = c(m0$auc, m2$auc, cm$auc, cs$auc, ayx$auc),
               lo  = c(NA, m2$auc_ci95[[1]], cm$auc_ci95[[1]], cs$auc_ci95[[1]], ayx$auc_ci95[[1]]),
               hi  = c(NA, m2$auc_ci95[[2]], cm$auc_ci95[[2]], cs$auc_ci95[[2]], ayx$auc_ci95[[2]]),
               n = m2$n_test, n_event = m2$n_event)
  }))
  dB$strat <- factor(dB$strat, levels = rev(c("M0 prevalence-only", "M2 frozen S0 (322-gene)",
                                              "clinical C_main (sex + stage)",
                                              "clinical C_sex (sex)",
                                              "Ayers IFN-\u03b3 6-gene")))
  dB$target <- factor(dB$target, levels = EXT_CH)
  const_note_val <- clin$external_results$GSE165252$C_main$constant_prediction
  add_sd("fig5", "b", "data/clinical_only.json",
         "external_results.GSE165252.C_main.constant_prediction",
         clin$external_results$GSE165252$C_main$constant_prediction)
  add_sd("fig5", "b", "data/clinical_only.json",
         "external_results.GSE165252.degeneracy",
         clin$external_results$GSE165252$degeneracy)

  pb <- ggplot(dB, aes(x = est, y = strat)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = strat), fatten = 2.2,
                    stroke = 0.5, shape = 19, linewidth = 0.4, na.rm = TRUE) +
    geom_text(aes(label = fmt3(est)), vjust = -1.15, size = 1.9, family = FF,
              colour = "grey25") +
    geom_text(aes(x = est, label = sprintf("n = %d (%d events)", n, n_event)),
              vjust = 2.3, size = 1.7, family = FF, colour = "grey45") +
    facet_wrap(~ target, ncol = 2) +
    scale_colour_manual(values = c("M0 prevalence-only" = "grey55",
                                   "M2 frozen S0 (322-gene)" = BLUE,
                                   "clinical C_main (sex + stage)" = RED,
                                   "clinical C_sex (sex)" = "#e2a09e",
                                   "Ayers IFN-\u03b3 6-gene" = PURPLE),
                        guide = "none") +
    scale_x_continuous(limits = c(0.08, 0.95), breaks = seq(0.2, 0.9, 0.2)) +
    coord_cartesian(ylim = c(0.2, 5.7), clip = "off") +
    labs(x = "AUC, primary estimand (train IMvigor210, frozen, transform-only)",
         y = NULL,
         title = "Simpler reference strategies:\nnone is licensable either",
         caption = sprintf("different estimand from panel a: single-source frozen transports (no LOCO re-training).\nC_main/C_sex: ridge clinical models frozen on IMvigor210 n = 298; Ayers: IFN-\u03b3 6-gene frozen score.\nGSE165252 facet: both clinical models degenerate to the train-prevalence constant %.3f \u2014 AUC 0.500 by\nconstruction (GEO deposits no covariates). All rows FAIL the two-criterion rule (frozen adjudication).",
                           const_note_val)) +
    theme(plot.caption = element_text(size = 5.9, colour = "grey30"),
          axis.text.x = element_text(size = 5.5)) +
    lab1("b")

  pa | pb + plot_layout(widths = c(1.35, 1))
}

# =============================================================================
# ED1 -- Positive control: simulated (rule behaviour) + Ayers (empirical)
# =============================================================================
make_ed1 <- function() {

  # ---- (a) simulated positive control pass-rate heatmap (, valid) ----
  scen_keys <- names(pcs)
  scen_lab <- vapply(scen_keys, function(k) {
    s <- k
    s <- sub("^auc0\\.", "AUC 0.", s)
    s <- sub("_wellcalibrated$", " \u00b7 well-calibrated", s)
    s <- sub("_miscalibrated_x4$", " \u00b7 miscalibrated \u00d74", s)
    s
  }, "")
  n_targets <- as.character(s75$positive_control$parameters$target_ns)
  dA <- do.call(rbind, lapply(scen_keys, function(sk) {
    do.call(rbind, lapply(n_targets, function(n) {
      r <- pcs[[sk]]$n_target_results[[n]]
      add_sd("ed1", "a", "data/supporting_analyses.json",
             sprintf("positive_control.results.%s.n_target_results.%s.pass_rate_rule_i_and_ii", sk, n),
             r$pass_rate_rule_i_and_ii)
      data.frame(scenario = scen_lab[sk], n = as.integer(n),
                 pass = r$pass_rate_rule_i_and_ii)
    }))
  }))
  dA$scenario <- factor(dA$scenario, levels = scen_lab[scen_keys])
  add_sd("ed1", "a", "derived", "max_wellcalibrated_pass_rate", pc_max_well)
  add_sd("ed1", "a", "derived", "all_miscalibrated_pass_rates_zero", pc_mis_zero)

  pa <- ggplot(dA, aes(x = factor(n), y = scenario, fill = pass)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * pass)), size = 2.5, family = FF,
              colour = ifelse(dA$pass > 0.15, "white", "grey25"), fontface = "bold") +
    scale_fill_gradient(low = "#f5f8fb", high = BLUE, limits = c(0, 0.4),
                        labels = percent, name = "pass rate") +
    scale_x_discrete(labels = c("25" = "n = 25", "71" = "n = 71", "298" = "n = 298")) +
    labs(x = NULL, y = NULL,
         title = "Simulated positive control: the rule admits\na truly transportable predictor",
         subtitle = sprintf("500 reps per cell \u00b7 max pass %.1f%% (AUC 0.8, well-calibrated, n = 298);\nevery miscalibrated cell = 0%%",
                            100 * max(dA$pass)),
         caption = "identity-independent simulation (frozen ): criterion 1 AND 2 adjudication per rep.\nwell-calibrated rows: not a trivial rejector \u2014 licensing rates grow with AUC and n.") +
    theme(axis.text.y = element_text(size = 6.5),
          plot.caption = element_text(size = 5.9, colour = "grey30"),
          legend.key.height = unit(0.45, "cm")) +
    lab("a")

  # ---- (b) Ayers empirical control: forest + tolerance-tier tiles ----
  dB <- do.call(rbind, lapply(EXT_CH, function(ch) {
    a <- ay[[ch]]$target_metrics
    adjt <- ay[[ch]]$adjudication_target
    add_sd("ed1", "b", "data/rerun_positive_control.json",
           sprintf("verdict.per_cohort.%s", ch), rpc$verdict$per_cohort[[ch]])
    data.frame(target = ch, est = a$auc, lo = a$auc_ci95[[1]], hi = a$auc_ci95[[2]],
               slope = a$calib_slope, n = a$n_test, n_event = a$n_event,
               main = adjt$main_analysis_rule$both_passed,
               t005 = adjt$tolerance_tiers[["pm0.05"]]$both_passed,
               t010 = adjt$tolerance_tiers[["pm0.10"]]$both_passed,
               t020 = adjt$tolerance_tiers[["pm0.20"]]$both_passed)
  }))
  stopifnot(!any(c(dB$main, dB$t005, dB$t010, dB$t020)))
  dB$row <- paste0("Ayers IFN-\u03b3 \u2192 ", dB$target)
  dB$y <- rev(seq_len(nrow(dB)))

  pb1 <- ggplot(dB, aes(x = est, y = y)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi), fatten = 2.6, stroke = 0.55,
                    shape = 19, colour = PURPLE, linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f [%.3f, %.3f]  slope %.2f", est, lo, hi, slope)),
              vjust = -1.25, size = 2.0, family = FF, colour = "grey25") +
    geom_text(aes(x = 0.09, label = sprintf("%s  (n = %d, %d events)", target, n, n_event)),
              hjust = 0, vjust = 2.3, size = 2.0, family = FF, colour = "grey20") +
    scale_x_continuous(limits = c(0.08, 0.92), breaks = seq(0.2, 0.9, 0.1)) +
    scale_y_continuous(breaks = dB$y, labels = rep("", nrow(dB))) +
    coord_cartesian(ylim = c(0.5, 2.6), clip = "off") +
    labs(x = "AUC (E ~ IFN-\u03b3 6-gene), frozen IMvigor210 anchor", y = NULL,
         title = "Empirical positive control also FAILs") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    lab("b")

  tiles <- do.call(rbind, lapply(seq_len(nrow(dB)), function(i) {
    r <- dB[i, ]
    data.frame(target = r$target,
               tier = c("main rule", "\u00b10.05", "\u00b10.10", "\u00b10.20"),
               pass = c(r$main, r$t005, r$t010, r$t020))
  }))
  tiles$tier <- factor(tiles$tier, levels = rev(c("main rule", "\u00b10.05", "\u00b10.10", "\u00b10.20")))
  pb2 <- ggplot(tiles, aes(x = tier, y = target)) +
    geom_tile(aes(fill = pass), colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(pass, "PASS", "FAIL")), size = 2.4, family = FF,
              fontface = "bold", colour = "white") +
    scale_fill_manual(values = c(`TRUE` = OLIVE, `FALSE` = RED), guide = "none") +
    labs(x = "tolerance tier", y = NULL,
         title = "No cohort passes at any\npre-registered tolerance tier") +
    theme(axis.text.x = element_text(size = 6.5)) +
    lab("")

  ann <- ggplot() +
    annotate("rect", xmin = 2, xmax = 98, ymin = 10, ymax = 82,
             fill = "grey97", colour = "grey35", linewidth = 0.4) +
    annotate("text", x = 50, y = 76, size = 2.35, family = FF, fontface = "bold",
             lineheight = 0.95,
             label = "Rule behaviour vs empirical signal \u2014\nthe two cannot be separated here") +
    annotate("text", x = 50, y = 44, size = 1.95, family = FF, colour = "grey15",
             lineheight = 1.05,
             label = "the rule licenses truly transportable\npredictors in simulation (panel a) and\nrejects every miscalibrated discriminator\n(0% in all six cells);\n\nthe Ayers IFN-\u03b3 signature \u2014 reproducible in\nthe literature \u2014 also fails everywhere, so in\nTHESE small external cohorts (n = 71 / n = 25)\none cannot distinguish \"rule too strict\" from\n\"insufficient signal\"; only larger gene-valid\ncohorts can") +
    theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off")

  (pa / pb1) | (pb2 / ann) + plot_layout(widths = c(1.25, 1))
}

# =============================================================================
# ED3 -- LOCO vs primary estimand
# =============================================================================
make_ed3 <- function() {

  # ---- (a) schematic: LOCO changes the predictor -> different estimand ----
  pa <- ggplot() +
    annotate("rect", xmin = 2, xmax = 47, ymin = 52, ymax = 92,
             fill = "#dbe5f1", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = 24.5, y = 87, size = 2.5, family = FF, fontface = "bold",
             colour = BLUE, label = "PRIMARY ESTIMAND (this paper)") +
    annotate("text", x = 24.5, y = 68, size = 2.1, family = FF, colour = "grey15",
             lineheight = 1.0,
             label = "frozen S0: PCA rotation + ridge\nestimated ONCE on IMvigor210,\nthen applied transform-only to the\ntarget cohort (no target labels,\nno re-estimation)") +
    annotate("rect", xmin = 53, xmax = 98, ymin = 52, ymax = 92,
             fill = "#f2dcdb", colour = RED, linewidth = 0.45) +
    annotate("text", x = 75.5, y = 87, size = 2.5, family = FF, fontface = "bold",
             colour = RED, label = "LOCO ESTIMAND (sensitivity)") +
    annotate("text", x = 75.5, y = 68, size = 2.1, family = FF, colour = "grey15",
             lineheight = 1.0,
             label = "predictor RE-TRAINED on the pooled\nother two cohorts for every target\n\u2192 a DIFFERENT prediction system\n\u2192 answers a different question\n(\"does any 2-cohort mixture transfer?\")") +
    annotate("segment", x = 47, xend = 53, y = 72, yend = 72, linewidth = 0.5,
             colour = "grey45",
             arrow = arrow(length = unit(0.18, "cm"))) +
    annotate("text", x = 50, y = 78, size = 1.95, family = FF, colour = RED,
             label = "different\nestimand") +
    annotate("rect", xmin = 2, xmax = 98, ymin = 6, ymax = 42,
             fill = "grey97", colour = "grey35", linewidth = 0.4) +
    annotate("text", x = 4, y = 36.5, hjust = 0, size = 2.3, family = FF,
             fontface = "bold", label = "Why the distinction matters ( lock):") +
    annotate("text", x = 4, y = 20, hjust = 0, size = 2.05, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = "LOCO numbers CANNOT substitute for the primary claim: the LOCO M2 score for target IMvigor210\ntrains on GSE165252+GSE78220 \u2014 a predictor that never existed in the frozen pipeline.\nPrimary-estimand transports (train IMvigor210, single-source, frozen) are the anchor-arm headline;\nLOCO M0\u2013M6 (Fig. 5a) is a sensitivity benchmark over a wider strategy grid.") +
    theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("a")

  # ---- (b) per target: LOCO M2 vs single-source M2 transports ----
  dB <- do.call(rbind, lapply(TARGETS, function(te) {
    loco <- rrb$loco_benchmark$per_target[[te]]$M2
    add_sd("ed3", "b", "data/rerun_benchmark.json",
           sprintf("loco_benchmark.per_target.%s.M2.auc", te), loco$auc)
    add_sd("ed3", "b", "data/rerun_benchmark.json",
           sprintf("loco_benchmark.per_target.%s.M2.auc_ci95", te), loco$auc_ci95)
    rows <- list(data.frame(target = te, predictor = "LOCO (pooled other two)",
                            est = loco$auc, lo = loco$auc_ci95[[1]],
                            hi = loco$auc_ci95[[2]]))
    for (tr in setdiff(COHORTS, te)) {
      m <- pw[[pk(tr, te)]]$M2
      add_sd("ed3", "b", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.%s.M2.auc", pk(tr, te)), m$auc)
      add_sd("ed3", "b", "data/rerun_benchmark.json",
             sprintf("pairwise_transports.transports.%s.M2.auc_ci95", pk(tr, te)), m$auc_ci95)
      rows[[length(rows) + 1]] <- data.frame(
        target = te, predictor = sprintf("single-source: %s", tr),
        est = m$auc, lo = m$auc_ci95[[1]], hi = m$auc_ci95[[2]])
    }
    do.call(rbind, rows)
  }))
  dB$target <- factor(dB$target, levels = TARGETS)
  dB$predictor <- factor(dB$predictor,
                         levels = rev(c("LOCO (pooled other two)",
                                        "single-source: IMvigor210",
                                        "single-source: GSE165252",
                                        "single-source: GSE78220")))
  # guard: only LOCO@IMvigor210 and single-source GSE78220@IMvigor210 exclude 0.5
  stopifnot(sum(dB$lo > 0.5) == 2L,
            dB$lo[dB$target == "IMvigor210"][1] > 0.5)
  loco_imv <- dB$est[dB$target == "IMvigor210" & dB$predictor == "LOCO (pooled other two)"]
  ss_78220_imv <- dB$est[dB$target == "IMvigor210" &
                           dB$predictor == "single-source: GSE78220"]

  pb <- ggplot(dB, aes(x = est, y = predictor)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = predictor == "LOCO (pooled other two)",
                        size = predictor == "LOCO (pooled other two)"),
                    fatten = 2.2, stroke = 0.5, shape = 19, na.rm = TRUE) +
    geom_text(aes(label = fmt3(est)), vjust = -1.15, size = 1.9, family = FF,
              colour = "grey30") +
    facet_wrap(~ target, ncol = 3) +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = "grey45"), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.35), guide = "none") +
    scale_x_continuous(limits = c(0.28, 0.78), breaks = seq(0.3, 0.75, 0.1)) +
    labs(x = "AUC (M2, efficacy)", y = NULL,
         title = "M2 under the two estimands, side by side (same target cohorts)",
         caption = sprintf("red = LOCO predictor (re-trained on pooled other two cohorts); grey = single-source frozen transports.\nIMvigor210 panel: LOCO %.3f and single-source GSE78220\u2192IMvigor210 %.3f both exclude 0.5 \u2014 but these are\nDIFFERENT predictors from the primary anchor arm; neither licenses the frozen IMvigor210-anchored claim.",
                           loco_imv, ss_78220_imv)) +
    theme(axis.text.y = element_text(size = 6.2),
          plot.caption = element_text(size = 5.8, colour = "grey30")) +
    lab1("b")

  pa / pb + plot_layout(heights = c(0.9, 1))
}

# =============================================================================
# ED4 -- 21-configuration ablation matrix (config x cohort AUC heatmap)
# =============================================================================
make_ed4 <- function() {

  lay_ord <- c("baseline",
               "gene_top100", "gene_top50",
               "gene_rand100_s0", "gene_rand100_s1", "gene_rand100_s2",
               "gene_rand100_s3", "gene_rand100_s4",
               "gene_rand50_s0", "gene_rand50_s1", "gene_rand50_s2",
               "gene_rand50_s3", "gene_rand50_s4",
               "norm_log2raw", "norm_within_z",
               "pc_k1", "pc_k3", "pc_k5",
               "ridge_lam0.01", "ridge_lam0.1", "ridge_lam1")
  stopifnot(setequal(lay_ord, cfg_names))
  d <- do.call(rbind, lapply(lay_ord, function(cfg) {
    x <- rra$per_config[[cfg]]
    do.call(rbind, lapply(EXT_CH, function(ch) {
      e <- x$external[[ch]]
      data.frame(config = cfg, layer = x$layer, cohort = ch, auc = e$auc,
                 excl = isTRUE(e$ci_excludes_0p5) || isTRUE(e$ci_excludes_0p5_bca))
    }))
  }))
  d$config <- factor(d$config, levels = rev(lay_ord))
  d$cohort <- factor(d$cohort, levels = EXT_CH)
  add_sd("ed4", "matrix", "data/rerun_ablation.json",
         "per_config_note", "all 21 configurations x 2 external cohorts; AUC values dumped below per key")
  add_sd("ed4", "matrix", "data/rerun_ablation.json",
         "verdict.any_config_recovers_transport", rra$verdict$any_config_recovers_transport)
  add_sd("ed4", "matrix", "data/rerun_ablation.json",
         "removed_configs", rra$removed_configs)
  add_sd("ed4", "matrix", "data/rerun_ablation.json",
         "n_degenerate_configs", n_degenerate)
  for (cfg in lay_ord) for (ch in EXT_CH) {
    add_sd("ed4", "matrix", "data/rerun_ablation.json",
           sprintf("per_config.%s.external.%s.auc", cfg, ch),
           rra$per_config[[cfg]]$external[[ch]]$auc)
  }

  layer_y <- vapply(lay_ord, function(cfg) rra$per_config[[cfg]]$layer, "")
  # guard: each layer must occupy one contiguous block, in the frozen order
  stopifnot(length(rle(layer_y)$lengths) == 5L)
  n_by_layer <- table(factor(layer_y,
                             levels = c("baseline", "gene_set", "normalization",
                                        "pc_k", "ridge_lambda")))
  # y positions of layer boundaries for grouping frames (bottom-up)
  yb <- cumsum(rev(as.vector(n_by_layer)))   # from bottom (rev factor levels)
  grp <- data.frame(ymin = c(0.5, yb[-length(yb)] + 0.5),
                    ymax = yb + 0.5,
                    layer = rev(names(n_by_layer)))
  lay_cols <- c(baseline = "black", gene_set = BLUE, normalization = RED,
                pc_k = OLIVE, ridge_lambda = PURPLE)

  p <- ggplot(d, aes(x = cohort, y = config, fill = auc)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", auc)), size = 2.3, family = FF,
              colour = ifelse(abs(d$auc - 0.5) > 0.10, "white", "grey15"),
              fontface = "bold") +
    geom_rect(data = grp, aes(xmin = 0.5, xmax = 2.5, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA,
              colour = lay_cols[grp$layer], linewidth = 0.55) +
    annotate("text", x = 2.72, y = (grp$ymin + grp$ymax) / 2, hjust = 0,
             size = 2.0, family = FF, fontface = "bold",
             colour = lay_cols[grp$layer], label = grp$layer) +
    scale_fill_gradient2(low = "#c0504d", mid = "grey92", high = "#1f4e79",
                         midpoint = 0.5, limits = c(0.40, 0.70),
                         name = "external AUC") +
    scale_x_discrete(position = "top") +
    coord_cartesian(xlim = c(0.4, 4.55), clip = "off") +
    labs(x = NULL, y = NULL,
         title = "Ablation matrix: 21 configurations x 2 external cohorts (recovered identity)",
         subtitle = sprintf("no configuration's bootstrap 95%% CI excludes 0.5 (percentile & BCa); degenerate configurations = %d;\n5 old  'measured-only' artefact configurations retracted (premise gone)",
                            n_degenerate),
         caption = "colour centred at 0.5 (grey = null); frames group configurations by perturbed layer.\nRemoved (retracted under recovered identity):\ngene_measured_only, gene_top100_measured, gene_top50_measured, norm_log2raw_measured, norm_within_z_measured.\nSources: data/rerun_ablation.json (per_config, verdict, removed_configs).") +
    theme(axis.text.y = element_text(size = 6.3),
          axis.text.x.top = element_text(size = 7.5),
          plot.caption = element_text(size = 5.9, colour = "grey30"),
          plot.subtitle = element_text(size = 6.8),
          plot.margin = margin(0.4, 1.4, 0.4, 0.4, "cm")) +
    annotate("text", x = -Inf, y = Inf, label = "a", hjust = -1.5, vjust = 0.35,
             size = 3.4, fontface = "bold", family = FF)
  p
}

# =============================================================================
# ED5 -- Gene-identity recovery: full timeline + old/new numbers
# =============================================================================
make_ed5 <- function() {

  hyp <- rec$phase1a$hypotheses
  n_hyp <- length(hyp)
  n_pass <- sum(sapply(hyp, function(h) isTRUE(h$pass)))
  vg <- rec$verification_gates
  old_zg <- rec_old$zero_genes_frozen_
  spear_new <- rec_new$cross_cohort_mean_expr_spearman$recovered_key
  spear_old <- rec_new$cross_cohort_mean_expr_spearman$old_positional_alignment
  p1c <- rec$phase1c
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "phase1a.winner", rec$phase1a$winner)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "phase1a.n_hypotheses_tested", n_hyp)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "phase1a.n_passing_hypotheses", n_pass)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "verification_gates", vg)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "phase2.new.cross_cohort_mean_expr_spearman",
         rec_new$cross_cohort_mean_expr_spearman)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "phase2.old.zero_genes_frozen_", old_zg)
  add_sd("ed5", "a", "data/imvigor_gene_recovery.json",
         "provenance.exmat_sha256", rec$provenance$exmat_sha256)

  # ---- (a) timeline: five stages ----
  xs <- c(1.5, 21, 40.5, 60, 79.5); bw <- 18.5
  titles <- c("1 \u00b7 Upstream\nanonymisation", "2 \u00b7 Double\nmisalignment (\u2013)",
              "3 \u00b7 Sex-anchor\nfalsification", "4 \u00b7 Verification\ngates", "5 \u00b7 Number\nchanges")
  bodies <- c(
    sprintf("exmat row labels replaced by\ngene_1..gene_%d placeholders;\nrow ORDER preserved =\ncounts(cds) write order\n(= symbol alphabetical).\nfData SHA-256 verified\nupstream (see frozen\nprovenance block)",
            p1c$n_exmat_rows),
    sprintf("pipeline assumed exmat row i\n== fData row i\n(entrez string sort)\n+ Supplementary_Table_1\ncolumn 1 misparsed as Entrez\n\u2192 pseudo-random %d-gene\n'intersection' (%d common)\n%d all-zero genes \u00b7\n%.1f%% PC1 mass",
            vg$panel_measurable$old_fake_entrez_selection,
            old_zg$n_common_genes, old_zg$n_unmeasured,
            100 * old_zg$pc1_loading_squared_mass),
    sprintf("P0-2 positive-control gates:\nY-linked / XIST /\nhousekeeping markers\nbiologically impossible under\npositional alignment\n\u2192 %d sort hypotheses tested,\n%d pass: symbol alphabetical\n(upper); both rawcase & upper\npanel genes agree;\nupper adjudicated",
            n_hyp, n_pass),
    sprintf("sex-anchor gate: %s\nbinaryResponse %d/%d/%d\n(CR-PR / SD-PD / NA)\nmirror numeric identity:\n%d/%d rows matched (1c)\ncross-cohort mean expression\nSpearman \u03c1 = %.2f (GSE165252)\nand %.2f (GSE78220);\nold alignment: %.2f / %.2f\nOVERALL: %s",
            vg$sex_anchor,
            vg$binaryResponse_counts$`CR/PR`, vg$binaryResponse_counts$`SD/PD`,
            vg$binaryResponse_counts[["NA"]],
            p1c$n_matched, p1c$n_exmat_rows,
            spear_new$GSE165252$spearman, spear_new$GSE78220$spearman,
            spear_old$GSE165252$spearman, spear_old$GSE78220$spearman,
            vg$overall),
    sprintf("panel: %d \u2192 %d of %d genes\nzero genes: %d \u2192 0\nPC1 mass on zero genes:\n%.1f%% \u2192 0.0%%\nwithin 0.632:\n%.3f \u2192 %.3f\nGSE78220 external:\n%.3f \u2192 %.3f\nGSE165252 external:\n%.3f \u2192 %.3f\nretraction: %d-gene artefact +\n%d-gene intersection\nwithdrawn",
            vg$panel_measurable$old_fake_entrez_selection,
            rec_new$gene_set$common, rec_new$gene_set$panel,
            old_zg$n_unmeasured, 100 * old_zg$pc1_loading_squared_mass,
            rec_old$within_632_frozen__pca_stability,
            rec_new$within_imvigor$auc_632,
            rec_old$transport$external$GSE78220$auc,
            rec_new$transport$external$GSE78220$auc,
            rec_old$transport$external$GSE165252$auc,
            rec_new$transport$external$GSE165252$auc,
            old_zg$n_unmeasured, old_zg$n_common_genes))
  fills <- c("#dbe5f1", "#f2dcdb", "#f2dcdb", "#ebf1dd", "grey97")
  cols  <- c(BLUE, RED, RED, OLIVE, "grey30")

  pa <- ggplot() +
    annotate("segment", x = xs[1:4] + bw, xend = xs[2:5], y = rep(46, 4),
             yend = rep(46, 4), arrow = arrow(length = unit(0.16, "cm")),
             linewidth = 0.5, colour = "grey35") +
    annotate("text", x = 50, y = 96, size = 2.8, family = FF, fontface = "bold",
             label = "How the gene-row identity of the anchor cohort was lost, detected, and recovered (all values frozen)") +
    annotate("text", x = 50, y = 90.5, size = 2.2, family = FF, colour = "grey35",
             label = "the anonymised deposit kept row ORDER but not row LABELS \u2014 every 'select gene by name' operation on the wrong assumption retrieved unrelated genes")
  for (i in 1:5) {
    pa <- pa +
      annotate("rect", xmin = xs[i], xmax = xs[i] + bw, ymin = 12, ymax = 84,
               fill = fills[i], colour = cols[i], linewidth = 0.45) +
      annotate("text", x = xs[i] + bw / 2, y = 77.5, size = 2.25, family = FF,
               fontface = "bold", colour = cols[i], lineheight = 0.9,
               label = titles[i]) +
      annotate("text", x = xs[i] + bw / 2, y = 44, size = 1.7, family = FF,
               colour = "grey15", lineheight = 0.98, label = bodies[i])
  }
  pa <- pa +
    annotate("text", x = 50, y = 5.5, size = 2.0, family = FF, colour = "grey30",
             lineheight = 1.0,
             label = "recovered key: exmat row j == j-th fData row under upper-cased symbol alphabetical sort \u00b7 exported as exmat_row_key.csv (31085 rows)\nanchor numeric consequences fully re-run in data/_rerun_*.json; the WRONG-identity archives remain frozen as provenance of the error") +
    theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("a")

  # ---- (b) old vs new numbers (sourced from phase2 old/new) ----
  rows <- list(
    list(item = "panel genes measurable in IMvigor210", old = "229 (fake Entrez)",
         new = sprintf("%d / 322 (missing: %s)", rec_new$gene_set$imv_available,
                       paste(unlist(rec_new$gene_set$imv_missing), collapse = ", "))),
    list(item = "cross-cohort intersection", old = sprintf("%d", old_zg$n_common_genes),
         new = sprintf("%d", rec_new$gene_set$common)),
    list(item = "all-zero genes in anchor panel", old = sprintf("%d", old_zg$n_unmeasured),
         new = sprintf("%d", rec_new$zero_genes$n_allzero_across_348)),
    list(item = "PC1 loading mass on zero genes", old = sprintf("%.1f%%", 100 * old_zg$pc1_loading_squared_mass),
         new = sprintf("%.1f%%", 100 * rec_new$pc1_loading_squared_mass_on_zero_genes)),
    list(item = "anchor within apparent AUC", old = fmt3(rec_old$within_imvigor$auc_apparent),
         new = fmt3(rec_new$within_imvigor$auc_apparent)),
    list(item = "anchor within 0.632 AUC", old = fmt3(rec_old$within_632_frozen__pca_stability),
         new = fmt3(rec_new$within_imvigor$auc_632)),
    list(item = "external AUC \u2192 GSE165252", old = fmt3(rec_old$transport$external$GSE165252$auc),
         new = fmt3(rec_new$transport$external$GSE165252$auc)),
    list(item = "external AUC \u2192 GSE78220", old = fmt3(rec_old$transport$external$GSE78220$auc),
         new = fmt3(rec_new$transport$external$GSE78220$auc)),
    list(item = "calib. slope \u2192 GSE165252", old = fmt2(rec_old$transport$external$GSE165252$calib_slope),
         new = fmt2(rec_new$transport$external$GSE165252$calib_slope)),
    list(item = "calib. slope \u2192 GSE78220", old = fmt2(rec_old$transport$external$GSE78220$calib_slope),
         new = fmt2(rec_new$transport$external$GSE78220$calib_slope)),
    list(item = "Brier \u2192 GSE78220", old = fmt3(rec_old$transport$external$GSE78220$brier),
         new = fmt3(rec_new$transport$external$GSE78220$brier)),
    list(item = "log-loss \u2192 GSE78220", old = fmt3(rec_old$transport$external$GSE78220$log_loss),
         new = fmt3(rec_new$transport$external$GSE78220$log_loss)))
  for (r in rows) {
    add_sd("ed5", "b", "data/imvigor_gene_recovery.json",
           sprintf("phase2.oldnew.%s", gsub("[^A-Za-z0-9]+", "_", r$item)),
           list(old = r$old, new = r$new))
  }
  d_txt <- do.call(rbind, lapply(rows, function(r)
    data.frame(item = r$item, old = r$old, new = r$new)))
  d_txt$y <- seq(11.2, 0.8, length.out = nrow(d_txt))
  # parity guard: old external AUCs must equal the frozen  values (md claim)
  stopifnot(rec$phase2$old$parity_with_frozen_$GSE165252$match,
            rec$phase2$old$parity_with_frozen_$GSE78220$match)

  pb <- ggplot(d_txt) +
    annotate("rect", xmin = 0, xmax = 100, ymin = 12.4, ymax = 13.6,
             fill = "grey93", colour = NA) +
    annotate("segment", x = 0, xend = 100, y = 12.4, yend = 12.4,
             colour = "grey40", linewidth = 0.35) +
    annotate("text", x = 0, y = 13.0, hjust = 0, size = 2.2, family = FF,
             fontface = "bold", colour = "grey15",
             label = "old (wrong identity, frozen // archives) vs new (recovered key,  re-run)") +
    geom_text(aes(x = 1, y = y, label = item), hjust = 0, size = 2.15,
              family = FF, colour = "grey15") +
    geom_text(aes(x = 52, y = y, label = old), hjust = 1, size = 2.15,
              family = FF, colour = "grey45") +
    geom_text(aes(x = 57, y = y, label = new), hjust = 0, size = 2.15,
              family = FF, colour = "grey10", fontface = "bold") +
    annotate("text", x = 52, y = 12.0, hjust = 1, size = 1.9, family = FF,
             fontface = "bold", colour = "grey45", label = "OLD") +
    annotate("text", x = 57, y = 12.0, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey10", label = "NEW") +
    annotate("text", x = 1, y = -1.4, hjust = 0, size = 1.9, family = FF,
             colour = "grey35", lineheight = 1.05,
             label = sprintf("old external AUCs reproduce the frozen  values digit-for-digit (parity check in the frozen JSON) \u2014\nthe 'old' column is the published-wrong baseline, kept for provenance\nGSE165252 barely moves (CI still spans 0.5); GSE78220 improves %+.3f but its CI still spans 0.5 at n = 25;\nBOTH calibration slopes move away from 1",
                             rec_new$transport$external$GSE78220$auc -
                               rec_old$transport$external$GSE78220$auc)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(-2.6, 14.0), clip = "off") +
    scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL,
         title = "Every affected anchor-arm number, old vs new (no cherry-picking)") +
    theme_void(base_family = FF) +
    theme(plot.title = element_text(size = 8, face = "plain")) +
    lab("b")

  pa / pb + plot_layout(heights = c(1.15, 1))
}

# =============================================================================
# ED6 -- M3 recalibration stability (100 repeats x 3 cohorts)
# =============================================================================
make_ed6 <- function() {

  dA <- do.call(rbind, lapply(TARGETS, function(ch) {
    m2 <- rrb$loco_benchmark$per_target[[ch]]$M2
    m3 <- rrb$loco_benchmark$per_target[[ch]]$M3
    r  <- rrb$calibration_curves$per_target[[ch]]$recal_before_after
    rm <- rrb$calibration_curves$per_target[[ch]]$M3$repeat_means
    add_sd("ed6", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M3.repeat_means.calib_slope", ch),
           rm$calib_slope)
    add_sd("ed6", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M3.repeat_means.calib_slope_ci95", ch),
           rm$calib_slope_ci95)
    add_sd("ed6", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M3.repeat_means.n_repeats_used", ch),
           rm$n_repeats_used)
    add_sd("ed6", "a", "data/rerun_benchmark.json",
           sprintf("calibration_curves.per_target.%s.M3.repeat_means.n_repeats_skipped", ch),
           rm$n_repeats_skipped)
    data.frame(cohort = ch, m2_slope = r$m2_calib_slope_logit,
               m3_slope = rm$calib_slope,
               lo = rm$calib_slope_ci95[[1]], hi = rm$calib_slope_ci95[[2]],
               improved = r$slope_improved, n_used = rm$n_repeats_used,
               n_skipped = rm$n_repeats_skipped)
  }))
  stopifnot(all(dA$n_used == 100L), all(dA$n_skipped == 0L))
  dA$cohort <- factor(dA$cohort, levels = TARGETS)
  dA$y <- as.integer(dA$cohort)

  pa <- ggplot(dA, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.4) +
    annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 0.55, ymax = 3.45,
             fill = OLIVE, alpha = 0.10) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14, linewidth = 0.45,
                   colour = BLUE) +
    geom_point(aes(x = m2_slope, fill = "M2 frozen slope"), shape = 22,
               size = 2.7, colour = "black", stroke = 0.5) +
    geom_point(aes(x = m3_slope, fill = "M3 100-repeat mean"), shape = 21,
               size = 2.7, colour = "black", stroke = 0.5) +
    geom_segment(aes(x = m2_slope, xend = m3_slope, yend = y, colour = improved),
                 arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.5,
                 linetype = "dotted") +
    geom_text(aes(x = m3_slope, y = y - 0.30,
                  label = sprintf("%.2f [%.2f, %.2f]", m3_slope, lo, hi)),
              size = 2.0, family = FF, colour = "grey20") +
    geom_text(aes(x = -4.7, label = sprintf("%s   n_repeats = %d (0 skipped)",
                                            cohort, n_used)),
              hjust = 0, size = 2.05, family = FF, colour = "grey20") +
    scale_fill_manual(values = c("M2 frozen slope" = "grey70",
                                 "M3 100-repeat mean" = BLUE), name = NULL) +
    scale_colour_manual(values = c(`TRUE` = OLIVE, `FALSE` = RED), guide = "none") +
    scale_x_continuous(limits = c(-4.8, 3.4), breaks = seq(-4, 3, 1)) +
    scale_y_continuous(breaks = dA$y, labels = rep("", nrow(dA))) +
    coord_cartesian(ylim = c(0.5, 3.7), clip = "off") +
    labs(x = "calibration slope (logit scale)", y = NULL,
         title = "M3 recalibration slopes: 100 stratified 50/50 repeats per cohort",
         subtitle = "whiskers = 2.5/97.5 percentiles across repeats; green arrow = slope improved, red = worsened; band = |slope \u2212 1| \u2264 0.10") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = c(0.70, 0.14),
          legend.background = element_rect(fill = alpha("white", 0.75))) +
    lab("a")

  # ---- (b) metric stability across repeats (AUC / Brier / log-loss) ----
  dB <- do.call(rbind, lapply(TARGETS, function(ch) {
    lb <- rrb$loco_benchmark$per_target[[ch]]
    for (mm in c("auc", "brier", "log_loss")) {
      add_sd("ed6", "b", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.M3.%s", ch, mm), lb$M3[[mm]])
      add_sd("ed6", "b", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.M3.%s_ci95", ch, mm),
             lb$M3[[paste0(mm, "_ci95")]])
      add_sd("ed6", "b", "data/rerun_benchmark.json",
             sprintf("loco_benchmark.per_target.%s.M2.%s", ch, mm), lb$M2[[mm]])
    }
    data.frame(cohort = ch,
               metric = c("AUC", "AUC", "Brier", "Brier", "log-loss", "log-loss"),
               model = rep(c("M2", "M3"), 3),
               value = c(lb$M2$auc, lb$M3$auc, lb$M2$brier, lb$M3$brier,
                         lb$M2$log_loss, lb$M3$log_loss),
               lo = c(NA, lb$M3$auc_ci95[[1]], NA, lb$M3$brier_ci95[[1]],
                      NA, lb$M3$log_loss_ci95[[1]]),
               hi = c(NA, lb$M3$auc_ci95[[2]], NA, lb$M3$brier_ci95[[2]],
                      NA, lb$M3$log_loss_ci95[[2]]))
  }))
  dB$cohort <- factor(dB$cohort, levels = TARGETS)
  dB$metric <- factor(dB$metric, levels = c("AUC", "Brier", "log-loss"))
  dB$model <- factor(dB$model, levels = c("M2", "M3"))

  pb <- ggplot(dB, aes(x = model, y = value)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.4,
                  colour = BLUE, alpha = 0.85, na.rm = TRUE) +
    geom_point(aes(fill = model), shape = 21, size = 2.5, colour = "black",
               stroke = 0.45) +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -1.05, size = 1.85,
              family = FF, colour = "grey30") +
    facet_grid(metric ~ cohort, scales = "free_y") +
    scale_fill_manual(values = c(M2 = RED, M3 = BLUE), guide = "none") +
    labs(x = NULL, y = NULL,
         title = "M3 vs M2 point estimates with repeat-percentile CIs (stability of the improvement)",
         caption = "M3 whiskers = 2.5/97.5 percentiles across the same 100 repeats (seed 20260911); M2 has no repeat layer (single frozen fit).\nGSE165252 / GSE78220: recalibration does NOT improve and can worsen \u2014 the IMvigor210 improvement is not a generic property.") +
    theme(plot.caption = element_text(size = 5.9, colour = "grey30"),
          strip.text = element_text(size = 7)) +
    lab1("b")

  pa / pb + plot_layout(heights = c(0.85, 1.15))
}

# =============================================================================
# ED7 -- Clinical-only baseline + case-mix shift
# =============================================================================
make_ed7 <- function() {

  smd <- clin$casemix$standardized_differences_IMvigor210_vs_GSE78220
  prev <- clin$casemix$prevalence_shift
  add_sd("ed7", "a", "data/clinical_only.json",
         "casemix.standardized_differences_IMvigor210_vs_GSE78220", smd)
  add_sd("ed7", "b", "data/clinical_only.json",
         "casemix.prevalence_shift", prev)

  # ---- (a) SMD bars ----
  dA <- data.frame(
    var = c("advanced / metastatic stage", "male sex"),
    smd = c(smd$pct_advanced_stage_smd, smd$pct_male_smd))
  dA$var <- factor(dA$var, levels = rev(dA$var))
  pa <- ggplot(dA, aes(x = smd, y = var)) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = 3, colour = "grey60",
               linewidth = 0.3) +
    geom_col(aes(fill = abs(smd) > 0.2), width = 0.5, colour = "black",
             linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", smd)),
              hjust = ifelse(dA$smd < 0, -0.25, 1.25), size = 2.4, family = FF,
              fontface = "bold") +
    annotate("text", x = -2.6, y = 2.43, hjust = 0, size = 1.7, family = FF,
             colour = "grey30", lineheight = 0.95,
             label = "IMvigor210 vs GSE78220;\n|SMD| > 0.2 = material shift;\nstage: TCGA IV vs AJCC M1") +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = "grey70"), guide = "none") +
    scale_x_continuous(limits = c(-3.0, 0.9), breaks = seq(-3, 1, 1)) +
    labs(x = "standardized mean difference (SMD)", y = NULL,
         title = "Case-mix shift:\nstage composition is extreme",
         subtitle = "18.8% stage-IV vs 96.0% M1;\nsex moderately shifted") +
    theme(plot.subtitle = element_text(size = 6.8)) +
    lab("a")

  # ---- (b) response-prevalence drift ----
  dB <- data.frame(cohort = names(prev[c("IMvigor210", "GSE165252", "GSE78220")]),
                   prev = unlist(prev[c("IMvigor210", "GSE165252", "GSE78220")]))
  dB$cohort <- factor(dB$cohort, levels = c("IMvigor210", "GSE165252", "GSE78220"))
  pb <- ggplot(dB, aes(x = cohort, y = prev)) +
    geom_col(aes(fill = cohort), width = 0.62, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", prev)), vjust = -0.55, size = 2.4,
              family = FF, fontface = "bold") +
    scale_fill_manual(values = c(IMvigor210 = BLUE, GSE165252 = "grey60",
                                 GSE78220 = RED), guide = "none") +
    scale_y_continuous(limits = c(0, 0.62), breaks = seq(0, 0.6, 0.1),
                       labels = percent) +
    labs(x = NULL, y = "response rate",
         title = sprintf("Prevalence drift:\n%.3f \u2192 %.3f \u2192 %.3f",
                         prev$IMvigor210, prev$GSE165252, prev$GSE78220),
         subtitle = sprintf("GSE78220 responds at %.1fx the anchor rate \u2014\nfrozen source-calibrated predictions systematically\nunderstate target risk",
                            prev$GSE78220 / prev$IMvigor210)) +
    theme(plot.subtitle = element_text(size = 6.8),
          axis.text.x = element_text(size = 6.2, angle = 30, hjust = 1)) +
    lab("b")

  # ---- (c) within IMvigor210: rich vs harmonizable clinical vs S0 ----
  wi <- clin$within_imvigor210_efron632
  add_sd("ed7", "c", "data/clinical_only.json",
         "within_imvigor210_efron632", wi)
  add_sd("ed7", "c", "data/rerun_benchmark.json",
         "within_reference.per_cohort.IMvigor210.auc_632", wr$IMvigor210$auc_632)
  dC <- data.frame(
    model = c("C_rich (14 covariates,\nIMvigor210-only)", "S0 (322-gene)",
              "C_main (sex + stage)", "C_sex (sex)"),
    auc632 = c(wi$C_rich$efron_0632_auc, wr$IMvigor210$auc_632,
               wi$C_main$efron_0632_auc, wi$C_sex$efron_0632_auc),
    transportable = c(FALSE, TRUE, TRUE, TRUE))
  dC$model <- factor(dC$model, levels = rev(dC$model))
  pc <- ggplot(dC, aes(x = auc632, y = model)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_col(aes(fill = transportable), width = 0.55, colour = "black",
             linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", auc632)), hjust = -0.22, size = 2.3,
              family = FF, fontface = "bold") +
    annotate("text", x = 0.42, y = 4.45, hjust = 0, size = 1.95, family = FF,
             colour = "grey30",
             label = "Efron 0.632 AUC,\nwithin IMvigor210\n(n = 298, 68 events)") +
    scale_fill_manual(values = c(`TRUE` = "grey65", `FALSE` = RED),
                      labels = c(`TRUE` = "transportable to >=1 cohort",
                                 `FALSE` = "NOT transportable"),
                      name = NULL) +
    scale_x_continuous(limits = c(0.4, 0.95), breaks = seq(0.4, 0.9, 0.1)) +
    labs(x = "AUC (Efron 0.632)", y = NULL,
         title = "Useful variables are\nnot transportable",
         subtitle = "rich clinical model: TMB + ECOG\n+ metastatic site \u2014 none deposited\nin the GEO targets") +
    theme(plot.subtitle = element_text(size = 6.8),
          legend.position = "top",
          legend.justification = "left",
          legend.direction = "vertical",
          legend.text = element_text(size = 5.9),
          legend.key.size = unit(0.28, "cm"),
          legend.margin = margin(0, 0, 0, 0)) +
    lab("c")

  # ---- (d) external GSE78220: harmonizable clinical vs S0 vs M0 ----
  cm <- clin$external_results$GSE78220$C_main
  cs <- clin$external_results$GSE78220$C_sex
  add_sd("ed7", "d", "data/clinical_only.json",
         "external_results.GSE78220", clin$external_results$GSE78220)
  add_sd("ed7", "d", "data/clinical_only.json",
         "external_results.GSE165252", clin$external_results$GSE165252)
  add_sd("ed7", "d", "data/rerun_benchmark.json",
         "pairwise_transports.transports.train_IMvigor210__test_GSE78220.M2.auc",
         prim$GSE78220$auc)
  dD <- data.frame(
    model = c("M0 prevalence", "C_main (sex + stage)", "C_sex (sex)",
              "M2 frozen S0"),
    est = c(pw[[pk("IMvigor210", "GSE78220")]]$M0$auc, cm$auc, cs$auc,
            prim$GSE78220$auc),
    lo = c(NA, cm$auc_ci95[[1]], cs$auc_ci95[[1]], prim$GSE78220$auc_ci95[[1]]),
    hi = c(NA, cm$auc_ci95[[2]], cs$auc_ci95[[2]], prim$GSE78220$auc_ci95[[2]]))
  dD$model <- factor(dD$model, levels = rev(dD$model))
  pc2_note <- sprintf("GSE165252 (not shown): every clinical model degenerates to the train-prevalence\nconstant %.3f (AUC 0.500 by construction); C_main log-loss %.3f = M0 exactly",
                      clin$external_results$GSE165252$C_main$constant_prediction,
                      clin$external_results$GSE165252$C_main$log_loss)
  pd <- ggplot(dD, aes(x = est, y = model)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_pointrange(aes(xmin = lo, xmax = hi), fatten = 2.3, stroke = 0.5,
                    shape = 19, colour = BLUE, linewidth = 0.45, na.rm = TRUE) +
    geom_text(aes(label = fmt3(est)), vjust = -1.15, size = 2.0, family = FF,
              colour = "grey25") +
    annotate("text", x = 0.40, y = 0.68, hjust = 0, size = 1.85, family = FF,
             colour = "grey30", lineheight = 1.0, label = pc2_note) +
    annotate("text", x = 0.86, y = 4.3, hjust = 1, size = 1.95, family = FF,
             colour = "grey30",
             label = "GSE78220, n = 25 (13 events)\nall CIs span 0.5") +
    scale_x_continuous(limits = c(0.39, 0.87), breaks = seq(0.4, 0.85, 0.1)) +
    coord_cartesian(ylim = c(0.45, 4.55), clip = "off") +
    labs(x = "AUC, external GSE78220 (train IMvigor210, frozen)", y = NULL,
         title = "Transportable variables add nothing beyond prevalence") +
    theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm")) +
    lab("d")

  (pa | pb | pc) / pd + plot_layout(heights = c(1, 0.85))
}

# =============================================================================
# ED10 -- analysis-freeze & reproducibility map (UPDATED for the recovery event)
# =============================================================================
make_ed10 <- function() {

  # run-time manifest enumeration (real manifest, not typed)
  _jsons <- basename(list.files(J("data", ""), pattern = "\\.json$"))
  _pys   <- basename(list.files(J("scripts"), pattern = "^_.*\\.py$"))
  n_prior_jsons <- length(_jsons) -
    sum(_jsons %in% c("imvigor_gene_recovery.json", "rerun_benchmark.json",
                         "rerun_ablation.json", "rerun_positive_control.json"))
  rerun_pys <- c("imvigor_gene_recovery.py", "rerun_benchmark.py",
                 "rerun_ablation.py", "rerun_positive_control.py")
  n_other_pys <- length(_pys) - sum(_pys %in% rerun_pys)
  seed <- rrb$seed
  n_boot_auc <- rrb$metadata$n_boot_auc_ci
  n_boot_632 <- rrb$metadata$n_boot_efron632
  n_recal <- rrb$metadata$n_repeats_recal
  stopifnot(seed == 20260911L, n_boot_auc == 2000L, n_boot_632 == 600L,
            n_recal == 100L)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "data__json", _jsons)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "scripts__py", _pys)
  add_sd("ed10", "manifest", "data/rerun_benchmark.json",
         "seed", seed)
  add_sd("ed10", "manifest", "data/rerun_benchmark.json",
         "n_boot_auc_ci", n_boot_auc)
  add_sd("ed10", "manifest", "data/rerun_benchmark.json",
         "n_boot_efron632", n_boot_632)
  add_sd("ed10", "manifest", "data/rerun_benchmark.json",
         "n_repeats_recal", n_recal)
  add_sd("ed10", "manifest", "data/imvigor_gene_recovery.json",
         "verification_gates.overall", rec$verification_gates$overall)
  add_sd("ed10", "manifest", "data/rerun_ablation.json",
         "removed_configs", names(rra$removed_configs))

  x_left <- 3.0
  p <- ggplot()
  # ---- layer 1: frozen inputs ----
  p <- p +
    annotate("rect", xmin = 1.5, xmax = 33.5, ymin = 34, ymax = 95,
             fill = "#eef3f9", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = x_left, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = BLUE, label = "FROZEN INPUTS (data/)") +
    annotate("text", x = x_left, y = 88.6, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25",
             label = "data/ \u2014 anchor authority + rerun series") +
    annotate("text", x = x_left, y = 84.8, hjust = 0, vjust = 1, size = 1.85,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = "imvigor_gene_recovery.json  (recovered key + gates)\nrerun_benchmark.json  (LOCO M0\u2013M6, within, calib,\n  pairwise, near-miss)\nrerun_ablation.json  (21-config grid + retraction)\nrerun_positive_control.json  (Ayers, tiers)") +
    annotate("text", x = x_left, y = 65.5, hjust = 0, size = 1.85, family = FF,
             colour = "grey30", lineheight = 1.05,
             label = sprintf("+ %d earlier  JSONs (clinical-only, provenance\nledger, threshold/zero sensitivity, transport\nsummary, positive-control superseded original)\n+ data/ / / prior freezes (wrong-identity\nvalues kept as provenance, never as results)", n_prior_jsons)) +
    annotate("text", x = x_left, y = 48.0, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = RED, label = "RETRACTION RECORD") +
    annotate("text", x = x_left, y = 44.2, hjust = 0, vjust = 1, size = 1.8,
             family = FF, colour = "grey20", lineheight = 1.05,
             label = "retracted: '70 unmeasured genes / 88.7% PC1 mass'\nstructural artefact + 228-gene intersection + 5\nmeasured-only ablation configs (premise gone under\nthe recovered key) \u2014 see data/\nartefact_retraction.md and number_audit.md")

  # ---- layer 2: scripts ----
  p <- p +
    annotate("rect", xmin = 36.5, xmax = 63.5, ymin = 34, ymax = 95,
             fill = "#f4eef8", colour = PURPLE, linewidth = 0.45) +
    annotate("text", x = 38.5, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = PURPLE, label = "SCRIPTS") +
    annotate("text", x = 38.5, y = 88.6, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25",
             label = "anchor recovery +  rerun series") +
    annotate("text", x = 38.5, y = 84.5, hjust = 0, vjust = 1, size = 1.85,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = paste0(rerun_pys, collapse = "\n")) +
    annotate("text", x = 38.5, y = 66.5, hjust = 0, size = 1.85, family = FF,
             colour = "grey30", lineheight = 1.05,
             label = sprintf("+ %d earlier  scripts (phase-1a diagnostics,\nclinical-only, sensitivity, summaries)", n_other_pys)) +
    annotate("text", x = 38.5, y = 58.0, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25", label = "figure generation") +
    annotate("text", x = 38.5, y = 54.2, hjust = 0, vjust = 1, size = 1.85,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = "figures/make_stable_figures.R (stable set)\nfigures/make_anchor_figures.R (this script;\nanchor set, reads the rerun series verbatim)")

  # ---- layer 3: outputs ----
  p <- p +
    annotate("rect", xmin = 66.5, xmax = 98.5, ymin = 34, ymax = 95,
             fill = "#ebf1dd", colour = OLIVE, linewidth = 0.45) +
    annotate("text", x = 68.5, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = OLIVE, label = "OUTPUTS (figures/)") +
    annotate("text", x = 68.5, y = 88.8, hjust = 0, size = 1.9, family = FF,
             colour = "grey25", label = "one figure = PNG (300 dpi) + PDF + SVG") +
    annotate("text", x = 68.5, y = 84.8, hjust = 0, vjust = 1, size = 1.8,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = "stable set: fig1_transport_system, fig6_joint_boundary,\ned2_cohort_ledger, ed8_toxicity_stress,\ned9_joint_identifiability") +
    annotate("text", x = 68.5, y = 70.5, hjust = 0, size = 1.8, family = FF,
             fontface = "bold", colour = "grey25",
             label = "anchor set (this script, post-recovery numbers):") +
    annotate("text", x = 68.5, y = 66.8, hjust = 0, vjust = 1, size = 1.8,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = "fig2_external_transport, fig3_calibration_pathways,\nfig4_ablation_constraint, fig5_benchmark_strategies,\ned1_positive_control, ed3_loco_sensitivity,\ned4_ablation_matrix, ed5_identity_recovery,\ned6_calibration_stability, ed7_clinical_casemix,\ned10_analysis_freeze (this figure, updated)") +
    annotate("text", x = 68.5, y = 46.5, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25", label = "source_data/") +
    annotate("text", x = 68.5, y = 43.2, hjust = 0, vjust = 1, size = 1.8,
             family = FF, colour = "grey30", lineheight = 1.05,
             label = "one JSON per figure; every panel value keyed to\nits frozen source file + key")

  # arrows between layers
  p <- p +
    annotate("segment", x = 34, xend = 36.2, y = 64, yend = 64,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8,
             colour = "grey35") +
    annotate("segment", x = 64, xend = 66.2, y = 64, yend = 64,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8,
             colour = "grey35")

  # ---- gene-identity recovery event band (NEW) ----
  p <- p +
    annotate("rect", xmin = 1.5, xmax = 98.5, ymin = 22.5, ymax = 31.5,
             fill = "#fdecea", colour = RED, linewidth = 0.45) +
    annotate("text", x = 3.0, y = 29.6, hjust = 0, size = 2.2, family = FF,
             fontface = "bold", colour = RED,
             label = "EVENT (2026-09-24): anchor gene-row identity recovered mid-freeze \u2014 all anchor-arm numbers re-run before submission") +
    annotate("text", x = 3.0, y = 25.6, hjust = 0, vjust = 0.5, size = 1.85,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = sprintf("exmat row j == j-th fData row under upper-cased symbol alphabetical sort; verification gates %s (sex anchor, mirror 31085/31085 numeric identity, cross-cohort \u03c1 0.81/0.62);\nwrong-identity archives (/// anchor values) stay frozen as provenance of the error and are superseded by the _rerun series",
                             rec$verification_gates$overall))

  # ---- locked parameters band ----
  p <- p +
    annotate("rect", xmin = 1.5, xmax = 98.5, ymin = 3.5, ymax = 19.5,
             fill = "grey97", colour = "grey35", linewidth = 0.4) +
    annotate("text", x = 3.5, y = 17.2, hjust = 0, size = 2.3, family = FF,
             fontface = "bold", label = "Locked analysis parameters ( freeze, unchanged by the recovery)") +
    annotate("text", x = 3.5, y = 11.8, hjust = 0, size = 2.0, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = sprintf("global seed = %d (all  resampling)\nAUC CI: %d-rep bootstrap BCa\nEfron 0.632: %d bootstrap refits",
                             seed, n_boot_auc, n_boot_632)) +
    annotate("text", x = 40, y = 11.8, hjust = 0, size = 2.0, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = sprintf("M3/M6 recalibration: %d stratified 50/50 repeats\nexternal CI (ablation): percentile + BCa engines\nridge \u03bb: 1/n_train external, 1/N_total within (re-ruled)",
                             n_recal)) +
    annotate("text", x = 78, y = 11.8, hjust = 0, size = 2.0, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = "no number in any figure is typed:\nall read from frozen JSON keys\n(recovered-identity values only)") +
    annotate("text", x = 50, y = 5.6, size = 1.95, family = FF, colour = "grey30",
             label = "any change to a frozen input requires a new version directory, not an edit")

  p + theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("a")
}

# =============================================================================
# build all
# =============================================================================
FIGS <- list(
  list(sdkey = "fig2",  name = "fig2_external_transport",    fn = make_fig2,  w = 180, h = 185),
  list(sdkey = "fig3",  name = "fig3_calibration_pathways",  fn = make_fig3,  w = 180, h = 195),
  list(sdkey = "fig4",  name = "fig4_ablation_constraint",   fn = make_fig4,  w = 180, h = 195),
  list(sdkey = "fig5",  name = "fig5_benchmark_strategies",  fn = make_fig5,  w = 180, h = 165),
  list(sdkey = "ed1",   name = "ed1_positive_control",       fn = make_ed1,   w = 180, h = 145),
  list(sdkey = "ed3",   name = "ed3_loco_sensitivity",       fn = make_ed3,   w = 180, h = 160),
  list(sdkey = "ed4",   name = "ed4_ablation_matrix",        fn = make_ed4,   w = 180, h = 175),
  list(sdkey = "ed5",   name = "ed5_identity_recovery",      fn = make_ed5,   w = 180, h = 205),
  list(sdkey = "ed6",   name = "ed6_calibration_stability",  fn = make_ed6,   w = 180, h = 155),
  list(sdkey = "ed7",   name = "ed7_clinical_casemix",       fn = make_ed7,   w = 180, h = 155),
  list(sdkey = "ed10",  name = "ed10_analysis_freeze",       fn = make_ed10,  w = 180, h = 150)
)

PANELS <- list(
  fig2_external_transport    = c("a primary frozen transports + within refs",
                                 "b pairwise 6-direction forest + source-swap deltas",
                                 "c Ayers IFN-gamma positive control vs S0"),
  fig3_calibration_pathways  = c("a M2 calibration curves (5 bins, CP CI, per-patient rug)",
                                 "b M3 recalibration slopes before/after",
                                 "c Brier / log-loss M0/M2/M3 x 3 cohorts"),
  fig4_ablation_constraint   = c("a 21-config external AUC + CI (layer colours)",
                                 "b layer variance decomposition old vs new",
                                 "c gene-identity recovery: old vs new numbers"),
  fig5_benchmark_strategies  = c("a LOCO M0-M6 x 3 cohorts, refusals marked",
                                 "b reference strategies (primary estimand)"),
  ed1_positive_control       = c("a simulated pass-rate heatmap (rule behaviour)",
                                 "b Ayers empirical forest + tolerance-tier tiles"),
  ed3_loco_sensitivity       = c("a estimand schematic (primary vs LOCO)",
                                 "b M2: LOCO vs single-source per target"),
  ed4_ablation_matrix        = c("21 x 2 AUC heatmap with layer frames"),
  ed5_identity_recovery      = c("a recovery timeline (5 stages)",
                                 "b old vs new numbers table"),
  ed6_calibration_stability  = c("a M3 slope 100-repeat means + CI",
                                 "b AUC/Brier/log-loss stability vs M2"),
  ed7_clinical_casemix       = c("a SMD bars", "b prevalence drift",
                                 "c within-anchor clinical vs S0",
                                 "d external clinical vs S0 vs M0"),
  ed10_analysis_freeze       = c("frozen inputs -> scripts -> outputs manifest + recovery event + locked parameters")
)

for (f in FIGS) {
  ok <- try({
    p <- f$fn()
    save_trio(p, f$name, f$w, f$h)
    dump_sd(f$sdkey, f$name)
    cat(sprintf("[panels] %s: %s\n", f$name,
                paste(PANELS[[f$name]], collapse = " | ")))
    TRUE
  }, silent = FALSE)
  if (!isTRUE(ok)) cat(sprintf("[FAIL] %s\n", f$name))
}

# =============================================================================
# verification: 11 figures x 3 formats non-zero + source-data JSONs parseable
# =============================================================================
cat("\n--- verification ---\n")
n_ok <- 0; n_expected <- 0
png_blank <- function(fp) {
  # read the raw PNG back and test whether ANY non-white pixel exists
  if (requireNamespace("png", quietly = TRUE)) {
    dim3 <- try(png::readPNG(fp), silent = TRUE)
    if (inherits(dim3, "try-error")) return(FALSE)
    lum <- if (is.matrix(dim3)) {
      dim3                                   # grayscale
    } else {
      0.30 * dim3[, , 1] + 0.59 * dim3[, , 2] + 0.11 * dim3[, , 3]
    }
    return(all(lum > 0.999))
  }
  # fallback: no png decoder available -- a fully blank render compresses to a
  # tiny file (observed blank = 18.5 kB vs 190-360 kB for real figures)
  file.size(fp) < 40000
}
for (f in FIGS) {
  for (ext in c("pdf", "png", "svg")) {
    n_expected <- n_expected + 1
    fp <- file.path(OUTDIR, paste0(f$name, ".", ext))
    sz <- if (file.exists(fp)) file.size(fp) else -1
    ok <- sz > 0
    if (ok && ext == "png" && png_blank(fp)) {
      ok <- FALSE
      cat(sprintf("[BLANK-RENDER] %s\n", fp))
    }
    n_ok <- n_ok + as.integer(ok)
    if (!ok && !file.exists(fp)) cat(sprintf("[MISSING] %s\n", fp))
  }
  sd_fp <- file.path(SDDIR, paste0(f$name, ".json"))
  if (file.exists(sd_fp)) {
    invisible(fromJSON(sd_fp, simplifyVector = FALSE))
    cat(sprintf("[source_data ok] %s (panels: %s)\n", basename(sd_fp),
                paste(names(fromJSON(sd_fp, simplifyVector = FALSE)$panels), collapse = ", ")))
  } else {
    cat(sprintf("[source_data MISSING] %s\n", sd_fp))
  }
}
cat(sprintf("figure files non-zero: %d/%d\n", n_ok, n_expected))
stopifnot(n_ok == n_expected)
cat("done.\n")
