#!/usr/bin/env Rscript
# =============================================================================
# make_stable_figures.R -- Nature Cancer Analysis  STABLE figure set
#
# Five figures that do NOT depend on the IMvigor210 gene-row identity (the
# anchor-arm gene identity is being restored separately; the remaining main
# figures will be produced after that settles). All statistical numbers are
# read from frozen JSON artifacts under data/ -- nothing data-derived is
# hard-coded; only geometry (box coordinates, dodge widths, annotation
# positions) is specified inline.
#
# Figures
#   Fig.1  fig1_transport_system    Transportability is a property of the full
#                                   prediction system
#          (a) construct -> measurement -> predictor -> outcome estimand ->
#              target population -> transport verdict chain (programmatic)
#          (b) primary estimand definition (IMvigor210 frozen S0 -> independent
#              target, no re-estimation; LOCO = different estimand boundary)
#          (c) cohort map, 9 cohorts (-corrected provenance)
#   Fig.6  fig6_joint_boundary      Joint efficacy-toxicity evidence boundary
#          (GSE186143, n = 47 joint-evaluable)
#          (a) four-quadrant occupancy + Clopper-Pearson 95% CI per cell
#          (b) 25-fold per-fold log-loss, joint multinomial vs independent
#          (c) identifiability (9 vs 6 effective parameters, ...)
#          (d) future paired-bridge design (programmatic schematic)
#   ED2    ed2_cohort_ledger        cohort provenance/endpoint ledger table
#   ED8    ed8_toxicity_stress      toxicity stress-test forest (0.5 ref line)
#   ED9    ed9_joint_identifiability fold occupancy / identifiability detail
#   ED10   ed10_analysis_freeze     analysis-freeze & reproducibility map
#
# Outputs: figures/{name}.{png,pdf,svg} (PNG 300 dpi, 180 mm wide)
#          figures/source_data/{name}.json (per-panel machine-readable)
#
# Frozen sources (read verbatim, asserted):
#   data/provenance_ledger.json        (cohort provenance, corrected)
#   data/transport_summary_joint.json  (p0_10 within refs; p0_11/ed9)
#   data/joint_stress_test.json        (per-fold deltas; tox recheck)
#   data/tox_meta.json                     (toxicity meta per-cohort)
#   data/disease_coverage.json         (cohort roster)
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

F_LEDGER <- J("data", "", "provenance_ledger.json")
F_TSJ    <- J("data", "", "transport_summary_joint.json")
F_JOINT  <- J("data", "", "joint_stress_test.json")
F_TOX60  <- J("data", "tox_meta.json")
F_DISC   <- J("data", "", "disease_coverage.json")

led    <- fromJSON(F_LEDGER, simplifyVector = FALSE)
tsj    <- fromJSON(F_TSJ,    simplifyVector = FALSE)
joint  <- fromJSON(F_JOINT,  simplifyVector = FALSE)
tox60  <- fromJSON(F_TOX60,  simplifyVector = FALSE)
disc   <- fromJSON(F_DISC,   simplifyVector = FALSE)

# ---------------------------------------------------------------- style ------
# mirrors make_figures.R so the  stable set is visually consistent
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

PAL  <- c("#1f4e79", "#c0504d", "#4f6228", "#7f5fa2")  # blue/red/olive/purple
BLUE <- "#1f4e79"; RED <- "#c0504d"; OLIVE <- "#4f6228"; PURPLE <- "#7f5fa2"

# bold lowercase panel tag, top-left, clipped-in ( convention)
lab <- function(txt) annotate("text", x = -Inf, y = Inf, label = txt,
                              hjust = -0.35, vjust = 2.2, size = 3.4,
                              fontface = "bold", family = FF)

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
  out <- list(figure = name, generated_by = "figures/make_stable_figures.R", panels = panels)
  write(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"),
        file.path(SDDIR, paste0(name, ".json")))
}

# ------------------------------------------------------------ save trio ------
SVG_OK <- requireNamespace("svglite", quietly = TRUE)
save_trio <- function(p, name, w_mm, h_mm) {
  png_file <- file.path(OUTDIR, paste0(name, ".png"))
  try({
    ggsave(png_file, p, width = w_mm, height = h_mm, units = "mm",
           dpi = 300, device = ragg::agg_png, bg = "white")
  }, silent = TRUE)
  if (!file.exists(png_file)) # fallback if ragg unavailable
    ggsave(png_file, p, width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white")
  pdf_file <- file.path(OUTDIR, paste0(name, ".pdf"))
  dev_lbl <- "quartz-pdf"
  # preferred on macOS: quartz PDF device with real Arial + unicode support
  try({
    ggsave(pdf_file, p, width = w_mm, height = h_mm, units = "mm",
           device = function(file, width, height, ...)
             grDevices::quartz(file = file, width = width, height = height,
                               type = "pdf", family = FF))
  }, silent = TRUE)
  if (!file.exists(pdf_file) || file.size(pdf_file) < 2000) {
    dev_lbl <- "pdf-Helvetica"
    try({
      ggsave(pdf_file, p, width = w_mm, height = h_mm, units = "mm",
             device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE)
    }, silent = TRUE)
  }
  svg_lbl <- "svglite"
  if (SVG_OK) {
    try({
      ggsave(file.path(OUTDIR, paste0(name, ".svg")), p, width = w_mm, height = h_mm,
             units = "mm", device = svglite::svglite, fix_text_size = FALSE, bg = "white")
    }, silent = TRUE)
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

# regex-extract helper: pull a sub-string out of a frozen JSON text field and
# stop hard if it is absent (numbers therefore remain sourced, not typed)
gx <- function(txt, pat, src) {
  m <- regmatches(txt, regexpr(pat, txt))
  if (length(m) != 1L || !nzchar(m))
    stop(sprintf("frozen-source assertion failed: pattern '%s' not found in %s", pat, src))
  m
}

# =============================================================================
# BLOCK 0 -- shared data preparation + assertions
# =============================================================================

# ---- 0.1 efficacy triad within-cohort references ( transport summary) ----
wr <- tsj$p0_10$within_reference
stopifnot(wr$IMvigor210$n == 298, wr$IMvigor210$n_event == 68,
          wr$GSE165252$n == 71,  wr$GSE165252$n_event == 23,
          wr$GSE78220$n == 25,   wr$GSE78220$n_event == 13)
for (ch in c("IMvigor210", "GSE165252", "GSE78220")) {
  add_sd("fig1", "b", "data/transport_summary_joint.json",
         sprintf("p0_10.within_reference.%s.n", ch), wr[[ch]]$n)
  add_sd("fig1", "b", "data/transport_summary_joint.json",
         sprintf("p0_10.within_reference.%s.n_event", ch), wr[[ch]]$n_event)
}

# ---- 0.2 ed9 joint block (GSE186143) ----------------------------------------
ed9 <- tsj$p0_11$ed9
stopifnot(ed9$n_joint_evaluable == 47L,
          ed9$cells$`DCB/severe` == 10L, ed9$cells$`DCB/no_severe` == 17L,
          ed9$cells$`NDB/severe` == 5L,  ed9$cells$`NDB/no_severe` == 15L,
          ed9$marginals$response_DCB == 27L, ed9$marginals$response_NDB == 20L,
          ed9$marginals$severe_irAE_yes == 15L, ed9$marginals$severe_irAE_no == 32L,
          ed9$n_folds_joint_better == 12L, ed9$n_folds_total == 25L,
          ed9$effective_params$joint_multinomial == 9L,
          ed9$effective_params$independent_pair == 6L,
          ed9$min_class_occupancy_per_fold == 1L,
          ed9$min_training_cell_across_folds == 4L,
          abs(ed9$pooled_log_loss$joint       - 1.3452) < 5e-4,
          abs(ed9$pooled_log_loss$independent - 1.3560) < 5e-4,
          abs(ed9$events_per_parameter - 5.222222) < 1e-4,
          abs(ed9$sign_permutation_p_on_frozen_deltas - 0.986260) < 1e-5,
          identical(ed9$fold_sizes_pattern, list(10L, 10L, 9L, 9L, 9L)) ||
            all(unlist(ed9$fold_sizes_pattern) == c(10, 10, 9, 9, 9)))

# per-fold deltas from the  frozen artifact (25 folds)
pf <- joint$joint_vs_independent$per_fold
stopifnot(length(pf) == 25L,
          joint$joint_vs_independent$n_folds_joint_better == 12L,
          joint$joint_vs_independent$n_folds_total == 25L,
          abs(joint$joint_vs_independent$log_loss_joint       - ed9$pooled_log_loss$joint) < 1e-9,
          abs(joint$joint_vs_independent$log_loss_independent - ed9$pooled_log_loss$independent) < 1e-9)
fold_df <- do.call(rbind, lapply(seq_along(pf), function(i) {
  add_sd("fig6", "b", "data/joint_stress_test.json",
         sprintf("joint_vs_independent.per_fold[%d].delta_joint_minus_independent", i - 1L),
         pf[[i]]$delta_joint_minus_independent)
  data.frame(fold = i, delta = pf[[i]]$delta_joint_minus_independent)
}))

# ---- 0.3 provenance ledger headline assertions ------------------------------
stopifnot(led$summary$corrected_items == 10L, led$summary$confirmed_items == 27L,
          led$summary$unverified_items == 3L)
f186 <- led$cohorts$GSE186143$fields
f166 <- led$cohorts$GSE166181$fields
f216 <- led$cohorts$GSE216329$fields
f253 <- led$cohorts$GSE253720$fields

# GSE186143 corrections quoted in the tasking, verified field-by-field:
n_chain_186 <- gx(f186$n_patients$verified_value, "61 melanoma patients",      "GSE186143 n_patients")
n_53_186    <- gx(f186$n_patients$verified_value, "the 53 = 26\\+27",          "GSE186143 n_patients")
n_47_186    <- gx(f186$n_patients$verified_value, "n=47 = 53 minus 6 NE",      "GSE186143 n_patients")
ici_186     <- paste(gx(f186$ici_regimen$verified_value, "anti-PD-1 monotherapy \\(29/53 labelled patients\\)", "GSE186143 ici_regimen"),
                     gx(f186$ici_regimen$verified_value, "anti-PD-1 \\+ anti-CTLA-4 combination \\(24/53\\)", "GSE186143 ici_regimen"),
                     sep = " / ")
stopifnot(f186$cancer_type$verdict == "CONFIRMED",
          grepl("Metastatic melanoma", f186$cancer_type$verified_value, fixed = TRUE),
          f186$n_patients$verdict == "CORRECTED",
          f186$ici_regimen$verdict == "CORRECTED",
          f186$baseline_sampling$verdict == "CONFIRMED")

# GSE216329 corrections
nsclc_216 <- gx(f216$cancer_type$verified_value, "18/24", "GSE216329 cancer_type")
irae_216  <- paste0(gx(f216$toxicity_endpoint$verified_value, "15 vs 9 patients", "GSE216329 toxicity_endpoint"),
                    " · ", gx(f216$toxicity_endpoint$verified_value, "grade 2-3", "GSE216329 toxicity_endpoint"))
stopifnot(f216$cancer_type$verdict == "CORRECTED",
          f216$modality$verdict == "CORRECTED",
          grepl("Single-cell RNA-seq", f216$modality$verified_value, fixed = TRUE),
          f216$toxicity_endpoint$verdict == "CORRECTED",
          grepl("NOT grade >= 3", f216$toxicity_endpoint$verified_value, fixed = TRUE))

# GSE253720 corrections
donors_253 <- gx(f253$n_patients$verified_value, "9 CPI-colitis",  "GSE253720 n_patients")
hc_253     <- gx(f253$n_patients$verified_value, "3 healthy controls", "GSE253720 n_patients")
uc_253     <- gx(f253$n_patients$verified_value, "2 ulcerative colitis", "GSE253720 n_patients")
n14_253    <- gx(f253$n_patients$verified_value, "14 donors", "GSE253720 n_patients")
stopifnot(f253$n_patients$verdict == "CORRECTED",
          f253$baseline_sampling$verdict == "CORRECTED",
          grepl("NOT pretreatment baseline", f253$baseline_sampling$verified_value, fixed = TRUE),
          grepl("steroid", f253$baseline_sampling$verified_value),
          f253$toxicity_endpoint$verdict == "CORRECTED")

# GSE166181 corrections
rn_166 <- gx(f166$n_patients$verified_value, "11 Responders, 9 Non-Responders", "GSE166181 n_patients")
stopifnot(f166$ici_regimen$verdict == "CORRECTED",
          grepl("nivolumab or pembrolizumab", f166$ici_regimen$verified_value),
          grepl("CD3\\+CD8\\+ T cells", f166$tissue_celltype$verified_value),
          f166$n_patients$verdict == "CONFIRMED")

# ---- 0.4 toxicity recheck + tox meta (ED8) -----------------------------------
rc <- joint$toxicity_cohort_recheck
stopifnot(rc$GSE189125$original_keys$cohort.n_patients == 13L,
          rc$GSE189125$original_keys$cohort.n_severe == 6L,
          rc$GSE189125$original_keys$cohort.n_no == 7L,
          rc$GSE189125$original_keys$verdict.n_closures == 5L,
          rc$GSE189125$original_keys$verdict.n_significant_p_mannwhitney == 2L,
          abs(rc$GSE189125$per_closure_aucs[[1]]$auc - 0.9048) < 1e-6)
g2 <- rc$GSE216329$original_keys__traceability
stopifnot(g2$n_total == 24L, g2$n_irAE == 15L, g2$n_no == 9L,
          abs(g2$auc_T_apparent - 0.6666667) < 1e-6,
          abs(g2$auc_T_632_corrected - 0.8042756) < 1e-6)
g3 <- rc$GSE253720$original_keys
stopifnot(g3$`auc_toxicity_CPI_colitis_vs_HC.n_patients` == 12L,
          g3$`auc_toxicity_CPI_colitis_vs_HC.n_CPI_colitis` == 9L,
          g3$`auc_toxicity_CPI_colitis_vs_HC.n_HC` == 3L,
          g3$`auc_toxicity_CPI_colitis_vs_HC.n_UC_excluded` == 2L,
          abs(g3$`auc_toxicity_CPI_colitis_vs_HC.point` - 0.918913) < 1e-6)
pc1 <- tox60$output$per_cohort[[1]]; pc2 <- tox60$output$per_cohort[[2]]
stopifnot(grepl("GSE216329", pc1$name), pc1$n == 24L, pc1$n_event == 15L, pc1$n_nonevent == 9L,
          grepl("GSE189125", pc2$name), pc2$n == 13L, pc2$n_event == 6L, pc2$n_nonevent == 7L,
          abs(tox60$output$corroborating_not_pooled$GSE189125_CD4Tem_univariate$auc - 0.905) < 1e-9,
          tox60$output$k == 2L)

# ---- 0.5 frozen locked parameters (ED10) -------------------------------------
SEED <- tsj$seed
stopifnot(SEED == 20260911L)
n_boot_outer <- tsj$p0_10$within_bootstrap$IMvigor210$n_boot_valid  # 2000
m_boot_inner <- gx(tsj$p0_10$delta_auc$IMvigor210$M2$method,
                   "600 internal OOB refits", "delta_auc.IMvigor210.M2.method")
sign_perm_note <- gx(tsj$p0_11$items[[12]]$note, "100,000 two-sided sign flips",
                     "p0_11.items[12].note")
stopifnot(n_boot_outer == 2000L)

# ---- 0.6 nine-cohort roster --------------------------------------------------
# roster basis: disease_coverage.json (13 local datasets) minus the four
# not-eligible / substrate-only series; the four audited cohorts are then
# overridden field-by-field with the  corrected provenance ledger.
elig <- Filter(function(c) !(c$id %in% c("GSE135222", "GSE249898", "GSE314698", "HCMI-CMDC")), disc$cohorts)
stopifnot(length(elig) == 9L)
roster_ids <- vapply(elig, function(c) c$id, "")
stopifnot(setequal(roster_ids, c("IMvigor210", "GSE165252", "GSE78220", "GSE166181",
                                 "GSE186143", "GSE189125", "GSE216329", "GSE253720",
                                 "GSE287541")))

AUDITED <- c("GSE186143", "GSE166181", "GSE216329", "GSE253720")
# per-cohort audit verdict counts, computed from the ledger (not typed)
verdict_counts <- lapply(AUDITED, function(id) {
  fl <- led$cohorts[[id]]$fields
  v <- vapply(fl, function(f) f$verdict, "")
  table(factor(v, levels = c("CONFIRMED", "CORRECTED", "UNVERIFIED")))
})
names(verdict_counts) <- AUDITED
stopifnot(sum(verdict_counts$GSE186143[["CORRECTED"]]) == 4L,
          sum(verdict_counts$GSE166181[["CORRECTED"]]) == 2L,
          sum(verdict_counts$GSE216329[["CORRECTED"]]) == 3L,
          sum(verdict_counts$GSE253720[["CORRECTED"]]) == 3L)

# cohort map data.frame (Fig.1c / ED2); every number traced above
cohort_map <- data.frame(
  id = roster_ids,
  arm = c("efficacy", "efficacy", "efficacy", "efficacy",
          "joint (E+T)", "toxicity", "toxicity", "toxicity", "available"),
  stringsAsFactors = FALSE
)
# order: efficacy block, toxicity block, joint, available
ord <- c("IMvigor210", "GSE165252", "GSE78220", "GSE166181",
         "GSE189125", "GSE216329", "GSE253720", "GSE186143", "GSE287541")
cohort_map <- cohort_map[match(ord, cohort_map$id), ]
rownames(cohort_map) <- NULL
# positional guard: arm MUST follow the id, never the list position
arm_of <- c(IMvigor210 = "efficacy", GSE165252 = "efficacy", GSE78220 = "efficacy",
            GSE166181 = "efficacy", GSE189125 = "toxicity", GSE216329 = "toxicity",
            GSE253720 = "toxicity", GSE186143 = "joint (E+T)",
            GSE287541 = "available")
cohort_map$arm <- unname(arm_of[cohort_map$id])
stopifnot(!any(is.na(cohort_map$arm)),
          cohort_map$arm[cohort_map$id == "GSE166181"] == "efficacy",
          cohort_map$arm[cohort_map$id == "GSE189125"] == "toxicity",
          cohort_map$arm[cohort_map$id == "GSE186143"] == "joint (E+T)",
          cohort_map$arm[cohort_map$id == "GSE287541"] == "available")

cm <- list(
  IMvigor210 = list(cancer = "urothelial", modality = "bulk RNA-seq (tumour)",
                    n = wr$IMvigor210$n, ev = wr$IMvigor210$n_event,
                    ici = "atezolizumab", sampling = "pretreatment baseline",
                    role = "source anchor (frozen S0)"),
  GSE165252  = list(cancer = "esophageal adenocarcinoma", modality = "bulk RNA-seq (tumour)",
                    n = wr$GSE165252$n, ev = wr$GSE165252$n_event,
                    ici = "atezolizumab", sampling = "pretreatment baseline",
                    role = "independent target"),
  GSE78220   = list(cancer = "melanoma", modality = "bulk RNA-seq (tumour)",
                    n = wr$GSE78220$n, ev = wr$GSE78220$n_event,
                    ici = "pembrolizumab", sampling = "pretreatment baseline",
                    role = "independent target"),
  GSE166181  = list(cancer = "metastatic melanoma", modality = "scRNA, FACS-sorted CD8 T",
                    n = 20L, ev = NA, ici = "nivolumab / pembrolizumab",
                    sampling = "T0 baseline (+T1/T2 on-tx)",
                    role = "efficacy (scRNA)"),
  GSE189125  = list(cancer = "melanoma", modality = "scRNA (5', PBMC)",
                    n = rc$GSE189125$original_keys$cohort.n_patients, ev = NA,
                    ici = "anti-PD-1 \u00b1 anti-CTLA-4", sampling = "pretreatment baseline",
                    role = "toxicity (CD4-Tem closure)"),
  GSE216329  = list(cancer = paste0("predominantly NSCLC (", nsclc_216, ")"),
                    modality = "scRNA (10x, PBMC T cells)",
                    n = pc1$n, ev = pc1$n_event,
                    ici = "anti-PD-1 (22/24 mono)", sampling = "baseline (_TN) + on-tx (_T)",
                    role = "toxicity"),
  GSE253720  = list(cancer = "CPI-colitis (mixed primaries)",
                    modality = "scRNA (colon biopsy)",
                    n = 14L, ev = NA,
                    ici = "anti-PD-1 6/9 \u00b7 combo 3/9",
                    sampling = "active steroid-exp./refractory colitis (NOT baseline)",
                    role = "toxicity (case-control probe)"),
  GSE186143  = list(cancer = "metastatic melanoma", modality = "bulk TPM (PBMC)",
                    n = ed9$n_joint_evaluable, ev = NA,
                    ici = "anti-PD-1 (29/53) \u00b7 combo (24/53)",
                    sampling = "pretreatment baseline PBMC",
                    role = "joint (E+T) stress test"),
  GSE287541  = list(cancer = "mixed (ChILI)", modality = "scRNA",
                    n = 28L, ev = NA, ici = "ICI", sampling = "longitudinal",
                    role = "available, not analysed")
)
for (id in cohort_map$id) {
  add_sd("fig1", "c", "data/provenance_ledger.json + data/transport_summary_joint.json + data/joint_stress_test.json + data/tox_meta.json + data/disease_coverage.json",
         sprintf("cohorts.%s.n", id), cm[[id]]$n)
  add_sd("fig1", "c", "same composite frozen sources",
         sprintf("cohorts.%s.cancer", id), cm[[id]]$cancer)
  add_sd("fig1", "c", "same composite frozen sources",
         sprintf("cohorts.%s.ici_regimen", id), cm[[id]]$ici)
}

# =============================================================================
# FIG 1 -- Transportability is a property of the full prediction system
# =============================================================================
make_fig1 <- function() {

  # ---- (a) system chain: six stages, programmatic boxes + arrows ----
  box_x <- seq(1, by = 16.6, length.out = 6)
  bw <- 14.6
  titles <- c("Construct", "Measurement", "Predictor", "Outcome\nestimand",
              "Target\npopulation", "Transport\nverdict")
  bodies <- c(
    "frozen immune-state\nconstruct S0\n(pre-specified\ngene-set definition)",
    "CLR \u2192 z \u2192 PCA\ntransport map\nbulk + scRNA",
    "ridge fit on source,\nthen frozen;\napplied\ntransform-only",
    "P(response | S0)\nin the target;\njoint (E,T) is a\nseparate estimand",
    "independent ICI\ncohorts\n(label-shifted:\ncase mix, platform)",
    "two-criterion gate:\nAUC CI excludes 0.5\nAND |slope \u2212 1|\n\u2264 0.10"
  )
  fills <- c("#dbe5f1", "#dbe5f1", "#f2dcdb", "#ebf1dd", "#ebf1dd", "grey97")
  cols  <- c(BLUE, BLUE, RED, OLIVE, OLIVE, "grey30")
  pa <- ggplot() +
    annotate("segment", x = 16.6 * 0:5 + 1, xend = 16.6 * 0:5 + 1 + 2.0,
             y = rep(46, 6), yend = rep(46, 6),
             arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.5,
             colour = "grey35") +
    annotate("text", x = 50, y = 97, size = 3.0, family = FF, fontface = "bold",
             label = "the transported object is the full prediction system, not a gene list") +
    annotate("text", x = 50, y = 91.5, size = 2.4, family = FF, colour = "grey35",
             label = "a weakness entering at ANY stage can break transport \u2014 the audit therefore scores the chain end-to-end")
  for (i in 1:6) {
    pa <- pa +
      annotate("rect", xmin = box_x[i], xmax = box_x[i] + bw, ymin = 22, ymax = 84,
               fill = fills[i], colour = cols[i], linewidth = 0.45) +
      annotate("text", x = box_x[i] + bw / 2, y = 77, size = 2.3, family = FF,
               fontface = "bold", colour = cols[i], lineheight = 0.9,
               label = titles[i]) +
      annotate("text", x = box_x[i] + bw / 2, y = 50, size = 1.9, family = FF,
               colour = "grey15", lineheight = 0.95, label = bodies[i])
  }
  pa <- pa +
    annotate("text", x = 50, y = 12, size = 2.4, family = FF, colour = "grey30",
             label = "Fig. 1a has no data content: it fixes the object whose transport is being claimed") +
    theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("a")

  # ---- (b) primary estimand definition ----
  pb <- ggplot() +
    annotate("rect", xmin = 2, xmax = 36, ymin = 30, ymax = 90,
             fill = "#dbe5f1", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = 19, y = 85, size = 2.7, family = FF, fontface = "bold",
             colour = BLUE, label = "SOURCE (frozen once)") +
    annotate("text", x = 19, y = 78.5, size = 2.15, family = FF, fontface = "bold",
             colour = "grey15",
             label = sprintf("IMvigor210 (n = %d, %d response events)",
                             wr$IMvigor210$n, wr$IMvigor210$n_event)) +
    annotate("text", x = 19, y = 55, size = 2.15, family = FF, colour = "grey15",
             lineheight = 1.0,
             label = "S0 construct defined\nPCA rotation + ridge fit\neverything estimated here,\nonce \u2014 then frozen\n(seed locked)") +
    annotate("segment", x = 37.5, xend = 56.5, y = 60, yend = 60,
             arrow = arrow(length = unit(0.30, "cm"), type = "closed"),
             linewidth = 1.5, colour = BLUE) +
    annotate("text", x = 47, y = 76, size = 2.0, family = FF, colour = "grey20",
             lineheight = 0.95,
             label = "frozen system\ntransported as-is\n(transform-only;\nno target labels\ntouched)") +
    annotate("rect", xmin = 58, xmax = 98, ymin = 30, ymax = 90,
             fill = "#f2dcdb", colour = RED, linewidth = 0.45) +
    annotate("text", x = 78, y = 85, size = 2.6, family = FF, fontface = "bold",
             colour = RED, label = "INDEPENDENT TARGETS (label-shifted)") +
    annotate("text", x = 78, y = 60, size = 2.15, family = FF, colour = "grey15",
             lineheight = 1.0, vjust = 0.5,
             label = sprintf("GSE165252 (n = %d, %d events)\nGSE78220 (n = %d, %d events)\n\nAUC \u00b7 calibration slope\nBrier \u00b7 log-loss\nscored with the frozen system",
                            wr$GSE165252$n, wr$GSE165252$n_event,
                            wr$GSE78220$n,  wr$GSE78220$n_event)) +
    annotate("rect", xmin = 2, xmax = 98, ymin = 4, ymax = 24,
             fill = "grey97", colour = "grey35", linewidth = 0.35) +
    annotate("text", x = 4, y = 20, hjust = 0, size = 2.45, family = FF,
             fontface = "bold", label = "Estimand boundary (pre-registered,  lock):") +
    annotate("text", x = 4, y = 8.2, hjust = 0, size = 2.1, family = FF,
             colour = "grey15", lineheight = 0.95, vjust = 0.5,
             label = "primary estimand = held-out performance of the frozen system in an independent target.\nLOCO / source-swapped refits answer a DIFFERENT question (demoted to sensitivity).\nM3-style recalibration consumes TARGET labels \u2014 an evaluation aid, not a transport.")
  pb <- pb + theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("b")

  # ---- (c) cohort map: 9 cohorts, size = n, fill = modality ----
  dC <- do.call(rbind, lapply(cohort_map$id, function(id) {
    x <- cm[[id]]
    ev_lab <- if (id %in% c("IMvigor210", "GSE165252", "GSE78220"))
      sprintf("%d events", x$ev)
    else if (id == "GSE166181") "11 R / 9 NR"
    else if (id == "GSE189125") "6 severe / 7 no"
    else if (id == "GSE216329") "15 irAE / 9 no"
    else if (id == "GSE253720") "9 CPI / 3 HC / 2 UC"
    else if (id == "GSE186143") "DCB 27 / NDB 20 (NE 6)"
    else "not analysed"
    data.frame(id = ifelse(id %in% AUDITED, paste0(id, " \u2020"), id),
               arm = cohort_map$arm[cohort_map$id == id],
               modality = x$modality,
               n = x$n, ev_lab = ev_lab)
  }))
  add_sd("fig1", "c", "data/provenance_ledger.json",
         "audit_marker", "dagger = fields re-verified against GEO SOFT + linked publications (see ED2)")
  dC$arm <- factor(dC$arm, levels = c("efficacy", "toxicity", "joint (E+T)", "available"))
  dC$modality <- factor(dC$modality, levels = c("bulk RNA-seq (tumour)", "bulk TPM (PBMC)",
                                                "scRNA (5', PBMC)", "scRNA (10x, PBMC T cells)",
                                                "scRNA, FACS-sorted CD8 T", "scRNA (colon biopsy)", "scRNA"))
  dC$y <- rev(seq_len(nrow(dC)))
  mod_cols <- c("bulk RNA-seq (tumour)" = BLUE, "bulk TPM (PBMC)" = BLUE,
                "scRNA (5', PBMC)" = "#3a78ae", "scRNA (10x, PBMC T cells)" = "#3a78ae",
                "scRNA, FACS-sorted CD8 T" = "#7fb2d6", "scRNA (colon biopsy)" = OLIVE,
                "scRNA" = "grey55")
  pc <- ggplot(dC, aes(x = arm, y = y)) +
    geom_point(aes(size = n, fill = modality), shape = 21, colour = "black",
               stroke = 0.5, alpha = 0.95) +
    geom_text(aes(label = sprintf("n = %d \u00b7 %s", n, ev_lab)),
              hjust = -0.08, size = 2.05, family = FF, colour = "grey15") +
    scale_size_continuous(range = c(2.5, 9), breaks = c(13, 25, 47, 71, 298), name = "n (patients/donors)") +
    scale_fill_manual(values = mod_cols, name = "modality") +
    scale_y_continuous(breaks = dC$y, labels = dC$id) +
    labs(x = NULL, y = NULL,
         title = "Cohort map: nine cohorts with endpoint roles (-corrected provenance)",
         caption = "\u2020 cohort fields re-verified against GEO SOFT records and linked publications in the  provenance audit; full field-level ledger in ED2.\nGSE186143 n = 47 joint-evaluable (53 labelled of 61 patients in series); GSE253720 = 14 donors (9 CPI-colitis + 3 HC + 2 UC).") +
    theme(axis.text.y = element_text(size = 6.8),
          plot.caption = element_text(size = 5.8, colour = "grey30"),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.direction = "vertical",
          legend.title.position = "top",
          legend.key.height = unit(0.26, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5)) +
    lab("c")

  pa / pb / pc + plot_layout(heights = c(1.05, 1.0, 1.35))
}

# =============================================================================
# FIG 6 -- Joint efficacy-toxicity evidence boundary (GSE186143, n = 47)
# =============================================================================
make_fig6 <- function() {

  cells <- ed9$cells
  cpc   <- ed9$per_cell_clopper_pearson_95

  # ---- (a) four-quadrant occupancy with Clopper-Pearson 95% CI ----
  dA <- data.frame(
    tox = factor(c("severe irAE", "no severe irAE", "severe irAE", "no severe irAE"),
                 levels = c("severe irAE", "no severe irAE")),
    eff = factor(c("DCB", "DCB", "NDB", "NDB"), levels = c("NDB", "DCB")),
    key = c("DCB/severe", "DCB/no_severe", "NDB/severe", "NDB/no_severe"))
  dA$n <- as.integer(c(cells$`DCB/severe`, cells$`DCB/no_severe`,
                       cells$`NDB/severe`, cells$`NDB/no_severe`))
  dA$prop <- c(cpc$`DCB/severe`$proportion, cpc$`DCB/no_severe`$proportion,
               cpc$`NDB/severe`$proportion, cpc$`NDB/no_severe`$proportion)
  dA$lo <- c(cpc$`DCB/severe`$cp95[[1]], cpc$`DCB/no_severe`$cp95[[1]],
             cpc$`NDB/severe`$cp95[[1]], cpc$`NDB/no_severe`$cp95[[1]])
  dA$hi <- c(cpc$`DCB/severe`$cp95[[2]], cpc$`DCB/no_severe`$cp95[[2]],
             cpc$`NDB/severe`$cp95[[2]], cpc$`NDB/no_severe`$cp95[[2]])
  for (k in dA$key) {
    add_sd("fig6", "a", "data/transport_summary_joint.json",
           sprintf("p0_11.ed9.per_cell_clopper_pearson_95.%s.n", k), cpc[[k]]$n)
    add_sd("fig6", "a", "data/transport_summary_joint.json",
           sprintf("p0_11.ed9.per_cell_clopper_pearson_95.%s.proportion", k), cpc[[k]]$proportion)
    add_sd("fig6", "a", "data/transport_summary_joint.json",
           sprintf("p0_11.ed9.per_cell_clopper_pearson_95.%s.cp95", k), cpc[[k]]$cp95)
  }
  add_sd("fig6", "a", "data/transport_summary_joint.json",
         "p0_11.ed9.n_joint_evaluable", ed9$n_joint_evaluable)
  add_sd("fig6", "a", "data/transport_summary_joint.json",
         "p0_11.ed9.marginals", ed9$marginals)

  dA$cell_lab <- sprintf("n = %d\n%.1f%%\n[%.1f, %.1f]", dA$n, 100 * dA$prop,
                         100 * dA$lo, 100 * dA$hi)
  pa <- ggplot(dA, aes(x = tox, y = eff, fill = prop)) +
    geom_tile(colour = "white", linewidth = 1.4) +
    geom_text(aes(label = cell_lab), size = 2.7, family = FF, fontface = "bold",
              colour = ifelse(dA$prop > 0.3, "white", "grey10"), lineheight = 0.95) +
    scale_fill_gradient(low = "#e8eff7", high = BLUE, labels = percent, name = "cell share") +
    annotate("text", x = 0.55, y = 3.02, hjust = 0, vjust = 1, size = 2.1, family = FF,
             colour = RED, lineheight = 0.95,
             label = "Clopper\u2013Pearson 95% CI per cell\nsmallest cell n = 5 (NDB & severe)") +
    scale_y_discrete(expand = expansion(add = c(0.25, 0.75))) +
    labs(x = "toxicity (T)", y = "efficacy (E)",
         title = sprintf("Joint (E,T) occupancy \u2014 GSE186143 (n = %d joint-evaluable)",
                         ed9$n_joint_evaluable),
         subtitle = sprintf("marginals: DCB %d / NDB %d \u00b7 severe %d / no severe %d \u00b7 OR (DCB | severe) = %.2f, Fisher exact p = %.2f",
                            ed9$marginals$response_DCB, ed9$marginals$response_NDB,
                            ed9$marginals$severe_irAE_yes, ed9$marginals$severe_irAE_no,
                            ed9$OR_E_T, ed9$fisher_exact_p_hand)) +
    theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.8)) +
    lab("a")
  # OR / Fisher values from the frozen ed9 block:
  stopifnot(abs(ed9$OR_E_T - 1.76470588) < 1e-6,
            abs(ed9$fisher_exact_p_hand - 0.5293323) < 1e-6)
  add_sd("fig6", "a", "data/transport_summary_joint.json",
         "p0_11.ed9.OR_E_T", ed9$OR_E_T)
  add_sd("fig6", "a", "data/transport_summary_joint.json",
         "p0_11.ed9.fisher_exact_p_hand", ed9$fisher_exact_p_hand)

  # ---- (b) per-fold log-loss, joint vs independent ----
  add_sd("fig6", "b", "data/transport_summary_joint.json",
         "p0_11.ed9.pooled_log_loss.joint", ed9$pooled_log_loss$joint)
  add_sd("fig6", "b", "data/transport_summary_joint.json",
         "p0_11.ed9.pooled_log_loss.independent", ed9$pooled_log_loss$independent)
  add_sd("fig6", "b", "data/transport_summary_joint.json",
         "p0_11.ed9.n_folds_joint_better", ed9$n_folds_joint_better)
  add_sd("fig6", "b", "data/transport_summary_joint.json",
         "p0_11.ed9.n_folds_total", ed9$n_folds_total)
  add_sd("fig6", "b", "data/transport_summary_joint.json",
         "p0_11.ed9.sign_permutation_p_on_frozen_deltas", ed9$sign_permutation_p_on_frozen_deltas)

  pb <- ggplot(fold_df, aes(x = fold, y = delta)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
    geom_col(aes(fill = delta < 0), width = 0.72, colour = NA) +
    scale_fill_manual(values = c(`TRUE` = OLIVE, `FALSE` = "grey70"), guide = "none") +
    annotate("text", x = 1, y = 0.46, hjust = 0, vjust = 1, size = 2.35,
             family = FF, colour = "grey20",
             label = sprintf("joint better in %d/%d folds", ed9$n_folds_joint_better, ed9$n_folds_total)) +
    annotate("text", x = 25.5, y = 0.40, hjust = 1, size = 2.25, family = FF,
             colour = "grey20", lineheight = 0.95,
             label = sprintf("pooled log-loss\njoint %.4f vs independent %.4f",
                             ed9$pooled_log_loss$joint, ed9$pooled_log_loss$independent)) +
    annotate("text", x = 25.5, y = 0.24, hjust = 1, size = 2.25, family = FF,
             colour = RED, fontface = "bold",
             label = sprintf("sign-permutation P = %.3f", ed9$sign_permutation_p_on_frozen_deltas)) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25)) +
    coord_cartesian(ylim = c(-0.3, 0.5)) +
    labs(x = "CV fold (5-fold \u00d7 5 repeats, seed 20260911)", y = "\u0394 log-loss (joint \u2212 independent)",
         title = "Joint multinomial buys nothing out-of-sample") +
    theme(plot.title = element_text(size = 8)) +
    lab("b")

  # ---- (c) identifiability ----
  add_sd("fig6", "c", "data/transport_summary_joint.json",
         "p0_11.ed9.effective_params", ed9$effective_params)
  add_sd("fig6", "c", "data/transport_summary_joint.json",
         "p0_11.ed9.events_per_parameter", ed9$events_per_parameter)
  add_sd("fig6", "c", "data/transport_summary_joint.json",
         "p0_11.ed9.min_training_cell_across_folds", ed9$min_training_cell_across_folds)
  add_sd("fig6", "c", "data/transport_summary_joint.json",
         "p0_11.ed9.min_class_occupancy_per_fold", ed9$min_class_occupancy_per_fold)

  dC <- data.frame(model = c("joint multinomial\n(4 classes)", "independent pair\n(2 logistics)"),
                   params = c(ed9$effective_params$joint_multinomial,
                              ed9$effective_params$independent_pair))
  dC$model <- factor(dC$model, levels = rev(dC$model))
  pc <- ggplot(dC, aes(x = model, y = params)) +
    geom_col(aes(fill = model), width = 0.52, colour = "black", linewidth = 0.35) +
    geom_text(aes(label = params), vjust = -0.6, size = 2.9, family = FF, fontface = "bold") +
    scale_fill_manual(values = c("independent pair\n(2 logistics)" = "grey70",
                                 "joint multinomial\n(4 classes)" = PURPLE),
                      guide = "none") +
    annotate("text", x = 1.5, y = 14.6, size = 2.35, family = FF, colour = "grey15",
             lineheight = 1.0,
             label = sprintf("events per parameter = %.1f (%d / %d)\n\u2014 far below the conventional 10",
                             ed9$events_per_parameter, ed9$n_joint_evaluable,
                             ed9$effective_params$joint_multinomial)) +
    annotate("text", x = 1.5, y = 12.0, size = 2.35, family = FF, colour = "grey15",
             lineheight = 1.0,
             label = sprintf("smallest cell n = %d (NDB & severe)\nsmallest training cell n = %d (every fold)",
                             ed9$cells$`NDB/severe`, ed9$min_training_cell_across_folds)) +
    coord_cartesian(ylim = c(0, 16)) +
    labs(x = NULL, y = "effective free parameters",
         title = "Identifiability of the joint estimand",
         subtitle = "no identifiable joint advantage at n = 47 \u2014 boundary, not discovery") +
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 6.8, colour = RED)) +
    lab("c")

  # ---- (d) future paired-bridge design (programmatic schematic) ----
  pd <- ggplot() +
    annotate("rect", xmin = 2, xmax = 98, ymin = 78, ymax = 97,
             fill = "grey93", colour = "grey45", linewidth = 0.35) +
    annotate("text", x = 50, y = 92.5, size = 2.25, family = FF,
             fontface = "bold", colour = "grey20",
             label = " boundary: E and T are measured in DIFFERENT cohorts") +
    annotate("text", x = 50, y = 84, size = 2.1, family = FF, colour = "grey30",
             lineheight = 0.95,
             label = "within-patient pairing unobservable \u2192 the joint (E,T) estimand\nis not estimable from any current retrospective cohort") +
    annotate("rect", xmin = 6, xmax = 30, ymin = 50, ymax = 72,
             fill = "#dbe5f1", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = 18, y = 65, size = 2.4, family = FF, fontface = "bold",
             colour = BLUE, label = "same patient") +
    annotate("text", x = 18, y = 56, size = 2.1, family = FF, colour = "grey15",
             lineheight = 0.95,
             label = "one baseline visit\n(prospective)") +
    annotate("rect", xmin = 40, xmax = 70, ymin = 66, ymax = 78,
             fill = "#ebf1dd", colour = OLIVE, linewidth = 0.4) +
    annotate("text", x = 55, y = 72, size = 2.15, family = FF, colour = "grey15",
             label = "tumour bulk RNA-seq \u2192 S0") +
    annotate("rect", xmin = 40, xmax = 70, ymin = 44, ymax = 56,
             fill = "#f2dcdb", colour = RED, linewidth = 0.4) +
    annotate("text", x = 55, y = 50, size = 2.15, family = FF, colour = "grey15",
             label = "PBMC scRNA \u2192 S0 (same draw)") +
    annotate("segment", x = 30, xend = 40, y = 64, yend = 72,
             arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.45, colour = "grey45") +
    annotate("segment", x = 30, xend = 40, y = 58, yend = 50,
             arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.45, colour = "grey45") +
    annotate("segment", x = 55, xend = 55, y = 44, yend = 36,
             arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.45, colour = "grey45") +
    annotate("rect", xmin = 26, xmax = 84, ymin = 22, ymax = 36,
             fill = "grey97", colour = "grey35", linewidth = 0.4) +
    annotate("text", x = 55, y = 29, size = 2.2, family = FF, colour = "grey15",
             lineheight = 0.95,
             label = "paired (E, T) outcomes per patient\n\u2192 within-patient E\u2013T coupling becomes estimable") +
    annotate("rect", xmin = 2, xmax = 98, ymin = 3, ymax = 16,
             fill = "#dbe5f1", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = 50, y = 11.5, size = 2.15, family = FF, colour = "grey10",
             lineheight = 0.95,
             label = "prospective multi-cancer cohort, paired n powered for the joint estimand") +
    annotate("text", x = 50, y = 6.5, size = 2.05, family = FF, colour = "grey15",
             label = "(the frozen prediction system is never re-estimated; the bridge supplies only the pairing)") +
    theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("d")

  (pa | pb) / (pc | pd) + plot_layout(heights = c(1.12, 1))
}

# =============================================================================
# ED2 -- cohort provenance / endpoint ledger (9 cohorts x fields)
# =============================================================================
make_ed2 <- function() {

  # field-level values; verdicts pulled from the ledger where audited
  # columns: cohort, cancer, modality, n(events), E endpoint, T endpoint,
  #          ICI regimen, sampling, role
  rows <- list(
    list(id = "IMvigor210", cancer = "urothelial", mod = "bulk (tumour)",
         n = sprintf("%d (%d ev)", wr$IMvigor210$n, wr$IMvigor210$n_event),
         ep = "response", tp = "\u2014", ici = "atezolizumab",
         smp = "baseline", role = "source anchor"),
    list(id = "GSE165252", cancer = "esoph. adeno.", mod = "bulk (tumour)",
         n = sprintf("%d (%d ev)", wr$GSE165252$n, wr$GSE165252$n_event),
         ep = "response", tp = "\u2014", ici = "atezolizumab",
         smp = "baseline", role = "target"),
    list(id = "GSE78220", cancer = "melanoma", mod = "bulk (tumour)",
         n = sprintf("%d (%d ev)", wr$GSE78220$n, wr$GSE78220$n_event),
         ep = "response", tp = "\u2014", ici = "pembrolizumab",
         smp = "baseline", role = "target"),
    list(id = "GSE166181", cancer = "m. melanoma", mod = "scRNA sort. CD8 T",
         n = "20 (11 R / 9 NR)", ep = "R vs NR (RECIST)", tp = "none annotated",
         ici = "nivo / pembro", smp = "T0 (+T1/T2)", role = "efficacy scRNA"),
    list(id = "GSE186143", cancer = "m. melanoma", mod = "bulk TPM (PBMC)",
         n = "61 / 53 / 47",
         ep = "DCB vs NDB (NE 6)\n[paper Suppl. Tab. 3]",
         tp = "severe irAE (grade \u22653)", ici = "anti-PD-1 29/53\n+ combo 24/53",
         smp = "baseline PBMC", role = "joint (E+T)"),
    list(id = "GSE189125", cancer = "melanoma", mod = "scRNA 5' (PBMC)",
         n = "13 (6/7)", ep = "\u2014", tp = "severe irAE",
         ici = "anti-PD-1 \u00b1 CTLA-4", smp = "baseline", role = "toxicity"),
    list(id = "GSE216329", cancer = "NSCLC-pred. (18/24)", mod = "scRNA 10x (T cells)",
         n = "24 (15/9)", ep = "none annotated",
         tp = "any irAE, grade 2-3,\norgan-specific", ici = "anti-PD-1 (22/24 mono)",
         smp = "baseline _TN +\non-tx _T", role = "toxicity"),
    list(id = "GSE253720", cancer = "CPI-colitis (mixed)", mod = "scRNA (colon biopsy)",
         n = "14 (9/3/2)",
         ep = "none annotated", tp = "CPI-colitis case def.\n(case-control)",
         ici = "anti-PD-1 6/9;\ncombo 3/9",
         smp = "active, steroid-exp.\n\u2014 NOT baseline", role = "toxicity (probe)"),
    list(id = "GSE287541", cancer = "mixed (ChILI)", mod = "scRNA",
         n = "28", ep = "\u2014", tp = "irAE (longitudinal)",
         ici = "ICI", smp = "longitudinal", role = "not yet analysed")
  )

  # ledger verdict mapping (which displayed field traces to which audited field)
  verdict_map <- list(
    GSE166181 = list(cancer = "cancer_type", mod = "modality", n = "n_patients",
                     ep = "efficacy_endpoint", tp = "toxicity_endpoint",
                     ici = "ici_regimen", smp = "baseline_sampling"),
    GSE186143 = list(cancer = "cancer_type", mod = "modality", n = "n_patients",
                     ep = "efficacy_endpoint", tp = "toxicity_endpoint",
                     ici = "ici_regimen", smp = "baseline_sampling"),
    GSE216329 = list(cancer = "cancer_type", mod = "modality", n = "n_patients_samples",
                     ep = "efficacy_endpoint", tp = "toxicity_endpoint",
                     ici = "ici_regimen", smp = "baseline_sampling"),
    GSE253720 = list(cancer = "cancer_type", mod = "modality", n = "n_patients",
                     ep = "efficacy_endpoint", tp = "toxicity_endpoint",
                     ici = "ici_regimen", smp = "baseline_sampling")
  )

  cols_df <- c(cohort = 1.0, cancer = 12.0, mod = 26.0, n = 39.5,
               ep = 51.0, tp = 66.5, ici = 82.5, smp = 98.0, role = 112.0)
  col_hjust <- c(cohort = 0, cancer = 0, mod = 0, n = 0, ep = 0, tp = 0,
                 ici = 0, smp = 0, role = 0)
  col_titles <- c(cohort = "Cohort", cancer = "Cancer", mod = "Modality",
                  n = "n (events)", ep = "Efficacy endpoint",
                  tp = "Toxicity endpoint", ici = "ICI regimen",
                  smp = "Sampling time", role = "Role")

  long <- do.call(rbind, lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    do.call(rbind, lapply(names(cols_df), function(f) {
      key <- if (f == "cohort") r$id else r[[f]]
      vd <- "NOT_AUDITED"
      if (r$id %in% AUDITED && f %in% names(verdict_map[[r$id]])) {
        fk <- verdict_map[[r$id]][[f]]
        if (!is.null(led$cohorts[[r$id]]$fields[[fk]]$verdict))
          vd <- led$cohorts[[r$id]]$fields[[fk]]$verdict
      }
      data.frame(row = i, field = f, x = cols_df[[f]], label = key, verdict = vd,
                 hjust = col_hjust[[f]])
    }))
  }))
  long$verdict <- factor(long$verdict,
                         levels = c("CONFIRMED", "CORRECTED", "NOT_AUDITED"))

  # audit-count column (right margin), computed from the ledger
  audit_lab <- sapply(rows, function(r) {
    if (r$id %in% AUDITED) {
      vc <- verdict_counts[[r$id]]
      sprintf("\u25cf %d / \u25b2 %d", vc[["CONFIRMED"]], vc[["CORRECTED"]])
    } else ""
  })
  add_sd("ed2", "table", "data/provenance_ledger.json",
         "summary.corrected_items", led$summary$corrected_items)
  add_sd("ed2", "table", "data/provenance_ledger.json",
         "summary.confirmed_items", led$summary$confirmed_items)
  add_sd("ed2", "table", "data/provenance_ledger.json",
         "summary.unverified_items", led$summary$unverified_items)
  for (id in AUDITED) {
    add_sd("ed2", "table", "data/provenance_ledger.json",
           sprintf("cohorts.%s.verdict_counts", id),
           list(CONFIRMED = unname(verdict_counts[[id]][["CONFIRMED"]]),
                CORRECTED = unname(verdict_counts[[id]][["CORRECTED"]]),
                UNVERIFIED = unname(verdict_counts[[id]][["UNVERIFIED"]])))
  }

  ytop <- length(rows) + 1.35
  # audit counts render as a small second line under the cohort name
  sub_df <- data.frame(x = 1.0, y = ytop - seq_along(rows) - 0.36,
                       label = audit_lab)
  p <- ggplot(long) +
    geom_text(aes(x = x, y = ytop - row, label = label, colour = verdict,
                  hjust = hjust),
              fontface = ifelse(long$verdict == "CORRECTED", "bold", "plain"),
              size = ifelse(long$field == "cohort", 2.2, 1.75), family = FF,
              lineheight = 0.9) +
    geom_text(data = sub_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 0.85, size = 1.55,
              family = FF, colour = "grey40") +
    geom_text(data = data.frame(x = cols_df, y = ytop + 0.45,
                                label = col_titles, field = names(cols_df)),
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 1.9, family = FF,
              fontface = "bold", colour = "grey15") +
    geom_hline(yintercept = ytop - 0.62, colour = "grey40", linewidth = 0.35) +
    annotate("segment", x = 0, xend = 124, y = 0.28, yend = 0.28,
             colour = "grey55", linewidth = 0.3) +
    scale_colour_manual(values = c(CONFIRMED = "grey10", CORRECTED = RED,
                                   NOT_AUDITED = "grey42"),
                        guide = "none") +
    annotate("text", x = 0, y = -1.5, hjust = 0, size = 1.95, family = FF,
             colour = "grey30", lineheight = 1.05,
             label = "black = CONFIRMED against GEO SOFT + linked publication \u00b7 red bold = CORRECTED in the  provenance audit \u00b7 grey = field outside audit scope\nsmall line under cohort name: \u25cf confirmed / \u25b2 corrected fields per audited cohort (ledger totals: 27 confirmed, 10 corrected, 3 unverified)\nn convention \u2014 GSE186143: 61 in series / 53 locally labelled / 47 joint-evaluable \u00b7 GSE253720: 14 donors (9 CPI-colitis / 3 HC / 2 UC)\nGSE166181: 20 (11 R / 9 NR) \u00b7 GSE216329: 24 (15 irAE / 9) \u00b7 GSE189125: 13 (6 severe / 7 no) \u00b7 m. = metastatic") +
    coord_cartesian(xlim = c(-1, 126), ylim = c(-2.8, ytop + 1.3), clip = "off") +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL,
         title = "Cohort provenance / endpoint ledger (nine endpoint-bearing cohorts, )") +
    theme_void(base_family = FF) +
    theme(plot.title = element_text(size = 9, face = "bold")) +
    lab("a")
  p
}

# =============================================================================
# ED8 -- toxicity stress-test forest (0.5 reference line)
# =============================================================================
make_ed8 <- function() {

  cd4tem_auc <- rc$GSE189125$per_closure_aucs[[1]]$auc   # closure a (PBMC)
  add_sd("ed8", "forest", "data/tox_meta.json",
         "corroborating_not_pooled.GSE189125_CD4Tem_univariate.auc",
         tox60$output$corroborating_not_pooled$GSE189125_CD4Tem_univariate$auc)
  add_sd("ed8", "forest", "data/tox_meta.json",
         "corroborating_not_pooled.GSE189125_CD4Tem_univariate.P",
         tox60$output$corroborating_not_pooled$GSE189125_CD4Tem_univariate$P)
  add_sd("ed8", "forest", "data/joint_stress_test.json",
         "toxicity_cohort_recheck.GSE189125.per_closure_aucs[0].auc (closure a, PBMC)", cd4tem_auc)
  add_sd("ed8", "forest", "data/joint_stress_test.json",
         "toxicity_cohort_recheck.GSE189125.original_keys.verdict.n_significant_p_mannwhitney",
         rc$GSE189125$original_keys$verdict.n_significant_p_mannwhitney)
  add_sd("ed8", "forest", "data/joint_stress_test.json",
         "toxicity_cohort_recheck.GSE189125.original_keys.verdict.n_closures",
         rc$GSE189125$original_keys$verdict.n_closures)
  add_sd("ed8", "forest", "data/tox_meta.json", "per_cohort[1] GSE189125 CD4-T S0",
         list(auc = pc2$auc, ci = pc2$ci, n = pc2$n, n_event = pc2$n_event, n_nonevent = pc2$n_nonevent))
  add_sd("ed8", "forest", "data/tox_meta.json", "per_cohort[0] GSE216329",
         list(auc = pc1$auc, ci = pc1$ci, n = pc1$n, n_event = pc1$n_event, n_nonevent = pc1$n_nonevent))
  add_sd("ed8", "forest", "data/joint_stress_test.json",
         "toxicity_cohort_recheck.GSE216329.original_keys__traceability",
         g2)
  add_sd("ed8", "forest", "data/joint_stress_test.json",
         "toxicity_cohort_recheck.GSE253720.original_keys", g3)
  add_sd("ed8", "forest", "data/tox_meta.json", "output.pooled",
         list(pooled_auc = tox60$output$pooled_auc, ci_low = tox60$output$ci_low,
              ci_high = tox60$output$ci_high, k = tox60$output$k))

  d <- data.frame(
    row = c(
      "GSE189125 \u00b7 CD4-Tem fraction (PBMC)",
      "GSE189125 \u00b7 CD4-T proportion S0",
      "GSE216329 \u00b7 S0 12-dim (0.632-corrected)",
      "GSE253720 \u00b7 CPI-colitis vs HC (probe)",
      "Pooled DL (k = 2, descriptive)"
    ),
    est = c(cd4tem_auc, pc2$auc, pc1$auc,
            g3$`auc_toxicity_CPI_colitis_vs_HC.point`, tox60$output$pooled_auc),
    lo  = c(NA, pc2$ci[[1]], pc1$ci[[1]],
            g3$`auc_toxicity_CPI_colitis_vs_HC.ci95`[[1]], tox60$output$ci_low),
    hi  = c(NA, pc2$ci[[2]], pc1$ci[[2]],
            g3$`auc_toxicity_CPI_colitis_vs_HC.ci95`[[2]], tox60$output$ci_high),
    n_lab = c("n = 13 (6 severe / 7 no)",
              "n = 13 (6 / 7)",
              "n = 24 (15 irAE / 9)",
              "n = 12 (9 CPI / 3 HC)",
              "k = 2 cohorts"),
    pooled = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  d$y <- rev(seq_len(nrow(d)))
  d$nudge <- c(0, 0.18, 0, 0, 0)   # row 2: lift n-label above its wide CI

  warn <- c(
    "no CI in frozen artifacts;\nauthors' Suppl. Table 3 value;\nsignificant in only 2/5 denominator closures",
    "wide CI spans nearly the\nfull AUC range at n = 13",
    "ANTICONSERVATIVE: 0.632-corrected 0.804\nEXCEEDS apparent 0.667; endpoint = any irAE\ngrade 2-3 (organ-specific), NOT grade \u22653 severe",
    "case-control confound (HC comparator\nwithout ICI exposure); biopsies at ACTIVE\nsteroid-experienced/refractory colitis \u2014\nnot baseline; upper CI unclipped (>1)",
    "inherits the anticonservative\n0.804 input; descriptive only"
  )

  p <- ggplot(d, aes(x = est, y = y)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_pointrange(data = subset(d, !is.na(lo)),
                    aes(xmin = lo, xmax = hi, colour = pooled, size = pooled),
                    fatten = 2.4, stroke = 0.6, shape = 19) +
    geom_point(data = subset(d, is.na(lo)), shape = 23, size = 2.8,
               fill = RED, colour = "black", stroke = 0.5) +
    geom_text(aes(label = sprintf("%.3f", est)), vjust = -1.05, size = 2.1,
              family = FF, colour = "grey20") +
    geom_text(aes(label = n_lab, y = y + nudge), x = 0.2, hjust = 0, size = 1.9,
              family = FF, colour = "grey30") +
    geom_text(data = data.frame(y = d$y, warn = warn),
              aes(x = 1.34, y = y, label = warn), inherit.aes = FALSE,
              hjust = 0, size = 1.8, family = FF, colour = "grey25",
              lineheight = 0.92) +
    scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = PURPLE), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.75), guide = "none") +
    scale_x_continuous(breaks = seq(0.2, 1.2, 0.2)) +
    scale_y_continuous(breaks = d$y, labels = d$row) +
    coord_cartesian(xlim = c(0.2, 1.32), clip = "off") +
    labs(x = "AUC (toxicity)", y = NULL,
         title = "Toxicity stress-test cohorts: every headline AUC carries a design caveat",
         caption = "diamond = point estimate without a frozen CI; pooled row = DerSimonian-Laird\nover the two pooled-eligible cohorts (descriptive only, k = 2).\nSources: data/tox_meta.json; data/joint_stress_test.json\n(toxicity_cohort_recheck); data/provenance_ledger.json (endpoint definitions).") +
    theme(axis.text.y = element_text(size = 6.6),
          plot.title = element_text(size = 8.5, face = "bold"),
          plot.caption = element_text(size = 5.8, colour = "grey30"),
          plot.margin = margin(0.4, 4.6, 0.4, 0.4, "cm")) +
    annotate("text", x = -Inf, y = Inf, label = "a", hjust = -0.35, vjust = 0.55,
             size = 3.4, fontface = "bold", family = FF)
  p
}

# =============================================================================
# ED9 -- fold occupancy / identifiability detail
# =============================================================================
make_ed9 <- function() {

  # ---- (a) per-fold class occupancy heatmap (25 folds) ----
  occ <- do.call(rbind, lapply(seq_along(ed9$fold_occupancy_by_repeat), function(r) {
    do.call(rbind, lapply(ed9$fold_occupancy_by_repeat[[r]], function(f) {
      data.frame(repeat_i = r, fold = f$fold,
                 `DCB & severe` = f$DCB_severe, `DCB & no severe` = f$DCB_no_severe,
                 `NDB & no severe` = f$NDB_no_severe, `NDB & severe` = f$NDB_severe,
                 check.names = FALSE)
    }))
  }))
  occ$fold_id <- seq_len(nrow(occ))
  stopifnot(nrow(occ) == 25L, all(occ$fold_sizes == ed9$fold_sizes_all ||
                                    identical(unlist(ed9$fold_sizes_all), rep(c(10L, 10L, 9L, 9L, 9L), 5))))
  long_occ <- do.call(rbind, lapply(c("DCB & severe", "DCB & no severe",
                                      "NDB & no severe", "NDB & severe"), function(cl) {
    data.frame(fold_id = occ$fold_id, repeat_i = occ$repeat_i,
               class = cl, count = occ[[cl]])
  }))
  long_occ$class <- factor(long_occ$class,
                           levels = rev(c("DCB & severe", "DCB & no severe",
                                          "NDB & no severe", "NDB & severe")))
  for (i in seq_len(nrow(occ))) {
    add_sd("ed9", "a", "data/transport_summary_joint.json",
           sprintf("p0_11.ed9.fold_occupancy_by_repeat[%d][%d]", occ$repeat_i[i] - 1L, occ$fold[i] - 1L),
           list(n_test = occ$`DCB & severe`[i] + occ$`DCB & no severe`[i] +
                  occ$`NDB & no severe`[i] + occ$`NDB & severe`[i],
                DCB_severe = occ$`DCB & severe`[i], DCB_no_severe = occ$`DCB & no severe`[i],
                NDB_no_severe = occ$`NDB & no severe`[i], NDB_severe = occ$`NDB & severe`[i]))
  }

  pa <- ggplot(long_occ, aes(x = fold_id, y = class, fill = count)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = count), size = 2.2, family = FF,
              colour = ifelse(long_occ$count >= 3, "white", "grey15")) +
    geom_vline(xintercept = c(5.5, 10.5, 15.5, 20.5), colour = "grey30",
               linetype = 3, linewidth = 0.3) +
    scale_fill_gradient(low = "#eef3f9", high = BLUE, name = "patients", breaks = 1:4) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25)) +
    labs(x = "CV fold (5 repeats, separated by dotted lines)", y = NULL,
         title = "Per-fold class occupancy: identical stratified folds for both model families",
         subtitle = sprintf("fold size pattern %s per repeat; smallest class held out in every fold = %d patient",
                            paste(ed9$fold_sizes_pattern, collapse = "/"),
                            ed9$min_class_occupancy_per_fold)) +
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 6.6),
          axis.text.y = element_text(size = 6.5)) +
    lab("a")

  # ---- (b) per-cell Clopper-Pearson 95% CI ----
  cp <- do.call(rbind, lapply(c("DCB/severe", "DCB/no_severe",
                                "NDB/no_severe", "NDB/severe"), function(k) {
    x <- ed9$per_cell_clopper_pearson_95[[k]]
    add_sd("ed9", "b", "data/transport_summary_joint.json",
           sprintf("p0_11.ed9.per_cell_clopper_pearson_95.%s", k),
           list(n = x$n, proportion = x$proportion, cp95 = x$cp95,
                full_width = x$full_width))
    data.frame(cell = gsub("/", " & ", k), prop = x$proportion,
               lo = x$cp95[[1]], hi = x$cp95[[2]], n = x$n)
  }))
  cp$cell <- factor(rev(paste0(cp$cell, "\n(n = ", cp$n, ")")),
                    levels = rev(paste0(cp$cell, "\n(n = ", cp$n, ")")))
  pb <- ggplot(cp, aes(x = prop, y = cell)) +
    geom_vline(xintercept = 1 / 4, linetype = 3, colour = "grey65", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), fatten = 2.6, linewidth = 0.5,
                    colour = BLUE, shape = 19) +
    annotate("text", x = 0.25, y = 4.55, size = 1.8, family = FF, colour = "grey45",
             label = "1/4 (uniform)") +
    scale_x_continuous(limits = c(0, 0.62), breaks = seq(0, 0.6, 0.1),
                       labels = function(x) sprintf("%.0f%%", 100 * x)) +
    scale_y_discrete(expand = expansion(add = c(0.55, 0.75))) +
    labs(x = "cell proportion (of n = 47 joint-evaluable), Clopper\u2013Pearson 95% CI",
         y = NULL, title = "Every cell proportion is weakly identified") +
    theme(plot.title = element_text(size = 8),
          axis.text.y = element_text(size = 6.3)) +
    lab("b")

  # ---- (c) identifiability ledger ----
  dC <- data.frame(metric = c("joint effective params", "independent effective params",
                              "events per parameter", "smallest cell n",
                              "smallest training cell n"),
                   value = c(ed9$effective_params$joint_multinomial,
                             ed9$effective_params$independent_pair,
                             ed9$events_per_parameter,
                             ed9$cells$`NDB/severe`,
                             ed9$min_training_cell_across_folds))
  dC$y <- rev(seq_len(nrow(dC)))
  dC$lab <- ifelse(dC$value == round(dC$value),
                   sprintf("%.0f", dC$value), sprintf("%.2f", dC$value))
  add_sd("ed9", "c", "data/transport_summary_joint.json",
         "p0_11.ed9.identifiability_ledger",
         list(effective_params_joint = ed9$effective_params$joint_multinomial,
              effective_params_independent = ed9$effective_params$independent_pair,
              events_per_parameter = ed9$events_per_parameter,
              smallest_cell_n = ed9$cells$`NDB/severe`,
              min_training_cell = ed9$min_training_cell_across_folds,
              sign_permutation_p = ed9$sign_permutation_p_on_frozen_deltas))
  pc <- ggplot(dC, aes(x = value, y = y)) +
    geom_col(aes(fill = metric), width = 0.58, colour = "black", linewidth = 0.3) +
    geom_vline(xintercept = 10, linetype = 2, colour = RED, linewidth = 0.35) +
    geom_text(aes(label = lab), hjust = -0.35, size = 2.3, family = FF,
              fontface = "bold") +
    annotate("text", x = 10, y = 0.42, size = 1.9, family = FF, colour = RED,
             hjust = 0.5, label = "conventional 10 events/parameter") +
    scale_fill_manual(values = c("joint effective params" = PURPLE,
                                 "independent effective params" = "grey65",
                                 "events per parameter" = OLIVE,
                                 "smallest cell n" = "#c0504d",
                                 "smallest training cell n" = "#e2a09e"),
                      guide = "none") +
    scale_y_continuous(breaks = 5:1, labels = dC$metric) +
    coord_cartesian(xlim = c(0, 12.5)) +
    labs(x = "value", y = NULL,
         title = "Identifiability ledger at n = 47") +
    theme(plot.title = element_text(size = 8),
          axis.text.y = element_text(size = 6.5)) +
    lab("c")

  pa / (pb | pc) + plot_layout(heights = c(1.15, 1))
}

# =============================================================================
# ED10 -- analysis-freeze & reproducibility map
# =============================================================================
make_ed10 <- function() {

  # enumerate the actual frozen artifacts at run time (real manifest, not typed)
  _jsons <- basename(list.files(J("data", ""), pattern = "\\.json$"))
  _jsons <- basename(list.files(J("data", ""), pattern = "\\.json$"))
  _jsons <- basename(list.files(J("data", ""), pattern = "\\.json$"))
  _jsons <- basename(list.files(J("data", ""), pattern = "\\.json$"))
  _pys   <- basename(list.files(J("scripts"), pattern = "^_.*\\.py$"))
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "data__json", _jsons)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "data__json", _jsons)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "data__json", _jsons)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "data__json", _jsons)
  add_sd("ed10", "manifest", "filesystem enumeration at run time",
         "scripts__py", _pys)
  add_sd("ed10", "manifest", "data/transport_summary_joint.json",
         "locked_seed", SEED)
  add_sd("ed10", "manifest", "data/transport_summary_joint.json",
         "bootstrap_outer_reps", n_boot_outer)
  add_sd("ed10", "manifest", "data/transport_summary_joint.json",
         "bootstrap_inner_OOB_refits", m_boot_inner)
  add_sd("ed10", "manifest", "data/transport_summary_joint.json",
         "sign_permutation_note", sign_perm_note)

  draw_file_list <- function(x, y_top, line_h = 3.1, size = 1.85,
                             prefix = "", col = "grey15") {
    txts <- paste0(prefix, x)
    n <- length(txts)
    annotate("text", x = x_left, y = y_top - (seq_len(n) - 1) * line_h,
             hjust = 0, vjust = 0.75, size = size, family = FF, colour = col,
             label = txts)
  }

  p <- ggplot()
  x_left <- 3.5
  # ---- layer 1: frozen inputs ----
  p <- p +
    annotate("rect", xmin = 1.5, xmax = 33.5, ymin = 30, ymax = 95,
             fill = "#eef3f9", colour = BLUE, linewidth = 0.45) +
    annotate("text", x = x_left, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = BLUE, label = "FROZEN INPUTS (data/)") +
    annotate("text", x = x_left, y = 88.6, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25", label = "data/") +
    draw_file_list(_jsons, 85.5, prefix = "  ") +
    annotate("text", x = x_left, y = 85.5 - 3.1 * length(_jsons) - 1.2,
             hjust = 0, size = 1.9, family = FF, fontface = "bold",
             colour = "grey25", label = "data/  \u00b7  /  \u00b7  /  (prior freezes)") +
    annotate("text", x = x_left, y = 85.5 - 3.1 * length(_jsons) - 4.2,
             hjust = 0, size = 1.8, family = FF, colour = "grey30", lineheight = 1.0,
             label = paste0("  : ", length(_jsons), " \u00b7 : ", length(_jsons),
                            " \u00b7 : ", length(_jsons), " JSON artifacts\n",
                            "  plus data/cd4tem_closure.json,\n",
                            "  data/tox_meta.json, data/disease_coverage.json"))

  # ---- layer 2: scripts ----
  p <- p +
    annotate("rect", xmin = 36.5, xmax = 63.5, ymin = 30, ymax = 95,
             fill = "#f4eef8", colour = PURPLE, linewidth = 0.45) +
    annotate("text", x = 38.5, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = PURPLE, label = "SCRIPTS") +
    annotate("text", x = 38.5, y = 88.6, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25", label = "scripts/ ( analyses)") +
    annotate("text", x = 38.5, y = 85.5 - (seq_along(_pys) - 1) * 3.1,
             hjust = 0, size = 1.85, family = FF, colour = "grey15", label = _pys) +
    annotate("text", x = 38.5, y = 85.5 - 3.1 * length(_pys) - 1.2,
             hjust = 0, size = 1.9, family = FF, fontface = "bold",
             colour = "grey25", label = "figure generation") +
    annotate("text", x = 38.5, y = 85.5 - 3.1 * length(_pys) - 4.2,
             hjust = 0, size = 1.85, family = FF, colour = "grey15", lineheight = 1.05,
             label = "  scripts/figures/make_figures.R ( set)\n  figures/make_stable_figures.R\n  (this script; stable set)")

  # ---- layer 3: outputs ----
  p <- p +
    annotate("rect", xmin = 66.5, xmax = 98.5, ymin = 30, ymax = 95,
             fill = "#ebf1dd", colour = OLIVE, linewidth = 0.45) +
    annotate("text", x = 68.5, y = 92.5, hjust = 0, size = 2.5, family = FF,
             fontface = "bold", colour = OLIVE, label = "OUTPUTS (figures/)") +
    annotate("text", x = 68.5, y = 88.8, hjust = 0, size = 1.9, family = FF,
             colour = "grey25", label = "one figure = PNG (300 dpi) + PDF + SVG") +
    annotate("text", x = 68.5, y = 84.5, hjust = 0, vjust = 1, size = 1.85,
             family = FF, colour = "grey15", lineheight = 1.05,
             label = "fig1_transport_system\nfig6_joint_boundary\ned2_cohort_ledger\ned8_toxicity_stress\ned9_joint_identifiability\ned10_analysis_freeze") +
    annotate("text", x = 68.5, y = 64.0, hjust = 0, size = 1.9, family = FF,
             fontface = "bold", colour = "grey25", label = "source_data/") +
    annotate("text", x = 68.5, y = 60.6, hjust = 0, vjust = 1, size = 1.85,
             family = FF, colour = "grey30", lineheight = 1.05,
             label = "one JSON per figure; every panel value\nkeyed to its frozen source file + key\n(machine-readable figure archive)")

  # arrows between layers
  p <- p +
    annotate("segment", x = 34, xend = 36.2, y = 62, yend = 62,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8,
             colour = "grey35") +
    annotate("segment", x = 64, xend = 66.2, y = 62, yend = 62,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8,
             colour = "grey35")

  # ---- locked parameters band ----
  p <- p +
    annotate("rect", xmin = 1.5, xmax = 98.5, ymin = 4, ymax = 24,
             fill = "grey97", colour = "grey35", linewidth = 0.4) +
    annotate("text", x = 3.5, y = 20.2, hjust = 0, size = 2.4, family = FF,
             fontface = "bold", label = "Locked analysis parameters ( freeze)") +
    annotate("text", x = 3.5, y = 15.0, hjust = 0, size = 2.05, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = sprintf("global seed = %d (all  resampling)\nbootstrap: %d outer reps \u00d7 %s\njoint stress test: 5-fold \u00d7 5 repeats, stratified on the 4 (E,T) cells",
                             SEED, n_boot_outer, m_boot_inner)) +
    annotate("text", x = 52, y = 15.0, hjust = 0, size = 2.05, family = FF,
             colour = "grey15", lineheight = 1.05,
             label = sprintf("sign-permutation test: %s\nridge \u03bb = 1/394 frozen project convention\nno number in any figure is typed: all read from frozen JSON keys",
                             sign_perm_note)) +
    annotate("text", x = 50, y = 6.8, size = 2.0, family = FF, colour = "grey30",
             label = "any change to a frozen input requires a new version directory, not an edit")

  p + theme_void(base_family = FF) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    lab("a")
}

# =============================================================================
# build all
# =============================================================================
FIGS <- list(
  list(sdkey = "fig1",  name = "fig1_transport_system",      fn = make_fig1,  w = 180, h = 215),
  list(sdkey = "fig6",  name = "fig6_joint_boundary",        fn = make_fig6,  w = 180, h = 175),
  list(sdkey = "ed2",   name = "ed2_cohort_ledger",          fn = make_ed2,   w = 180, h = 130),
  list(sdkey = "ed8",   name = "ed8_toxicity_stress",        fn = make_ed8,   w = 180, h = 115),
  list(sdkey = "ed9",   name = "ed9_joint_identifiability",  fn = make_ed9,   w = 180, h = 185),
  list(sdkey = "ed10",  name = "ed10_analysis_freeze",       fn = make_ed10,  w = 180, h = 135)
)

PANELS <- list(
  fig1_transport_system     = c("a system chain (no data)", "b primary estimand", "c cohort map (9 cohorts)"),
  fig6_joint_boundary       = c("a four-quadrant occupancy + CP 95% CI", "b per-fold log-loss (25 folds)",
                                "c identifiability", "d future paired-bridge design"),
  ed2_cohort_ledger         = c("cohort x field ledger table (9 cohorts, verdict-coloured)"),
  ed8_toxicity_stress       = c("toxicity stress-test forest (5 rows, 0.5 reference)"),
  ed9_joint_identifiability = c("a fold occupancy heatmap (25 folds x 4 classes)",
                                "b per-cell Clopper-Pearson CI", "c identifiability ledger"),
  ed10_analysis_freeze      = c("frozen inputs -> scripts -> outputs manifest + locked parameters")
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
# verification: 5 figures x 3 formats non-zero + 5 source-data JSONs parseable
# =============================================================================
cat("\n--- verification ---\n")
n_ok <- 0; n_expected <- 0
for (f in FIGS) {
  for (ext in c("pdf", "png", "svg")) {
    n_expected <- n_expected + 1
    fp <- file.path(OUTDIR, paste0(f$name, ".", ext))
    sz <- if (file.exists(fp)) file.size(fp) else -1
    ok <- sz > 0
    n_ok <- n_ok + as.integer(ok)
    if (!ok) cat(sprintf("[MISSING/EMPTY] %s\n", fp))
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
