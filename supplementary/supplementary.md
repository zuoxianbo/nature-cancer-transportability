# Supplementary Information: Decomposing the transportability of immune-state biomarkers across checkpoint immunotherapy

**Supplementary Notes 1–6 and Tables S1–S8.** All numbers are read programmatically from the frozen JSON archives (`data/*.json`, `data/*.json`; seed 20260911) and are reproducible with `scripts/analysis/*.py`, `scripts/analysis/*.py` (Python 3.13, NumPy/SciPy/scikit-learn) and `scripts/figures/make_anchor_figures.R`, `scripts/figures/make_figures.R`, `scripts/figures/make_stable_figures.R` (R 4.6.1, ggplot2 4.0.3, 300 dpi, theme_classic, Arial). No value in this document is transcribed by hand.

## Supplementary Note 1 | Cohort provenance and endpoint harmonisation

All cohorts are public GEO deposits re-verified field-by-field against GEO SOFT records and the linked primary publications during the  provenance audit. Endpoint roles follow the primary estimand of the main text: IMvigor210 is the single source anchor; GSE165252 and GSE78220 are the external targets; all source-swapped (leave-one-cohort-out, LOCO) refits are sensitivity analyses and never define the primary estimand. GSE287541 (ChILI) is locked as withheld and contributes to no estimate; the single-cell toxicity cohorts (GSE189125, GSE216329, GSE253720) are used only in the stress tests of Extended Data Fig. 8 with the design caveats recorded there.

**Table S1. Cohort ledger.** n and event counts are the analysis-locked labelled sets (`rerun_benchmark.json:cohorts`; casemix block of `clinical_only.json`).

| Cohort | GSE accession | Assay | n (labelled) | Events | Endpoint | Role |
|---|---|---|---|---|---|---|
| IMvigor210 | GSE139551 | bulk RNA-seq (tumour) | 298 | 68 | binary durable clinical benefit / response | source anchor (within reference only) |
| GSE165252 | GSE165252 | bulk RNA-seq (tumour) | 71 | 23 | binary durable clinical benefit / response | external target (primary transport) |
| GSE78220 | GSE78220 | bulk RNA-seq (tumour) | 25 | 13 | binary durable clinical benefit / response | external target (primary transport) |
| GSE186143 | GSE186143 | scRNA-seq (PBMC, bulk TPM) | 47 joint-evaluable (53 labelled of 61 in series) | 27 DCB / 20 NDB | efficacy × severe irAE (4-cell) | joint identifiability analysis |
| GSE189125 | GSE189125 | scRNA-seq (PBMC) | 13 | 6 severe / 7 no | severe irAE | toxicity stress test |
| GSE216329 | GSE216329 | scRNA-seq 10x (T cells) | 24 | 15 irAE / 9 no | any irAE grade 2–3 (organ-specific) | toxicity stress test |
| GSE253720 | GSE253720 | scRNA-seq (colon biopsy) | 14 | 9 CPI / 3 HC | CPI-colitis (case-control) | toxicity stress test (probe) |
| GSE287541 | GSE287541 | scRNA-seq | 28 (61 in series) | — | irAE (longitudinal) | withheld; not used in any estimate |

GSE186143 contributes 47 joint-evaluable patients (response × severe-irAE 4-cell design; 53 labelled patients in the series). Toxicity-cohort endpoint definitions are NOT harmonisable across cohorts (severe grade ≥3 vs any grade 2–3 organ-specific vs case-control colitis) and are labelled non-equivalent rather than pooled.

## Supplementary Note 2 | Statistical estimands and transport criteria

**Primary estimand.** The primary transport estimand is the held-out performance of the frozen IMvigor210-anchored system (S0: library-size normalisation → CLR transform → 322-gene panel → per-cohort uncentred PCA (k = 2) → Procrustes alignment with IMvigor210 as reference (rotation + translation only) → ridge logistic fit in the anchor), applied transform-only to the target cohort. No target label enters any fitting step. Source-swapped LOCO refits answer a different question (does any two-cohort mixture transfer?) and are reported as sensitivity analyses only.

**Licensing convention.** A transport claim requires both (C1) the external AUC bootstrap 95% CI (BCa, 2,000 resamples, seed 20260911) to exclude 0.5 in the beneficial direction, and (C2) |calibration slope − 1| ≤ 0.10 on the logit scale. The ±0.10 tolerance is a declared reporting convention, not an estimated quantity.

**Threshold sensitivity (P0-3).** Re-running the complete adjudication under AUC-only, calibration-only and joint criteria, with calibration-slope tolerances of ±0.05, ±0.10 and ±0.20 and under both slope conventions (probability-scale and logit-scale), changes no verdict: flip count = 0 (probability-scale convention) and 0 (logit-scale convention); transport permitted at any tolerance = false. The central conclusion does not depend on the 0.10 tolerance.

**Continuous metrics.** AUC, calibration intercept and slope, Brier score and log-loss are reported alongside every licensing decision (Tables S3, S4, S6), because a dual threshold can only summarise, not replace, the continuous transport behaviour.

## Supplementary Note 3 | Representation audit (21-configuration grid)

The controlled perturbation grid spans five layers: the gene set (baseline 322-gene panel; data-driven top-100 / top-50 subsets; five random 100-gene and five random 50-gene draws), the normalisation (log2-CPM-raw; within-cohort z-score), the PCA dimension (k = 1, 3, 5) and the ridge penalty (λ = 0.01, 0.1, 1 × 1/n). Five former “measured-only” configurations were retracted: under the recovered gene identity the anchor panel has zero all-zero genes, so every measured-only subset coincides with the full panel (`rerun_ablation.json:removed_configs`, 5 entries). The 21 non-degenerate configurations are the final grid.

**Table S5. 21-configuration representation audit** (external AUCs of the frozen primary transport; BCa 95% CI; `rerun_ablation.json:per_config`).

| Configuration | Layer | n genes | k | AUC GSE165252 [95% CI] | AUC GSE78220 [95% CI] |
|---|---|---|---|---|---|
| baseline | baseline | 320 | 2 | 0.504 [0.354, 0.658] | 0.660 [0.404, 0.862] |
| gene_top100 | gene_set | 100 | 2 | 0.535 [0.387, 0.677] | 0.513 [0.260, 0.758] |
| gene_top50 | gene_set | 50 | 2 | 0.466 [0.330, 0.611] | 0.596 [0.327, 0.800] |
| gene_rand100_s0 | gene_set | 100 | 2 | 0.538 [0.387, 0.689] | 0.526 [0.276, 0.750] |
| gene_rand100_s1 | gene_set | 100 | 2 | 0.573 [0.413, 0.725] | 0.571 [0.331, 0.792] |
| gene_rand100_s2 | gene_set | 100 | 2 | 0.486 [0.339, 0.642] | 0.615 [0.364, 0.827] |
| gene_rand100_s3 | gene_set | 100 | 2 | 0.549 [0.402, 0.697] | 0.590 [0.340, 0.806] |
| gene_rand100_s4 | gene_set | 100 | 2 | 0.481 [0.332, 0.637] | 0.615 [0.364, 0.827] |
| gene_rand50_s0 | gene_set | 50 | 2 | 0.467 [0.326, 0.606] | 0.647 [0.390, 0.844] |
| gene_rand50_s1 | gene_set | 50 | 2 | 0.495 [0.351, 0.648] | 0.622 [0.332, 0.833] |
| gene_rand50_s2 | gene_set | 50 | 2 | 0.418 [0.273, 0.570] | 0.654 [0.388, 0.851] |
| gene_rand50_s3 | gene_set | 50 | 2 | 0.602 [0.456, 0.732] | 0.513 [0.254, 0.747] |
| gene_rand50_s4 | gene_set | 50 | 2 | 0.565 [0.416, 0.700] | 0.538 [0.301, 0.792] |
| norm_log2raw | normalization | 320 | 2 | 0.476 [0.328, 0.626] | 0.654 [0.407, 0.853] |
| norm_within_z | normalization | 320 | 2 | 0.509 [0.353, 0.661] | 0.654 [0.390, 0.853] |
| pc_k1 | pc_k | 320 | 1 | 0.500 [0.356, 0.644] | 0.622 [0.360, 0.846] |
| pc_k3 | pc_k | 320 | 3 | 0.485 [0.337, 0.635] | 0.679 [0.417, 0.872] |
| pc_k5 | pc_k | 320 | 5 | 0.464 [0.328, 0.611] | 0.583 [0.330, 0.798] |
| ridge_lam0.01 | ridge_lambda | 320 | 2 | 0.504 [0.354, 0.658] | 0.660 [0.404, 0.862] |
| ridge_lam0.1 | ridge_lambda | 320 | 2 | 0.504 [0.354, 0.658] | 0.660 [0.404, 0.862] |
| ridge_lam1 | ridge_lambda | 320 | 2 | 0.504 [0.354, 0.658] | 0.660 [0.404, 0.862] |

No configuration establishes transport in either target cohort (any = false); the observed maximum external AUC over the grid is 0.679 and is not distinguishable from the selection-adjusted permutation null (Table S8b).

## Supplementary Note 4 | Gene-row identity recovery

The anonymised IMvigor210 deposit preserved row ORDER but not row labels: the expression matrix carried placeholder identifiers (gene_1 … gene_31085). All – pipeline versions had assumed exmat row i == fData row i under an Entrez string sort, and one Supplementary Table column was misparsed as Entrez identifiers, producing a pseudo-random 228-gene intersection with 70 all-zero genes carrying 88.7% of PC1 loading mass. The recovered key (exmat row j == the j-th fData row under upper-cased symbol alphabetical sort) was adjudicated among candidate sorts by biologically impossible-marker gates (Y-linked / XIST / housekeeping), mirror numeric identity and cross-cohort mean-expression rank correlation; all verification gates pass (`imvigor_gene_recovery.json:verification_gates`, OVERALL PASS). The corrected row key is released as `exmat_row_key.csv` (31,085 rows, SHA-256 recorded in the provenance ledger).

**Table S2. Gene-row key and panel accounting** (`rerun_benchmark.json:anchor`; `imvigor_gene_recovery.json`).

| Quantity | Value |
|---|---|
| Gene-set definition size | 322 |
| Measurable in IMvigor210 (recovered key) | 320 / 322 (missing: IGHG1, PRKN) |
| Measurable in GSE165252 | 322 |
| Measurable in GSE78220 | 320 |
| Cross-cohort projection panel (common) | 320 |
| All-zero genes across the 348-gene union | 0 |
| Within-cohort 0.632 AUC, anchor (old → new) | 0.567 → 0.657 |
| External AUC → GSE78220 (old → new) | 0.410 → 0.660 |
| Two-criterion LOCO verdict (old → new) | refused → refused (unchanged) |

The former 70-zero-gene pattern was a data-mapping error, not biological evidence; all anchor-arm numbers were re-run under the recovered key before submission (P0-1).

## Supplementary Note 5 | Model competition and positive controls

Eight strategies were compared under the identical LOCO design (identical held-out patients): M0 prevalence-only intercept; M1 T-cell-receptor-only; M2 frozen S0; M3 supervised recalibration of M2 (100 stratified 50/50 repeats); M4 correlation-aligned (CORAL) variant; M5a/M5b pooled with platform- and disease-intercept shifts; M6 ridge recalibration. Because a single score cannot summarise a transport claim, all four probability metrics are reported (Table S6) and no method is ranked as a winner.

**Table S6. Model competition under the LOCO sensitivity estimand** (`rerun_benchmark.json:loco_benchmark` + `:adjudication`). Verdict = two-criterion rule (C1 ∧ C2); IMvigor210 as target is served by the pooled other two cohorts.

| Model | Target | AUC [95% CI] | Slope (logit) | Brier | Log-loss | Verdict |
|---|---|---|---|---|---|---|
| M0 | IMvigor210 | 0.500 [0.500, 0.500] | NA (constant) | 0.198 | 0.587 | refused |
| M0 | GSE165252 | 0.500 [0.500, 0.500] | NA (constant) | 0.224 | 0.643 | refused |
| M0 | GSE78220 | 0.500 [0.500, 0.500] | NA (constant) | 0.324 | 0.864 | refused |
| M1 | IMvigor210 | 0.537 [0.463, 0.615] | 0.53 | 0.214 | 0.619 | refused |
| M1 | GSE165252 | 0.500 [0.356, 0.644] | 0.13 | 0.226 | 0.649 | refused |
| M1 | GSE78220 | 0.622 [0.360, 0.846] | 20.73 | 0.325 | 0.865 | refused |
| M2 | IMvigor210 | 0.663 [0.585, 0.731] | 2.93 | 0.196 | 0.583 | refused |
| M2 | GSE165252 | 0.511 [0.361, 0.666] | -0.06 | 0.264 | 0.824 | refused |
| M2 | GSE78220 | 0.615 [0.365, 0.833] | 1.87 | 0.248 | 0.688 | refused |
| M3 | IMvigor210 | 0.665 [0.604, 0.730] | 1.15 | 0.168 | 0.511 | refused |
| M3 | GSE165252 | 0.445 [0.356, 0.522] | -1.61 | 0.229 | 0.654 | refused |
| M3 | GSE78220 | 0.512 [0.154, 0.731] | -1.61 | 0.281 | 0.832 | refused |
| M4 | IMvigor210 | 0.660 [0.581, 0.729] | 3.10 | 0.197 | 0.584 | refused |
| M4 | GSE165252 | 0.511 [0.361, 0.666] | -0.06 | 0.264 | 0.824 | refused |
| M4 | GSE78220 | 0.628 [0.378, 0.840] | 1.88 | 0.242 | 0.675 | refused |
| M5a | IMvigor210 | 0.665 [0.586, 0.733] | 1.69 | 0.188 | 0.562 | refused |
| M5a | GSE165252 | 0.491 [0.340, 0.639] | -0.07 | 0.241 | 0.705 | refused |
| M5a | GSE78220 | 0.660 [0.404, 0.862] | 2.04 | 0.312 | 0.832 | refused |
| M5b | IMvigor210 | 0.667 [0.589, 0.735] | 1.78 | 0.188 | 0.563 | refused |
| M5b | GSE165252 | 0.514 [0.363, 0.665] | -0.05 | 0.243 | 0.715 | refused |
| M5b | GSE78220 | 0.615 [0.365, 0.833] | 1.47 | 0.315 | 0.842 | refused |
| M6 | IMvigor210 | 0.668 [0.602, 0.731] | 1.13 | 0.167 | 0.510 | refused |
| M6 | GSE165252 | 0.445 [0.357, 0.515] | -0.73 | 0.229 | 0.654 | refused |
| M6 | GSE78220 | 0.512 [0.154, 0.731] | -0.75 | 0.281 | 0.834 | refused |

All 24 cells are refused. Seven LOCO cells exceed the AUC-only criterion in the beneficial direction and fail only on calibration (`rerun_benchmark.json:near_miss:auc_passes_but_slope_fails`); the closest overall is M3 → IMvigor210 (AUC 0.665, slope 1.15).

**Positive controls (P0-2).** The published 18-gene T-cell-inflamed GEP (Cristescu et al. 2018) was pushed through the same frozen transport pipeline and the same two-criterion evaluation under the recovered anchor identity; the 6-gene IFN-γ panel (Ayers et al. 2017) is sensitivity only. Results (`data/positive_control_gep18.json`):

| Signature | Cohort | AUC [95% CI] | Slope | Two-criterion verdict |
|---|---|---|---|---|
| 18-gene T-cell-inflamed GEP (primary) | GSE165252 | 0.546 [0.403, 0.697] | 0.52 | fail |
| 18-gene T-cell-inflamed GEP (primary) | GSE78220 | 0.571 [0.318, 0.788] | 0.70 | fail |
| 6-gene IFN-γ (sensitivity) | GSE165252 | 0.561 [0.422, 0.712] | 0.49 | fail |
| 6-gene IFN-γ (sensitivity) | GSE78220 | 0.494 [0.247, 0.721] | 0.25 | fail |

Within the anchor (IMvigor210, Efron 0.632), the 18-gene GEP reaches 0.573 and the 6-gene panel 0.605: known reproducible signatures attain within-cohort signal in these cohorts, yet both fail the two-criterion rule in every external cohort and at every slope tolerance. The evaluation environment is therefore not a trivial rejector; the refusal of the 322-gene representation is not explained by an over-strict rule alone.

## Supplementary Note 6 | Joint efficacy–toxicity analysis (GSE186143)

GSE186143 is the only cohort in this study carrying both an efficacy endpoint (durable clinical benefit) and a severe immune-related toxicity endpoint in the same patients. With 47 joint-evaluable patients the 4-cell occupancy is:

| | severe irAE yes | severe irAE no |
|---|---|---|
| DCB | 10 | 17 |
| NDB | 5 | 15 |

Marginals: 27 DCB / 20 NDB; 15 severe / 32 no severe irAE. The efficacy–toxicity cross-tabulation is not significant (odds ratio 1.76; Fisher's exact P = 0.529; per-cell Clopper–Pearson 95% CIs all span their marginal rates).

A joint multinomial model of the 4-cell outcome was compared with a pair of independent binomial models by 25-fold cross-validation: pooled log-loss 1.345 (joint) vs 1.356 (independent); the joint model wins in 12 of 25 folds, and the fold-level difference is not distinguishable from zero (sign-permutation P = 0.986).

Identifiability assessment: the joint multinomial uses 9 effective parameters against 6 for the independent pair, i.e. 5.2 events per parameter, with a minimum training cell of 4 patients (minimum class occupancy 1 per fold). The cohort cannot identify a joint efficacy–toxicity representation: this is reported as an identifiability boundary and a design requirement for paired prospective studies, not as evidence about a shared mechanism.

## Table S3. Primary and sensitivity AUCs

**S3a. Primary frozen transports (the primary estimand).**

| Train → target | Model | AUC [BCa 95% CI] | n target (events) |
|---|---|---|---|
| IMvigor210 → GSE165252 | M2 (frozen S0) | 0.504 [0.354, 0.658] | 71 (23) |
| IMvigor210 → GSE78220 | M2 (frozen S0) | 0.660 [0.404, 0.862] | 25 (13) |

**S3b. True external pairwise transports, all six ordered directions (sensitivity; single-source frozen, transform-only).**

| Train → target | AUC [95% CI] | n (events) |
|---|---|---|
| IMvigor210 → GSE165252 | 0.504 [0.354, 0.658] | 71 (23) |
| IMvigor210 → GSE78220 | 0.660 [0.404, 0.862] | 25 (13) |
| GSE165252 → IMvigor210 | 0.391 [0.321, 0.468] | 298 (68) |
| GSE165252 → GSE78220 | 0.442 [0.224, 0.699] | 25 (13) |
| GSE78220 → IMvigor210 | 0.640 [0.560, 0.710] | 298 (68) |
| GSE78220 → GSE165252 | 0.515 [0.367, 0.672] | 71 (23) |

**S3c. Descriptive pooled external AUC (k = 3 cohorts).** 0.6146 [0.5147, 0.7146], I² = 0.387, Q P = 0.195. Reported for completeness only; with three cohorts of very different sizes it is not used to license transport (Extended Data Fig. 3).

## Table S4. Calibration and probability-accuracy metrics

**S4a. LOCO M2 (frozen S0) and M3 (100-repeat recalibration mean) per target.**

| Target | Model | Calib. intercept | Calib. slope (logit) | Brier | Log-loss | Mean predicted | Observed rate |
|---|---|---|---|---|---|---|---|
| IMvigor210 | M2 | -0.04 | 2.93 | 0.196 | 0.583 | nan | nan |
| IMvigor210 | M3 | 0.17 | 1.15 | 0.168 | 0.511 | nan | nan |
| GSE165252 | M2 | -0.86 | -0.06 | 0.264 | 0.824 | nan | nan |
| GSE165252 | M3 | -1.95 | -1.61 | 0.229 | 0.654 | nan | nan |
| GSE78220 | M2 | 0.66 | 1.87 | 0.248 | 0.688 | nan | nan |
| GSE78220 | M3 | 0.18 | -1.61 | 0.281 | 0.832 | nan | nan |

**S4b. Supervised recalibration before/after (M2 → M3, logit-scale slope; 100 stratified 50/50 repeats).**

| Target | M2 slope | M3 repeat-mean slope [95% CI] | Slope improved |
|---|---|---|---|
| IMvigor210 | 2.93 | 1.15 [0.46, 2.70] | yes |
| GSE165252 | -0.06 | -1.61 [-6.78, 3.13] | no |
| GSE78220 | 1.87 | -1.61 [-26.06, 11.95] | no |

**S4c. Within-cohort 0.632 references (Efron optimism-corrected, 600 refits).**

| Cohort | AUC (0.632) [95% CI] |
|---|---|
| IMvigor210 | 0.657 [0.586, 0.728] |
| GSE165252 | 0.428 [0.326, 0.530] |
| GSE78220 | 0.622 [0.419, 0.825] |

## Table S7. Case-mix and outcome-shift audit

| Quantity | IMvigor210 | GSE165252 | GSE78220 |
|---|---|---|---|
| n (labelled) | 298 | 71 | 25 |
| Response rate | 0.228 | 0.324 | 0.520 |
| Male sex | 0.782 | n/a | 0.680 |
| Stage III–IV / M1 | 0.188 | n/a | n/a |
| SMD stage (IMvigor210 vs GSE78220) | — | — | -2.50 |
| SMD male sex (IMvigor210 vs GSE78220) | — | — | 0.23 |

Stage constructs are not equivalent between cohorts (IMvigor210 AJCC stage I–IV vs GSE78220 M1 disease); the standardized differences are reported only for the male-sex contrast and as a descriptive stage gap, and no causal case-mix attribution is made. In GSE165252 no covariates are deposited; every clinical model degenerates to the train-prevalence constant 0.228 (Extended Data Fig. 7).

## Table S8. Simulation and power/precision scenarios (P0-4, P0-5)

**S8a. External-validation precision simulation** (binormal equal-variance latent score; 5,000 simulations × 2,000-bootstrap CIs per cell; outcome prevalence fixed at the observed target rates; `precision_simulation.json`). P(C1) = probability that the 95% CI excludes 0.5 in the beneficial direction.

| True AUC | n = 71, 23 events: mean CI width | P(C1) | n = 25, 13 events: mean CI width | P(C1) |
|---|---|---|---|---|
| 0.55 | 0.290 | 9.9% | 0.474 | 5.5% |
| 0.60 | 0.284 | 28.3% | 0.464 | 12.3% |
| 0.65 | 0.275 | 54.1% | 0.446 | 24.1% |
| 0.70 | 0.262 | 79.2% | 0.426 | 38.6% |
| 0.75 | 0.243 | 94.4% | 0.393 | 59.1% |
| 0.80 | 0.219 | 99.3% | 0.353 | 76.9% |
| 0.85 | 0.188 | 100.0% | 0.300 | 91.5% |

At the observed sizes, only large true effects are detectable reliably: for n = 71 (GSE165252) the probability of excluding 0.5 exceeds 90% only for true AUC ≥ 0.75; for n = 25 (GSE78220) even a true AUC of 0.80 yields 77% power. Conversely, the observed primary estimates cannot establish transport, but they also cannot rule out moderate true effects; this is the precision boundary of the present data.

**S8b. Selection-adjusted permutation null for the 21-configuration search (P0-4).** Anchor outcome labels permuted 2,000 times (target labels fixed; the same permuted labelling applied to all 21 configurations within a permutation; `permutation_null.json`). The observed maximum external AUC over the grid is compared with the permutation max-AUC distribution.

| Target cohort | Observed max AUC (argmax config) | Permuted max: mean | q95 | max | One-sided P |
|---|---|---|---|---|---|
| GSE165252 | 0.602 (gene_rand50_s3) | 0.582 | 0.644 | 0.652 | 0.394 |
| GSE78220 | 0.679 (pc_k3) | 0.637 | 0.744 | 0.827 | 0.398 |

The observed configuration-selection maxima lie inside the permutation null distributions: the 21-configuration search does not contain a signal beyond what selection alone produces, and the uniform refusal of the grid is not an artefact of multiple testing against a favourable direction.

---
*Reproducibility: every value above is extracted programmatically from the frozen JSON archives cited in place; SHA-256 checksums of the frozen inputs are recorded in `data/analysis_manifest.json` and the Extended Data Fig. 10 manifest. Seed 20260911 throughout.*
