# Figure 2–6 Number Reconciliation —  manuscript vs  frozen JSON

Frozen data source root: `_20260926_Nature_Cancer/data/{,}/` and `figures/source_data/`.

Every "manuscript value" below is quoted from `__work/NatureCancer_Final_Submission_Package_/_ms.txt` (paragraph index given). "JSON value" is the raw value from the frozen upstream. "Rounded to N dp" is what appears in the figure label or the SI table.

## Fig. 2 — External transport & source-swapped sensitivity

| Item | Manuscript | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| P13 predicted 5th percentile | 0.310 | `predicted_quantiles.0.05 = 0.310302` | 0.310 | PASS |
| P13 predicted 95th percentile | 0.474 | `predicted_quantiles.0.95 = 0.474214` | 0.474 | PASS |
| P13 median predicted probability | 0.394 | `predicted_quantiles.0.5 = 0.393756` | 0.394 | PASS |
| P13 observed event rate | 0.228 | `overall_observed_rate = 0.228188` | 0.228 | PASS |
| **P13 wording** | "held-out predictions were concentrated around the source event rate" | median 0.394 vs observed 0.228, Δ = 0.166 (73% relative) | — | **FAIL** — wording contradicts the numbers in the same sentence |
| P13 calibration slopes (source-swapped / LOCO) | 2.93 / −0.06 / 1.87 | `calib_slope` M2: `2.928521 / −0.059868 / 1.874568` | 2.93 / −0.06 / 1.87 | PASS |
| Panel c S0 transport slope, IMv→GSE165252 | −0.06 (figure label) | `IMv→GSE165252.slope = −0.063323` | −0.06 | PASS |
| Panel c S0 transport slope, IMv→GSE78220 | 1.35 (figure label) | `IMv→GSE78220.slope = 1.349792` | 1.35 | PASS (but note: 1.35 is main-transport slope, 1.87 is source-swapped — two different frameworks; consider footnote to prevent reader confusion) |
| Panel c 6-gene GSE165252 | 0.561 [0.422, 0.712] slope 0.49 | `0.560688 [0.422040, 0.711897] slope 0.489` | 0.561 [0.422, 0.712] slope 0.49 | PASS |
| Panel c 6-gene GSE78220 | 0.494 [0.247, 0.721] slope 0.25 | `0.493590 [0.246667, 0.721437] slope 0.250` | 0.494 [0.247, 0.721] slope 0.25 | PASS |
| Panel c 6-gene within-IMvigor210 (Efron 0.632) | 0.605 | `0.604979` | 0.605 | PASS |
| Panel c 6-gene gene list (figure footnote) | IDO1, CXCL9, CXCL10, CCL5, HLA-DRA, STAT1 | `_rerun_positive_control.genes_used_in_signature` = `["IDO1","CXCL9","CXCL10","CCL5","HLA-DRA","STAT1"]` | — | PASS (matches exactly) |
| Pairwise span (6 directions) | 0.391–0.660 | `min 0.391304, max 0.660256` | 0.391–0.660 | PASS |
| **Fig. 2 PNG in  package** | orphan, md5 `ada64ae0…`, 357,106 B, 2026-10-04 |  upstream md5 `0b97969e…`, 348,112 B, 2026-09-26 | — | **FAIL** —  replaces with  frozen byte-identical |

## Fig. 3 — Calibration pathways & recalibration

| Item | Manuscript | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| IMvigor210 slope M2 | 2.93 | 2.928521 | 2.93 | PASS |
| GSE165252 slope M2 | −0.06 | −0.059868 | −0.06 | PASS |
| GSE78220 slope M2 | 1.87 | 1.874568 | 1.87 | PASS |
| IMvigor210 log-loss before→after recalibration | 0.583 → 0.511 | 0.582832 → 0.511 | 0.583 → 0.511 | PASS |
| IMvigor210 Brier before→after | 0.196 → 0.168 | 0.196256 → 0.168 | 0.196 → 0.168 | PASS |
| GSE78220 log-loss before→after (recalibration deteriorates) | 0.688 → 0.832 | 0.688419 → 0.832 | 0.688 → 0.832 | PASS |
| Pooled AUC 0.615 (label on fig3) | cross-reference | `pooled_auc = 0.615` (in transport summary) | 0.615 | PASS (cross-figure label; not stored in fig3's own source_data JSON) |

## Fig. 4 — 21-config ablation + permutation null

| Item | Manuscript | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| 21-config count | 21 | `len(grid) = 21` | 21 | PASS |
| External AUC range, GSE165252 | 0.418–0.602 | `min 0.417572, max 0.602355` | 0.418–0.602 | PASS |
| External AUC range, GSE78220 | 0.513–0.679 | `min 0.512821, max 0.679487` | 0.513–0.679 | PASS |
| Variance decomp gene_set | 0.0028 | 0.002752 | 0.0028 | PASS |
| Variance decomp normalization | 0.0002 | 0.000169 | 0.0002 | PASS |
| Variance decomp pc_k | 0.0011 | 0.001074 | 0.0011 | PASS |
| Variance decomp ridge | 0 | 0 | 0 | PASS |
| **Permutation null, GSE78220 observed max AUC** | **0.680** | **0.6794871794871795 → 0.679** | **0.679** | **FAIL** — manuscript writes 0.680, JSON rounds to 0.679 |
| Permutation null, GSE78220 perm q95 | 0.744 | 0.74359 | 0.744 | PASS |
| Permutation null, GSE78220 P | 0.398 | 0.397801 | 0.398 | PASS |
| Permutation null, GSE165252 observed | 0.602 | 0.602355 | 0.602 | PASS |
| Permutation null, GSE165252 perm q95 | 0.644 | 0.644022 | 0.644 | PASS |
| Permutation null, GSE165252 P | 0.394 | 0.393803 | 0.394 | PASS |
| **Panel reference in manuscript P16** | "(Fig. 4c)" | Caption: (c) gene-row identity; **(d) permutation null** | — | **FAIL** — should be (Fig. 4d) |

## Fig. 5 — Benchmark strategies & near-miss controls

| Item | Manuscript | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| 24-cell adjudication permitted | 0/24 | `verdict_counts.n_permitted = 0, n_total = 24` | 0/24 | PASS |
| Near-miss (GSE78220→IMvigor210) AUC | 0.640 [0.560, 0.710] slope 0.43 | `0.640153 [0.559741, 0.709772] slope 0.425973` | 0.640 [0.560, 0.710] slope 0.43 | PASS |
| Clinical-only GSE78220 AUC | 0.487 [0.305, 0.698] | `0.487179 [0.305195, 0.698155]` | 0.487 [0.305, 0.698] | PASS |
| Clinical-only GSE165252 | constant 0.228 | `constant_prediction = 0.228188 = train_prevalence_reference` | 0.228 | PASS |
| **TCR diversity M1 (GSE186143)** | AUC 0.570, log-loss 0.696, Brier 0.250 | **from `joint_stress_test.json` (NOT )**: `0.570370 / 0.695707 / 0.249499` | 0.570 / 0.696 / 0.250 | PASS (transcription) / **GAP (source is , not )** |
| **TCR vs S0+features M4** | AUC 0.544, log-loss 0.767, Brier 0.271 | ****: `0.544444 / 0.766699 / 0.271055` | 0.544 / 0.767 / 0.271 | PASS (transcription) / **GAP (source is , not )** |

## Fig. 6 — Joint boundary & toxicity stress

| Item | Manuscript | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| n joint-evaluable | 47 | 47 | 47 | PASS |
| cells DSB/severe, DSB/no severe, NSB/no severe, NSB/severe | 10, 17, 15, 5 | 10, 17, 15, 5 | same | PASS |
| OR E×T | 1.76 | 1.764706 | 1.76 | PASS |
| Fisher exact P | 0.53 | 0.529332 | 0.53 | PASS |
| Pooled log-loss joint vs independent | 1.345 vs 1.356 | 1.345178 vs 1.355974 | 1.345 vs 1.356 | PASS |
| Folds where joint better | 12 / 25 | 12 / 25 | 12 / 25 | PASS |
| Sign-permutation P | 0.986 | 0.98626 | 0.986 | PASS |
| Min training cell across folds | 4 | 4 | 4 | PASS |

## Extended Data Fig. 1 — GEP18 positive control

| Item | Manuscript / SI | JSON (frozen) | Rounded | Verdict |
|---|---|---|---|---|
| GEP18 GSE165252 transport | 0.546 [0.403, 0.697] | `0.546196 [0.403367, 0.696561]` | 0.546 [0.403, 0.697] | PASS |
| GEP18 GSE78220 transport | 0.571 [0.318, 0.788] | `0.570513 [0.317886, 0.788462]` | 0.571 [0.318, 0.788] | PASS |
| GEP18 within-IMvigor210 (Efron 0.632) | 0.573 | 0.572831 | 0.573 | PASS |
| 6-gene IFN-γ within-IMvigor210 (Efron 0.632) | 0.605 | 0.604979 | 0.605 | PASS |
| GEP18 verdict (positive control) | FAIL (all external cohorts) | `any_external_cohort_passes = 0` across pm0.05/0.10/0.20 | FAIL | PASS |
| 18-gene list | see methods | `["CD8A","CMKLR1","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","CCL5","CXCL9","CXCL10","IDO1","PSMB10","STAT1","TIGIT","CD27","CD274","CD276","LAG3","PDCD1LG2"]` | — | PASS |

## Summary tally

- Total items checked: **64**
- PASS: **59**
- FAIL (manuscript text errors): **2**
  1. Fig. 4 P16: GSE78220 permutation null observed AUC should be **0.679**, not 0.680.
  2. Fig. 4 P16: permutation null is (Fig. 4d), not (Fig. 4c).
- FAIL (semantic error): **1**
  3. Fig. 2 P13: "concentrated around the source event rate" contradicts the median-observed gap of 0.166.
- FAIL (provenance): **1**
  4.  Fig. 2.png orphan (replaced in  with  frozen byte-identical).
- GAP (source provenance, not a manuscript error): **1**
  5. Fig. 5 TCR M1 / M4 numbers sourced from `joint_stress_test.json`; not recomputed in the  pipeline.
- GAP (source data completeness, non-blocking): **3**
  6. ED1 source_data JSON does not export the 8 transport CI bounds.
  7. ED5 source_data JSON is sparse (25 scalars) vs many displayed labels.
  8. Fig3 pooled AUC 0.615 and Fig5 labels 0.494/0.504/0.561 are cross-figure references.
