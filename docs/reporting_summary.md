# Reporting Summary —  (Nature Cancer, Analysis)

Manuscript: "Decomposing the transportability of immune-state biomarkers across checkpoint immunotherapy"
Package: _20260926_Nature_Cancer

## Statistics & general

- n of independently replicated analyses: one primary frozen-transport analysis (2 external targets), one LOCO sensitivity benchmark (8 models x 3 targets = 24 cells), one 21-configuration representation audit, one permutation null (2,000 permutations x 21 configurations x 2 targets), one precision simulation (5,000 simulations x 2,000 bootstrap CIs x 7 true-AUC values x 2 target sizes), two positive controls (18-gene GEP primary; 6-gene IFN-gamma sensitivity), one joint efficacy-toxicity identifiability analysis (GSE186143, n = 47 joint-evaluable).
- All analyses were pre-specified in the versioned analysis scripts before benchmark execution on the recovered gene identity; they were NOT prospectively preregistered. The two-criterion licensing rule (AUC CI excludes 0.5 beneficially AND |slope-1| <= 0.10 logit) is a declared reporting convention; threshold sensitivity (ED6, Supplementary Note 2) shows no verdict depends on the 0.10 tolerance.
- Exclusion criteria and exclusions: GSE287541 withheld from all estimates (longitudinal design, no usable endpoint crosswalk); five former "measured-only" ablation configurations voided under the recovered gene identity (documented retraction); GSE135222 excluded (no response labels). All exclusions are stated in the manuscript and Supplementary Note 1/3.

## Data flow & reproducibility

- Every number in the manuscript, figures, Extended Data and Supplementary Tables is read programmatically from frozen JSON archives (data/, data/; seed 20260911 throughout; SHA-256 checksums in analysis_manifest.json). No value is hand-typed.
- Figures: R 4.6.1, ggplot2 4.0.3, theme_classic, Arial, 300 dpi; PNG + PDF + SVG per figure; code in scripts/figures/make_figures.R, make_anchor_figures.R, make_stable_figures.R (helper: scripts/ggbiomed.R). No AI-generated figures.
- One end-to-end build from the corrected source-of-truth to final figures was executed before submission (P0-10). The gene-row identity recovery event and the retraction of the affected configurations are documented (Supplementary Note 4, ED5, ED10).
- Python 3.13 (NumPy/SciPy/scikit-learn); environment pinned in requirements.txt / requirements_R.txt.

## Prespecified endpoints vs reported

- Primary endpoint: external AUC of the frozen IMvigor210-anchored S0 representation with BCa 95% CI. Reported: 0.504 [0.354, 0.658] (GSE165252), 0.660 [0.404, 0.862] (GSE78220).
- Secondary endpoints: calibration intercept/slope, Brier score, log-loss, recalibration behaviour, representation-configuration grid, positive controls, joint efficacy-toxicity identifiability. All reported in Results, Extended Data and Tables S3-S8.
- Deviations: none from the declared analysis plan. A data-mapping error affecting the anchor cohort was discovered mid-analysis, corrected, and every affected number re-run (P0-1); the retraction is disclosed rather than silently superseded.

## Sample sizes & power

- Observed: IMvigor210 n = 298 (68 events); GSE165252 n = 71 (23 events); GSE78220 n = 25 (13 events); GSE186143 n = 47 joint-evaluable (27 DCB / 20 NDB; 15 severe irAE).
- Precision simulation (ED9, Table S8a): at n = 71, probability of excluding AUC = 0.5 exceeds 90% only for true AUC >= 0.75; at n = 25, even true AUC 0.80 gives 77% power. The study cannot establish transport at these sizes, nor rule out moderate true effects; this is stated in the manuscript.
- Joint analysis identifiability: 5.2 events per effective parameter; 25-fold CV log-loss difference between joint and independent models not distinguishable from zero (sign-permutation P = 0.986).

## Multiple comparisons & selection

- The 21-configuration search is selection-adjusted: observed maximum external AUC vs permutation max-AUC distribution (2,000 anchor-label permutations, paired across configurations): one-sided P = 0.394 (GSE165252) and 0.398 (GSE78220). All configurations are retained in Supplementary Table S5; none is selected post hoc.

## LLM use disclosure

- Large-language-model assistants were used for code drafting support, language editing and internal review checklists under human direction. They were not used to generate scientific content, analyses, numbers or figures; all results derive from the archived deterministic code and frozen data.

## Conflicts, funding, data

- Funding: [TO BE CONFIRMED AT SUBMISSION]. Competing interests: none declared. All data public (see data_availability_statement.txt); no restrictions.
