# Decomposing the transportability of immune-state biomarkers across checkpoint immunotherapy

Reproduction repository for the *Nature Cancer* Analysis submission. All numbers, figures, and source data referenced by the manuscript are archived here in a frozen, byte-checksummed state.

**Target:** *Nature Cancer*, Analysis track
**Abstract:** 150 words, unreferenced
**Main text:** 1,909 / 4,000 words
**Methods:** 627 / 3,000 words
**Display items:** 6 main figures + 10 Extended Data figures
**Random seed (all analyses):** 20260911

## Scientific claim

Under a prespecified frozen transport estimand, the immune-state representation did not establish reproducible external transport in the tested independent checkpoint-immunotherapy cohorts. No stronger claim (biology absent / non-transportable in principle / single-cause attribution) is made anywhere in the repository.

## Repository layout

| Path | Content |
|---|---|
| `submission/` | Full *Nature Cancer* submission package: manuscript, supplementary information, cover letter, reporting summary, checklist, revision notes, figure production specification, all 16 display items in PNG/PDF/SVG, per-figure source data in JSON and XLSX, and a frozen SHA-256 manifest. |
| `manuscript/manuscript.md` | Markdown source of truth for the main manuscript. |
| `supplementary/supplementary.md` | Markdown source for Supplementary Notes 1–6 and Tables S1–S8. |
| `scripts/analysis/` | Python analysis: 18-gene GEP positive control, permutation null, precision simulation, clinical-only, gene-row recovery, transport summary, benchmark, ablation, sensitivity, CDS extract, phase-1a panel checks. |
| `scripts/figures/` | R 4.6.1 + ggplot2 plotting scripts (6 main + 10 ED figures) and the local `ggbiomed.R` helper. |
| `scripts/doc/` | Builders that render `manuscript.md` and `supplementary.md` into `.docx`. |
| `scripts/compliance/compliance_check.py` | Automated pre-submission compliance checker. |
| `data/` | Frozen JSON results from every rerun (benchmark, ablation, permutation null, positive control, precision simulation, transport summary, sensitivity, provenance ledger, number audit, analysis manifest). |
| `figures/` | 16 display items × 3 formats (PDF, PNG 300 dpi, SVG), plus `source_data/` — one per-figure JSON archive. |
| `docs/` | Cover letter, reporting summary, data availability statement (plaintext). |
| `requirements.txt`, `requirements_R.txt` | Python and R dependency pin. |

## Honest-result notes (retained intentionally)

1. The 18-gene T-cell-inflamed GEP positive control **fails** the two-criterion rule in both external cohorts (AUC 0.546 and 0.571) while passing within the anchor (0.654). This is presented as evidence that the evaluation environment is not a trivial rejector — not as a validated positive transport.
2. The selection-adjusted permutation null shows the 21-configuration maximum is **indistinguishable from selection alone** (P = 0.394 / 0.398).
3. The precision simulation shows n = 25 / n = 71 cannot establish transport, and also cannot rule out moderate true effects. Both directions are stated.
4. A gene-row identity mapping error from an earlier pipeline generation was retracted and corrected before submission; the wrong-identity archives are retained only as provenance.

## Reproduction

All results are reproducible from the frozen JSONs under `data/`. Patient-level S0 representation matrices used by the rerun series are derived from public GEO accessions (see `docs/data_availability_statement.txt`); only the derived representation, not raw counts, is required to re-run the analyses.

```bash
# 1. Python analyses (seed 20260911)
python scripts/analysis/positive_control_gep18.py
python scripts/analysis/permutation_null.py
python scripts/analysis/precision_simulation.py
python scripts/analysis/rerun_benchmark.py
python scripts/analysis/rerun_ablation.py
python scripts/analysis/rerun_positive_control.py
python scripts/analysis/rerun_sensitivity_manifest.py
python scripts/analysis/imvigor_gene_recovery.py
python scripts/analysis/transport_summary_joint.py
python scripts/analysis/clinical_only.py
python scripts/analysis/threshold_zero_sensitivity.py
python scripts/analysis/cds_extract.py

# 2. Figures (R 4.6.1 + ggplot2 4.0.3)
Rscript scripts/figures/make_anchor_figures.R
Rscript scripts/figures/make_stable_figures.R
Rscript scripts/figures/make_final_figures.R

# 3. Documents
python scripts/doc/make_manuscript_docx.py
python scripts/doc/make_supplementary_docx.py

# 4. Pre-submission compliance
python scripts/compliance/compliance_check.py
```

## Data availability

All primary data are public GEO accessions (IMvigor210 GSE139551; external targets GSE165252, GSE78220; joint efficacy–toxicity GSE186143/GSE186144). The corrected gene-row identity key for the IMvigor210 deposit is included with the analysis manifest. All GEO accessions were re-verified field-by-field against GEO SOFT records during the provenance audit. Full statement in `docs/data_availability_statement.txt`.

## Code availability

This repository. Every number quoted in the manuscript is reproducible from the frozen JSONs under `data/` and the analysis scripts under `scripts/`.

## License

Analysis code: Apache License 2.0.
Manuscript, figures, and results: © the authors, available under the *Nature* publishing terms.
