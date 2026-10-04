#  final-submission build notes

Generated at (UTC): 2026-10-04T08:36:20Z

## What changed from 

1. **Fig_2.png replaced.**  shipped an orphan PNG (md5 `ada64ae0d6ff1c8ec77340986c0238ca`, 357,106 B, 2026-10-04 15:17) that is not byte-identical to any frozen upstream. Replaced with  frozen PNG (md5 `0b97969e8fb0cf1ea9122a3b90e170ed`, 348,112 B, 2026-09-26 18:01).

2. **Source_Data/ added.**  contained no Source Data folder.  ships 16 JSON + 16 xlsx pairs, all copied byte-for-byte from `_20260926_Nature_Cancer/figures/source_data/`.

3. **frozen_manifest.json added.** SHA256 of every shipped file, plus the figure→source_data pairing table.

## Figure 2–6 number reconciliation — PASS/FAIL

See `RECONCILIATION_TABLE.md` for the full table. Headline:

- 21-config ranges (Fig. 4a): PASS (0.418–0.602 / 0.513–0.679)
- Variance decomposition (Fig. 4b): PASS (0.0028 / 0.0002 / 0.0011 / 0)
- Permutation null GSE165252 (Fig. 4d): PASS (0.602 vs 0.644, P 0.394)
- Permutation null GSE78220 (Fig. 4d): **FAIL** — manuscript writes `0.680`, JSON is `0.679487…` = **0.679**
- Panel reference for permutation null: **FAIL** — manuscript cites (Fig. 4c); caption defines (Fig. 4d)
- Fig. 2 calibration quantiles: PASS (0.310 / 0.474 / 0.394 / 0.228)
- Fig. 2 P13 wording 'concentrated around the source event rate': **FAIL** — median 0.394 vs observed 0.228 (Δ 0.166); wording contradicts data
- Fig. 3 calibration slopes + recalibration: PASS
- Fig. 5 24-cell adjudication + near-miss + clinical-only: PASS
- Fig. 6 joint boundary + toxicity stress: PASS
- GEP18 positive control (Extended Data Fig. 1): PASS

## Source data completeness gaps (do not block submission, but disclose)

- ED1 source_data JSON does not export the 8 CI bounds visible in the figure.
- ED5 source_data JSON is sparse (25 scalars) vs the many AUC-like labels on the figure.
- Fig3 pooled AUC 0.615 and Fig5 labels 0.494/0.504/0.561 are cross-figure references.

## Blocking item — Code availability DOI

The  manuscript (P52), Reporting Summary, Revision Notes, Cover Letter, and README all still carry a DOI placeholder. See `DOI_BLOCKER.md`.
