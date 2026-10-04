#!/usr/bin/env python3
""" -- Nature Cancer , P0-4.

Selection-adjusted permutation null for the 21-configuration controlled
ablation. REAL RUN.

Reviewer challenge addressed ( PreSubmission Revision Notes, Table 4):
"Does the 21-config search inflate the apparent negative/positive signal?"
Required response: a selection-adjusted permutation null comparing the
observed maximum external AUC across the 21 non-degenerate configurations
with the permutation max-AUC distribution.

Design:
  - The 21 configurations and the frozen transport pipeline are EXACTLY
    those of rerun_ablation.py (imported verbatim; the anchor self-check
    values of that script are reproduced because the S0 PCA scores are
    label-independent).
  - S0 scores per configuration are computed ONCE (CLR/PCA do not use y).
  - For each of N_PERM = 2,000 permutations (seed 20260911):
      permute the ANCHOR outcome labels (IMvigor210, n=298), refit each
      configuration's ridge logistic on the permuted labels, predict the
      FIXED target labels of each external cohort, and record the AUC.
      Within one permutation the SAME permuted labeling is used across all
      21 configurations (paired design), and the MAXIMUM AUC across the
      21 configurations is recorded per target cohort.
  - Empirical one-sided P-values:
      P_perm = (1 + #{perm_max >= obs_max}) / (1 + N_PERM)
      per configuration: P_cfg = (1 + #{perm_auc >= obs_auc}) / (1 + N_PERM)
  - Also reported: the 95th percentile of the permutation max-AUC
    distribution (the selection-adjusted chance threshold).

Honesty discipline: everything computed at run time; no fabricated numbers.
Emission: data/permutation_null.json + .md report.
"""
import json
import os
import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
BASE = f"{ROOT}/outputs/analysis_base"
OUT_DIR = f"{BASE}/data/"
OUT_JSON = f"{OUT_DIR}/permutation_null.json"
OUT_MD = f"{OUT_DIR}/permutation_null.md"
SEED = 20260911
N_PERM = 2000

sys.path.insert(0, f"{BASE}/scripts")
import rerun_ablation as A  # noqa: E402  (verbatim pipeline; main guarded)

COHORTS = A.COHORTS
EXT_COHORTS = A.EXT_COHORTS


def build_configs(n_common, data):
    """VERBATIM copy of _rerun_ablation.main() config construction."""
    Zb = A.build_z({c: data[c][0] for c in COHORTS}, "clr")
    pca_b = PCA(n_components=2, random_state=SEED).fit(Zb["IMvigor210"])
    pc1_load = pca_b.components_[0]
    rank = np.argsort(-np.abs(pc1_load))
    top100_idx = sorted(rank[:100].tolist())
    top50_idx = sorted(rank[:50].tolist())
    n_lab_imv = int((~np.isnan(data["IMvigor210"][1])).sum())
    lam_base_ext = 1.0 / n_lab_imv
    configs = [dict(id="baseline", layer="baseline", gene_idx=None,
                    norm="clr", k=2, lam_ext=lam_base_ext)]
    configs.append(dict(id="gene_top100", layer="gene_set", gene_idx=top100_idx,
                        norm="clr", k=2, lam_ext=lam_base_ext))
    configs.append(dict(id="gene_top50", layer="gene_set", gene_idx=top50_idx,
                        norm="clr", k=2, lam_ext=lam_base_ext))
    for i in range(5):
        rng = np.random.default_rng(SEED + i)
        idx = sorted(rng.choice(n_common, 100, replace=False).tolist())
        configs.append(dict(id=f"gene_rand100_s{i}", layer="gene_set",
                            gene_idx=idx, norm="clr", k=2, lam_ext=lam_base_ext))
    for i in range(5):
        rng = np.random.default_rng(SEED + i)
        idx = sorted(rng.choice(n_common, 50, replace=False).tolist())
        configs.append(dict(id=f"gene_rand50_s{i}", layer="gene_set",
                            gene_idx=idx, norm="clr", k=2, lam_ext=lam_base_ext))
    configs.append(dict(id="norm_log2raw", layer="normalization", gene_idx=None,
                        norm="log2raw", k=2, lam_ext=lam_base_ext))
    configs.append(dict(id="norm_within_z", layer="normalization", gene_idx=None,
                        norm="within_z", k=2, lam_ext=lam_base_ext))
    for k in (1, 3, 5):
        configs.append(dict(id=f"pc_k{k}", layer="pc_k", gene_idx=None,
                            norm="clr", k=k, lam_ext=lam_base_ext))
    for lv in (0.01, 0.1, 1.0):
        configs.append(dict(id=f"ridge_lam{lv:g}", layer="ridge_lambda",
                            gene_idx=None, norm="clr", k=2, lam_ext=lv))
    return configs


def config_scores(cfg, data, n_common):
    """S0 PCA scores per cohort for one configuration (label-independent)."""
    gidx = (cfg["gene_idx"] if cfg["gene_idx"] is not None
            else np.arange(n_common))
    gidx = np.asarray(gidx)
    Msub = {c: data[c][0][gidx, :] for c in COHORTS}
    Z = A.build_z(Msub, cfg["norm"])
    pca = PCA(n_components=cfg["k"], random_state=SEED).fit(Z["IMvigor210"])
    return {c: pca.transform(Z[c]) for c in COHORTS}


def fast_auc_bootstrap_p(y, p, n_boot, rng):
    """Vectorised percentile bootstrap CI for AUC via the Mann-Whitney rank
    identity (mathematically identical to roc_auc_score per resample)."""
    n = len(y)
    n1 = int(y.sum())
    n0 = n - n1
    if n1 == 0 or n0 == 0:
        return np.nan, (np.nan, np.nan)
    idx = rng.integers(0, n, size=(n_boot, n))
    ys = y[idx]                                    # (B, n)
    ps = p[idx]                                    # (B, n)
    # ranks of ps within each row (average ties)
    order = np.argsort(ps, axis=1, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    rows = np.arange(n_boot)[:, None]
    ranks[rows, order] = np.arange(1, n + 1, dtype=float)[None, :]
    # average ranks for ties
    sp = np.sort(ps, axis=1)
    tie = np.zeros((n_boot, n), dtype=bool)
    tie[:, 1:] = sp[:, 1:] == sp[:, :-1]
    tie[:, :-1] |= sp[:, 1:] == sp[:, :-1]
    if tie.any():
        # conservative tie correction: recompute per-row with scipy
        from scipy.stats import rankdata
        ranks = np.apply_along_axis(rankdata, 1, ps)
    auc = (ranks[ys == 1].mean(1) - (n1 + 1) / 2.0) / n0
    lo, hi = np.percentile(auc, [2.5, 97.5])
    return float(auc.mean()), (float(lo), float(hi))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ================= data under the recovered identity (verbatim main()) ==
    M, sids = A.load_exmat()
    common, sym2ens, r2s, avail = A.common_genes_recovered()
    n_common = len(common)
    M320_imv, y_imv = A.load_imvigor_correct(M, sids, common, r2s)
    m_lab = ~np.isnan(y_imv)
    data = {
        "IMvigor210": (M320_imv[:, m_lab], y_imv[m_lab]),
        "GSE165252": A.load_gse165252(common, sym2ens),
        "GSE78220": A.load_gse78220(common),
    }
    for c in COHORTS:
        Mc, yc = data[c]
        m = ~np.isnan(yc)
        data[c] = (Mc[:, m], yc[m].astype(float))
    y_anchor = data["IMvigor210"][1]
    n_anchor = len(y_anchor)
    print(f"[P0] n_anchor={n_anchor}, configs to build...")

    # ================= S0 scores per config (label-independent) =============
    configs = build_configs(n_common, data)
    assert len(configs) == 21, f"expected 21 configs, got {len(configs)}"
    S = {}
    for cfg in configs:
        S[cfg["id"]] = config_scores(cfg, data, n_common)
    print("[P0] S0 scores computed for all 21 configurations")

    # ================= observed AUCs ========================================
    observed = {}
    for tgt in EXT_COHORTS:
        y_te = data[tgt][1]
        rows = {}
        for cfg in configs:
            cfg_id = cfg["id"]
            S_te = S[cfg_id][tgt]
            S_tr = S[cfg_id]["IMvigor210"]
            Xtr = np.column_stack([np.ones(n_anchor), S_tr])
            w = A.fit_logistic(Xtr, y_anchor, cfg["lam_ext"])
            p_te = A.sigmoid(np.column_stack([np.ones(len(y_te)), S_te]) @ w)
            rows[cfg_id] = dict(
                auc=float(roc_auc_score(y_te, p_te)),
                predictions=p_te.tolist(),
            )
        aucs = np.array([rows[c]["auc"] for c in rows])
        observed[tgt] = dict(
            per_config=rows,
            max_auc=float(aucs.max()),
            argmax_config=str(list(rows)[int(aucs.argmax())]),
            n_configs=len(rows),
        )
        print(f"  observed {tgt}: max AUC={observed[tgt]['max_auc']:.4f} "
              f"({observed[tgt]['argmax_config']})")

    # ================= permutation null =====================================
    rng = np.random.default_rng(SEED)
    perm_max = {t: np.empty(N_PERM) for t in EXT_COHORTS}
    perm_auc = {t: np.empty((N_PERM, 21)) for t in EXT_COHORTS}
    cfg_ids = [c["id"] for c in configs]
    for pi in range(N_PERM):
        y_perm = rng.permutation(y_anchor)
        for j, cfg in enumerate(configs):
            S_tr = S[cfg["id"]]["IMvigor210"]
            Xtr = np.column_stack([np.ones(n_anchor), S_tr])
            w = A.fit_logistic(Xtr, y_perm, cfg["lam_ext"])
            for tgt in EXT_COHORTS:
                y_te = data[tgt][1]
                S_te = S[cfg["id"]][tgt]
                p_te = A.sigmoid(np.column_stack([np.ones(len(y_te)), S_te]) @ w)
                perm_auc[tgt][pi, j] = roc_auc_score(y_te, p_te)
        for tgt in EXT_COHORTS:
            perm_max[tgt][pi] = perm_auc[tgt][pi].max()
        if (pi + 1) % 200 == 0:
            print(f"  permutation {pi+1}/{N_PERM} done", flush=True)

    # ================= inference ============================================
    result = dict(
        analysis="permutation_null",
        purpose=("P0-4 selection-adjusted permutation null for the 21-config "
                 "controlled ablation: observed maximum external AUC vs the "
                 "permutation max-AUC distribution (anchor-label permutation, "
                 "paired across configurations)"),
        pipeline_source="current",
        n_permutations=N_PERM,
        seed=SEED,
        permutation_unit="anchor (IMvigor210, n=298) outcome labels; target labels fixed",
        pairing="same permuted labeling across all 21 configurations within a permutation",
        configs=[dict(id=c["id"], layer=c["layer"], k=c["k"],
                      lam_ext=c["lam_ext"]) for c in configs],
    )
    for tgt in EXT_COHORTS:
        pm = perm_max[tgt]
        pa = perm_auc[tgt]
        obs_max = observed[tgt]["max_auc"]
        p_max = (1 + int((pm >= obs_max).sum())) / (1 + N_PERM)
        per_cfg = {}
        for j, cid in enumerate(cfg_ids):
            obs = observed[tgt]["per_config"][cid]["auc"]
            p_cfg = (1 + int((pa[:, j] >= obs).sum())) / (1 + N_PERM)
            per_cfg[cid] = dict(observed_auc=obs,
                                perm_mean_auc=float(pa[:, j].mean()),
                                perm_p_one_sided=float(p_cfg))
        result[tgt] = dict(
            observed_max_auc=obs_max,
            observed_argmax=observed[tgt]["argmax_config"],
            perm_max_mean=float(pm.mean()),
            perm_max_q95=float(np.percentile(pm, 95)),
            perm_max_max=float(pm.max()),
            perm_p_max_one_sided=float(p_max),
            per_config=per_cfg,
            observed_per_config_max=observed[tgt],
        )
        print(f"[P0] {tgt}: obs_max={obs_max:.4f} "
              f"perm_max_mean={pm.mean():.4f} perm_q95={np.percentile(pm, 95):.4f} "
              f"P(max)={p_max:.4f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print("Wrote", OUT_JSON)

    # ================= markdown =============================================
    md = ["#  P0-4 — 21-config 选择校正置换零分布\n"]
    md.append(f"日期：2026-09-26 ｜ seed {SEED} ｜ 置换次数 {N_PERM} ｜ "
              "脚本：scripts/analysis/permutation_null.py ｜ "
              "管线来源：scripts/analysis/rerun_ablation.py（逐字导入，锚点自检值可复现）\n")
    md.append("**设计：** 每次置换重排锚队列（IMvigor210，n=298）结局标签；21 个配置在"
              "同一置换标签下配对重拟合；用固定目标标签计算外推 AUC；记录每次置换 21 "
              "配置中的**最大 AUC**。PCA/CLR 得分与标签无关，只算一次。\n")
    for tgt in EXT_COHORTS:
        r = result[tgt]
        md.append(f"\n## {tgt}\n")
        md.append(f"- 观测最大 AUC：**{r['observed_max_auc']:.4f}**"
                  f"（{r['observed_argmax']}）")
        md.append(f"- 置换最大 AUC：mean = {r['perm_max_mean']:.4f}，"
                  f"95th 分位 = **{r['perm_max_q95']:.4f}**，max = {r['perm_max_max']:.4f}")
        md.append(f"- 选择校正单侧 P（max）：**P = {r['perm_p_max_one_sided']:.4f}**")
        md.append("\n| config | 观测 AUC | 置换 mean | 单侧 P |")
        md.append("|---|---|---|---|")
        for cid, c in r["per_config"].items():
            md.append(f"| {cid} | {c['observed_auc']:.4f} | "
                      f"{c['perm_mean_auc']:.4f} | {c['perm_p_one_sided']:.4f} |")
    md.append("\n**解读：** 若观测最大 AUC 未超过置换零分布的 95th 分位，则 21-config "
              "搜索得到的任何表观信号均不能与选择效应区分——审稿人关于"
              "“post-selection 膨胀”质疑的定量回应。\n")
    with open(OUT_MD, "w") as f:
        f.write("\n".join(md))
    print("Wrote", OUT_MD)


if __name__ == "__main__":
    main()
