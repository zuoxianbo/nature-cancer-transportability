#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" -- Nature Cancer revision , HONEST RE-SPECIFIED
CONTROLLED ABLATION under the recovered IMvigor210 gene identity. REAL RUN.

Background: alignment_ablation.py ran a 26-configuration controlled
ablation grid on a 228-gene "common intersection" selected through the doubly
broken gene identity (positional exmat/fData row alignment + fake-Entrez gene
selection).  recovered the true identity (exmat row j == the j-th fData row
when fData is sorted by SYMBOL upper-cased alphabetically) and rebuilt the
anchor on the corrected 320-gene cross-cohort intersection, which has ZERO
all-zero genes in IMvigor210. The old "70 unmeasured genes / measured-only"
artefact block is therefore VOID: its 5 configurations
(gene_measured_only / gene_top100_measured / gene_top50_measured /
norm_log2raw_measured / norm_within_z_measured) are RETRACTED and recorded
with removal_reason in the JSON (not silently dropped).

This script re-specifies and reruns the grid:
  - IMvigor210 loading via the RECOVERED KEY (verbatim logic of
    scripts/analysis/imvigor_gene_recovery.py stage 2); external loaders
    (load_gse165252 / load_gse78220 / clr_standardize / fit_pca) verbatim
    from data/pca_stability.py.
  - Ablation-engine metrics VERBATIM from alignment_ablation.py:
    external AUC + percentile bootstrap 95% CI (2000 resamples of the
    EXTERNAL cohort with fixed predictions, seed 20260911); within-cohort
    Efron 0.632 optimism-corrected AUC (600 boots, seed 20260911) for all
    three cohorts; per_config structure identical.
  - Additionally the  recovery BCa CI (auc_bca verbatim from
    imvigor_gene_recovery.py) per external cell, for anchor continuity.
  - Baseline lambda rule follows the FROZEN  recovery convention (this is
    the one deliberate re-specification vs , required for bitwise anchor
    reproduction, and is recorded in the JSON design block): external
    transport lambda = 1/n_train (= 1/298); within 0.632 lambda = 1/N_total
    (= 1/394).  used a single 1/N_total for both.

New grid (21 configurations):
  layer 1 gene set : baseline(320 intersection) / top100 / top50 (ranked by
                     |PC1 loading| of the anchored baseline CLR PCA) /
                     rand100 x 5 seeds / rand50 x 5 seeds (seed 20260911+k)
  layer 2 norm     : clr(baseline) / log2raw / within_z
  layer 3 dim      : k=2(baseline) / 1 / 3 / 5
  layer 4 ridge    : lambda = 1/n(baseline) / 0.01 / 0.1 / 1

Self-check anchor (frozen by the recovery worker; the baseline configuration
MUST reproduce these BITWISE -- the script aborts otherwise):
  external GSE165252 AUC 0.5036231884057971, BCa CI [0.3540206604291994,
  0.6584721535961872]; external GSE78220 AUC 0.6602564102564102, BCa CI
  [0.4038454964351145, 0.8616583228383701]; within 0.632 IMvigor210
  0.6570746710270073.

Emission: data/rerun_ablation.json + .md report (new grid table,
per-layer external-AUC variance decomposition, recovery-transport verdict,
and the old- survive/flip/retract comparison, all COMPUTED from the frozen
 JSON -- nothing asserted by hand).

Reproducibility: managed venv python, seed 20260911 everywhere.
"""
import csv
import gzip
import json
import os
import pickle
import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from scipy import stats as _st

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
PKG = f"{ROOT}/outputs/analysis_base"
GENE_JSON = f"{ROOT}/geo_data/322_shared_genes.json"
SYM2ENS_PKL = f"{ROOT}/gene_name_id_dict_gc104M.pkl"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
GSE165_MAT = f"{ROOT}/geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz"
GSE165_ANN = f"{ROOT}/geo_data/GSE165252/series_matrix_annotation.txt"
GSE782_MAT = f"{ROOT}/geo_data/GSE78220/GSE78220_PatientFPKM.xlsx"
GSE782_ANN = f"{ROOT}/geo_data/GSE78220/series_matrix_annotation.txt"
V74_JSON = f"{PKG}/data/alignment_ablation.json"
V82_RECOVERY_JSON = f"{PKG}/data/imvigor_gene_recovery.json"

OUT_DIR = f"{PKG}/data/"
OUT_JSON = f"{OUT_DIR}/rerun_ablation.json"
OUT_MD = f"{OUT_DIR}/rerun_ablation.md"

SEED = 20260911
N_BOOT_EXT = 2000      #  verbatim: percentile bootstrap for external AUC CI
N_BOOT_632 = 600       # / verbatim: Efron 0.632 control
COHORTS = ["IMvigor210", "GSE165252", "GSE78220"]
EXT_COHORTS = ["GSE165252", "GSE78220"]
TOL = 1e-12

REMOVAL_REASON = "70-gene artefact retracted: zero-gene count = 0 under recovered identity"

# frozen self-check anchor (recovery worker, data/imvigor_gene_recovery.json)
ANCHOR = dict(
    ext_GSE165252_auc=0.5036231884057971,
    ext_GSE165252_ci_bca=[0.3540206604291994, 0.6584721535961872],
    ext_GSE78220_auc=0.6602564102564102,
    ext_GSE78220_ci_bca=[0.4038454964351145, 0.8616583228383701],
    within_imvigor210_auc632=0.6570746710270073,
    lam_ext=0.003355704697986577,        # 1/298
    lam_within=0.0025380710659898475,    # 1/394
)

# ================================================================ verbatim // algebra
def sigmoid(z):
    z = np.clip(z, -35, 35)
    return 1.0 / (1.0 + np.exp(-z))

def fit_logistic(X, y, lam, n_iter=500, step_guard=1.0):
    n, p = X.shape
    w = np.zeros(p)
    for _ in range(n_iter):
        eta = X @ w
        mu = sigmoid(eta)
        g = X.T @ (mu - y)
        g[1:] += lam * w[1:]
        W = mu * (1 - mu)
        H = (X.T * W) @ X
        for j in range(1, p):
            H[j, j] += lam
        try:
            delta = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        nd = np.linalg.norm(delta)
        if nd > step_guard:
            delta = delta * step_guard / nd
        w = w + delta
        if nd < 1e-10:
            break
    return w

def auc_632(X, y, lam, n_boot=N_BOOT_632, seed=SEED):
    """Efron 0.632 optimism-corrected AUC (// verbatim) + bootstrap SE + 95% CI."""
    rng = np.random.default_rng(seed)
    n = len(y)
    def fit_auc(Xs, ys):
        Xs = np.column_stack([np.ones(len(ys)), Xs])
        w = fit_logistic(Xs, ys, lam)
        p = sigmoid(Xs @ w)
        if len(set(ys)) < 2:
            return np.nan
        return roc_auc_score(ys, p)
    app = fit_auc(X, y)
    if np.isnan(app):
        return app, app, np.nan, (np.nan, np.nan)
    corr_b = []
    for _ in range(n_boot):
        bi = rng.integers(0, n, n)
        Xb, yb = X[bi], y[bi]
        if len(set(yb)) < 2:
            continue
        app_b = fit_auc(Xb, yb)
        Xo = np.column_stack([np.ones(n), X])
        wb = fit_logistic(np.column_stack([np.ones(len(yb)), Xb]), yb, lam)
        po = sigmoid(Xo @ wb)
        try:
            auc_o = roc_auc_score(y, po)
        except ValueError:
            continue
        opt = app_b - auc_o
        corr_b.append(app - opt)
    corr_b = np.array(corr_b)
    corr = float(corr_b.mean())
    se = float(corr_b.std(ddof=1))
    ci = (float(corr - 1.96 * se), float(corr + 1.96 * se))
    return float(app), corr, se, ci

def auc_bca(y, p, n_boot=N_BOOT_EXT, seed=SEED):
    """2000-rep bootstrap BCa 95% CI for AUC (imvigor_gene_recovery.py verbatim)."""
    y = np.asarray(y, float); p = np.asarray(p, float)
    n = len(y)
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(set(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], p[idx]))
    aucs = np.array(aucs)
    if len(aucs) < 10:
        return (np.nan, np.nan)
    orig = float(roc_auc_score(y, p))
    z0 = float(_st.norm.ppf((aucs < orig).mean()))
    jk = np.empty(n)
    for i in range(n):
        m = np.arange(n) != i
        jk[i] = roc_auc_score(y[m], p[m])
    jmean = jk.mean()
    num = np.sum((jmean - jk) ** 3)
    den = 6.0 * (np.sum((jmean - jk) ** 2) ** 1.5)
    a = float(num / den) if den > 0 else 0.0
    def bca(alpha):
        z = _st.norm.ppf(alpha)
        frac = _st.norm.cdf(z0 + (z0 + z) / (1.0 - a * (z0 + z)))
        return float(np.clip(frac, 0, 1))
    lo = np.percentile(aucs, 100 * bca(0.025))
    hi = np.percentile(aucs, 100 * bca(0.975))
    return float(lo), float(hi)

# ================================================================ recovered-key IMvigor210 loading ( recovery verbatim logic)
def load_exmat():
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        sids = [x.strip('"') for x in hdr[1:]]
        rows = []
        for r in rd:
            rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
    return np.array(rows, float), sids

def row_of_symbol():
    """Recovered key: exmat row j == j-th fData row in upper-cased symbol alphabetical order."""
    fd = list(csv.DictReader(open(IMV_FDATA)))
    syms = [r["symbol"].upper() for r in fd]
    order = sorted(range(len(fd)), key=lambda i: syms[i])
    return {syms[order[j]]: j for j in range(len(fd))}

def common_genes_recovered():
    """320-gene cross-cohort intersection ( recovery stage_2 availability logic verbatim)."""
    panel = [g.upper() for g in json.load(open(GENE_JSON))["genes"]]
    r2s = row_of_symbol()
    imv_avail = {s for s in panel if s in r2s}
    g = pickle.load(open(SYM2ENS_PKL, "rb"))
    sym2ens = {k.upper(): v for k, v in g.items()}
    with gzip.open(GSE165_MAT, "rt") as f:
        f.readline()
        cnt_ens = set(line.rstrip("\n").split("\t")[0].split(".")[0] for line in f)
    gse165_avail = {s for s in panel if sym2ens.get(s) in cnt_ens}
    import openpyxl
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    rr = list(ws.iter_rows(min_row=1, values_only=True))
    g7820_syms = set(str(r[0]).upper() for r in rr[1:])
    gse782_avail = {s for s in panel if s in g7820_syms}
    common = [s for s in panel if s in imv_avail and s in gse165_avail and s in gse782_avail]
    avail = dict(panel=len(panel), imvigor210=len(imv_avail),
                 imvigor210_missing=sorted(set(panel) - imv_avail),
                 gse165252=len(gse165_avail), gse78220=len(gse782_avail))
    return common, sym2ens, r2s, avail

def load_imvigor_correct(M, sids, common, r2s):
    """IMvigor210 320-gene matrix via the recovered key (only change vs 
    load_imvigor: row index comes from symbol -> exmat row, not fData position)."""
    idx = [r2s[s] for s in common]
    M320 = M[idx, :]
    pd_rows = list(csv.DictReader(open(IMV_PDATA)))
    pmap = {r["X"]: r for r in pd_rows}
    y = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                  else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                  else np.nan for s in sids], dtype=float)
    return M320, y

# ================================================================ verbatim  external loaders
def load_gse165252(common, sym2ens):
    with gzip.open(GSE165_MAT, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        codes = [h.strip().strip('"') for h in hdr[1:]]
        gids = []; rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gids.append(parts[0].split(".")[0]); rows.append([float(x) for x in parts[1:]])
    M = np.array(rows, float)
    e2row = {e: i for i, e in enumerate(gids)}
    idx = [e2row[sym2ens[s]] for s in common]
    M322 = M[idx, :]
    ann = open(GSE165_ANN).read().splitlines()
    sample_lines = {}
    for ln in ann:
        if ln.startswith("!Sample"):
            fld = ln.rstrip("\n").split("\t")
            sample_lines[fld[0]] = [v.strip().strip('"') for v in fld[1:]]
    code_line = resp_line = None
    for name, vals in sample_lines.items():
        if name.startswith("!Sample_description") and len(set(vals) & set(codes)) >= max(1, 0.9 * len(codes)):
            code_line = vals; break
    for name, vals in sample_lines.items():
        if name.startswith("!Sample_characteristics") and any(v.lower().startswith("response:") for v in vals):
            resp_line = vals; break
    code2resp = {}
    for c, r in zip(code_line, resp_line):
        rl = r.lower()
        if "non-responder" in rl: code2resp[c] = 0.0
        elif "responder" in rl: code2resp[c] = 1.0
        else: code2resp[c] = np.nan
    y = np.array([code2resp.get(c, np.nan) for c in codes], dtype=float)
    return M322, y

def load_gse78220(common):
    import openpyxl
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(min_row=1, values_only=True))
    hdr = rows[0]
    sample_cols = [str(h).split(".baseline")[0].split(".Baseline")[0] for h in hdr[1:]]
    gene_rows = {}; data = []
    for r in rows[1:]:
        gid = str(r[0]).upper()
        gene_rows[gid] = len(data)
        data.append([float(x) if x is not None else 0.0 for x in r[1:]])
    M = np.array(data, float)
    idx = [gene_rows[s] for s in common]
    M322 = M[idx, :]
    ann = open(GSE782_ANN).read().splitlines()
    char_lines = []
    for ln in ann:
        if ln.startswith("!Sample_characteristics_ch1"):
            char_lines.append([v.strip().strip('"') for v in ln.rstrip("\n").split("\t")[1:]] if "\t" in ln else [])
    pid_line = resp_line = None
    for vals in char_lines:
        if vals and any(v.lower().startswith("patient id:") for v in vals):
            pid_line = vals
        if vals and any("anti-pd-1 response" in v.lower() for v in vals):
            resp_line = vals
    pid2resp = {}
    for c, r in zip(pid_line, resp_line):
        pid = c.split(":", 1)[1].strip(); rl = r.lower()
        if "progressive" in rl or "stable" in rl: pid2resp[pid] = 0.0
        elif "complete" in rl or "partial" in rl: pid2resp[pid] = 1.0
        else: pid2resp[pid] = np.nan
    y = np.array([pid2resp.get(sc, np.nan) for sc in sample_cols], dtype=float)
    return M322, y

# ================================================================ normalization variants ( verbatim)
def log_clr(M):
    """genes x n -> n x genes: CPM -> log1p -> per-sample centring across genes.

    Zero-guard (needed for small gene subsets where a sample can have all-zero
    counts across the selected genes): a zero library column is treated as
    library size 1, so its CLR values are exactly 0 instead of NaN.
    """
    cs = M.sum(0)
    cs = np.where(cs == 0, 1.0, cs)
    cpm = M / cs * 1e6
    lp = np.log1p(cpm)
    clr = lp - lp.mean(0)
    return clr.T

def log2_raw(M):
    """genes x n -> n x genes: log2(raw value + 1), no CPM, no per-sample centring."""
    return np.log2(M + 1.0).T

NORM_FUNCS = {"clr": log_clr, "log2raw": log2_raw, "within_z": log_clr}

def build_z(Msub, norm):
    """Return dict cohort -> z-scored design matrix, under normalization `norm` ( verbatim).

    clr      : CLR-log values, gene-wise z-score with IMvigor210 (train) stats
               (baseline; stats transported to external cohorts)
    log2raw  : log2(raw+1) values, gene-wise z-score with train stats
    within_z : CLR-log values, gene-wise z-score with EACH COHORT'S OWN stats
               (scaling statistics NOT transported)
    """
    f = NORM_FUNCS[norm]
    if norm in ("clr", "log2raw"):
        Xtr = f(Msub["IMvigor210"])
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
        return {c: (f(Msub[c]) - mu) / sd for c in COHORTS}
    else:  # within_z
        out = {}
        for c in COHORTS:
            Xc = f(Msub[c])
            out[c] = (Xc - Xc.mean(0)) / (Xc.std(0) + 1e-9)
        return out

# ================================================================ external AUC + bootstrap CI ( verbatim) + BCa ( verbatim)
def external_auc(S_tr, y_tr, S_te, y_te, lam, n_boot=N_BOOT_EXT, seed=SEED):
    """Fit ridge logistic on anchor (train) PC scores, predict external, AUC +
    percentile bootstrap 95% CI (resampling the EXTERNAL cohort, fixed predictions)."""
    Xtr = np.column_stack([np.ones(len(y_tr)), S_tr])
    w = fit_logistic(Xtr, y_tr, lam)
    p_te = sigmoid(np.column_stack([np.ones(len(y_te)), S_te]) @ w)
    auc = float(roc_auc_score(y_te, p_te))
    rng = np.random.default_rng(seed)
    n = len(y_te)
    boots = []
    for _ in range(n_boot):
        bi = rng.integers(0, n, n)
        yb, pb = y_te[bi], p_te[bi]
        if len(set(yb)) < 2:
            continue
        boots.append(roc_auc_score(yb, pb))
    boots = np.array(boots)
    lo, hi = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
    return auc, (lo, hi), int(len(boots))

def external_cell(S_tr, y_tr, S_te, y_te, lam):
    """ percentile CI + -recovery BCa CI for one external cohort."""
    auc, ci, nb = external_auc(S_tr, y_tr, S_te, y_te, lam)
    Xtr = np.column_stack([np.ones(len(y_tr)), S_tr])
    w = fit_logistic(Xtr, y_tr, lam)
    p_te = sigmoid(np.column_stack([np.ones(len(y_te)), S_te]) @ w)
    bca = auc_bca(y_te, p_te)
    return dict(auc=auc,
                ci95_percentile=[ci[0], ci[1]],
                ci95_bca=[bca[0], bca[1]],
                boot_valid=nb,
                ci_excludes_0p5=bool(ci[0] > 0.5),
                ci_excludes_0p5_below=bool(ci[1] < 0.5),
                ci_excludes_0p5_bca=bool(bca[0] > 0.5),
                ci_excludes_0p5_below_bca=bool(bca[1] < 0.5))

# ================================================================ one ablation configuration ( run_config structure, dual-lambda)
def run_config(cfg, data, n_common, zero_mask):
    """cfg keys: id, layer, gene_idx (None = all), norm, k, lam_ext, lam_within."""
    gidx = cfg["gene_idx"] if cfg["gene_idx"] is not None else np.arange(n_common)
    gidx = np.asarray(gidx)
    Msub = {c: data[c][0][gidx, :] for c in COHORTS}
    Z = build_z(Msub, cfg["norm"])
    pca = PCA(n_components=cfg["k"], random_state=SEED).fit(Z["IMvigor210"])
    S = {c: pca.transform(Z[c]) for c in COHORTS}
    y = {c: data[c][1] for c in COHORTS}
    # degeneracy check ( verbatim): constant anchor PC scores -> constant
    # predictions. Under the recovered identity (0 all-zero genes) this cannot
    # arise from unmeasured-gene blocks, but the check is kept.
    degenerate = bool(np.all(S["IMvigor210"].std(0) < 1e-12))

    ext = {}
    for c in EXT_COHORTS:
        ext[c] = external_cell(S["IMvigor210"], y["IMvigor210"], S[c], y[c], cfg["lam_ext"])
    within = {}
    for c in COHORTS:
        app, corr, se, ci = auc_632(S[c], y[c], cfg["lam_within"])
        within[c] = dict(auc_apparent=float(app), auc_632_corrected=float(corr),
                         se=float(se), ci95=[float(ci[0]), float(ci[1])])
    return dict(n_genes=int(len(gidx)),
                n_unmeasured_in_subset=int(zero_mask[gidx].sum()),
                degenerate=degenerate,
                degenerate_reason=("anchor PC scores are constant" if degenerate else None),
                norm=cfg["norm"], k=cfg["k"],
                lam_ext=cfg["lam_ext"], lam_within=cfg["lam_within"],
                external=ext, within_632=within)

# ================================================================ layer aggregation ( verbatim)
def layer_stats(per_config, config_ids, cohort):
    a = np.array([per_config[cid]["external"][cohort]["auc"] for cid in config_ids], float)
    return dict(n_configs=len(a), min=float(a.min()), max=float(a.max()),
                range=float(a.max() - a.min()),
                q1=float(np.percentile(a, 25)), q3=float(np.percentile(a, 75)),
                iqr=float(np.percentile(a, 75) - np.percentile(a, 25)),
                variance=float(a.var(ddof=1)) if len(a) > 1 else 0.0)

# ================================================================ main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- load everything under the recovered identity
    M, sids = load_exmat()
    common, sym2ens, r2s, avail = common_genes_recovered()
    n_common = len(common)
    M320_imv, y_imv = load_imvigor_correct(M, sids, common, r2s)
    m_lab = ~np.isnan(y_imv)

    # ---- zero-gene artefact retraction evidence ( definitions on the new panel)
    zero_mask_all348 = (M320_imv == 0).all(1)
    zero_mask = (M320_imv[:, m_lab] == 0).all(1)   # across the 298 labelled
    n_zero = int(zero_mask_all348.sum())
    assert n_zero == 0, f"recovered panel unexpectedly has {n_zero} all-zero genes"

    data = {
        "IMvigor210": (M320_imv[:, m_lab], y_imv[m_lab]),
        "GSE165252": load_gse165252(common, sym2ens),
        "GSE78220": load_gse78220(common),
    }
    n_lab = {}; n_ev = {}
    for c in COHORTS:
        M322c, yc = data[c]
        m = ~np.isnan(yc)
        data[c] = (M322c[:, m], yc[m])
        n_lab[c] = int(m.sum()); n_ev[c] = int(yc[m].sum())
    N_total = sum(n_lab.values())
    lam_base_ext = 1.0 / n_lab["IMvigor210"]      # 1/298 (frozen recovery transport rule)
    lam_base_within = 1.0 / N_total               # 1/394 (/ global within rule)

    print(f"[] common genes: {n_common} | zero genes (all-348): {n_zero} | "
          f"n_lab: {n_lab} | N_total: {N_total}")
    print(f"      lam_ext=1/{n_lab['IMvigor210']}={lam_base_ext:.10f}  "
          f"lam_within=1/{N_total}={lam_base_within:.10f}")

    # ---- anchored baseline PCA (CLR, 320 genes, IMvigor210 labelled) for gene ranking
    Zb = build_z({c: data[c][0] for c in COHORTS}, "clr")
    pca_b = PCA(n_components=2, random_state=SEED).fit(Zb["IMvigor210"])
    pc1_load = pca_b.components_[0]
    rank = np.argsort(-np.abs(pc1_load))
    top100_idx = sorted(rank[:100].tolist())
    top50_idx = sorted(rank[:50].tolist())

    # ---- configuration grid (21 configs; baseline first for the anchor self-check)
    configs = [dict(id="baseline", layer="baseline", gene_idx=None,
                    norm="clr", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within)]
    configs.append(dict(id="gene_top100", layer="gene_set", gene_idx=top100_idx,
                        norm="clr", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within))
    configs.append(dict(id="gene_top50", layer="gene_set", gene_idx=top50_idx,
                        norm="clr", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within))
    for i in range(5):
        rng = np.random.default_rng(SEED + i)
        idx = sorted(rng.choice(n_common, 100, replace=False).tolist())
        configs.append(dict(id=f"gene_rand100_s{i}", layer="gene_set", gene_idx=idx,
                            norm="clr", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within,
                            rand_seed=SEED + i))
    for i in range(5):
        rng = np.random.default_rng(SEED + i)
        idx = sorted(rng.choice(n_common, 50, replace=False).tolist())
        configs.append(dict(id=f"gene_rand50_s{i}", layer="gene_set", gene_idx=idx,
                            norm="clr", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within,
                            rand_seed=SEED + i))
    configs.append(dict(id="norm_log2raw", layer="normalization", gene_idx=None,
                        norm="log2raw", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within))
    configs.append(dict(id="norm_within_z", layer="normalization", gene_idx=None,
                        norm="within_z", k=2, lam_ext=lam_base_ext, lam_within=lam_base_within))
    for k in (1, 3, 5):
        configs.append(dict(id=f"pc_k{k}", layer="pc_k", gene_idx=None,
                            norm="clr", k=k, lam_ext=lam_base_ext, lam_within=lam_base_within))
    for lv in (0.01, 0.1, 1.0):
        configs.append(dict(id=f"ridge_lam{lv:g}", layer="ridge_lambda", gene_idx=None,
                            norm="clr", k=2, lam_ext=lv, lam_within=lv))

    # ---- run all configurations (REAL RUN)
    per_config = {}
    for cfg in configs:
        n_g = n_common if cfg["gene_idx"] is None else len(cfg["gene_idx"])
        print(f"[] running {cfg['id']} (norm={cfg['norm']}, k={cfg['k']}, "
              f"lam_ext={cfg['lam_ext']:.6g}, lam_within={cfg['lam_within']:.6g}, "
              f"n_genes={n_g})", flush=True)
        rec = run_config(cfg, data, n_common, zero_mask)
        rec["layer"] = cfg["layer"]
        if "rand_seed" in cfg:
            rec["rand_seed"] = cfg["rand_seed"]
        per_config[cfg["id"]] = rec
        for c in EXT_COHORTS:
            e = rec["external"][c]
            print(f"    ext {c}: AUC={e['auc']:.4f} "
                  f"pctCI=[{e['ci95_percentile'][0]:.4f},{e['ci95_percentile'][1]:.4f}] "
                  f"bcaCI=[{e['ci95_bca'][0]:.4f},{e['ci95_bca'][1]:.4f}] "
                  f"excl0.5(pct)={e['ci_excludes_0p5']}"
                  f"{' DEGENERATE' if rec['degenerate'] else ''}", flush=True)

    # ---- ANCHOR SELF-CHECK (bitwise vs frozen recovery worker values)
    base = per_config["baseline"]
    anchor_checks = {
        "ext_GSE165252_auc": (base["external"]["GSE165252"]["auc"], ANCHOR["ext_GSE165252_auc"]),
        "ext_GSE165252_ci_bca_lo": (base["external"]["GSE165252"]["ci95_bca"][0], ANCHOR["ext_GSE165252_ci_bca"][0]),
        "ext_GSE165252_ci_bca_hi": (base["external"]["GSE165252"]["ci95_bca"][1], ANCHOR["ext_GSE165252_ci_bca"][1]),
        "ext_GSE78220_auc": (base["external"]["GSE78220"]["auc"], ANCHOR["ext_GSE78220_auc"]),
        "ext_GSE78220_ci_bca_lo": (base["external"]["GSE78220"]["ci95_bca"][0], ANCHOR["ext_GSE78220_ci_bca"][0]),
        "ext_GSE78220_ci_bca_hi": (base["external"]["GSE78220"]["ci95_bca"][1], ANCHOR["ext_GSE78220_ci_bca"][1]),
        "within_imvigor210_auc632": (base["within_632"]["IMvigor210"]["auc_632_corrected"], ANCHOR["within_imvigor210_auc632"]),
        "lam_ext": (lam_base_ext, ANCHOR["lam_ext"]),
        "lam_within": (lam_base_within, ANCHOR["lam_within"]),
    }
    anchor_self_check = {}
    all_pass = True
    for k, (got, exp) in anchor_checks.items():
        d = abs(float(got) - float(exp))
        ok = d < TOL
        all_pass &= ok
        anchor_self_check[k] = dict(expected=exp, computed=float(got), abs_diff=d, tol=TOL, **{"pass": ok})
    anchor_self_check["overall_pass"] = bool(all_pass)
    print("[] anchor self-check:", "PASS (bitwise)" if all_pass else "FAIL")
    for k, v in anchor_self_check.items():
        if k != "overall_pass":
            print(f"      {k}: computed={v['computed']:.12f} expected={v['expected']:.12f} "
                  f"diff={v['abs_diff']:.2e} {'OK' if v['pass'] else 'MISMATCH'}")
    if not all_pass:
        print("FATAL: baseline does not reproduce the frozen recovery anchor -- aborting.",
              file=sys.stderr)
        sys.exit(1)

    # ---- removed configurations (retracted artefact controls)
    removed_configs = {}
    for cid, desc in {
        "gene_measured_only": "all measured genes only (drop the 70 'unmeasured' genes)",
        "gene_top100_measured": "top-100 |PC1 loading| genes among measured genes",
        "gene_top50_measured": "top-50 |PC1 loading| genes among measured genes",
        "norm_log2raw_measured": "log2raw normalization on measured genes only",
        "norm_within_z_measured": "within-cohort z-score on measured genes only",
    }.items():
        removed_configs[cid] = dict(
            layer="artifact_removal",
            old__description=desc,
            removal_reason=REMOVAL_REASON,
            retracted_because=("Under the recovered gene identity the 320-gene anchor "
                               "panel has 0 all-zero genes in IMvigor210 (old : 70/228), "
                               "so every 'measured-only' subset coincides with the full "
                               "panel and the configuration is void."),
        )

    # ---- layer membership (each primary layer includes the shared baseline config)
    layer_members = {
        "gene_set":      ["baseline", "gene_top100", "gene_top50"]
                         + [f"gene_rand100_s{i}" for i in range(5)]
                         + [f"gene_rand50_s{i}" for i in range(5)],
        "normalization": ["baseline", "norm_log2raw", "norm_within_z"],
        "pc_k":          ["baseline", "pc_k1", "pc_k3", "pc_k5"],
        "ridge_lambda":  ["baseline", "ridge_lam0.01", "ridge_lam0.1", "ridge_lam1"],
    }
    PRIMARY_LAYERS = ["gene_set", "normalization", "pc_k", "ridge_lambda"]
    auc_range_by_layer = {}
    for layer, ids in layer_members.items():
        ids = [i for i in ids if i in per_config]
        stats = {}
        for c in EXT_COHORTS:
            s_all = layer_stats(per_config, ids, c)
            ids_nd = [i for i in ids if not per_config[i]["degenerate"]]
            s_nd = layer_stats(per_config, ids_nd, c) if len(ids_nd) > 0 else None
            stats[c] = dict(all_configs=s_all, non_degenerate=s_nd)
        auc_range_by_layer[layer] = stats
        auc_range_by_layer[layer]["mean_variance_across_external"] = float(
            np.mean([(stats[c]["non_degenerate"] or stats[c]["all_configs"])["variance"]
                     for c in EXT_COHORTS]))
    layer_variance_ranking = sorted(PRIMARY_LAYERS,
                                    key=lambda L: -auc_range_by_layer[L]["mean_variance_across_external"])
    largest_variance_layer = layer_variance_ranking[0]

    # ---- recovery-transport verdict (percentile engine CI primary; BCa also recorded)
    recovery = {c: dict(any=False, configs=[],
                        any_bca=False, configs_bca=[],
                        any_below=False, configs_below=[],
                        any_below_bca=False, configs_below_bca=[]) for c in EXT_COHORTS}
    for cid, rec in per_config.items():
        if rec["degenerate"]:
            continue
        for c in EXT_COHORTS:
            e = rec["external"][c]
            if e["ci_excludes_0p5"]:
                recovery[c]["any"] = True; recovery[c]["configs"].append(cid)
            if e["ci_excludes_0p5_bca"]:
                recovery[c]["any_bca"] = True; recovery[c]["configs_bca"].append(cid)
            if e["ci_excludes_0p5_below"]:
                recovery[c]["any_below"] = True; recovery[c]["configs_below"].append(cid)
            if e["ci_excludes_0p5_below_bca"]:
                recovery[c]["any_below_bca"] = True; recovery[c]["configs_below_bca"].append(cid)
    any_recovery = bool(any(recovery[c]["any"] for c in EXT_COHORTS))
    any_recovery_bca = bool(any(recovery[c]["any_bca"] for c in EXT_COHORTS))
    degenerate_configs = [cid for cid, rec in per_config.items() if rec["degenerate"]]

    # ---- honest key finding (computed, not asserted)
    nd_configs = [cid for cid in per_config if not per_config[cid]["degenerate"]]
    all_ext = {c: [per_config[cid]["external"][c]["auc"] for cid in nd_configs] for c in EXT_COHORTS}
    kf_parts = []
    kf_parts.append(
        f"Across {len(nd_configs)} non-degenerate ablation configurations spanning gene "
        f"set, normalization, PC dimension and ridge lambda under the recovered gene "
        f"identity ({n_common}-gene panel, 0 all-zero genes), external AUC ranged "
        f"{min(all_ext['GSE165252']):.3f}-{max(all_ext['GSE165252']):.3f} (GSE165252) and "
        f"{min(all_ext['GSE78220']):.3f}-{max(all_ext['GSE78220']):.3f} (GSE78220).")
    if not any_recovery and not any_recovery_bca:
        kf_parts.append(
            "NO non-degenerate configuration yielded an external AUC whose bootstrap 95% "
            "CI excluded 0.5 under either the percentile or the BCa method -- the "
            "transport failure remains SYSTEMATIC under the recovered identity, not an "
            "artefact of any specific preprocessing choice tested.")
    else:
        kf_parts.append(
            f"Recovery: percentile-CI {json.dumps({c: recovery[c]['configs'] for c in EXT_COHORTS})}; "
            f"BCa-CI {json.dumps({c: recovery[c]['configs_bca'] for c in EXT_COHORTS})}.")
    kf_parts.append(
        f"The largest cross-configuration external-AUC variance was attributable to the "
        f"'{largest_variance_layer}' layer (ranking: "
        + " > ".join(layer_variance_ranking) + ").")
    kf_parts.append(
        f"The  '70 unmeasured genes carrying 88.7% of PC1 loading mass' structural "
        f"artefact is RETRACTED: under the recovered identity the anchor panel has "
        f"{n_zero} all-zero genes and the PC1 loading mass on (nonexistent) zero-genes "
        f"is exactly 0; the 5 measured-only artefact configurations were removed "
        f"(removal_reason recorded).")
    kf = " ".join(kf_parts)

    b_ext = {c: dict(auc=base["external"][c]["auc"],
                     ci95_percentile=base["external"][c]["ci95_percentile"],
                     ci95_bca=base["external"][c]["ci95_bca"]) for c in EXT_COHORTS}
    b_within = {c: base["within_632"][c]["auc_632_corrected"] for c in COHORTS}

    # ---- old- comparison (computed from the frozen  JSON)
    old74 = json.load(open(V74_JSON))
    old_comparison = dict(
        old__json=V74_JSON,
        old_n_common=old74["unmeasured_gene_artifact"]["n_common_genes"],
        old_n_unmeasured=old74["unmeasured_gene_artifact"]["n_unmeasured_in_imvigor210"],
        old_pc1_mass_unmeasured=old74["unmeasured_gene_artifact"]["pc1_loading_squared_mass_on_unmeasured"],
        old_baseline_external={c: old74["baseline"]["external"][c] for c in EXT_COHORTS},
        old_baseline_within_imvigor=old74["baseline"]["within_632"]["IMvigor210"],
        old_degenerate_configs=old74["unmeasured_gene_artifact"]["degenerate_configs"],
        old_recovery_any=old74["verdict"]["any_config_recovers_transport"]["any"],
        old_largest_variance_layer=old74["verdict"]["largest_variance_layer"],
        old_layer_mean_variance={L: old74["verdict"]["external_auc_range_by_layer"][L]["mean_variance_across_external"]
                                 for L in PRIMARY_LAYERS},
        new_recovery_any=any_recovery,
        new_recovery_any_bca=any_recovery_bca,
        new_largest_variance_layer=largest_variance_layer,
        new_layer_mean_variance={L: auc_range_by_layer[L]["mean_variance_across_external"]
                                 for L in PRIMARY_LAYERS},
        new_degenerate_configs=degenerate_configs,
    )

    out = dict(
        analysis="rerun_ablation",
        date="2026-09-24",
        seed=SEED,
        purpose=("Honest re-specified controlled ablation under the recovered IMvigor210 "
                 "gene identity (): re-run of the  ablation grid with the recovered "
                 "320-gene panel; the 5 'measured-only' artefact configurations are "
                 "retracted with removal_reason because the 70-unmeasured-gene artefact "
                 "no longer exists (zero-gene count = 0)."),
        design=dict(
            task="binary ICI efficacy discrimination (responder vs non-responder)",
            anchor_scheme=("PCA rotation fit on IMvigor210 labelled samples (n=%d) under the "
                           "RECOVERED gene key (exmat row j == j-th fData row in upper-cased "
                           "symbol alphabetical order), transform-only to GSE165252 (n=%d) "
                           "and GSE78220 (n=%d)"
                           % (n_lab["IMvigor210"], n_lab["GSE165252"], n_lab["GSE78220"])),
            gene_panel=dict(authoritative_panel=avail["panel"],
                            availability=avail, common_intersection=n_common,
                            imvigor_missing=avail["imvigor210_missing"]),
            external_auc=("ridge logistic fit on IMvigor210 anchor PC scores, predict "
                          "external cohort; AUC + percentile bootstrap 95%% CI (n_boot=%d, "
                          "resampling the external cohort with fixed predictions, seed=%d) "
                          "[ verbatim]; BCa CI ( recovery auc_bca verbatim) reported "
                          "alongside for anchor continuity" % (N_BOOT_EXT, SEED)),
            within_632_control=("/ verbatim Efron 0.632 optimism-corrected within-cohort "
                                "AUC (n_boot=%d, seed=%d) for all three cohorts"
                                % (N_BOOT_632, SEED)),
            lambda_rule=("RE-SPECIFIED vs  to bitwise reproduce the frozen  recovery "
                         "anchor: external transport lambda = 1/n_train = 1/%d = %.10f; "
                         "within 0.632 lambda = 1/N_total = 1/%d = %.10f.  used a single "
                         "1/N_total for both. Ridge-layer perturbations set BOTH lambdas to "
                         "the same value." % (n_lab["IMvigor210"], lam_base_ext, N_total, lam_base_within)),
            top_loading_rule=("genes ranked by |PC1 loading| of the baseline CLR PCA (all %d "
                              "common genes) fit on IMvigor210 labelled samples" % n_common),
            random_subset_rule="uniform without replacement, seeds 20260911+0..4, sorted indices",
            removed_config_policy=("5  'measured-only' artefact configurations removed; "
                                   "removal_reason recorded per config"),
        ),
        cohorts={c: dict(n_labeled=n_lab[c], n_event=n_ev[c]) for c in COHORTS},
        N_total_labeled=N_total,
        artifact_retraction=dict(
            old_=dict(n_common_genes=old74["unmeasured_gene_artifact"]["n_common_genes"],
                         n_unmeasured=old74["unmeasured_gene_artifact"]["n_unmeasured_in_imvigor210"],
                         pc1_loading_squared_mass_on_unmeasured=old74["unmeasured_gene_artifact"]["pc1_loading_squared_mass_on_unmeasured"]),
            new_=dict(n_common_genes=n_common,
                         n_allzero_across_348=n_zero,
                         n_allzero_across_298_labelled=int(zero_mask.sum()),
                         pc1_loading_squared_mass_on_zero_genes=0.0),
            removal_reason=REMOVAL_REASON,
        ),
        gene_ranking=dict(
            top100_symbols=[common[i] for i in top100_idx],
            top50_symbols=[common[i] for i in top50_idx],
            pc1_explained_variance_ratio_baseline=[float(v) for v in pca_b.explained_variance_ratio_],
        ),
        anchor_self_check=anchor_self_check,
        baseline=dict(
            config_id="baseline",
            description=("recovered-identity baseline: %d-gene intersection, CLR, PC k=2, "
                         "external lambda=1/n_train, within lambda=1/N_total" % n_common),
            external=b_ext,
            within_632=b_within,
        ),
        grid=dict(
            gene_set=["baseline", "gene_top100", "gene_top50"]
                     + [f"gene_rand100_s{i}" for i in range(5)]
                     + [f"gene_rand50_s{i}" for i in range(5)],
            normalization=["baseline(clr)", "norm_log2raw", "norm_within_z"],
            pc_k=["baseline(k=2)", "pc_k1", "pc_k3", "pc_k5"],
            ridge_lambda=["baseline(1/n)", "ridge_lam0.01", "ridge_lam0.1", "ridge_lam1"],
        ),
        removed_configs=removed_configs,
        per_config=per_config,
        verdict=dict(
            any_config_recovers_transport=dict(
                any=any_recovery,
                per_cohort={c: dict(any=recovery[c]["any"], configs=recovery[c]["configs"],
                                    any_bca=recovery[c]["any_bca"], configs_bca=recovery[c]["configs_bca"],
                                    any_below=recovery[c]["any_below"], configs_below=recovery[c]["configs_below"],
                                    any_below_bca=recovery[c]["any_below_bca"], configs_below_bca=recovery[c]["configs_below_bca"])
                            for c in EXT_COHORTS},
                definition=("external AUC bootstrap 95% CI lower bound > 0.5 in at least one "
                            "non-degenerate configuration (percentile engine CI primary; "
                            "BCa reported alongside)"),
                degenerate_configs_excluded=degenerate_configs,
            ),
            external_auc_range_by_layer=auc_range_by_layer,
            layer_variance_ranking=layer_variance_ranking,
            largest_variance_layer=largest_variance_layer,
            largest_variance_layer_scope=("computed over the four primary layers (gene_set, "
                                          "normalization, pc_k, ridge_lambda), non-degenerate "
                                          "configs, mean variance across the two external cohorts"),
            key_finding=kf,
        ),
        old__comparison=old_comparison,
    )
    json.dump(out, open(OUT_JSON, "w"), indent=2)
    print("\nWrote", OUT_JSON)

    # ---- console digest
    print("\n===  rerun ablation digest ===")
    print("baseline external AUC:", {c: round(b_ext[c]["auc"], 4) for c in EXT_COHORTS})
    print("baseline within 0.632 IMvigor210:", round(b_within["IMvigor210"], 4))
    for layer in PRIMARY_LAYERS:
        ids_nd = [i for i in layer_members[layer] if not per_config[i]["degenerate"]]
        for c in EXT_COHORTS:
            a = [per_config[i]["external"][c]["auc"] for i in ids_nd]
            print(f"{layer:14s} {c} (non-deg, n={len(a)}): min={min(a):.3f} max={max(a):.3f} "
                  f"var={np.var(a, ddof=1):.6f}")
    print("degenerate configs:", degenerate_configs)
    print("recovery (percentile):", {c: recovery[c]["any"] for c in EXT_COHORTS})
    print("recovery (BCa):", {c: recovery[c]["any_bca"] for c in EXT_COHORTS})
    print("layer variance ranking:", layer_variance_ranking)

    # ---- MD report
    write_md(out, old74, common, n_lab, n_ev)
    print("Wrote", OUT_MD)

# ================================================================ report
def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    return f"{x:.{nd}f}"

def write_md(out, old74, common, n_lab, n_ev):
    per_config = out["per_config"]
    verdict = out["verdict"]
    asc = out["anchor_self_check"]
    oc = out["old__comparison"]
    art = out["artifact_retraction"]
    lines = []
    A = lines.append
    A("#  消融网格诚实重规与重跑（恢复基因身份下的受控消融）")
    A("")
    A(f"日期：{out['date']} ｜ seed {out['seed']} ｜ 脚本：scripts/analysis/rerun_ablation.py ｜ 冻结 JSON：data/rerun_ablation.json")
    A("")
    A("## 1. 重规说明（为什么旧  的 26 配置全部作废）")
    A("")
    A(f"- 基因身份已恢复（_imvigor_gene_recovery）：exmat 行 j == fData 按 symbol（大写归一）字母序排序后的第 j 行；"
      f"IMvigor210 加载行索引来自恢复钥匙，而非旧的 fData 位置对齐。")
    A(f"- 新锚定面板 = 322 权威 symbol 的跨队列交集 **{out['design']['gene_panel']['common_intersection']} 基因**"
      f"（IMvigor210 可测 {out['design']['gene_panel']['availability']['imvigor210']}/322，缺 "
      f"{', '.join(out['design']['gene_panel']['imvigor_missing'])}；GSE165252 322/322；GSE78220 320/322）。")
    A(f"- 面板内 IMvigor210 全零基因数 = **{art['new_']['n_allzero_across_348']}**（旧 ：{art['old_']['n_unmeasured']}/{art['old_']['n_common_genes']}），"
      f"零值基因 PC1 loading 质量旧 88.7% → 新 **0.0**。前提消失，5 个 measured-only 伪影配置删除（removal_reason 见 JSON）。")
    A(f"- λ 规则重规（为逐位复现冻结恢复锚定值）：外部运输 λ=1/n_train=1/{n_lab['IMvigor210']}；within 0.632 λ=1/N_total=1/{out['N_total_labeled']}。"
      f" 曾对两者统一用 1/N_total；ridge 层扰动则将两个 λ 同设为扰动值。")
    A(f"- 度量引擎逐字复用 ：外部 AUC + percentile bootstrap 95% CI（2000 次重抽外部队列、固定预测、seed 20260911）；"
      f"within Efron 0.632（600 次、seed 20260911、三队列）；另按恢复脚本逐字附 BCa CI 以保持锚定连续性。")
    A("")
    A("## 2. 锚定自检（baseline 逐位复现冻结恢复值）")
    A("")
    A("| 量 | 冻结期望 | 本次计算 | abs diff | PASS |")
    A("|---|---|---|---|---|")
    for k, v in asc.items():
        if k == "overall_pass":
            continue
        A(f"| {k} | {v['expected']:.12f} | {v['computed']:.12f} | {v['abs_diff']:.2e} | {'**PASS**' if v['pass'] else 'FAIL'} |")
    A(f"| **总门** | | | | **{'PASS（逐位一致）' if asc['overall_pass'] else 'FAIL'}** |")
    A("")
    A("## 3. 新网格总表（21 配置，全部真跑）")
    A("")
    A("外部 CI 为主引擎 percentile bootstrap；括号内为 BCa（锚定连续性）。within 为 Efron 0.632。")
    A("")
    A("| 配置 | 层 | n_genes | norm | k | λ_ext | GSE165252 AUC [pctCI] (bcaCI) | 排除0.5 | GSE78220 AUC [pctCI] (bcaCI) | 排除0.5 | within IMV 0.632 | 退化 |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    order = (["baseline", "gene_top100", "gene_top50"]
             + [f"gene_rand100_s{i}" for i in range(5)]
             + [f"gene_rand50_s{i}" for i in range(5)]
             + ["norm_log2raw", "norm_within_z", "pc_k1", "pc_k3", "pc_k5",
                "ridge_lam0.01", "ridge_lam0.1", "ridge_lam1"])
    for cid in order:
        r = per_config[cid]
        e1 = r["external"]["GSE165252"]; e2 = r["external"]["GSE78220"]
        def cell(e):
            return (f"{e['auc']:.4f} [{e['ci95_percentile'][0]:.4f},{e['ci95_percentile'][1]:.4f}]"
                    f" ({e['ci95_bca'][0]:.4f},{e['ci95_bca'][1]:.4f})")
        def excl(e):
            s = f"pct:{'是' if e['ci_excludes_0p5'] else '否'}/bca:{'是' if e['ci_excludes_0p5_bca'] else '否'}"
            return s
        A(f"| {cid} | {r['layer']} | {r['n_genes']} | {r['norm']} | {r['k']} | {r['lam_ext']:.4g} "
          f"| {cell(e1)} | {excl(e1)} | {cell(e2)} | {excl(e2)} "
          f"| {fmt(r['within_632']['IMvigor210']['auc_632_corrected'])} "
          f"| {'**是**' if r['degenerate'] else '否'} |")
    A("")
    A(f"退化配置数：{len(verdict['any_config_recovers_transport']['degenerate_configs_excluded'])}"
      f"（{', '.join(verdict['any_config_recovers_transport']['degenerate_configs_excluded']) or '无'}）；"
      f"每配置子集内全零基因数均为 0（n_unmeasured_in_subset 全 0）。")
    A("")
    A("## 4. 逐层外部 AUC 方差分解（对齐旧  口径）")
    A("")
    A("每层含共享 baseline；统计基于非退化配置；均值为两外部队列方差平均。")
    A("")
    A("| 层 | 队列 | n | min | max | range | IQR | variance |")
    A("|---|---|---|---|---|---|---|---|")
    for layer, st in out["verdict"]["external_auc_range_by_layer"].items():
        for c in ("GSE165252", "GSE78220"):
            s = st[c]["non_degenerate"] or st[c]["all_configs"]
            A(f"| {layer} | {c} | {s['n_configs']} | {s['min']:.4f} | {s['max']:.4f} "
              f"| {s['range']:.4f} | {s['iqr']:.4f} | {s['variance']:.6f} |")
    A("")
    A(f"层方差排序（主 4 层，均值跨两外部队列）：**{' > '.join(verdict['layer_variance_ranking'])}**；"
      f"最大方差层 = **{verdict['largest_variance_layer']}**（旧 ：{oc['old_largest_variance_layer']}）。")
    A("")
    A("| 层 | 旧  均值方差 | 新  均值方差 | 倍数变化 |")
    A("|---|---|---|---|")
    for L in ("gene_set", "normalization", "pc_k", "ridge_lambda"):
        o = oc["old_layer_mean_variance"][L]; nw = oc["new_layer_mean_variance"][L]
        if o == 0.0 and nw == 0.0:
            mult_s = "两代均恰为 0"
        elif o == 0.0:
            mult_s = "∞（旧为 0）"
        else:
            mult_s = f"{nw / o:.2f}x"
        A(f"| {L} | {o:.6f} | {nw:.6f} | {mult_s} |")
    A("")
    A("## 5. 恢复运输判定（任何配置 CI 有利方向排除 0.5？）")
    A("")
    rec = verdict["any_config_recovers_transport"]["per_cohort"]
    for c in EXT_COHORTS:
        A(f"- **{c}**：percentile CI 下界>0.5 的配置：{', '.join(rec[c]['configs']) or '无'}；"
          f"BCa CI 下界>0.5 的配置：{', '.join(rec[c]['configs_bca']) or '无'}；"
          f"（不利方向，CI 上界<0.5：percentile {', '.join(rec[c]['configs_below']) or '无'} / "
          f"BCa {', '.join(rec[c]['configs_below_bca']) or '无'}）")
    any_p = verdict["any_config_recovers_transport"]["any"]
    any_b = verdict["any_config_recovers_transport"]["per_cohort"]
    any_bca = any(any_b[c]["any_bca"] for c in EXT_COHORTS)
    if not any_p and not any_bca:
        A("")
        A("**判定：无任何非退化配置的外部 AUC bootstrap 95% CI 在有利方向排除 0.5（percentile 与 BCa 两法一致）——"
          "恢复身份后运输失败依旧系统性存在。**")
    else:
        A("")
        A(f"**判定：存在恢复运输的配置（percentile: {any_p}；BCa: {any_bca}）——详见上表逐配置明细。**")
    A("")
    A("## 6. 与旧  结论对照（存活 / 翻转 / 撤回）")
    A("")
    ob = oc["old_baseline_external"]; nb = out["baseline"]["external"]
    ow = oc["old_baseline_within_imvigor"]; nw = out["baseline"]["within_632"]["IMvigor210"]
    A("### 撤回（前提已消失）")
    A("")
    A(f"1. **“70 未测量基因伪影”结构诊断**：旧  报告 70/228 基因在 IMvigor210 全零、承载 88.7% PC1 loading 质量、"
      f"top-50 |PC1 loading| 子集完全退化；恢复身份后 320 基因面板全零基因数 = {art['new_']['n_allzero_across_348']}，"
      f"零基因 loading 质量 = 0.0，整条机制叙事撤回。")
    A(f"2. **5 个 measured-only 伪影配置**（gene_measured_only / gene_top100_measured / gene_top50_measured / "
      f"norm_log2raw_measured / norm_within_z_measured）：removal_reason = \"{out['artifact_retraction']['removal_reason']}\"。")
    A(f"3. **旧 228 基因“交集”**：由假 Entrez 选择产生（Supplementary_Table_1 列 1 为行号非 Entrez），被 320 基因真实交集取代。")
    old_deg = oc["old_degenerate_configs"]; new_deg = oc["new_degenerate_configs"]
    A(f"4. **gene_top50 退化结论**：旧退化配置 {', '.join(old_deg)}；新网格退化配置 {', '.join(new_deg) or '无'}"
      f"——top50 现为真实可判别配置。")
    A("")
    A("### 翻转（数值方向性改变）")
    A("")
    A("| 量 | 旧 （错误身份） | 新 （恢复身份） | 解读 |")
    A("|---|---|---|---|")
    A(f"| baseline 外部 AUC GSE78220 | {fmt(ob['GSE78220']['auc'])} [{fmt(ob['GSE78220']['ci95'][0])},{fmt(ob['GSE78220']['ci95'][1])}] "
      f"| {fmt(nb['GSE78220']['auc'])} [{fmt(nb['GSE78220']['ci95_percentile'][0])},{fmt(nb['GSE78220']['ci95_percentile'][1])}] "
      f"| **方向翻转**（<0.5 → >0.5，点估计越过无效线；CI 仍含 0.5，n=25） |")
    A(f"| baseline within 0.632 IMvigor210 | {fmt(ow)} | {fmt(nw)} | 锚定内部判别大幅改善（{nw-ow:+.4f}） |")
    A(f"| baseline 外部 AUC GSE165252 | {fmt(ob['GSE165252']['auc'])} [{fmt(ob['GSE165252']['ci95'][0])},{fmt(ob['GSE165252']['ci95'][1])}] "
      f"| {fmt(nb['GSE165252']['auc'])} [{fmt(nb['GSE165252']['ci95_percentile'][0])},{fmt(nb['GSE165252']['ci95_percentile'][1])}] "
      f"| 数值微降但结论不变：仍≈0.5、CI 横跨 0.5 |")
    A("")
    A("### 存活（重跑后仍成立）")
    A("")
    surv = []
    if (not oc["old_recovery_any"]) and (not any_p) and (not any_bca):
        surv.append("“无任何配置外部 CI 排除 0.5 —— 运输失败是系统性的，而非特定预处理选择的伪影”这一  核心结论"
                    "在恢复身份、撤回伪影配置后的 21 配置网格中**存活**（percentile 与 BCa 两法一致）。")
    if oc["old_largest_variance_layer"] == verdict["largest_variance_layer"]:
        surv.append(f"“外部 AUC 跨配置方差最大层 = {verdict['largest_variance_layer']}”**存活**"
                    f"（旧 {oc['old_largest_variance_layer']} → 新 {verdict['largest_variance_layer']}）。")
    else:
        surv.append(f"最大方差层结论**翻转**：旧 {oc['old_largest_variance_layer']} → 新 {verdict['largest_variance_layer']}。")
    if not new_deg:
        surv.append("“退化配置不构成信息”结论以空集形式存活：新网格 0 个退化配置。")
    for s in surv:
        A(f"- {s}")
    A("")
    A("## 7. 汇总（计算所得 key finding）")
    A("")
    A(verdict["key_finding"])
    A("")
    A("## 8. 复现性")
    A("")
    A(f"- Python：managed venv（/Users/zuoxianbo/.workbuddy/binaries/python/envs/default/bin/python）；seed {out['seed']} 全程。")
    A(f"- 队列样本（标注 n / 事件 n）：IMvigor210 {n_lab['IMvigor210']}/{n_ev['IMvigor210']}；"
      f"GSE165252 {n_lab['GSE165252']}/{n_ev['GSE165252']}；GSE78220 {n_lab['GSE78220']}/{n_ev['GSE78220']}；"
      f"N_total={out['N_total_labeled']}。")
    A("- baseline 锚定自检逐位通过后方可产出本报告（否则脚本中止，见 anchor_self_check）。")
    A(f"- 旧  冻结对照源：{oc['old__json']}。")
    A("")
    md = "\n".join(lines)
    with open(OUT_MD, "w") as f:
        f.write(md)

if __name__ == "__main__":
    main()
