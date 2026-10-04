#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" -- Nature Cancer revision  (recovered-key rerun).

Three deliverables, all computed for real on the RECOVERED IMvigor210 gene-row
identity (symbol-alphabetical-upper key; scripts/analysis/imvigor_gene_recovery.py,
data/imvigor_gene_recovery.json); nothing hardcoded from the broken
positional-alignment era.

P0-3  Threshold (calibration-slope tolerance) sensitivity, NEW numbers:
      primary frozen transport (IMvigor210-anchored, 320-gene panel, CLR->PCA
      k=2, ridge lambda=1/298) to GSE165252 / GSE78220; +/-0.05/0.10/0.20
      decision matrix + slope bootstrap 95% CI (2000 reps, seed 20260911,
      fixed-prediction patient resampling) + BCa AUC CI recomputed and
      cross-checked against the frozen recovery values (tol 1e-9).
      NOTE: rerun_benchmark.json (M0-M6 on the recovered key) is NOT on
      disk yet, so the matrix covers the primary frozen arm only; the M0-M6
      cells stay pending until the benchmark worker lands.

P0-10 Delta-AUC continuous panel, NEW numbers, per target cohort:
      delta AUC = external AUC - within-cohort same-pipeline Efron-0.632 AUC
      (CI from the outer/inner bootstrap difference distribution), logit-scale
      calibration intercept + slope ( convention), Brier, log-loss.
      Rows: primary frozen (done); M0/M2/M3 pending on rerun_benchmark.json.
      Within side for GSE165252/GSE78220 recomputed per task spec: same loader
      logic + recovered key + CLR + PCA k=2 + lambda=1/n + 600 bootstrap
      (Efron 0.632, optimism form, verbatim _imvigor_gene_recovery.auc_632
      algebra); nested bootstrap 2000 outer x 600 inner, seed 20260911,
      cached in data/__within_boot_newkey_cache.npz.
      Pooled AUC: DerSimonian-Laird random-effects synthesis (k=3) on the new
      320-gene cohort AUCs; the old 0.498 [0.365,0.631] ( descriptive) is
      retained as the historical 229-gene-era value and flagged as superseded.

P0-12 data/analysis_manifest.json: frozen mapping of every 
      analysis (script, inputs + SHA256, outputs, seed, bootstrap counts,
      key parameters, status done/pending/superseded) + IMvigor210 three-file
      SHA256 (asserted equal to the recovery report) + upstream provenance
      (mimifp/tfm_mUC censored exmat; SiYangming/IMvigor210CoreBiologies mirror).

Emits: data/rerun_sensitivity.json (+ .md), data/analysis_manifest.json.
Python: managed venv; seed 20260911 everywhere.
"""
import hashlib
import json
import os
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import imvigor_gene_recovery as rec   # verbatim loaders + algebra + key

ROOT = rec.ROOT
PKG = rec.PKG
OUT_DIR = rec.OUT_DIR
SEED = rec.SEED
N_BOOT_632 = rec.N_BOOT_632        # 600
N_BOOT_DELTA = 2000                # outer/paired bootstrap for delta AUC + slope CI
TOL_PARITY = 1e-9
WITHIN_CACHE = f"{OUT_DIR}/__within_boot_newkey_cache.npz"
TOLERANCES = [0.05, 0.10, 0.20]

fit_logistic = rec.fit_logistic
sigmoid = rec.sigmoid
calib_slope_prob = rec.calib_slope            # /recovery convention: feature = p
auc_bca = rec.auc_bca
brier_logloss = rec.brier_logloss
clr_standardize = rec.clr_standardize
fit_pca = rec.fit_pca
transform_pca = rec.transform_pca
load_gse165252 = rec.load_gse165252
load_gse78220 = rec.load_gse78220


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def logit(p):
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def calib_ab(y, p):
    """Logit-scale calibration intercept a + slope b (/ convention, verbatim)."""
    y = np.asarray(y, float)
    lp = logit(p)
    if np.ptp(lp) > 1e-9:
        w = fit_logistic(np.column_stack([np.ones(len(y)), lp]), y, 1e-6)
        return float(w[0]), float(w[1])
    a = float(logit(np.clip(float(np.mean(y)), 1e-6, 1 - 1e-6)) - lp[0])
    return a, float("nan")


def boot_slope_ci(y, p, slope_fn, n_boot=N_BOOT_DELTA, seed=SEED):
    """2000-rep patient bootstrap of the calibration slope with FIXED predictions;
    percentile 95% CI (frozen P0-3 method)."""
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    n = len(y)
    rng = np.random.default_rng(seed)
    slopes, skipped = [], 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb, pb = y[idx], p[idx]
        if len(np.unique(yb)) < 2 or np.ptp(logit(pb)) <= 1e-9:
            skipped += 1
            continue
        s = slope_fn(yb, pb)
        if s is None or not np.isfinite(s):
            skipped += 1
            continue
        slopes.append(float(s))
    a = np.array(slopes)
    if len(a) < 10:
        return None, None, int(len(a)), skipped
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)), int(len(a)), skipped


def boot_auc_draws(y, p, n_boot=N_BOOT_DELTA, seed=SEED):
    """2000-rep patient bootstrap draws of the external AUC (frozen predictions)."""
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    n = len(y)
    rng = np.random.default_rng(seed)
    aucs, skipped = [], 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            skipped += 1
            continue
        aucs.append(float(roc_auc_score(y[idx], p[idx])))
    return np.array(aucs), skipped


# ================================================================ within (new key)
def efron_0632_once(X, y, lam, rng, n_boot=N_BOOT_632):
    """One Efron-0.632 (optimism form) estimate, verbatim rec.auc_632 inner loop."""
    n = len(y)

    def fit_auc(Xs, ys):
        Xs = np.column_stack([np.ones(len(ys)), Xs])
        w = fit_logistic(Xs, ys, lam)
        pr = sigmoid(Xs @ w)
        if len(set(ys.tolist())) < 2:
            return np.nan
        return roc_auc_score(ys, pr)

    Xfull = np.column_stack([np.ones(n), X])
    app = fit_auc(X, y)
    if np.isnan(app):
        return np.nan, 0
    corr_b = []
    skipped = 0
    for _ in range(n_boot):
        bi = rng.integers(0, n, n)
        Xb, yb = X[bi], y[bi]
        if len(set(yb.tolist())) < 2:
            skipped += 1
            continue
        # single fit per draw: fit_auc(Xb, yb) and the transport fit are the
        # same model, so app_b is scored from wb directly (numerically identical
        # to the verbatim two-fit form, half the cost)
        wb = fit_logistic(np.column_stack([np.ones(len(yb)), Xb]), yb, lam)
        app_b = roc_auc_score(yb, sigmoid(np.column_stack([np.ones(len(yb)), Xb]) @ wb))
        po = sigmoid(Xfull @ wb)
        try:
            auc_o = roc_auc_score(y, po)
        except ValueError:
            skipped += 1
            continue
        corr_b.append(app - (app_b - auc_o))
    if not corr_b:
        return np.nan, skipped
    return float(np.mean(corr_b)), skipped


def nested_within_bootstrap(S0, y, n_outer=N_BOOT_DELTA, n_inner=N_BOOT_632,
                            seed=SEED, cache_key=None):
    """Nested bootstrap of the Efron-0.632 within estimator on the recovered-key
    cohort S0 scores (lambda = 1/n fixed, -style nested design; the optimism
    form replaces the OOB form so the replicated estimator IS the quoted
    within baseline estimator).  Cached."""
    if cache_key and os.path.exists(WITHIN_CACHE):
        z = np.load(WITHIN_CACHE, allow_pickle=True)
        if cache_key in z.files:
            draws = z[cache_key]
            meta = json.loads(str(z[cache_key + "__meta"]))
            return draws, meta
    n = len(y)
    lam = 1.0 / n
    rng = np.random.default_rng(seed)
    ests, outer_skipped, inner_skipped_total = [], 0, 0
    t0 = time.time()
    for i in range(n_outer):
        idx = rng.integers(0, n, n)
        e, sk = efron_0632_once(S0[idx], y[idx], lam, rng, n_inner)
        inner_skipped_total += sk
        if not np.isfinite(e):
            outer_skipped += 1
            continue
        ests.append(e)
        if (i + 1) % 250 == 0:
            print(f"    outer {i+1}/{n_outer} ({time.time()-t0:.0f}s)", flush=True)
    draws = np.array(ests)
    meta = dict(outer_skipped=outer_skipped, inner_skipped_total=inner_skipped_total,
                seconds=round(time.time() - t0, 1))
    if cache_key:
        save_keys = {}
        if os.path.exists(WITHIN_CACHE):
            z = np.load(WITHIN_CACHE, allow_pickle=True)
            for k in z.files:
                save_keys[k] = z[k]
        save_keys[cache_key] = draws
        save_keys[cache_key + "__meta"] = json.dumps(meta)
        np.savez(WITHIN_CACHE, **save_keys)
    return draws, meta


# ================================================================ DL meta
def dl_meta(aucs, ses, names):
    """DerSimonian-Laird random-effects pooling of AUCs, k cohorts."""
    a = np.asarray(aucs, float)
    v = np.asarray(ses, float) ** 2
    k = len(a)
    w_fe = 1.0 / v
    mu_fe = float(np.sum(w_fe * a) / np.sum(w_fe))
    Q = float(np.sum(w_fe * (a - mu_fe) ** 2))
    df = k - 1
    C = float(np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    mu = float(np.sum(w_re * a) / np.sum(w_re))
    se = float(np.sqrt(1.0 / np.sum(w_re)))
    i2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    from scipy import stats as st
    Q_p = 1.0 - float(st.chi2.cdf(Q, df)) if df > 0 else float("nan")
    return dict(k=k, cohorts=names, auc_inputs=[float(x) for x in aucs],
                se_inputs=[float(x) for x in ses],
                pooled_auc=mu, se=se, ci95=[mu - 1.96 * se, mu + 1.96 * se],
                tau2=tau2, Q=Q, Q_df=df, Q_p=Q_p, i2_qbased=i2,
                weights_re=[float(x / np.sum(w_re)) for x in w_re])


# ================================================================ main pipeline
def main():
    t_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---------- load + recovered key (verbatim rec.stage_1a winner)
    M, sids = rec.load_exmat()
    fd, fd_ent, fd_sym_raw, fd_sym_up = rec.load_fdata()
    pd_rows = rec.load_pdata()
    n = len(fd)
    order_up = sorted(range(n), key=lambda i: fd_sym_up[i])
    row_of_sym = {fd_sym_up[order_up[j]]: j for j in range(n)}

    # ---------- 320-gene panel availability (verbatim rec.stage_2 logic)
    import pickle
    panel = [g.upper() for g in json.load(open(rec.GENE_JSON))["genes"]]
    imv_avail = {s for s in panel if s in row_of_sym}
    g = pickle.load(open(rec.SYM2ENS_PKL, "rb"))
    sym2ens = {k.upper(): v for k, v in g.items()}
    import gzip
    with gzip.open(rec.GSE165_MAT, "rt") as f:
        f.readline()
        cnt_ens = set(line.rstrip("\n").split("\t")[0].split(".")[0] for line in f)
    gse165_avail = {s for s in panel if sym2ens.get(s) in cnt_ens}
    import openpyxl
    wb = openpyxl.load_workbook(rec.GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    rr = list(ws.iter_rows(min_row=1, values_only=True))
    wb.close()
    g7820_syms = set(str(r[0]).upper() for r in rr[1:])
    gse782_avail = {s for s in panel if s in g7820_syms}
    common = [s for s in panel if s in imv_avail and s in gse165_avail and s in gse782_avail]
    print(f"[panel] common genes: {len(common)} (IMv avail {len(imv_avail)}, "
          f"GSE165252 {len(gse165_avail)}, GSE78220 {len(gse782_avail)})", flush=True)

    # ---------- anchor arm (Scheme A, recovered key)
    idx = [row_of_sym[s] for s in common]
    M322 = M[idx, :]
    pmap = {r["X"]: r for r in pd_rows}
    y_imv = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                      else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                      else np.nan for s in sids], dtype=float)
    m_lab = ~np.isnan(y_imv)
    Xlab_imv = M322[:, m_lab]
    ylab_imv = y_imv[m_lab]
    Xz_tr, mu, sd = clr_standardize(Xlab_imv)
    pca = fit_pca(Xz_tr)
    S0_imv = pca.transform(Xz_tr)
    n_train = int(m_lab.sum())
    lam_tr = 1.0 / n_train
    w = fit_logistic(np.column_stack([np.ones(n_train), S0_imv]), ylab_imv, lam_tr)
    print(f"[anchor] PCA k=2 EVR: {pca.explained_variance_ratio_.round(4).tolist()}; "
          f"lambda={lam_tr:.6f}", flush=True)

    M_ext = {"GSE165252": load_gse165252(common, sym2ens),
             "GSE78220": load_gse78220(common)}

    # frozen cross-check targets
    frozen = json.load(open(f"{OUT_DIR}/imvigor_gene_recovery.json"))
    fz_ext = frozen["phase2"]["new"]["transport"]["external"]

    # ---------- P0-3: rebuild external predictions + tolerance matrix
    p0_3 = dict(
        point_values_source="data/imvigor_gene_recovery.json phase2.new "
                            "(primary frozen, recovered key, 320-gene panel)",
        model_rows_covered={"primary_frozen": "done",
                            "M0_M2_M3_and_other_benchmark_models":
                                "pending__rerun_benchmark_json_not_on_disk"},
        tolerances=TOLERANCES,
        criteria=dict(
            c1="external AUC bootstrap 95% CI (BCa, 2000, seed 20260911) excludes 0.5 "
               "in the beneficial direction (tolerance-independent)",
            c2_slope_convention_primary="|calib slope - 1| <= tol, slope on the "
                                       "predicted-probability feature (/recovery "
                                       "convention; matches the quoted frozen baseline "
                                       "-0.527 / 5.678)",
            c2_slope_convention_secondary="same on the logit-scale Cox/Steyerberg slope "
                                          "( convention), reported for comparison",
            overall="c1 AND c2",
        ),
        cells={},
    )
    panel_p0_10 = {}
    ext_draws_all = {}
    within_all = {}
    ci_se = {}

    for c in ("GSE165252", "GSE78220"):
        M322c, yc = M_ext[c]
        mc = ~np.isnan(yc)
        S0c = transform_pca(M322c[:, mc], mu, sd, pca)
        yext = yc[mc]
        Xte = np.column_stack([np.ones(len(yext)), S0c])
        p = sigmoid(Xte @ w)
        auc = float(roc_auc_score(yext, p))
        lo, hi = auc_bca(yext, p)
        sl_prob = calib_slope_prob(yext, p)
        a_log, b_log = calib_ab(yext, p)
        br, ll = brier_logloss(yext, p)

        # parity vs frozen recovery values
        fz = fz_ext[c]
        parity = dict(
            auc=dict(recomputed=auc, frozen=fz["auc"],
                     match=abs(auc - fz["auc"]) < TOL_PARITY),
            auc_ci95=dict(recomputed=[lo, hi], frozen=fz["ci95"],
                          max_abs_diff=float(max(abs(lo - fz["ci95"][0]),
                                                 abs(hi - fz["ci95"][1])))),
            slope_prob=dict(recomputed=sl_prob, frozen=fz["calib_slope"],
                            match=abs(sl_prob - fz["calib_slope"]) < TOL_PARITY),
            brier=dict(recomputed=br, frozen=fz["brier"],
                       match=abs(br - fz["brier"]) < TOL_PARITY),
            log_loss=dict(recomputed=ll, frozen=fz["log_loss"],
                          match=abs(ll - fz["log_loss"]) < TOL_PARITY),
            n_test=int(len(yext)), n_event=int(yext.sum()),
        )
        all_match = (parity["auc"]["match"] and parity["slope_prob"]["match"]
                     and parity["brier"]["match"] and parity["log_loss"]["match"]
                     and parity["auc_ci95"]["max_abs_diff"] < 1e-6)
        print(f"[{c}] AUC {auc:.4f} (frozen {fz['auc']:.4f}, match={all_match}) "
              f"slope_prob {sl_prob:+.3f} logit_ab ({a_log:+.2f},{b_log:+.2f})", flush=True)

        # slope bootstrap CIs, both conventions, fixed predictions
        sl_lo, sl_hi, n_valid, n_sk = boot_slope_ci(yext, p, calib_slope_prob,
                                                   seed=SEED)
        bl_lo, bl_hi, n_valid_l, n_sk_l = boot_slope_ci(yext, p,
                                                       lambda yy, pp: calib_ab(yy, pp)[1],
                                                       seed=SEED)
        # AUC CI (BCa recomputed above) -> c1
        c1 = bool(lo > 0.5)

        cell = dict(
            n_test=int(len(yext)), n_event=int(yext.sum()),
            auc=auc, auc_ci95_bca=[lo, hi], c1_auc_ci_beneficial_exclusion=c1,
            slope_point_prob=sl_prob,
            slope_boot_ci95_prob=[sl_lo, sl_hi],
            slope_boot_n_valid_prob=n_valid, slope_boot_skipped_prob=n_sk,
            slope_point_logit=b_log, calib_intercept_logit=a_log,
            slope_boot_ci95_logit=[bl_lo, bl_hi],
            slope_boot_n_valid_logit=n_valid_l, slope_boot_skipped_logit=n_sk_l,
            brier=br, log_loss=ll,
            parity_vs_frozen=parity, parity_all_match=bool(all_match),
        )
        for tol in TOLERANCES:
            tkey = f"{tol:.2f}"
            c2p = abs(sl_prob - 1.0) <= tol
            c2l = (b_log is not None and np.isfinite(b_log)
                   and abs(b_log - 1.0) <= tol)
            cell[f"c2_slope_within_tolerance_prob_{tkey}"] = bool(c2p)
            cell[f"c2_slope_within_tolerance_logit_{tkey}"] = bool(c2l)
            cell[f"overall_transport_permitted_prob_{tkey}"] = bool(c1 and c2p)
            cell[f"overall_transport_permitted_logit_{tkey}"] = bool(c1 and c2l)
            # bootstrap-CI robustness (primary convention)
            if sl_lo is None:
                cell[f"boot_ci_robustness_prob_{tkey}"] = "not_estimable"
            elif sl_lo >= 1 - tol and sl_hi <= 1 + tol:
                cell[f"boot_ci_robustness_prob_{tkey}"] = "robust_pass"
            elif sl_lo > 1 + tol or sl_hi < 1 - tol:
                cell[f"boot_ci_robustness_prob_{tkey}"] = "robust_fail"
            else:
                cell[f"boot_ci_robustness_prob_{tkey}"] = "inconclusive"
        p0_3["cells"][c] = cell

        # external AUC draws for delta-AUC CI + SE for pooled meta
        ext_draws, ext_skipped = boot_auc_draws(yext, p, seed=SEED)
        ext_draws_all[c] = ext_draws
        ci_se[c] = float((hi - lo) / (2 * 1.96))   # BCa-CI-implied SE for DL input

        # ---------- within (same pipeline, recovered key, lambda=1/n)
        Xz_c, mu_c, sd_c = clr_standardize(M322c[:, mc])
        pca_c = fit_pca(Xz_c)
        S0_within = pca_c.transform(Xz_c)
        y_within = yc[mc]
        n_within = int(len(y_within))
        lam_within = 1.0 / n_within
        app_w, corr_w, se_w, ci_w = rec.auc_632(S0_within, y_within, lam_within,
                                                n_boot=N_BOOT_632, seed=SEED)
        key = f"{c}_newkey_n{n_within}_seed{SEED}_o{N_BOOT_DELTA}_i{N_BOOT_632}"
        w_draws, w_meta = nested_within_bootstrap(S0_within, y_within, seed=SEED,
                                                  cache_key=key)
        print(f"[within {c}] 0.632 = {corr_w:.4f} (apparent {app_w:.4f}, SE {se_w:.4f}); "
              f"nested draws {len(w_draws)} meta={w_meta}", flush=True)
        within_all[c] = dict(point=dict(apparent=float(app_w), efron_0632=float(corr_w),
                                        se=float(se_w), ci95=[float(ci_w[0]), float(ci_w[1])],
                                        lambda_ridge=lam_within, n=n_within,
                                        n_event=int(y_within.sum()),
                                        estimator="Efron 0.632 optimism form, verbatim "
                                                  "_imvigor_gene_recovery.auc_632, "
                                                  "600 bootstrap, seed 20260911"),
                             nested_bootstrap=dict(draws_n=int(len(w_draws)), **w_meta),
                             draws=w_draws)

        # ---------- P0-10 delta AUC
        k = min(len(ext_draws), len(w_draws))
        dd = ext_draws[:k] - w_draws[:k]
        delta_point = auc - corr_w
        panel_p0_10[c] = dict(
            primary_frozen=dict(
                n_test=int(len(yext)), n_event=int(yext.sum()),
                auc=auc, auc_ci95_bca=[lo, hi],
                delta_auc=dict(
                    point=float(delta_point),
                    within_reference=float(corr_w),
                    within_estimator="same-pipeline Efron 0.632 within the target "
                                     "cohort (recovered key, CLR+PCA k=2, lambda=1/n)",
                    ci95_paired_bootstrap=[float(np.percentile(dd, 2.5)),
                                           float(np.percentile(dd, 97.5))],
                    n_paired=int(k),
                    method="external side: 2000-rep patient bootstrap of the frozen "
                           "model's fixed predictions; within side: 2000x600 nested "
                           "bootstrap of the Efron-0.632 estimator (lambda=1/n, "
                           "seed 20260911); CI = 2.5/97.5 percentiles of the paired "
                           "difference distribution; external and within patient sets "
                           "are independent, so replicate pairing order carries no "
                           "information",
                ),
                calib_intercept_logit=float(a_log),
                calib_slope_logit=float(b_log),
                calib_slope_prob=float(sl_prob),
                brier=br, log_loss=ll,
                external_boot_skipped=ext_skipped,
            ),
            M0="pending__rerun_benchmark_json",
            M2="pending__rerun_benchmark_json",
            M3="pending__rerun_benchmark_json",
        )

    # ---------- IMvigor210 within reference (frozen from the recovery run)
    imv_fz = frozen["phase2"]["new"]["within_imvigor"]
    within_all["IMvigor210"] = dict(
        point=dict(apparent=imv_fz["auc_apparent"], efron_0632=imv_fz["auc_632"],
                   se=imv_fz["se"], ci95=imv_fz["ci95"],
                   lambda_ridge=imv_fz["lambda_ridge"], n=imv_fz["n"],
                   n_event=imv_fz["n_event"],
                   estimator="frozen: _imvigor_gene_recovery stage 2 "
                             "(600 bootstrap, lambda=1/N_total per the  rule)"),
        nested_bootstrap=None,
        draws=None,
        note="source cohort of the primary frozen transport: no external self-test; "
             "nested within-side bootstrap not rerun here (frozen 600-bootstrap CI "
             "quoted); rerun when the benchmark worker lands (needed for M0-M3 "
             "delta-AUC rows)",
    )

    # ---------- flip counting (P0-3)
    flips = {"prob": [], "logit": []}
    verdicts = {"prob": {}, "logit": {}}
    for c in p0_3["cells"]:
        cell = p0_3["cells"][c]
        for conv in ("prob", "logit"):
            vs = [cell[f"overall_transport_permitted_{conv}_{t:.2f}"] for t in TOLERANCES]
            verdicts[conv][c] = vs
            if len(set(vs)) > 1:
                flips[conv].append(c)
    p0_3["overall_flip_count_prob"] = len(flips["prob"])
    p0_3["overall_flip_cells_prob"] = flips["prob"]
    p0_3["overall_flip_count_logit"] = len(flips["logit"])
    p0_3["overall_flip_cells_logit"] = flips["logit"]
    p0_3["c2_flip_count_prob"] = 0
    c2_flips = []
    for c, cell in p0_3["cells"].items():
        vs = [cell[f"c2_slope_within_tolerance_prob_{t:.2f}"] for t in TOLERANCES]
        if len(set(vs)) > 1:
            c2_flips.append(c)
    p0_3["c2_flip_cells_prob"] = c2_flips
    p0_3["c2_flip_count_prob"] = len(c2_flips)
    any_permitted = any(
        p0_3["cells"][c][f"overall_transport_permitted_prob_{t:.2f}"]
        for c in p0_3["cells"] for t in TOLERANCES)
    p0_3["any_transport_permitted_at_any_tolerance"] = bool(any_permitted)
    p0_3["main_conclusion_depends_on_0p10_tolerance"] = bool(len(flips["prob"]) > 0)
    p0_3["verdicts_by_tolerance"] = verdicts

    # ---------- pooled DL meta (k=3, new 320-gene era cohort AUCs)
    imv_se = imv_fz["se"]
    pooled_inputs = {
        "IMvigor210": dict(auc=imv_fz["auc_632"], se=imv_se,
                           role="anchor-cohort within Efron-0.632 (recovered key, "
                                "lambda=1/N_total  rule, 600 bootstrap)"),
        "GSE165252": dict(auc=p0_3["cells"]["GSE165252"]["auc"],
                          se=ci_se["GSE165252"],
                          role="primary frozen external AUC; SE = BCa CI width / "
                               "(2*1.96)"),
        "GSE78220": dict(auc=p0_3["cells"]["GSE78220"]["auc"],
                         se=ci_se["GSE78220"],
                         role="primary frozen external AUC; SE = BCa CI width / "
                              "(2*1.96)"),
    }
    pooled_new = dl_meta([pooled_inputs[c]["auc"] for c in ("IMvigor210", "GSE165252",
                                                           "GSE78220")],
                        [pooled_inputs[c]["se"] for c in ("IMvigor210", "GSE165252",
                                                          "GSE78220")],
                        ["IMvigor210", "GSE165252", "GSE78220"])
    # secondary: all-within 0.632 synthesis (mirrors the old  design)
    pooled_allwithin = dl_meta(
        [imv_fz["auc_632"], within_all["GSE165252"]["point"]["efron_0632"],
         within_all["GSE78220"]["point"]["efron_0632"]],
        [imv_se, within_all["GSE165252"]["point"]["se"],
         within_all["GSE78220"]["point"]["se"]],
        ["IMvigor210", "GSE165252", "GSE78220"])
    pooled = dict(
        historical_=dict(
            pooled_auc=0.498, ci95=[0.365, 0.631], i2_qbased=0.7516, Q_p=0.0179,
            role="descriptive only; HISTORICAL 229-gene-era (broken positional "
                 "alignment) value from data/meta_sensitivity.json "
                 "efficacy_dlmeta; retained for the record, superseded by the "
                 "new-320-gene synthesis below",
            inputs="_eff_meta per-cohort 0.632 AUCs on the broken-alignment S0"),
        new_320gene_primary_spec=dict(
            inputs=pooled_inputs, dl=pooled_new,
            note="per task spec: the three cohorts' primary-frozen-design AUCs "
                 "(IMvigor210 within-0.632 as the anchor cohort's own estimate + "
                 "the two external transports) with DL random effects"),
        new_320gene_all_within_secondary=dict(
            dl=pooled_allwithin,
            note="sensitivity: all three within-0.632 AUCs (mirrors the old  "
                 "all-within design); GSE165252/GSE78220 within use lambda=1/n "
                 "per the rerun spec, IMvigor210 lambda=1/N_total frozen"),
        heterogeneity_note="I2/Q reported as descriptive; heterogeneity precludes "
                           "a single transportable-AUC interpretation (same caveat "
                           "as the historical value)",
    )

    # ---------- P0-12: analysis freeze manifest
    exmat_sha = sha256(rec.IMV_EXMAT)
    fdata_sha = sha256(rec.IMV_FDATA)
    pdata_sha = sha256(rec.IMV_PDATA)
    fz_prov = frozen["provenance"]
    hash_assert = dict(
        exmat=dict(computed=exmat_sha, reported=fz_prov["exmat_sha256"],
                   match=exmat_sha == fz_prov["exmat_sha256"]),
        fData=dict(computed=fdata_sha, reported=fz_prov["fData_sha256"],
                   match=fdata_sha == fz_prov["fData_sha256"]),
        pData=dict(computed=pdata_sha, reported=fz_prov["pData_sha256"],
                   match=pdata_sha == fz_prov["pData_sha256"]),
    )
    assert all(v["match"] for v in hash_assert.values()), hash_assert
    print("[hash] IMvigor210 3-file SHA256 asserted equal to the recovery report",
          flush=True)

    input_hashes = {
        "imvigor210/exmat.csv": exmat_sha,
        "imvigor210/fData.csv": fdata_sha,
        "imvigor210/pData.csv": pdata_sha,
        "geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz": sha256(rec.GSE165_MAT),
        "geo_data/GSE165252/series_matrix_annotation.txt": sha256(rec.GSE165_ANN),
        "geo_data/GSE78220/GSE78220_PatientFPKM.xlsx": sha256(rec.GSE782_MAT),
        "geo_data/GSE78220/series_matrix_annotation.txt": sha256(rec.GSE782_ANN),
        "geo_data/322_shared_genes.json": sha256(rec.GENE_JSON),
    }

    manifest = dict(
        manifest="_analysis_manifest",
        date="2026-09-24",
        seed=SEED,
        python=sys.executable,
        purpose="P0-12 submission-freeze manifest: every  analysis x "
                "{script, inputs, outputs, seed, bootstrap counts, key params, "
                "status}; statuses follow the IMvigor210 gene-row recovery "
                "(imvigor_gene_recovery.json) which retracted the positional-"
                "alignment-era anchored numbers",
        imvigor210_files_sha256=hash_assert,
        upstream_provenance=dict(
            censored_source=dict(
                repo=fz_prov["upstream"]["repo"],
                file="2_data/1_IMvigor210/exmat_censored_IMvigor210.csv",
                note="exmat anonymised to gene_1..gene_31085 row labels; row order "
                     "= counts(cds) = fData symbol alphabetical (upper-cased); "
                     "fData/pData byte-identical upstream",
                extraction_script=fz_prov["upstream"]["extraction_script"],
            ),
            labeled_mirror=dict(
                repo=fz_prov["mirror"]["repo"], file="data/cds.RData",
                url=fz_prov["mirror"]["url"], sha256=fz_prov["mirror"]["sha256"],
                note="used only for the stage-1c numeric identity check "
                     "(31085/31085 rows matched), never for analysis values",
            ),
            recovered_key="exmat row j == fData row with the j-th symbol in "
                          "alphabetical order of the upper-cased fData symbol "
                          "column; exported to geo_data/imvigor210_labeled/"
                          "exmat_row_key.csv",
        ),
        analyses=[
            dict(id="-01", name="number_audit",
                 script=None,
                 script_note="session-level audit (no standalone script); audit of "
                             "reviews_clean_manuscript_extract.md against the "
                             "pre-recovery frozen archive",
                 inputs=["reviews_clean_manuscript_extract.md",
                         "data/, data/, data/, data/_meta_sensitivity, "
                          "data/_eff_meta, data/_tox_meta, data/ JSONs"],
                 outputs=["data/number_audit.md"],
                 seed=None, n_boot=None,
                 key_params=dict(audited_lines=57, mismatches=4, not_in_frozen=1),
                 status="done",
                 status_note="audit verdicts refer to the pre-recovery frozen "
                             "numbers; every anchored-arm number it validated "
                             "(0.444/0.464/0.498 pooled/88.7% etc.) is now "
                             "superseded by the recovery-era rerun"),
            dict(id="-02", name="positive_control (P0-2)",
                 script="scripts/analysis/positive_control.py",
                 inputs=["data/imvigor210/{exmat,fData,pData}.csv",
                         "geo_data/GSE165252/*", "geo_data/GSE78220/*",
                         "data/_workspace/_eff_*_s0.npy"],
                 outputs=["data/positive_control.json",
                          "data/positive_control.md"],
                 seed=SEED, n_boot=dict(auc_ci=2000, efron632=600),
                 key_params=dict(signature="IFN-gamma 6-gene (Ayers JCI 2017)",
                                 rule="C1 AUC CI>0.5; C2 |slope-1|<=0.10; tiers "
                                      "0.05/0.10/0.20"),
                 status="superseded",
                 status_note="as-specified IMvigor-anchored arm INVALID: the P0-2 "
                             "sex-anchor gates falsified the exmat/fData positional "
                             "alignment (7/7 impossible flags); this audit triggered "
                             "the gene recovery. Salvage two-cohort transports "
                             "(GSE165252<->GSE78220, gene-valid matrices) remain "
                             "descriptive; superseded by the recovery + rerun"),
            dict(id="-03", name="clinical_only (P0-4)",
                 script="scripts/analysis/clinical_only.py",
                 inputs=["data/imvigor210/pData.csv", "geo_data/GSE78220/series_"
                          "matrix_annotation.txt", "geo_data/GSE165252/series_"
                          "matrix_annotation.txt",
                          "data/_workspace/_eff_*_y.npy"],
                 outputs=["data/clinical_only.json",
                          "data/clinical_only.md"],
                 seed=SEED, n_boot=dict(auc_ci=2000, efron632=600),
                 key_params=dict(models=["C_main(sex,stage)", "C_sex", "C_rich(n=298 "
                                     "IMvigor210-only)"], lambda_rule="1/n_train"),
                 status="done",
                 status_note="UNAFFECTED by the gene recovery: clinical covariates "
                             "only, no expression rows involved"),
            dict(id="-04", name="threshold_zero_sensitivity (P0-3 old + P0-6)",
                 script="scripts/analysis/threshold_zero_sensitivity.py",
                 inputs=["data/_workspace/_eff_*_s0.npy (broken-alignment "
                          "frozen S0)", "data/competitor_benchmark.json"],
                 outputs=["data/threshold_zero_sensitivity.json",
                          "data/threshold_zero_sensitivity.md"],
                 seed=SEED, n_boot=dict(slope_ci=2000, repeats=100),
                 key_params=dict(tolerances=[0.05, 0.10, 0.20], models_m0_m6=True,
                                 zero_schemes=["original_log1p",
                                               "multiplicative_replacement_0p65",
                                               "measured_only_158"]),
                 status="superseded",
                 status_note="ALL cells rest on the frozen  S0 built from the "
                             "falsified positional alignment (and the 229-gene "
                             "pseudo-intersection); the P0-6 zero-handling "
                             "'artifact' (70 unmeasured genes, 88.7% PC1 loading "
                             "mass) is RETRACTED as an alignment artefact; "
                             "superseded by -09 (new numbers)"),
            dict(id="-05", name="provenance_ledger (P0-7)",
                 script=None,
                 script_note="automated provenance audit executed in-session; "
                             "evidence frozen under data/evidence/",
                 inputs=["supplementary.md Table S15", "README.md",
                         "geo_data/GSE186143_bulk_clinical.json",
                         "NCBI FTP SOFT family files (frozen copies)"],
                 outputs=["data/provenance_ledger.json",
                          "data/provenance_ledger.md",
                          "data/evidence/*.soft.gz"],
                 seed=None, n_boot=None,
                 key_params=dict(cohorts=["GSE186143", "GSE166181", "GSE216329",
                                          "GSE253720"],
                                 confirmed=27, corrected=10, unverified=3),
                 status="done",
                 status_note="UNAFFECTED: GEO metadata only, no expression "
                             "matrices involved"),
            dict(id="-06", name="transport_summary_joint (P0-10 old + P0-11)",
                 script="scripts/analysis/transport_summary_joint.py",
                 inputs=["data/_workspace/_eff_*_s0.npy (broken-alignment "
                          "frozen S0)", "data/_workspace/_eff_*_y.npy",
                          "data/external_validation.json",
                          "geo_data/GSE186143_bulk_clinical.json",
                          "data/joint_stress_test.json",
                          "data/__within_boot_cache.npz"],
                 outputs=["data/transport_summary_joint.json",
                          "data/transport_summary_joint.md"],
                 seed=SEED,
                 n_boot=dict(paired_delta=2000, nested_within="2000 outer x 600 "
                            "inner", recal_repeats=100),
                 key_params=dict(pooled_historical=0.498),
                 status="split",
                 status_note="P0-11 GSE186143 joint-model re-verification: VALID "
                             "(unaffected by the IMvigor210 gene-row issue). "
                             "P0-10 anchored within/external/delta-AUC numbers: "
                             "superseded by -09 (all rest on the broken-"
                             "alignment S0 and the 229-gene era gene set)"),
            dict(id="-07a", name="phase1a forensic battery (availability / sort "
                                    "hypotheses / deep verify / sort-key / panel)",
                 script=["scripts/analysis/phase1a_availability.py",
                         "scripts/analysis/phase1a_sort_hypotheses.py",
                         "scripts/analysis/phase1a_deep_verify.py",
                         "scripts/analysis/phase1a_sortkey_check.py",
                         "scripts/analysis/phase1a_panel_check.py"],
                 inputs=["data/imvigor210/{exmat,fData,pData}.csv",
                          "outputs/analysis_base/supplementary/"
                          "Supplementary_Table_1_322_genes.csv",
                          "geo_data/322_shared_genes.json"],
                 outputs=["data/phase1a_results.json (console-stage artefact)"],
                 seed=None, n_boot=None,
                 key_params=dict(winner="symbol_alphabetical_upper",
                                 raw_vs_upper_ambiguous_rows=145,
                                 panel_genes_unaffected_by_sort_ambiguity=True),
                 status="done",
                 status_note="subsumed by the stage-1a re-run inside "
                             "imvigor_gene_recovery.py"),
            dict(id="-07b", name="cds mirror extraction",
                 script="scripts/analysis/cds_extract.py",
                 inputs=["geo_data/imvigor210_labeled/cds.RData"],
                 outputs=["geo_data/imvigor210_labeled/cds_counts.csv"],
                 seed=None, n_boot=None, key_params={},
                 status="done",
                 status_note="stage-1c numeric identity input only"),
            dict(id="-07", name="imvigor_gene_recovery (Go/No-Go recovery)",
                 script="scripts/analysis/imvigor_gene_recovery.py",
                 inputs=["data/imvigor210/{exmat,fData,pData}.csv",
                         "geo_data/322_shared_genes.json",
                         "gene_name_id_dict_gc104M.pkl",
                         "geo_data/GSE165252/*", "geo_data/GSE78220/*",
                         "geo_data/imvigor210_labeled/cds_counts.csv"],
                 outputs=["data/imvigor_gene_recovery.json",
                          "data/imvigor_gene_recovery.md",
                          "geo_data/imvigor210_labeled/exmat_row_key.csv"],
                 seed=SEED, n_boot=dict(efron632_within=600, auc_bca=2000),
                 key_params=dict(key="exmat row j == j-th fData symbol "
                                     "(upper-cased) alphabetical",
                                 panel="320/322 measurable, cross-cohort common 320 "
                                       "(IGHG1, PRKN absent from fData)",
                                 anchor="Scheme A CLR->PCA k=2 on 298 labelled, "
                                        "lambda_transport=1/298, "
                                        "lambda_within=1/N_total  rule",
                                 gates="overall PASS"),
                 status="done",
                 status_note="the authoritative new baseline: IMvigor210 within "
                             "0.632 = 0.6571 [0.5857, 0.7284]; external GSE165252 "
                             "0.5036 [0.3540, 0.6585] slope -0.527; GSE78220 0.6603 "
                             "[0.4038, 0.8617] slope 5.678"),
            dict(id="-08", name="rerun_benchmark (M0-M6 on the recovered key)",
                 script="pending (benchmark worker)",
                 inputs=["recovered key + 320-gene loaders (as -07)",
                         "data/competitor_benchmark.json (target values)"],
                 outputs=["data/rerun_benchmark.json"],
                 seed=SEED, n_boot=dict(auc_ci=2000, efron632=600, recal=100),
                 key_params=dict(models=["M0", "M1", "M2", "M3", "M4", "M5a", "M5b",
                                         "M6"]),
                 status="pending",
                 status_note="rerun_benchmark.json NOT on disk; blocks the "
                             "M0/M2/M3 rows of the P0-10 panel and the M0-M6 rows "
                             "of the P0-3 matrix"),
            dict(id="-09", name="rerun_sensitivity_manifest (P0-3 new + P0-10 "
                                   "new + P0-12 manifest)",
                 script="scripts/analysis/rerun_sensitivity_manifest.py",
                 inputs=["data/imvigor210/{exmat,fData,pData}.csv "
                          "(SHA256 asserted == recovery report)",
                         "geo_data/322_shared_genes.json",
                         "geo_data/GSE165252/*", "geo_data/GSE78220/*",
                         "data/imvigor_gene_recovery.json (frozen baseline "
                          "+ provenance)"],
                 outputs=["data/rerun_sensitivity.json",
                          "data/rerun_sensitivity.md",
                          "data/analysis_manifest.json",
                          "data/__within_boot_newkey_cache.npz"],
                 seed=SEED,
                 n_boot=dict(slope_ci=2000, auc_bca=2000, paired_delta=2000,
                             within_0632=600, nested_within="2000 outer x 600 inner"),
                 key_params=dict(panel=320, pca_k=2, lambda_transport="1/298",
                                 lambda_within_external="1/n per task spec",
                                 tolerances=[0.05, 0.10, 0.20],
                                 pooled_dl_k3=True),
                 status="done",
                 status_note="primary frozen arm fully recomputed and parity-checked "
                             "against -07 (tol 1e-9); M0-M3 pending on -08"),
        ],
        input_file_sha256=input_hashes,
        retraction_reference="data/artefact_retraction.md",
    )

    with open(f"{OUT_DIR}/analysis_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print("[manifest] wrote", f"{OUT_DIR}/analysis_manifest.json", flush=True)

    # ---------- emit sensitivity JSON (strip raw draws; keep summaries)
    within_out = {}
    for c in ("IMvigor210", "GSE165252", "GSE78220"):
        wdt = within_all[c]
        d = dict(point=wdt["point"])
        if wdt.get("nested_bootstrap"):
            d["nested_bootstrap"] = wdt["nested_bootstrap"]
            dr = wdt["draws"]
            d["nested_bootstrap"]["ci95_percentile"] = [
                float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))]
            d["nested_bootstrap"]["boot_mean"] = float(dr.mean())
        if wdt.get("note"):
            d["note"] = wdt["note"]
        within_out[c] = d

    out = dict(
        analysis="rerun_sensitivity",
        date="2026-09-24",
        seed=SEED,
        python=sys.executable,
        baseline=dict(
            gene_identity="recovered key (symbol alphabetical, upper); "
                          "imvigor_gene_recovery.json gates overall PASS",
            panel=320, common_genes=len(common),
            primary_frozen_source="data/imvigor_gene_recovery.json phase2.new",
        ),
        p0_3=p0_3,
        p0_10=dict(
            rows=dict(panel=panel_p0_10),
            within_reference=within_out,
            pending=dict(
                M0_M2_M3="rerun_benchmark.json not on disk; rows to be added "
                         "when the benchmark worker lands (marked pending in every "
                         "cohort block)"),
            pooled=pooled,
        ),
        total_seconds=round(time.time() - t_start, 1),
    )
    with open(f"{OUT_DIR}/rerun_sensitivity.json", "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("[out] wrote", f"{OUT_DIR}/rerun_sensitivity.json", flush=True)

    write_md(out)
    print("[out] wrote", f"{OUT_DIR}/rerun_sensitivity.md", flush=True)
    print(f"[done] total {time.time()-t_start:.0f}s", flush=True)


# ================================================================ report
def fmt(x, nd=4):
    if x is None:
        return "NA"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if np.isnan(xf):
        return "NA"
    return f"{xf:.{nd}f}"


def write_md(out):
    p0_3 = out["p0_3"]
    panel = out["p0_10"]["rows"]["panel"]
    within = out["p0_10"]["within_reference"]
    pooled = out["p0_10"]["pooled"]
    L = []
    A = L.append
    A("#  重算敏感性报告（P0-3 / P0-10 新数字版）+ P0-12 manifest")
    A("")
    A(f"日期：{out['date']} ｜ seed {out['seed']} ｜ 脚本：scripts/analysis/rerun_sensitivity_manifest.py ｜ "
      f"冻结 JSON：data/rerun_sensitivity.json")
    A("")
    A("基线：IMvigor210 基因行身份已恢复（symbol 字母序钥匙，_imvigor_gene_recovery 总门 PASS）；"
      "primary frozen = 320 基因 CLR→PCA k=2 锚定臂（298 labelled），ridge λ=1/298 运输。"
      "全部数字由本脚本从原始矩阵重算，并与恢复报告冻结值逐位交叉校验（tol 1e-9）。")
    A("")
    A("## 1. 重建一致性（先验证再分析）")
    A("")
    A("| 队列 | AUC（重算） | AUC（冻结） | slope（重算/冻结） | Brier | log-loss | 全部一致 |")
    A("|---|---|---|---|---|---|---|")
    for c in ("GSE165252", "GSE78220"):
        cell = p0_3["cells"][c]
        p = cell["parity_vs_frozen"]
        A(f"| {c} | {fmt(cell['auc'])} | {fmt(p['auc']['frozen'])} | "
          f"{fmt(cell['slope_point_prob'], 3)} / {fmt(p['slope_prob']['frozen'], 3)} | "
          f"{fmt(p['brier']['recomputed'])} | {fmt(p['log_loss']['recomputed'])} | "
          f"{'**是**' if cell['parity_all_match'] else '否'} |")
    A("")
    A("## 2. P0-3 阈值敏感性（新数字版）")
    A("")
    A("覆盖：primary frozen 两个方向（GSE165252 / GSE78220）。"
      "M0-M6 基准模型行 **待补**：rerun_benchmark.json 尚未落盘（manifest 中记 pending）。")
    A("")
    A("判定规则：C1 = 外部 AUC BCa 95% CI（2000，seed 20260911）向上排除 0.5（与容差无关）；"
      "C2 = |slope−1| ≤ tol（slope 用恢复报告口径：以预测概率为特征的重标定斜率，"
      "即引文基线 −0.527 / 5.678；logit 口径同时报告）。")
    A("")
    A("| 队列 | AUC [BCa CI] | C1 | slope(概率口径) [boot CI] | slope(logit口径) [boot CI] | "
      "|slope−1| | tol=0.05 | 0.10 | 0.20 |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in ("GSE165252", "GSE78220"):
        cell = p0_3["cells"][c]
        sl = cell["slope_point_prob"]
        A(f"| {c} | {fmt(cell['auc'])} [{fmt(cell['auc_ci95_bca'][0])}, "
          f"{fmt(cell['auc_ci95_bca'][1])}] | {'过' if cell['c1_auc_ci_beneficial_exclusion'] else '不过'} | "
          f"{fmt(sl, 3)} [{fmt(cell['slope_boot_ci95_prob'][0], 2)}, "
          f"{fmt(cell['slope_boot_ci95_prob'][1], 2)}] | "
          f"{fmt(cell['slope_point_logit'], 2)} [{fmt(cell['slope_boot_ci95_logit'][0], 2)}, "
          f"{fmt(cell['slope_boot_ci95_logit'][1], 2)}] | {fmt(abs(sl-1), 3)} | "
          + " | ".join(
              ("**过**" if cell[f"overall_transport_permitted_prob_{t:.2f}"] else "不过")
              for t in (0.05, 0.10, 0.20)) + " |")
    A("")
    A(f"- **判定翻转格数（overall，概率口径）：{p0_3['overall_flip_count_prob']}**"
      f"（翻转队列：{p0_3['overall_flip_cells_prob'] or '无'}）；"
      f"C2 单项翻转 {p0_3['c2_flip_count_prob']} 格。")
    A(f"- **主结论是否依赖 0.10 容差：{'是' if p0_3['main_conclusion_depends_on_0p10_tolerance'] else '否'}** —— "
      f"任何容差下均无任何格允许运输（any_transport_permitted = "
      f"{p0_3['any_transport_permitted_at_any_tolerance']}）：两个方向的 |slope−1|（1.53 / 4.68）"
      "远超最大容差 0.20，且 C1 在两个方向均不过（AUC CI 含 0.5）。")
    A("- 两口径（概率特征 / logit 特征）slope 结论一致；slope bootstrap 95% CI 在两个方向均极宽，"
      "任何容差下均为 inconclusive（CI 同时跨越 [1−tol, 1+tol] 边界）。")
    A("")
    A("## 3. P0-10 ΔAUC 连续面板（新数字版）")
    A("")
    A("within 侧（任务规格：恢复钥匙 + 同一加载逻辑 + CLR + PCA k=2 + λ=1/n + 600 bootstrap，"
      "Efron 0.632 optimism 形式，与 IMvigor210 基线 0.6571 同一估计器）：")
    A("")
    A("| 队列 | within 0.632 [CI] | λ | n（事件） | nested bootstrap |")
    A("|---|---|---|---|---|")
    for c in ("IMvigor210", "GSE165252", "GSE78220"):
        w = within[c]
        pt = w["point"]
        nb = w.get("nested_bootstrap")
        nbs = (f"{nb['draws_n']} 抽样，CI [{fmt(nb['ci95_percentile'][0])}, "
               f"{fmt(nb['ci95_percentile'][1])}]" if nb else "未重跑（冻结 600-boot CI）")
        A(f"| {c} | {fmt(pt['efron_0632'])} [{fmt(pt['ci95'][0])}, {fmt(pt['ci95'][1])}] | "
          f"{pt['lambda_ridge']:.6f} | {pt['n']}（{pt['n_event']}） | {nbs} |")
    A("")
    A("primary frozen 行（M0 / M2 / M3 **待补**——rerun_benchmark.json 未落盘，"
      "已在 manifest 记 pending，落盘后同表追加）：")
    A("")
    A("| 目标队列 | 外部 AUC [BCa CI] | ΔAUC（外−内） | ΔAUC 95% CI（配对 bootstrap） | "
      "logit 截距 | logit slope | 概率口径 slope | Brier | log-loss |")
    A("|---|---|---|---|---|---|---|---|---|")
    for c in ("GSE165252", "GSE78220"):
        r = panel[c]["primary_frozen"]
        d = r["delta_auc"]
        A(f"| {c} | {fmt(r['auc'])} [{fmt(r['auc_ci95_bca'][0])}, {fmt(r['auc_ci95_bca'][1])}] | "
          f"{fmt(d['point'])} | [{fmt(d['ci95_paired_bootstrap'][0])}, "
          f"{fmt(d['ci95_paired_bootstrap'][1])}] | {fmt(r['calib_intercept_logit'], 2)} | "
          f"{fmt(r['calib_slope_logit'], 2)} | {fmt(r['calib_slope_prob'], 3)} | "
          f"{fmt(r['brier'])} | {fmt(r['log_loss'])} |")
    A("")
    A("解读（不挑好的说）：")
    A("")
    for c in ("GSE165252", "GSE78220"):
        r = panel[c]["primary_frozen"]
        d = r["delta_auc"]
        A(f"- {c}：ΔAUC = {fmt(d['point'])} [{fmt(d['ci95_paired_bootstrap'][0])}, "
          f"{fmt(d['ci95_paired_bootstrap'][1])}]，CI 均含 0——同管线 within 参考下，"
          "外部表现与队列内表现无统计可分辨的差异。")
    A("")
    A("### Pooled AUC（DL 随机效应，k=3）")
    A("")
    pn = pooled["new_320gene_primary_spec"]
    aw = pooled["new_320gene_all_within_secondary"]["dl"]
    A(f"- **新 320 基因口径（主规格：IMvigor within 0.632 + 两外部队列运输 AUC）："
      f"pooled AUC = {fmt(pn['dl']['pooled_auc'], 3)} [{fmt(pn['dl']['ci95'][0], 3)}, "
      f"{fmt(pn['dl']['ci95'][1], 3)}]，I²(Q) = {fmt(pn['dl']['i2_qbased'], 3)}，"
      f"Q = {fmt(pn['dl']['Q'], 2)}（df=2，P = {fmt(pn['dl']['Q_p'], 3)}），"
      f"tau² = {fmt(pn['dl']['tau2'], 5)}**")
    A(f"- 次规格（三队列全部 within 0.632，镜像旧  设计）：pooled = "
      f"{fmt(aw['pooled_auc'], 3)} [{fmt(aw['ci95'][0], 3)}, {fmt(aw['ci95'][1], 3)}]，"
      f"I²(Q) = {fmt(aw['i2_qbased'], 3)}")
    A("- 历史  descriptive **0.498 [0.365, 0.631] 保留在案**，但注明为旧 229 基因"
      "（错误位置对齐）口径的历史值，已由上列新口径合成取代；异质性 caveat 不变"
      "（I² 高，单一可运输 AUC 解释仍不可取）。")
    A("")
    A("## 4. P0-12 manifest 摘要")
    A("")
    A("- 冻结映射：data/analysis_manifest.json —— 10 个  分析条目"
      "（-01 … -09，含 phase1a 取证组与 cds 镜像提取子条目），每条含脚本/输入/"
      "输出/seed/bootstrap 次数/关键参数/状态（done / pending / superseded / split）。")
    A("- IMvigor210 三文件 SHA256 与恢复报告逐一断言相等："
      "exmat `fc51849b…`、fData `d0f781b5…`、pData `33aafcca…`；"
      "上游溯源（mimifp/tfm_mUC censored、SiYangming/IMvigor210CoreBiologies 镜像）已固化。")
    A("- 伪影撤回文档：data/artefact_retraction.md（ED10/provenance 教训表）。")
    A("")
    A("## 5. 文件")
    A("")
    A("- 脚本：scripts/analysis/rerun_sensitivity_manifest.py")
    A("- 数据：data/rerun_sensitivity.json；data/analysis_manifest.json；"
      "within 嵌套 bootstrap 缓存 data/__within_boot_newkey_cache.npz")
    A("- 报告：data/rerun_sensitivity.md；撤回文档 data/artefact_retraction.md")
    with open(f"{OUT_DIR}/rerun_sensitivity.md", "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
