#!/usr/bin/env python3
""" -- Nature Cancer revision , anchored-arm full rerun.

The IMvigor210 gene-row identity has been RECOVERED (_imvigor_gene_recovery:
exmat row j == the j-th fData row under upper-cased SYMBOL alphabetical sort;
the old "positional alignment" is falsified), so every anchored analysis must be
rerun.  This script reruns, with REAL computation only:

  1. M0-M6 (incl. M5a/M5b) x 3 target cohorts LOCO competitor benchmark, 
     conventions verbatim (pooled-two training, ridge lambda=1/n_train,
     AUC + 2000-rep BCa CI, logit-scale calibration intercept/slope, Brier,
     log-loss, two-criterion licensing rule: AUC CI excludes 0.5 in the
     beneficial direction AND |slope-1| <= 0.10);
  2. within-cohort same-pipeline reference for the three cohorts (Efron 0.632,
     600 bootstrap, seed 20260911) in the /-recovery convention (global
     lambda 1/N_total, optimism-corrected), plus the  OOB convention as a
     secondary field;
  3. calibration-curve block (_supporting_analyses analysis-1 machinery):
     M2 frozen per-patient prediction distribution, target event rate contrast,
     5 equal-frequency bins, M3 100-repeat recalibration before/after slopes,
     GSE165252 M3 slope CI;
  4. true-external pairwise transports, all 6 ordered directions (F1-F3 and the
     source-swap sensitivity directions),  pairwise block verbatim;
  5. near-miss check: does ANY direction x model have an AUC CI excluding 0.5
     in the beneficial direction under the new numbers?  Full-table confirmation.

Anchor pipeline (sole authority): scripts/analysis/imvigor_gene_recovery.py --
panel = 322 authoritative symbols, IMvigor availability by recovered
symbol->row key (320/322; IGHG1 and PRKN absent), cross-cohort common set =
320 genes, IMvigor in-panel all-zero genes = 0.  Everything else is imported
VERBATIM: _competitor_benchmark (ridge algebra, metrics_block, auc_bca,
recal_repeats, CORAL, location-scale), _competitor_benchmark_
(prob_scale_slope, efron632_within), _supporting_analyses
(equal_freq_bins, r_arrays, recal_repeats_).

HARD SELF-CHECKS against the frozen new-anchor baselines (tolerance 1e-9;
the script ABORTS on any mismatch -- numbers are never fabricated):
  within IMvigor210 0.632 = 0.6570746710270073 [0.5857004961602726, 0.728448845893742]
  F1-style GSE165252: AUC 0.5036231884057971 [0.3540206604291994, 0.6584721535961872]
                      prob-scale slope -0.526789606103434
  F1-style GSE78220:  AUC 0.6602564102564102 [0.4038454964351145, 0.8616583228383701]
                      prob-scale slope 5.678169924426832

Emits: data/rerun_benchmark.json (+ .md report with the old-vs-new
comparison against the frozen / artifacts).
Reproducibility: managed venv python, seed 20260911 everywhere.
"""
import json
import os
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
PKG = f"{ROOT}/outputs/analysis_base"
SCRIPTS = f"{PKG}/scripts"
sys.path.insert(0, SCRIPTS)

# ---- sole-authority anchor pipeline + verbatim project algebra (imports) ----
import imvigor_gene_recovery as rec          # recovered key, loaders,  auc_632
from competitor_benchmark import (
    sigmoid, fit_logistic, fit_logistic_gen, calib_ab, auc_bca, metrics_block,
    recal_repeats, coral_align, locscale_align, _clean,
    COHORTS, SEED, N_BOOT, SLOPE_TOL, M4_GAMMA)
from competitor_benchmark_ import prob_scale_slope, efron632_within
import supporting_analyses as sup            # equal_freq_bins, r_arrays, recal_repeats_

OUTDIR = f"{PKG}/data/"
OUT_JSON = f"{OUTDIR}/rerun_benchmark.json"
OUT_MD = f"{OUTDIR}/rerun_benchmark.md"
V75_JSON = f"{PKG}/data/competitor_benchmark.json"
V71_PS_JSON = f"{PKG}/data/pca_stability.json"
V82_RECOVERY_JSON = f"{OUTDIR}/imvigor_gene_recovery.json"

MODEL_ORDER = ["M0", "M1", "M2", "M3", "M4", "M5a", "M5b", "M6"]
N_BOOT_632 = 600
TOL = 1e-9

# ---- frozen new-anchor baselines (recovery worker, imvigor_gene_recovery.json) ----
FROZEN = dict(
    within_imvigor210=dict(
        auc_apparent=0.6658567774936061,
        auc_632=0.6570746710270073,
        se=0.036415395340170784,
        ci95=[0.5857004961602726, 0.728448845893742]),
    f1_style={
        "GSE165252": dict(
            auc=0.5036231884057971,
            ci95=[0.3540206604291994, 0.6584721535961872],
            calib_slope_probability_scale=-0.526789606103434,
            brier=0.262702511059918,
            log_loss=0.8141936892958811),
        "GSE78220": dict(
            auc=0.6602564102564102,
            ci95=[0.4038454964351145, 0.8616583228383701],
            calib_slope_probability_scale=5.678169924426832,
            brier=0.22802786705931374,
            log_loss=0.6480391960560606)})


# ================================================================ anchor (recovered key)
def build_anchor():
    """Rebuild the anchored S0 for all three cohorts with the RECOVERED gene-row
    identity, verbatim _imvigor_gene_recovery stage-2 logic: panel symbols
    matched directly against the recovered symbol->row key; Scheme A
    CLR -> z -> PCA k=2 fit on the 298 labelled IMvigor210 samples;
    transform-only to the external cohorts."""
    import gzip
    import pickle
    import openpyxl

    panel = [g.upper() for g in json.load(open(rec.GENE_JSON))["genes"]]
    fd, fd_ent, fd_sym_raw, fd_sym_up = rec.load_fdata()
    n_fd = len(fd)
    order_up = sorted(range(n_fd), key=lambda i: fd_sym_up[i])
    row_of_sym = {fd_sym_up[order_up[j]]: j for j in range(n_fd)}

    imv_avail = {s for s in panel if s in row_of_sym}
    g = pickle.load(open(rec.SYM2ENS_PKL, "rb"))
    sym2ens = {k.upper(): v for k, v in g.items()}
    with gzip.open(rec.GSE165_MAT, "rt") as f:
        f.readline()
        cnt_ens = set(line.rstrip("\n").split("\t")[0].split(".")[0] for line in f)
    gse165_avail = {s for s in panel if sym2ens.get(s) in cnt_ens}
    wb = openpyxl.load_workbook(rec.GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    rr = list(ws.iter_rows(min_row=1, values_only=True))
    g7820_syms = set(str(r[0]).upper() for r in rr[1:])
    gse782_avail = {s for s in panel if s in g7820_syms}
    common = [s for s in panel
              if s in imv_avail and s in gse165_avail and s in gse782_avail]

    M, sids = rec.load_exmat()
    idx = [row_of_sym[s] for s in common]
    M322 = M[idx, :]
    pd_rows = rec.load_pdata()
    pmap = {r["X"]: r for r in pd_rows}
    y_imv = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                      else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                      else np.nan for s in sids], dtype=float)
    m_lab = ~np.isnan(y_imv)
    Xlab_imv = M322[:, m_lab]
    ylab_imv = y_imv[m_lab]

    Xz_tr, mu, sd = rec.clr_standardize(Xlab_imv)
    pca = rec.fit_pca(Xz_tr)
    S0_imv = pca.transform(Xz_tr)

    zero_all = (M322 == 0).all(1)
    anchor_info = dict(
        key=("exmat row j == the j-th fData row under upper-cased symbol "
             "alphabetical sort (_imvigor_gene_recovery recovered key)"),
        panel=322,
        imv_available=len(imv_avail),
        imv_missing=sorted(set(panel) - imv_avail),
        gse165252_available=len(gse165_avail),
        gse78220_available=len(gse782_avail),
        common=len(common),
        n_allzero_genes_across_348=int(zero_all.sum()),
        pc_explained_variance_ratio=[float(pca.explained_variance_ratio_[0]),
                                     float(pca.explained_variance_ratio_[1])],
    )

    D = {"IMvigor210": dict(s0=S0_imv, y=ylab_imv,
                            n=int(m_lab.sum()), n_event=int(ylab_imv.sum()))}
    M_ext = {"GSE165252": rec.load_gse165252(common, sym2ens),
             "GSE78220": rec.load_gse78220(common)}
    for c, (M322c, yc) in M_ext.items():
        mc = ~np.isnan(yc)
        S0c = rec.transform_pca(M322c[:, mc], mu, sd, pca)
        D[c] = dict(s0=S0c, y=yc[mc], n=int(mc.sum()), n_event=int(yc[mc].sum()))
    return D, common, anchor_info


# ================================================================ within references
def within_reference(D):
    """Three-cohort within reference, /-recovery convention verbatim:
    Efron 0.632 optimism-corrected AUC (rec.auc_632, n_boot=600, seed 20260911)
    with the GLOBAL ridge lambda 1/N_total -- the exact estimator behind the
    frozen new-anchor within baseline."""
    N_total = sum(D[c]["n"] for c in COHORTS)
    lam = 1.0 / N_total
    out = {}
    for c in COHORTS:
        app, corr, se, ci = rec.auc_632(D[c]["s0"], D[c]["y"], lam)
        out[c] = dict(
            n=D[c]["n"], n_event=D[c]["n_event"],
            lambda_ridge=lam,
            auc_apparent=float(app), auc_632=float(corr),
            se=float(se), ci95=[float(ci[0]), float(ci[1])],
            estimator=("/-recovery verbatim: Efron 0.632 optimism-corrected "
                       "AUC (app - mean(app_boot - auc_boot-on-full)), n_boot=600, "
                       "seed 20260911, global lambda 1/N_total"),
            n_boot=N_BOOT_632, seed=SEED)
    return out, N_total, lam


# ================================================================ LOCO benchmark ( verbatim)
def loco_benchmark(D):
    res_all, adjud_all, m2_pred = {}, {}, {}
    for tgt in COHORTS:
        src = [c for c in COHORTS if c != tgt]
        Xs = np.concatenate([D[c]["s0"] for c in src], 0)
        ys = np.concatenate([D[c]["y"] for c in src], 0)
        coh = np.concatenate([[c] * D[c]["n"] for c in src])
        Xt, yt = D[tgt]["s0"], D[tgt]["y"]
        n_src = len(ys)
        lam = 1.0 / n_src
        res = {}

        # ---- M0 case-mix-only (degenerates to prevalence-only; disclosed) ----
        p0 = float(ys.mean())
        res["M0"] = metrics_block(yt, np.full(len(yt), p0))

        # ---- M1 single-feature (PC1 only) ----
        w1 = fit_logistic(np.column_stack([np.ones(n_src), Xs[:, 0]]), ys, lam)
        res["M1"] = metrics_block(yt, sigmoid(w1[0] + Xt[:, 0] * w1[1]))

        # ---- M2 frozen S0 (verbatim  F4 /  baseline algebra) ----
        w2 = fit_logistic(np.column_stack([np.ones(n_src), Xs]), ys, lam)
        eta2_t = w2[0] + Xt @ w2[1:]
        res["M2"] = metrics_block(yt, sigmoid(eta2_t))
        m2_pred[tgt] = dict(eta=eta2_t, p=sigmoid(eta2_t),
                            train="+".join(src), ridge_lambda=lam)

        # ---- M3 supervised target-domain recalibration ----
        res["M3"] = recal_repeats(eta2_t, yt)

        # ---- M4 hierarchical (shared slopes + shrunken cohort intercepts) ----
        a_c = {}
        for c in src:
            m = coh == c
            a_c[c] = float(fit_logistic_gen(
                np.ones((int(m.sum()), 1)), ys[m], 1e-6,
                offset=Xs[m] @ w2[1:], penalize=[False])[0])
        a_cs = {c: w2[0] + M4_GAMMA * (a_c[c] - w2[0]) for c in src}
        off_slope = np.array([a_cs[c] for c in coh])
        b_sh = fit_logistic_gen(Xs, ys, lam, offset=off_slope, penalize=[True, True])
        a_t = float(fit_logistic_gen(
            np.ones((n_src, 1)), ys, 1e-6, offset=Xs @ b_sh, penalize=[False])[0])
        res["M4"] = metrics_block(yt, sigmoid(a_t + Xt @ b_sh))
        res["M4"]["m4_internals"] = dict(
            pooled_intercept=float(w2[0]), source_intercepts_raw=a_c,
            source_intercepts_shrunk=a_cs, shared_slopes=[float(x) for x in b_sh],
            target_intercept=a_t, shrinkage_gamma=M4_GAMMA)

        # ---- M5a CORAL (unsupervised) ----
        res["M5a"] = metrics_block(yt, sigmoid(w2[0] + coral_align(Xt, Xs) @ w2[1:]))

        # ---- M5b location-scale (ComBat-style, unsupervised) ----
        Xt_ls = locscale_align(Xt, Xs)
        res["M5b"] = metrics_block(yt, sigmoid(w2[0] + Xt_ls @ w2[1:]))

        # ---- M6 proposed audit composite ----
        eta6 = w2[0] + Xt_ls @ w2[1:]
        res["M6"] = recal_repeats(eta6, yt)

        res_all[tgt] = res

        # ---- corrected two-criterion adjudication (direction-aware,  verbatim) ----
        for mk, mres in res.items():
            if mres.get("not_executable"):
                adjud_all.setdefault(mk, {})[tgt] = dict(
                    transport_permitted=False, reason="not_executable")
                continue
            lo, hi = mres["auc_ci95"]
            auc_ok = bool(lo is not None and not np.isnan(lo) and lo > 0.5)
            auc_adv = bool(hi is not None and not np.isnan(hi) and hi < 0.5)
            sl = mres["calib_slope"]
            slope_ok = bool(sl is not None and not np.isnan(sl)
                            and abs(sl - 1.0) <= SLOPE_TOL)
            failed = []
            if not auc_ok:
                failed.append("auc_ci_excludes_0p5_in_beneficial_direction")
            if not slope_ok:
                failed.append("calibration_slope_within_tolerance")
            adjud_all.setdefault(mk, {})[tgt] = dict(
                auc=mres["auc"], auc_ci95=[lo, hi],
                auc_ci_excludes_0p5_beneficial_direction=auc_ok,
                auc_ci_excludes_0p5_adverse_direction=auc_adv,
                calib_slope=sl, slope_within_tolerance=slope_ok,
                criteria_failed=failed,
                transport_permitted=bool(auc_ok and slope_ok))
    return res_all, adjud_all, m2_pred


# ================================================================ calibration curves ( supporting analysis-1 machinery)
def calibration_block(D, m2_pred, res_all):
    results = {}
    xchecks = {}
    for tgt in COHORTS:
        yt = D[tgt]["y"]
        eta = m2_pred[tgt]["eta"]
        p2 = m2_pred[tgt]["p"]
        rm2 = res_all[tgt]["M2"]

        m2_bins = sup.equal_freq_bins(yt, p2, 5)
        m2 = dict(
            model=("M2 = frozen S0 ridge (both PCs) trained on the pooled other two "
                   "cohorts, verbatim  F4 /  algebra (lambda=1/n_train); "
                   "new  recovered-key anchor"),
            train=m2_pred[tgt]["train"], ridge_lambda=m2_pred[tgt]["ridge_lambda"],
            n_test=int(len(yt)), n_event=int(yt.sum()),
            auc=rm2["auc"], auc_ci95=rm2["auc_ci95"],
            calib_intercept=rm2["calib_intercept"], calib_slope=rm2["calib_slope"],
            log_loss=rm2["log_loss"], brier=rm2["brier"],
            overall_mean_predicted=float(np.mean(p2)),
            overall_observed_rate=float(np.mean(yt)),
            predicted_quantiles={str(q): float(np.quantile(p2, q))
                                 for q in (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0)},
            per_patient_predicted=p2.tolist(),
            per_patient_y=yt.tolist(),
            bins=m2_bins, r_arrays=sup.r_arrays(m2_bins),
        )

        # ---- M3 (100 stratified 50/50 recalibration repeats; representative fold) ----
        m3r = sup.recal_repeats_(eta, yt)
        lls = np.array(m3r["_per_ll"], float)
        med = float(np.median(lls))
        r_star = int(np.argmin(np.abs(lls - med)))
        ye_s, pe_s = m3r["_per_eval"][r_star]
        fa, fb = m3r["_per_fit_ab"][r_star]
        a3, b3 = calib_ab(ye_s, pe_s)
        m3_bins = sup.equal_freq_bins(ye_s, pe_s, 5)
        m3 = dict(
            model=("M3 = M2 linear predictor eta + target-internal stratified 50/50 "
                   "supervised logistic recalibration logit(p)=a+b*eta, x 100 repeats "
                   "(seed 20260911), verbatim / protocol"),
            uses_target_labels=True,
            repeat_means=dict(
                log_loss=m3r["log_loss"], log_loss_ci95=m3r["log_loss_ci95"],
                brier=m3r["brier"], brier_ci95=m3r["brier_ci95"],
                calib_intercept=m3r["calib_intercept"],
                calib_intercept_ci95=m3r["calib_intercept_ci95"],
                calib_slope=m3r["calib_slope"], calib_slope_ci95=m3r["calib_slope_ci95"],
                auc=m3r["auc"], auc_ci95=m3r["auc_ci95"],
                mean_recal_a=m3r["mean_recal_a"], mean_recal_b=m3r["mean_recal_b"],
                n_repeats_used=m3r["n_repeats_used"],
                n_repeats_skipped=m3r["n_repeats_skipped"],
                n_eval_half_mean=m3r["n_eval_half_mean"]),
            representative_fold=dict(
                rule=("repeat whose eval-half log-loss is closest to the median of "
                      "the 100 eval-half log-losses (ties -> smallest index)"),
                repeat_index=r_star + 1,
                eval_half_log_loss=float(lls[r_star]),
                median_eval_half_log_loss=med,
                fit_half_recal_a=fa, fit_half_recal_b=fb,
                n_eval=int(len(ye_s)), n_event_eval=int(ye_s.sum()),
                auc=float(roc_auc_score(ye_s, pe_s)),
                calib_intercept=a3, calib_slope=b3),
            bins=m3_bins, r_arrays=sup.r_arrays(m3_bins),
            overall_mean_predicted=float(np.mean(pe_s)),
            overall_observed_rate=float(np.mean(ye_s)),
            note=("The representative fold supplies the plotted curve; inference "
                  "rests on the repeat means (/ frozen convention)."),
        )

        # ---- internal cross-check: sup.recal_repeats_ must reproduce the LOCO
        #      recal_repeats repeat means bit-for-bit (identical rng consumption) ----
        lo_m3 = res_all[tgt]["M3"]
        xc = dict(
            calib_slope=abs(m3r["calib_slope"] - lo_m3["calib_slope"]),
            calib_intercept=abs(m3r["calib_intercept"] - lo_m3["calib_intercept"]),
            log_loss=abs(m3r["log_loss"] - lo_m3["log_loss"]),
            auc=abs(m3r["auc"] - lo_m3["auc"]),
            mean_recal_a=abs(m3r["mean_recal_a"] - lo_m3["mean_recal_a"]),
            mean_recal_b=abs(m3r["mean_recal_b"] - lo_m3["mean_recal_b"]))
        worst = max(xc.values())
        assert worst < 1e-12, f"[calib {tgt}] M3 dual-path mismatch: {xc}"
        xchecks[tgt] = dict(max_abs_diff=worst, fields=xc)

        results[tgt] = dict(
            M2=m2, M3=m3,
            recal_before_after=dict(
                m2_calib_slope_logit=rm2["calib_slope"],
                m3_repeat_mean_calib_slope_logit=m3r["calib_slope"],
                m3_calib_slope_ci95=m3r["calib_slope_ci95"],
                slope_improved=bool(abs(m3r["calib_slope"] - 1.0)
                                    < abs(rm2["calib_slope"] - 1.0))))
        print(f"[calib] {tgt}: M2 slope(logit)={rm2['calib_slope']:.3f} -> "
              f"M3 repeat-mean slope={m3r['calib_slope']:.3f} "
              f"CI=[{m3r['calib_slope_ci95'][0]:.3f},{m3r['calib_slope_ci95'][1]:.3f}] "
              f"| M3 dual-path xcheck max|d|={worst:.1e}")
    return results, xchecks


# ================================================================ pairwise transports ( verbatim)
def pairwise_transports(D):
    pairwise = {}
    for s in COHORTS:
        for t in COHORTS:
            if s == t:
                continue
            Xs, ys = D[s]["s0"], D[s]["y"]
            Xt, yt = D[t]["s0"], D[t]["y"]
            w = fit_logistic(np.column_stack([np.ones(len(ys)), Xs]), ys, 1.0 / len(ys))
            p2t = sigmoid(w[0] + Xt @ w[1:])
            pairwise[f"train_{s}__test_{t}"] = dict(
                train_cohort=s, test_cohort=t,
                n_train=int(len(ys)), n_test=int(len(yt)),
                train_prevalence=float(ys.mean()),
                lambda_rule="1/n_train (single-source training; identical algebra)",
                M0=metrics_block(yt, np.full(len(yt), float(ys.mean()))),
                M2=metrics_block(yt, p2t),
                M2_calib_slope_probability_scale__convention=prob_scale_slope(yt, p2t))
    return pairwise


# ================================================================ near-miss analysis
def near_miss_analysis(res_all, pairwise):
    cases = []
    for t in COHORTS:
        for mk in MODEL_ORDER:
            if mk == "M0":
                continue
            r = res_all[t][mk]
            if r.get("not_executable"):
                continue
            lo, hi = r["auc_ci95"]
            sl = r["calib_slope"]
            auc_ok = bool(lo is not None and not np.isnan(lo) and lo > 0.5)
            slope_ok = bool(sl is not None and not np.isnan(sl)
                            and abs(sl - 1.0) <= SLOPE_TOL)
            cases.append(dict(
                scope="LOCO_benchmark", model=mk, target=t,
                train="+".join([c for c in COHORTS if c != t]),
                auc=r["auc"], auc_ci95=[lo, hi], calib_slope=sl,
                auc_pass=auc_ok, slope_pass=slope_ok,
                permitted=bool(auc_ok and slope_ok)))
    for key, v in pairwise.items():
        r = v["M2"]
        lo, hi = r["auc_ci95"]
        sl = r["calib_slope"]
        slp = v["M2_calib_slope_probability_scale__convention"]
        auc_ok = bool(lo is not None and not np.isnan(lo) and lo > 0.5)
        slope_ok = bool(sl is not None and not np.isnan(sl)
                        and abs(sl - 1.0) <= SLOPE_TOL)
        cases.append(dict(
            scope="true_external_pairwise", model="M2",
            target=v["test_cohort"], train=v["train_cohort"],
            auc=r["auc"], auc_ci95=[lo, hi], calib_slope=sl,
            calib_slope_probability_scale__convention=slp,
            auc_pass=auc_ok, slope_pass=slope_ok,
            permitted=bool(auc_ok and slope_ok)))

    # ---- full confirmation table (INCLUDING M0, the degenerate constant model):
    #      every direction x model, is the AUC CI excluding 0.5 in the beneficial
    #      direction? ----
    confirmation = []
    for t in COHORTS:
        for mk in MODEL_ORDER:
            r = res_all[t][mk]
            if r.get("not_executable"):
                continue
            lo, hi = r["auc_ci95"]
            confirmation.append(dict(
                scope="LOCO_benchmark", model=mk, target=t,
                train="+".join([c for c in COHORTS if c != t]),
                auc=r["auc"], auc_ci95=[lo, hi],
                ci_excludes_0p5_beneficial=bool(
                    lo is not None and not np.isnan(lo) and lo > 0.5),
                ci_excludes_0p5_adverse=bool(
                    hi is not None and not np.isnan(hi) and hi < 0.5)))
    for key, v in pairwise.items():
        for mk in ("M0", "M2"):
            r = v[mk]
            lo, hi = r["auc_ci95"]
            confirmation.append(dict(
                scope="true_external_pairwise", model=mk,
                target=v["test_cohort"], train=v["train_cohort"],
                auc=r["auc"], auc_ci95=[lo, hi],
                ci_excludes_0p5_beneficial=bool(
                    lo is not None and not np.isnan(lo) and lo > 0.5),
                ci_excludes_0p5_adverse=bool(
                    hi is not None and not np.isnan(hi) and hi < 0.5)))

    auc_pass_slope_fail = [c for c in cases if c["auc_pass"] and not c["slope_pass"]]
    slope_pass_auc_fail = [c for c in cases if c["slope_pass"] and not c["auc_pass"]]
    adverse = [c for c in cases
               if c["auc_ci95"][1] is not None and c["auc_ci95"][1] < 0.5]
    failing = [c for c in cases if not c["permitted"]]
    closest = None
    if failing:
        closest = max(failing,
                      key=lambda c: (c["auc_ci95"][0] if c["auc_ci95"][0] is not None
                                     else -1.0))
    closest_by_slope = None
    if failing:
        def _slope_dist(c):
            sl = c.get("calib_slope")
            return abs(sl - 1.0) if sl is not None else float("inf")
        closest_by_slope = min(failing, key=_slope_dist)
    beneficial = [c for c in confirmation if c["ci_excludes_0p5_beneficial"]]
    return dict(
        note=("near-miss = configuration that clears one licensing criterion but not "
              "the other; the full confirmation table includes M0 (constant "
              "predictions, degenerate CI [0.5,0.5]) and answers the question: does "
              "ANY direction x model have an AUC CI excluding 0.5 in the beneficial "
              "direction under the new anchor numbers?"),
        any_auc_ci_excludes_0p5_in_beneficial_direction=bool(beneficial),
        n_full_confirmation_entries=len(confirmation),
        beneficial_exclusion_entries=beneficial,
        adverse_exclusion_entries=[c for c in confirmation
                                   if c["ci_excludes_0p5_adverse"]],
        closest_to_passing_overall=closest,
        closest_to_passing_by_slope_distance=(
            dict(closest_by_slope, slope_distance_from_1=abs(closest_by_slope["calib_slope"] - 1.0))
            if closest_by_slope is not None and closest_by_slope.get("calib_slope") is not None
            else None),
        auc_passes_but_slope_fails=auc_pass_slope_fail,
        slope_passes_but_auc_fails=slope_pass_auc_fail,
        full_confirmation_table=confirmation)


# ================================================================ old vs new comparison
def old_vs_new(res_all, within_new, within_conv, pairwise):
     = json.load(open(V75_JSON))
    ps = json.load(open(V71_PS_JSON))

    loco_cmp = {}
    max_auc = dict(abs_delta=0.0, where=None, old=None, new=None)
    max_slope = dict(abs_delta=0.0, where=None, old=None, new=None)
    for t in COHORTS:
        for mk in MODEL_ORDER:
            old = ["results"]["per_target_cohort"][t][mk]
            new = res_all[t][mk]
            d_auc = abs(new["auc"] - old["auc"])
            if d_auc > max_auc["abs_delta"]:
                max_auc = dict(abs_delta=d_auc, where=f"LOCO {mk} @ {t}",
                               old=old["auc"], new=new["auc"])
            o_sl, n_sl = old["calib_slope"], new["calib_slope"]
            if o_sl is not None and n_sl is not None and not np.isnan(n_sl):
                d_sl = abs(n_sl - o_sl)
                if d_sl > max_slope["abs_delta"]:
                    max_slope = dict(abs_delta=d_sl, where=f"LOCO {mk} @ {t}",
                                     old=o_sl, new=n_sl)
            loco_cmp.setdefault(t, {})[mk] = dict(
                old_auc=old["auc"], new_auc=new["auc"],
                delta_auc=new["auc"] - old["auc"],
                old_auc_ci95=old["auc_ci95"], new_auc_ci95=new["auc_ci95"],
                old_calib_slope=o_sl, new_calib_slope=n_sl,
                old_log_loss=old["log_loss"], new_log_loss=new["log_loss"],
                old_brier=old["brier"], new_brier=new["brier"])

    oldA = ps["per_scheme_per_cohort"]["A_anchor_IMV"]
    oldW = ["within_cohort_same_pipeline"]["per_cohort"]
    within_cmp = {}
    max_w = dict(abs_delta=0.0, where=None, old=None, new=None)
    for c in COHORTS:
        d = within_new[c]["auc_632"] - oldA[c]["auc_632_corrected"]
        if abs(d) > max_w["abs_delta"]:
            max_w = dict(abs_delta=abs(d), where=f"within 0.632 @ {c}",
                         old=oldA[c]["auc_632_corrected"], new=within_new[c]["auc_632"])
        within_cmp[c] = dict(
            old__schemeA_convention=dict(
                auc_apparent=oldA[c]["auc_apparent"],
                auc_632=oldA[c]["auc_632_corrected"], ci95=oldA[c]["auc_ci95"]),
            new__recovery_convention=dict(
                auc_apparent=within_new[c]["auc_apparent"],
                auc_632=within_new[c]["auc_632"], ci95=within_new[c]["ci95"]),
            delta_auc_632=d,
            old__efron_oob_convention=dict(
                efron_0632_auc=oldW[c]["efron_0632_auc"],
                apparent_auc=oldW[c]["apparent_auc"]),
            new__efron_oob_convention=dict(
                efron_0632_auc=within_conv[c]["efron_0632_auc"],
                apparent_auc=within_conv[c]["apparent_auc"]),
            delta_efron_oob=within_conv[c]["efron_0632_auc"]
                            - oldW[c]["efron_0632_auc"])

    pw_cmp = {}
    max_p = dict(abs_delta=0.0, where=None, old=None, new=None)
    for k, v in pairwise.items():
        old = ["true_external_pairwise_transports"]["transports"][k]
        d = abs(v["M2"]["auc"] - old["M2"]["auc"])
        if d > max_p["abs_delta"]:
            max_p = dict(abs_delta=d, where=f"pairwise {k} M2",
                         old=old["M2"]["auc"], new=v["M2"]["auc"])
        pw_cmp[k] = dict(
            old_auc=old["M2"]["auc"], new_auc=v["M2"]["auc"],
            delta_auc=v["M2"]["auc"] - old["M2"]["auc"],
            old_auc_ci95=old["M2"]["auc_ci95"], new_auc_ci95=v["M2"]["auc_ci95"],
            old_prob_scale_slope=old["M2_calib_slope_probability_scale__convention"],
            new_prob_scale_slope=v["M2_calib_slope_probability_scale__convention"])

    all_max = [max_auc, max_slope, max_w, max_p]
    overall = max(all_max, key=lambda x: x["abs_delta"])
    overall = dict(what=dict(loco_auc=max_auc, loco_calib_slope=max_slope,
                             within_auc632=max_w, pairwise_auc=max_p),
                   largest=overall)
    return dict(loco=loco_cmp, within=within_cmp, pairwise=pw_cmp,
                max_changes=overall,
                old_sources=dict(_benchmark=V75_JSON, _pca_stability=V71_PS_JSON))


# ================================================================ self-checks
def run_self_checks(anchor_info, within_new, pairwise):
    """Hard checks vs the frozen new-anchor baselines (tolerance 1e-9).
    Returns the check block; ABORTS (sys.exit) on any mismatch."""
    issues = []
    checks = {}

    # ---- anchor identity cross-check vs the frozen recovery JSON ----
    recj = json.load(open(V82_RECOVERY_JSON))["phase2"]["new"]
    gs = recj["gene_set"]
    anchor_cks = dict(
        common=anchor_info["common"] == gs["common"],
        imv_available=anchor_info["imv_available"] == gs["imv_available"],
        imv_missing=anchor_info["imv_missing"] == gs["imv_missing"],
        n_allzero=anchor_info["n_allzero_genes_across_348"]
                  == recj["zero_genes"]["n_allzero_across_348"],
        pc_evr_1=abs(anchor_info["pc_explained_variance_ratio"][0]
                     - recj["pc_explained_variance_ratio"][0]) < TOL,
        pc_evr_2=abs(anchor_info["pc_explained_variance_ratio"][1]
                     - recj["pc_explained_variance_ratio"][1]) < TOL)
    checks["anchor_vs_recovery_json"] = dict(
        fields=anchor_cks, all_pass=all(anchor_cks.values()))
    if not all(anchor_cks.values()):
        issues.append(f"anchor mismatch vs recovery JSON: {anchor_cks}")

    # ---- within IMvigor210 (frozen) ----
    fw = FROZEN["within_imvigor210"]
    wn = within_new["IMvigor210"]
    w_cks = dict(
        apparent=abs(wn["auc_apparent"] - fw["auc_apparent"]),
        auc_632=abs(wn["auc_632"] - fw["auc_632"]),
        se=abs(wn["se"] - fw["se"]),
        ci_lo=abs(wn["ci95"][0] - fw["ci95"][0]),
        ci_hi=abs(wn["ci95"][1] - fw["ci95"][1]))
    checks["within_imvigor210"] = dict(
        frozen=fw, recomputed=dict(
            apparent=wn["auc_apparent"], auc_632=wn["auc_632"],
            se=wn["se"], ci95=wn["ci95"]),
        abs_diff=w_cks, max_abs_diff=max(w_cks.values()),
        all_pass=max(w_cks.values()) < TOL)
    if max(w_cks.values()) >= TOL:
        issues.append(f"within IMvigor210 mismatch: {w_cks}")

    # ---- F1-style (train IMvigor210 -> external), frozen ----
    f1 = {}
    for c in ("GSE165252", "GSE78220"):
        ff = FROZEN["f1_style"][c]
        v = pairwise[f"train_IMvigor210__test_{c}"]
        m2 = v["M2"]
        cks = dict(
            auc=abs(m2["auc"] - ff["auc"]),
            ci_lo=abs(m2["auc_ci95"][0] - ff["ci95"][0]),
            ci_hi=abs(m2["auc_ci95"][1] - ff["ci95"][1]),
            prob_scale_slope=abs(
                v["M2_calib_slope_probability_scale__convention"]
                - ff["calib_slope_probability_scale"]),
            brier=abs(m2["brier"] - ff["brier"]),
            log_loss=abs(m2["log_loss"] - ff["log_loss"]))
        f1[c] = dict(frozen=ff,
                     recomputed=dict(
                         auc=m2["auc"], auc_ci95=m2["auc_ci95"],
                         calib_slope_probability_scale=(
                             v["M2_calib_slope_probability_scale__convention"]),
                         brier=m2["brier"], log_loss=m2["log_loss"]),
                     abs_diff=cks, max_abs_diff=max(cks.values()),
                     all_pass=max(cks.values()) < TOL)
        if max(cks.values()) >= TOL:
            issues.append(f"F1-style {c} mismatch: {cks}")
    checks["f1_style_train_IMvigor210"] = f1

    checks["all_pass"] = not issues
    if issues:
        print("\n!! SELF-CHECK FAILURE -- aborting, numbers NOT written:")
        for s in issues:
            print("   -", s)
        sys.exit(1)
    return checks


# ================================================================ markdown report
def _f(x, nd=3):
    if x is None:
        return "NA"
    if isinstance(x, float) and (np.isnan(x) or np.isinf(x)):
        return "NA"
    return f"{x:.{nd}f}"


def write_md(out):
    A = []
    A.append("#  锚定臂全量重跑：M0–M6 竞争基准 + within 参照 + 校准曲线 + M3 再校准 + pairwise 运输\n")
    A.append(f"日期：{out['date']} ｜ seed {out['seed']} ｜ 脚本：scripts/analysis/rerun_benchmark.py ｜ "
             f"冻结 JSON：data/rerun_benchmark.json")
    A.append("")
    A.append("**锚定管线（唯一权威）**：IMvigor210 基因行身份已恢复——exmat 第 j 行 == fData 按 symbol（大写归一）"
             "字母序排序后的第 j 个基因（scripts/analysis/imvigor_gene_recovery.py）；面板 = 322 清单按恢复钥匙直配 "
             f"（IMv 可得 {out['anchor']['imv_available']}/322，缺 {', '.join(out['anchor']['imv_missing'])}；"
             f"三队列交集 {out['anchor']['common']} 基因）；IMvigor210 面板内全零基因数 = "
             f"{out['anchor']['n_allzero_genes_across_348']}（旧管线为 70）。"
             "旧“位置对齐 + 假 entrez 选择”结果全部作废，本文件为重跑后的新数字。\n")

    # ---- self-checks ----
    sc = out["self_check"]
    A.append("## 1. 自检对齐（冻结新锚定基线，容差 1e-9）\n")
    A.append("| 检查项 | 冻结值 | 重算值 | max |Δ| | 判定 |")
    A.append("|---|---|---|---|---|")
    w = sc["within_imvigor210"]
    A.append(f"| within IMvigor210 0.632 | {w['frozen']['auc_632']:.10f} | "
             f"{w['recomputed']['auc_632']:.10f} | {w['max_abs_diff']:.2e} | "
             f"{'PASS' if w['all_pass'] else 'FAIL'} |")
    for c, v in sc["f1_style_train_IMvigor210"].items():
        A.append(f"| F1-style {c}（AUC/CI/prob-slope/Brier/log-loss 全字段） | "
                 f"AUC {v['frozen']['auc']:.10f} | AUC {v['recomputed']['auc']:.10f} | "
                 f"{v['max_abs_diff']:.2e} | {'PASS' if v['all_pass'] else 'FAIL'} |")
    A.append(f"| 锚定身份（common/missing/零基因/PC EVR vs 恢复 JSON） | — | — | — | "
             f"{'PASS' if sc['anchor_vs_recovery_json']['all_pass'] else 'FAIL'} |")
    A.append(f"\n**全部自检 {'PASS' if sc['all_pass'] else 'FAIL'}**（不一致即中止，未产出任何数字）。\n")

    # ---- anchor diagnostics ----
    an = out["anchor"]
    A.append("## 2. 锚定诊断\n")
    A.append(f"- 面板 322 → IMv 可得 {an['imv_available']}（缺 {', '.join(an['imv_missing'])}）；"
             f"GSE165252 可得 {an['gse165252_available']}；GSE78220 可得 {an['gse78220_available']}；"
             f"三队列交集 **{an['common']}** 基因（旧假 entrez 交集 228）。")
    A.append(f"- 面板内全零基因（348 样本）= {an['n_allzero_genes_across_348']}（旧 70）；"
             f"PC1/PC2 解释方差比 = {100*an['pc_explained_variance_ratio'][0]:.1f}% / "
             f"{100*an['pc_explained_variance_ratio'][1]:.1f}%。")
    A.append(f"- 样本：IMvigor210 {out['cohorts']['IMvigor210']['n']}（事件 "
             f"{out['cohorts']['IMvigor210']['n_event']}）、GSE165252 {out['cohorts']['GSE165252']['n']}"
             f"（{out['cohorts']['GSE165252']['n_event']}）、GSE78220 {out['cohorts']['GSE78220']['n']}"
             f"（{out['cohorts']['GSE78220']['n_event']}）；N_total = {out['within_reference']['N_total']}，"
             f"within λ = 1/N_total = {out['within_reference']['lambda_ridge']:.6f}。\n")

    # ---- LOCO tables ----
    A.append("## 3. M0–M6 LOCO 竞争基准（新锚定， 口径）\n")
    for t in COHORTS:
        A.append(f"### target {t}（train = {'+'.join(c for c in COHORTS if c != t)}）\n")
        A.append("| model | AUC [95% CI] | calib slope b | calib intercept a | Brier | log-loss |")
        A.append("|---|---|---|---|---|---|")
        for mk in MODEL_ORDER:
            r = out["loco_benchmark"]["per_target"][t][mk]
            if r.get("not_executable"):
                A.append(f"| {mk} | NOT EXECUTABLE | | | | |")
                continue
            A.append(f"| {mk} | {_f(r['auc'])} [{_f(r['auc_ci95'][0])},{_f(r['auc_ci95'][1])}] "
                     f"| {_f(r['calib_slope'], 2)} | {_f(r['calib_intercept'], 2)} "
                     f"| {_f(r['brier'], 4)} | {_f(r['log_loss'], 4)} |")
        A.append("")

    # ---- adjudication ----
    A.append("## 4. 两判据判定总表（AUC CI 有利方向排除 0.5 + |slope−1|≤0.10）\n")
    A.append("| model | target | AUC [95% CI] | c1: CI 下界>0.5 | c2: \\|slope−1\\|≤0.10 | transport permitted |")
    A.append("|---|---|---|---|---|---|")
    adj = out["adjudication"]["per_model"]
    for mk in MODEL_ORDER:
        for t in COHORTS:
            e = adj[mk][t]
            A.append(f"| {mk} | {t} | {_f(e['auc'])} [{_f(e['auc_ci95'][0])},{_f(e['auc_ci95'][1])}] "
                     f"| {'YES' if e['auc_ci_excludes_0p5_beneficial_direction'] else 'no'} "
                     f"| {'YES' if e['slope_within_tolerance'] else 'no'} "
                     f"| **{'YES' if e['transport_permitted'] else 'no'}** |")
    vc = out["adjudication"]["verdict_counts"]
    A.append(f"\n**汇总：{vc['n_refused']}/{vc['n_total']} 判定 refused 运输"
             f"（{vc['n_permitted']} permitted）**；按模型："
             + "，".join(f"{mk} {vc['per_model_n_permitted'][mk]}/3" for mk in MODEL_ORDER) + "。\n")

    # ---- within ----
    A.append("## 5. within 同管线参照（三队列，Efron 0.632，600 bootstrap，seed 20260911）\n")
    A.append("主口径 = /-recovery 约定（全局 λ=1/N_total，optimism 修正 0.632）——与冻结基线逐位一致；"
             "附  OOB 约定（λ=1/n，0.368·apparent+0.632·OOB）供与旧  within 对照。\n")
    A.append("| cohort | n (ev) | apparent | **0.632 AUC [95% CI]**（recovery 口径） | 旧  Scheme A | 0.632（ OOB 口径） | 旧  OOB |")
    A.append("|---|---|---|---|---|---|---|")
    wc = out["within_reference"]["per_cohort"]
    oc = out["old_vs_new"]["within"]
    for c in COHORTS:
        w = wc[c]
        o = oc[c]
        A.append(f"| {c} | {w['n']} ({w['n_event']}) | {_f(w['auc_apparent'])} | "
                 f"**{_f(w['auc_632'])} [{_f(w['ci95'][0])},{_f(w['ci95'][1])}]** | "
                 f"{_f(o['old__schemeA_convention']['auc_632'])} | "
                 f"{_f(o['new__efron_oob_convention']['efron_0632_auc'])} | "
                 f"{_f(o['old__efron_oob_convention']['efron_0632_auc'])} |")
    A.append("")

    # ---- calibration curves ----
    A.append("## 6. 校准曲线块（M2 冻结预测 per-patient 分布 + M3 再校准前后）\n")
    cal = out["calibration_curves"]["per_target"]
    for t in COHORTS:
        m2 = cal[t]["M2"]
        m3 = cal[t]["M3"]
        rm = m3["repeat_means"]
        q = m2["predicted_quantiles"]
        A.append(f"### target {t}（M2 train = {m2['train']}）\n")
        A.append(f"- **M2 per-patient 预测分布**：min {_f(q['0.0'])} / 中位 {_f(q['0.5'])} / "
                 f"max {_f(q['1.0'])}（5%–95%：{_f(q['0.05'])}–{_f(q['0.95'])}）；"
                 f"总体均值预测 {_f(m2['overall_mean_predicted'])} vs **目标事件率 "
                 f"{_f(m2['overall_observed_rate'])}**（n={m2['n_test']}，事件 {m2['n_event']}）。")
        A.append(f"- M2 logit 尺度 calib：intercept {_f(m2['calib_intercept'], 2)}，"
                 f"slope {_f(m2['calib_slope'], 2)}。")
        A.append(f"- **M3 再校准前后 slope（logit 尺度）**：{_f(cal[t]['recal_before_after']['m2_calib_slope_logit'], 2)}"
                 f" → **{_f(cal[t]['recal_before_after']['m3_repeat_mean_calib_slope_logit'], 2)}**"
                 f"（100 重复均值，CI [{_f(rm['calib_slope_ci95'][0], 2)}, "
                 f"{_f(rm['calib_slope_ci95'][1], 2)}]；斜率改善："
                 f"{'是' if cal[t]['recal_before_after']['slope_improved'] else '否'}）。")
        A.append("- 5 等频 bin（预测→观测，n）：" + "；".join(
            f"{b['mean_predicted']:.3f}→{b['observed_rate']:.3f} (n={b['n']})"
            for b in m2["bins"]))
        A.append("")

    # ---- pairwise ----
    A.append("## 7. pairwise 真外部运输（全部 6 方向：F1–F3 + 源互换敏感性）\n")
    A.append("| 方向（train→test） | M2 AUC [95% CI]（新） | prob-scale slope（新） | 旧 AUC [95% CI] | 旧 prob-slope |")
    A.append("|---|---|---|---|---|")
    pc = out["old_vs_new"]["pairwise"]
    for k in [f"train_{s}__test_{t}" for s in COHORTS for t in COHORTS if s != t]:
        v = pc[k]
        direction = k.replace("train_", "").replace("__", " → ").replace("test_", "")
        A.append(f"| {direction} "
                 f"| {_f(v['new_auc'])} [{_f(v['new_auc_ci95'][0])},{_f(v['new_auc_ci95'][1])}] "
                 f"| {_f(v['new_prob_scale_slope'], 2)} "
                 f"| {_f(v['old_auc'])} [{_f(v['old_auc_ci95'][0])},{_f(v['old_auc_ci95'][1])}] "
                 f"| {_f(v['old_prob_scale_slope'], 2)} |")
    A.append("")

    # ---- near-miss ----
    nm = out["near_miss"]
    A.append("## 8. near-miss 检查（新数字下有无任何方向×模型的 AUC CI 在有利方向排除 0.5）\n")
    if nm["any_auc_ci_excludes_0p5_in_beneficial_direction"]:
        A.append(f"**有 {len(nm['beneficial_exclusion_entries'])} 项**在有利方向排除 0.5：")
        for c0 in nm["beneficial_exclusion_entries"]:
            A.append(f"- {c0['scope']} {c0['model']} {c0['train']}→{c0['target']}: "
                     f"AUC {_f(c0['auc'])} [{_f(c0['auc_ci95'][0])},{_f(c0['auc_ci95'][1])}]")
    else:
        A.append(f"**无**。全表 {nm['n_full_confirmation_entries']} 项（LOCO 8 模型×3 目标 + pairwise 6 方向×M0/M2，"
                 "含退化 M0）逐项确认：没有任何方向×模型的 AUC CI 下界 > 0.5。")
    A.append(f"- GSE78220 F1 方向（train IMvigor210）：AUC 0.660 [0.404, 0.862]——CI 下界 0.404 不过（与冻结基线一致）。")
    if nm["adverse_exclusion_entries"]:
        A.append("- **不利方向**（CI 上界 < 0.5）排除 0.5 的项（同样不满足判据 1 的有利方向要求）：" + "；".join(
            f"{c0['scope']} {c0['model']} {c0['train']}→{c0['target']} "
            f"[{_f(c0['auc_ci95'][0])},{_f(c0['auc_ci95'][1])}]"
            for c0 in nm["adverse_exclusion_entries"]))
    else:
        A.append("- 不利方向（CI 上界 < 0.5）无任何排除。")
    cl = nm["closest_to_passing_overall"]
    if cl:
        A.append(f"- 按 CI 下界最近接通过者：{cl['scope']} {cl['model']} {cl['train']}→{cl['target']}，"
                 f"CI 下界 {_f(cl['auc_ci95'][0], 4)}。")
    cs = nm.get("closest_to_passing_by_slope_distance")
    if cs:
        A.append(f"- 按约束距离最近接通过者（本表瓶颈在 slope 判据）：{cs['scope']} {cs['model']} "
                 f"{cs['train']}→{cs['target']}，slope {_f(cs['calib_slope'], 3)}，"
                 f"|slope−1| = {_f(cs['slope_distance_from_1'], 3)}。")
    if nm["auc_passes_but_slope_fails"]:
        A.append("- **AUC 过、slope 挂**（证明判据 2 必要）：" + "；".join(
            f"{c0['scope']} {c0['model']} {c0['train']}→{c0['target']} "
            f"AUC {_f(c0['auc'])} slope {_f(c0['calib_slope'], 1)}"
            for c0 in nm["auc_passes_but_slope_fails"]))
    else:
        A.append("- 无“AUC 过、slope 挂”的 near-miss 项（AUC 判据在全表均未通过）。")
    if nm["slope_passes_but_auc_fails"]:
        A.append("- **slope 过、AUC 挂**（证明判据 1 必要）：" + "；".join(
            f"{c0['scope']} {c0['model']} {c0['train']}→{c0['target']} "
            f"slope {_f(c0['calib_slope'], 2)} AUC {_f(c0['auc'])}"
            for c0 in nm["slope_passes_but_auc_fails"]))
    else:
        A.append("- 无“slope 过、AUC 挂”的 near-miss 项。")
    A.append("")

    # ---- old vs new ----
    A.append("## 9. 新旧对照总表（旧 = 错误对齐冻结档案；新 = 恢复钥匙重跑）\n")
    A.append("| model | target | 旧 AUC [95% CI] | 新 AUC [95% CI] | ΔAUC | 旧 slope | 新 slope |")
    A.append("|---|---|---|---|---|---|---|")
    lc = out["old_vs_new"]["loco"]
    for t in COHORTS:
        for mk in MODEL_ORDER:
            e = lc[t][mk]
            A.append(f"| {mk} | {t} | {_f(e['old_auc'])} [{_f(e['old_auc_ci95'][0])},{_f(e['old_auc_ci95'][1])}] "
                     f"| {_f(e['new_auc'])} [{_f(e['new_auc_ci95'][0])},{_f(e['new_auc_ci95'][1])}] "
                     f"| {e['delta_auc']:+.4f} | {_f(e['old_calib_slope'], 2)} | {_f(e['new_calib_slope'], 2)} |")
    mc = out["old_vs_new"]["max_changes"]["what"]
    mc_largest = out["old_vs_new"]["max_changes"]["largest"]
    A.append("")
    A.append("### 新旧最大变化\n")
    A.append(f"- LOCO AUC 最大 |Δ| = {mc['loco_auc']['abs_delta']:.4f}"
             f"（{mc['loco_auc']['where']}：{mc['loco_auc']['old']:.4f} → {mc['loco_auc']['new']:.4f}）。")
    A.append(f"- LOCO calib slope 最大 |Δ| = {mc['loco_calib_slope']['abs_delta']:.2f}"
             f"（{mc['loco_calib_slope']['where']}）。")
    A.append(f"- within 0.632 最大 |Δ| = {mc['within_auc632']['abs_delta']:.4f}"
             f"（{mc['within_auc632']['where']}）。")
    A.append(f"- pairwise AUC 最大 |Δ| = {mc['pairwise_auc']['abs_delta']:.4f}"
             f"（{mc['pairwise_auc']['where']}）。")
    A.append(f"- **全表最大变化 = {mc_largest['abs_delta']:.4f}（{mc_largest['where']}）**。\n")

    md = "\n".join(A)
    with open(OUT_MD, "w") as f:
        f.write(md)
    print("Wrote", OUT_MD)


# ================================================================ main
def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ---- 1. anchor with the recovered key ----
    D, common, anchor_info = build_anchor()
    print("[anchor] common:", anchor_info["common"],
          "| missing:", anchor_info["imv_missing"],
          "| zero genes:", anchor_info["n_allzero_genes_across_348"],
          "| PC EVR:", [round(x, 4) for x in anchor_info["pc_explained_variance_ratio"]])
    print("[anchor] cohorts:", {c: (D[c]["n"], D[c]["n_event"]) for c in COHORTS})

    # ---- 2. within references (frozen-checked convention first: fail fast) ----
    within_new, N_total, lam_within = within_reference(D)
    print("[within] recovery-convention 0.632:",
          {c: round(within_new[c]["auc_632"], 4) for c in COHORTS})
    within_conv = efron632_within(D)
    print("[within] -OOB-convention 0.632:",
          {c: round(within_conv[c]["efron_0632_auc"], 4) for c in COHORTS})

    # ---- 3. pairwise transports (needed for the F1-style frozen check) ----
    pairwise = pairwise_transports(D)
    print("[pairwise] M2 AUCs:",
          {k: round(v["M2"]["auc"], 4) for k, v in pairwise.items()})

    # ---- 4. HARD self-checks vs frozen baselines (abort on mismatch) ----
    self_check = run_self_checks(anchor_info, within_new, pairwise)
    print("[self-check] ALL PASS:",
          self_check["all_pass"],
          "| within max|d| =", f"{self_check['within_imvigor210']['max_abs_diff']:.2e}",
          "| F1 max|d| =", max(v["max_abs_diff"]
                               for v in self_check["f1_style_train_IMvigor210"].values()))

    # ---- 5. LOCO benchmark M0-M6 ----
    res_all, adjud_all, m2_pred = loco_benchmark(D)
    print("[loco] done")

    # ---- 6. calibration curves + M3 recalibration ----
    calib, calib_xchecks = calibration_block(D, m2_pred, res_all)

    # ---- 7. near-miss ----
    near = near_miss_analysis(res_all, pairwise)

    # ---- 8. old vs new ----
    ovn = old_vs_new(res_all, within_new, within_conv, pairwise)

    # ---- verdict counts ----
    n_perm = {mk: int(sum(1 for t in COHORTS
                           if adjud_all[mk][t].get("transport_permitted")))
              for mk in MODEL_ORDER}
    n_total = len(MODEL_ORDER) * len(COHORTS)
    n_permitted = int(sum(n_perm.values()))
    verdict_counts = dict(n_total=n_total, n_permitted=n_permitted,
                          n_refused=n_total - n_permitted,
                          per_model_n_permitted=n_perm)

    out = dict(
        analysis="rerun_benchmark",
        date="2026-09-24",
        seed=SEED,
        metadata=dict(
            n_boot_auc_ci=N_BOOT,
            n_boot_efron632=N_BOOT_632,
            n_repeats_recal=100,
            slope_tolerance=SLOPE_TOL,
            anchor_pipeline="scripts/analysis/imvigor_gene_recovery.py (recovered key; sole authority)",
            algebra_imported_verbatim=[
                "_competitor_benchmark (sigmoid/fit_logistic/fit_logistic_gen/calib_ab/"
                "auc_bca/metrics_block/recal_repeats/coral_align/locscale_align)",
                "_competitor_benchmark_ (prob_scale_slope, efron632_within, "
                "two-criterion direction-aware adjudication)",
                "_supporting_analyses (equal_freq_bins, r_arrays, recal_repeats_)",
                "_imvigor_gene_recovery (clr_standardize/fit_pca/transform_pca/"
                "load_gse165252/load_gse78220/auc_632, recovered-key loading)"],
            calibration_convention=("logit-scale logistic recalibration "
                                    "logit(P(y=1)) = a + b*logit(p_hat) for the LOCO "
                                    "benchmark/adjudication (/); the pairwise "
                                    "blocks additionally report the  probability-"
                                    "scale slope, which is what the frozen F1-style "
                                    "baseline uses"),
        ),
        anchor=anchor_info,
        cohorts={c: dict(n=D[c]["n"], n_event=D[c]["n_event"]) for c in COHORTS},
        self_check=self_check,
        within_reference=dict(
            per_cohort=within_new,
            N_total=N_total, lambda_ridge=lam_within,
            _oob_convention=dict(
                role=("secondary reference:  efron632_within convention "
                      "(lambda=1/n per cohort, 0.368*apparent + 0.632*mean OOB)"),
                per_cohort=within_conv)),
        loco_benchmark=dict(
            design=("leave-one-cohort-out (LOCO): for each held-out target cohort, "
                    "train = pooled remaining 2 cohorts, test = ALL labeled patients "
                    "of the target; identical test sets for every model"),
            per_target=res_all),
        adjudication=dict(
            rule=dict(
                c1="external AUC bootstrap 95% CI excludes 0.5 IN THE BENEFICIAL "
                   "direction (CI lower bound > 0.5)",
                c2=f"|calibration slope - 1| <= {SLOPE_TOL} (logit scale)",
                transport_permitted="both criteria met"),
            per_model={mk: adjud_all[mk] for mk in MODEL_ORDER},
            verdict_counts=verdict_counts),
        calibration_curves=dict(
            method=("LOCO M2 frozen predictions, 5 equal-frequency bins with exact "
                    "Clopper-Pearson 95% CI; M3 = 100 stratified 50/50 supervised "
                    "recalibration repeats (repeat means + representative fold); "
                    "_supporting_analyses machinery verbatim"),
            per_target=calib,
            m3_dual_path_cross_check=calib_xchecks),
        pairwise_transports=dict(
            note=("descriptive single-source (true external) transports over all 6 "
                  "ordered directions (F1-F3 and the source-swap sensitivity "
                  "directions); train_IMvigor210__test_* entries are the F1-style "
                  "frozen-checked numbers"),
            transports=pairwise),
        near_miss=near,
        old_vs_new=ovn,
    )

    with open(OUT_JSON, "w") as f:
        json.dump(_clean(out), f, indent=2)
    print("Wrote", OUT_JSON)

    write_md(out)

    # ---- console summary ----
    print("\n================ FINAL SUMMARY ================")
    print(f"adjudication: {verdict_counts['n_refused']}/{verdict_counts['n_total']} "
          f"refused (permitted {n_permitted}); per model:",
          {mk: f"{n_perm[mk]}/3" for mk in MODEL_ORDER})
    mc = ovn["max_changes"]["what"]
    print(f"max old->new change: LOCO AUC |d|={mc['loco_auc']['abs_delta']:.4f} "
          f"({mc['loco_auc']['where']}); within |d|={mc['within_auc632']['abs_delta']:.4f} "
          f"({mc['within_auc632']['where']}); pairwise |d|={mc['pairwise_auc']['abs_delta']:.4f} "
          f"({mc['pairwise_auc']['where']})")
    print(f"largest overall: {ovn['max_changes']['largest']['abs_delta']:.4f} "
          f"({ovn['max_changes']['largest']['where']})")
    print("near-miss: any beneficial-direction CI exclusion:",
          near["any_auc_ci_excludes_0p5_in_beneficial_direction"],
          f"({len(near['beneficial_exclusion_entries'])} entries; "
          f"{len(near['adverse_exclusion_entries'])} adverse-direction)")
    if near["closest_to_passing_overall"]:
        cl = near["closest_to_passing_overall"]
        print(f"  closest to passing: {cl['scope']} {cl['model']} {cl['train']}->{cl['target']} "
              f"CI_lo={cl['auc_ci95'][0]:.4f}")
    print("files:")
    print("  ", OUT_JSON)
    print("  ", OUT_MD)
    print("  ", f"{SCRIPTS}/rerun_benchmark.py")


if __name__ == "__main__":
    main()
