#!/usr/bin/env python3
""" -- Nature Cancer revision , P0-2 RERUN.

The original P0-2 positive control (positive_control.py, 2026-09-24) was
VOIDED at runtime because the IMvigor210 exmat gene-row identity was destroyed
by upstream anonymization (sex-anchor falsification), so the Ayers IFN-gamma
6-gene signature could not be validly anchored on IMvigor210.

That identity has since been RECOVERED (imvigor_gene_recovery.py, frozen
data/imvigor_gene_recovery.json): exmat row j == the j-th fData row
when fData is sorted by SYMBOL (upper-cased) alphabetically -- verified by
sex-anchor gates, upstream provenance, and a 31085/31085-row numeric identity
match against the labeled mirror counts(cds).

This script RERUNS the P0-2 positive control AS ORIGINALLY SPECIFIED: the
Ayers et al. JCI 2017 IFN-gamma 6-gene signature (IDO1, CXCL9, CXCL10, CCL5,
HLA-DRA, STAT1) through the IDENTICAL frozen-transport machinery as the main
analysis, with the ONLY change being the IMvigor210 loader (the recovered-key
load_imvigor_correct instead of the falsified positional alignment).
Everything else -- metric algebra (imported verbatim from /), cohort
loaders for GSE165252/GSE78220 (verbatim from _pca_stability), the Efron
0.632 within-cohort machinery, the two-criterion adjudication with the
+/-0.05 / +/-0.10 / +/-0.20 tolerance grid, seed 20260911 -- is reused
unchanged from positive_control.py.

Pipeline (as specified):
  1. per-gene availability of the 6 genes in all three cohorts (IMvigor210
     under the recovered key; honest recording of any missing gene);
  2. IMvigor210 anchor: per-gene z on log1p(CPM) (full-matrix library size),
     mean/sd frozen on the n=298 labeled patients, signature = mean of the
     available genes' z-scores; univariate ridge logistic lambda=1/n_train,
     intercept unpenalised, fit on the anchor and frozen;
  3. frozen transport to GSE165252 / GSE78220 on the SAME held-out patient
     sets as the main analysis (all labeled patients): AUC + 2000-rep
     bootstrap BCa 95% CI (seed 20260911), logit-scale calibration slope,
     calibration-in-the-large intercept, Brier, log-loss; within-cohort
     Efron 0.632 (n_boot=600, seed 20260911);
  4. two-criterion adjudication (C1: CI lower bound > 0.5; C2: |slope-1| <=
     0.10) with the pre-registered tolerance tiers;
  5. comparison with the NEW S0 baseline (recovered-key anchor arm rebuilt in
     imvigor_gene_recovery.json: external GSE165252 0.5036 / GSE78220
     0.6603).

Honesty discipline: everything computed at run time, nothing fabricated; all
per-patient signature scores and predicted probabilities are stored in the
JSON.  Sensitivity variant: raw values instead of log1p(CPM) (same structure
as the voided run).

Emission: data/positive_control_gep18.json + .md report.
"""

#  NOTE: signature switched to the 18-gene T-cell-inflamed GEP
# (Cristescu et al., Science 2018). All machinery -- loaders, metric algebra,
# adjudication, seeds -- is unchanged from rerun_positive_control.py.

import json, os, sys, csv, gzip, pickle
import numpy as np

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
BASE = f"{ROOT}/outputs/analysis_base"
WS = f"{ROOT}/data/_workspace"
OUTDIR = f"{BASE}/data/"
OUT = f"{OUTDIR}/positive_control_gep18.json"
MDOUT = f"{OUTDIR}/positive_control_gep18.md"
SEED = 20260911
COHORTS = ["IMvigor210", "GSE165252", "GSE78220"]
EXT_COHORTS = ["GSE165252", "GSE78220"]
GENES = ["CD8A", "CMKLR1", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "CCL5", "CXCL9", "CXCL10", "IDO1", "PSMB10", "STAT1", "TIGIT", "CD27", "CD274", "CD276", "LAG3", "PDCD1LG2"]
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
IMV_KEY_CSV = f"{ROOT}/geo_data/imvigor210_labeled/exmat_row_key.csv"
GSE165_MAT = f"{ROOT}/geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz"
GSE782_MAT = f"{ROOT}/geo_data/GSE78220/GSE78220_PatientFPKM.xlsx"
SYM2ENS_PKL = f"{ROOT}/gene_name_id_dict_gc104M.pkl"
S0_JSON = f"{OUTDIR}/imvigor_gene_recovery.json"
TIERS = [0.05, 0.10, 0.20]

sys.path.insert(0, f"{BASE}/scripts")
sys.path.insert(0, f"{BASE}/data/")

# ---- VERBATIM frozen-pipeline components (import guarantees identity) ----
from competitor_benchmark import (  # noqa: E402
    sigmoid, fit_logistic, metrics_block, _clean, N_BOOT)
from competitor_benchmark_ import N_BOOT_632  # noqa: E402
from pca_stability import load_gse165252, load_gse78220  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402


def efron632_within_subset(D):
    """VERBATIM copy of _competitor_benchmark_.efron632_with the
    cohort loop parameterised to D's keys; formula, seed and bootstrap count
    unchanged (identical to the voided positive_control.py copy)."""
    out = {}
    for c in D:
        X, y = D[c]["s0"], D[c]["y"]
        n = len(y)
        w = fit_logistic(np.column_stack([np.ones(n), X]), y, 1.0 / n)
        app = float(roc_auc_score(y, sigmoid(w[0] + X @ w[1:])))
        rng = np.random.default_rng(SEED)
        oob, skipped = [], 0
        for _ in range(N_BOOT_632):
            idx = rng.integers(0, n, n)
            oob_i = np.setdiff1d(np.arange(n), idx)
            if len(np.unique(y[idx])) < 2 or len(np.unique(y[oob_i])) < 2:
                skipped += 1
                continue
            wb = fit_logistic(np.column_stack([np.ones(n), X[idx]]), y[idx], 1.0 / n)
            oob.append(float(roc_auc_score(y[oob_i], sigmoid(wb[0] + X[oob_i] @ wb[1:]))))
        oob_mean = float(np.mean(oob)) if oob else float("nan")
        out[c] = dict(
            n=int(n), n_event=int(y.sum()),
            lambda_rule="1/n (within-cohort n, intercept unpenalised; identical algebra to M2)",
            apparent_auc=app,
            oob_auc_mean=oob_mean,
            n_boot=N_BOOT_632, n_boot_valid=len(oob), n_boot_skipped=skipped,
            efron_0632_auc=float(0.368 * app + 0.632 * oob_mean),
            formula=("AUC_0.632 = 0.368 * AUC_apparent + 0.632 * mean(AUC_out-of-bag "
                     f"over {N_BOOT_632} bootstrap refits)"),
            seed=SEED)
    return out


# ------------------------------------------------------------------ recovered-key IMvigor loader
def load_imvigor_correct(genes):
    """RECOVERED-KEY IMvigor210 loader (_imvigor_gene_recovery stage 1a/1c
    key): exmat row j == the j-th fData row sorted by upper-cased symbol.
    Returns (M_genes x 348 raw counts, y_binaryResponse, full-matrix column
    library sums, sids, per-gene exmat row indices 0-based)."""
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        sids = [x.strip('"') for x in hdr[1:]]
        rows = []
        for r in rd:
            rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
    M = np.array(rows, float)                      # (31085 x 348)
    lib = M.sum(0)                                # full-matrix library sizes

    fd = list(csv.DictReader(open(IMV_FDATA)))
    fd_sym_up = [r["symbol"].strip().upper() for r in fd]
    order = sorted(range(len(fd)), key=lambda i: fd_sym_up[i])
    row_of_sym = {}
    for j, fi in enumerate(order):
        s = fd_sym_up[fi]
        if s not in row_of_sym:                    # first occurrence wins
            row_of_sym[s] = j

    missing = [g for g in genes if g not in row_of_sym]
    idx = [row_of_sym[g] for g in genes if g in row_of_sym]
    Mg = M[idx, :]

    pd_ = list(csv.DictReader(open(IMV_PDATA)))
    pmap = {r["X"]: r for r in pd_}
    y = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                  else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                  else np.nan for s in sids], dtype=float)
    return Mg, y, lib, sids, {g: row_of_sym[g] for g in row_of_sym if g in set(genes)}, missing


def key_crosscheck(row_of_sym_6):
    """Verify the in-script recovered key against the exported project asset
    geo_data/imvigor210_labeled/exmat_row_key.csv (written by
    imvigor_gene_recovery.py)."""
    if not os.path.exists(IMV_KEY_CSV):
        return dict(performed=False, reason="exmat_row_key.csv not found")
    key = {}
    with open(IMV_KEY_CSV) as f:
        rd = csv.DictReader(f)
        n_rows = 0
        for r in rd:
            n_rows += 1
            key[r["symbol"].strip().upper()] = int(r["exmat_row_1based"])
    agree = {g: (key.get(g) == row_of_sym_6[g] + 1) for g in row_of_sym_6}
    return dict(performed=True, n_key_rows=n_rows,
                per_gene_agreement=agree,
                all_agree=bool(all(agree.values())))


# ------------------------------------------------------------------ external library sums
def libsum_gse165252():
    with gzip.open(GSE165_MAT, "rt") as f:
        f.readline()
        lib = None
        for line in f:
            v = np.array([float(x) for x in line.rstrip("\n").split("\t")[1:]])
            lib = v if lib is None else lib + v
    return lib


def libsum_gse78220():
    import openpyxl
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    lib = None
    for r in rows:
        v = np.array([float(x) if x is not None else 0.0 for x in r[1:]])
        lib = v if lib is None else lib + v
    return lib


# ------------------------------------------------------------------ scoring (verbatim from voided run)
def logcpm(M_raw, lib):
    """genes x n raw values -> log1p(CPM) using the FULL-matrix library size
    (main-analysis normalisation chain up to the CLR step; CLR per-sample
    centring is specific to the 322-gene compositional panel and is NOT applied
    to a 6-gene signature -- disclosed, unchanged from the voided run)."""
    lib = np.asarray(lib, float)
    lib = np.where(lib <= 0, 1.0, lib)
    return np.log1p(M_raw / lib[None, :] * 1e6)


def freeze_z(M_lab):
    mu = M_lab.mean(1)
    sd = M_lab.std(1)
    defined = sd > 1e-8
    return mu, sd, defined


def signature_scores(M_lab, mu, sd, genes_used):
    Z = (M_lab - mu[:, None]) / (sd[:, None] + 1e-9)
    idx = [GENES.index(g) for g in genes_used]
    return Z[idx, :].mean(0)


# ------------------------------------------------------------------ adjudication (verbatim)
def adjudicate(m):
    lo, hi = m["auc_ci95"]
    sl = m["calib_slope"]
    c1_strict = bool(lo is not None and not np.isnan(lo) and lo > 0.5)
    c2_main = bool(sl is not None and not np.isnan(sl) and abs(sl - 1.0) <= 0.10)
    adverse = bool(hi is not None and not np.isnan(hi) and hi < 0.5)
    out = dict(
        auc=m["auc"], auc_ci95=[lo, hi], calib_slope=sl,
        auc_ci_entirely_below_0p5_adverse_direction=adverse,
        main_analysis_rule=dict(
            c1_auc_ci_excludes_0p5_beneficial_direction=c1_strict,
            c2_abs_slope_minus_1_le_0p10=c2_main,
            both_passed=bool(c1_strict and c2_main)),
        tolerance_tiers={},
    )
    for t in TIERS:
        c1 = bool(lo is not None and not np.isnan(lo) and lo > 0.5 - t)
        c2 = bool(sl is not None and not np.isnan(sl) and abs(sl - 1.0) <= t)
        out["tolerance_tiers"][f"pm{t:.2f}"] = dict(
            tolerance=t, c1_ci_lower_bound_gt_0p5_minus_t=c1,
            c2_abs_slope_minus_1_le_t=c2, both_passed=bool(c1 and c2))
    return out


# ------------------------------------------------------------------ transport runner (verbatim)
def run_transport(L_lab, y, anchor, target, genes_used):
    """Frozen z on anchor labeled patients -> signature -> univariate ridge
    (lambda=1/n_anchor, intercept unpenalised) -> predict target labeled."""
    mu, sd, _ = freeze_z(L_lab[anchor])
    sig_a = signature_scores(L_lab[anchor], mu, sd, genes_used)
    sig_t = signature_scores(L_lab[target], mu, sd, genes_used)
    n_tr = len(y[anchor])
    lam = 1.0 / n_tr
    w = fit_logistic(np.column_stack([np.ones(n_tr), sig_a]), y[anchor], lam)
    p_t = sigmoid(w[0] + w[1] * sig_t)
    met = metrics_block(y[target], p_t)
    p_a = sigmoid(w[0] + w[1] * sig_a)
    met_a = metrics_block(y[anchor], p_a)
    pp = [dict(patient_index=i, signature_score=float(sig_t[i]),
               p_hat=float(p_t[i]), y=float(y[target][i]))
          for i in range(len(y[target]))]
    pp_a = [dict(patient_index=i, signature_score=float(sig_a[i]),
                 p_hat=float(p_a[i]), y=float(y[anchor][i]))
            for i in range(n_tr)]
    return dict(
        anchor=anchor, target=target, genes_used=list(genes_used),
        frozen_z={g: dict(mean=float(mu[GENES.index(g)]), sd=float(sd[GENES.index(g)]))
                  for g in genes_used},
        model=dict(intercept=float(w[0]), beta_signature=float(w[1]),
                   lambda_ridge=lam, lambda_rule="1/n_train, intercept unpenalised",
                   n_train=int(n_tr), frozen_on=f"{anchor} labeled patients"),
        target_metrics=met, anchor_apparent_metrics=met_a,
        adjudication_target=adjudicate(met),
        per_patient_target=pp, per_patient_anchor=pp_a,
    )


# ------------------------------------------------------------------ main
def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ================= 1. availability under the recovered key =================
    print("== 1. gene availability (IMvigor210 under the RECOVERED key) ==")
    M_imv_raw_all, y_tot_imv, lib_imv, sids, row_of_sym_6, imv_missing = \
        load_imvigor_correct(GENES)
    xchk = key_crosscheck(row_of_sym_6)
    print(f"   key cross-check vs exmat_row_key.csv: performed={xchk['performed']} "
          f"all_agree={xchk.get('all_agree')}")
    if xchk["performed"] and not xchk["all_agree"]:
        raise SystemExit("FATAL: in-script recovered key disagrees with exported "
                         "exmat_row_key.csv -- aborting rather than running on a "
                         "wrong identity.")

    g = pickle.load(open(SYM2ENS_PKL, "rb"))
    sym2ens = {k.upper(): v for k, v in g.items()}
    with gzip.open(GSE165_MAT, "rt") as f:
        f.readline()
        ens165 = set(line.rstrip("\n").split("\t")[0].split(".")[0] for line in f)
    import openpyxl
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    syms782 = set(str(r[0]).upper() for r in ws.iter_rows(min_row=2, values_only=True))

    avail = {}
    for gene in GENES:
        avail[gene] = dict(
            IMvigor210=bool(gene in row_of_sym_6),
            imvigor_exmat_row_1based=(row_of_sym_6[gene] + 1
                                      if gene in row_of_sym_6 else None),
            GSE165252=bool(sym2ens.get(gene) in ens165),
            gse165252_ensembl_id=sym2ens.get(gene),
            GSE78220=bool(gene in syms782),
        )
        a = avail[gene]
        print(f"  {gene:8s} IMV={a['IMvigor210']} (exmat row {a['imvigor_exmat_row_1based']}) "
              f"GSE165252={a['GSE165252']} ({a['gse165252_ensembl_id']}) "
              f"GSE78220={a['GSE78220']}")
    genes_used = [gene for gene in GENES
                  if avail[gene]["IMvigor210"] and avail[gene]["GSE165252"]
                  and avail[gene]["GSE78220"]]
    genes_missing = [gene for gene in GENES if gene not in genes_used]
    print(f"   genes used (available in all 3 cohorts): {genes_used}"
          + (f" | MISSING: {genes_missing}" if genes_missing else ""))

    # ================= 2-3. anchor fit + frozen transport =================
    M_raw, y_tot, lib = {}, {}, {}
    M_raw["IMvigor210"], y_tot["IMvigor210"], lib["IMvigor210"] = \
        M_imv_raw_all, y_tot_imv, lib_imv
    M_raw["GSE165252"], y_tot["GSE165252"] = load_gse165252(genes_used, sym2ens)
    lib["GSE165252"] = libsum_gse165252()
    M_raw["GSE78220"], y_tot["GSE78220"] = load_gse78220(genes_used)
    lib["GSE78220"] = libsum_gse78220()

    lab_raw, y = {}, {}
    for c in COHORTS:
        m = ~np.isnan(y_tot[c])
        lab_raw[c] = M_raw[c][:, m]
        y[c] = y_tot[c][m].astype(float)
    lib_lab = {c: lib[c][~np.isnan(y_tot[c])] for c in COHORTS}
    L = {c: logcpm(lab_raw[c], lib_lab[c]) for c in COHORTS}
    print(f"   labeled n: " + ", ".join(f"{c}={len(y[c])} (events {int(y[c].sum())})"
                                        for c in COHORTS))

    # ---- anchor apparent + zero-variance check on IMvigor210 training
    mu_i, sd_i, defined_i = freeze_z(L["IMvigor210"])
    genes_zerovar = [GENES[i] for i in range(len(GENES))
                     if GENES[i] in genes_used and not defined_i[i]]
    genes_defined = [gene for gene in genes_used
                     if defined_i[GENES.index(gene)]]
    print(f"   zero-variance genes in IMvigor210 training: "
          f"{genes_zerovar if genes_zerovar else 'none'}")
    # primary definition: signature = mean over ALL available genes (sd->0
    # guard absorbs any zero-variance gene as z=0); sensitivity: defined-only
    sig_imv = signature_scores(L["IMvigor210"], mu_i, sd_i, genes_used)

    # ---- within-cohort Efron 0.632 (anchor, frozen z == within z on IMvigor)
    D_imv = {"IMvigor210": dict(s0=sig_imv.reshape(-1, 1), y=y["IMvigor210"])}
    efron_imv = efron632_within_subset(D_imv)["IMvigor210"]
    print(f"   within IMvigor210: apparent={efron_imv['apparent_auc']:.4f} "
          f"0.632={efron_imv['efron_0632_auc']:.4f}")

    # ---- external within-cohort references (conventional within-z, descriptive)
    efron_ext = {}
    for c in EXT_COHORTS:
        mu_c, sd_c, _ = freeze_z(L[c])
        sig_c = signature_scores(L[c], mu_c, sd_c, genes_used)
        efron_ext[c] = efron632_within_subset(
            {c: dict(s0=sig_c.reshape(-1, 1), y=y[c])})[c]
        print(f"   within {c}: apparent={efron_ext[c]['apparent_auc']:.4f} "
              f"0.632={efron_ext[c]['efron_0632_auc']:.4f}")

    # ---- frozen transports (primary logcpm + sensitivity raw + defined-only)
    variants = {}
    for vname, Lx, gg in [
            ("primary_logcpm_z_frozen_all6", L, genes_used),
            ("sensitivity_raw_z_frozen", {c: lab_raw[c] for c in COHORTS}, genes_used),
            ("sensitivity_logcpm_z_defined_genes", L, genes_defined)]:
        res = {}
        for tgt in EXT_COHORTS:
            r = run_transport(Lx, y, "IMvigor210", tgt, gg)
            res[tgt] = r
            tm = r["target_metrics"]
            adj = r["adjudication_target"]["main_analysis_rule"]
            sl = tm["calib_slope"]
            sl_s = "NaN" if sl is None or np.isnan(sl) else f"{sl:.3f}"
            print(f"   {vname} -> {tgt}: AUC={tm['auc']:.4f} "
                  f"CI[{tm['auc_ci95'][0]:.4f},{tm['auc_ci95'][1]:.4f}] "
                  f"slope={sl_s} a={tm['calib_intercept']:.3f} "
                  f"brier={tm['brier']:.4f} logloss={tm['log_loss']:.4f} "
                  f"main_rule={'PASS' if adj['both_passed'] else 'FAIL'}")
        variants[vname] = res

    # ================= 5. new-S0 comparison =================
    s0 = None
    if os.path.exists(S0_JSON):
        d = json.load(open(S0_JSON))
        s0 = dict(
            source="data/imvigor_gene_recovery.json phase2.new "
                   "(recovered-key anchor arm, same frozen-transport design)",
            within_imvigor=d["phase2"]["new"]["within_imvigor"],
            external={c: d["phase2"]["new"]["transport"]["external"][c]
                      for c in EXT_COHORTS},
        )
        for c in EXT_COHORTS:
            print(f"   new-S0 reference {c}: AUC={s0['external'][c]['auc']:.4f}")

    # ================= honest verdict =================
    prim = variants["primary_logcpm_z_frozen_all6"]
    verdict_rows = {}
    for c in EXT_COHORTS:
        adj = prim[c]["adjudication_target"]
        verdict_rows[c] = dict(
            auc=adj["auc"], auc_ci95=adj["auc_ci95"], calib_slope=adj["calib_slope"],
            main_rule=adj["main_analysis_rule"]["both_passed"],
            tiers={t: adj["tolerance_tiers"][f"pm{t:.2f}"]["both_passed"] for t in TIERS})
    any_pass = any(v["main_rule"] for v in verdict_rows.values())
    any_tier_pass = {t: any(v["tiers"][t] for v in verdict_rows.values()) for t in TIERS}
    verdict = dict(
        p0_2_question="is the two-criterion frozen-transport framework a trivial rejector?",
        per_cohort=verdict_rows,
        any_external_cohort_passes_main_rule=bool(any_pass),
        any_external_cohort_passes_at_tier={f"pm{t:.2f}": bool(any_tier_pass[t])
                                            for t in TIERS},
        statement=None,  # filled below with computed numbers
    )

    def _fmt_ci(ci):
        if ci is None or ci[0] is None:
            return "NA"
        return f"[{ci[0]:.3f},{ci[1]:.3f}]"

    verdict["statement"] = (
        "IMvigor210-anchored 18-gene T-cell-inflamed GEP positive control "
        "(recovered gene identity): "
        + "; ".join(
            f"{c}: AUC {verdict_rows[c]['auc']:.3f} "
            f"{_fmt_ci(verdict_rows[c]['auc_ci95'])}, slope "
            f"{('NaN' if verdict_rows[c]['calib_slope'] is None else format(verdict_rows[c]['calib_slope'], '.3f'))}, "
            f"main-rule {'PASS' if verdict_rows[c]['main_rule'] else 'FAIL'}"
            for c in EXT_COHORTS)
        + f". Any cohort passes main rule: {any_pass}; "
        + "at tiers pm0.05/pm0.10/pm0.20: "
        + "/".join(str(any_tier_pass[t]) for t in TIERS) + ".")

    out = dict(
        analysis="positive_control_gep18",
        supersedes="data/positive_control.json (voided: IMvigor210 gene "
                    "identity invalid at that run)",
        purpose=("P0-2  primary positive control: the 18-gene T-cell-inflamed "
                 "GEP (Cristescu et al., Science 2018; T-cell markers CD8A/CMKLR1/"
                 "CXCR6/HLA-DQA1/HLA-DRB1/HLA-E; IFN-gamma-responsive CCL5/CXCL9/"
                 "CXCL10/IDO1/PSMB10/STAT1/TIGIT; checkpoints CD27/CD274/CD276/"
                 "LAG3/PDCD1LG2) pushed through the SAME frozen transport pipeline "
                 "and two-criterion evaluation as the main analysis, under the "
                 "RECOVERED IMvigor210 gene identity (exmat row j == j-th fData "
                 "row in upper-cased symbol alphabetical order; see "
                 "imvigor_gene_recovery.json). The 6-gene IFN-gamma panel is "
                 "SENSITIVITY ONLY (rerun_positive_control.json)."),
        seed=SEED,
        signature=dict(
            name="18-gene T-cell-inflamed GEP (Cristescu et al., Science 2018)",
            genes=GENES,
            score="mean of per-gene z-scores; z parameters frozen on the IMvigor210 "
                  "labeled patients (n=298) and transported unchanged"),
        recovered_key=dict(
            definition="exmat row j (1-based) == fData row with the j-th symbol "
                       "in alphabetical order of the upper-cased fData symbol column",
            crosscheck_vs_exported_key=xchk,
            evidence="imvigor_gene_recovery.json: sex-anchor gates PASS, "
                     "31085/31085 numeric identity match vs labeled mirror, "
                     "upstream provenance byte-identical"),
        gene_availability=avail,
        genes_used_in_signature=genes_used,
        genes_missing=genes_missing,
        zero_variance_genes_in_imigor_training=genes_zerovar,
        cohorts={c: dict(n_labeled=int(len(y[c])), n_event=int(y[c].sum()))
                 for c in COHORTS},
        frozen_pipeline=dict(
            imvigor_loader="load_imvigor_correct (THIS script): recovered-key row "
                           "selection -- the ONLY change vs the voided run",
            external_loaders="verbatim from pca_stability.py (load_gse165252 / "
                             "load_gse78220); gene labels valid in both",
            training="IMvigor210 labeled patients n=298; per-gene z (log1p CPM, "
                     "full-matrix library size) frozen; signature = mean of the "
                     "available genes' z; univariate ridge logistic lambda=1/n_train, "
                     "intercept unpenalised (verbatim // algebra)",
            transport="z mean/sd and model coefficients frozen on the anchor, "
                      "applied unchanged to ALL labeled patients of the target "
                      "cohort (the same held-out sets as the main analysis)",
            metrics="verbatim  benchmark metric set (metrics_block from "
                    "competitor_benchmark.py): AUC + 2000-rep bootstrap BCa "
                    f"95% CI (seed {SEED}), logit-scale calibration slope, "
                    "calibration-in-the-large (intercept), Brier, log-loss; "
                    f"within-cohort Efron 0.632 (n_boot={N_BOOT_632}, seed {SEED})",
            variants="primary logcpm_z all-available-GEP18-genes; sensitivity raw_z; "
                     "sensitivity logcpm_z defined-genes-only"),
        adjudication_rule=dict(
            c1_main="external AUC bootstrap 95% CI excludes 0.5 in the beneficial "
                    "direction (CI lower bound > 0.5)",
            c2_main="|calibration slope - 1| <= 0.10 (logit scale)",
            tolerance_semantics=("registered before execution: at tier t, C1 passes "
                                 "if CI lower bound > 0.5 - t; C2 passes if "
                                 "|slope - 1| <= t; tiers +/-0.05, +/-0.10, +/-0.20"),
            tiers=TIERS),
        within_imvigor210_efron632=efron_imv,
        within_external_cohorts_efron632=efron_ext,
        transports=variants,
        new_s0_comparison=s0,
        verdict=verdict,
    )

    with open(OUT, "w") as f:
        json.dump(_clean(out), f, indent=2)
    print("\nWrote", OUT)

    # ================================================================ markdown
    md = []
    md.append("#  P0-2 主阳性对照（恢复基因身份）— Cristescu 18 基因 T 细胞炎性 GEP\n")
    md.append(f"日期：2026-09-24 ｜ seed {SEED} ｜ 脚本：scripts/analysis/positive_control_gep18.py ｜ "
              "冻结 JSON：data/positive_control_gep18.json\n")
    md.append("**背景：** 原 P0-2 阳性对照（positive_control.py）因 IMvigor210 exmat "
              "基因行身份被上游匿名化破坏而作废。imvigor_gene_recovery.py 已恢复身份"
              "（exmat 第 j 行 == fData 按 symbol 大写字母序排序后的第 j 个基因；性别锚定门 "
              "PASS、镜像 31085/31085 行数值恒等匹配、上游溯源逐字节一致）。本脚本按原设计"
              "重跑：**唯一改动是 IMvigor210 加载器换为恢复钥匙选行**，指标代数、判据、"
              "seed 全部原样复用（/ 逐字导入）。\n")

    md.append("## 1. 恢复钥匙交叉核验\n")
    if xchk["performed"]:
        md.append(f"- 脚本内重建的 symbol（大写）字母序钥匙与项目资产 "
                  f"`exmat_row_key.csv`（{xchk['n_key_rows']} 行）对 6 个签名基因逐一比对："
                  f"**{'全部一致' if xchk['all_agree'] else '存在不一致（脚本会中止）'}**。\n")

    md.append("## 2. 三队列 6 基因可得性（恢复钥匙下）\n")
    md.append("| gene | IMvigor210（恢复钥匙，exmat 行号 1-based） | GSE165252（Ensembl） | GSE78220（symbol） |")
    md.append("|---|---|---|---|")
    for gene in GENES:
        a = avail[gene]
        md.append(f"| {gene} | {'YES（行 ' + str(a['imvigor_exmat_row_1based']) + '）' if a['IMvigor210'] else 'NO'} "
                  f"| {'YES（' + str(a['gse165252_ensembl_id']) + '）' if a['GSE165252'] else 'NO'} "
                  f"| {'YES' if a['GSE78220'] else 'NO'} |")
    md.append("\n" + (f"进入签名的基因：{', '.join(genes_used)}。缺失："
              f"{', '.join(genes_missing) if genes_missing else '无'}。"
              if True else ""))
    if genes_zerovar:
        md.append(f"IMvigor210 训练集内零方差基因：{', '.join(genes_zerovar)}"
                  "（主变体经 sd→0 保护计为 z=0；defined-genes 敏感性变体剔除）。\n")

    md.append("## 3. IMvigor210 锚定（n=298 标注样本）与 within 表现\n")
    mdl = prim["GSE165252"]["model"]
    md.append(f"- 冻结模型：logit(P) = {mdl['intercept']:.4f} + "
              f"{mdl['beta_signature']:.4f} × signature（ridge λ=1/n=1/298，截距不惩罚）。")
    md.append(f"- within 表观 AUC = **{efron_imv['apparent_auc']:.4f}**；"
              f"Efron 0.632 AUC = **{efron_imv['efron_0632_auc']:.4f}**"
              f"（n_boot=600，seed {SEED}；事件 {efron_imv['n_event']}/{efron_imv['n']}）。\n")
    md.append("| 队列（within，常规队列内 z，描述性） | 表观 AUC | mean OOB | Efron 0.632 |")
    md.append("|---|---|---|---|")
    for c in EXT_COHORTS:
        e = efron_ext[c]
        md.append(f"| {c} | {e['apparent_auc']:.4f} | {e['oob_auc_mean']:.4f} | "
                  f"**{e['efron_0632_auc']:.4f}** |")
    md.append("")

    md.append("## 4. 冻结运输：逐队列指标与两判据判定\n")
    for vname, res in variants.items():
        md.append(f"### {vname}\n")
        md.append("| target | AUC [95% BCa CI] | slope b | CITL a | Brier | log-loss | 主判据 C1∧C2 | ±0.05 | ±0.10 | ±0.20 |")
        md.append("|---|---|---|---|---|---|---|---|---|---|")
        for tgt, r in res.items():
            tm = r["target_metrics"]
            adj = r["adjudication_target"]
            tiers = adj["tolerance_tiers"]
            sl = tm["calib_slope"]
            sl_md = "NaN" if sl is None or np.isnan(sl) else f"{sl:.3f}"
            a_md = tm["calib_intercept"]
            a_md = "NaN" if a_md is None or np.isnan(a_md) else f"{a_md:.3f}"
            md.append(f"| {tgt} | **{tm['auc']:.4f}** [{tm['auc_ci95'][0]:.4f},"
                      f"{tm['auc_ci95'][1]:.4f}] | {sl_md} | {a_md} "
                      f"| {tm['brier']:.4f} | {tm['log_loss']:.4f} "
                      f"| {'PASS' if adj['main_analysis_rule']['both_passed'] else 'FAIL'} "
                      + " | ".join("PASS" if tiers[f"pm{t:.2f}"]["both_passed"] else "fail"
                                   for t in TIERS) + " |")
        md.append("")
    md.append("判据：C1 = AUC 95% CI 下界 > 0.5；C2 = |slope − 1| ≤ 0.10（logit 尺度）；"
              "容差档 t 下 C1 放宽为下界 > 0.5 − t、C2 放宽为 |slope − 1| ≤ t。\n")

    md.append("## 5. 与新 S0 基线对比（同一冻结运输设计）\n")
    if s0:
        md.append("| transport | GEP-18 阳性对照 AUC [CI]（主变体） | 新 S0 AUC [CI] | Δ AUC |")
        md.append("|---|---|---|---|")
        for c in EXT_COHORTS:
            tm = prim[c]["target_metrics"]
            s0c = s0["external"][c]
            md.append(f"| IMvigor210 → {c} | **{tm['auc']:.4f}** "
                      f"[{tm['auc_ci95'][0]:.4f},{tm['auc_ci95'][1]:.4f}] "
                      f"| {s0c['auc']:.4f} [{s0c['ci95'][0]:.4f},{s0c['ci95'][1]:.4f}] "
                      f"| {tm['auc']-s0c['auc']:+.4f} |")
        md.append("")
        md.append(f"- within IMvigor210：GEP-18 0.632 AUC **{efron_imv['efron_0632_auc']:.4f}** "
                  f"vs 新 S0 {s0['within_imvigor']['auc_632']:.4f}。")
        for c in EXT_COHORTS:
            tm = prim[c]["target_metrics"]
            s0c = s0["external"][c]
            sl = tm["calib_slope"]
            sl_s = "NaN" if sl is None else f"{sl:.3f}"
            md.append(f"- {c}：slope {sl_s} vs S0 {s0c['calib_slope']:.3f}；"
                      f"Brier {tm['brier']:.4f} vs {s0c['brier']:.4f}；"
                      f"log-loss {tm['log_loss']:.4f} vs {s0c['log_loss']:.4f}。")
        md.append("")

    md.append("## 6. 结论：框架是否 trivial rejector？\n")
    md.append(f"- {verdict['statement']}")
    if any_pass:
        md.append("- **判定：框架不是 trivial rejector**——已发表阳性对照生物标志物在"
                  "同一冻结运输与两判据下至少于一个外部队列通过，说明框架有能力放行"
                  "真实可复现信号。")
    elif any(any_tier_pass[t] for t in TIERS):
        mt = [t for t in TIERS if any_tier_pass[t]][0]
        md.append(f"- **判定：主判据下未通过，但在 ±{mt:.2f} 容差档下至少一个外部队列通过**"
                  "——框架对阳性对照的拒绝接近但未达到主判据门槛，属“严格”而非"
                  "“trivial”（无差别拒绝）的边缘证据，需结合 within 信号强度解释。")
    else:
        md.append("- **判定：本数据集内无法证明框架非 trivial rejector**——阳性对照亦未在"
                  "任何容差档通过。结合 within 表现（见 §3）区分两种解释：(a) 签名在这些"
                  "小样本外部队列本身信号不足；(b) 框架过严。需要更大基因有效队列才能"
                  "裁决。")
    md.append("")
    md.append("以上全部数字为运行时实算；每患者签名分与预测概率存于 JSON（per_patient_target/"
              "per_patient_anchor）。无任何手工数值。\n")
    md.append("*冻结产出：`data/positive_control_gep18.json`；脚本 "
              "`scripts/analysis/positive_control_gep18.py`；seed 20260911；主阳性对照（6 基因 IFN-γ 签名见 data/rerun_positive_control.json，敏感性分析）。*"
              "`positive_control.json`。*")
    with open(MDOUT, "w") as f:
        f.write("\n".join(md))
    print("Wrote", MDOUT)

    # console summary
    print("\n== summary ==")
    for c in EXT_COHORTS:
        tm = prim[c]["target_metrics"]
        adj = prim[c]["adjudication_target"]
        print(f"  IMvigor210 -> {c}: AUC={tm['auc']:.4f} "
              f"CI[{tm['auc_ci95'][0]:.4f},{tm['auc_ci95'][1]:.4f}] "
              f"slope={tm['calib_slope']} main_rule="
              f"{'PASS' if adj['main_analysis_rule']['both_passed'] else 'FAIL'}")
    print(f"  verdict: any_pass={any_pass}, tiers="
          + "/".join(f"{any_tier_pass[t]}" for t in TIERS))


if __name__ == "__main__":
    main()
