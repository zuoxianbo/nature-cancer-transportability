#!/usr/bin/env python3
""" -- Nature Cancer revision , P0-4 submission gate.

Establish a TRUE clinical-only competitive baseline (not just prevalence-only)
for the frozen-S0 transport analysis, on the IDENTICAL held-out patients as
.

Design (mirrors the main analysis transport direction):
    train  IMvigor210 labeled n=298  (ridge logistic, lambda=1/n_train,
                                      intercept unpenalised, / verbatim
                                      algebra, one-hot categoricals, honest
                                      missingness handling + recording)
    frozen -> test  GSE165252 (n=71) and GSE78220 (n=25), the SAME labeled
    patients as the frozen S0 /  benchmark (rebuilt label vectors are
    asserted equal to the frozen _eff_*_y.npy at run time).

Honest variable availability audit (from the actual source files):
  - IMvigor210 pData: Sex, Baseline ECOG, FMOne TMB (64/298 missing),
    Met.Disease.Status (27/298 missing), TCGA.Subtype (stage group I-IV),
    Tissue, Received.platinum, Intravesical.BCG, Tobacco.  NOTE: the column
    "Sample.age" is a SAMPLE-TIMING category ("(less than) 1 year", "1-2
    years", "more than 2 years"), NOT patient age; IMvigor210 carries NO
    patient age variable in the deposited pData.
  - GSE78220 series matrix: gender, age (yrs), disease status (AJCC M
    category M0/M1a/M1b/M1c), previous MAPKi, vital status, per-gene
    mutations (NOT TMB).
  - GSE165252 series matrix: ONLY tissue / treatment / response.  No
    patient-level covariate of any kind is deposited in GEO.
Consequence: the common covariate set across ALL THREE cohorts is EMPTY.
Maximal harmonizable set (IMvigor210 -> GSE78220 transportable):
  {sex, advanced/metastatic-stage binary}.  Patient age cannot enter any
  trained model (absent in the training cohort).  On GSE165252 the
  clinical-only baseline necessarily degenerates to a constant (covariate-
  marginal) prediction -- recorded honestly, not skipped.

Models:
  C_main : ridge logistic on {sex_M, stage_adv}             (harmonizable set)
  C_sex  : ridge logistic on {sex_M} only                   (minimal fallback)
  C_rich : IMvigor210-only sensitivity model with the full local covariate
           set (sex, stage, ECOG, log10 TMB imputed + missingness indicator,
           met-status one-hot + missing indicator, platinum, BCG, tobacco),
           evaluated WITHIN IMvigor210 by apparent + Efron 0.632 (600
           bootstrap refits, seed 20260911) -- never transported.

Metrics: IDENTICAL to  (imported verbatim from _competitor_benchmark):
log-loss, Brier, calibration-in-the-large (a), logit-scale calibration slope
(b), AUC + 2000-rep bootstrap BCa 95% CI (seed 20260911).  External AUC 0.632
is not defined (no target refits); within-cohort 0.632 uses the 
efron632_within formula on the clinical design matrices.

Outputs: data/clinical_only.json + clinical_only.md
Seed 20260911 everywhere.  Nothing hardcoded; all numbers computed at run
time.  Any degeneracy is reported as such, never fabricated.
"""
import json, os, sys, csv
import numpy as np

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
BASE = f"{ROOT}/outputs/analysis_base"
SCRIPTS = f"{BASE}/scripts"
sys.path.insert(0, SCRIPTS)
WS = f"{ROOT}/data/_workspace"

# ---- / VERBATIM algebra and constants (import guarantees identity) ----
from competitor_benchmark import (
    sigmoid, fit_logistic, metrics_block, logit, SEED, N_BOOT)

OUTDIR = f"{BASE}/data/"
OUT = f"{OUTDIR}/clinical_only.json"
MDOUT = f"{OUTDIR}/clinical_only.md"
V75_JSON = f"{BASE}/data/competitor_benchmark.json"
N_BOOT_632 = 600
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
GSE165_ANN = f"{ROOT}/geo_data/GSE165252/series_matrix_annotation.txt"
GSE165_MAT = f"{ROOT}/geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz"
GSE782_ANN = f"{ROOT}/geo_data/GSE78220/series_matrix_annotation.txt"
GSE782_MAT = f"{ROOT}/geo_data/GSE78220/GSE78220_PatientFPKM.xlsx"
COHORTS = ["IMvigor210", "GSE165252", "GSE78220"]

# ================================================================== loaders
def load_imvigor():
    """Rebuild the IMvigor210 labeled subset (n=298) EXACTLY as the frozen
    pipeline: sample order = exmat.csv header; labels from pData binaryResponse
    (CR/PR=1, SD/PD=0); covariates joined by sample id.  Asserted equal to
    _eff_IMvigor210_y.npy."""
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f); h = next(rd)
        sids = [x.strip('"') for x in h[1:]]
    pd_rows = list(csv.DictReader(open(IMV_PDATA)))
    pmap = {r["X"]: r for r in pd_rows}
    y = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                  else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                  else np.nan for s in sids], dtype=float)
    yf = np.load(f"{WS}/_eff_IMvigor210_y.npy")
    assert np.array_equal(yf, y, equal_nan=True), "IMvigor210 y mismatch vs frozen"
    m = ~np.isnan(y)
    sids_lab = [s for s, k in zip(sids, m) if k]
    cov = []
    for s in sids_lab:
        r = pmap[s]
        tmb = r["FMOne.mutation.burden.per.MB"]
        cov.append(dict(
            sid=s,
            sex_male=1.0 if r["Sex"] == "M" else 0.0,
            stage4=1.0 if r["TCGA.Subtype"] == "IV" else 0.0,
            stage_group=r["TCGA.Subtype"],
            ecog=r["Baseline.ECOG.Score"],
            tmb=float(tmb) if tmb not in ("", "NA") else np.nan,
            met_status=r["Met.Disease.Status"],           # may be "NA"
            platinum=1.0 if r["Received.platinum"] == "Y" else 0.0,
            bcg=1.0 if r["Intravesical.BCG.administered"] == "Y" else 0.0,
            tobacco=r["Tobacco.Use.History"],
        ))
    return cov, y[m]

def load_gse78220():
    """Rebuild the GSE78220 labeled subset (n=25) EXACTLY as : xlsx columns
    split at '.baseline'; labels from the annotation 'anti-pd-1 response'
    mapped by patient id; clinical covariates (gender / age / AJCC M category)
    parsed per-position from the characteristics lines (which are NOT
    uniformly keyed -- each position carries its own 'key: value' tag).
    Pt16.OnTx / Pt27A / Pt27B columns fall outside the patient-id crosswalk
    exactly as in the frozen pipeline; Pt27 therefore contributes no labeled
    sample.  Asserted equal to _eff_GSE78220_y.npy (labeled subset)."""
    import openpyxl
    from collections import defaultdict
    ann = open(GSE782_ANN).read().splitlines()
    lines = [[v.strip().strip('"') for v in ln.rstrip('\n').split('\t')[1:]]
             for ln in ann if ln.startswith("!Sample_characteristics_ch1")]
    npos = len(lines[0])
    per_pos = [dict() for _ in range(npos)]
    for vals in lines:
        for i, v in enumerate(vals):
            if ":" in v:
                k, val = v.split(":", 1)
                per_pos[i][k.strip()] = val.strip()
    pid_cov = {}
    for pos in per_pos:
        pid = pos.get("patient id")
        if not pid:
            continue
        d = pid_cov.setdefault(pid, {})
        for ann_key, out_key in [("gender", "gender"), ("age (yrs)", "age"),
                                 ("disease status", "disease"),
                                 ("anti-pd-1 response", "response")]:
            if pos.get(ann_key) is not None:
                d[out_key] = pos[ann_key]
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    hdr = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
    raw_cols = [str(h) for h in hdr[1:]]
    sample_cols = [h.split(".baseline")[0].split(".Baseline")[0] for h in raw_cols]
    y = []
    for sc in sample_cols:
        r = pid_cov.get(sc, {}).get("response")
        rl = r.lower() if r else ""
        if "progressive" in rl or "stable" in rl:
            y.append(0.0)
        elif "complete" in rl or "partial" in rl:
            y.append(1.0)
        else:
            y.append(np.nan)
    y = np.array(y, dtype=float)
    yf = np.load(f"{WS}/_eff_GSE78220_y.npy")
    m = ~np.isnan(y)
    assert np.array_equal(yf, y[m]), "GSE78220 y mismatch vs frozen"
    cov = []
    for sc, k in zip(sample_cols, m):
        if not k:
            continue
        d = pid_cov[sc]
        cov.append(dict(
            pid=sc,
            sex_male=1.0 if d.get("gender") == "M" else 0.0,
            age=float(d["age"]) if d.get("age") not in (None, "NA") else np.nan,
            m_stage=d.get("disease"),
            stage_adv=1.0 if str(d.get("disease", "")).startswith("M1") else 0.0,
        ))
    excluded = [rc for rc, k in zip(raw_cols, m) if not k]
    return cov, y[m], excluded

def load_gse165252():
    """Rebuild the GSE165252 labeled subset (n=71) EXACTLY as : norm.cnt
    columns; response from the characteristics 'response:' line mapped through
    the !Sample_description TE codes.  AUDIT RESULT (recorded, not assumed):
    the deposited GEO metadata carries ONLY tissue / treatment / response --
    no patient-level covariate of any kind.  Asserted equal to the frozen y."""
    import gzip
    with gzip.open(GSE165_MAT, 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        codes = [h.strip().strip('"') for h in hdr[1:]]
    ann = open(GSE165_ANN).read().splitlines()
    sample_lines = {}
    char_lines = []
    for ln in ann:
        if ln.startswith("!Sample"):
            fld = ln.rstrip('\n').split('\t')
            vals = [v.strip().strip('"') for v in fld[1:]]
            if fld[0].startswith("!Sample_characteristics"):
                char_lines.append(vals)
            sample_lines.setdefault(fld[0], []).append(vals)
    code_line = None
    for name, vlist in sample_lines.items():
        for vals in vlist:
            if name.startswith("!Sample_description") and \
               len(set(vals) & set(codes)) >= max(1, int(0.9 * len(codes))):
                code_line = vals
                break
        if code_line:
            break
    resp_line = None
    for vals in char_lines:
        if any(v.lower().startswith("response:") for v in vals):
            resp_line = vals
            break
    code2resp = {}
    for c, r in zip(code_line, resp_line):
        rl = r.lower()
        if "non-responder" in rl:
            code2resp[c] = 0.0
        elif "responder" in rl:
            code2resp[c] = 1.0
        else:
            code2resp[c] = np.nan
    y = np.array([code2resp.get(c, np.nan) for c in codes], dtype=float)
    yf = np.load(f"{WS}/_eff_GSE165252_y.npy")
    m = ~np.isnan(y)
    assert np.array_equal(yf, y[m]), "GSE165252 y mismatch vs frozen"
    keys = set()
    for vals in char_lines:
        for v in vals:
            if ":" in v:
                keys.add(v.split(":", 1)[0].strip().lower())
    return [c for c, k in zip(codes, m) if k], y[m], sorted(keys)

# ================================================================== helpers
def efron632_clinical(X, y, lam, n_boot=N_BOOT_632, seed=SEED):
    """ efron632_within applied to a CLINICAL design matrix: apparent AUC +
    Efron 0.632 optimism-corrected AUC (0.368*apparent + 0.632*mean OOB AUC
    over bootstrap refits, lambda fixed at the full-cohort rule)."""
    from sklearn.metrics import roc_auc_score
    n = len(y)
    def fit_pred(Xa, ya, Xe):
        w = fit_logistic(np.column_stack([np.ones(len(ya)), Xa]), ya, lam)
        return sigmoid(np.column_stack([np.ones(len(Xe)), Xe]) @ w)
    app = float(roc_auc_score(y, fit_pred(X, y, X)))
    rng = np.random.default_rng(seed)
    oob, skipped = [], 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        oob_i = np.setdiff1d(np.arange(n), idx)
        if len(np.unique(y[idx])) < 2 or len(np.unique(y[oob_i])) < 2:
            skipped += 1
            continue
        p_oob = fit_pred(X[idx], y[idx], X[oob_i])
        oob.append(float(roc_auc_score(y[oob_i], p_oob)))
    oob_mean = float(np.mean(oob)) if oob else float("nan")
    return dict(
        n=int(n), n_event=int(y.sum()),
        lambda_rule=f"1/n (n={n}, intercept unpenalised;  fit_logistic verbatim)",
        apparent_auc=app, oob_auc_mean=oob_mean,
        n_boot=n_boot, n_boot_valid=len(oob), n_boot_skipped=skipped,
        efron_0632_auc=float(0.368 * app + 0.632 * oob_mean),
        formula=("AUC_0.632 = 0.368 * AUC_apparent + 0.632 * mean(AUC_out-of-bag "
                 f"over {n_boot} bootstrap refits)"),
        seed=seed)

def smd_prop(p1, p2):
    denom = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2.0)
    return float((p1 - p2) / denom) if denom > 0 else float("nan")

def adjudicate(mres):
    lo, hi = mres["auc_ci95"]
    auc_ok = bool(lo is not None and not np.isnan(lo) and lo > 0.5)
    sl = mres["calib_slope"]
    slope_ok = bool(sl is not None and not np.isnan(sl) and abs(sl - 1.0) <= 0.1)
    return dict(auc_ci_excludes_0p5_beneficial_direction=auc_ok,
                slope_within_tolerance=slope_ok,
                transport_permitted=bool(auc_ok and slope_ok))

# ================================================================== main
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    cov_imv, y_imv = load_imvigor()
    cov_782, y_782, excluded_782 = load_gse78220()
    ids_165, y_165, keys_165 = load_gse165252()
    n_imv, n_782, n_165 = len(y_imv), len(y_782), len(y_165)
    assert (n_imv, int(y_imv.sum())) == (298, 68)
    assert (n_165, int(y_165.sum())) == (71, 23)
    assert (n_782, int(y_782.sum())) == (25, 13)

    # ---------------- variable availability matrix (from the actual files) ----
    tmb_missing = sum(1 for c in cov_imv if np.isnan(c["tmb"]))
    met_missing = sum(1 for c in cov_imv if c["met_status"] == "NA")
    availability = {
        "sex": {
            "IMvigor210": {"status": "available", "detail": "pData Sex (M/F), 0/298 missing"},
            "GSE165252": {"status": "unavailable",
                          "detail": "GEO metadata carries only tissue/treatment/response"},
            "GSE78220": {"status": "available", "detail": "characteristics 'gender' (M/F), 25/25"},
            "harmonizable": "IMvigor210 + GSE78220 (same construct, same coding)"},
        "patient_age": {
            "IMvigor210": {"status": "unavailable",
                           "detail": ("pData has NO patient age; the column 'Sample.age' is a "
                                      "sample-timing category ('(less than) 1 year' / '1-2 years' / "
                                      "'more than 2 years'), NOT age")},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "available", "detail": "characteristics 'age (yrs)', numeric"},
            "harmonizable": ("NO -- age exists in only one cohort and NOT in the training "
                             "cohort, so it cannot enter any trained clinical model")},
        "stage": {
            "IMvigor210": {"status": "available",
                           "detail": "TCGA.Subtype stage group I/II/III/IV (0/298 missing)"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "available",
                         "detail": "characteristics 'disease status' = AJCC M category M0/M1a/M1b/M1c"},
            "harmonizable": ("crude binary only: IMvigor210 stage group IV vs I-III; GSE78220 "
                             "M1(any) vs M0. Constructs differ (TNM stage group vs M category) "
                             "-- flagged as a pre-declared crude harmonization")},
        "ecog_performance_status": {
            "IMvigor210": {"status": "available", "detail": "Baseline.ECOG.Score 0/1/2, 0/298 missing"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "unavailable", "detail": "not deposited"},
            "harmonizable": "NO (single cohort)"},
        "tumor_mutation_burden": {
            "IMvigor210": {"status": "available",
                           "detail": f"FMOne.mutation.burden.per.MB, {tmb_missing}/298 missing"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "unavailable",
                         "detail": "only per-gene mutation calls (BRAF/NRAS/NF1), no TMB"},
            "harmonizable": "NO -- IMvigor210-only; used in the IMvigor210-only sensitivity model"},
        "metastatic_disease_status": {
            "IMvigor210": {"status": "available",
                           "detail": f"Met.Disease.Status (Visceral/Liver/LN Only), {met_missing}/298 missing"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "partially_available",
                         "detail": "folded into the AJCC M category, not the same construct"},
            "harmonizable": "NO -- IMvigor210-only sensitivity"},
        "prior_therapy": {
            "IMvigor210": {"status": "available",
                           "detail": "Received.platinum (Y/N), Intravesical.BCG.administered (Y/N)"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "available", "detail": "previous MAPKi (Y/N)"},
            "harmonizable": ("NO -- different drugs/classes (platinum+BCG in bladder cancer vs "
                             "MAPK inhibitors in melanoma); not a common construct")},
        "smoking": {
            "IMvigor210": {"status": "available", "detail": "Tobacco.Use.History NEVER/PREVIOUS/CURRENT"},
            "GSE165252": {"status": "unavailable", "detail": "not deposited"},
            "GSE78220": {"status": "unavailable", "detail": "not deposited"},
            "harmonizable": "NO (single cohort)"},
        "response_outcome": {
            "IMvigor210": {"status": "available", "detail": "binaryResponse CR/PR=1, SD/PD=0 (n=298)"},
            "GSE165252": {"status": "available", "detail": "responder/non-responder (n=71)"},
            "GSE78220": {"status": "available", "detail": "anti-PD-1 response CR/PR vs SD/PD (n=25)"},
            "harmonizable": "label, not a predictor"},
    }
    common_all3 = []   # honestly EMPTY
    harmonizable_set = ["sex", "advanced_stage_binary"]

    # ---------------- design matrices ----------------
    X_imv_main = np.array([[c["sex_male"], c["stage4"]] for c in cov_imv])
    X_imv_sex = np.array([[c["sex_male"]] for c in cov_imv])
    # IMvigor210-only rich model: honest missingness handling (imputation +
    # missingness indicators), every rule recorded.
    tmb_vals = np.array([c["tmb"] for c in cov_imv])
    tmb_med = float(np.nanmedian(tmb_vals))
    met_cats = ["Visceral", "Liver", "LN Only"]
    tob_cats = ["PREVIOUS", "CURRENT"]
    rows_rich = []
    for c in cov_imv:
        tmb_missing_i = 1.0 if np.isnan(c["tmb"]) else 0.0
        logtmb = np.log10(tmb_med + 1.0) if np.isnan(c["tmb"]) else np.log10(c["tmb"] + 1.0)
        met_missing_i = 1.0 if c["met_status"] == "NA" else 0.0
        row = [c["sex_male"], c["stage4"],
               1.0 if c["ecog"] == "1" else 0.0, 1.0 if c["ecog"] == "2" else 0.0,
               logtmb, tmb_missing_i] \
              + [1.0 if c["met_status"] == k else 0.0 for k in met_cats] + [met_missing_i] \
              + [c["platinum"], c["bcg"]] \
              + [1.0 if c["tobacco"] == k else 0.0 for k in tob_cats]
        rows_rich.append(row)
    X_imv_rich = np.array(rows_rich)
    rich_feat_names = ["sex_male", "stage4", "ecog1", "ecog2", "log10_tmb_plus1",
                       "tmb_missing"] + [f"met_{k.replace(' ', '_').lower()}" for k in met_cats] \
                      + ["met_missing", "platinum", "bcg", "tob_previous", "tob_current"]

    lam_imv = 1.0 / n_imv

    def fit_frozen(X, y):
        w = fit_logistic(np.column_stack([np.ones(len(y)), X]), y, lam_imv)
        return w

    w_main = fit_frozen(X_imv_main, y_imv)
    w_sex = fit_frozen(X_imv_sex, y_imv)

    # ---------------- external application ----------------
    X_782_main = np.array([[c["sex_male"], c["stage_adv"]] for c in cov_782])
    X_782_sex = np.array([[c["sex_male"]] for c in cov_782])
    p_782_main = sigmoid(w_main[0] + X_782_main @ w_main[1:])
    p_782_sex = sigmoid(w_sex[0] + X_782_sex @ w_sex[1:])
    # GSE165252: no covariates deposited -> covariate-marginal constant
    # prediction = mean training-set predicted probability (documented degeneracy)
    p_imv_main = sigmoid(w_main[0] + X_imv_main @ w_main[1:])
    p_imv_sex = sigmoid(w_sex[0] + X_imv_sex @ w_sex[1:])
    const_165_main = float(p_imv_main.mean())
    const_165_sex = float(p_imv_sex.mean())

    ext = {}
    ext["GSE78220"] = {
        "note": ("identical 25 held-out patients as the frozen S0 /  benchmark "
                 "(label vector asserted equal to _eff_GSE78220_y.npy)"),
        "C_main": metrics_block(y_782, p_782_main),
        "C_sex": metrics_block(y_782, p_782_sex),
    }
    ext["GSE78220"]["C_main"]["features_used"] = ["sex_male", "stage_adv"]
    ext["GSE78220"]["C_sex"]["features_used"] = ["sex_male"]
    ext["GSE165252"] = {
        "note": ("identical 71 held-out patients as the frozen S0 /  benchmark "
                 "(label vector asserted equal to _eff_GSE165252_y.npy)"),
        "degeneracy": ("NO patient-level covariate is deposited for GSE165252 in GEO "
                       "(verified: the only characteristics keys are "
                       f"{keys_165}); the frozen clinical model therefore CANNOT produce "
                       "patient-varying predictions here and degenerates to a constant "
                       "covariate-marginal prediction (mean training-set predicted "
                       "probability). AUC is exactly 0.5 by construction and the "
                       "calibration slope is not estimable -- this is a data constraint, "
                       "recorded honestly, not a modelling shortcut."),
        "C_main": metrics_block(y_165, np.full(n_165, const_165_main)),
        "C_sex": metrics_block(y_165, np.full(n_165, const_165_sex)),
    }
    ext["GSE165252"]["C_main"]["constant_prediction"] = const_165_main
    ext["GSE165252"]["C_sex"]["constant_prediction"] = const_165_sex
    ext["GSE165252"]["train_prevalence_reference"] = float(y_imv.mean())

    adjud = {c: {m: adjudicate(ext[c][m]) for m in ("C_main", "C_sex")}
             for c in ("GSE78220", "GSE165252")}

    # ---------------- within-IMvigor210 0.632 (clinical models) ----------------
    within = {
        "C_main": efron632_clinical(X_imv_main, y_imv, lam_imv),
        "C_sex": efron632_clinical(X_imv_sex, y_imv, lam_imv),
        "C_rich": efron632_clinical(X_imv_rich, y_imv, lam_imv),
    }
    # apparent calibration of the rich model in-cohort (descriptive)
    w_rich = fit_frozen(X_imv_rich, y_imv)
    p_rich = sigmoid(w_rich[0] + X_imv_rich @ w_rich[1:])

    # ---------------- comparison vs  (frozen S0 transports, same patients) ----
     = json.load(open(V75_JSON))
    pw = ["true_external_pairwise_transports"]["transports"]
    cmp_out = {}
    for tgt, cmet in (("GSE165252", "C_main"), ("GSE78220", "C_main")):
        key = f"train_IMvigor210__test_{tgt}"
        m2 = pw[key]["M2"]; m0 = pw[key]["M0"]
        cm = ext[tgt][cmet]
        cmp_out[tgt] = dict(
            comparison=("C_main clinical-only vs M2 frozen-S0 vs M0 prevalence-only, "
                        "all trained on IMvigor210 (n=298) and applied to the SAME "
                        f"held-out patients of {tgt}"),
            C_main=dict(auc=cm["auc"], auc_ci95=cm["auc_ci95"],
                        calib_slope=cm["calib_slope"], log_loss=cm["log_loss"],
                        brier=cm["brier"]),
            M2_frozen_s0=dict(auc=m2["auc"], auc_ci95=m2["auc_ci95"],
                              calib_slope=m2["calib_slope"], log_loss=m2["log_loss"],
                              brier=m2["brier"]),
            M0_prevalence=dict(auc=m0["auc"], log_loss=m0["log_loss"], brier=m0["brier"]),
            delta_auc_C_minus_M2=float(cm["auc"] - m2["auc"]),
            delta_logloss_C_minus_M2=float(cm["log_loss"] - m2["log_loss"]),
            delta_brier_C_minus_M2=float(cm["brier"] - m2["brier"]),
        )

    # ---------------- case-mix descriptive table (P0-9) ----------------
    sex_m_imv = float(np.mean([c["sex_male"] for c in cov_imv]))
    sex_m_782 = float(np.mean([c["sex_male"] for c in cov_782]))
    adv_imv = float(np.mean([c["stage4"] for c in cov_imv]))
    adv_782 = float(np.mean([c["stage_adv"] for c in cov_782]))
    age_782 = np.array([c["age"] for c in cov_782], float)
    stage_dist_imv = {k: int(sum(1 for c in cov_imv if c["stage_group"] == k))
                      for k in ["I", "II", "III", "IV"]}
    mstage_dist_782 = {}
    for c in cov_782:
        mstage_dist_782[c["m_stage"]] = mstage_dist_782.get(c["m_stage"], 0) + 1
    casemix = {
        "IMvigor210": dict(n=n_imv, n_event=int(y_imv.sum()),
                           response_rate=float(y_imv.mean()),
                           pct_male=sex_m_imv,
                           patient_age="not available (Sample.age is a sample-timing category)",
                           stage_group_distribution=stage_dist_imv,
                           pct_stageIV=adv_imv,
                           ecog_distribution={"0": int(sum(1 for c in cov_imv if c["ecog"] == "0")),
                                              "1": int(sum(1 for c in cov_imv if c["ecog"] == "1")),
                                              "2": int(sum(1 for c in cov_imv if c["ecog"] == "2"))},
                           tmb_missing=tmb_missing),
        "GSE165252": dict(n=n_165, n_event=int(y_165.sum()),
                          response_rate=float(y_165.mean()),
                          pct_male="not available",
                          patient_age="not available",
                          stage="not available",
                          deposited_characteristics_keys=keys_165),
        "GSE78220": dict(n=n_782, n_event=int(y_782.sum()),
                         response_rate=float(y_782.mean()),
                         pct_male=sex_m_782,
                         patient_age=dict(mean=float(age_782.mean()),
                                          sd=float(age_782.std(ddof=1)),
                                          min=float(age_782.min()),
                                          max=float(age_782.max())),
                         m_stage_distribution=mstage_dist_782,
                         pct_M1=adv_782),
        "standardized_differences_IMvigor210_vs_GSE78220": dict(
            pct_male_smd=smd_prop(sex_m_imv, sex_m_782),
            pct_advanced_stage_smd=smd_prop(adv_imv, adv_782),
            note=("positive SMD = IMvigor210 higher; |SMD|>0.2 conventionally indicates "
                  "material case-mix shift. Sex: %.1f%% vs %.1f%% male; advanced stage "
                  "(IMvigor210 stage-group IV vs GSE78220 AJCC M1): %.1f%% vs %.1f%%."
                  % (100 * sex_m_imv, 100 * sex_m_782, 100 * adv_imv, 100 * adv_782))),
        "prevalence_shift": dict(
            IMvigor210=float(y_imv.mean()), GSE165252=float(y_165.mean()),
            GSE78220=float(y_782.mean()),
            note="response-rate case-mix shift across the three cohorts"),
    }

    # ---------------- limitation statement draft ----------------
    limitation = (
        "Clinical-only baseline limitation (P0-4 audit, computed not assumed): the "
        "patient-level covariates deposited with the three efficacy cohorts are not "
        "mutually harmonizable. GSE165252 (PERFECT) deposits no patient-level covariate "
        "of any kind (only tissue, treatment arm and response), the IMvigor210 pData "
        "contains no patient age (its 'Sample.age' column is a sample-timing category), "
        "and the staging constructs differ (IMvigor210 TCGA stage group vs GSE78220 AJCC "
        "M category). The covariate set common to all three cohorts is therefore empty, "
        "and the maximal set harmonizable between the training cohort (IMvigor210) and "
        "any external cohort is {sex, a crude advanced/metastatic-stage binary}, "
        "transportable only to GSE78220. On GSE165252 every clinical-only model "
        "necessarily degenerates to a constant (covariate-marginal) prediction with AUC "
        "exactly 0.5, so for that cohort prevalence-only is the strongest clinical "
        "baseline the deposited data permit. These constraints are properties of the "
        "deposited data, not modelling choices; a richer clinical baseline would require "
        "patient-level covariate tables from the trial sponsors/publications that are not "
        "part of the GEO deposits."
    )

    out = dict(
        analysis="clinical_only",
        purpose=("P0-4 submission gate: true clinical-only competitive baseline on the "
                 "identical held-out patients as the frozen S0 /  benchmark"),
        seed=SEED,
        metadata=dict(n_boot_auc_ci=N_BOOT, n_boot_efron632=N_BOOT_632,
                      ridge_lambda_rule="lambda = 1/n_train (n_train=298), intercept unpenalised",
                      algebra="_competitor_benchmark.fit_logistic imported verbatim",
                      label_reconstruction_check=("rebuilt label vectors asserted equal to "
                                                  "_eff_{IMvigor210,GSE165252,GSE78220}_y.npy "
                                                  "(frozen S0 pipeline subsets: 298/71/25)")),
        variable_availability_matrix=availability,
        harmonization_rules=dict(
            sex="male=1, female=0 (identical construct and coding in IMvigor210 pData 'Sex' and GSE78220 'gender')",
            advanced_stage_binary=("IMvigor210: TCGA.Subtype=='IV' -> 1 else 0; GSE78220: AJCC M "
                                   "category startswith 'M1' -> 1, 'M0' -> 0. PRE-DECLARED CRUDE "
                                   "harmonization across different staging constructs; flagged"),
            patient_age="EXCLUDED from all trained models: not present in IMvigor210 (training cohort)",
            missingness=("C_main/C_sex: no missing values in their 1-2 features (verified). "
                         "C_rich (IMvigor210-only sensitivity): TMB missing in "
                         f"{tmb_missing}/298 -> median imputation of log10(TMB+1) + missingness "
                         f"indicator; Met.Disease.Status missing in {met_missing}/298 -> "
                         "'missing' indicator with all three category dummies set to 0"),
        ),
        common_set_all_three=common_all3,
        common_set_all_three_statement=("EMPTY: no patient-level covariate is available in all "
                                        "three cohorts; recorded honestly per task P0-4 step 4"),
        maximal_harmonizable_set=harmonizable_set,
        cohorts=dict(
            IMvigor210=dict(n=n_imv, n_event=int(y_imv.sum()),
                            excluded_columns=[]),
            GSE165252=dict(n=n_165, n_event=int(y_165.sum()),
                           sample_ids=ids_165),
            GSE78220=dict(n=n_782, n_event=int(y_782.sum()),
                          patient_ids=[c["pid"] for c in cov_782],
                          excluded_columns=excluded_782,
                          exclusion_note=("Pt16.OnTx / Pt27A / Pt27B columns fall outside the "
                                          "patient-id crosswalk exactly as in the frozen / "
                                          "pipeline; Pt27 therefore contributes no labeled sample")),
        ),
        models=dict(
            C_main=dict(
                definition="ridge logistic on {sex_male, stage_adv} trained on IMvigor210 n=298",
                features=["sex_male", "stage_adv"],
                intercept=float(w_main[0]),
                coefficients=dict(sex_male=float(w_main[1]), stage_adv=float(w_main[2])),
                frozen=True),
            C_sex=dict(
                definition="ridge logistic on {sex_male} only (minimal feasible fallback)",
                features=["sex_male"],
                intercept=float(w_sex[0]),
                coefficients=dict(sex_male=float(w_sex[1])),
                frozen=True),
            C_rich=dict(
                definition=("IMvigor210-ONLY sensitivity model (NOT transportable): full local "
                            "covariate set with honest missingness handling"),
                features=rich_feat_names,
                intercept=float(w_rich[0]),
                coefficients={k: float(v) for k, v in zip(rich_feat_names, w_rich[1:])},
                transportable=False),
        ),
        external_results=ext,
        adjudication_two_criterion=adjud,
        within_imvigor210_efron632=within,
        comparison_vs_=cmp_out,
        casemix=casemix,
        limitation_statement_draft=limitation,
    )
    json.dump(out, open(OUT, "w"), indent=2, allow_nan=True)
    print("Wrote", OUT)

    # ------------------------------------------------------------ console summary
    print("\n=== harmonizable set ===")
    print("common all-3:", common_all3, "| maximal (IMV->GSE78220):", harmonizable_set)
    print("\n=== external metrics ===")
    for c in ("GSE78220", "GSE165252"):
        for m in ("C_main", "C_sex"):
            r = ext[c][m]
            print(f"{c:10s} {m:7s} AUC={r['auc']:.4f} CI={r['auc_ci95']} "
                  f"slope={r['calib_slope']} logloss={r['log_loss']:.4f} brier={r['brier']:.4f}")
    print("\n=== within-IMvigor210 (Efron 0.632) ===")
    for m, r in within.items():
        print(f"{m:7s} apparent={r['apparent_auc']:.4f} oob={r['oob_auc_mean']:.4f} "
              f"0.632={r['efron_0632_auc']:.4f}")
    print("\n=== vs  (train IMvigor210, same patients) ===")
    for tgt, d in cmp_out.items():
        print(f"{tgt:10s} C_main AUC={d['C_main']['auc']:.4f}  M2 AUC={d['M2_frozen_s0']['auc']:.4f}  "
              f"dAUC(C-M2)={d['delta_auc_C_minus_M2']:+.4f}")

if __name__ == "__main__":
    main()
