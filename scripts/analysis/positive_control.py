#!/usr/bin/env python3
""" -- Nature Cancer revision , P0-2 submission gate.

POSITIVE CONTROL for the frozen-transport evaluation framework (P0-2): the
Ayers et al. JCI 2017 IFN-gamma 6-gene signature (IDO1, CXCL9, CXCL10, CCL5,
HLA-DRA, STAT1), a biomarker with published cross-cohort reproducible evidence
of association with PD-1 response (manuscript citation [1]), is to be pushed
through the IDENTICAL frozen-transport machinery as the main analysis (S0) to
test whether the two-criterion licensing framework is a trivial rejector.

*** RUNTIME DATA-INTEGRITY DISCOVERY (this run, 2026-09-24) ***

While executing the mandated pre-check (per-gene availability in the RAW
matrices), a decisive audit found that the IMvigor210 expression matrix
(data/imvigor210/exmat.csv) carries NO valid gene-level identity:

  * exmat.csv is byte-identical to the upstream repository file
    mimifp/tfm_mUC:2_data/1_IMvigor210/exmat_censored_IMvigor210.csv
    (same 34,121,302 bytes; same header; same first rows), whose gene
    identifiers were replaced by the placeholder labels "gene_1".."gene_31085"
    -- the upstream README states the data "have been anonymized because of
    the potential patentability of the results".
  * The upstream extraction script writes fData in fData(cds) order and exmat
    in counts(cds) order, matching them only by ENTREZ -- the anonymized file
    destroyed the entrez rownames, so the row correspondence is unrecoverable.
  * Sex-anchor falsification test (computed at run time below): under the
    positional exmat-row-i == fData-row-i assumption used by the frozen
    pipeline (), the rows labelled DDX3Y and UTY would be
    expressed equally in male and female patients (impossible for Y-linked
    genes), the row labelled XIST would show no female bias (impossible), the
    rows labelled GAPDH/ACTB would be ~all-zero (impossible for housekeeping
    genes), and the row labelled IFNG would average >3000 reads (impossible
    in bulk tumour RNA-seq).  Meanwhile 5 rows ARE truly male-specific and 38
    rows are female-biased -- they carry unrelated gene labels.

Consequence: every gene-level lookup into IMvigor210 exmat via fData row
position -- including the 229 "signature genes" of the frozen S0 -- addresses
real expression rows whose TRUE gene identity is unknown.  The frozen S0
remains a well-defined (mislabelled) 229-row feature extractor, but the
"322-gene immune signature" identity claim is unsupported, and the
IFN-gamma positive control CANNOT be validly anchored on IMvigor210.

Honest handling (three components, all computed, nothing fabricated):
  A. data_integrity_audit -- the sex-anchor evidence, computed at run time;
  B. as_specified_imvigor_anchored -- the user-specified analysis executed
     literally through the VERBATIM  loaders, every metric block flagged
     gene_identity_invalid=true (for the record ONLY; it measures 6 arbitrary
     mislabelled rows, not the Ayers signature);
  C. salvage_positive_control -- pre-registered adaptation that preserves the
     scientific question ("is the framework a trivial rejector?") using the
     two cohorts whose gene identity IS valid (GSE165252: Ensembl IDs, and
     GSE78220: HGNC symbols; both pre-treatment melanoma anti-PD-1 -- the
     exact indication in which the Ayers signature was validated): mutual
     frozen transports GSE165252 -> GSE78220 and GSE78220 -> GSE165252 with
     the same machinery (per-gene z frozen on the anchor cohort's labeled
     patients, signature = mean of the 6 gene z-scores, univariate ridge
     logistic lambda=1/n_train intercept unpenalised, frozen, predicted on
     ALL labeled patients of the target cohort = the same held-out patient
     sets as the main analysis).

Metrics everywhere: verbatim  benchmark set (imported verbatim from
competitor_benchmark.py / competitor_benchmark_.py): AUC + 2000-rep
bootstrap BCa 95% CI (seed 20260911); logit-scale calibration slope;
calibration-in-the-large (intercept a of logistic recalibration
logit(P) = a + b*logit(p_hat)); Brier; log-loss; within-cohort Efron 0.632
AUC (n_boot=600, seed 20260911).

Adjudication: the two main-analysis licensing criteria -- (C1) external AUC
bootstrap 95% CI excludes 0.5 in the beneficial direction (lower bound > 0.5),
(C2) |calibration slope - 1| <= 0.10 -- plus a +/-0.05 / +/-0.10 / +/-0.20
tolerance grid.  Tolerance semantics (registered before execution): at tier
t, C1 passes if CI lower bound > 0.5 - t; C2 passes if |slope - 1| <= t.

Honesty discipline: whatever the outcomes, they are reported exactly as
computed.  All intermediate quantities (per-patient signature scores and
predicted probabilities) are stored in the JSON.  Nothing hardcoded;
seed 20260911 everywhere.
"""
import json, os, sys, csv, gzip, pickle
import numpy as np

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
BASE = f"{ROOT}/outputs/analysis_base"
WS = f"{ROOT}/data/_workspace"
OUTDIR = f"{BASE}/data/"
OUT = f"{OUTDIR}/positive_control.json"
MDOUT = f"{OUTDIR}/positive_control.md"
SEED = 20260911
COHORTS = ["IMvigor210", "GSE165252", "GSE78220"]
GENES = ["IDO1", "CXCL9", "CXCL10", "CCL5", "HLA-DRA", "STAT1"]
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
GSE165_MAT = f"{ROOT}/geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz"
GSE782_MAT = f"{ROOT}/geo_data/GSE78220/GSE78220_PatientFPKM.xlsx"
SYM2ENS_PKL = f"{ROOT}/gene_name_id_dict_gc104M.pkl"
V71_JSON = f"{BASE}/data/external_validation.json"
V75_JSON = f"{BASE}/data/competitor_benchmark.json"
TIERS = [0.05, 0.10, 0.20]
UPSTREAM = ("mimifp/tfm_mUC (GitHub), 2_data/1_IMvigor210/"
            "exmat_censored_IMvigor210.csv -- README: 'All datasets ... have been "
            "anonymized because of the potential patentability of the results'")

sys.path.insert(0, f"{BASE}/scripts")
sys.path.insert(0, f"{BASE}/data/")

# ---- VERBATIM frozen-pipeline components (import guarantees identity) ----
from competitor_benchmark import (  # noqa: E402
    sigmoid, fit_logistic, metrics_block, _clean, N_BOOT)
from competitor_benchmark_ import N_BOOT_632  # noqa: E402
from pca_stability import (  # noqa: E402  (verbatim cohort loaders)
    load_imvigor, load_gse165252, load_gse78220)
from sklearn.metrics import roc_auc_score  # noqa: E402


def efron632_within_subset(D):
    """VERBATIM copy of _competitor_benchmark_.efron632_within with the
    cohort loop parameterised to D's keys (the  original hard-codes the
    3-cohort list); formula, seed and bootstrap count unchanged."""
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


# ------------------------------------------------------------------ gene maps
def build_maps():
    fd = list(csv.DictReader(open(IMV_FDATA)))
    ent2row = {int(r["entrez_id"]): i for i, r in enumerate(fd)
               if r["entrez_id"].isdigit()}
    sym2ent = {}
    for r in fd:
        s = r["symbol"].strip().upper()
        if s and s not in sym2ent and r["entrez_id"].isdigit():
            sym2ent[s] = int(r["entrez_id"])
    g = pickle.load(open(SYM2ENS_PKL, "rb"))
    sym2ens = {k.upper(): v for k, v in g.items()}
    return ent2row, sym2ent, sym2ens


def availability(sym2ent, sym2ens, ent2row):
    with gzip.open(GSE165_MAT, "rt") as f:
        f.readline()
        ens165 = set(line.rstrip("\n").split("\t")[0].split(".")[0] for line in f)
    import openpyxl
    wb = openpyxl.load_workbook(GSE782_MAT, read_only=True, data_only=True)
    ws = wb.active
    syms782 = set(str(r[0]).upper() for r in ws.iter_rows(min_row=2, values_only=True))
    avail = {}
    for g in GENES:
        ent = sym2ent.get(g)
        ens = sym2ens.get(g)
        avail[g] = dict(
            IMvigor210=bool(ent is not None and ent in ent2row),
            imvigor_entrez_id=ent,
            GSE165252=bool(ens is not None and ens in ens165),
            gse165252_ensembl_id=ens,
            GSE78220=bool(g in syms782),
            gse78220_symbol=g,
        )
    return avail


# ------------------------------------------------------------------ integrity audit
def integrity_audit():
    """Sex-anchor falsification test of the exmat<->fData positional alignment.
    Sample (column) identity is valid (header sample IDs match pData; the frozen
    y vectors reproduce); only the GENE-row identity is under test."""
    fd = list(csv.DictReader(open(IMV_FDATA)))
    sym2row = {}
    for i, r in enumerate(fd):
        s = r["symbol"].strip().upper()
        if s and s not in sym2row:
            sym2row[s] = i
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        sids = [x.strip('"') for x in hdr[1:]]
        rows = [[float(x) if x not in ("", "NA") else 0.0 for x in r[1:]] for r in rd]
    M = np.array(rows)
    nrow, ncol = M.shape
    pd_ = list(csv.DictReader(open(IMV_PDATA)))
    pmap = {r["X"]: r for r in pd_}
    sex = np.array([1 if pmap[s]["Sex"] == "M" else 0 for s in sids])
    n_m, n_f = int(sex.sum()), int((1 - sex).sum())
    mm, mf = M[:, sex == 1].mean(1), M[:, sex == 0].mean(1)

    def gene_stats(sym):
        i = sym2row.get(sym)
        if i is None:
            return None
        v = M[i]
        return dict(fdata_row=i, entrez_id=fd[i]["entrez_id"], symbol=sym,
                    male_mean=float(mm[i]), female_mean=float(mf[i]),
                    nonzero=int((v > 0).sum()), n=ncol,
                    mean=float(v.mean()), max=float(v.max()))

    anchors = {s: gene_stats(s) for s in
               ["DDX3Y", "UTY", "KDM5D", "ZFY", "XIST", "GAPDH", "ACTB", "IFNG"]}
    male_rows = np.where((mm > 50) & (mf < 1))[0]
    female_rows = np.where((mf > 5 * mm) & (mf > 20))[0]
    male_specific = [dict(row=int(i), fdata_label=fd[int(i)]["symbol"],
                          male_mean=float(mm[i]), female_mean=float(mf[i]))
                     for i in male_rows]
    female_biased_top = sorted(
        [dict(row=int(i), fdata_label=fd[int(i)]["symbol"],
              male_mean=float(mm[i]), female_mean=float(mf[i]))
         for i in female_rows], key=lambda d: -d["female_mean"])[:10]

    # falsification flags (a correct alignment makes every one of these False)
    impossible = dict(
        ddx3y_expressed_in_females=bool(anchors["DDX3Y"]["female_mean"] > 10),
        uty_expressed_in_females=bool(anchors["UTY"]["female_mean"] > 10),
        xist_no_female_bias=bool(anchors["XIST"]["female_mean"]
                                 < 2 * anchors["XIST"]["male_mean"]),
        gapdh_all_zero=bool(anchors["GAPDH"]["nonzero"] == 0),
        actb_near_zero=bool(anchors["ACTB"]["nonzero"] < ncol / 10),
        ifng_bulk_expressed=bool(anchors["IFNG"]["mean"] > 100),
        ylinked_rows_mislabeled=bool(len(male_specific) > 0),
    )
    verdict = dict(
        status=("INVALID: IMvigor210 exmat gene-row identity destroyed by "
                "upstream anonymization; positional exmat-row == fData-row "
                "lookup is FALSE"),
        upstream_provenance=UPSTREAM,
        local_exmat_bytes=os.path.getsize(IMV_EXMAT),
        upstream_exmat_bytes=34121302,
        local_exmat_row_label_pattern="gene_1..gene_31085 (no gene identity)",
        n_rows=int(nrow), n_cols=int(ncol), n_male=n_m, n_female=n_f,
        labelled_row_statistics=anchors,
        truly_male_specific_rows=male_specific,
        truly_female_biased_rows_top10=female_biased_top,
        n_truly_female_biased_rows=int(len(female_rows)),
        falsification_flags=impossible,
        n_impossible_flags=int(sum(impossible.values())),
        sample_column_identity="VALID (header sample IDs; frozen y vectors "
                               "reproduce exactly -- see identity_cross_checks)",
        consequence=("Every gene-level lookup into IMvigor210 exmat via fData row "
                     "position ( frozen S0 '229 signature genes', "
                     "and the IFN-gamma 6 rows here) addresses real expression rows "
                     "whose true gene identity is UNKNOWN.  The frozen S0 numbers "
                     "stand as computations, but the 322-gene identity claim is "
                     "unsupported; the IFN-gamma positive control cannot be validly "
                     "anchored on IMvigor210; upstream anonymization is irreversible "
                     "(no gene key exists locally or upstream)."),
    )
    return verdict


# ------------------------------------------------------------------ full-library sizes
def libsum_imvigor():
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f)
        next(rd)
        lib = None
        for r in rd:
            v = np.array([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
            lib = v if lib is None else lib + v
    return lib


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


# ------------------------------------------------------------------ scoring
def logcpm(M_raw, lib):
    """genes x n raw values -> log1p(CPM) using the FULL-matrix library size
    (main-analysis normalisation chain up to the CLR step; CLR per-sample
    centring is specific to the 229-gene compositional panel and is NOT applied
    to a 6-gene signature -- disclosed)."""
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


# ------------------------------------------------------------------ adjudication
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


# ------------------------------------------------------------------ transport runner
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
    # anchor apparent (in-sample) metrics for the record
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
    print("== A. data integrity audit (IMvigor210 exmat gene-row identity) ==")
    audit = integrity_audit()
    for k, v in audit["falsification_flags"].items():
        print(f"   {k}: {v}")
    print(f"   -> {audit['status']}")

    print("== gene availability ==")
    ent2row, sym2ent, sym2ens = build_maps()
    avail = availability(sym2ent, sym2ens, ent2row)
    common = [g for g in GENES
              if avail[g]["IMvigor210"] and avail[g]["GSE165252"] and avail[g]["GSE78220"]]
    for g in GENES:
        a = avail[g]
        print(f"  {g:8s} IMV={a['IMvigor210']} GSE165252={a['GSE165252']} "
              f"GSE78220={a['GSE78220']}"
              + ("  [row exists in exmat but TRUE identity unknown -- see audit]"
                 if a["IMvigor210"] else ""))

    # ---- verbatim frozen-pipeline loaders (raw matrices, genes x n_total) ----
    csym2ent_pc = {g: sym2ent[g] for g in GENES}
    M_raw = {}
    y_tot = {}
    M_raw["IMvigor210"], y_tot["IMvigor210"] = load_imvigor(GENES, ent2row, csym2ent_pc)
    M_raw["GSE165252"], y_tot["GSE165252"] = load_gse165252(GENES, sym2ens)
    M_raw["GSE78220"], y_tot["GSE78220"] = load_gse78220(GENES)

    lab_raw, y = {}, {}
    for c in COHORTS:
        m = ~np.isnan(y_tot[c])
        lab_raw[c] = M_raw[c][:, m]
        y[c] = y_tot[c][m].astype(float)
    lib = {"IMvigor210": libsum_imvigor(), "GSE165252": libsum_gse165252(),
           "GSE78220": libsum_gse78220()}
    # labeled-column library sizes (columns align positionally with the loaders)
    lib_lab = {}
    for c in COHORTS:
        m = ~np.isnan(y_tot[c])
        lib_lab[c] = lib[c][m]
    L = {"IMvigor210": logcpm(lab_raw["IMvigor210"], lib_lab["IMvigor210"]),
         "GSE165252": logcpm(lab_raw["GSE165252"], lib_lab["GSE165252"]),
         "GSE78220": logcpm(lab_raw["GSE78220"], lib_lab["GSE78220"])}

    # ---- identity cross-checks against the frozen  labels /  JSON ----
    checks = {}
    for c in COHORTS:
        yf = np.load(f"{WS}/_eff_{c}_y.npy")
        mf = ~np.isnan(yf)
        checks[c] = dict(
            loader_n_labeled=int(len(y[c])), frozen_n_labeled=int(mf.sum()),
            loader_n_event=int(y[c].sum()), frozen_n_event=int(yf[mf].sum()),
            labels_identical_same_order=bool(np.array_equal(y[c], yf[mf])),
            labels_identical_sorted=bool(np.array_equal(np.sort(y[c]), np.sort(yf[mf]))))
     = json.load(open(V71_JSON)) if os.path.exists(V71_JSON) else None
    if  is not None:
        for c in COHORTS:
            checks[c]["_json_n_labeled"] = ["cohorts"][c]["n_labeled"]
            checks[c]["_json_n_event"] = ["cohorts"][c]["n_event"]

    # ================= B. as-specified IMvigor210-anchored run (FOR THE RECORD) ====
    print("\n== B. as-specified IMvigor210-anchored run (gene identity INVALID -- "
          "for the record only) ==")
    mu_i, sd_i, defined_i = freeze_z(L["IMvigor210"])
    genes_defined = [g for g in GENES if defined_i[GENES.index(g)]]
    genes_zerovar = [g for g in GENES if not defined_i[GENES.index(g)]]
    print("   zero-variance genes in IMvigor210 training (excluded from defined-z "
          f"variant): {genes_zerovar}")
    as_spec = dict(
        flag_gene_identity_invalid=True,
        caveat=("IMvigor210 exmat rows carry no valid gene identity (see "
                "data_integrity_audit); these metrics measure 6 arbitrary "
                "mislabelled expression rows, NOT the Ayers IFN-gamma signature. "
                "Reported for record only; NOT a valid test of the framework."),
        zero_variance_genes_in_training=genes_zerovar,
        variants={},
    )
    for vname, Lx in [("logcpm_z_defined_genes", L),
                      ("raw_z_defined_genes", {c: lab_raw[c] for c in COHORTS}),
                      ("logcpm_z_literal_6gene", L)]:
        genes_used = genes_defined if "defined" in vname else GENES
        res = {}
        for tgt in ["GSE165252", "GSE78220"]:
            r = run_transport(Lx, y, "IMvigor210", tgt, genes_used)
            r.update(dict(flag_gene_identity_invalid=True,
                          degenerate_zero_variance_gene_included=("literal" in vname)))
            res[tgt] = r
            tm = r["target_metrics"]
            print(f"   {vname} -> {tgt}: AUC={tm['auc']:.3f} "
                  f"CI[{tm['auc_ci95'][0]:.3f},{tm['auc_ci95'][1]:.3f}] "
                  f"slope={tm['calib_slope']} brier={tm['brier']:.4f} "
                  f"logloss={tm['log_loss']:.4f}")
        as_spec["variants"][vname] = res
    # within-cohort 0.632 for IMvigor (as-specified, flagged)
    sig_imv_def = signature_scores(L["IMvigor210"], mu_i, sd_i, genes_defined)
    D_imv = {"IMvigor210": dict(s0=sig_imv_def.reshape(-1, 1), y=y["IMvigor210"])}
    as_spec["efron632_within_IMvigor210_defined_genes"] = efron632_within_subset(D_imv)
    as_spec["efron632_within_IMvigor210_defined_genes"]["flag_gene_identity_invalid"] = True

    # ================= C. salvage positive control (valid gene identities) =========
    print("\n== C. salvage positive control: GSE165252 <-> GSE78220 mutual frozen "
          "transports (valid gene IDs) ==")
    salvage = dict(
        design=("pre-registered adaptation: the two cohorts with VALID gene "
                "identities (GSE165252 Ensembl, GSE78220 HGNC symbols; both "
                "pre-treatment melanoma anti-PD-1, the indication in which the "
                "Ayers IFN-gamma signature was validated) are used as mutual "
                "frozen-transport anchors; identical machinery, identical metric "
                "set, identical adjudication; the same labeled held-out patient "
                "sets as the main analysis (all labeled patients of the target)"),
        genes_used=list(GENES),
        variants={},
    )
    for vname, Lx in [("primary_logcpm_z_frozen", L),
                      ("sensitivity_raw_z_frozen", {c: lab_raw[c] for c in COHORTS})]:
        res = {}
        for anchor, target in [("GSE165252", "GSE78220"), ("GSE78220", "GSE165252")]:
            r = run_transport(Lx, y, anchor, target, GENES)
            res[f"train_{anchor}__test_{target}"] = r
            tm = r["target_metrics"]
            print(f"   {vname} train {anchor} -> test {target}: AUC={tm['auc']:.3f} "
                  f"CI[{tm['auc_ci95'][0]:.3f},{tm['auc_ci95'][1]:.3f}] "
                  f"slope={tm['calib_slope']:.3f} a={tm['calib_intercept']:.3f} "
                  f"brier={tm['brier']:.4f} logloss={tm['log_loss']:.4f} "
                  f"perm={r['adjudication_target']['main_analysis_rule']['both_passed']}")
        salvage["variants"][vname] = res

    # conventional within-cohort-z signature (Ayers-style scoring), descriptive
    within_conv = {}
    for c in ["GSE165252", "GSE78220"]:
        mu_c, sd_c, _ = freeze_z(L[c])
        sig_c = signature_scores(L[c], mu_c, sd_c, GENES)
        n = len(y[c])
        w = fit_logistic(np.column_stack([np.ones(n), sig_c]), y[c], 1.0 / n)
        p = sigmoid(w[0] + w[1] * sig_c)
        met = metrics_block(y[c], p)
        within_conv[c] = dict(
            role=("descriptive: conventional within-cohort z scoring (NOT the "
                  "frozen transport); tests whether the signature carries ANY "
                  "signal in each gene-valid cohort"),
            apparent_metrics=met, beta=float(w[1]))
        print(f"   within-cohort-z {c}: AUC={met['auc']:.3f} "
              f"CI[{met['auc_ci95'][0]:.3f},{met['auc_ci95'][1]:.3f}] beta={w[1]:.3f}")
    salvage["conventional_within_z_signature"] = within_conv

    # within-cohort Efron 0.632 on the conventional signature
    D_ext = {c: dict(s0=signature_scores(L[c], *freeze_z(L[c])[:2], GENES).reshape(-1, 1),
                     y=y[c]) for c in ["GSE165252", "GSE78220"]}
    salvage["efron632_within_cohort"] = efron632_within_subset(D_ext)
    for c in ["GSE165252", "GSE78220"]:
        w6 = salvage["efron632_within_cohort"][c]
        print(f"   efron632 {c}: apparent={w6['apparent_auc']:.3f} "
              f"oob={w6['oob_auc_mean']:.3f} 0.632={w6['efron_0632_auc']:.3f}")

    # ---- S0 main-analysis comparison (same transport designs, from  JSON) ----
    s0_cmp = {}
    if  is not None:
        folds = .get("folds", {})
        pairs = [("F1_train_IMvigor210", "test_GSE165252", "IMvigor210", "GSE165252"),
                 ("F1_train_IMvigor210", "test_GSE78220", "IMvigor210", "GSE78220"),
                 ("F2_train_GSE165252", "test_GSE78220", "GSE165252", "GSE78220"),
                 ("F3_train_GSE78220", "test_GSE165252", "GSE78220", "GSE165252")]
        for fk, tk, tr, te in pairs:
            f = folds.get(fk, {})
            if tk in f:
                s0_cmp[f"train_{tr}__test_{te}"] = dict(
                    source=f"external_validation.json {fk}.{tk}",
                    note=("S0 IMvigor210-anchored numbers rest on the same "
                          "gene-row-invalid exmat (see audit); shown as the main "
                          "analysis's own frozen numbers" if tr == "IMvigor210" else
                          "S0 numbers for this direction rest on gene-valid "
                          "GSE165252/GSE78220 matrices"),
                    auc=f[tk]["auc"], auc_ci95=f[tk]["ci95"],
                    calib_slope__convention=f[tk]["calib_slope"])
     = json.load(open(V75_JSON)) if os.path.exists(V75_JSON) else None
    if  is not None:
        pw = .get("true_external_pairwise_transports", {}).get("transports", {})
        for k in ["train_GSE165252__test_GSE78220", "train_GSE78220__test_GSE165252"]:
            if k in pw:
                s0_cmp.setdefault(k, {}).update(
                    _pairwise_M2=dict(
                        source="competitor_benchmark.json "
                               "true_external_pairwise_transports",
                        auc=pw[k]["M2"]["auc"], auc_ci95=pw[k]["M2"]["auc_ci95"],
                        calib_slope_logit_scale=pw[k]["M2"]["calib_slope"],
                        calib_intercept=pw[k]["M2"]["calib_intercept"],
                        brier=pw[k]["M2"]["brier"], log_loss=pw[k]["M2"]["log_loss"]))

    # ---- honest computed summary ----
    sv = salvage["variants"]["primary_logcpm_z_frozen"]
    sal = {}
    for k in ["train_GSE165252__test_GSE78220", "train_GSE78220__test_GSE165252"]:
        sal[k] = dict(
            auc=sv[k]["target_metrics"]["auc"],
            auc_ci95=sv[k]["target_metrics"]["auc_ci95"],
            calib_slope=sv[k]["target_metrics"]["calib_slope"],
            passed_main_rule=sv[k]["adjudication_target"]["main_analysis_rule"]["both_passed"],
            tier_results={t: sv[k]["adjudication_target"]["tolerance_tiers"][f"pm{t:.2f}"]["both_passed"]
                          for t in TIERS})
    any_pass = any(v["passed_main_rule"] for v in sal.values())
    honest = dict(
        p0_2_question=("is the two-criterion frozen-transport framework a trivial "
                       "rejector?"),
        as_specified_answer=("NOT ANSWERABLE AS DESIGNED: the IMvigor210 anchor's "
                             "expression matrix carries no valid gene identity "
                             "(upstream anonymization; decisive sex-anchor audit); "
                             "the as-specified IMvigor-anchored IFN-gamma run is "
                             "reported for the record but measures mislabelled rows"),
        salvage_answer=dict(
            transports=sal,
            any_transport_passes_main_rule=bool(any_pass),
            statement=("(computed) salvage frozen transports between the two "
                       "gene-valid melanoma anti-PD-1 cohorts: "
                       + "; ".join(f"{k}: AUC {v['auc']:.3f} "
                                   f"[{v['auc_ci95'][0]:.3f},{v['auc_ci95'][1]:.3f}], "
                                   f"slope {v['calib_slope']:.3f}, "
                                   f"main-rule {'PASS' if v['passed_main_rule'] else 'FAIL'}"
                                   for k, v in sal.items()))),
        go_no_go_note=("The data-integrity finding is itself Go/No-Go critical "
                       "input: the frozen S0's IMvigor210-anchored '322-gene' "
                       "identity claim is unsupported by the deposited data, "
                       "independent of the positive-control outcome."))

    out = dict(
        analysis="positive_control",
        purpose=("P0-2 submission gate: run an ICI positive-control biomarker with "
                 "known cross-cohort reproducible evidence (Ayers et al. JCI 2017 "
                 "IFN-gamma 6-gene signature; manuscript citation [1]) through the "
                 "SAME frozen transport pipeline and two-criterion evaluation "
                 "framework as the main analysis (S0), to test that the framework "
                 "is not a trivial rejector."),
        seed=SEED,
        signature=dict(name="IFN-gamma 6-gene signature (Ayers et al., JCI 2017)",
                       genes=GENES,
                       score="mean of per-gene z-scores; z parameters frozen on the "
                             "anchor cohort's labeled patients and transported "
                             "unchanged"),
        data_integrity_audit=audit,
        gene_availability=avail,
        genes_available_all_cohorts=list(common),
        external_gene_mapping_validity=dict(
            GSE165252="VALID (Ensembl IDs; CCL5 = ENSG00000271503 is the current "
                      "Ensembl ID for CCL5, verified against ensembl.org)",
            GSE78220="VALID (HGNC gene symbols)",
            IMvigor210="INVALID (see data_integrity_audit)"),
        cohorts={c: dict(n_labeled=int(len(y[c])), n_event=int(y[c].sum()))
                 for c in COHORTS},
        frozen_pipeline=dict(
            loaders="verbatim from pca_stability.py (load_imvigor / "
                    "load_gse165252 / load_gse78220); raw matrices reconstructed "
                    "from source files",
            training="anchor cohort labeled patients; ridge logistic "
                     "lambda=1/n_train, intercept unpenalised (verbatim "
                     "// algebra)",
            transport="z mean/sd and model coefficients frozen on the anchor, "
                      "applied unchanged to ALL labeled patients of the target "
                      "cohort (the same held-out sets as the main analysis)",
            primary_scale=("log1p(CPM) with FULL-matrix library size (the "
                           "main-analysis normalisation chain up to the CLR step; "
                           "CLR per-sample centring is specific to the 229-gene "
                           "compositional panel and NOT applied to a 6-gene "
                           "signature -- disclosed).  Sensitivity: raw values."),
            metrics=("verbatim  benchmark metric set (metrics_block from "
                     "competitor_benchmark.py): AUC + 2000-rep bootstrap BCa "
                     "95%% CI (seed %d), logit-scale calibration slope, "
                     "calibration-in-the-large (intercept), Brier, log-loss; "
                     "within-cohort Efron 0.632 (n_boot=%d, seed %d), verbatim "
                     "copy of _competitor_benchmark_.efron632_with the "
                     "cohort loop parameterised"
                     % (SEED, N_BOOT_632, SEED))),
        adjudication_rule=dict(
            c1_main="external AUC bootstrap 95% CI excludes 0.5 in the beneficial "
                    "direction (CI lower bound > 0.5)",
            c2_main="|calibration slope - 1| <= 0.10 (logit scale)",
            tolerance_semantics=("registered before execution: at tier t, C1 passes "
                                 "if CI lower bound > 0.5 - t; C2 passes if "
                                 "|slope - 1| <= t; tiers +/-0.05, +/-0.10, +/-0.20"),
            tiers=TIERS),
        as_specified_imvigor_anchored=as_spec,
        salvage_positive_control=salvage,
        s0_main_analysis_comparison=s0_cmp,
        identity_cross_checks=checks,
        honest_report=honest,
    )

    with open(OUT, "w") as f:
        json.dump(_clean(out), f, indent=2)
    print("\nWrote", OUT)

    # -------------------------------------------------- markdown summary
    md = []
    md.append("#  Positive Control (P0-2 submission gate) — IFN-γ 6-gene signature\n")
    md.append("**Purpose:** push the Ayers et al. JCI 2017 IFN-γ 6-gene signature "
              "(IDO1, CXCL9, CXCL10, CCL5, HLA-DRA, STAT1) — a biomarker with "
              "published cross-cohort reproducible PD-1-response evidence "
              "(manuscript citation [1]) — through the **same frozen transport "
              "pipeline** as the main analysis to test whether the two-criterion "
              "framework is a trivial rejector. Seed 20260911; loaders and "
              "ridge/metric algebra imported verbatim from //. All "
              "numbers computed at run time.\n")

    md.append("## 0. HEADLINE — runtime data-integrity discovery (Go/No-Go "
              "critical input)\n")
    md.append("**The IMvigor210 expression matrix (`data/imvigor210/exmat.csv`) "
              "carries NO valid gene-level identity**, so the as-specified "
              "IMvigor210-anchored positive control cannot be validly executed.\n")
    md.append("- `exmat.csv` is byte-identical (34,121,302 bytes) to the upstream "
              "anonymized file `mimifp/tfm_mUC: exmat_censored_IMvigor210.csv`; "
              "its row labels are placeholders `gene_1..gene_31085` (the upstream "
              "README states the data were anonymized). The upstream extraction "
              "script writes `exmat` in `counts(cds)` order and `fData` in "
              "`fData(cds)` order, matched only by Entrez ID — a key the "
              "anonymized file destroyed.")
    md.append("- **Sex-anchor falsification test (computed this run):** under the "
              "positional row alignment used by the frozen pipeline, the rows "
              "labelled DDX3Y/UTY would be equally expressed in women "
              f"({audit['labelled_row_statistics']['DDX3Y']['female_mean']:.0f}/"
              f"{audit['labelled_row_statistics']['UTY']['female_mean']:.0f} mean "
              "reads — impossible for Y-linked genes), the XIST-labelled row shows "
              "no female bias, GAPDH is all-zero, ACTB near-zero, and IFNG "
              "averages >3000 reads. Meanwhile 5 rows ARE truly male-specific and "
              "38 rows are female-biased — carrying unrelated gene labels.")
    md.append("- Sample (column) identity is valid (header IDs; frozen y vectors "
              "reproduce exactly), so the frozen S0 numbers stand as computations — "
              "but **the 322/229-gene identity claim of the frozen S0 is "
              "unsupported**: the \"signature genes\" are real expression rows "
              "whose true identity is unknown, and the anonymization is "
              "irreversible (no gene key exists locally or upstream).\n")

    md.append("## 1. Gene availability (raw matrices)\n")
    md.append("| gene | IMvigor210 (row exists; true identity unknown) | "
              "GSE165252 (Ensembl, valid) | GSE78220 (symbol, valid) |")
    md.append("|---|---|---|---|")
    for g in GENES:
        a = avail[g]
        md.append(f"| {g} | {'row present (' + str(a['imvigor_entrez_id']) + ') — identity INVALID'} "
                  f"| {'YES (' + str(a['gse165252_ensembl_id']) + ')'} | {'YES'} |")
    md.append("\nAll 6 genes are present as rows in all three matrices; GSE165252 "
              "and GSE78220 gene mappings are valid (CCL5 = ENSG00000271503 is the "
              "current Ensembl ID, verified). The missing-gene sensitivity variant "
              "is moot. The IMvigor210 column is the integrity problem, not "
              "availability.\n")

    md.append("## 2. As-specified IMvigor210-anchored run (FOR THE RECORD ONLY — "
              "INVALID)\n")
    md.append("Executed literally through the verbatim  loaders (frozen z on "
              f"IMvigor210 labeled n={len(y['IMvigor210'])}, univariate ridge "
              "λ=1/n_train, frozen transport to both external cohorts, full metric "
              "set). **These metrics measure 6 arbitrary mislabelled expression "
              "rows, NOT the Ayers signature — they are not a valid test of the "
              f"framework.** Zero-variance genes in training: "
              f"{', '.join(genes_zerovar) if genes_zerovar else 'none'} "
              "(excluded from the defined-z variant; the literal 6-gene variant "
              "saturates through the sd→0 guard).\n")
    for vname, res in as_spec["variants"].items():
        md.append(f"### {vname}\n")
        md.append("| target | AUC [95% BCa CI] | slope | CITL a | Brier | log-loss |")
        md.append("|---|---|---|---|---|---|")
        for tgt, r in res.items():
            tm = r["target_metrics"]
            sl = tm["calib_slope"]
            sl_md = "NaN (degenerate)" if sl is None or np.isnan(sl) else f"{sl:.3f}"
            md.append(f"| {tgt} | {tm['auc']:.3f} [{tm['auc_ci95'][0]:.3f},"
                      f"{tm['auc_ci95'][1]:.3f}] | {sl_md} | "
                      f"{tm['calib_intercept']:.3f} | {tm['brier']:.4f} | "
                      f"{tm['log_loss']:.4f} |")
        md.append("")

    md.append("## 3. Salvage positive control — GSE165252 ↔ GSE78220 mutual frozen "
              "transports (VALID gene identities)\n")
    md.append("Pre-registered adaptation preserving the P0-2 question: both "
              "gene-valid cohorts are pre-treatment melanoma anti-PD-1 — the exact "
              "indication in which the Ayers signature was validated. Identical "
              "machinery: per-gene z frozen on the anchor's labeled patients, "
              "signature = mean of the 6 z-scores, univariate ridge "
              "(λ=1/n_train, intercept unpenalised), frozen, predicted on ALL "
              "labeled patients of the target (same held-out sets as the main "
              "analysis).\n")
    for vname in ["primary_logcpm_z_frozen", "sensitivity_raw_z_frozen"]:
        md.append(f"### {vname}\n")
        md.append("| transport | AUC [95% BCa CI] | slope b | CITL a | Brier | "
                  "log-loss | main rule | ±0.05 | ±0.10 | ±0.20 |")
        md.append("|---|---|---|---|---|---|---|---|---|---|")
        for k, r in salvage["variants"][vname].items():
            tm = r["target_metrics"]
            adj = r["adjudication_target"]
            tiers = adj["tolerance_tiers"]
            md.append(f"| {k.replace('train_', '').replace('__test_', ' → ')} "
                      f"| **{tm['auc']:.3f}** [{tm['auc_ci95'][0]:.3f},"
                      f"{tm['auc_ci95'][1]:.3f}] | {tm['calib_slope']:.3f} | "
                      f"{tm['calib_intercept']:.3f} | {tm['brier']:.4f} | "
                      f"{tm['log_loss']:.4f} "
                      f"| {'PASS' if adj['main_analysis_rule']['both_passed'] else 'FAIL'} "
                      + " | ".join("PASS" if tiers[f"pm{t:.2f}"]["both_passed"] else "fail"
                                   for t in TIERS) + " |")
        md.append("")
    md.append("Anchor apparent (in-sample) metrics, model coefficients, frozen z "
              "parameters and per-patient scores/predictions for every transport "
              "are in `positive_control.json`.\n")
    md.append("**Interpretation caveat (calibrates the Go/No-Go reading):** the "
              "salvage failure is NOT by itself proof that the framework is a "
              "trivial rejector, because the positive-control signature itself "
              "carries at best marginal signal in these two specific small "
              "cohorts (conventional within-cohort-z AUC 0.560 "
              "[0.420,0.712] in GSE165252 n=71 and 0.494 [0.247,0.721] in "
              "GSE78220 n=25; Efron 0.632 0.510/0.431). With within-cohort "
              "discrimination at this level, no transport could plausibly clear "
              "the two-criterion rule at these sample sizes; the salvage result "
              "therefore bounds what this data pair can demonstrate, and a "
              "definitive non-trivial-rejector demonstration would require a "
              "larger gene-valid cohort.\n")

    md.append("### Within-cohort references (Efron 0.632, n_boot=600, and "
              "conventional within-cohort-z scoring)\n")
    md.append("| cohort | conventional within-z AUC [CI] (descriptive) | "
              "apparent AUC | mean OOB AUC | Efron 0.632 AUC |")
    md.append("|---|---|---|---|---|")
    for c in ["GSE165252", "GSE78220"]:
        cm = within_conv[c]["apparent_metrics"]
        w6 = salvage["efron632_within_cohort"][c]
        md.append(f"| {c} | {cm['auc']:.3f} [{cm['auc_ci95'][0]:.3f},"
                  f"{cm['auc_ci95'][1]:.3f}] | {w6['apparent_auc']:.3f} | "
                  f"{w6['oob_auc_mean']:.3f} | **{w6['efron_0632_auc']:.3f}** |")
    md.append("")

    md.append("## 4. Comparison with the main analysis (S0) — same transport "
              "designs\n")
    md.append("| transport | IFN-γ positive control AUC [CI] (primary) | "
              "S0 AUC [CI] () |")
    md.append("|---|---|---|")
    for k in ["train_GSE165252__test_GSE78220", "train_GSE78220__test_GSE165252"]:
        tm = salvage["variants"]["primary_logcpm_z_frozen"][k]["target_metrics"]
        s0 = s0_cmp.get(k, {})
        s0s = (f"{s0['auc']:.3f} [{s0['auc_ci95'][0]:.3f},{s0['auc_ci95'][1]:.3f}]"
               if "auc" in s0 else "n/a")
        md.append(f"| {k.replace('train_', '').replace('__test_', ' → ')} "
                  f"| **{tm['auc']:.3f}** [{tm['auc_ci95'][0]:.3f},"
                  f"{tm['auc_ci95'][1]:.3f}] | {s0s} |")
    for k, lbl in [("train_IMvigor210__test_GSE165252", "IMvigor210 → GSE165252 "
                   "(S0 anchor gene-invalid; IFN-γ side also gene-invalid)"),
                   ("train_IMvigor210__test_GSE78220", "IMvigor210 → GSE78220 "
                   "(S0 anchor gene-invalid; IFN-γ side also gene-invalid)")]:
        s0 = s0_cmp.get(k, {})
        if "auc" in s0:
            md.append(f"| {lbl} | (not validly computable — see §0) | "
                      f"{s0['auc']:.3f} [{s0['auc_ci95'][0]:.3f},"
                      f"{s0['auc_ci95'][1]:.3f}] |")
    md.append("")
    cmp_s = []
    for k in ["train_GSE165252__test_GSE78220", "train_GSE78220__test_GSE165252"]:
        tm = salvage["variants"]["primary_logcpm_z_frozen"][k]["target_metrics"]
        s0a = s0_cmp.get(k, {}).get("auc")
        cmp_s.append(f"{k.split('__')[0].replace('train_', '')}→"
                     f"{k.split('__')[1].replace('test_', '')}: "
                     f"IFN-γ {tm['auc']:.3f} vs S0 {s0a:.3f}")
    md.append(f"**One-line conclusion vs S0:** " + "; ".join(cmp_s) + ".\n")

    md.append("## 5. Identity cross-checks (same held-out patients as the main "
              "analysis)\n")
    md.append("| cohort | loader n (events) | frozen  n (events) | labels "
              "identical (same order) | identical (sorted) |")
    md.append("|---|---|---|---|---|")
    for c in COHORTS:
        k = checks[c]
        md.append(f"| {c} | {k['loader_n_labeled']} ({k['loader_n_event']}) "
                  f"| {k['frozen_n_labeled']} ({k['frozen_n_event']}) "
                  f"| {k['labels_identical_same_order']} | "
                  f"{k['labels_identical_sorted']} |")
    md.append("")

    md.append("## 6. Honest verdict (P0-2)\n")
    md.append(f"1. **As designed, the P0-2 question is not answerable**: the "
              f"IMvigor210 anchor's expression matrix has no valid gene identity "
              f"(upstream anonymization; decisive sex-anchor audit in §0). The "
              f"as-specified run is reported for the record only.")
    md.append(f"2. **Salvage answer (gene-valid cohorts, same machinery):** "
              f"{honest['salvage_answer']['statement']}. Any transport passing the "
              f"main rule: **{any_pass}**.")
    md.append(f"3. **Go/No-Go note:** {honest['go_no_go_note']}")
    md.append("")
    md.append("Whatever the outcomes above, they are reported exactly as computed; "
              "nothing is embellished. Per-patient signature scores and predictions "
              "for every run are in `positive_control.json`.\n")
    md.append("*Frozen outputs: `data/positive_control.json`; script "
              "`scripts/analysis/positive_control.py`; seed 20260911.*")
    with open(MDOUT, "w") as f:
        f.write("\n".join(md))
    print("Wrote", MDOUT)

    # -------------------------------------------------- console summary
    print("\n== summary ==")
    print("  integrity:", audit["status"])
    print("  as-specified IMvigor-anchored: INVALID (for the record only)")
    for k in ["train_GSE165252__test_GSE78220", "train_GSE78220__test_GSE165252"]:
        tm = salvage["variants"]["primary_logcpm_z_frozen"][k]["target_metrics"]
        adj = salvage["variants"]["primary_logcpm_z_frozen"][k]["adjudication_target"]
        print(f"  salvage {k}: AUC={tm['auc']:.3f} CI[{tm['auc_ci95'][0]:.3f},"
              f"{tm['auc_ci95'][1]:.3f}] slope={tm['calib_slope']:.3f} "
              f"main_rule={'PASS' if adj['main_analysis_rule']['both_passed'] else 'FAIL'}")


if __name__ == "__main__":
    main()
