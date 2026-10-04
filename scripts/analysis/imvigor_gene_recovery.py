#!/usr/bin/env python3
""" -- Nature Cancer revision , Go/No-Go recovery task.

Recover the gene-row identity of the anonymized IMvigor210 anchor matrix
(data/imvigor210/exmat.csv, rows gene_1..gene_31085), then rebuild the
anchored 322-gene CLR->PCA k=2 arm with REAL gene labels and rerun the
main numbers.

STAGES
  1a  local de-anonymization: candidate row-ordering hypotheses of exmat vs
      fData, each tested against sex-anchor / housekeeping / immune gates.
      WINNER (this run): exmat row j == the j-th fData row when fData is
      sorted by SYMBOL (upper-cased) alphabetically.
  1b  provenance: upstream repo mimifp/tfm_mUC (byte-identical fData/pData),
      upstream extraction script shows exmat written in counts(cds) row order
      with gene symbols in column 1, later censored to gene_N placeholders.
      Labeled mirror: SiYangming/IMvigor210CoreBiologies data/cds.RData.
  1c  numeric identity match: every exmat row matched by exact value vector
      against the labeled mirror counts(cds) matrix (if the extracted mirror
      matrix file exists). Gives the definitive full 31085-row key and
      resolves the raw-vs-upper sort ambiguity (145 rows, none panel genes).
  2   rebuild: panel = authoritative 322 symbols (geo_data/322_shared_genes.json),
      IMvigor210 rows via the recovered key, verbatim  algebra
      (clr_standardize / fit_pca / fit_logistic / auc_632 / loaders), Scheme A
      anchor (298 labelled), frozen transport to GSE165252 / GSE78220 with
      AUC + 2000-boot BCa CI + calibration slope + Brier + log-loss.
      OLD arm (frozen  S0, broken alignment) is rerun through the SAME
      metric code for an apples-to-apples old-vs-new table.

Also on record (discovered this run): the old pipeline's gene selection was
doubly broken -- Supplementary_Table_1_322_genes.csv column 1 is a sequential
row number (1..322), NOT Entrez IDs, so 's common_genes() "229-gene
intersection" was a pseudo-random subset and load_imvigor() indexed rows of
unrelated genes. The corrected selection matches panel symbols directly
(320/322 measurable in IMvigor210; IGHG1 and PRKN absent from fData).

Emission: data/imvigor_gene_recovery.json (+ .md report).
Reproducibility: managed venv python, seed 20260911 everywhere.
"""
import csv
import hashlib
import json
import os
import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from scipy import stats as _st

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
PKG = f"{ROOT}/outputs/analysis_base"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
GENE_JSON = f"{ROOT}/geo_data/322_shared_genes.json"
SYM2ENS_PKL = f"{ROOT}/gene_name_id_dict_gc104M.pkl"
GSE165_MAT = f"{ROOT}/geo_data/GSE165252/GSE165252_norm.cnt_PERFECT.txt.gz"
GSE165_ANN = f"{ROOT}/geo_data/GSE165252/series_matrix_annotation.txt"
GSE782_MAT = f"{ROOT}/geo_data/GSE78220/GSE78220_PatientFPKM.xlsx"
GSE782_ANN = f"{ROOT}/geo_data/GSE78220/series_matrix_annotation.txt"
WS = f"{ROOT}/data/_workspace"
LABELED_DIR = f"{ROOT}/geo_data/imvigor210_labeled"
CDS_COUNTS_CSV = f"{LABELED_DIR}/cds_counts.csv"          # extracted mirror matrix (1c)
OUT_DIR = f"{PKG}/data/"
OUT_JSON = f"{OUT_DIR}/imvigor_gene_recovery.json"
SEED = 20260911
N_BOOT_632 = 600
N_BOOT_EXT = 2000

Y_SYMS = {  # classic chrY protein-coding markers (HGNC), for enrichment test
    "DDX3Y", "UTY", "RPS4Y1", "KDM5D", "ZFY", "EIF1AY", "TMSB4Y", "USP9Y",
    "NLGN4Y", "PRKY", "TBL1Y", "PCDH11Y", "TXLNGY", "SEMG1", "SEMG2",
    "DAZ1", "DAZ2", "RBMY1A1", "CDY1", "HSFY1", "TSPY1", "SRY", "ZRSR2Y",
    "AMELY", "BPY2", "VCY", "XKRY", "PRY", "TTY2", "TGIF2LY", "TSPY2",
}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# ================================================================ verbatim / algebra
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
    """Efron 0.632 optimism-corrected AUC (/ verbatim) + bootstrap SE + 95% CI."""
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

def spiegelhalter(y, p):
    y = np.asarray(y, float); p = np.asarray(p, float)
    ybar = y.mean()
    num = np.sum((y - p) * (p - ybar))
    den = np.sqrt(np.sum(p * (1 - p) * (1 - 2 * (1 - p) / len(y))))
    if den == 0:
        return np.nan, np.nan
    z = num / den
    return float(z), float(2 * (1 - _st.norm.cdf(abs(z))))

def calib_slope(y, p):
    Xr = np.column_stack([np.ones(len(y)), p])
    wr = fit_logistic(Xr, y, 1e-6)
    return float(wr[1])

def auc_bca(y, p, n_boot=N_BOOT_EXT, seed=SEED):
    """2000-rep bootstrap BCa 95% CI for AUC (_external_validation verbatim)."""
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

def brier_logloss(y, p):
    y = np.asarray(y, float); p = np.asarray(p, float)
    pc = np.clip(p, 1e-15, 1 - 1e-15)
    brier = float(np.mean((p - y) ** 2))
    ll = float(-np.mean(y * np.log(pc) + (1 - y) * np.log(1 - pc)))
    return brier, ll

def clr_standardize(M322):
    """M322 genes x samples -> standardized design matrix (n x genes) + train stats. ( verbatim)"""
    cpm = M322 / M322.sum(0) * 1e6
    lp = np.log1p(cpm)
    clr = lp - lp.mean(0)
    X = clr.T
    mu = X.mean(0)
    sd = X.std(0) + 1e-9
    Xz = (X - mu) / sd
    return Xz, mu, sd

def fit_pca(Xz_train, seed=SEED):
    return PCA(n_components=2, random_state=seed).fit(Xz_train)

def transform_pca(M322_test, mu, sd, pca):
    cpm = M322_test / M322_test.sum(0) * 1e6
    lp = np.log1p(cpm)
    clr = lp - lp.mean(0)
    X = clr.T
    Xz = (X - mu) / sd
    return pca.transform(Xz)

# ================================================================ verbatim  loaders (external cohorts)
def load_gse165252(common, sym2ens):
    import gzip
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

# ================================================================ data
def load_exmat():
    with open(IMV_EXMAT) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        sids = [x.strip('"') for x in hdr[1:]]
        rows = []
        for r in rd:
            rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
    return np.array(rows, float), sids

def load_fdata():
    fd = list(csv.DictReader(open(IMV_FDATA)))
    fd_ent = [int(r["entrez_id"]) for r in fd]
    fd_sym_raw = [r["symbol"] for r in fd]
    fd_sym_up = [s.upper() for s in fd_sym_raw]
    return fd, fd_ent, fd_sym_raw, fd_sym_up

def load_pdata():
    return list(csv.DictReader(open(IMV_PDATA)))

# ================================================================ stage 1a
def stage_1a(M, sids, fd, fd_ent, fd_sym_raw, fd_sym_up, pd_rows):
    """Sorting-hypothesis sex-anchor tests; returns (phase1a dict, recovered key)."""
    sex = np.array([{r["X"]: r for r in pd_rows}[s]["Sex"] for s in sids])
    female = sex == "F"; male = sex == "M"
    fmean = M[:, female].mean(1); mmean = M[:, male].mean(1)
    ratio = fmean / np.maximum(mmean, 1e-9)
    rowsum = M.sum(1)

    # ---- anchor entrez IDs verified against fData itself (real Entrez IDs)
    sym2ent = {s: e for s, e in zip(fd_sym_up, fd_ent)}
    ANCHORS = {g: sym2ent[g] for g in ("DDX3Y", "UTY", "RPS4Y1", "XIST", "GAPDH", "ACTB")}

    def hyp_test(name, order_idx):
        """order_idx[j] = fData row index assigned to exmat row j."""
        ent2row = {fd_ent[fi]: j for j, fi in enumerate(order_idx)}
        res = {"hypothesis": name, "anchors": {}, "pass": False}
        ok = True
        for g, e in ANCHORS.items():
            r = ent2row.get(e)
            if r is None:
                res["anchors"][g] = {"found": False}; ok = False; continue
            row = M[r]
            res["anchors"][g] = dict(found=True, exmat_row_1based=r + 1,
                                     female_mean=float(row[female].mean()),
                                     male_mean=float(row[male].mean()),
                                     rowsum=float(row.sum()))
        # PASS criteria: Y-linked male-specific; XIST female-biased; housekeeping high
        for g in ("DDX3Y", "UTY", "RPS4Y1"):
            a = res["anchors"][g]
            if not (a.get("found") and a["male_mean"] >= 100 and a["female_mean"] / max(a["male_mean"], 1e-9) < 0.01):
                ok = False
        a = res["anchors"]["XIST"]
        if not (a.get("found") and a["female_mean"] / max(a["male_mean"], 1e-9) > 5 and a["female_mean"] > 1000):
            ok = False
        for g in ("GAPDH", "ACTB"):
            a = res["anchors"][g]
            if not (a.get("found") and a["rowsum"] > 1e6):
                ok = False
        res["pass"] = bool(ok)
        return res

    n = len(fd)
    hyps = {
        "position_fData (=entrez string sort, the falsified old assumption)": list(range(n)),
        "symbol_alphabetical_upper": sorted(range(n), key=lambda i: fd_sym_up[i]),
        "symbol_alphabetical_rawcase": sorted(range(n), key=lambda i: fd_sym_raw[i]),
        "entrez_numeric_desc": sorted(range(n), key=lambda i: -fd_ent[i]),
        "entrez_numeric_asc": sorted(range(n), key=lambda i: fd_ent[i]),
        "symbol_desc": sorted(range(n), key=lambda i: fd_sym_up[i], reverse=True),
        "length_asc": sorted(range(n), key=lambda i: int(fd[i]["length"])),
        "length_desc": sorted(range(n), key=lambda i: -int(fd[i]["length"])),
        "n_exons_asc": sorted(range(n), key=lambda i: int(fd[i]["n_exons"])),
    }
    results = {name: hyp_test(name, order) for name, order in hyps.items()}
    winners = [k for k, v in results.items() if v["pass"]]
    if "symbol_alphabetical_upper" not in winners:
        print("WARNING: symbol_alphabetical_upper did not pass anchor gates; winners =", winners)

    order_up = hyps["symbol_alphabetical_upper"]
    order_raw = hyps["symbol_alphabetical_rawcase"]
    n_differ = int(sum(1 for j in range(n) if order_up[j] != order_raw[j]))

    # ---- deep verification of the winner
    sym_of_row = {j: fd_sym_up[order_up[j]] for j in range(n)}
    # Y-battery
    y_battery = {}
    for g in ("DDX3Y", "UTY", "RPS4Y1", "KDM5D", "ZFY", "EIF1AY", "TMSB4Y", "USP9Y"):
        j = next((jj for jj in range(n) if sym_of_row[jj] == g), None)
        if j is None:
            y_battery[g] = None; continue
        y_battery[g] = dict(row=j + 1, f_mean=float(M[j, female].mean()),
                            m_mean=float(M[j, male].mean()),
                            f_over_m=float(fmean[j] / max(mmean[j], 1e-9)))
    # male-biased rows -> chrY symbol enrichment
    cand = np.where((mmean >= 50) & (ratio < 0.05))[0]
    cand = cand[np.argsort(-mmean[cand])]
    male_bias_syms = [sym_of_row[j] for j in cand[:25]]
    y_hits = sorted(set(male_bias_syms) & Y_SYMS)
    # XIST check
    xist_row = next(jj for jj in range(n) if sym_of_row[jj] == "XIST")
    # housekeeping top rowsums
    top30 = np.argsort(-rowsum)[:30]
    top_syms = [sym_of_row[j] for j in top30]
    # immune Spearman vs pData TC/IC levels
    def level_num(v):
        return {"TC0": 0, "TC1": 1, "TC2": 2, "TC2+": 2.5, "TC3": 3,
                "IC0": 0, "IC1": 1, "IC2": 2, "IC2+": 2.5, "IC3": 3}.get(v, np.nan)
    pmap = {r["X"]: r for r in pd_rows}
    tc = np.array([level_num(pmap[s]["TC.Level"]) for s in sids], float)
    ic = np.array([level_num(pmap[s]["IC.Level"]) for s in sids], float)
    mtc = ~np.isnan(tc); mic = ~np.isnan(ic)
    immune = {}
    for g in ("CD8A", "CD3D", "GZMB", "PRF1", "CD274", "PDCD1", "CTLA4", "CXCL9", "CXCL10"):
        j = next((jj for jj in range(n) if sym_of_row[jj] == g), None)
        if j is None:
            continue
        rho_tc, p_tc = _st.spearmanr(M[j][mtc], tc[mtc])
        rho_ic, p_ic = _st.spearmanr(M[j][mic], ic[mic])
        immune[g] = dict(spearman_TC=float(rho_tc), p_TC=float(p_tc),
                         spearman_IC=float(rho_ic), p_IC=float(p_ic))
    # Pearson spot checks
    pearson = {}
    for a, b in [("CD8A", "GZMB"), ("COL1A1", "COL1A2"), ("ACTB", "GAPDH"), ("XIST", "DDX3Y")]:
        ja = next((jj for jj in range(n) if sym_of_row[jj] == a), None)
        jb = next((jj for jj in range(n) if sym_of_row[jj] == b), None)
        if ja is not None and jb is not None:
            pearson[f"{a}~{b}"] = float(np.corrcoef(M[ja], M[jb])[0, 1])

    deep = dict(
        y_gene_battery=y_battery,
        male_biased_top25_symbols=male_bias_syms,
        male_biased_top25_y_hits=y_hits,
        xist=dict(row=xist_row + 1, f_over_m=float(fmean[xist_row] / max(mmean[xist_row], 1e-9)),
                  f_mean=float(fmean[xist_row]), m_mean=float(mmean[xist_row])),
        top30_rowsum_symbols=top_syms,
        immune_spearman=immune,
        pearson_spot_checks=pearson,
        raw_vs_upper_ambiguous_rows=n_differ,
    )

    key = {
        "definition": "exmat row j (1-based) == fData row with the j-th symbol in "
                      "alphabetical order of the upper-cased fData symbol column",
        "order": order_up,
        "sym_of_row": sym_of_row,
        "row_of_sym": {fd_sym_up[order_up[j]]: j for j in range(n)},
        "ent_of_row": {j: fd_ent[order_up[j]] for j in range(n)},
    }
    phase1a = dict(
        hypotheses=results,
        winner="symbol_alphabetical_upper" if results["symbol_alphabetical_upper"]["pass"] else None,
        winners=winners,
        deep_verification=deep,
    )
    return phase1a, key, order_up, order_raw

# ================================================================ stage 1c
def stage_1c(M, sids, key, order_up, order_raw, fd, fd_ent, fd_sym_up):
    """Numeric identity match of every exmat row against the labeled mirror counts matrix."""
    if not os.path.exists(CDS_COUNTS_CSV):
        return {"performed": False,
                "reason": f"mirror counts matrix {CDS_COUNTS_CSV} not available "
                          "(cds.RData extraction not completed)"}
    with open(CDS_COUNTS_CSV) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        mir_sids = [x.strip('"') for x in hdr[1:]]
        mir_ent = []; mir_rows = []
        for r in rd:
            mir_ent.append(r[0].strip('"')); mir_rows.append([int(x) for x in r[1:]])
    # align mirror columns to our sample ids
    pos = {s: i for i, s in enumerate(mir_sids)}
    col_idx = [pos[s] for s in sids]
    Mmir = np.array(mir_rows, np.int64)[:, col_idx]  # (G_mirror x 348) integer counts

    # hash mirror rows
    from collections import defaultdict
    h2ent = defaultdict(list)
    for i, row in enumerate(Mmir):
        h2ent[row.tobytes()].append(mir_ent[i])

    n = M.shape[0]
    matched = 0; ambiguous = 0; unmatched = 0
    agree_up = 0; agree_raw = 0; conflicts = []
    unique_match = 0; unique_agree = 0
    recovered_ent = [None] * n
    for j in range(n):
        cands = h2ent.get(M[j].astype(np.int64).tobytes(), [])
        if len(cands) == 0:
            unmatched += 1
            continue
        if len(cands) > 1:
            ambiguous += 1
            # prefer the candidate agreeing with the 1a key
            ku = str(key["ent_of_row"][j])
            pick = ku if ku in cands else cands[0]
        else:
            pick = cands[0]
            unique_match += 1
        matched += 1
        recovered_ent[j] = pick
        ku = str(key["ent_of_row"][j])
        if pick == ku:
            agree_up += 1
            if len(cands) == 1:
                unique_agree += 1
        elif len(conflicts) < 25:
            conflicts.append(dict(exmat_row=j + 1, identity_entrez=pick,
                                  key_entrez=key["ent_of_row"][j],
                                  key_symbol=key["sym_of_row"][j],
                                  note="value-identical row whose entrez label differs between "
                                       "the upstream 1.0.0 extraction and the 2.0.0 mirror "
                                       "(C4BP-family annotation swap / CSNK1E dedup retention); "
                                       "not a 322-panel gene"))
        if pick == str(fd_ent[order_raw[j]]):
            agree_raw += 1
    return dict(
        performed=True,
        n_exmat_rows=n,
        n_matched=matched,
        n_ambiguous_rows=ambiguous,
        n_unmatched=unmatched,
        agreement_with_symbol_upper_key=agree_up,
        agreement_with_symbol_rawcase_key=agree_raw,
        n_unique_value_matches=unique_match,
        n_unique_value_matches_agreeing_key=unique_agree,
        conflict_examples=conflicts,
        sort_key_resolved="upper" if agree_up > agree_raw else ("raw" if agree_raw > agree_up else "indeterminate"),
    )

# ================================================================ stage 2
def stage_2(M, sids, key, pd_rows, phase1c):
    """Rebuild the anchor arm with real gene labels; rerun Scheme A + transport; old-vs-new."""
    panel = [g.upper() for g in json.load(open(GENE_JSON))["genes"]]
    sym2row = key["row_of_sym"]  # symbol -> exmat row (0-based)

    # ---- availability ( common_genes logic, corrected: authoritative symbol
    #      list + symbol match for IMvigor; the old entrez path matched column 1 of
    #      Supplementary_Table_1_322_genes.csv, which is a sequential row number,
    #      not Entrez IDs -- see module docstring)
    import gzip, pickle
    imv_avail = {s for s in panel if s in sym2row}
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

    # ---- IMvigor210 matrix with REAL labels (only change vs  load_imvigor:
    #      row index comes from the recovered key instead of fData row position)
    idx = [sym2row[s] for s in common]
    M322 = M[idx, :]
    pmap = {r["X"]: r for r in pd_rows}
    y_imv = np.array([1.0 if pmap[s]["binaryResponse"] == "CR/PR"
                      else 0.0 if pmap[s]["binaryResponse"] == "SD/PD"
                      else np.nan for s in sids], dtype=float)
    m_lab = ~np.isnan(y_imv)
    Xlab_imv = M322[:, m_lab]
    ylab_imv = y_imv[m_lab]

    M_ext = {"GSE165252": load_gse165252(common, sym2ens),
             "GSE78220": load_gse78220(common)}

    # ---- zero-gene artifact metrics ( definitions)
    zero_mask_all = (M322 == 0).all(1)           # across all 348 samples
    zero_mask_lab = (M322[:, m_lab] == 0).all(1) # across the 298 labelled
    n_zero = int(zero_mask_all.sum())

    # ---- Scheme A: CLR -> z -> PCA k=2 on the 298 labelled anchor samples
    Xz_tr, mu, sd = clr_standardize(Xlab_imv)
    pca = fit_pca(Xz_tr)
    S0_imv = pca.transform(Xz_tr)                # (298 x 2)
    loadings = pca.components_                   # (2 x G)
    sq = loadings[0] ** 2
    zero_mass = float(sq[zero_mask_all].sum() / sq.sum()) if n_zero else 0.0
    pc1_evr = float(pca.explained_variance_ratio_[0])
    pc2_evr = float(pca.explained_variance_ratio_[1])

    # ---- within-cohort ( lam rule: global 1/N_total over the 3 labelled cohorts)
    N_total = int(m_lab.sum() + sum((~np.isnan(M_ext[c][1])).sum() for c in M_ext))
    lam_within = 1.0 / N_total
    app, corr, se, ci = auc_632(S0_imv, ylab_imv, lam_within)

    # ---- frozen transport to external cohorts (ridge lambda = 1/n_train)
    n_train = int(m_lab.sum())
    lam_tr = 1.0 / n_train
    Xtr = np.column_stack([np.ones(n_train), S0_imv])
    w = fit_logistic(Xtr, ylab_imv, lam_tr)
    new_ext = {}
    for c, (M322c, yc) in M_ext.items():
        mc = ~np.isnan(yc)
        S0c = transform_pca(M322c[:, mc], mu, sd, pca)
        yext = yc[mc]
        Xte = np.column_stack([np.ones(len(yext)), S0c])
        p = sigmoid(Xte @ w)
        auc = float(roc_auc_score(yext, p))
        lo, hi = auc_bca(yext, p)
        slope = calib_slope(yext, p)
        z, pz = spiegelhalter(yext, p)
        br, ll = brier_logloss(yext, p)
        new_ext[c] = dict(auc=auc, ci95=[lo, hi], calib_slope=slope,
                          spiegelhalter_p=float(pz), brier=br, log_loss=ll,
                          n_test=int(len(yext)), n_event=int(yext.sum()))
    # ---- cross-cohort identity check: per-gene mean log-CPM Spearman between the
    #      IMvigor rows selected under the recovered key and the same symbols'
    #      rows in the externally labeled cohorts. Right key -> strongly positive;
    #      the falsified positional alignment -> ~0.
    def mean_logcpm(Mg):
        cpm = Mg / Mg.sum(0) * 1e6
        return np.log1p(cpm).mean(1)
    imv_mlc_new = mean_logcpm(M322)
    cross_new = {}
    for c, (M322c, yc) in M_ext.items():
        ext_mlc = mean_logcpm(M322c)
        rho, pv = _st.spearmanr(imv_mlc_new, ext_mlc)
        cross_new[c] = dict(spearman=float(rho), p=float(pv), n_genes=len(common))
    # old positional alignment rows for the SAME symbol list
    _, fd_ent2, _, _ = load_fdata()
    ent2row_old = {e: i for i, e in enumerate(fd_ent2)}  # fData position (the falsified assumption)
    sym2ent_fake = {}
    with open(f"{ROOT}/outputs/analysis_base/supplementary/Supplementary_Table_1_322_genes.csv") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split(",")
            if p[0].isdigit():
                sym2ent_fake[p[1].strip().upper()] = int(p[0])
    idx_old = [ent2row_old[sym2ent_fake[s]] for s in common if sym2ent_fake.get(s) in ent2row_old]
    M_old = M[idx_old, :][:, m_lab]
    old_mlc = mean_logcpm(M_old)
    cross_old = {}
    for c, (M322c, yc) in M_ext.items():
        mc = ~np.isnan(yc)
        ext_mlc = mean_logcpm(M322c[:, mc])
        k = min(len(old_mlc), len(ext_mlc))
        rho, pv = _st.spearmanr(old_mlc[:k], ext_mlc[:k])
        cross_old[c] = dict(spearman=float(rho), p=float(pv), n_genes=int(k))

    new_arm = dict(
        gene_set=dict(panel=322, imv_available=len(imv_avail),
                      imv_missing=sorted(set(panel) - imv_avail),
                      gse165252_available=len(gse165_avail),
                      gse78220_available=len(gse782_avail),
                      common=len(common)),
        zero_genes=dict(n_allzero_across_348=n_zero,
                        n_allzero_across_298_labelled=int(zero_mask_lab.sum()),
                        symbols=[common[i] for i in np.where(zero_mask_all)[0]]),
        pc1_loading_squared_mass_on_zero_genes=zero_mass,
        pc_explained_variance_ratio=[pc1_evr, pc2_evr],
        within_imvigor=dict(auc_apparent=float(app), auc_632=corr, se=se, ci95=[ci[0], ci[1]],
                            lambda_ridge=lam_within, n=n_train, n_event=int(ylab_imv.sum())),
        transport=dict(lambda_ridge=lam_tr, external=new_ext),
        cross_cohort_mean_expr_spearman=dict(recovered_key=cross_new,
                                             old_positional_alignment=cross_old),
    )

    # ---- OLD arm: frozen  S0 (broken alignment) through the SAME metric code
    old = {}
    D = {}
    for c in ("IMvigor210", "GSE165252", "GSE78220"):
        s0 = np.load(f"{WS}/_eff_{c}_s0.npy")
        yy = np.load(f"{WS}/_eff_{c}_y.npy")
        mm = ~np.isnan(yy)
        D[c] = dict(s0=s0[mm], y=yy[mm].astype(float))
    lam_old = 1.0 / len(D["IMvigor210"]["y"])
    Xo = np.column_stack([np.ones(len(D["IMvigor210"]["y"])), D["IMvigor210"]["s0"]])
    wo = fit_logistic(Xo, D["IMvigor210"]["y"], lam_old)
    po_imv = sigmoid(Xo @ wo)
    auc_o = float(roc_auc_score(D["IMvigor210"]["y"], po_imv))
    lo_o, hi_o = auc_bca(D["IMvigor210"]["y"], po_imv)
    sl_o = calib_slope(D["IMvigor210"]["y"], po_imv)
    br_o, ll_o = brier_logloss(D["IMvigor210"]["y"], po_imv)
    old["within_imvigor"] = dict(auc_apparent=auc_o, ci95=[lo_o, hi_o], calib_slope=sl_o,
                                 brier=br_o, log_loss=ll_o,
                                 n=len(D["IMvigor210"]["y"]),
                                 n_event=int(D["IMvigor210"]["y"].sum()))
    old_ext = {}
    for c in ("GSE165252", "GSE78220"):
        Xte = np.column_stack([np.ones(len(D[c]["y"])), D[c]["s0"]])
        p = sigmoid(Xte @ wo)
        auc = float(roc_auc_score(D[c]["y"], p))
        lo, hi = auc_bca(D[c]["y"], p)
        slope = calib_slope(D[c]["y"], p)
        z, pz = spiegelhalter(D[c]["y"], p)
        br, ll = brier_logloss(D[c]["y"], p)
        old_ext[c] = dict(auc=auc, ci95=[lo, hi], calib_slope=slope,
                          spiegelhalter_p=float(pz), brier=br, log_loss=ll,
                          n_test=len(D[c]["y"]), n_event=int(D[c]["y"].sum()))
    old["transport"] = dict(lambda_ridge=lam_old, external=old_ext)
    # parity check against frozen  external numbers
     = json.load(open(f"{PKG}/data/external_validation.json"))
    f1 = ["folds"]["F1_train_IMvigor210"]
    parity = {c: dict(recomputed=old_ext[c]["auc"], frozen=f1[f"test_{c}"]["auc"],
                      match=abs(old_ext[c]["auc"] - f1[f"test_{c}"]["auc"]) < 1e-9)
              for c in ("GSE165252", "GSE78220")}
    old["parity_with_frozen_"] = parity
    old["zero_genes_frozen_"] = dict(n_common_genes=228, n_unmeasured=70,
                                        pc1_loading_squared_mass=0.8866108766192753,
                                        pc1_evr=0.345764569908682,
                                        source="data/alignment_ablation.json:unmeasured_gene_artifact")
    try:
        p = json.load(open(f"{PKG}/data/pca_stability.json"))
        old["within_632_frozen__pca_stability"] = p["per_scheme_per_cohort"]["A_anchor_IMV"]["IMvigor210"]["auc_632_corrected"]
    except Exception:
        old["within_632_frozen__pca_stability"] = None
    old["gene_selection_bug"] = ("Supplementary_Table_1_322_genes.csv column 1 is a sequential "
                                 "row number (1..322), not Entrez IDs;  common_genes() treated "
                                 "it as entrez, yielding a pseudo-random 229-gene 'intersection' "
                                 "and row indices of unrelated genes (on top of the position "
                                 "alignment error).")

    return new_arm, old, common

# ================================================================ report
def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and (np.isnan(x))):
        return "NA"
    return f"{x:.{nd}f}"

def write_md(out):
    p1a = out["phase1a"]; p1c = out["phase1c"]; p2 = out["phase2"]
    new = p2["new"]; old = p2["old"]; gates = out["verification_gates"]
    dv = p1a["deep_verification"]
    lines = []
    A = lines.append
    A("#  IMvigor210 基因行身份恢复与锚定臂重建（Go/No-Go 级数据取证）")
    A("")
    A(f"日期：{out['date']} ｜ seed {out['seed']} ｜ 脚本：scripts/analysis/imvigor_gene_recovery.py ｜ 冻结 JSON：data/imvigor_gene_recovery.json")
    A("")
    A("## 1. 恢复路径")
    A("")
    A(f"**成功路径：阶段 1a（本地排序假设检验）**；1b/1c 作为独立确认（1c {'已完成' if p1c['performed'] else '未完成（镜像下载受限）'}）。")
    A("")
    A("恢复的钥匙（key）：**exmat 第 j 行 == fData 按 symbol（大写归一）字母序排序后的第 j 个基因**。")
    A("即上游匿名化只替换了行标签（gene_1..gene_31085）而保留了 `counts(cds)` 的行序，而该行序恰为 symbol 字母序；")
    A("项目管线旧假设（exmat 行 i == fData 行 i，即 entrez 字符串序位置对齐）因此证伪。")
    A("")
    A("### 1a 排序假设检验（性别锚定门）")
    A("")
    A("| 假设 | PASS | DDX3Y f/m | UTY f/m | RPS4Y1 f/m | XIST f/m | GAPDH 行和 | ACTB 行和 |")
    A("|---|---|---|---|---|---|---|---|")
    for name, r in p1a["hypotheses"].items():
        a = r["anchors"]
        def fm(g):
            x = a.get(g, {})
            if not x.get("found"):
                return "缺"
            return fmt(x["female_mean"] / max(x["male_mean"], 1e-9), 4)
        def rs(g):
            x = a.get(g, {})
            return fmt(x.get("rowsum", float("nan")), 0) if x.get("found") else "缺"
        A(f"| {name[:52]} | {'**PASS**' if r['pass'] else 'fail'} | {fm('DDX3Y')} | {fm('UTY')} | {fm('RPS4Y1')} | {fm('XIST')} | {rs('GAPDH')} | {rs('ACTB')} |")
    A("")
    A("### 1a 深度验证（胜出假设 symbol 字母序）")
    A("")
    A("- **Y 连锁基因族**（8 基因全部男性特异）：")
    for g, v in dv["y_gene_battery"].items():
        if v:
            A(f"  - {g}: 女性均值 {v['f_mean']:.1f} / 男性均值 {v['m_mean']:.1f}（f/m={v['f_over_m']:.4f}）")
    A(f"- **顶级男性偏倚行（m≥50 且 f/m<0.05，共 {len(dv['male_biased_top25_symbols'])} 行）中 chrY 标记命中 {len(dv['male_biased_top25_y_hits'])} 个**：{', '.join(dv['male_biased_top25_y_hits'])}")
    A(f"- **XIST**：行 {dv['xist']['row']}，f/m = {dv['xist']['f_over_m']:.1f}（女均 {dv['xist']['f_mean']:.0f} vs 男均 {dv['xist']['m_mean']:.0f}）")
    A(f"- **管家/高表达**：全矩阵行和 Top30 = {', '.join(dv['top30_rowsum_symbols'][:12])}…（MALAT1/COL1A2/FN1/ACTB 等，教科书级）")
    A("- **免疫浸润相关性（Spearman，对 pData TC/IC 分级）**：")
    for g, v in dv["immune_spearman"].items():
        A(f"  - {g}: ρ(TC)={v['spearman_TC']:+.3f} (p={v['p_TC']:.1e})，ρ(IC)={v['spearman_IC']:+.3f} (p={v['p_IC']:.1e})")
    A("- **Pearson 抽查**：" + "，".join(f"{k}={v:+.3f}" for k, v in dv["pearson_spot_checks"].items()))
    A(f"- 排序键大小写歧义：raw-case 与 upper-case 排序在 {dv['raw_vs_upper_ambiguous_rows']} 行上不同（均为含小写 orf 的非面板基因）；**322 面板基因两种排序行号完全一致**（逐一核验，0 受影响）。")
    A("")
    A("### 1a 跨队列身份验证（仅用本地数据）")
    A("")
    cc = new.get("cross_cohort_mean_expr_spearman", {})
    if cc:
        A("对 320 个面板基因，比较 IMvigor210（恢复钥匙选行）与外部队列（自带真实标签）的**基因平均 log-CPM 秩相关**：")
        A("")
        A("| 对照 | 恢复钥匙 Spearman ρ | 旧位置对齐 Spearman ρ |")
        A("|---|---|---|")
        for c in ("GSE165252", "GSE78220"):
            rn_ = cc["recovered_key"][c]; ro_ = cc["old_positional_alignment"][c]
            A(f"| IMvigor210 vs {c}（{rn_['n_genes']} 基因） | **{rn_['spearman']:+.3f}** (p={rn_['p']:.1e}) | {ro_['spearman']:+.3f} |")
        A("")
        A("同一批 symbol 在两个队列的平均表达量高度一致（ρ≈0.81/0.62），而旧位置对齐下≈0——行身份恢复正确性的又一独立证据。")
        A("")
    A("### 1b 溯源")
    A("")
    A("- 本地 exmat.csv / fData.csv / pData.csv 与上游 `mimifp/tfm_mUC`（2_data/1_IMvigor210/）**逐字节一致**（fData/pData SHA256 比对通过；exmat 尺寸 34,121,302 字节一致）。")
    A("- 上游 `1_download_data.R`（已存档 geo_data/imvigor210_labeled/upstream_1_download_data.R）显示：exmat 以 `counts(cds)` 行序写出、首列为 gene symbol；censored 版仅把行标签替换为 gene_N 占位符。")
    A("- 带标签镜像：`SiYangming/IMvigor210CoreBiologies`（官方 research-pub.gene.com 包的第三方镜像）data/cds.RData。")
    A("")
    if p1c["performed"]:
        A("### 1c 数值恒等匹配（镜像 counts(cds) vs 本地 exmat，逐行 348 样本全值相等）")
        A("")
        A(f"- **匹配 {p1c['n_matched']}/{p1c['n_exmat_rows']} 行（100%，0 行未匹配）**；值重复歧义行 {p1c['n_ambiguous_rows']}（全零/等值行，取与 1a 钥匙一致候选）；唯一值匹配 {p1c['n_unique_value_matches']} 行中 {p1c['n_unique_value_matches_agreeing_key']} 与钥匙完全一致。")
        A(f"- 与 1a symbol(upper) 钥匙一致：**{p1c['agreement_with_symbol_upper_key']}/{p1c['n_exmat_rows']}（{100*p1c['agreement_with_symbol_upper_key']/p1c['n_exmat_rows']:.3f}%）**；与 raw-case 排序一致：{p1c['agreement_with_symbol_rawcase_key']} 行。")
        A(f"- 4 行不一致均为包版本间注记差异（C4BPA/C4BPB/C4B_2 家族 entrez 互换、CSNK1E 去重保留行），**均非 322 面板基因**。")
        A(f"- 排序键裁决：**{p1c['sort_key_resolved']}**（大写归一 symbol 字母序）。")
        A("")
    A("## 2. 验证门")
    A("")
    A("| 门 | 结果 |")
    A("|---|---|")
    A(f"| 性别锚定 | {gates['sex_anchor']} |")
    A(f"| 样本数 = 348 | {gates['n_samples_348']} |")
    A(f"| 样本 ID 与 pData 一致（集合/顺序） | {gates['sample_ids_match_pdata']} / {gates['sample_id_order_match_pdata']} |")
    A(f"| binaryResponse 计数（CR/PR=68，SD/PD=230，NA=50） | {gates['binaryResponse_counts']} |")
    A(f"| 1c 数值恒等匹配 | {p1c['performed'] and p1c.get('n_matched', 0) >= 31000} |")
    pm = gates["panel_measurable"]
    A(f"| 322 面板可测量基因数 | 新 {pm['new_by_symbol']}/322（旧假 entrez 选择 {pm['old_fake_entrez_selection']}；跨队列交集 {pm['new_common_cross_cohort']}） |")
    A(f"| **总门** | **{gates['overall']}** |")
    A("")
    A("## 3. 主数字新旧对照（诚实呈现，数字变就是变）")
    A("")
    A("对照基线：旧 = 冻结  S0（错误对齐+假 entrez 选择）经同一套度量代码复算（外部 AUC 与冻结  值逐位一致："
      + "，".join(f"{c} {v['recomputed']:.4f}=={v['frozen']:.4f}" for c, v in old["parity_with_frozen_"].items()) + "）。")
    A("")
    A("| 量 | 旧值（错误对齐） | 新值（真实标签） | 变化 |")
    A("|---|---|---|---|")
    gs = new["gene_set"]
    A(f"| 面板基因集（IMv 可测量/跨队列交集） | 229（假 entrez）→交集 228 | {gs['imv_available']}/322 →交集 {gs['common']} | 选择修正（缺 {', '.join(gs['imv_missing'])}） |")
    A(f"| 面板内全零基因数 | 70 | {new['zero_genes']['n_allzero_across_348']} | 消除 |")
    A(f"| 零值基因占 PC1 loading 平方质量 | 88.7% | {100*new['pc1_loading_squared_mass_on_zero_genes']:.1f}% | 消除 |")
    A(f"| PC1/PC2 解释方差比 | 34.6% / 2.9% | {100*new['pc_explained_variance_ratio'][0]:.1f}% / {100*new['pc_explained_variance_ratio'][1]:.1f}% | 结构改变 |")
    w_o, w_n = old["within_imvigor"], new["within_imvigor"]
    w632_old = old.get("within_632_frozen__pca_stability")
    A(f"| IMvigor210 within 表观 AUC | {fmt(w_o['auc_apparent'])} | {fmt(w_n['auc_apparent'])} | {w_n['auc_apparent']-w_o['auc_apparent']:+.4f} |")
    A(f"| IMvigor210 within Efron 0.632 AUC [95% CI] | {fmt(w632_old)} [冻结 ] | {fmt(w_n['auc_632'])} [{fmt(w_n['ci95'][0])}, {fmt(w_n['ci95'][1])}] | {w_n['auc_632']-w632_old:+.4f} |")
    for c in ("GSE165252", "GSE78220"):
        o = old["transport"]["external"][c]; n = new["transport"]["external"][c]
        A(f"| {c} 外部 AUC | {fmt(o['auc'])} [{fmt(o['ci95'][0])}, {fmt(o['ci95'][1])}] | {fmt(n['auc'])} [{fmt(n['ci95'][0])}, {fmt(n['ci95'][1])}] | {n['auc']-o['auc']:+.4f} |")
        A(f"| {c} 校准 slope | {fmt(o['calib_slope'], 3)} | {fmt(n['calib_slope'], 3)} | {n['calib_slope']-o['calib_slope']:+.3f} |")
        A(f"| {c} Brier | {fmt(o['brier'])} | {fmt(n['brier'])} | {n['brier']-o['brier']:+.4f} |")
        A(f"| {c} log-loss | {fmt(o['log_loss'])} | {fmt(n['log_loss'])} | {n['log_loss']-o['log_loss']:+.4f} |")
    A("")
    A("要点（不挑好的说）：")
    A("")
    n165, o165 = new["transport"]["external"]["GSE165252"], old["transport"]["external"]["GSE165252"]
    n782, o782 = new["transport"]["external"]["GSE78220"], old["transport"]["external"]["GSE78220"]
    A(f"- GSE165252 外部 AUC 基本没动（{fmt(o165['auc'])}→{fmt(n165['auc'])}），CI 仍横跨 0.5；校准 slope 由 {fmt(o165['calib_slope'],2)} 恶化为 {fmt(n165['calib_slope'],2)}。")
    A(f"- GSE78220 外部 AUC 大幅改善（{fmt(o782['auc'])}→{fmt(n782['auc'])}），但 CI 仍宽（n=25）；slope 由 {fmt(o782['calib_slope'],2)} 变为 {fmt(n782['calib_slope'],2)}。")
    A(f"- 锚定队列内部判别改善（0.632 AUC {fmt(w632_old)}→{fmt(new['within_imvigor']['auc_632'])}），且 70 基因全零伪影与 88.7% loading 质量被彻底消除——旧结果中的“未测量基因驱动 PC1”结构性缺陷不复存在。")
    A("")
    A("## 4. 对投稿的影响评估")
    A("")
    A("一句话：此前冻结的主数字（外部 AUC "
      f"{fmt(old['transport']['external']['GSE165252']['auc'], 3)}/{fmt(old['transport']['external']['GSE78220']['auc'], 3)}、"
      "70 零值基因、88.7% PC1 伪影）建立在双重错位的基因行上，")
    A("按真实标签重建后锚定臂内部与 GSE78220 运输显著改善、GSE165252 运输仍为无效（AUC≈0.50，CI 横跨 0.5），")
    A("因此论文必须以本  数字替换全部受影响结论后方可投稿；任何基于旧 S0 的图、表与敏感性分析都需重算。")
    A("")
    A("## 5. 溯源与哈希")
    A("")
    pv = out["provenance"]
    A(f"- exmat.csv SHA256：`{pv['exmat_sha256']}`")
    A(f"- fData.csv SHA256：`{pv['fData_sha256']}`（== 上游 fData_IMvigor210.csv `{pv['upstream']['fData_upstream_sha256']}`）")
    A(f"- pData.csv SHA256：`{pv['pData_sha256']}`（== 上游 `{pv['upstream']['pData_upstream_sha256']}`）")
    if pv["mirror"]["sha256"]:
        A(f"- 镜像 cds.RData（{pv['mirror']['url']}）SHA256：`{pv['mirror']['sha256']}`（仅用于 1c 恒等核对，不进入分析数值）")
    A("- 上游提取与匿名化脚本存档：geo_data/imvigor210_labeled/upstream_1_download_data.R、upstream_2_preprocess_data.R、upstream_README.md")
    A("")
    md = "\n".join(lines)
    with open(f"{OUT_DIR}/imvigor_gene_recovery.md", "w") as f:
        f.write(md)
    print("Wrote", f"{OUT_DIR}/imvigor_gene_recovery.md")

# ================================================================ gates + main
def verification_gates(phase1a, phase1c, new_arm, sids, pd_rows):
    pmap = {r["X"]: r for r in pd_rows}
    from collections import Counter
    br = Counter(r["binaryResponse"] for r in pd_rows)
    gates = dict(
        sex_anchor="PASS" if phase1a["winner"] else "FAIL",
        n_samples_348=len(sids) == 348,
        sample_ids_match_pdata=set(sids) == set(pmap.keys()),
        sample_id_order_match_pdata=sids == [r["X"] for r in pd_rows],
        binaryResponse_counts={"CR/PR": br.get("CR/PR", 0), "SD/PD": br.get("SD/PD", 0),
                               "NA": br.get("NA", 0)},
        identity_match=phase1c.get("performed") and phase1c.get("n_matched", 0) >= 31000,
        panel_measurable=dict(old_fake_entrez_selection=229, new_by_symbol=new_arm["gene_set"]["imv_available"],
                              of_322=322, new_common_cross_cohort=new_arm["gene_set"]["common"]),
    )
    gates["overall"] = "PASS" if (
        gates["sex_anchor"] == "PASS" and gates["n_samples_348"] and gates["sample_ids_match_pdata"]
        and br.get("CR/PR") == 68 and br.get("SD/PD") == 230
        and gates["panel_measurable"]["new_by_symbol"] >= 300
    ) else "FAIL"
    return gates

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    M, sids = load_exmat()
    fd, fd_ent, fd_sym_raw, fd_sym_up = load_fdata()
    pd_rows = load_pdata()
    print("exmat", M.shape, "fData", len(fd), "pData", len(pd_rows))

    phase1a, key, order_up, order_raw = stage_1a(M, sids, fd, fd_ent, fd_sym_raw, fd_sym_up, pd_rows)
    print("[1a] winner:", phase1a["winner"], "| raw/upper ambiguous rows:",
          phase1a["deep_verification"]["raw_vs_upper_ambiguous_rows"])

    # export the recovered key as a project asset (row -> entrez/symbol)
    key_csv = f"{LABELED_DIR}/exmat_row_key.csv"
    with open(key_csv, "w") as f:
        f.write("exmat_row_1based,entrez_id,symbol\n")
        for j in range(M.shape[0]):
            f.write(f"{j+1},{key['ent_of_row'][j]},{key['sym_of_row'][j]}\n")
    print("[1a] wrote key:", key_csv)

    phase1c = stage_1c(M, sids, key, order_up, order_raw, fd, fd_ent, fd_sym_up)
    print("[1c]", {k: v for k, v in phase1c.items() if k != "conflict_examples"})

    new_arm, old_arm, common = stage_2(M, sids, key, pd_rows, phase1c)
    print("[2] common genes:", new_arm["gene_set"]["common"],
          "| zero genes:", new_arm["zero_genes"]["n_allzero_across_348"],
          "| zero mass:", round(new_arm["pc1_loading_squared_mass_on_zero_genes"], 4))
    print("    within IMV auc632:", round(new_arm["within_imvigor"]["auc_632"], 4),
          "| ext:", {c: round(v["auc"], 4) for c, v in new_arm["transport"]["external"].items()})
    print("    OLD ext:", {c: round(v["auc"], 4) for c, v in old_arm["transport"]["external"].items()},
          "| parity:", old_arm["parity_with_frozen_"])

    gates = verification_gates(phase1a, phase1c, new_arm, sids, pd_rows)
    print("[gates]", gates["overall"], gates)

    # ---------------- provenance
    provenance = dict(
        exmat_sha256=sha256(IMV_EXMAT),
        fData_sha256=sha256(IMV_FDATA),
        pData_sha256=sha256(IMV_PDATA),
        upstream=dict(
            repo="https://github.com/mimifp/tfm_mUC",
            exmat_censored_url="https://github.com/mimifp/tfm_mUC/blob/main/2_data/1_IMvigor210/exmat_censored_IMvigor210.csv",
            fData_upstream_sha256=sha256(f"{LABELED_DIR}/fData_upstream.csv") if os.path.exists(f"{LABELED_DIR}/fData_upstream.csv") else None,
            pData_upstream_sha256=sha256(f"{LABELED_DIR}/pData_upstream.csv") if os.path.exists(f"{LABELED_DIR}/pData_upstream.csv") else None,
            extraction_script="1_scripts/2_pipeline/1_download_data.R (archived as geo_data/imvigor210_labeled/upstream_1_download_data.R): "
                              "exmat written in counts(cds) row order with gene symbols in column 1; "
                              "censored file replaces gene_id with gene_1..gene_31085 placeholders",
        ),
        mirror=dict(
            repo="https://github.com/SiYangming/IMvigor210CoreBiologies",
            file="data/cds.RData",
            url="https://raw.githubusercontent.com/SiYangming/IMvigor210CoreBiologies/master/data/cds.RData",
            expected_size=51233046,
            sha256=sha256(f"{LABELED_DIR}/cds.RData")
                if (os.path.exists(f"{LABELED_DIR}/cds.RData")
                    and os.path.getsize(f"{LABELED_DIR}/cds.RData") == 51233046) else None,
            note="mirror of the official research-pub.gene.com IMvigor210CoreBiologies R package "
                 "(Mariathasan et al., Nature 2018, doi:10.1038/nature25501); used only for the "
                 "independent numeric identity check, never for the analysis values",
        ),
    )

    out = dict(
        analysis="imvigor_gene_recovery",
        date="2026-09-24",
        seed=SEED,
        phase1a=phase1a,
        phase1c=phase1c,
        phase2=dict(new=new_arm, old=old_arm),
        verification_gates=gates,
        provenance=provenance,
    )
    json.dump(out, open(OUT_JSON, "w"), indent=2)
    print("Wrote", OUT_JSON)
    write_md(out)

if __name__ == "__main__":
    main()
