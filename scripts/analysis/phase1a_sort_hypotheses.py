#!/usr/bin/env python3
""" phase 1a -- local de-anonymization: sorting-hypothesis sex-anchor tests.

For each candidate ordering hypothesis of exmat rows vs fData genes, build
entrez->row mapping and test:
  - DDX3Y/UTY/RPS4Y1 (Y-linked): female mean ~ 0, male mean > 0
  - XIST: female/male mean ratio > 5
  - GAPDH/ACTB (housekeeping): high nonzero expression
  - global: rows that are male-specific in exmat should map to Y-linked symbols
"""
import csv
import json
import numpy as np

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
OUT = f"{ROOT}/outputs/analysis_base/data/phase1a_results.json"

# ---------------------------------------------------------------- load
with open(IMV_EXMAT) as f:
    rd = csv.reader(f)
    hdr = next(rd)
    sids = [x.strip('"') for x in hdr[1:]]
    rows = []
    for r in rd:
        rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
M = np.array(rows, float)  # (31085 x 348)
print("exmat:", M.shape)

fd = list(csv.DictReader(open(IMV_FDATA)))
fd_ent = [r["entrez_id"] for r in fd]
fd_sym = [r["symbol"].upper() for r in fd]
fd_ent_int = [int(e) for e in fd_ent]

pd = list(csv.DictReader(open(IMV_PDATA)))
pmap = {r["X"]: r for r in pd}
sex = np.array([pmap[s]["Sex"] for s in sids])
female = sex == "F"
male = sex == "M"
print("F:", female.sum(), "M:", male.sum())

# ---------------------------------------------------------------- anchors
ANCHORS = {
    "DDX3Y": 8653, "UTY": 7404, "RPS4Y1": 6255,
    "XIST": 7503, "GAPDH": 2597, "ACTB": 60,
}
# a broader Y-linked symbol panel (HGNC chrY protein-coding classics)
Y_SYMS = {
    "DDX3Y", "UTY", "RPS4Y1", "KDM5D", "ZFY", "EIF1AY", "TMSB4Y", "USP9Y",
    "DAZ1", "DAZ2", "DAZ3", "DAZ4", "RBMY1A1", "RBMY1B", "RBMY1D", "RBMY1E",
    "RBMY1F", "RBMY1J", "PRKY", "CDY1", "CDY1B", "CDY2A", "CDY2B", "HSFY1",
    "HSFY2", "TSPY1", "TSPY2", "TSPY3", "TSPY4", "TSPY8", "TSPY10", "TSPY12",
    "NLGN4Y", "AMELY", "TXLNG2P", "BPY2", "VCY", "VCY3", "XKRY", "HSFY",
    "PRY", "PRY2", "PRY3", "PRY4", "TTY1", "TTY2", "TTY3", "TTY4", "TTY5",
    "SRY", "TGIF2LY", "UTY2", "ZFY2", "GYG2P", "CSPG4P1Y", "DDX3Y2",
    "SEPTIN4", "SEPTIN5", "ZRSR2Y", "GLYCTK", "PCDH11Y", "GTP2A3P", "WWC3",
    "L1CAP2", "TBL1Y", "TSPY13P", "CHDCH1", "ENOX1", "AC002444.2",
}

# ---------------------------------------------------------------- diagnostics independent of hypothesis
rowsum = M.sum(1)
print("exmat rowsum monotonic nondecreasing:", bool(np.all(np.diff(rowsum) >= 0)),
      "nonincreasing:", bool(np.all(np.diff(rowsum) <= 0)))
fmean = M[:, female].mean(1)
mmean = M[:, male].mean(1)
male_spec = (fmean == 0) & (mmean >= 5)
print("rows female==0 & male>=5:", int(male_spec.sum()))
fem_bias = fmean / np.maximum(mmean, 1e-9)
fem_rows = np.where((fem_bias > 5) & (fmean > 1))[0]
print("rows f/m>5 & fmean>1:", len(fem_rows))

# ---------------------------------------------------------------- hypothesis machinery
def mapping_from_order(order_idx):
    """order_idx: list of fData row indices in exmat-row order (exmat row j == fData row order_idx[j]).
    Returns entrez_int -> exmat row index."""
    ent2row = {}
    for j, fi in enumerate(order_idx):
        ent2row[fd_ent_int[fi]] = j
    return ent2row

def test_hypothesis(name, ent2row):
    res = {"hypothesis": name, "anchors": {}, "pass": False}
    ok = True
    for g, e in ANCHORS.items():
        if e not in ent2row:
            res["anchors"][g] = {"found": False}
            ok = False
            continue
        r = ent2row[e]
        row = M[r]
        mf = float(row[female].mean())
        mm = float(row[male].mean())
        res["anchors"][g] = {
            "found": True, "exmat_row_1based": r + 1,
            "female_mean": mf, "male_mean": mm,
            "rowsum": float(row.sum()),
        }
    # PASS criteria
    for g in ("DDX3Y", "UTY", "RPS4Y1"):
        a = res["anchors"][g]
        if not a.get("found") or not (a["female_mean"] <= 1.0 and a["male_mean"] >= 5.0):
            ok = False
    a = res["anchors"].get("XIST", {})
    if not a.get("found") or not (a["female_mean"] / max(a["male_mean"], 1e-9) > 5):
        ok = False
    for g in ("GAPDH", "ACTB"):
        a = res["anchors"][g]
        if not a.get("found") or not (a["rowsum"] > 1e5):
            ok = False
    # global Y-enrichment: male-specific rows -> symbols under this mapping
    row2sym = {j: fd_sym[fi] for j, fi in enumerate(order_idx)} if ent2row and "_order" in res else None
    res["pass"] = bool(ok)
    return res

def full_test(name, order_idx):
    """order_idx[j] = fData row index for exmat row j."""
    ent2row = mapping_from_order(order_idx)
    res = test_hypothesis(name, ent2row)
    # global Y-enrichment
    rows_ms = np.where(male_spec)[0]
    syms_ms = [fd_sym[order_idx[j]] for j in rows_ms]
    hits = sorted(set(syms_ms) & Y_SYMS)
    res["male_specific_rows"] = int(len(rows_ms))
    res["male_specific_y_hits"] = hits
    res["y_hit_rate"] = len(hits) / max(len(syms_ms), 1)
    # XIST-global: top fem-bias rows mapped
    top = fem_rows[np.argsort(-fem_bias[fem_rows])][:10]
    res["top_fem_bias_symbols"] = [
        {"symbol": fd_sym[order_idx[j]], "row": int(j + 1),
         "f_over_m": float(fem_bias[j])} for j in top
    ]
    # housekeeping sanity: what maps to top-50 rowsum rows
    top50 = np.argsort(-rowsum)[:20]
    res["top_rowsum_symbols"] = [
        {"symbol": fd_sym[order_idx[j]], "rowsum": float(rowsum[j])} for j in top50
    ]
    return res

# ---------------------------------------------------------------- hypotheses
n = len(fd)
hyps = {}

# (0) position: exmat row j == fData row j  (== entrez string sort, the falsified one)
hyps["position_fData"] = list(range(n))

# (i) == fData sorted by entrez string  (identical to position since fData IS string-sorted)
# already covered by position_fData; keep explicit label for the record

# (ii) == fData sorted by symbol alphabetical
hyps["symbol_asc"] = sorted(range(n), key=lambda i: fd_sym[i])

# (iii) == fData sorted by entrez numeric descending
hyps["entrez_numeric_desc"] = sorted(range(n), key=lambda i: -fd_ent_int[i])

# (iv) other systematic orderings
hyps["entrez_numeric_asc"] = sorted(range(n), key=lambda i: fd_ent_int[i])
hyps["symbol_desc"] = sorted(range(n), key=lambda i: fd_sym[i], reverse=True)
fd_len = [int(r["length"]) for r in fd]
hyps["length_asc"] = sorted(range(n), key=lambda i: fd_len[i])
hyps["length_desc"] = sorted(range(n), key=lambda i: -fd_len[i])
fd_nex = [int(r["n_exons"]) for r in fd]
hyps["n_exons_asc"] = sorted(range(n), key=lambda i: fd_nex[i])
fd_sym2 = [r["Symbol"].upper() for r in fd]
hyps["Symbol_col_asc"] = sorted(range(n), key=lambda i: fd_sym2[i])

results = {}
for name, order in hyps.items():
    r = full_test(name, order)
    results[name] = r
    print("=" * 70)
    print("HYPOTHESIS:", name, "PASS:", r["pass"])
    for g, a in r["anchors"].items():
        if a.get("found"):
            print(f"  {g:8s} row={a['exmat_row_1based']:6d} f={a['female_mean']:10.2f} "
                  f"m={a['male_mean']:10.2f} rowsum={a['rowsum']:12.0f}")
        else:
            print(f"  {g:8s} NOT FOUND")
    print("  male-specific rows:", r["male_specific_rows"],
          "Y-hits:", r["male_specific_y_hits"][:15],
          "hit-rate:", round(r["y_hit_rate"], 3))

json.dump(results, open(OUT, "w"), indent=2)
print("\nWrote", OUT)
