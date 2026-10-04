#!/usr/bin/env python3
""" phase 1a -- does the raw-vs-upper sort ambiguity affect the 322 panel genes?
Plus immune-gene vs TC/IC correlation with fixed level parsing."""
import csv
import json
import numpy as np
from scipy.stats import spearmanr

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
G322_CSV = f"{ROOT}/outputs/analysis_base/supplementary/Supplementary_Table_1_322_genes.csv"
G322_JSON = f"{ROOT}/geo_data/322_shared_genes.json"

fd = list(csv.DictReader(open(IMV_FDATA)))
fd_sym_raw = [r["symbol"] for r in fd]
fd_sym_up = [s.upper() for s in fd_sym_raw]
fd_ent = [int(r["entrez_id"]) for r in fd]

order_raw = sorted(range(len(fd)), key=lambda i: fd_sym_raw[i])   # exmat row j == fData row order_raw[j]
order_up = sorted(range(len(fd)), key=lambda i: fd_sym_up[i])
# row index (exmat) of a gene, under each ordering:
row_raw = {fd_sym_up[fi]: j for j, fi in enumerate(order_raw)}
row_up = {fd_sym_up[fi]: j for j, fi in enumerate(order_up)}
ent_row_raw = {fd_ent[fi]: j for j, fi in enumerate(order_raw)}
ent_row_up = {fd_ent[fi]: j for j, fi in enumerate(order_up)}

# 322 panel
panel = json.load(open(G322_JSON))
print("322 json type:", type(panel), "len:", len(panel))
print("first entries:", panel[:5] if isinstance(panel, list) else list(panel.items())[:5])

csyms, entrz = [], []
with open(G322_CSV) as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split(",")
        if p[0].isdigit():
            entrz.append(int(p[0]))
            csyms.append(p[1].strip().upper())
print("panel from CSV:", len(csyms))

# availability + row consistency
n_found = 0
n_diff = 0
diff_genes = []
for s, e in zip(csyms, entrz):
    r_raw = ent_row_raw.get(e)
    r_up = ent_row_up.get(e)
    if r_raw is not None and r_up is not None:
        n_found += 1
        if r_raw != r_up:
            n_diff += 1
            diff_genes.append((s, e, r_raw, r_up))
    elif s in row_raw or s in row_up:
        # entrez mismatch but symbol present
        rr = row_raw.get(s)
        ru = row_up.get(s)
        n_found += 1
        if rr != ru:
            n_diff += 1
            diff_genes.append((s, e, rr, ru))
print("panel genes found in IMvigor fData (by entrez or symbol):", n_found, "/", len(csyms))
print("panel genes whose row differs raw vs upper:", n_diff)
print(diff_genes)

# ---- immune correlation with fixed level parsing
with open(IMV_EXMAT) as f:
    rd = csv.reader(f)
    hdr = next(rd)
    sids = [x.strip('"') for x in hdr[1:]]
    rows = []
    for r in rd:
        rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
M = np.array(rows, float)
pd = list(csv.DictReader(open(IMV_PDATA)))
pmap = {r["X"]: r for r in pd}

def level_num(v):
    m = {"TC0": 0, "TC1": 1, "TC2": 2, "TC2+": 2.5, "TC3": 3,
         "IC0": 0, "IC1": 1, "IC2": 2, "IC2+": 2.5, "IC3": 3}
    return m.get(v, np.nan)

tc = np.array([level_num(pmap[s]["TC.Level"]) for s in sids], float)
ic = np.array([level_num(pmap[s]["IC.Level"]) for s in sids], float)
mtc = ~np.isnan(tc)
mic = ~np.isnan(ic)

sym_of_row_up = {j: fd_sym_up[order_up[j]] for j in range(len(fd))}
print("\nimmune correlation (Spearman, upper-sort key):")
for g in ("CD8A", "CD3D", "CD3E", "GZMB", "PRF1", "CD274", "PDCD1", "CTLA4", "MS4A1",
          "KRT5", "KRT14", "KRT17", "CD8B", "CXCL9", "CXCL10"):
    j = None
    for jj in range(len(fd)):
        if sym_of_row_up[jj] == g:
            j = jj
            break
    if j is None:
        print(f"  {g}: not in IMvigor fData")
        continue
    rho_tc, p_tc = spearmanr(M[j][mtc], tc[mtc])
    rho_ic, p_ic = spearmanr(M[j][mic], ic[mic])
    print(f"  {g:8s} row={j+1:6d} rho(TC)={rho_tc:+.3f} (p={p_tc:.1e})  rho(IC)={rho_ic:+.3f} (p={p_ic:.1e})")
