#!/usr/bin/env python3
""" phase 1a refinement -- deep verification of the symbol_asc hypothesis."""
import csv
import json
import numpy as np

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
IMV_EXMAT = f"{ROOT}/data/imvigor210/exmat.csv"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"

with open(IMV_EXMAT) as f:
    rd = csv.reader(f)
    hdr = next(rd)
    sids = [x.strip('"') for x in hdr[1:]]
    rows = []
    for r in rd:
        rows.append([float(x) if x not in ("", "NA") else 0.0 for x in r[1:]])
M = np.array(rows, float)

fd = list(csv.DictReader(open(IMV_FDATA)))
fd_ent = [int(r["entrez_id"]) for r in fd]
fd_sym_raw = [r["symbol"] for r in fd]
fd_sym = [s.upper() for s in fd_sym_raw]

pd = list(csv.DictReader(open(IMV_PDATA)))
pmap = {r["X"]: r for r in pd}
sex = np.array([pmap[s]["Sex"] for s in sids])
female = sex == "F"
male = sex == "M"

# symbol_asc order: exmat row j == j-th fData row in symbol-sorted order
order = sorted(range(len(fd)), key=lambda i: fd_sym[i])
sym_of_row = {j: fd_sym[order[j]] for j in range(len(fd))}
ent_of_row = {j: fd_ent[order[j]] for j in range(len(fd))}

# check for duplicate symbols (order ambiguity)
from collections import Counter
cdup = Counter(fd_sym)
dups = {s: c for s, c in cdup.items() if c > 1}
print("duplicate symbols in fData:", len(dups), list(dups.items())[:10])

# RPS4Y1 present?
for i, s in enumerate(fd_sym):
    if s == "RPS4Y1":
        print("RPS4Y1 found in fData: entrez", fd_ent[i])

# ---- 1) female nonzero pattern for DDX3Y / UTY rows
for g in ("DDX3Y", "UTY", "RPS4Y1", "KDM5D", "ZFY", "EIF1AY", "TMSB4Y", "USP9Y"):
    r = None
    for j in range(len(fd)):
        if sym_of_row[j] == g:
            r = j
            break
    if r is None:
        print(f"{g}: not in symbol set")
        continue
    row = M[r]
    fnz = int((row[female] > 0).sum())
    mnz = int((row[male] > 0).sum())
    print(f"{g:8s} row={r+1:6d} f_mean={row[female].mean():10.2f} m_mean={row[male].mean():10.2f} "
          f"f_nonzero={fnz}/76 m_nonzero={mnz}/272 f/m={row[female].mean()/max(row[male].mean(),1e-9):.4f}")

# ---- 2) top male-biased rows (m high, f/m tiny) -> should be Y genes
fmean = M[:, female].mean(1)
mmean = M[:, male].mean(1)
ratio = fmean / np.maximum(mmean, 1e-9)
cand = np.where((mmean >= 50) & (ratio < 0.05))[0]
cand = cand[np.argsort(-mmean[cand])]
print("\ntop male-biased rows (mmean>=50, f/m<0.05):", len(cand))
for j in cand[:25]:
    print(f"  row={j+1:6d} sym={sym_of_row[j]:12s} m={mmean[j]:10.1f} f={fmean[j]:8.2f}")

# ---- 3) top female-biased rows -> XIST should be #1
candf = np.where((fmean >= 10) & (ratio > 5))[0]
candf = candf[np.argsort(-ratio[candf])]
print("\ntop female-biased rows (fmean>=10, f/m>5):", len(candf))
for j in candf[:15]:
    print(f"  row={j+1:6d} sym={sym_of_row[j]:12s} f={fmean[j]:12.1f} m={mmean[j]:10.1f} f/m={ratio[j]:.1f}")

# ---- 4) immune-gene vs pData TC.Level / IC.Level correlation check
def level_num(v):
    return {"TC0": 0, "TC1": 1, "TC2": 2, "TC3": 3, "IC0": 0, "IC1": 1, "IC2": 2, "IC3": 3}.get(v, np.nan)

tc = np.array([level_num(pmap[s]["TC.Level"]) for s in sids], float)
ic = np.array([level_num(pmap[s]["IC.Level"]) for s in sids], float)
from scipy.stats import spearmanr
print("\nimmune correlation spot-checks (Spearman):")
for g in ("CD8A", "CD3D", "CD3E", "GZMB", "PRF1", "CD274", "PDCD1", "CTLA4", "MS4A1", "KRT5", "KRT14", "KRT17"):
    r = None
    for j in range(len(fd)):
        if sym_of_row[j] == g:
            r = j
            break
    if r is None:
        print(f"  {g}: not found")
        continue
    rho_tc, p_tc = spearmanr(M[r], tc)
    rho_ic, p_ic = spearmanr(M[r], ic)
    print(f"  {g:8s} rho(TC)={rho_tc:+.3f} (p={p_tc:.1e})  rho(IC)={rho_ic:+.3f} (p={p_ic:.1e})")

# ---- 5) housekeeping top-rowsum symbols
rowsum = M.sum(1)
top50 = np.argsort(-rowsum)[:30]
print("\ntop-30 rowsum symbols:")
for j in top50:
    print(f"  row={j+1:6d} sym={sym_of_row[j]:12s} rowsum={rowsum[j]:14.0f}")

# ---- 6) symbol case check: is sort stable wrt original case?
order_raw = sorted(range(len(fd)), key=lambda i: fd_sym_raw[i])
print("\nsymbol_asc(raw-case) == symbol_asc(upper)?", order_raw == order)
