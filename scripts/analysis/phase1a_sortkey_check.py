#!/usr/bin/env python3
""" phase 1a -- resolve sort-key ambiguity (raw vs upper) + immune checks + panel overlap."""
import csv
import json
import numpy as np
from scipy.stats import spearmanr

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
IMV_PDATA = f"{ROOT}/data/imvigor210/pData.csv"
G322 = f"{ROOT}/geo_data/322_shared_genes.json"

fd = list(csv.DictReader(open(IMV_FDATA)))
fd_sym_raw = [r["symbol"] for r in fd]
fd_sym_up = [s.upper() for s in fd_sym_raw]
fd_Symbol = [r["Symbol"] for r in fd]

# lowercase-containing symbols?
lower_syms = [(i, s) for i, s in enumerate(fd_sym_raw) if s != s.upper()]
print("symbols containing lowercase:", len(lower_syms))
for i, s in lower_syms[:40]:
    print("  fData row", i, repr(s), "entrez", fd[i]["entrez_id"], "Symbol col:", repr(fd_Symbol[i]))

# symbol col == Symbol col?
print("symbol col == Symbol col (raw):", fd_sym_raw == fd_Symbol)

# orderings
order_raw = sorted(range(len(fd)), key=lambda i: fd_sym_raw[i])
order_up = sorted(range(len(fd)), key=lambda i: fd_sym_up[i])
diff = [j for j in range(len(fd)) if order_raw[j] != order_up[j]]
print("rows where raw-sort position != upper-sort position:", len(diff))

# among lowercase symbols, their neighbors in upper sort vs raw sort
# the only ambiguity: relative order among symbols equal ignoring case? check collisions
from collections import Counter
cu = Counter(fd_sym_up)
print("collisions after upper():", {s: c for s, c in cu.items() if c > 1})

# TC/IC levels distribution
pd = list(csv.DictReader(open(IMV_PDATA)))
from collections import Counter as C2
print("TC.Level values:", C2(r["TC.Level"] for r in pd))
print("IC.Level values:", C2(r["IC.Level"] for r in pd))
