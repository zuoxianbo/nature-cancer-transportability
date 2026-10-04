#!/usr/bin/env python3
""" phase 1a -- panel availability breakdown: entrez vs symbol matching."""
import csv
import json

ROOT = "/Users/zuoxianbo/Desktop/SCI论文/Sandy_Nature"
IMV_FDATA = f"{ROOT}/data/imvigor210/fData.csv"
G322_CSV = f"{ROOT}/outputs/analysis_base/supplementary/Supplementary_Table_1_322_genes.csv"

fd = list(csv.DictReader(open(IMV_FDATA)))
fd_ent = {int(r["entrez_id"]) for r in fd}
fd_sym = {r["symbol"].upper() for r in fd}
fd_ent_by_sym = {r["symbol"].upper(): int(r["entrez_id"]) for r in fd}

csyms, entrz = [], []
with open(G322_CSV) as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split(",")
        if p[0].isdigit():
            entrz.append(int(p[0]))
            csyms.append(p[1].strip().upper())

n_ent = sum(1 for e in entrz if e in fd_ent)
n_sym = sum(1 for s in csyms if s in fd_sym)
print("panel size:", len(csyms))
print("panel entrez found in fData entrez set:", n_ent)
print("panel symbol found in fData symbol set:", n_sym)

# mismatch breakdown
ent_missing = [(s, e) for s, e in zip(csyms, entrz) if e not in fd_ent]
sym_missing = [s for s in csyms if s not in fd_sym]
print("\npanel genes with entrez NOT in fData:", len(ent_missing))
print("panel genes with symbol NOT in fData:", len(sym_missing), sym_missing)
# among entrez-missing: is symbol present, and what entrez does fData give?
both = []
for s, e in ent_missing:
    if s in fd_sym:
        fe = fd_ent_by_sym[s]
        both.append((s, e, fe))
print("entrez-missing but symbol-present:", len(both))
for s, e, fe in both[:20]:
    print(f"  {s:12s} panel_entrez={e:10d} fData_entrez={fe:10d}")
# entrez present but pointing at a DIFFERENT symbol in fData?
fd_sym_by_ent = {int(r["entrez_id"]): r["symbol"].upper() for r in fd}
conflict = [(s, e, fd_sym_by_ent[e]) for s, e in zip(csyms, entrz)
            if e in fd_ent and fd_sym_by_ent[e] != s]
print("\nentrez present but fData symbol differs:", len(conflict))
for s, e, fs in conflict[:20]:
    print(f"  panel {s:12s} entrez={e:10d} -> fData symbol {fs}")
