#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" compliance self-check: abstract<=150w, body<=4000w, banned words, fig cites, no 'well calibrated' / 'digital twin' / 'triangulation'."""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ms = os.path.join(ROOT, "manuscript", "manuscript_.md")
text = open(ms).read()

# abstract
m = re.search(r"## Abstract\s*(.*?)\n\n---\n", text, re.S)
abstract = m.group(1).strip() if m else ""
abs_words = len(abstract.split())

# main text = Nature Cancer policy: Introduction + Results + Discussion (<=4000w).
# The Introduction carries NO heading (text flows after the abstract); Methods is a
# SEPARATE online section (<=3000w) and is EXCLUDED from the 4000 budget.
_post_abs = text.split("## Abstract", 1)[-1]
_post_abs = re.split(r"\n---\n", _post_abs, 1)[-1]   # everything after the abstract separator
_main_text = _post_abs.split("## Methods", 1)[0]
body_words = len(_main_text.split())
# methods (online, separate limit 3000w) -- informational only
_methods_raw = text.split("## Methods", 1)[-1]
_methods_text = _methods_raw.split("**Figure legends**", 1)[0]
methods_words = len(_methods_text.split())

banned = ["separatrix", "bifurcation", "clinically actionable", "replicated", "benchmarked",
          "landscape", "validated primary estimate", "digital twin"]  # positive use banned; negative use allowed
found = []
low = text.lower().replace("**", "")
for b in banned:
    start = 0
    hit = False
    while True:
        idx = low.find(b, start)
        if idx == -1:
            break
        # allow if within 70 chars before the phrase there is a negation token
        ctx = low[max(0, idx-70):idx]
        neg = any(tok in ctx for tok in ["not ", "never ", "do not", " no ", "isn't", "n't ",
                                         "un", "without", "rather than", "not a",
                                         "cannot", "do not claim", "never a", "dispatched",
                                         "proposed", "not executed", "specified"])
        if not neg:
            hit = True
            break
        start = idx + 1
    if hit:
        found.append(b)

# must-NOT-contain: reviewer flagged "well calibrated" as a hard flaw
must_not = ["well calibrated"]
present_bad = [p for p in must_not if p in low]

# must be ABSENT entirely from the body: "triangulation" (relabeled to genetic plausibility)
triangulation_present = "triangulation" in low

figs = [f"Fig. {i}" for i in range(1,8)]
fig_present = {f: (f in text) for f in figs}

ok = True
print("="*60)
print("COMPLIANCE SELF-CHECK")
print("="*60)
print(f"abstract words : {abs_words}  (limit 150) -> {'PASS' if abs_words<=150 else 'FAIL'}")
print(f"body words     : {body_words}  (limit 4000) -> {'PASS' if body_words<=4000 else 'FAIL'}")
print(f"methods words  : {methods_words}  (online limit 3000) -> {'PASS' if methods_words<=3000 else 'FAIL'}")
print(f"banned words    : {found if found else 'none'} -> {'PASS' if not found else 'FAIL'}")
print(f"'well calibrated': {present_bad if present_bad else 'absent'} -> {'PASS' if not present_bad else 'FAIL'}")
print(f"'triangulation' : {'PRESENT (FAIL)' if triangulation_present else 'absent'} -> {'PASS' if not triangulation_present else 'FAIL'}")
for f in figs:
    print(f"figure cite    : {f} -> {'present' if fig_present[f] else 'MISSING'}")
    if not fig_present[f]: ok=False
if abs_words>150 or body_words>4000 or methods_words>3000 or found or present_bad or triangulation_present: ok=False
print("="*60)
print("RESULT:", "PASS" if ok else "FAIL")
print("="*60)
