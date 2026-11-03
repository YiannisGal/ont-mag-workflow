#!/usr/bin/env python3
"""Tip tables for the SP05 mcrA / dsrA / dsrB trees (S_mcrA_tree, S_dsr_trees).

Categories follow fig14_15_mcra.py's own subunit rule (description text of the reference
FASTA, not tree position). The two restriction-enzyme "McrA" (WP_484489188.1,
WP_471601478.1) are flagged DROP and pruned by the plotting script.
"""
import os, re, sys
WS2 = "${WS2}"
T = os.path.join(WS2, "99_thesis_audit/SP05_dsr_mcr/trees")
OUT = sys.argv[1]
DROP = {"WP_484489188.1", "WP_471601478.1"}

desc = {}
for fn in ("acr.faa", "mcr.faa"):
    for ln in open(os.path.join(WS2, "15_marker_trees/mcrA/refs", fn)):
        if ln.startswith(">"):
            acc, _, d = ln[1:].strip().partition(" ")
            desc[acc] = d

def subunit(d):
    d = d.lower()
    if "component a2" in d or re.search(r"\ba2\b", d): return "A2"
    for s in ("alpha", "beta", "gamma"):
        if re.search(r"\b%s\b" % s, d): return s
    return "other"

def species(d):
    m = re.search(r"\[([^\]]+)\]", d)
    return m.group(1) if m else ""

tips = [l[1:].strip().split()[0] for l in open(os.path.join(T, "mcrA_in.faa")) if l.startswith(">")]
with open(os.path.join(OUT, "mcrA_tips.tsv"), "w") as f:
    f.write("tip\tcat\tlabel\n")
    for t in tips:
        p = t.split("|")
        if p[0] == "QUERY":
            c, lab = "query", " ".join(p[1:2] + p[3:]) if len(p) > 3 else " ".join(p[1:])
        else:
            acc = p[1]; d = desc.get(acc, "")
            if acc in DROP: c = "DROP"
            elif p[0] == "MCR": c = "MCR alpha"
            else:
                su = subunit(d)
                c = {"alpha": "ACR alpha", "A2": "ACR component A2"}.get(su, "ACR other subunits")
            lab = "%s %s" % (species(d) or "", acc)
        f.write("%s\t%s\t%s\n" % (t, c, lab.strip()))

for g in ("dsrA", "dsrB"):
    tips = [l[1:].strip().split()[0] for l in open(os.path.join(T, "%s_in.faa" % g)) if l.startswith(">")]
    with open(os.path.join(OUT, "%s_tips.tsv" % g), "w") as f:
        f.write("tip\tcat\tlabel\n")
        for t in tips:
            p = t.split("|")
            c = {"QUERY": "query", "REF_RED": "reductive", "REF_OXI": "oxidative (rDsr)"}[p[0]]
            lab = " ".join(p[1:]).replace("_", " ") if p[0] != "QUERY" else " ".join(p[1:])
            f.write("%s\t%s\t%s\n" % (t, c, lab))
print("tip tables written")
