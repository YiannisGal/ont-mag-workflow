#!/usr/bin/env bash
# 92_mcra_panel_composition.sh -- the decoy panel is not what the figure says
# it is, and the short alignment is a symptom rather than bad luck.
#
# WHAT STEP 91 PRINTED. Two things, and the second is the bigger one.
#
# (1) TWO SEQUENCES IN THE MCR PANEL ARE BACTERIAL RESTRICTION ENZYMES.
#
#       WP_484489188.1  5-methylcytosine-specific endonuclease McrA [Escherichia coli]
#       WP_471601478.1  type IV restriction endonuclease McrA [Streptomyces sp.]
#
#     "McrA" is also the name of a bacterial restriction endonuclease that has
#     nothing whatever to do with methanogenesis. Neither of the two search
#     terms visible in step 91's grep can produce them, so a THIRD term must
#     exist in the fetch script, and section A prints that block verbatim
#     instead of through a filter. iqtree.log already convicts these two on
#     other grounds: the E. coli sequence is 62,59% gaps in the alignment,
#     against 0,00% for essentially every genuine MCR.
#
# (2) THE 36-SEQUENCE "ACR DECOY PANEL" IS MOSTLY NOT THE ALPHA SUBUNIT.
#
#     Reading refs/acr.faa's own deflines, the panel contains
#
#       "methyl coenzyme M reductase system, component A2"   -- a DIFFERENT
#            protein (McrC/AtwA-like), not a subunit of the enzyme at all
#       "... subunit beta", "... subunit gamma"              -- paralogues of
#            each other, not of the alpha subunit
#       "... subunit alpha"                                  -- the only
#            sequences that are actually comparable to our queries
#
#     Section B counts them. If the alpha count is small, then the control the
#     figure advertises -- "can you tell MCR alpha from ACR alpha?" -- is being
#     done by a handful of sequences while the rest of the panel answers a much
#     easier question that nobody asked: can you tell an alpha subunit from a
#     gamma subunit. THE CAPTION OF fig14 CURRENTLY SAYS 36. That number must
#     become whatever section B returns.
#
#     It also explains the alignment. 105 sequences over 139 amino-acid sites
#     looked like an unlucky trim. If the alignment contains alpha, beta, gamma
#     and a non-subunit protein, then 139 columns is what survives aligning
#     things that are not homologous over most of their length -- a symptom of
#     the reference set, not of the trimming.
#
# WHY IT MATTERS AND WHY IT MAY STILL BE FINE. Every query's nearest ACR
# neighbour was reported as Candidatus Syntropharchaeum, which is where the
# genuine alpha sequences are. If the nearest ACR is an ALPHA in every case,
# the measured margin stands and only the wording changes. If it is a component
# A2 or a gamma, the margin is measured against the wrong thing and the
# comparison has to be redone against alpha only. Section D decides that from
# the distance matrix, per query, by name.
#
# Nothing is fetched, edited or deleted here.
set -u
WS2="${WS2}"
OUT="$WS2/92_panel_composition"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. the search terms, verbatim and unfiltered"
echo "=============================================================="
echo "   Step 91 grepped this block and the grep could not match a term that"
echo "   does not contain the words it was looking for. Printed whole here."
echo
for f in "$WS2"/scripts/15_marker_trees.sh; do
    [ -f "$f" ] || { echo "   MISSING $f"; continue; }
    echo "   --- ${f#$WS2/}, lines 90-135 ---"
    sed -n '90,135p' "$f" | cat -n | sed 's/^/      /'
done

echo
echo "=============================================================="
echo "B. what each panel actually contains, by subunit"
echo "=============================================================="
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections
WS2, OUT = sys.argv[1], sys.argv[2]
RD = os.path.join(WS2, "15_marker_trees", "mcrA", "refs")


def classify(desc):
    d = desc.lower()
    # order matters: "component A2" is not a subunit at all, and must be
    # tested before the subunit words, because its defline also contains
    # "methyl coenzyme M reductase"
    if "component a2" in d or "component, a2" in d:
        return "component A2 (not a subunit)"
    if "endonuclease" in d or "restriction" in d:
        return "RESTRICTION ENZYME (not MCR at all)"
    if re.search(r"\b(subunit )?alpha\b", d) or d.endswith(" alpha"):
        return "alpha"
    if re.search(r"\b(subunit )?beta\b", d) or d.endswith(" beta"):
        return "beta"
    if re.search(r"\b(subunit )?gamma\b", d) or d.endswith(" gamma"):
        return "gamma"
    if "sulfoethylthiotransferase subunit alpha" in d:
        return "alpha"
    return "UNCLASSIFIED"


def org_of(desc):
    m = re.search(r"\[([^\]]+)\]", desc)
    return m.group(1) if m else ""


ROWS = {}
for fn, panel in (("mcr.faa", "MCR"), ("acr.faa", "ACR")):
    p = os.path.join(RD, fn)
    if not os.path.exists(p):
        print("   MISSING %s" % p)
        continue
    rows = []
    for ln in open(p, errors="replace"):
        if not ln.startswith(">"):
            continue
        h = ln[1:].rstrip("\n")
        accn = h.split()[0]
        desc = h[len(accn):].strip()
        rows.append((accn, classify(desc), org_of(desc), desc))
    ROWS[panel] = rows
    c = collections.Counter(r[1] for r in rows)
    print()
    print("   %s panel -- refs/%s, %d sequence(s)" % (panel, fn, len(rows)))
    for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
        print("      %-36s %3d" % (k, v))
    # the alpha members, named, because they are the only comparable ones
    al = [r for r in rows if r[1] == "alpha"]
    print("      -> %d alpha-subunit sequence(s)" % len(al))
    if panel == "ACR":
        g = collections.Counter(" ".join(r[2].split()[:2]) for r in al)
        print("         from %d organism name(s): %s"
              % (len(g), ", ".join("%s x%d" % (k, v) for k, v in g.most_common())))
        for accn, kind, org, desc in al:
            print("         %-18s %s" % (accn, org))
    bad = [r for r in rows if r[1].startswith("RESTRICTION")]
    for accn, kind, org, desc in bad:
        print("      !! %-18s %s" % (accn, desc[:110]))

with open(os.path.join(OUT, "panel_composition.tsv"), "w") as fh:
    fh.write("panel\taccession\tclass\torganism\tdescription\n")
    for panel, rows in ROWS.items():
        for accn, kind, org, desc in rows:
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (panel, accn, kind, org, desc))
print()
print("   -> 92_panel_composition/panel_composition.tsv")

nal = len([r for r in ROWS.get("ACR", []) if r[1] == "alpha"])
print()
print("   THE NUMBER fig14's CAPTION MUST USE: %d alpha-subunit decoy(s)," % nal)
print("   not 36. The other %d ACR-panel entries are other proteins."
      % (len(ROWS.get("ACR", [])) - nal))
PYEOF

echo
echo "=============================================================="
echo "C. how much of each sequence is gap in the alignment"
echo "=============================================================="
echo "   IQ-TREE prints a gap percentage per sequence. A sequence that is"
echo "   mostly gaps is not homologous to the rest over the aligned region,"
echo "   and it is dragging the alignment: every column it cannot fill is a"
echo "   column the trimmer is likely to drop for everyone."
echo
"$PY" - "$WS2" <<'PYEOF'
import sys, os, re
WS2 = sys.argv[1]
lg = os.path.join(WS2, "15_marker_trees", "mcrA", "mcrA.log")
if not os.path.exists(lg):
    lg = os.path.join(WS2, "15_marker_trees", "mcrA", "iqtree.log")
if not os.path.exists(lg):
    print("   MISSING the IQ-TREE log"); sys.exit(0)
rows = []
for ln in open(lg, errors="replace"):
    m = re.match(r"\s*\d+\s+(\S+)\s+([\d.]+)%\s+(\w+)\s+([\d.]+)%", ln)
    if m:
        rows.append((float(m.group(2)), m.group(1), m.group(3), m.group(4)))
if not rows:
    print("   no per-sequence gap table found in %s" % os.path.basename(lg))
    sys.exit(0)
rows.sort(reverse=True)
print("   %d sequence(s) with a gap figure; the 20 worst:" % len(rows))
print("   %8s  %-46s %s" % ("gap %", "sequence", "composition test"))
for gap, nm, verdict, pval in rows[:20]:
    print("   %7.2f%%  %-46s %s" % (gap, nm, verdict))
import statistics
print()
print("   median gap %.2f%%, mean %.2f%%"
      % (statistics.median(r[0] for r in rows),
         sum(r[0] for r in rows) / len(rows)))
for panel in ("QUERY", "MCR", "ACR"):
    sub = [r[0] for r in rows if r[1].startswith(panel + "|")]
    if sub:
        print("   %-6s n=%3d  gap %% min %.2f  median %.2f  max %.2f"
              % (panel, len(sub), min(sub),
                 statistics.median(sub), max(sub)))
PYEOF

echo
echo "=============================================================="
echo "D. per query: nearest ACR ALPHA, by name, from the distance matrix"
echo "=============================================================="
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re
WS2, OUT = sys.argv[1], sys.argv[2]
D = os.path.join(WS2, "15_marker_trees", "mcrA")
RD = os.path.join(D, "refs")


def classify(desc):
    d = desc.lower()
    if "component a2" in d or "component, a2" in d:
        return "A2"
    if "endonuclease" in d or "restriction" in d:
        return "RESTRICTION"
    if re.search(r"\balpha\b", d):
        return "alpha"
    if re.search(r"\bbeta\b", d):
        return "beta"
    if re.search(r"\bgamma\b", d):
        return "gamma"
    return "?"


KIND, ORG = {}, {}
for fn in ("mcr.faa", "acr.faa"):
    p = os.path.join(RD, fn)
    if not os.path.exists(p):
        continue
    for ln in open(p, errors="replace"):
        if not ln.startswith(">"):
            continue
        h = ln[1:].rstrip("\n")
        a = h.split()[0]
        desc = h[len(a):].strip()
        KIND[a] = classify(desc)
        m = re.search(r"\[([^\]]+)\]", desc)
        ORG[a] = m.group(1) if m else ""

md = os.path.join(D, "mcrA.mldist")
lines = [l.rstrip("\n") for l in open(md, errors="replace") if l.strip()]
names, M = [], []
for l in lines[1:]:
    f = l.split()
    if len(f) < 2:
        continue
    names.append(f[0]); M.append([float(x) for x in f[1:]])
idx = {n: i for i, n in enumerate(names)}


def panel(t):
    return t.split("|")[0]


def accn(t):
    return t.split("|")[-1]


qs = [n for n in names if panel(n) == "QUERY"]
print("   %-20s %-34s %8s   %-34s %8s   %s"
      % ("query", "nearest MCR alpha (archaeal)", "dist",
         "nearest ACR ALPHA", "dist", "ratio"))
rat = []
for q in qs:
    i = idx[q]

    def best(pred):
        c = [(M[i][j], names[j]) for j in range(len(names))
             if j != i and pred(names[j])]
        return min(c) if c else None

    m_ok = best(lambda t: panel(t) == "MCR" and KIND.get(accn(t)) == "alpha")
    a_al = best(lambda t: panel(t) == "ACR" and KIND.get(accn(t)) == "alpha")
    a_any = best(lambda t: panel(t) == "ACR")
    f = q.split("|")
    if m_ok and a_al:
        r = a_al[0] / m_ok[0] if m_ok[0] else float("inf")
        rat.append(r)
        print("   %-20s %-34s %8.4f   %-34s %8.4f   x%.1f"
              % (f[1].replace("ETO_", ""),
                 ("%s %s" % (ORG.get(accn(m_ok[1]), "?"), accn(m_ok[1])))[:34],
                 m_ok[0],
                 ("%s %s" % (ORG.get(accn(a_al[1]), "?"), accn(a_al[1])))[:34],
                 a_al[0], r))
    if a_any and a_al and a_any[1] != a_al[1]:
        print("        (the nearest ACR of ANY kind is %s, a %s, at %.4f -- so the"
              % (accn(a_any[1]), KIND.get(accn(a_any[1]), "?"), a_any[0]))
        print("         margin quoted before was measured against that, not against alpha)")
if rat:
    print()
    print("   against ALPHA decoys only: %.1f to %.1f times further"
          % (min(rat), max(rat)))

print()
print("   and the ten nearest neighbours again, now labelled by what each")
print("   reference actually is:")
for q in qs:
    i = idx[q]
    d = sorted(((M[i][j], names[j]) for j in range(len(names)) if j != i))[:10]
    f = q.split("|")
    comp = {}
    for _, nm in d:
        k = "QUERY" if panel(nm) == "QUERY" else "%s %s" % (panel(nm), KIND.get(accn(nm), "?"))
        comp[k] = comp.get(k, 0) + 1
    print("   %-20s %s" % (f[1].replace("ETO_", ""),
                           ", ".join("%s x%d" % (k, v) for k, v in sorted(comp.items()))))
PYEOF

echo
echo "=============================================================="
echo "E. the same question of the mttB panels"
echo "=============================================================="
echo "   mttB_verdicts.tsv reports every query CONFIRMED with its closest"
echo "   reference at 1,22 to 1,28 substitutions per site, and the SAME"
echo "   accession closest to all five. A closest relative that far away, and"
echo "   identical for every query, is what you see when the reference set is"
echo "   not close to the queries at all."
echo
for f in "$WS2"/15_marker_trees/mttB/refs/*.faa; do
    [ -f "$f" ] || continue
    echo "   --- ${f#$WS2/}: $(grep -c '^>' "$f") sequence(s) ---"
    grep "^>" "$f" | sed 's/^/      /' | head -50
done

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   B gives fig14's caption its real decoy count, and tells us whether"
echo "   the mcrA reference set has to be rebuilt before the figure is final."
echo "   D decides whether the margin already measured is sound: if every"
echo "   query's nearest ACR is an ALPHA subunit, only the wording changes."
echo "   C tells us whether the 139-site alignment is a trimming accident or"
echo "   a consequence of aligning proteins that are not homologous."
