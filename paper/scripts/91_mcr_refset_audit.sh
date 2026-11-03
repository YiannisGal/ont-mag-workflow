#!/usr/bin/env bash
# 91_mcr_refset_audit.sh -- the MCR reference panel may not be clean.
#
# WHAT fig14 FOUND WHILE DRAWING THE TREE.
#
# The 36 ACR sequences were put into the mcrA alignment as a DECOY PANEL, and
# the figure roots the tree on them. They are not a clade. Inside the smallest
# group that contains all of them there are also SEVEN sequences carrying the
# MCR label, and the organism names the script recovered for them include
#
#     Methanopyrus kandleri, Methanopyrus sp.   -- archaea, deep-branching
#     Candidatus Methanoglobus                  -- archaeon
#     Streptomyces sp., Escherichia coli        -- BACTERIA
#
# Methyl-coenzyme M reductase is an archaeal enzyme. A Streptomyces or an
# E. coli protein in an MCR reference set is not a phylogenetic result, it is a
# reference set that was assembled by matching protein NAMES. This project has
# already been bitten by exactly that: step 63 found seven bacterial proteins
# that an NCBI protein-name query had imported into the mtmB TARGET set, and
# removing them took 67 reads off that count.
#
# BUT THE ORGANISM NAMES ABOVE ARE MY OWN REGEX'S OUTPUT, NOT THE FILE'S.
# fig14 strips a defline down to what is in square brackets and then trims
# enzyme words out of it. A mangled description would produce exactly this
# symptom without anything being wrong in the data. So section A prints the
# RAW deflines, untouched, and every conclusion waits for them.
#
# WHAT TURNS ON IT.
#
#   * FINDINGS_REPORT section 5: "mcrA phylogenetic placement: CONFIRMED for 9
#     proteins across 7 genomes, ten nearest neighbours all MCR in every case."
#     If some MCR-labelled references are not MCR, "all MCR" is a statement
#     about labels rather than about enzymes. Section D recomputes the ten
#     nearest neighbours of every query FROM THE DISTANCE MATRIX and prints
#     which references they are, so the claim can be restated in terms of named
#     organisms instead of a panel tag.
#   * fig14's caption, which currently names the intruders.
#   * the mttB tree, built the same way, and fig15 if it is drawn.
#
# Nothing is deleted and no reference set is edited here. Rule 121.
set -u
WS2="${WS2}"
OUT="$WS2/91_mcr_refset_audit"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. the raw deflines of every reference in the mcrA alignment"
echo "=============================================================="
echo "   Printed verbatim. My own name extraction is not trusted here."
echo
for f in "$WS2"/15_marker_trees/mcrA/mcrA_all.faa; do
    [ -f "$f" ] || { echo "   MISSING $f"; continue; }
    echo "   ${f#$WS2/}: $(grep -c '^>' "$f") sequence(s)"
    grep "^>" "$f" | sed 's/^/      /'
done
echo
echo "   and whatever is in the refs/ subdirectory (rule 125: name the"
echo "   extension, a directory is not one kind of file):"
ls -la "$WS2/15_marker_trees/mcrA/refs" 2>/dev/null | awk 'NR>1 {printf "      %10s  %s\n", $5, $9}'
for f in "$WS2"/15_marker_trees/mcrA/refs/*.faa "$WS2"/15_marker_trees/mcrA/refs/*.fa \
         "$WS2"/15_marker_trees/mcrA/refs/*.tsv "$WS2"/15_marker_trees/mcrA/refs/*.txt; do
    [ -f "$f" ] || continue
    echo "   --- ${f#$WS2/} ---"
    if grep -q "^>" "$f" 2>/dev/null; then
        grep "^>" "$f" | head -80 | sed 's/^/      /'
    else
        head -40 "$f" | sed 's/^/      /'
    fi
done

echo
echo "=============================================================="
echo "B. how this reference set was built"
echo "=============================================================="
echo "   A set assembled by matching protein NAMES will contain whatever else"
echo "   carries that name. A set assembled from a curated profile or from"
echo "   named genomes will not. The command that made it decides how much"
echo "   the panel label is worth."
echo
for f in "$WS2"/scripts/*marker_tree* "$WS2"/scripts/*15_* "$WS2"/15_marker_trees/*.sh \
         "$WS2"/15_marker_trees/mcrA/*.log; do
    [ -f "$f" ] || continue
    echo "   --- ${f#$WS2/} ---"
    grep -nE "esearch|efetch|eutils|entrez|datasets|wget|curl|blastp|hmmsearch|\
term=|query=|MCR|ACR|alkyl|methyl-coenzyme" "$f" 2>/dev/null | head -20 | sed 's/^/      /'
done

echo
echo "=============================================================="
echo "C. is any MCR-labelled reference non-archaeal"
echo "=============================================================="
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections
WS2, OUT = sys.argv[1], sys.argv[2]
D = os.path.join(WS2, "15_marker_trees", "mcrA")

HEAD = {}
p = os.path.join(D, "mcrA_all.faa")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        if ln.startswith(">"):
            h = ln[1:].rstrip("\n")
            HEAD[h.split()[0]] = h
print("   %d defline(s) read" % len(HEAD))

# Archaeal methanogen/ANME genus stems. This list is a FILTER FOR ATTENTION,
# not a classifier: anything it does not recognise is printed for inspection
# rather than being called bacterial.
ARC = ("methano", "methermicoccus", "archaeoglobus", "ferroglobus",
       "candidatus methano", "anme", "syntrophoarchaeum", "argoarchaeum",
       "methanoliparum", "bathyarchae", "verstraetearchae", "korarchae",
       "thermoplasma", "halo", "nitrososphaer", "sulfolobus", "thermococc",
       "pyrococc", "methanopyrus", "candidatus syntroph", "candidatus bathy",
       "candidatus methanophagales", "candidatus alkanophaga")

rows = []
for key, h in sorted(HEAD.items()):
    panel = key.split("|")[0] if "|" in key else "?"
    if panel not in ("MCR", "ACR"):
        continue
    m = re.search(r"\[([^\]]+)\]", h)
    org = m.group(1) if m else ""
    low = org.lower()
    known = any(low.startswith(a) or a in low for a in ARC)
    rows.append((panel, key, org, known))

for panel in ("MCR", "ACR"):
    sub = [r for r in rows if r[0] == panel]
    unk = [r for r in sub if not r[3]]
    print()
    print("   %s panel: %d sequence(s), %d with an organism this step does not"
          " recognise as archaeal" % (panel, len(sub), len(unk)))
    if not sub:
        continue
    gen = collections.Counter((r[2].split()[0] if r[2] else "(no organism)")
                              for r in sub)
    print("      genera present: %s"
          % ", ".join("%s %d" % (k, v) for k, v in gen.most_common()))
    for _, key, org, known in unk:
        print("      >> %-28s %s" % (key.split("|")[-1], org or "(no [organism] in defline)"))
        print("         raw: %s" % HEAD[key][:150])

with open(os.path.join(OUT, "refset_organisms.tsv"), "w") as fh:
    fh.write("panel\ttip\torganism\trecognised_archaeal\n")
    for panel, key, org, known in rows:
        fh.write("%s\t%s\t%s\t%s\n" % (panel, key, org, "yes" if known else "no"))
print()
print("   -> 91_mcr_refset_audit/refset_organisms.tsv")
PYEOF

echo
echo "=============================================================="
echo "D. the ten nearest neighbours of each query, from the DISTANCE MATRIX"
echo "=============================================================="
echo "   mcrA.mldist is the maximum-likelihood distance matrix, computed before"
echo "   the tree and independent of its topology. The report's CONFIRMED"
echo "   verdicts came from it. Recomputed here so the claim can be stated as"
echo "   named organisms rather than as a panel tag."
echo
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re
WS2, OUT = sys.argv[1], sys.argv[2]
D = os.path.join(WS2, "15_marker_trees", "mcrA")
md = os.path.join(D, "mcrA.mldist")
if not os.path.exists(md):
    print("   MISSING %s" % md); sys.exit(0)

lines = [l.rstrip("\n") for l in open(md, errors="replace") if l.strip()]
try:
    n = int(lines[0].split()[0])
except (ValueError, IndexError):
    print("   first line is not a count: %r" % lines[0][:80]); sys.exit(0)
names, M = [], []
for l in lines[1:]:
    f = l.split()
    if len(f) < 2:
        continue
    names.append(f[0])
    M.append([float(x) for x in f[1:]])
print("   %d x %d matrix (header says %d)" % (len(names), len(M[0]) if M else 0, n))
if len(names) != n:
    print("   !! row count does not match the header; not using this file")
    sys.exit(0)

HEAD = {}
p = os.path.join(D, "mcrA_all.faa")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        if ln.startswith(">"):
            h = ln[1:].rstrip("\n")
            HEAD[h.split()[0]] = h


def org(t):
    h = HEAD.get(t, "")
    m = re.search(r"\[([^\]]+)\]", h)
    return m.group(1) if m else t.split("|")[-1]


def panel(t):
    return t.split("|")[0] if "|" in t else "?"


idx = {nm: i for i, nm in enumerate(names)}
qs = [nm for nm in names if panel(nm) == "QUERY"]
print("   %d query sequence(s) in the matrix" % len(qs))
print()
for q in qs:
    i = idx[q]
    d = sorted(((M[i][j], names[j]) for j in range(len(names)) if j != i))[:10]
    comp = {}
    for _, nm in d:
        comp[panel(nm)] = comp.get(panel(nm), 0) + 1
    f = q.split("|")
    print("   %s  %s" % (f[1] if len(f) > 1 else q, f[3] if len(f) > 3 else ""))
    print("      ten nearest: %s"
          % ", ".join("%s %d" % (k, v) for k, v in sorted(comp.items())))
    for dist, nm in d:
        print("        %7.4f  %-6s %-30s %s"
              % (dist, panel(nm), org(nm)[:30], nm.split("|")[-1]))
    print()
PYEOF

echo
echo "=============================================================="
echo "E. the same two questions for the mttB tree"
echo "=============================================================="
echo "   Built the same way, and fig15 would inherit whatever is wrong here."
echo
for f in "$WS2"/15_marker_trees/mttB/mttB_all.faa; do
    [ -f "$f" ] || { echo "   MISSING $f"; continue; }
    echo "   ${f#$WS2/}: $(grep -c '^>' "$f") sequence(s); first 30 deflines:"
    grep "^>" "$f" | head -30 | sed 's/^/      /'
done
echo
echo "   mttB verdicts:"
[ -f "$WS2/15_marker_trees/mttB/mttB_verdicts.tsv" ] && \
    sed 's/^/      /' "$WS2/15_marker_trees/mttB/mttB_verdicts.tsv"

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   If the raw deflines confirm non-archaeal organisms in the MCR panel,"
echo "   then the panel tag 'MCR' means 'matched a protein name', section 5's"
echo "   'ten nearest neighbours all MCR' has to be restated as named"
echo "   organisms, and fig14's caption says so. Section D already gives the"
echo "   replacement wording."
echo
echo "   If instead the deflines are clean and my extraction mangled them, the"
echo "   only thing wrong is fig14's label code, and the note comes out of the"
echo "   caption."
echo
echo "   Either way the placement itself is unaffected: every query's nearest"
echo "   relative is 11 to 102 times closer than the nearest ACR sequence, and"
echo "   section D prints the organisms behind that number."
