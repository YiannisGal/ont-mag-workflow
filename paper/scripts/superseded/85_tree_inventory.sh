#!/usr/bin/env bash
# 85_tree_inventory.sh -- what phylogenies exist, and are they figure-ready?
#
# WHY. Three claims in the findings report rest on trees that no figure shows:
#   mcrA placement       "CONFIRMED for 9 proteins across 7 genomes"
#   MttB placement       "confirmed in MTP4 and both Methanosarcina, nowhere else"
#   DsrAB placement      "reductive-type, not reverse-type" -- the competition
#                        premise of the whole chapter depends on this one
# and a catalogue phylogenomic tree is standard in a MAG chapter and absent.
#
# Reads only. Nothing is drawn until this says what is actually on disk: a
# figure script written against an assumed Newick file would fail on the first
# real one, and a tree drawn without knowing whether its labels are support
# values or branch lengths is a picture, not evidence.
set -u
WS2="${WS2}"
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. every tree file in the workspace"
echo "=============================================================="
find "$WS2" -type f \( -name "*.treefile" -o -name "*.tree" -o -name "*.nwk" \
     -o -name "*.newick" -o -name "*.contree" -o -name "*.nex" \) \
     ! -path "*/envs/*" ! -path "*/.pixi/*" 2>/dev/null \
  | sort | while read -r f; do
      printf '   %8s  %s\n' "$(du -h "$f" | cut -f1)" "${f#$WS2/}"
    done
echo

echo "=============================================================="
echo "B. can each one be parsed, and what is in it?"
echo "=============================================================="
# A tree is figure-ready only if we can say how many tips it has, whether the
# internal labels are support values, and which tips are OURS. All three decide
# how the panel is drawn.
find "$WS2" -type f \( -name "*.treefile" -o -name "*.tree" -o -name "*.contree" \) \
     ! -path "*/envs/*" ! -path "*/.pixi/*" 2>/dev/null | sort > /tmp/trees.txt
"$PY" - "$WS2" < /tmp/trees.txt <<'PYEOF'
import sys, re, os
WS2 = sys.argv[1]
paths = [l.strip() for l in sys.stdin if l.strip()]
if not paths:
    print("   no parseable tree files found")
    sys.exit(0)

TOK = re.compile(r"[(),;]|[^(),:;]+")


def tips_and_internals(nw):
    """minimal Newick walk -- no dependency, because the figure scripts will
    have to do the same and adding a parser library to this machine's envs is
    exactly what the project rules forbid"""
    tips, internals, depth, maxdepth = [], [], 0, 0
    i, n = 0, len(nw)
    buf = ""
    prev = ""
    while i < n:
        c = nw[i]
        if c == "(":
            depth += 1; maxdepth = max(maxdepth, depth)
            if buf.strip():
                pass
            buf = ""; prev = "("
        elif c in ",);":
            lab = buf.strip()
            if lab:
                name = lab.split(":")[0].strip().strip("'\"")
                if prev == ")":
                    internals.append(name)
                else:
                    tips.append(name)
            buf = ""
            if c == ")":
                depth -= 1; prev = ")"
            else:
                prev = c
        else:
            buf += c
        i += 1
    return tips, internals, maxdepth


for p in paths:
    try:
        nw = open(p, errors="replace").read().strip()
    except OSError as e:
        print("   %s -- unreadable (%s)" % (p[len(WS2) + 1:], e)); continue
    if not nw.endswith(";"):
        nw += ";"
    try:
        tips, internals, d = tips_and_internals(nw)
    except Exception as e:
        print("   %s -- parse failed (%s)" % (p[len(WS2) + 1:], e)); continue

    sup = [x for x in internals if x]
    numeric = [x for x in sup if re.match(r"^[\d.]+(/[\d.]+)*$", x)]
    ours = [t for t in tips if re.search(r"(ETO_|LIQ_|SED_|MAG_|contig_)", t)]
    print("   %s" % p[len(WS2) + 1:])
    print("      %d tips, %d labelled internal node(s), max depth %d"
          % (len(tips), len(sup), d))
    print("      internal labels numeric (support values): %d of %d%s"
          % (len(numeric), len(sup),
             "   e.g. " + sup[0] if sup else ""))
    print("      tips that look like OURS: %d%s"
          % (len(ours), ("   e.g. " + ours[0]) if ours else ""))
    if tips:
        print("      example reference tip: %s"
              % next((t for t in tips if t not in ours), "(none)"))
PYEOF
echo

echo "=============================================================="
echo "C. the verdict tables the trees were built to produce"
echo "=============================================================="
for f in "$WS2"/15_marker_trees/*/*verdicts*.tsv "$WS2"/58_dsr_finish/*.tsv \
         "$WS2"/12_phylogeny/derep_clusters.tsv; do
  [ -s "$f" ] || continue
  echo "--- ${f#$WS2/}  ($(( $(wc -l < "$f") - 1 )) rows)"
  head -4 "$f" | cut -c1-170 | sed 's/^/     /'
done
echo

echo "=============================================================="
echo "D. alignments and models, for the caption"
echo "=============================================================="
# A tree figure's caption must say what was aligned, with what, under which
# model, and how support was computed. Those are in the IQ-TREE log, not in
# anyone's memory.
for f in "$WS2"/12_phylogeny/trees/*.iqtree "$WS2"/12_phylogeny/ref_trees/*.iqtree \
         "$WS2"/15_marker_trees/*/*.iqtree; do
  [ -s "$f" ] || continue
  echo "--- ${f#$WS2/}"
  grep -m1 -E "Model of substitution|Best-fit model" "$f" | sed 's/^/     /'
  grep -m1 -E "Input data:|number of taxa" "$f" | sed 's/^/     /'
  grep -m1 -E "ultrafast bootstrap|Replicates" "$f" | sed 's/^/     /'
done
for f in "$WS2"/15_marker_trees/*/*.log "$WS2"/58_dsr_finish/*.log; do
  [ -s "$f" ] || continue
  echo "--- ${f#$WS2/}"
  grep -m2 -i -E "fasttree|mafft|trimal|iqtree|model" "$f" | cut -c1-140 | sed 's/^/     /'
done
echo

echo "=============================================================="
echo "WHAT I NEED FROM THIS"
echo "=============================================================="
echo "   B tells me which trees can be drawn and whether their internal labels"
echo "   are support values (drawn) or empty (not drawn, and the caption must"
echo "   say support was not computed rather than leaving the reader to assume)."
echo "   C tells me which tips to highlight. D gives the caption its method"
echo "   sentence. Then fig12 (catalogue phylogenomics), fig13 (DsrAB,"
echo "   reductive vs reverse) and fig14 (mcrA) can be written against reality."
