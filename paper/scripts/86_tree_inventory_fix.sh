#!/usr/bin/env bash
# 86_tree_inventory_fix.sh -- parse the trees, and check who is in them.
#
# STEP 85 SECTION B WAS A SHELL BUG, NOT AN ABSENCE. It did
#     "$PY" - "$WS2" < /tmp/trees.txt <<'PYEOF'
# where the heredoc supplies the PROGRAM on stdin. The last redirect wins, so
# the file list was never read and the step announced "no parseable tree files
# found" while section A was listing twenty-six of them. The list has to be an
# argument, because stdin is already taken.
#
# AND SECTION D RAISED SOMETHING THAT MUST BE SETTLED BEFORE ANY TREE IS DRAWN.
#
#     12_phylogeny/trees/ar53.iqtree     8 sequences
#     12_phylogeny/trees/bac120.iqtree  28 sequences
#
# The catalogue is 21 bacteria and 8 archaea. Eight in ar53 is exactly right.
# TWENTY-EIGHT in bac120 is not 21, and 8 + 28 = 36 against a catalogue of 29,
# so at least seven genomes are in BOTH trees. Either bac120 was run on the
# 28 dereplicated genomes irrespective of domain -- in which case archaea were
# aligned against bacterial marker profiles, and a tool trained on one domain
# applied to another fails silently and in one direction (rule 66) -- or the
# count means something else.
#
# A figure cannot be drawn from a tree whose membership is unexplained, so
# section C joins every tip to its domain in MASTER.tsv and counts.
set -u
WS2="${WS2}"
OUT="$WS2/86_tree_inventory"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

find "$WS2" -type f \( -name "*.treefile" -o -name "*.tree" -o -name "*.contree" \) \
     ! -path "*/envs/*" ! -path "*/.pixi/*" 2>/dev/null | sort > "$OUT/trees.txt"
echo "   $(wc -l < "$OUT/trees.txt") tree file(s) to parse"
echo

"$PY" - "$WS2" "$OUT/trees.txt" "$OUT" <<'PYEOF'
import sys, os, re, csv, collections
WS2, LIST, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
paths = [l.strip() for l in open(LIST) if l.strip()]


def parse(nw):
    """Minimal Newick reader. Written here rather than imported: adding a
    parser library to this machine's environments is what the project rules
    forbid, and the figure scripts will need the same thirty lines."""
    nw = nw.strip()
    if not nw.endswith(";"):
        nw += ";"
    tips, internal = [], []
    buf, prev, depth, maxd = "", "", 0, 0
    for c in nw:
        if c == "(":
            depth += 1; maxd = max(maxd, depth); buf = ""; prev = "("
        elif c in ",);":
            lab = buf.strip()
            if lab:
                name = lab.split(":")[0].strip().strip("'\"")
                if name:
                    (internal if prev == ")" else tips).append(name)
            buf = ""
            if c == ")":
                depth -= 1; prev = ")"
            else:
                prev = c
        else:
            buf += c
    return tips, internal, maxd


# ---- domain per genome, from MASTER
dom, gen = {}, {}
m = os.path.join(WS2, "47_master", "MASTER.tsv")
if os.path.exists(m):
    for r in csv.DictReader(open(m, errors="replace"), delimiter="\t"):
        d = (r.get("domain") or "?").strip()
        g = (r.get("genus") or "").strip()
        for k in ("mag_id", "working_name", "qc_bin_name", "bin"):
            v = (r.get(k) or "").strip()
            if v:
                dom[v] = d; gen[v] = g
    print("   MASTER.tsv: %d genome name(s) with a domain" % len(dom))
else:
    print("   !! MASTER.tsv missing -- domain composition cannot be checked")
print()

NUM = re.compile(r"^[\d.]+(/[\d.]+)*$")
rows = []
for p in paths:
    try:
        tips, internal, d = parse(open(p, errors="replace").read())
    except Exception as e:
        print("   %-52s PARSE FAILED (%s)" % (p[len(WS2) + 1:], e)); continue
    sup = [x for x in internal if x]
    numeric = sum(1 for x in sup if NUM.match(x))
    rows.append((p, tips, sup, numeric, d))

print("=" * 66)
print("A. every tree, parsed")
print("=" * 66)
print("   %-46s %6s %8s %8s" % ("file", "tips", "labels", "numeric"))
for p, tips, sup, numeric, d in rows:
    print("   %-46s %6d %8d %8d"
          % (p[len(WS2) + 1:][-46:], len(tips), len(sup), numeric))
print()
print("   'numeric' = internal labels that look like support values. A tree")
print("   with 0 of them has no support to draw, and its caption must say so")
print("   rather than let the reader assume the topology is supported.")
print()

print("=" * 66)
print("B. who is in the catalogue trees")
print("=" * 66)


def own(t):
    """our genome names appear as bare bin names or ETO_ ids"""
    for k in (t, t.split("|")[0], t.split("~")[0]):
        if k in dom:
            return k
    return None


for p, tips, sup, numeric, d in rows:
    rel = p[len(WS2) + 1:]
    if not rel.startswith("12_phylogeny"):
        continue
    c = collections.Counter()
    unknown = []
    for t in tips:
        k = own(t)
        if k:
            c[dom[k]] += 1
        else:
            c["reference/unmatched"] += 1
            unknown.append(t)
    print("   %s  (%d tips)" % (rel, len(tips)))
    for k, v in c.most_common():
        print("      %-24s %d" % (k, v))
    if unknown:
        print("      unmatched examples: %s" % ", ".join(unknown[:4]))
    if rel.endswith(("bac120.treefile", "bac120.contree")) and c.get("Archaea"):
        print("      >> %d ARCHAEAL genome(s) in a bac120 tree. Bacterial marker"
              % c["Archaea"])
        print("         profiles applied to archaea fail silently and in one")
        print("         direction. This tree cannot be published as it stands.")
    print()

print("=" * 66)
print("C. the dsr tree -- the one the chapter's premise rests on")
print("=" * 66)
for p, tips, sup, numeric, d in rows:
    if "58_dsr_finish" not in p:
        continue
    print("   %s: %d tips, %d support labels" % (p[len(WS2) + 1:], len(tips), numeric))
    q = [t for t in tips if t.startswith("QUERY") or "contig_" in t]
    print("      our sequences (%d): %s" % (len(q), ", ".join(q[:6])))
    # the reverse-type references were put in on purpose; the figure must show
    # them as a named clade or the "not rDsr" claim is invisible
    REV = ("Allochromatium", "Thiobacillus", "Chlorobium", "Desulfurivibrio",
           "Chlorobaculum", "Halorhodospira", "Thermochromatium")
    rv = [t for t in tips if any(r.lower() in t.lower() for r in REV)]
    print("      reverse-type (rDsr) references present (%d): %s"
          % (len(rv), ", ".join(rv[:6]) if rv else "NONE FOUND"))
    if not rv:
        print("      >> the 'neither goes near rDsr' claim cannot be drawn from")
        print("         this tree if the rDsr references are not tips in it.")
    with open(os.path.join(OUT, "dsr_tips.txt"), "w") as fh:
        fh.write("\n".join(tips) + "\n")
    print("      all tips written to 86_tree_inventory/dsr_tips.txt")
print()

print("=" * 66)
print("D. marker trees -- how many of the tips are ours")
print("=" * 66)
for p, tips, sup, numeric, d in rows:
    rel = p[len(WS2) + 1:]
    if not rel.startswith("15_marker_trees") or not rel.endswith(".treefile"):
        continue
    ours = [t for t in tips if t.startswith("ETO_") or "|ETO" in t]
    pan = collections.Counter(t.split("|")[0] for t in tips if not t.startswith("ETO_"))
    print("   %s: %d tips, %d ours" % (rel, len(tips), len(ours)))
    print("      reference panels: %s"
          % ", ".join("%s %d" % (k, v) for k, v in pan.most_common(6)))
PYEOF
echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   B settles whether the catalogue phylogenomic tree is drawable at all."
echo "   If archaea are sitting in the bac120 tree, fig12 must be built from"
echo "   ar53 and bac120 as TWO panels over the genomes that belong in each,"
echo "   and step 12 has a real error to record -- not a figure to draw around."
echo "   C settles whether fig13 can show the reductive-vs-reverse contrast,"
echo "   which is the only reason that panel is worth a page."
