#!/usr/bin/env bash
# 89_fig_inputs.sh -- the inputs for the two phylogenetic figures that are still
# missing, printed in full before either script is written.
#
# WHY THIS STEP EXISTS AND WHY IT IS ONLY PRINTING.
#
# fig13 was written twice. The first version drew two genes as one tree, showed
# no support values under a caption that said support was shown, and named two
# reference genera that are not tips in the tree at all. Every one of those
# three failures was a guess about a file I had not looked at. The inventory in
# steps 85-86 is what caught them, and it caught them only because it printed.
#
# So before fig14 (mcrA) and fig12 (catalogue phylogenomics) are written, this
# step prints the things those scripts will have to assume:
#
#   * the exact tip names of the mcrA tree, grouped by reference panel, so the
#     script can tell OUR nine proteins from the 60 MCR and 36 ACR references
#     by a rule taken from the data rather than from a naming habit;
#   * whether the file is .treefile or .contree, and whether the two differ;
#   * the substitution model, the number of ALIGNMENT SITES and the number of
#     bootstrap replicates, which belong in the caption (rule 130);
#   * the nine nominated mcrA proteins and which genome each came from;
#   * the curated KOfam thresholds and the per-genome scores, because panel (a)
#     of fig14 has to show why SEVEN of eight archaea are in the tree and the
#     eighth is not -- and 41_env_checks/archaea_kofam_markers.tsv must not be
#     used for that, it is a presence screen with no threshold applied at all;
#   * the same for the ar53/bac120 reference trees, since fig12 is next.
#
# NOTHING IS COMPUTED HERE AND NOTHING IS WRITTEN OUTSIDE ITS OWN DIRECTORY.
#
# Rule 125 is the one to watch in this step: a directory is not one kind of
# file. Every loop below names the extension it wants.
set -u
WS2="${WS2}"
OUT="$WS2/89_fig_inputs"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. what is actually in the marker-tree and phylogeny directories"
echo "=============================================================="
for d in "$WS2"/15_marker_trees/mcrA "$WS2"/15_marker_trees/mttB \
         "$WS2"/12_phylogeny/trees "$WS2"/12_phylogeny/ref_trees; do
    if [ ! -d "$d" ]; then echo "   MISSING $d"; continue; fi
    echo "   ${d#$WS2/}"
    ls -la "$d" | awk 'NR>1 {printf "      %10s  %s\n", $5, $9}'
    echo
done

echo "=============================================================="
echo "B. the mcrA tree, tip by tip"
echo "=============================================================="
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections, hashlib
WS2, OUT = sys.argv[1], sys.argv[2]


def parse_tips(nw):
    """tips only; the same minimal reader as step 86"""
    nw = nw.strip()
    if not nw.endswith(";"):
        nw += ";"
    tips, buf, prev = [], "", ""
    for c in nw:
        if c == "(":
            buf, prev = "", "("
        elif c in ",);":
            lab = buf.strip()
            if lab and prev != ")":
                nm = lab.split(":")[0].strip().strip("'\"")
                if nm:
                    tips.append(nm)
            buf = ""
            prev = ")" if c == ")" else c
        else:
            buf += c
    return tips


def show(d, stems):
    """print every tree file under d, and say whether they are the same tree"""
    if not os.path.isdir(d):
        print("   MISSING %s" % d)
        return {}
    seen = {}
    for fn in sorted(os.listdir(d)):
        if not fn.endswith((".treefile", ".contree", ".tree", ".nwk")):
            continue
        p = os.path.join(d, fn)
        try:
            txt = open(p, errors="replace").read()
        except OSError as e:
            print("   %-28s UNREADABLE (%s)" % (fn, e))
            continue
        tips = parse_tips(txt)
        seen[fn] = tips
        print("   %-28s %4d tips   md5 %s"
              % (fn, len(tips), hashlib.md5(txt.encode()).hexdigest()[:7]))
    if len(seen) > 1:
        ks = sorted(seen)
        for i in range(len(ks) - 1):
            a, b = ks[i], ks[i + 1]
            same = set(seen[a]) == set(seen[b])
            print("      %s vs %s: same tip set = %s" % (a, b, same))
    return seen


for marker in ("mcrA", "mttB"):
    d = os.path.join(WS2, "15_marker_trees", marker)
    print("   --- %s ---" % marker)
    seen = show(d, None)
    if not seen:
        print()
        continue
    # the panel prefix is whatever precedes the first '|'. Print EVERY distinct
    # one with its count and three examples, so the figure script's rule for
    # "ours" can be read off this output instead of being invented.
    fn = sorted(seen, key=lambda k: (not k.endswith(".contree"), k))[0]
    tips = seen[fn]
    pan = collections.OrderedDict()
    for t in tips:
        key = t.split("|")[0] if "|" in t else "(no pipe)"
        pan.setdefault(key, []).append(t)
    print("   tip-name prefixes in %s:" % fn)
    for k, v in sorted(pan.items(), key=lambda kv: -len(kv[1])):
        print("      %-16s %4d   e.g. %s" % (k, len(v), " | ".join(v[:3])))
    with open(os.path.join(OUT, "%s_tips.txt" % marker), "w") as fh:
        fh.write("\n".join(tips) + "\n")
    print("   all %d tip names -> 89_fig_inputs/%s_tips.txt" % (len(tips), marker))
    # and print in full the ones that are not in the biggest two panels --
    # ours should be among them, and if they are not, I need to see that now
    big = [k for k, v in sorted(pan.items(), key=lambda kv: -len(kv[1]))[:2]]
    rest = [t for t in tips if (t.split("|")[0] if "|" in t else "(no pipe)") not in big]
    print("   tips outside the two largest panels (%d), in full:" % len(rest))
    for t in rest:
        print("      %s" % t)
    print()

print("   --- the catalogue trees ---")
for sub in ("trees", "ref_trees"):
    d = os.path.join(WS2, "12_phylogeny", sub)
    print("   12_phylogeny/%s" % sub)
    seen = show(d, None)
    for fn, tips in seen.items():
        with open(os.path.join(OUT, "%s_%s_tips.txt" % (sub, fn)), "w") as fh:
            fh.write("\n".join(tips) + "\n")
        if "ref" in fn or sub == "ref_trees":
            ours = [t for t in tips if t.startswith("ETO_") or "MAG" in t
                    or t.startswith(("das_", "bin", "metabat", "vamb", "semibin"))]
            print("      %s: %d tip(s) look like ours; first few: %s"
                  % (fn, len(ours), ", ".join(ours[:5])))
    print()
PYEOF

echo "=============================================================="
echo "C. models, alignment length, replicates -- for the captions"
echo "=============================================================="
for f in "$WS2"/15_marker_trees/mcrA/*.iqtree "$WS2"/15_marker_trees/mttB/*.iqtree \
         "$WS2"/12_phylogeny/ref_trees/*.iqtree "$WS2"/12_phylogeny/trees/*.iqtree; do
    [ -f "$f" ] || continue
    echo "   ${f#$WS2/}"
    grep -E "Input data:|Model of substitution:|Number of (constant|distinct|parsimony|invariant)|ultrafast bootstrap|replicate|Best-fit model|Log-likelihood of the tree" "$f" \
        | head -12 | sed 's/^/      /'
    echo
done

echo "=============================================================="
echo "D. the nine nominated mcrA proteins, and where each came from"
echo "=============================================================="
CAND="$WS2/11_kofam/verify/mcrA_K00399_candidates.faa"
if [ -f "$CAND" ]; then
    echo "   $CAND"
    grep "^>" "$CAND" | sed 's/^/      /'
    echo "      ---"
    echo "      $(grep -c '^>' "$CAND") sequence(s)"
else
    echo "   MISSING $CAND"
    echo "   looking for anything similar:"
    find "$WS2/11_kofam" -name "*K00399*" -o -name "*mcrA*" 2>/dev/null | sed 's/^/      /'
fi
echo
echo "   and the verdict table the report quotes:"
for f in "$WS2"/15_marker_trees/mcrA/mcrA_verdicts.tsv "$WS2"/15_marker_trees/mcrA/*verdict*; do
    [ -f "$f" ] || continue
    echo "   ${f#$WS2/}"
    cat "$f" | sed 's/^/      /'
    break
done

echo
echo "=============================================================="
echo "E. the curated thresholds, read from ko_list"
echo "=============================================================="
KL="$WS2/36_gapclose/ko_list"
if [ -f "$KL" ]; then
    head -1 "$KL" | sed 's/^/      /'
    grep -E "^(K00399|K00401|K00402|K14083|K16176|K16178|K04480)\b" "$KL" | sed 's/^/      /'
else
    echo "   MISSING $KL"
fi

echo
echo "=============================================================="
echo "F. per-genome K00399 scores -- EVERY row, above threshold or not"
echo "=============================================================="
echo "   41_env_checks/archaea_kofam_markers.tsv is NOT used here. It applies"
echo "   no score threshold at all (41f line 124) and the report already"
echo "   records that its 'K00399 = 4' is a raw hit count, not a gene count."
echo
"$PY" - "$WS2" <<'PYEOF'
import sys, os, re, collections
WS2 = sys.argv[1]

# --- thresholds
THR = {}
kl = os.path.join(WS2, "36_gapclose", "ko_list")
if os.path.exists(kl):
    for ln in open(kl, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 1 and re.match(r"^K\d{5}$", f[0]):
            try:
                THR[f[0]] = float(f[1])
            except ValueError:
                pass
print("   thresholds read: K00399 %s, K00401 %s, K00402 %s"
      % (THR.get("K00399"), THR.get("K00401"), THR.get("K00402")))

# --- the hmmsearch table. Field 0 target, field 2 query, field 5 score
#     (rule 114: identify the column before comparing its contents to anything)
tb = os.path.join(WS2, "45_kofam_full", "all_hits.tbl")
if not os.path.exists(tb):
    print("   MISSING %s" % tb)
    cand = []
    for root, dirs, files in os.walk(os.path.join(WS2, "45_kofam_full")):
        for f in files:
            if f.endswith((".tbl", ".tsv")):
                cand.append(os.path.join(root, f)[len(WS2) + 1:])
    print("   what is there instead: %s" % ", ".join(sorted(cand)[:12]))
    sys.exit(0)

print("   %s -- first two non-comment lines verbatim, so the columns are visible:"
      % tb[len(WS2) + 1:])
n = 0
for ln in open(tb, errors="replace"):
    if ln.startswith("#"):
        continue
    print("      %s" % ln.rstrip("\n")[:200])
    n += 1
    if n >= 2:
        break
print()

rows = []
for ln in open(tb, errors="replace"):
    if ln.startswith("#"):
        continue
    f = ln.split()
    if len(f) < 6:
        continue
    if f[2] not in ("K00399", "K00401", "K00402"):
        continue
    try:
        rows.append((f[0], f[2], float(f[5])))
    except ValueError:
        pass
print("   %d row(s) for the three mcr KOs" % len(rows))

# --- protein id -> genome. Print the join before believing it (rule 111).
prot = sorted({r[0] for r in rows})
print("   %d distinct protein id(s); the first five verbatim:" % len(prot))
for p in prot[:5]:
    print("      %s" % p)


def genome_of(pid):
    if "~" in pid:
        return pid.rpartition("~")[0]
    m = re.match(r"^(ETO_(?:SED|LIQ)_MAG_\d+)", pid)
    if m:
        return m.group(1)
    return None


gs = [genome_of(p) for p in prot]
ok = sum(1 for g in gs if g)
print("   %d of %d protein id(s) carry a genome name" % (ok, len(prot)))
if ok == 0:
    print("   !! the join fails. The figure cannot claim a per-genome score")
    print("      from this file and the next step has to find the mapping.")
    sys.exit(0)

best = collections.defaultdict(dict)
for pid, ko, sc in rows:
    g = genome_of(pid)
    if not g:
        continue
    cur = best[g].get(ko)
    if cur is None or sc > cur[1]:
        best[g][ko] = (pid, sc)

print()
print("   %-22s %-10s %10s %10s  %s" % ("genome", "KO", "best", "threshold", "verdict"))
for g in sorted(best):
    for ko in ("K00399", "K00401", "K00402"):
        if ko not in best[g]:
            continue
        pid, sc = best[g][ko]
        t = THR.get(ko)
        v = "--" if t is None else ("ABOVE" if sc >= t else "below")
        print("   %-22s %-10s %10.1f %10s  %s" % (g[:22], ko, sc,
                                                  ("%.2f" % t) if t else "?", v))

print()
print("   and every K00399 hit of every archaeal genome, in full (rule 112):")
arc = [g for g in best if "K00399" in best[g]]
for g in sorted(arc):
    hs = sorted([(sc, pid) for pid, ko, sc in rows
                 if ko == "K00399" and genome_of(pid) == g], reverse=True)
    print("      %s: %d hit(s)  %s" % (g, len(hs),
                                       ", ".join("%.1f" % s for s, _ in hs[:6])))
PYEOF

echo
echo "=============================================================="
echo "G. the fragmented genome, for panel (a) of fig14"
echo "=============================================================="
HC="$WS2/54_mcr_fragmentation/hits_with_coords.tsv"
if [ -f "$HC" ]; then
    echo "   ${HC#$WS2/}  ($(wc -l < "$HC") lines)"
    head -1 "$HC" | sed 's/^/      /'
    grep -E "K00399|K00401|K00402" "$HC" | sed 's/^/      /'
else
    echo "   MISSING $HC"
    ls -la "$WS2/54_mcr_fragmentation" 2>/dev/null | sed 's/^/      /'
fi

echo
echo "=============================================================="
echo "H. archaeal genomes: name, genus, sediment breadth"
echo "=============================================================="
"$PY" - "$WS2" <<'PYEOF'
import sys, os, csv
WS2 = sys.argv[1]
m = os.path.join(WS2, "47_master", "MASTER.tsv")
if not os.path.exists(m):
    print("   MISSING %s" % m); sys.exit(0)
rd = list(csv.DictReader(open(m, errors="replace"), delimiter="\t"))
print("   MASTER.tsv columns: %s" % ", ".join(rd[0].keys()) if rd else "(empty)")
print()
print("   %-20s %-22s %-10s %s" % ("mag_id", "genus", "domain", "sed_breadth"))
for r in rd:
    if (r.get("domain") or "").strip() != "Archaea":
        continue
    br = ""
    for k in r:
        if "breadth" in k.lower() and ("sed" in k.lower() or "SED" in k):
            br = r[k]
            break
    print("   %-20s %-22s %-10s %s" % (r.get("mag_id", "?"), r.get("genus", "?"),
                                       r.get("domain", "?"), br or "(no column)"))
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   B decides how fig14 tells our nine sequences from the 96 references,"
echo "     and whether the ACR decoy panel is named in the tip labels at all."
echo "     If it is not, the decoy panel cannot be drawn as a named clade and"
echo "     the figure loses its control -- which is the only reason to draw it."
echo "   C gives the caption its alignment length. mcrA is 105 sequences over"
echo "     139 sites; a caption that hides that is hiding the weakness."
echo "   D+E+F decide whether panel (a) can show WHY seven of eight archaea are"
echo "     in the tree. If F's join fails, panel (a) gets no scores and the"
echo "     figure will say so rather than print numbers from memory."
