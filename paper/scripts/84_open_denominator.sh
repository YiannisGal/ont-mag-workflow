#!/usr/bin/env bash
# 84_open_denominator.sh -- the ceiling compared an open reference set against
# a closed one.
#
# WHAT STEP 83 SHOWED. Both controls passed: mcrA is 100% archaeal in both
# samples, and every AAM* accession is Methanosarcina acetivorans C2A, which is
# what Volume XV called by hand. So the carrier split is sound -- and it kills
# the bacterial explanation, because only 65 of 283 sediment mtbB reads are
# bacterial and the archaeal-only ratio is still 9.08x.
#
# THE ANSWER IS IN THE ORGANISM NAMES, NOT THE COUNTS.
#
#   mcrA references   9, ALL of them our own contigs: Methanosarcina,
#                     Methanofollis, MTP4. A CLOSED set -- it can only detect
#                     the methanogens already in this catalogue.
#
#   mtbB references   51, including Candidatus Methanomassiliicoccus
#                     intestinalis (110 sediment reads, 39% of the label),
#                     Methanomassiliicoccus sp. (17) and Methanonatronarchaeia
#                     archaeon (18). An OPEN set spanning methanogen orders
#                     with NO genome in this catalogue.
#
# A ratio between an open numerator and a closed denominator is not a ceiling
# test. Sediment reads from a Methanomassiliicoccales can match that order's
# mtbB, because it is in the panel, while their mcrA has nothing to match --
# our nine mcrA references are all Methanosarcinales and Methanomicrobiales.
# The excess is where the reference sets differ, not where the biology is.
#
# AND IT PREDICTS SOMETHING CHECKABLE. If those 110 reads are real, a
# methyl-reducing methanogen order is present in Etoliko sediment and absent
# from the 29 MAGs. That is a claim about the catalogue's coverage, and this
# project already has four independent taxonomic profiles of the same sediment
# that can be asked about it WITHOUT a new search.
#
# So: section A quantifies the asymmetry, section B asks the profilers, and
# only section C says whether a new search is warranted.
set -u
WS2="${WS2}"
P83="$WS2/83_who_really"
OUT="$WS2/84_open_denominator"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. taxonomic breadth of each label's reference set"
echo "=============================================================="
"$PY" - "$WS2/82_mtbB_who/reference_taxonomy.tsv" "$P83" <<'PYEOF'
import sys, os, csv, collections
src, P83 = sys.argv[1], sys.argv[2]

own, ext = {}, {}
p = os.path.join(P83, "own_tax.tsv")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 5:
            own[f[0]] = f[4] or "?"          # genus, from our own MAGs
p = os.path.join(P83, "ext_tax.tsv")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) == 2 and f[1]:
            ext[f[0]] = f[1]

# the genera this catalogue actually contains
OURS = {"Methanosarcina", "Methanofollis", "MTP4", "Methanogenium",
        "Methanobacterium_C", "Methanobacterium"}


def genus(ref):
    r = ref.split("|")[-1]
    if r in own:
        return own[r], True
    o = ext.get(r, "")
    if not o:
        return "(unresolved)", False
    g = o.split()[1] if o.startswith("Candidatus ") else o.split()[0]
    return g, g in OURS


rows = collections.defaultdict(lambda: collections.defaultdict(int))
for r in csv.DictReader(open(src), delimiter="\t"):
    if r["sample"] != "SED":
        continue
    g, mine = genus(r["reference"])
    rows[r["label"]][(g, mine)] += int(r["reads"])

for lab in ("mcrA", "mtmB", "DECOYmtaB", "DECOYmtbB"):
    d = rows.get(lab)
    if not d:
        continue
    tot = sum(d.values())
    inside = sum(n for (g, mine), n in d.items() if mine)
    print("   %s -- %d sediment reads, %d genus/genera"
          % (lab, tot, len({g for g, _ in d})))
    for (g, mine), n in sorted(d.items(), key=lambda kv: -kv[1])[:8]:
        print("      %5d  %-44s %s"
              % (n, g[:44], "in this catalogue" if mine else "NOT in it"))
    print("      reads on genera present in the catalogue: %d of %d (%.0f%%)"
          % (inside, tot, 100.0 * inside / tot if tot else 0))
    print()
print("   A label whose reads sit on genera the catalogue does not contain is")
print("   detecting organisms the catalogue missed. A denominator built only")
print("   from the catalogue's own contigs cannot see those organisms at all,")
print("   so the ratio between them measures the gap between two reference")
print("   sets, not a biological ceiling.")
PYEOF
echo

echo "=============================================================="
echo "B. do the sediment profilers see these lineages independently?"
echo "=============================================================="
# Four profiles of the same sediment already exist. If a methyl-reducing
# methanogen order is really there, at least one should say so -- and if none
# does, the 110 reads need a different explanation.
for pat in Methanomassiliicoccus Methanomassiliicoccales Thermoplasmatota \
           Methanonatronarchaeia Methanofastidiosa Verstraetearchaeota; do
  echo "--- $pat"
  hit=0
  for f in "$WS2"/22_singlem/SED_profile.tsv \
           "$WS2"/24_deep_taxonomy/singlem_ranks.tsv \
           "$WS2"/20_amplicon/*SED*.tsv "$WS2"/21_kaiju/*SED*.tsv \
           "$WS2"/19_silva/*SED*.tsv; do
    [ -s "$f" ] || continue
    n=$(grep -c -i "$pat" "$f" 2>/dev/null || true)
    if [ "${n:-0}" -gt 0 ]; then
      hit=1
      echo "     ${f#$WS2/}: $n line(s)"
      grep -i -m3 "$pat" "$f" | cut -c1-150 | sed 's/^/        /'
    fi
  done
  [ "$hit" -eq 0 ] && echo "     not found in any sediment profile on disk"
done
echo
echo "   (a profiler that never saw the lineage is not evidence of absence if"
echo "    its database lacks it -- but SingleM and GTDB-Tk both carry"
echo "    Thermoplasmatota, so a zero there is informative)"
echo

echo "=============================================================="
echo "C. is a new search warranted?"
echo "=============================================================="
echo "   If B found the lineage in an independent sediment profile:"
echo "     the 110 reads are corroborated, mtbB's excess is explained as an"
echo "     OPEN-vs-CLOSED reference asymmetry, and the reportable finding is"
echo "     about catalogue coverage -- a methyl-reducing methanogen order in"
echo "     the sediment with no MAG. No new search needed for that sentence."
echo
echo "   If B found nothing:"
echo "     the reads rest on one reference and one method. Then it is worth"
echo "     adding divergent-lineage mcrA references to the denominator and"
echo "     re-running step 61's OWN search command -- extracted from its"
echo "     script, not rewritten -- so numerator and denominator finally span"
echo "     the same taxa. That step is written from B's output."
echo
echo "   Either way mtbB stays out of the chapter as a per-gene abundance."
echo "   What changes is the reason given, and the reason is now a real one."
