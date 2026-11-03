#!/usr/bin/env bash
# 69_auxotrophy_pairing.sh -- the corroboration test, done properly this time.
#
#   bash scripts/run.sh scripts/69_auxotrophy_pairing_XXXXXXX.sh
#
# 🔴 SECTION C OF STEP 68 CANNOT BE READ THE WAY IT INVITES YOU TO READ IT.
#
# It printed a column headed "19b auxo" and a sentence saying that a genome
# receiving many compounds while 19b calls it prototrophic is an artefact. Then it
# printed zeros for LIQ_MAG_015, LIQ_MAG_012, LIQ_MAG_020 and SED_MAG_004.
#
# Those zeros do not mean prototrophic. 19b makes a FIRM auxotrophy call only for
# genomes at least 90% complete -- that is the whole point of 19b, and rule 97's
# ancestor: a missing marker in an incomplete genome is missing evidence, not a
# missing pathway. A blank in that column means "19b declined to call", and I
# printed it next to a sentence that reads it as "19b says no". Two different
# things in one column, with no way for the reader to tell them apart.
#
#   > Rule 113. Never put "absent" and "not tested" in the same column. A blank
#   > that could mean either is worse than no column, because it invites the
#   > reader to choose the convenient one.
#
# AND MTP4 IS MISSING FROM THAT TABLE ENTIRELY. It receives nothing, so it had no
# row -- and it is the genome whose zero the whole figure turns on. The question
# that matters is precisely: does 19b think MTP4 can make its own amino acids? If
# it does, "receives nothing" is corroborated. If 19b says MTP4 is auxotrophic
# while the model says it needs nobody, then its model was gap-filled to
# self-sufficiency and the empty column is an artefact.
#
# THE RIGHT INSTRUMENT WAS ALSO IN THE OUTPUT. pathway_completeness.tsv carries a
# CONTINUOUS completeness per pathway for ALL 29 genomes, with no 90% gate. It
# needs no firm call, so it has no blanks to misread.
#
# THE TEST. Within each genome, compare the 19b completeness of the pathways whose
# products it RECEIVES against the pathways it does not receive. Paired inside the
# genome, so genome quality cancels: a poorly annotated genome looks incomplete
# everywhere and contributes nothing to the difference. If received pathways are
# systematically less complete, two unrelated methods agree about WHICH compounds,
# not merely about which genomes.
#
# HONEST LIMIT, STATED UP FRONT. 19b covers 13 amino acid pathways and 7 vitamins.
# It has no Ala, Asp, Glu, Asn, Gln or Phe, and aspartate and glutamate are two of
# the four most-transferred compounds. Those receipts are untestable here and are
# counted and reported as untestable rather than dropped.
#
# READ-ONLY except 69_pairing/. Seconds. No environment touched, no network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/69_pairing"; mkdir -p "$OUT"
exec > >(tee "$OUT/69.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "69_auxotrophy_pairing -- do the two methods agree about WHICH compounds"
echo "WS2 = $WS2"
echo "started $(date -Is)"

DET="$WS2/43_smetana/run_43q/pairs/pairs_detailed.tsv"
NAMES="$WS2/67_network/compound_names.tsv"
PC="$WS2/19b_auxotrophy/pathway_completeness.tsv"
AUX="$WS2/19b_auxotrophy/auxotrophies.tsv"
MASTER="$WS2/47_master/MASTER.tsv"
for f in "$DET" "$NAMES" "$PC" "$MASTER"; do
  [ -s "$f" ] || { no "missing $f"; echo "STEP_STATUS: FAILED"; exit 1; }
done

python3 - "$DET" "$NAMES" "$PC" "$AUX" "$MASTER" "$OUT" <<'PYEOF'
import sys, os, csv, re, collections
det, namesf, pcf, auxf, master, out = sys.argv[1:7]

# ---- the compound -> 19b pathway map, written out so it can be argued with
MAP = {
    "cysteine": "Cys", "methionine": "Met", "arginine": "Arg",
    "threonine": "Thr", "tryptophan": "Trp", "tyrosine": "Tyr",
    "glycine": "Gly", "serine": "Ser", "proline": "Pro",
    "histidine": "His", "lysine": "Lys", "leucine": "Leu",
    "isoleucine": "Ile_Val", "valine": "Ile_Val",
    "cobalamin": "B12_cobalamin", "thiamin": "B1_thiamine",
    "riboflavin": "B2_riboflavin", "pantothen": "B5_pantothenate",
    "pan": "B5_pantothenate", "pyridox": "B6_pyridoxal",
    "biotin": "B7_biotin", "folate": "B9_folate",
}
NOCOVER = ("alanine", "aspart", "glutam", "asparagine", "phenylalanine")

M = {}
for r in csv.DictReader(open(master), delimiter="\t"):
    if r.get("working_name"):
        M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")

nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(namesf), delimiter="\t")}


def bare(c):
    m = re.search(r"(cpd\d+)", c or "")
    return m.group(1) if m else c


def path_of(name):
    s = (name or "").lower()
    if s.startswith("d-"):            # D-serine is not the L pathway's product
        return None
    for k, v in MAP.items():
        if k in s:
            return v
    return None


# ---- what each genome receives
rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
gets = collections.defaultdict(set)
untest = collections.Counter()
unmapped = collections.Counter()
for r in rows:
    try:
        if float(r[low["smetana"]]) <= 0:
            continue
    except ValueError:
        continue
    g = M.get(r[low["receiver"]], r[low["receiver"]])
    name = nm.get(bare(r[low["compound"]]), "")
    p = path_of(name)
    if p:
        gets[g].add(p)
    elif any(k in (name or "").lower() for k in NOCOVER):
        untest[g] += 1
    else:
        unmapped[g] += 1

# ---- 19b completeness, continuous, all genomes, no 90% gate
PC = list(csv.DictReader(open(pcf), delimiter="\t"))
pcols = [c for c in PC[0] if c not in
         ("mag_id", "working_name", "sample", "domain", "genus", "completeness")]
AACOLS = [c for c in pcols if not c.startswith("B")]
comp = {}
gcomp = {}
for r in PC:
    g = (r.get("mag_id", "") or "").replace("ETO_", "")
    comp[g] = {}
    for c in pcols:
        try: comp[g][c] = float(r[c])
        except (TypeError, ValueError): pass
    try: gcomp[g] = float(r.get("completeness", "nan"))
    except ValueError: gcomp[g] = float("nan")

print("    19b covers %d pathway(s): %s" % (len(pcols), ", ".join(pcols)))
print("    of those, %d are amino acids." % len(AACOLS))
print("    NOT covered, and therefore untestable: %s" % ", ".join(NOCOVER))

# ================================================================ A
print()
print("================================================================")
print("A. GENOME COMPLETENESS -- so a blank is never mistaken for a zero")
print("================================================================")
firm = collections.defaultdict(set)
if os.path.exists(auxf):
    for r in csv.DictReader(open(auxf), delimiter="\t"):
        g = (r.get("mag_id", "") or "").replace("ETO_", "")
        firm[g].add(r.get("pathway", ""))
print("      %-16s %12s %-14s %8s" % ("mag", "completeness", "firm calls?", "n firm"))
for g in sorted(gets, key=lambda x: -len(gets[x])):
    c = gcomp.get(g, float("nan"))
    elig = "eligible" if c >= 90.0 else "NOT eligible"
    print("      %-16s %12.2f %-14s %8s"
          % (g, c, elig, len(firm.get(g, ())) if c >= 90.0 else "-"))
print()
print("    A genome below 90% gets no firm call at all. Step 68 printed those as")
print("    0 and invited them to be read as prototrophic. They are not.")

# ================================================================ B
print()
print("================================================================")
print("B. THE PAIRED TEST -- received pathways against the rest, within genome")
print("================================================================")
print("      %-16s %6s %8s %6s %8s %9s  %s"
      % ("mag", "n_rec", "mean_rec", "n_not", "mean_not", "diff", "direction"))
diffs = []
for g in sorted(gets, key=lambda x: -len(gets[x])):
    if g not in comp:
        print("      %-16s  no 19b row" % g); continue
    got = [comp[g][p] for p in gets[g] if p in comp[g]]
    rest = [comp[g][p] for p in comp[g] if p not in gets[g]]
    if not got or not rest:
        print("      %-16s %6d %8s %6d %8s %9s  %s"
              % (g, len(got), "-", len(rest), "-", "-", "not testable"))
        continue
    a = sum(got) / len(got); b = sum(rest) / len(rest)
    diffs.append((g, a - b))
    print("      %-16s %6d %8.3f %6d %8.3f %+9.3f  %s"
          % (g, len(got), a, len(rest), b, a - b,
             "received are LESS complete" if a < b else "received are MORE complete"))

if diffs:
    neg = sum(1 for _, d in diffs if d < 0)
    mean = sum(d for _, d in diffs) / len(diffs)
    print()
    print("    %d of %d genomes have received pathways LESS complete than the rest;"
          % (neg, len(diffs)))
    print("    mean within-genome difference %+0.3f." % mean)
    print()
    print("    A sign test on %d genomes: %d negative. Under the null of no relation"
          % (len(diffs), neg))
    print("    each genome is a coin flip, so %d/%d is %s." %
          (neg, len(diffs),
           "what chance gives you" if abs(neg - len(diffs) / 2.0) <= 1.0
           else "a lopsided result worth reporting"))
    print("    This is PAIRED INSIDE EACH GENOME, so a badly annotated genome looks")
    print("    incomplete everywhere and cancels. It tests agreement about WHICH")
    print("    compounds, which is a much harder thing to get right by accident")
    print("    than agreement about which genomes.")

# ================================================================ C
print()
print("================================================================")
print("C. THE TWO METHANOGENS, SIDE BY SIDE -- including the one with no row")
print("================================================================")
arc = [g for g in comp if g in ("LIQ_MAG_003", "LIQ_MAG_011")]
if not arc:
    print("    neither archaeon found in 19b's table")
else:
    print("      %-18s %s" % ("pathway", "  ".join("%-14s" % a for a in arc)))
    for p in pcols:
        line = "      %-18s" % p
        for a in arc:
            v = comp[a].get(p)
            mark = ""
            if v is not None and p in gets.get(a, ()):
                mark = " <-recv"
            line += "  %-14s" % (("%.2f%s" % (v, mark)) if v is not None else "-")
        print(line)
    print()
    for a in arc:
        aa = [comp[a][c] for c in AACOLS if c in comp[a]]
        print("      %s: mean amino-acid pathway completeness %.3f over %d pathway(s);"
              % (a, sum(aa) / len(aa) if aa else float("nan"), len(aa)))
        print("        genome completeness %.2f; receives %d mapped pathway(s), "
              "%d untestable transfer(s)"
              % (gcomp.get(a, float("nan")), len(gets.get(a, ())), untest.get(a, 0)))
    print()
    print("    THE DECISIVE COMPARISON. MTP4 receives nothing. If its amino-acid")
    print("    completeness is HIGH, the model and the annotation agree that it can")
    print("    make its own, and the empty column is a result. If it is LOW -- as low")
    print("    as MAG_011's, say -- then two methods disagree, the model was filled")
    print("    to self-sufficiency, and the empty column is the gap-filler talking.")

# ================================================================ D
print()
print("================================================================")
print("D. WHAT WAS UNTESTABLE")
print("================================================================")
print("      %-16s %12s %12s" % ("mag", "untestable", "unmapped"))
for g in sorted(set(list(untest) + list(unmapped) + list(gets))):
    print("      %-16s %12d %12d" % (g, untest.get(g, 0), unmapped.get(g, 0)))
print()
print("    'untestable' = aspartate, glutamate, alanine, asparagine, phenylalanine:")
print("    real transfers with no 19b pathway to check them against. 'unmapped' =")
print("    everything else -- nucleosides, organic acids, PPi. Neither is evidence")
print("    of anything; both are counted so the test's coverage is visible.")

with open(os.path.join(out, "paired_test.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["mag", "genome_completeness", "n_received_pathways",
                "mean_completeness_received", "mean_completeness_other", "difference"])
    for g, d in diffs:
        got = [comp[g][p] for p in gets[g] if p in comp[g]]
        rest = [comp[g][p] for p in comp[g] if p not in gets[g]]
        w.writerow([g, "%.2f" % gcomp.get(g, float("nan")), len(got),
                    "%.4f" % (sum(got) / len(got)), "%.4f" % (sum(rest) / len(rest)),
                    "%+.4f" % d])
print()
print("    wrote %s" % os.path.join(out, "paired_test.tsv"))
PYEOF

say "E. WHAT THIS DECIDES, IN ONE LINE EACH"
cat <<'EOT'
  MTP4's amino-acid completeness HIGH   -> the empty column is a result. The
                                           chapter can say the reconstruction
                                           predicts MTP4 depends on nobody for
                                           amino acids -- consistent with fig04,
                                           where its route is the non-competitive
                                           substrate rather than a partner
  MTP4's amino-acid completeness LOW    -> two methods disagree; the empty column
                                           is gap-filling and must be reported as
                                           a model limitation, not a finding

  The paired test lopsidedly negative   -> the two methods agree about WHICH
                                           compounds move, not just which genomes.
                                           That is the strongest statement this
                                           reconstruction can support
  The paired test near even             -> agreement is at genome level only. The
                                           figure stands as capacity; the compound
                                           identities stay in the model's voice
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
