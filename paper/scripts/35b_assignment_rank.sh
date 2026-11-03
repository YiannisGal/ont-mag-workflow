#!/usr/bin/env bash
# =============================================================================
# 35b_assignment_rank.sh  --  where do Kaiju's reads ACTUALLY land, by rank?
#
# THE PROBLEM WITH STEP 35's SECTION 4
#   It reported this, and I told you to cite it as the novelty measure:
#
#     rank        SED taxa   SED % reads
#     phylum           139        37.27%
#     class            152        26.54%
#     order            305        19.38%
#     family           725        14.90%
#     genus           3763        12.49%
#     species        21508        38.48%     <-- HIGHER THAN GENUS
#
#   A rank-assignment curve MUST be non-increasing as rank deepens: a read that
#   can be placed in a genus can always be placed in a family. 38.48% at species
#   against 12.49% at genus is impossible as a decay curve, so the quantity is
#   not what I said it was.
#
#   The cause is not Kaiju. `kaiju2table -r <rank>` walks each read's assigned
#   taxon UP the NCBI lineage to the requested rank and reports it only if a
#   node of exactly that rank exists. Large parts of the progenomes taxonomy —
#   candidate divisions, "unclassified <phylum>", MAG-derived lineages — have a
#   species-ranked leaf and NO genus-ranked ancestor. Those reads count at
#   species and vanish at genus. The intermediate ranks are therefore measuring
#   LINEAGE COMPLETENESS IN NCBI, not novelty in Etoliko, and the difference
#   phylum% - genus% = 24.8 points that step 35 printed is not interpretable.
#
# WHAT THIS STEP MEASURES INSTEAD
#   The rank of the taxon each read was ACTUALLY ASSIGNED TO, read straight out
#   of the .kaiju output and nodes.dmp. That needs no lineage walk, so it cannot
#   be distorted by missing intermediate nodes. It answers the real question:
#   how deep does the evidence go?
#
#   It then quantifies the artefact directly: of the reads assigned at species
#   level, how many have no genus-ranked ancestor at all. That number IS the
#   size of the distortion in step 35's section 4.
#
# CONTROLS
#   1. The classified-read total recomputed here MUST equal step 35's section 1
#      (2,856,701 sediment). Two independent passes over the same file.
#   2. The per-rank assignment counts MUST sum to the classified total. A
#      partition that does not sum is a parsing failure, not a result.
#   3. Every taxid in the .kaiju output must resolve in nodes.dmp. Unresolvable
#      taxids are counted and reported, never silently dropped.
#
# OUTPUT  35_kaiju/assignment_rank.tsv · ASSIGNMENT_RANK.txt
# RUNTIME 3-8 min (one pass over each .kaiju file)
# USAGE   ws2 && bash scripts/35b_assignment_rank.sh 2>&1 | tee logs/35b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="35b_assignrank"
D="$WS2/35_kaiju"

echo "######## 35b — WHERE DO THE READS ACTUALLY LAND? ########"
echo "  started: $(date -Is)"

[ -s "$D/SED.kaiju" ] || { echo "  !! no $D/SED.kaiju — run step 35 first."; exit 2; }

DBDIR="${KAIJU_DBDIR:-${DB_DIR}/kaiju}"
NODES=""; NAMES=""
for n in "$DBDIR/nodes.dmp" "$D/nodes.dmp"; do [ -s "$n" ] && { NODES="$n"; break; }; done
for n in "$DBDIR/names.dmp" "$D/names.dmp"; do [ -s "$n" ] && { NAMES="$n"; break; }; done
if [ -z "$NODES" ]; then
  echo "  !! nodes.dmp not found under $DBDIR — the rank of a taxon cannot be"
  echo "     looked up without it, and guessing is not an option."
  exit 3
fi
echo "  nodes : $NODES"
echo "  names : ${NAMES:-<not found — taxon names will be omitted>}"

echo
python3 - "$D" "$NODES" "$NAMES" <<'PY'
import sys, os
from collections import Counter, defaultdict

D, NODES, NAMES = sys.argv[1], sys.argv[2], sys.argv[3]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

# ---------------------------------------------------------- taxonomy
say("  reading nodes.dmp ...")
parent, rank = {}, {}
with open(NODES, errors='replace') as fh:
    for line in fh:
        f = [x.strip() for x in line.split('|')]
        if len(f) < 3:
            continue
        try:
            t = int(f[0]); p = int(f[1])
        except ValueError:
            continue
        parent[t] = p
        rank[t] = f[2]
say(f"  {len(rank):,} taxa")

name = {}
if NAMES and os.path.exists(NAMES):
    with open(NAMES, errors='replace') as fh:
        for line in fh:
            f = [x.strip() for x in line.split('|')]
            if len(f) >= 4 and f[3] == 'scientific name':
                try:
                    name[int(f[0])] = f[1]
                except ValueError:
                    pass

MAJOR = ('superkingdom', 'domain', 'phylum', 'class', 'order',
         'family', 'genus', 'species')

def lineage_ranks(t):
    """set of ranks present on the path from t to the root"""
    out, seen, hops = set(), set(), 0
    while t in parent and t not in seen and hops < 100:
        seen.add(t); hops += 1
        r = rank.get(t)
        if r:
            out.add(r)
        p = parent[t]
        if p == t:
            break
        t = p
    return out

# ---------------------------------------------------------- the reads
for code in ('SED', 'LIQ'):
    p = os.path.join(D, f'{code}.kaiju')
    if not os.path.exists(p):
        continue
    rule(f'{code} — RANK OF THE TAXON EACH READ WAS ASSIGNED TO')
    tot = cls = unresolved = 0
    byrank = Counter()
    toptax = Counter()
    # of species-level assignments, how many have no genus ancestor
    sp_no_genus = sp_total = 0
    with open(p, errors='replace') as fh:
        for line in fh:
            tot += 1
            if line[:1] != 'C':
                continue
            cls += 1
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                unresolved += 1; continue
            try:
                t = int(f[2])
            except ValueError:
                unresolved += 1; continue
            r = rank.get(t)
            if r is None:
                unresolved += 1; continue
            byrank[r] += 1
            toptax[t] += 1
            if r == 'species':
                sp_total += 1
                if 'genus' not in lineage_ranks(t):
                    sp_no_genus += 1

    # ---- CONTROL 2: the partition must sum to the classified total
    partsum = sum(byrank.values()) + unresolved
    say(f"  reads total      : {tot:,}")
    say(f"  classified       : {cls:,}   ({100.0 * cls / max(1, tot):.2f}%)")
    say(f"  rank-resolved    : {sum(byrank.values()):,}")
    say(f"  taxid unresolved : {unresolved:,}")
    if partsum != cls:
        say(f"  !! PARTITION DOES NOT SUM: {partsum:,} vs {cls:,} classified.")
        say("     That is a parsing failure, not a result. Nothing below is usable.")
        continue
    say("  [PASS] the per-rank counts partition the classified reads exactly.")
    say()
    say(f"  {'assigned at rank':<22}{'reads':>14}{'% of classified':>18}{'% of all reads':>17}")
    say('-' * W)
    for r, n in byrank.most_common(14):
        say(f"  {r:<22}{n:>14,}{100.0 * n / max(1, cls):>17.2f}%"
            f"{100.0 * n / max(1, tot):>16.2f}%")
    say()

    # ---- the artefact, quantified
    rule(f'{code} — WHY STEP 35 SHOWED SPECIES ABOVE GENUS')
    say("  `kaiju2table -r genus` walks each read's taxon UP to a genus-ranked")
    say("  ancestor and drops the read if there is not one. A species-ranked leaf")
    say("  with no genus above it therefore counts at species and disappears at")
    say("  genus, which inverts the curve.")
    say()
    if sp_total:
        say(f"  reads assigned at species rank        : {sp_total:,}")
        say(f"  ...of those, with NO genus ancestor   : {sp_no_genus:,}"
            f"   ({100.0 * sp_no_genus / sp_total:.1f}%)")
        say()
        if sp_no_genus > 0.2 * sp_total:
            say("  >> That is the artefact, measured. Step 35's section 4 intermediate")
            say("     ranks are reporting NCBI lineage completeness as much as they are")
            say("     reporting this community's novelty, and the phylum-minus-genus")
            say("     difference it printed is NOT a novelty statistic. Use the")
            say("     assignment-rank table above instead.")
        else:
            say("  >> Too few to explain the inversion. Something else is going on —")
            say("     do not use either table until it is understood.")
    else:
        say("  no species-rank assignments at all, so this is not the explanation.")
    say()
    say(f"  {'most-assigned taxa':<46}{'reads':>14}{'rank':>14}")
    say('-' * W)
    for t, n in toptax.most_common(10):
        say(f"  {(name.get(t, str(t)))[:45]:<46}{n:>14,}{rank.get(t, '?'):>14}")

rule('WHAT TO REPORT')
say("  1. The citable depth measure is the ASSIGNMENT-RANK table: the rank of")
say("     the taxon each read was actually assigned to. It needs no lineage")
say("     walk, so missing intermediate nodes cannot distort it.")
say("  2. Step 35's section 4 percentages at class/order/family/genus are")
say("     CONFOUNDED by NCBI lineage completeness. Do not quote the")
say("     phylum-minus-genus difference as a novelty statistic — I was wrong to")
say("     suggest it, and the species-above-genus inversion is the proof.")
say("  3. Step 35's sections 1, 2, 3 and the taxon COUNTS are unaffected: they")
say("     do not depend on walking a lineage to a named rank.")
say()
say("  LIMITS:")
say("   1. Assignment rank is a property of the DATABASE's taxonomy as much as")
say("      of the read. A database whose leaves are mostly species will assign")
say("      mostly at species. Name the database with the number.")
say("   2. Kaiju's LCA means a read matching two divergent references is pushed")
say("      up; a high-rank assignment is an admission, not a measurement of the")
say("      organism's actual novelty.")

with open(os.path.join(D, 'ASSIGNMENT_RANK.txt'), 'w') as fh:
    fh.write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'ASSIGNMENT_RANK.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "report written" "$([ -s "$D/ASSIGNMENT_RANK.txt" ] && echo 0 || echo 1)" || OK=1
if grep -q 'PARTITION DOES NOT SUM' "$D/ASSIGNMENT_RANK.txt" 2>/dev/null; then
  echo "  [STOP] the rank partition did not sum — nothing here is usable."; OK=1
fi
[ "$OK" -eq 0 ] && runlog "35b_assignrank complete" || runlog "**35b_assignrank FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
