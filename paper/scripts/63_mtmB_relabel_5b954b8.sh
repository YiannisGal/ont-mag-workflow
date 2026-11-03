#!/usr/bin/env bash
# 63_mtmB_relabel.sh -- the mtmB target set is mostly bacteria. Relabel, recount.
#
#   bash scripts/run.sh scripts/63_mtmB_relabel_XXXXXXX.sh
#
# WHAT THE ACCESSION LISTING SHOWED
#   AAF68952.2    Methanosarcina barkeri          genuine methanogen mtmB
#   AAM03597.1    Methanosarcina acetivorans C2A  genuine
#   AAM06345.1    Methanosarcina acetivorans C2A  genuine
#   CAD2231239.1  Desulfosporosinus meridiei      BACTERIUM
#   CAD2231240.1  Desulfosporosinus meridiei      BACTERIUM
#   GFP26094.1    Ca. Hakubella thermalkaliphila  BACTERIUM (partial)
#   GFP26864.1    Ca. Hakubella thermalkaliphila  BACTERIUM
#   GFP32384.1    Ca. Hakubella thermalkaliphila  BACTERIUM (partial)
#   GFP43010.1    Ca. Hakubella thermalkaliphila  BACTERIUM
#   CAB1368335.1  Denitratisoma oestradiolicum    BACTERIUM
#
# MY FETCH PUT THEM THERE. The term
#   "methylamine--corrinoid protein Co-methyltransferase"[Protein Name]
# is what NCBI calls the MttB SUPERFAMILY, not specifically methanogen mtmB, so
# it imported seven bacterial proteins straight into the target set. 62 flagged
# only two of them, because the other five each sat under the 4x/20% threshold.
# Every read they hold is still counted as mtmB.
#
# NO RE-SEARCH IS NEEDED. Which reference a read hit best does not change when
# the reference is renamed, so this is a pure relabel of 61's existing hits:
# the seven bacterial entries move to a DECOYbact bucket and the mtmB count
# falls out correctly. Exact, and it takes seconds.
#
# THE DESULFOSPOROSINUS ENTRIES ARE NOT A COINCIDENCE. That is the same genus as
# LIQ_contig_12880, the unbinned sulfate reducer identified in steps 55-59 at
# 88-90% DsrAB identity to D. meridiei. A Desulfosporosinus in this sediment
# carrying MttB-superfamily methylamine methyltransferases is a plausible source
# of exactly these reads -- which is a reason to move them out of mtmB, and an
# observation worth a sentence of its own.
#
# Read-only except 63_mtmB_relabel/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/63_mtmB_relabel"; mkdir -p "$OUT"
P61="$WS2/61_mtmB"
exec > >(tee "$OUT/63.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
ok(){ printf '  [OK]   %s\n' "$*"; }; bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "63 -- mtmB relabelled by organism"
echo "date : $(date -Is)"

[ -s "$P61/SED_hits.tsv" ] || { bad "run 61 first"; echo "STEP_STATUS: FAILED"; exit 1; }

# the seven bacterial entries, by accession
cat > "$OUT/bacterial.txt" <<'EOT'
CAD2231239.1
CAD2231240.1
GFP26094.1
GFP26864.1
GFP32384.1
GFP43010.1
CAB1368335.1
EOT
ok "$(wc -l < "$OUT/bacterial.txt") bacterial entries to move out of mtmB"

NSED=$(awk '$1=="N"{print $2; exit}' "$WS2/34_unbinned/read_len_hist.tsv" 2>/dev/null)
NSED="${NSED:-6827767}"

say "COUNTS BEFORE AND AFTER THE RELABEL"
python3 - "$P61" "$OUT" "$NSED" <<'PY'
import sys, os
from collections import defaultdict
P61, OUT, NSED = sys.argv[1], sys.argv[2], int(sys.argv[3])
MINLEN = 30
bact = {l.strip() for l in open(os.path.join(OUT, 'bacterial.txt')) if l.strip()}

nref = defaultdict(int)
for line in open(os.path.join(P61, 'panel.faa')):
    if line.startswith('>'):
        nref[line[1:].split('|')[0]] += 1

def relabel(ref):
    p = ref.split('|')
    return 'DECOYbact' if p[-1] in bact else p[0]

for tag in ('SED', 'LIQ'):
    p = os.path.join(P61, '%s_hits.tsv' % tag)
    if not os.path.exists(p): continue
    best = {}
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: ln, bs = int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        cur = best.get(f[0])
        if cur is None or bs > cur[1]:
            best[f[0]] = (f[1], bs)

    before = defaultdict(int); after = defaultdict(int)
    whichbact = defaultdict(int)
    for _q, (ref, _bs) in best.items():
        before[ref.split('|')[0]] += 1
        lab = relabel(ref)
        after[lab] += 1
        if lab == 'DECOYbact':
            whichbact[ref.split('|')[-1]] += 1

    print('\n  %s' % tag)
    print('    %-14s %9s %11s %12s' % ('marker', 'before', 'after', 'per 1M'))
    for mk in sorted(set(before) | set(after), key=lambda k: -before.get(k, 0)):
        pm = ('%12.2f' % (1e6 * after.get(mk, 0) / NSED)) if tag == 'SED' else '%12s' % '-'
        print('    %-14s %9d %11d %s' % (mk, before.get(mk, 0), after.get(mk, 0), pm))

    if whichbact:
        print('\n    reads held by each bacterial entry:')
        for a, c in sorted(whichbact.items(), key=lambda kv: -kv[1]):
            print('      %-16s %5d' % (a, c))

    m, c = after.get('mtmB', 0), after.get('mcrA', 0)
    if c:
        print('\n    mtmB / mcrA = %.2f' % (m / c))
PY

say "HOW TO READ IT"
cat <<'EOT'
  The mtmB row "after" is the count that belongs in the chapter. Everything the
  seven bacterial entries were holding is superfamily sequence from organisms
  that are not methanogens, and it was being reported as a methanogen marker.

  The mtmB/mcrA ratio is denominator-free, so it needs no assumption about
  library size. Compare the two samples: the cultures were enriched with no
  sulfate under H2/CO2, which selects hydrogenotrophs, so a LOWER mtmB/mcrA in
  the cultures than in the sediment is the expected direction and is the result
  worth reporting.

  Calibration for the sediment value: our own catalogue holds 11 mtmB proteins
  and 9 mcrA above their curated thresholds. Corrected for read+gene span
  (mtmB 1,864 bp, mcrA 2,914 bp) that predicts a read ratio near 0.78.

  AND ONE THING TO CARRY FORWARD. Two of the seven were Desulfosporosinus
  meridiei -- the same genus as LIQ_contig_12880, the unbinned sulfate reducer
  from steps 55-59. Bacterial MttB-superfamily methyltransferases in this
  sediment are not hypothetical; one of their likely carriers is already
  identified in this dataset.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
