#!/usr/bin/env bash
# 65_refset_rarefy.sh -- is mtbB's excess a gene or a reference count?
#
#   bash scripts/run.sh scripts/65_refset_rarefy_XXXXXXX.sh
#
# THE CROSS-CHECK IN 64 FAILED, AND IT WAS RIGHT TO. mtmB came back at 13 where
# step 63 measured 18. The panels differ -- 63 gave mtmB 28 references and no
# bacterial decoy bucket; 64 gave it 19 archaeal references and 79 DECOYbact --
# so five reads moved. At counts of 13 to 18, a read-level number is sensitive to
# panel composition at the ~30% level. That is a fact about the measurement and
# it belongs in the Methods, not a bug to be tuned away.
#
# WHICH EXPOSES A CONFOUND IN 64's CEILING TEST. The panel is not balanced:
#
#     marker   references   reads   reads per reference
#     mtbB         61        211           3.46
#     mcrA          9         24           2.67
#     mtmB         19         13           0.68
#     mtaB         19         13           0.68
#
# mtbB is 8.79x mcrA in raw reads and 1.3x mcrA PER REFERENCE. A marker with 61
# references samples far more sequence space than one with 9, so it recruits more
# divergent reads -- that is sensitivity, not abundance, and the expected-ratio
# formula in 64 accounts for gene copies and gene length but NOT for how much of
# the family each panel can see. Declaring mtbB "still measuring a family" on an
# unbalanced panel would repeat in a new form the mistake this whole sequence has
# been about.
#
# THE TEST. Rarefy every marker's reference set to the SAME SIZE, repeatedly, and
# recompute. No re-searching is needed: which reference a read hit best does not
# change, so each replicate simply drops the unsampled references and reassigns
# every read to its best surviving hit (diamond kept up to five per read). Reads
# with nothing left are counted and reported rather than silently lost.
#
# IF mtbB's excess survives rarefaction, it is the gene family and the marker
# stays unreportable. IF it collapses toward its expectation, the excess was the
# panel and mtbB is reportable at its rarefied ratio, with the panel stated.
#
# Seconds, local, no network. Read-only except 65_rarefy/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/65_rarefy"; mkdir -p "$OUT"
P64="$WS2/64_mtbB_mtaB"
exec > >(tee "$OUT/65.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
ok(){ printf '  [OK]   %s\n' "$*"; }; bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "65 -- reference-set rarefaction"
echo "date : $(date -Is)"

for f in "$P64/SED_hits.tsv" "$P64/panel.faa"; do
  [ -s "$f" ] || { bad "run 64 first ($f missing)"; echo "STEP_STATUS: FAILED"; exit 1; }
done
ok "using 64's hits and panel"

NSED=$(awk '$1=="N"{print $2; exit}' "$WS2/34_unbinned/read_len_hist.tsv" 2>/dev/null)
NSED="${NSED:-6827767}"
REPS="${REPS:-200}"

python3 - "$P64" "$OUT" "$NSED" "$REPS" <<'PY'
import sys, os, random, statistics as st
from collections import defaultdict
P64, OUT, NSED, REPS = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
MINLEN = 30
SPAN = {'mcrA': 2914, 'mtmB': 1864, 'mtbB': 2278, 'mtaB': 2593}
COPIES = {}
cp = os.path.join(P64, 'own_copies.tsv')
if os.path.exists(cp):
    for line in open(cp):
        k, v = line.split('\t'); COPIES[k] = int(v)

refs = defaultdict(list)
for line in open(os.path.join(P64, 'panel.faa')):
    if line.startswith('>'):
        name = line[1:].strip()
        refs[name.split('|')[0]].append(name)

def hits(tag):
    """read -> list of (bitscore, reference), best first"""
    h = defaultdict(list)
    p = os.path.join(P64, '%s_hits.tsv' % tag)
    if not os.path.exists(p): return h
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: ln, bs = int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        h[f[0]].append((bs, f[1]))
    for k in h: h[k].sort(reverse=True)
    return h

MARK = ('mcrA', 'mtmB', 'mtbB', 'mtaB')

for tag in ('SED', 'LIQ'):
    H = hits(tag)
    if not H: continue
    print('\n' + '=' * 70)
    print('  %s' % tag)
    print('=' * 70)

    full = defaultdict(int)
    for q, hl in H.items(): full[hl[0][1].split('|')[0]] += 1
    print('\n  as measured in 64 (unbalanced panel)')
    print('    %-12s %8s %8s %12s' % ('marker', 'reads', 'n refs', 'reads/ref'))
    for m in sorted(full, key=lambda k: -full[k]):
        n = len(refs.get(m, []))
        print('    %-12s %8d %8d %12.2f' % (m, full[m], n, full[m] / max(1, n)))

    # rarefy to the smallest TARGET/CEILING reference set
    N = min(len(refs[m]) for m in MARK if refs.get(m))
    print('\n  rarefied: every reference set cut to %d, %d replicates' % (N, REPS))
    res = defaultdict(list); lost = []
    rng = random.Random(20260913)
    for _ in range(REPS):
        keep = set()
        for m, rl in refs.items():
            keep |= set(rl if len(rl) <= N else rng.sample(rl, N))
        c = defaultdict(int); nl = 0
        for q, hl in H.items():
            s = next((r for _bs, r in hl if r in keep), None)
            if s is None: nl += 1
            else: c[s.split('|')[0]] += 1
        for m in MARK: res[m].append(c.get(m, 0))
        lost.append(nl)
    print('    reads with no surviving hit: %.0f of %d (%.1f%%)'
          % (st.mean(lost), len(H), 100 * st.mean(lost) / max(1, len(H))))

    print('\n    %-12s %16s %16s' % ('marker', 'reads (mean sd)', 'per 1M'))
    for m in MARK:
        v = res[m]
        pm = ('%16s' % ('%.2f' % (1e6 * st.mean(v) / NSED))) if tag == 'SED' else '%16s' % '-'
        print('    %-12s %10.1f %4.1f %s'
              % (m, st.mean(v), st.pstdev(v), pm))

    if tag == 'SED':
        mc = res['mcrA']
        print('\n    RATIO TO mcrA, rarefied, against each marker own expectation')
        print('    %-8s %14s %14s %12s  %s' %
              ('marker', 'expected r', 'observed r', 'obs/exp', 'verdict'))
        for m in ('mtmB', 'mtbB', 'mtaB'):
            if not COPIES.get(m) or not COPIES.get('mcrA'): continue
            exp = (COPIES[m] * SPAN[m]) / (COPIES['mcrA'] * SPAN['mcrA'])
            rr = [a / b for a, b in zip(res[m], mc) if b]
            if not rr: continue
            obs = st.mean(rr)
            f = obs / exp
            verdict = ('consistent with the gene' if f <= 2.5 else
                       'STILL a family-wide signal')
            print('    %-8s %14.2f %14.2f %12.2f  %s' % (m, exp, obs, f, verdict))

        print('''
    Compare with 64 on the unbalanced panel: mtbB observed 8.79, obs/exp 9.2.
    If the rarefied obs/exp is near 1, that excess was reference count -- the
    marker with 61 references was simply able to see more of the family than
    mcrA with 9. If it is still far above, it is the gene family and mtbB
    stays out of the chapter.''')
PY

say "AND THE mtmB STABILITY QUESTION"
cat <<'EOT'
  63 measured 18, 64 measured 13, on panels that differ in reference count and
  in whether a bacterial decoy bucket exists. Both are honest numbers from
  honest panels. What that means for the chapter:

    DO NOT write a single read count for mtmB. Write the range and its cause:
    "13-18 sediment reads depending on reference-panel composition (1.9-2.6 per
    million)".

    DO write the direction, which is stable across both panels: mtmB per mcrA is
    higher in the sediment than in the enrichment cultures -- 0.75 vs 0.32 on
    63's panel, 0.54 vs 0.31 on 64's. Between 1.7x and 2.3x, same sign, same
    conclusion. The enrichment selected hydrogenotrophs; the sediment carries
    relatively more methylamine methyltransferase.

  A ratio of two markers measured in the SAME run cancels most of the panel
  sensitivity. A raw count does not. That is the general lesson and it applies
  to every read-level number in this chapter.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
