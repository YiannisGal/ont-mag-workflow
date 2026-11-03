#!/usr/bin/env bash
# 62_mtmB_flagged.sh -- what are the two references taking 61% of mtmB?
#
#   bash scripts/run.sh scripts/62_mtmB_flagged_XXXXXXX.sh
#
# WHAT 61 SETTLED, AND WHAT IT LEFT
# The decoys worked. 34_unbinned's 1,978 sediment mtmB reads became 85 once
# mtmC, mtbB, mttB, mtaB and the MttB superfamily were all able to win a read:
# 568 reads went to siblings instead. So 96% of 34's mtmB count was the family,
# not the gene, and 34's row must not be quoted.
#
# But 85 is still 3.5x the mcrA count of 24, and the ceiling argument is if
# anything UNDERSTATED by gene length: 34's own length model gives mcrA a
# read+gene span of 2,914 bp against mtmB's 1,864, so per gene copy mcrA should
# collect 1.56x MORE reads. Corrected for that, mtmB is running ~5.5x above a
# marker that every methanogen carries and it does not.
#
# AND 61 PRINTED THE LIKELY REASON WITHOUT ACTING ON IT
#   mtmB|REF|GFP32384.1    31.8% of 85 reads   8.9x its expected share
#   mtmB|REF|CAB1368335.1  29.4% of 85 reads   8.2x its expected share
# Two of twenty-eight references take 61% of the count. That is the same shape
# as the mtaB reference that took 98% and the two dsrB references that took 55%
# and 40% -- and 34 REMOVED those reads. 61 only flagged them.
#
# WHAT THIS DOES, all locally, in seconds:
#   A  names those references from refs_raw.faa, which kept the NCBI description
#      and the source organism in brackets
#   B  asks what the reads they captured hit SECOND. If a read's runner-up is a
#      decoy, that read is superfamily sequence that happened to land on a
#      mis-annotated entry -- which is evidence, not an assumption
#   C  recomputes every count with flagged references excluded, 34's convention
#   D  re-tests the ceiling on the corrected numbers
#
# No searching, no network. Read-only except 62_mtmB_flagged/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/62_mtmB_flagged"; mkdir -p "$OUT"
P61="$WS2/61_mtmB"
exec > >(tee "$OUT/62.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
ok(){ printf '  [OK]   %s\n' "$*"; }; bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "62 -- the flagged mtmB references"
echo "date : $(date -Is)"
echo "out  : $OUT"

for f in "$P61/SED_hits.tsv" "$P61/panel.faa"; do
  [ -s "$f" ] || { bad "missing $f -- run 61 first"; echo "STEP_STATUS: FAILED"; exit 1; }
done
ok "using 61's hits and panel"

NSED=$(awk '$1=="N"{print $2; exit}' "$WS2/34_unbinned/read_len_hist.tsv" 2>/dev/null)
NSED="${NSED:-6827767}"

python3 - "$P61" "$OUT" "$NSED" <<'PY'
import sys, os
from collections import defaultdict
P61, OUT, NSED = sys.argv[1], sys.argv[2], int(sys.argv[3])
MINLEN = 30

# ---- panel composition
nref = defaultdict(int)
for line in open(os.path.join(P61, 'panel.faa')):
    if line.startswith('>'):
        nref[line[1:].split('|')[0]] += 1

# ---- NCBI descriptions, kept verbatim by 61 in refs_raw.faa
desc = {}
raw = os.path.join(P61, 'refs_raw.faa')
if os.path.exists(raw):
    for line in open(raw, errors='replace'):
        if line.startswith('>'):
            acc = line[1:].split()[0]
            desc[acc] = line[1:].strip()

def load(tag):
    """best hit per read, plus the full ranked hit list for the runner-up test"""
    best, allh = {}, defaultdict(list)
    p = os.path.join(P61, '%s_hits.tsv' % tag)
    if not os.path.exists(p): return best, allh
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: ln, bs = int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        allh[f[0]].append((bs, f[1]))
        cur = best.get(f[0])
        if cur is None or bs > cur[1]:
            best[f[0]] = (f[1], bs)
    return best, allh

for tag in ('SED', 'LIQ'):
    best, allh = load(tag)
    if not best:
        continue
    cnt = defaultdict(int); perref = defaultdict(lambda: defaultdict(int))
    for q, (ref, _bs) in best.items():
        mk = ref.split('|')[0]
        cnt[mk] += 1; perref[mk][ref] += 1

    print('\n' + '=' * 72)
    print('  %s' % tag)
    print('=' * 72)

    # ---------------------------------------------------------- A  who are they
    flagged = set()
    print('\n  A. REFERENCES TAKING MORE THAN 4x THEIR EXPECTED SHARE')
    for mk in sorted(perref):
        if mk == 'CONTROL_arc':
            continue                      # a control is not a count; see note below
        n = max(1, nref.get(mk, 1)); exp = 1.0 / n; tot = cnt[mk]
        if tot < 20: continue
        for ref, c in sorted(perref[mk].items(), key=lambda kv: -kv[1]):
            share = c / tot
            if share >= 4 * exp and share > 0.20:
                flagged.add(ref)
                parts = ref.split('|')
                acc = parts[-1]
                d = desc.get(acc, '')
                org = ''
                if d.endswith(']') and '[' in d:
                    org = d[d.rfind('[') + 1:-1]
                    d = d[:d.rfind('[')].strip()
                print('    %-12s %5.1f%% of %-5d reads  %5.1fx expected   %s'
                      % (mk, 100 * share, tot, share / exp, ref))
                if d:
                    print('        description : %s' % d.split(' ', 1)[-1][:88])
                    print('        organism    : %s' % (org or '(none given)'))
                elif 'OWN' in parts:
                    print('        one of OUR OWN proteins, above its KOfam threshold')
                else:
                    print('        (no description recorded)')
    if not flagged:
        print('    none')

    # ---------------------------------------------------------- B  runner-up test
    print('\n  B. WHAT DID THOSE READS HIT SECOND?')
    print('     If a read captured by a flagged mtmB reference has a DECOY as its')
    print('     runner-up, the read is superfamily sequence and the flagged entry')
    print('     is very likely mis-annotated. This is measured, not assumed.')
    for mk in ('mtmB',):
        fl = [r for r in flagged if r.startswith(mk + '|')]
        if not fl:
            print('     no flagged %s reference' % mk); continue
        runner = defaultdict(int); nsecond = 0
        for q, (ref, _bs) in best.items():
            if ref not in fl: continue
            hits = sorted(allh[q], reverse=True)
            second = next((h[1] for h in hits if h[1] != ref), None)
            if second is None:
                runner['(no second hit)'] += 1
            else:
                runner[second.split('|')[0]] += 1; nsecond += 1
        tot = sum(runner.values())
        for k, v in sorted(runner.items(), key=lambda kv: -kv[1]):
            print('     %-16s %5d  (%.0f%%)' % (k, v, 100 * v / max(1, tot)))
        dec = sum(v for k, v in runner.items() if k.startswith('DECOY'))
        if nsecond:
            print('     -> %.0f%% of the reads with a runner-up have a DECOY second'
                  % (100 * dec / nsecond))

    # ---------------------------------------------------------- C  corrected
    print('\n  C. COUNTS WITH FLAGGED REFERENCES EXCLUDED (34\'s convention)')
    cnt2 = defaultdict(int)
    for q, (ref, _bs) in best.items():
        if ref in flagged: continue
        cnt2[ref.split('|')[0]] += 1
    print('    %-14s %10s %10s %12s' % ('marker', 'raw', 'corrected', 'per 1M'))
    for mk in sorted(set(cnt) | set(cnt2), key=lambda k: -cnt.get(k, 0)):
        pm = ('%12.2f' % (1e6 * cnt2.get(mk, 0) / NSED)) if tag == 'SED' else '%12s' % '-'
        print('    %-14s %10d %10d %s' % (mk, cnt.get(mk, 0), cnt2.get(mk, 0), pm))

    # ---------------------------------------------------------- D  ceiling again
    if tag == 'SED':
        print('\n  D. THE CEILING, ON THE CORRECTED NUMBERS')
        m, c = cnt2.get('mtmB', 0), cnt2.get('mcrA', 0)
        print('     mcrA %d   mtmB %d' % (c, m))
        if c == 0:
            print('     >> no mcrA left to test against')
        else:
            r = m / c
            print('     raw ratio            : %.2f' % r)
            print('     length-corrected     : %.2f   (mcrA read+gene 2,914 bp vs'
                  ' mtmB 1,864 bp, so mcrA should collect 1.56x more per copy)' % (r * 1.56))
            if r <= 1.0:
                print('     >> WITHIN THE CEILING. mtmB is reportable at this count,'
                      ' with the exclusions stated.')
            else:
                print('     >> STILL ABOVE. Report mtmB as present but UNQUANTIFIED.')

print('''
  A NOTE ON THE CONTROL, WHICH THE SAME RULE FLAGS
  CONTROL_arc|ETO_LIQ_MAG_003|ETOLIQMAG003_01125 took 48% of the control reads,
  6.3x its expected share. That is NOT promiscuity and is not excluded here: the
  control is 13 ribosomal/EF-1a sequences from 8 genomes, and the sediment
  archaea are simply closer to one of them than to the others. The 4x/1-in-n
  rule detects an UNEVEN DISTRIBUTION, which is evidence of promiscuity only
  when the references are supposed to be interchangeable. For a marker panel
  they are; for a control set drawn from different genomes they are not.
  Applying the rule blindly to the control would have thrown away the proof that
  the search works.''')
PY

say "WHAT THIS DECIDES"
cat <<'EOT'
  If section B shows most of the captured reads with a DECOY runner-up, the two
  flagged entries are mis-annotated superfamily proteins, excluding them is
  correct, and section C's corrected mtmB is the number to consider.

  If the runner-ups are mostly mtmB, the reads really are mtmB-like and the
  excess over mcrA is not a reference artefact. Then the honest conclusion is
  the harder one: a 50%-identity blastx on ONT reads cannot resolve this family
  at all, and no read-level mtmB number belongs in the chapter.

  EITHER WAY, two corrections stand:
    1. 34_unbinned's mtmB row -- 1,978 sediment reads, sed_over_rpoB 0.3253 --
       is 96% superfamily and must not be quoted.
    2. 34's promiscuity rule of "one reference over 25%" is blind in small
       reference sets: with 7 references the expected share is already 14%.
       Restate it as a multiple of 1/n.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
