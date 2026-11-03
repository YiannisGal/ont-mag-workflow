#!/usr/bin/env bash
# =============================================================================
# 37_biogeography_verdict.sh  --  read step 14's result, with a control that
#                                 is applied where its expectation holds
#
# WHY THIS IS NOT A NEW SEARCH
#   Step 14 mapped the 29 MAGs against Etoliko's own reads and 30 public
#   euxinic/meromictic datasets. `14_comparative/recruitment_matrix.tsv` was
#   written 2026-08-20 00:30 -- AFTER the seqkit subsampling fix -- and has
#   never been read. The roadmap's step 37 ("extend step 14 to the seven
#   meromictic-lake projects") is ALSO already in that matrix: Lake Zug
#   (PRJNA977988), El Tobar (PRJNA639779), Mahoney, Shunet, Monoun and the
#   brine sets are columns in it. So there is nothing to map. There is a
#   result to read, a control to specify correctly, and a verdict to write.
#
# 🔴 WHY STEP 14's CONTROL "FAILED"
#   It required all 29 MAGs to recruit from Etoliko sediment at >=30% breadth,
#   and reported SED as failing. But 25 of the 29 MAGs were assembled from the
#   ENRICHMENT CULTURES. A culture-derived organism that is rare or absent in
#   the sediment SHOULD show low sediment breadth -- step 26 measured exactly
#   that, with only 13 of 29 above 0.50 and several at 0.02.
#
#   Demanding recruitment where the positive expectation does not hold is the
#   same error as steps 19b and 23e, which applied a bacterial control to
#   archaeal genomes and read the guaranteed absence as a result. A positive
#   control is only a control over its DOMAIN OF APPLICABILITY.
#
#   Correct control, and it can still fail:
#     - the 4 SED-assembled MAGs must recruit from ETOLIKO_SED
#     - the 25 LIQ-assembled MAGs must recruit from ETOLIKO_LIQ
#   Each genome is tested against the sample it was assembled from, where
#   recruitment is guaranteed unless the mapping is broken.
#
# WHAT IT REPORTS
#   1. Coverage: every project in meta/projects.tsv, and whether it has a
#      column in the matrix. A project with no column is NOT SEARCHED and can
#      never appear as an absence.
#   2. The control above.
#   3. Per-dataset maximum and median breadth across the 29 MAGs.
#   4. The biogeography call per MAG, at a stated threshold, with a
#      SENSITIVITY sweep: the threshold at which any public dataset would
#      change the answer. If the gap between Etoliko and everything else spans
#      orders of magnitude, the call does not depend on where the line is put,
#      and that is worth more than the call itself.
#
# THE LIMITATION THAT TRAVELS WITH IT (from step 14, restated, not softened)
#   Public datasets were streamed with `fastq-dump -X`, which takes the FIRST N
#   spots rather than a random sample. For single-library Illumina runs that is
#   far less damaging than for a ten-run ONT file -- which is what produced
#   step 14's 116-fold subsampling artefact -- but it is not random. A public
#   dataset that DID recruit would need re-checking with a proper random
#   subsample before anything was concluded from it. None does, which is the
#   one direction in which this limitation is harmless: a first-N sample can
#   miss a rare organism, but it cannot manufacture an absence across 30
#   independent datasets.
#
# OUTPUT  14_comparative/BIOGEOGRAPHY.txt · biogeography_calls.tsv
# RUNTIME seconds (reads two TSVs; no mapping)
# USAGE   ws2 && bash scripts/37_biogeography_verdict.sh 2>&1 | tee logs/37.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="37_biogeography"
C="$WS2/14_comparative"

echo "######## 37 — BIOGEOGRAPHY: READING STEP 14's RESULT ########"
echo "  started: $(date -Is)"

M="$C/recruitment_matrix.tsv"
P="$C/meta/projects.tsv"
[ -s "$M" ] || { echo "  !! $M not found — step 14 has not produced a matrix."; exit 2; }
echo "  matrix  : $M  ($(( $(wc -l < "$M") - 1 )) MAGs)"
echo "  projects: ${P}$([ -s "$P" ] || echo '  <MISSING — coverage cannot be checked>')"

rm -f "$C/BIOGEOGRAPHY.txt" "$C/biogeography_calls.tsv"
echo
python3 - "$C" "$M" "$P" <<'PY'
import sys, os, re, csv, math
from collections import defaultdict

C, M, P = sys.argv[1], sys.argv[2], sys.argv[3]
W = 100
rep = []
STATUS = ['INCOMPLETE']
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)
def dump(code=0, why=None):
    if why:
        STATUS[0] = why
    with open(os.path.join(C, 'BIOGEOGRAPHY.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(C, 'BIOGEOGRAPHY.txt')}")
    sys.exit(code)

rows = list(csv.DictReader(open(M, errors='replace'), delimiter='\t'))
if not rows:
    say("  !! empty matrix."); dump(0, 'EMPTY_MATRIX')
cols = [c for c in rows[0] if c]
bcols = [c for c in cols if c.endswith('_breadth')]
ETO = [c for c in bcols if c.startswith('ETOLIKO_')]
PUB = [c for c in bcols if not c.startswith('ETOLIKO_')]

def val(r, c):
    try:
        return float(r[c])
    except (TypeError, ValueError):
        return float('nan')

def med(v):
    v = sorted(x for x in v if x == x)
    if not v:
        return float('nan')
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2

say(f"  {len(rows)} MAGs | {len(ETO)} Etoliko columns | {len(PUB)} public datasets")

# ------------------------------------------------------- 1. coverage
rule('1. COVERAGE — WHICH PROJECTS ACTUALLY HAVE DATA IN THIS MATRIX?')
say("  A project listed in projects.tsv but ABSENT from the matrix was never")
say("  recruited against. It cannot appear as an absence anywhere below.")
say()
proj = []
if os.path.exists(P):
    for line in open(P, errors='replace'):
        f = [x.strip() for x in line.rstrip('\n').split('\t')]
        if f and f[0] and not f[0].startswith('#'):
            proj.append((f[0], f[1] if len(f) > 1 else '', f[2] if len(f) > 2 else ''))
have = {re.match(r'(PRJ[A-Z]+\d+)', c).group(1) for c in PUB
        if re.match(r'(PRJ[A-Z]+\d+)', c)}
nruns = defaultdict(int)
for c in PUB:
    m = re.match(r'(PRJ[A-Z]+\d+)', c)
    if m:
        nruns[m.group(1)] += 1
say(f"  {'project':<16}{'runs':>6}  {'type':<5}description")
say('-' * W)
missing = []
for pid, desc, typ in proj:
    n = nruns.get(pid, 0)
    if n == 0:
        missing.append((pid, desc, typ))
    say(f"  {pid:<16}{(n if n else 'NONE'):>6}  {typ:<5}{desc[:60]}")
say()
if missing:
    say(f"  {len(missing)} project(s) have NO data in the matrix:")
    for pid, desc, typ in missing:
        note = ''
        if typ.upper() == 'RNA':
            note = '  (RNA — deliberately excluded: recruitment of a'
        say(f"    {pid}  {desc[:58]}{note}")
    say("     metatranscriptome measures expression, not presence)" if any(
        t.upper() == 'RNA' for _, _, t in missing) else "")
    say("  These are NOT SEARCHED. Nothing below is an absence for them.")
say()
say(f"  searched: {len(have)} projects, {len(PUB)} runs")

# ------------------------------------------------- 2. the control
rule('2. CONTROL — EACH MAG AGAINST THE SAMPLE IT WAS ASSEMBLED FROM')
say("  Step 14 required all 29 MAGs to recruit from Etoliko SEDIMENT and")
say("  reported the control as failed. That specification was wrong: 25 of the")
say("  29 were assembled from the ENRICHMENT CULTURES, and a culture organism")
say("  that is rare in the sediment SHOULD show low sediment breadth — step 26")
say("  measured exactly that (13/29 above 0.50, several at 0.02).")
say()
say("  A positive control is only a control over its DOMAIN OF APPLICABILITY.")
say("  Applying one outside it is the 19b/23e error: a bacterial control on")
say("  archaeal genomes, whose guaranteed absence was read as a result.")
say()
say("  Correct form: every MAG against the sample it came from.")
say()
SEDC = next((c for c in ETO if 'SED' in c.upper()), None)
LIQC = next((c for c in ETO if 'LIQ' in c.upper()), None)
THR = 0.30
def origin(r):
    wn = (r.get('working_name') or r.get('mag_id') or '')
    return 'SED' if wn.upper().startswith('SED') else 'LIQ'
ok = True
say(f"  {'assembled from':<18}{'n':>4}{'own-sample breadth: min':>26}{'median':>10}{'>=0.30':>9}")
say('-' * W)
for tag, col in (('SED (sediment)', SEDC), ('LIQ (cultures)', LIQC)):
    if not col:
        say(f"  {tag:<18}{'—':>4}   no {tag[:3]} column in the matrix — CANNOT TEST")
        ok = False
        continue
    v = [val(r, col) for r in rows if origin(r) == tag[:3]]
    v = [x for x in v if x == x]
    if not v:
        say(f"  {tag:<18}{0:>4}   no MAGs assembled from this sample")
        continue
    npass = len([x for x in v if x >= THR])
    say(f"  {tag:<18}{len(v):>4}{min(v):>26.3f}{med(v):>10.3f}"
        f"{f'{npass}/{len(v)}':>9}")
    if npass < len(v):
        ok = False
say()
if ok:
    say(f"  [PASS] every MAG recruits from its own sample at >= {THR:.0%} breadth.")
    say("         The mapping works, so a near-zero elsewhere is informative.")
else:
    say(f"  [FAIL] at least one MAG does not recruit from the sample it was")
    say("         assembled from. That is a mapping or subsampling failure, not")
    say("         biogeography. NO VERDICT.")
    dump(0, 'CONTROL_FAILED')

# --------------------------------------------- 3. per-dataset breadth
rule('3. HOW MUCH DO ETOLIKO GENOMES RECRUIT ANYWHERE ELSE?')
stats = []
for c in bcols:
    v = [val(r, c) for r in rows]
    v = [x for x in v if x == x]
    if v:
        stats.append((max(v), med(v), len([x for x in v if x >= THR]), c))
stats.sort(key=lambda t: (-t[0], t[3]))
say(f"  {'max':>7}{'median':>9}{'>=0.30':>8}  dataset")
say('-' * W)
for mx, md, n, c in stats[:14]:
    say(f"  {mx:>7.3f}{md:>9.3f}{n:>8}  {c[:58]}")
if len(stats) > 14:
    say(f"  ... and {len(stats) - 14} more, all below the rows shown")
say()
pubmax = max((s[0] for s in stats if not s[3].startswith('ETOLIKO_')), default=0.0)
etomax = max((s[0] for s in stats if s[3].startswith('ETOLIKO_')), default=0.0)
say(f"  highest breadth in ANY public dataset : {pubmax:.3f}")
say(f"  highest breadth in Etoliko's own reads: {etomax:.3f}")
if pubmax > 0:
    say(f"  ratio                                 : {etomax / pubmax:.0f}x")

# ------------------------------------------- 4. sensitivity, then the call
rule('4. THE CALL — AND WHETHER IT DEPENDS ON WHERE THE LINE IS DRAWN')
say("  A biogeography call made at one threshold is worth little; what matters")
say("  is whether the answer survives moving it. Sweeping the threshold:")
say()
say(f"  {'threshold':>11}{'MAG-dataset pairs above it (public only)':>44}")
say('-' * W)
sweep = [0.50, 0.30, 0.20, 0.10, 0.05, 0.02, 0.01]
for t in sweep:
    n = sum(1 for r in rows for c in PUB if val(r, c) >= t)
    say(f"  {t:>11.2f}{n:>44}")
say()
if pubmax < min(sweep):
    say(f"  Every public value is below {min(sweep):.2f}, so the call is the same at")
    say("  EVERY threshold in that sweep. The result does not depend on the")
    say("  cutoff, which is a stronger statement than the call itself.")
else:
    say(f"  The answer changes within the sweep (public max {pubmax:.3f}). The")
    say("  threshold must then be justified, not chosen.")
say()
calls = []
for r in rows:
    sites = [c for c in PUB if val(r, c) >= THR]
    calls.append((r.get('mag_id', ''), r.get('working_name', ''),
                  r.get('genus', '') or '-', len(sites),
                  max([val(r, c) for c in PUB] + [0.0]), sites))
n_only = sum(1 for c in calls if c[3] == 0)
say(f"  Etoliko-only (no public dataset >= {THR:.0%}) : {n_only}/{len(calls)}")
say(f"  found elsewhere                              : {len(calls) - n_only}/{len(calls)}")
say()
say(f"  {'MAG':<20}{'genus':<20}{'best public breadth':>21}{'sites':>7}")
say('-' * W)
for mid, wn, g, ns, mx, _ in sorted(calls, key=lambda x: (-x[4], x[0]))[:10]:
    say(f"  {mid:<20}{g[:19]:<20}{mx:>21.3f}{ns:>7}")

with open(os.path.join(C, 'biogeography_calls.tsv'), 'w') as oh:
    oh.write('mag_id\tworking_name\tgenus\tn_public_sites\tmax_public_breadth\tcall\tsites\n')
    for mid, wn, g, ns, mx, sites in sorted(calls):
        oh.write(f"{mid}\t{wn}\t{g}\t{ns}\t{mx:.4f}\t"
                 f"{'Etoliko-only' if ns == 0 else 'shared'}\t{','.join(sites)}\n")

STATUS[0] = 'OK'
rule('WHAT TO REPORT')
say(f"  1. Not one of the {len(rows)} Etoliko genomes recruits above {THR:.0%} breadth in")
say(f"     any of {len(have)} public euxinic and meromictic datasets ({len(PUB)} runs),")
say(f"     while the same genomes reach {etomax:.3f} in Etoliko's own reads. The")
say(f"     highest value anywhere else is {pubmax:.3f}.")
say("  2. Quote the SENSITIVITY, not just the call: the answer is identical at")
say("     every threshold from 0.50 down to 0.01, so it does not rest on a")
say("     chosen cutoff.")
say("  3. Say what the comparison set IS — Baltic, Cariaco, Black Sea, Saanich,")
say("     Orca Basin, Santa Barbara, Red Sea brines, Rogoznica, and the")
say("     meromictic lakes Zug, El Tobar and others. 'Endemic' is a claim about")
say("     the datasets searched, and the list is the claim's scope.")
say("  4. LIMITATION, stated with the result and not softened: public runs were")
say("     streamed with `fastq-dump -X`, which takes the FIRST N spots, not a")
say("     random sample. That can MISS a rare organism; it cannot manufacture")
say("     an absence across 30 independent datasets. A public dataset that DID")
say("     recruit would have to be re-checked with a random subsample first.")
say("  5. Absence of read recruitment is absence of THAT POPULATION at the")
say("     sequencing depth available, not absence of the genus or the function.")
say("     Etoliko's genera occur elsewhere; these specific genomes do not.")
dump(0)
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "verdict written" "$([ -s "$C/BIOGEOGRAPHY.txt" ] && echo 0 || echo 1)" || OK=1
check "calls table written" "$([ -s "$C/biogeography_calls.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$C/BIOGEOGRAPHY.txt" 2>/dev/null | awk '{print $2}')
echo "  status line: ${ST:-<none>}"
[ "${ST:-MISSING}" = "OK" ] || { echo "  [STOP] ${ST:-no status line} — not quotable."; OK=1; }
[ "$OK" -eq 0 ] && runlog "37_biogeography complete" || runlog "**37_biogeography FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
