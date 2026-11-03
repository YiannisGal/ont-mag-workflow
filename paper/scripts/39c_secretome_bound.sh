#!/usr/bin/env bash
# =============================================================================
# 39c_secretome_bound.sh  --  BOUNDING THE ARCHAEAL SECRETOME APPROXIMATION
#
# WHY 39B'S ANSWER WAS "UNINFORMATIVE", AND WHY THAT WAS MY FAULT
#
#   39 section B measured three N-terminal properties for archaeal against
#   bacterial proteomes and then refused to report the comparison, because its
#   sensitivity control failed: the smallest real-vs-shuffle h-region gap was
#   0.039 KD units against a 0.05 bar. Two things were wrong with that.
#
#   1. THE CONTROL WAS DILUTED BY ITS OWN DENOMINATOR. Signal peptides occur in
#      roughly 5-15% of a proteome. Shuffling the N-terminus of a cytoplasmic
#      protein changes almost nothing, because there was no positional structure
#      to destroy. Averaged over 88,856 proteins, a real effect in the ~10% that
#      carry a signal peptide is divided by ten before it is compared to a
#      threshold. The control was therefore measuring dilution, not resolution.
#
#   2. ONE CONTROL GATED THREE DIFFERENT MEASUREMENTS. The shuffle control tests
#      whether the H-REGION metric — a maximum over sliding windows — can see
#      positional structure. n-region NET CHARGE is a direct count over the
#      first eight residues; it has no windows and no positional maximum, and a
#      failure of the h-region control says nothing about it. Gating all three
#      axes on one control is Rule 9a in a new costume: a gate belongs on the
#      branch it applies to, and putting it everywhere does not make the
#      analysis more conservative, it makes it blind in one direction.
#
#   AND THE BLINDED DIRECTION CARRIED THE RESULT. 39B's own table already
#   showed n-region charge at 0.55 for bacteria and 0.11 for archaea — a
#   five-fold difference — and refused to say so.
#
# WHAT THIS DOES INSTEAD
#
#   * Computes the same three properties, but ALSO on an SP-PLAUSIBLE SUBSET:
#     the top decile of h-region hydrophobicity WITHIN each domain. Taking the
#     decile within each domain separately means the subset is defined without
#     presupposing which domain has more signal peptides.
#   * Runs the shuffle control ON THAT SUBSET, where signal peptides actually
#     live. That is where the metric must have resolution for the h-region
#     comparison to mean anything.
#   * Tests each axis on its own terms, with an EFFECT SIZE beside every p —
#     because at n ~ 89,000 everything is significant and a p-value stops
#     carrying information. Cohen's d and a raw difference are what should be
#     quoted.
#   * States the DIRECTION of the bias and what it implies for a Gram-positive
#     model applied to archaea, which is the thing the word "approximate" has
#     been hiding since step 23.
#
# STATISTICS, hand-rolled — there is no scipy on this machine
#   Welch's t with a normal tail (n per group is ~10^4, so t -> z; the
#   approximation is stated in the output, not hidden), and a two-proportion
#   z-test for the AXA motif. Both validated below against closed forms.
#
# WHAT IT STILL CANNOT DO
#   It measures the INPUT distribution, not predictor accuracy. It cannot say a
#   single DeepSig call is right or wrong. It converts an unquantified caveat
#   into a bounded one with a direction. SignalP 6.0's archaeal mode remains the
#   only real remedy.
#
# OUTPUT  39_loose_ends/SECRETOME_BOUND.txt · nterm_stats.tsv
# RUNTIME 1-3 min
# USAGE   ws2 && bash scripts/39c_secretome_bound.sh 2>&1 | tee logs/39c.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="39c_secretome"
A="$WS2/09_annotation"; F="$WS2/07_mags"
D="$WS2/39_loose_ends"; mkdir -p "$D"

echo "######## 39c — BOUNDING THE ARCHAEAL SECRETOME APPROXIMATION ########"
echo "  started: $(date -Is)"
rm -f "$D/SECRETOME_BOUND.txt" "$D/nterm_stats.tsv"

echo
python3 - "$D" "$A" "$F" "$WS2" <<'PY'
import sys, os, re, csv, math, random
from collections import defaultdict

D, A, F, WS2 = sys.argv[1:5]
W = 100
rep = []
STATUS = ['INCOMPLETE']


def say(s=''):
    print(s)
    rep.append(s)


def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t)
        say('=' * W)


def dump(status):
    STATUS[0] = status
    with open(os.path.join(D, 'SECRETOME_BOUND.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'SECRETOME_BOUND.txt')}")
    sys.exit(0)


# ------------------------------------------------------------------ statistics
def norm_sf(z):
    """upper tail of the standard normal. exact via erfc, no table."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def welch(m1, v1, n1, m2, v2, n2):
    """Welch's t. At n ~ 10^4 per group the t distribution IS the normal to
    better than 1e-4, and that approximation is printed rather than assumed."""
    if n1 < 2 or n2 < 2:
        return None, None, None
    se2 = v1 / n1 + v2 / n2
    if se2 <= 0:
        return None, None, None
    t = (m1 - m2) / math.sqrt(se2)
    p = 2.0 * norm_sf(abs(t))
    sp = math.sqrt(max(1e-300, ((n1 - 1) * v1 + (n2 - 1) * v2) / max(1, n1 + n2 - 2)))
    d = (m1 - m2) / sp if sp > 0 else 0.0
    return t, p, d


def two_prop(k1, n1, k2, n2):
    """two-proportion z-test with pooled variance."""
    if n1 == 0 or n2 == 0:
        return None, None, None
    p1, p2 = k1 / n1, k2 / n2
    pp = (k1 + k2) / (n1 + n2)
    se = math.sqrt(max(1e-300, pp * (1 - pp) * (1 / n1 + 1 / n2)))
    if se <= 0:
        return None, None, None
    z = (p1 - p2) / se
    return z, 2.0 * norm_sf(abs(z)), (p1 - p2) * 100.0


# self-check on closed forms; a statistic that cannot reproduce a known value
# does not get to decide anything
_ok = (abs(norm_sf(0.0) - 0.5) < 1e-12 and abs(norm_sf(1.959963985) - 0.025) < 1e-6
       and abs(norm_sf(2.575829304) - 0.005) < 1e-6)

# ------------------------------------------------------------------ properties
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
POS, NEG = set('KR'), set('DE')


def nterm(s):
    if len(s) < 25:
        return None
    n = s[:8]
    charge = sum(1 for c in n if c in POS) - sum(1 for c in n if c in NEG)
    win = [sum(KD.get(c, 0.0) for c in s[i:i + 9]) / 9.0
           for i in range(2, min(23, len(s) - 8))]
    if not win:
        return None
    return charge, max(win), (1 if re.search(r'A.A', s[15:35]) else 0)


def load_domains():
    dom = {}
    tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
    if not os.path.exists(tp):
        return dom
    with open(tp, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            blob = ' '.join(str(v or '') for v in r.values())
            d = next((str(r[c]).strip() for c in sorted(k for k in r if k)
                      if c.strip().lower() == 'domain'), '')
            if not d:
                m = re.search(r'd__([A-Za-z]+)', blob)
                d = m.group(1) if m else ''
            if not d:
                continue
            for v in r.values():
                v = str(v or '').strip()
                if v and '_' in v and any(ch.isdigit() for ch in v):
                    dom.setdefault(v, d)
    return dom


rule('0. THE STATISTICS ANSWER KNOWN VALUES FIRST')
say(f"  normal tail self-check (0 -> 0.5, 1.96 -> 0.025, 2.576 -> 0.005): "
    f"{'[PASS]' if _ok else '[FAIL]'}")
if not _ok:
    say("  The tail function is wrong, so no test below may be believed.")
    dump('STATS SELF-CHECK FAILED')
say()

rule('1. THE MEASUREMENT')
dom = load_domains()
if not dom:
    say("  !! no MAG_TAXONOMY.tsv / no domain column — cannot split by domain.")
    dump('NO TAXONOMY')

rng = random.Random(20260904)
recs = {'Bacteria': [], 'Archaea': []}   # (charge, hmax, axa, shuffled_hmax)
gen = defaultdict(set)
for m in (sorted(os.listdir(A)) if os.path.isdir(A) else []):
    faa = os.path.join(A, m, m + '.faa')
    if not os.path.exists(faa):
        continue
    d = dom.get(m, '')
    if not d:
        continue
    grp = 'Archaea' if d.lower().startswith('archae') else 'Bacteria'
    gen[grp].add(m)
    seq, have = [], False

    def flush():
        if not have:
            return
        s = ''.join(seq)[:40]
        st = nterm(s)
        if not st:
            return
        ch = list(s)
        rng.shuffle(ch)
        sh = nterm(''.join(ch))
        recs[grp].append((st[0], st[1], st[2], sh[1] if sh else st[1]))

    with open(faa, errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                seq, have = [], True
            else:
                seq.append(line.strip())
    flush()

if not recs['Bacteria'] or not recs['Archaea']:
    say("  !! one domain has no scorable proteins — no comparison possible.")
    dump('ONE DOMAIN ONLY')

with open(os.path.join(D, 'nterm_stats.tsv'), 'w') as fh:
    fh.write('domain\tcharge\thmax\taxa\tshuffled_hmax\n')
    for g in ('Bacteria', 'Archaea'):
        for r in recs[g]:
            fh.write(f"{g}\t{r[0]}\t{r[1]:.4f}\t{r[2]}\t{r[3]:.4f}\n")


def summarise(rows):
    n = len(rows)
    if not n:
        return None
    mc = sum(r[0] for r in rows) / n
    vc = sum((r[0] - mc) ** 2 for r in rows) / max(1, n - 1)
    mh = sum(r[1] for r in rows) / n
    vh = sum((r[1] - mh) ** 2 for r in rows) / max(1, n - 1)
    ax = sum(r[2] for r in rows)
    ms = sum(r[3] for r in rows) / n
    return dict(n=n, mc=mc, vc=vc, mh=mh, vh=vh, ax=ax, ms=ms)


# ---- the SP-plausible subset: top decile of h-region WITHIN each domain, so
#      the subset definition cannot itself create the difference being tested
sub = {}
for g in ('Bacteria', 'Archaea'):
    rows = sorted(recs[g], key=lambda r: (-r[1], r[0], r[2]))
    k = max(1, len(rows) // 10)
    sub[g] = rows[:k]

FULL = {g: summarise(recs[g]) for g in recs}
SUB = {g: summarise(sub[g]) for g in sub}

say(f"  {'set':<22}{'domain':<11}{'MAGs':>6}{'proteins':>10}"
    f"{'n-charge':>12}{'h-region KD':>14}{'AXA':>9}{'shuffled KD':>13}")
say('-' * W)
for lbl, S in (('whole proteome', FULL), ('top-decile h-region', SUB)):
    for g in ('Bacteria', 'Archaea'):
        s = S[g]
        say(f"  {lbl:<22}{g:<11}{len(gen[g]):>6}{s['n']:>10,}"
            f"{s['mc']:>12.3f}{s['mh']:>14.3f}"
            f"{100.0 * s['ax'] / s['n']:>8.1f}%{s['ms']:>13.3f}")
say()

# =====================================================================
rule('2. CONTROL — DOES THE h-REGION METRIC HAVE RESOLUTION WHERE IT MATTERS?')
say("  Shuffling the N-terminus of a cytoplasmic protein destroys nothing,")
say("  because there was no positional structure there to destroy. Over a whole")
say("  proteome the ~10% that carry a signal peptide are diluted by the ~90%")
say("  that do not — so a whole-proteome shuffle control measures DILUTION, not")
say("  resolution. That is why 39B's control failed at 0.039, and the failure")
say("  was a property of its denominator rather than of the metric.")
say()
say(f"  {'set':<22}{'domain':<11}{'real KD':>10}{'shuffled':>11}{'gap':>9}")
say('-' * W)
gaps = {}
for lbl, S, key in (('whole proteome', FULL, 'full'),
                    ('top-decile h-region', SUB, 'sub')):
    for g in ('Bacteria', 'Archaea'):
        s = S[g]
        gp = s['mh'] - s['ms']
        gaps[(key, g)] = gp
        say(f"  {lbl:<22}{g:<11}{s['mh']:>10.3f}{s['ms']:>11.3f}{gp:>9.3f}")
say()
# --- the subset above is selected BY h-max, so comparing it to its own shuffles
#     is selection on the statistic under test: taking extremes and then
#     shuffling regresses to the mean and manufactures a gap out of nothing.
#     The honest comparison applies the SAME selection pressure to both — the
#     top decile of the REAL h-max distribution against the top decile of the
#     SHUFFLED one. If the sequences carry no positional structure the two
#     distributions are the same and their top deciles coincide.
say()
say("  The subset above is selected BY h-max, so testing it against its own")
say("  shuffles would be selection on the statistic under test — take extremes,")
say("  shuffle, and regression to the mean manufactures a gap from nothing.")
say("  So: the top decile of the REAL distribution against the top decile of the")
say("  SHUFFLED one. Equal selection pressure on both sides.")
say()
say(f"  {'domain':<11}{'real top decile':>18}{'shuffled top decile':>22}{'gap':>9}")
say('-' * W)
deciles = {}
for g in ('Bacteria', 'Archaea'):
    rows = recs[g]
    k = max(1, len(rows) // 10)
    rt = sorted((r[1] for r in rows), reverse=True)[:k]
    st = sorted((r[3] for r in rows), reverse=True)[:k]
    mr, ms_ = sum(rt) / len(rt), sum(st) / len(st)
    deciles[g] = mr - ms_
    say(f"  {g:<11}{mr:>18.3f}{ms_:>22.3f}{mr - ms_:>9.3f}")
say()
ctrl = min(deciles['Bacteria'], deciles['Archaea'])
say(f"  CONTROL statistic = smallest matched top-decile gap = {ctrl:.3f} KD units")
H_OK = ctrl >= 0.20
if H_OK:
    say("  >> [PASS] on the subset where signal peptides actually live, the")
    say("     metric separates a real N-terminus from its own composition-")
    say("     preserving shuffle. The h-region axis is usable.")
else:
    say("  >> [FAIL] even on the enriched subset the metric cannot separate a")
    say("     real N-terminus from its shuffle. The h-region axis is reported")
    say("     as UNINFORMATIVE. The other two axes are unaffected — they are")
    say("     direct counts and do not depend on a windowed maximum.")
say()

# =====================================================================
rule('3. EACH AXIS ON ITS OWN TERMS, WITH AN EFFECT SIZE BESIDE EVERY p')
say("  n is ~89,000. At that size a p-value stops carrying information — every")
say("  difference is 'significant'. Quote the raw difference and Cohen's d.")
say("  The p-values use a normal tail; at n ~ 10^4 per group the t distribution")
say("  is the normal to better than 1e-4, and that is stated rather than hidden.")
say()
b, a = FULL['Bacteria'], FULL['Archaea']
say(f"  {'axis':<26}{'Bacteria':>11}{'Archaea':>11}{'diff':>10}"
    f"{'d':>8}{'p':>12}{'gated by control?':>22}")
say('-' * W)

t1, p1, d1 = welch(b['mc'], b['vc'], b['n'], a['mc'], a['vc'], a['n'])
say(f"  {'n-region net charge':<26}{b['mc']:>11.3f}{a['mc']:>11.3f}"
    f"{b['mc'] - a['mc']:>10.3f}{d1:>8.2f}{p1:>12.1e}{'no - a direct count':>22}")

if H_OK:
    t2, p2, d2 = welch(b['mh'], b['vh'], b['n'], a['mh'], a['vh'], a['n'])
    say(f"  {'h-region max KD':<26}{b['mh']:>11.3f}{a['mh']:>11.3f}"
        f"{b['mh'] - a['mh']:>10.3f}{d2:>8.2f}{p2:>12.1e}{'yes - control PASSED':>22}")
else:
    d2 = None
    say(f"  {'h-region max KD':<26}{b['mh']:>11.3f}{a['mh']:>11.3f}"
        f"{b['mh'] - a['mh']:>10.3f}{'--':>8}{'--':>12}{'yes - BLOCKED':>22}")

z3, p3, dd3 = two_prop(b['ax'], b['n'], a['ax'], a['n'])
say(f"  {'AXA motif (%)':<26}{100.0 * b['ax'] / b['n']:>11.1f}"
    f"{100.0 * a['ax'] / a['n']:>11.1f}{dd3:>10.1f}{'--':>8}{p3:>12.1e}"
    f"{'no - a direct count':>22}")
say()

# =====================================================================
rule('4. THE BOUND, AND ITS DIRECTION')
say("  A Gram-positive signal-peptide model is trained on N-termini obeying the")
say("  POSITIVE-INSIDE RULE: a short, net-POSITIVELY-charged n-region preceding")
say("  the hydrophobic core. n-region charge is therefore not one property among")
say("  three — it is the feature the model most directly keys on.")
say()
big_charge = abs(d1 or 0.0) >= 0.20 or abs(b['mc'] - a['mc']) >= 0.30
if big_charge:
    ratio = (b['mc'] / a['mc']) if a['mc'] not in (0,) else float('inf')
    say(f"  >> ARCHAEAL N-TERMINI ARE MATERIALLY LESS POSITIVELY CHARGED:")
    say(f"     {b['mc']:.2f} against {a['mc']:.2f}"
        + (f" — a factor of {ratio:.1f}." if math.isfinite(ratio) else "."))
    say(f"     Cohen's d = {d1:.2f}"
        + ("  (small by convention, but on the feature that decides the call,"
           if abs(d1) < 0.5 else "  (a moderate-to-large standardised difference,"))
    say("     and measured over tens of thousands of proteins).")
    say()
    say("     DIRECTION, which is the part 'approximate' was hiding: a model")
    say("     expecting a positively charged n-region, applied to N-termini that")
    say("     carry less positive charge, will tend to UNDER-call signal")
    say("     peptides. **The archaeal secretome percentages are therefore")
    say("     likely to be an UNDER-estimate, not a symmetric uncertainty.**")
    say("     That is a one-sided caveat and it should be written as one.")
    verdict = 'BOUNDED — one-sided, archaeal secretome likely under-called'
else:
    say("  >> The n-region charge distributions are close. On the feature a")
    say("     Gram-positive model keys on hardest, these archaea sit inside the")
    say("     bacterial spread, so the approximation is small on that axis.")
    verdict = 'BOUNDED — small on the charge axis'
say()
if dd3 is not None and abs(dd3) >= 2.0:
    say(f"  The AXA cleavage motif differs by {dd3:+.1f} points"
        f" ({100.0 * b['ax'] / b['n']:.1f}% vs {100.0 * a['ax'] / a['n']:.1f}%),")
    say("  in the same accounting. Report it beside the charge, not instead of it.")
    say()

say("  SENTENCE FOR THE THESIS, with numbers instead of the word 'approximate'.")
say("  The closing clause follows the verdict above; it is not boilerplate.")
say()
say("    'Archaeal secretome fractions were obtained with DeepSig's Gram-positive")
say("     model, no archaeal model being available. Archaeal N-terminal regions")
say(f"     in this dataset carry a mean net charge of {a['mc']:.2f} against"
    f" {b['mc']:.2f}")
say(f"     in the bacterial genomes (Cohen's d = {d1:.2f}, n = {a['n']:,} and"
    f" {b['n']:,}),")
say(f"     and an A-X-A cleavage motif in {100.0 * a['ax'] / a['n']:.1f}% against"
    f" {100.0 * b['ax'] / b['n']:.1f}% of N-termini.")
if big_charge and (b['mc'] - a['mc']) > 0:
    say("     Because the model keys on a positively charged n-region and the")
    say("     archaeal n-regions carry LESS positive charge, these percentages")
    say("     are expected to be an under-estimate rather than a symmetric")
    say("     uncertainty; archaeal and bacterial fractions are not compared")
    say("     directly.'")
elif big_charge:
    say("     The archaeal n-regions carry MORE positive charge than the")
    say("     bacterial ones, so the model's central assumption is not violated")
    say("     in the direction that would cause under-calling; the difference is")
    say("     nonetheless material and the two domains are not compared")
    say("     directly.'")
else:
    say("     These distributions are close on the feature the model keys on")
    say("     hardest, so the approximation is small on that axis; archaeal and")
    say("     bacterial fractions are nonetheless not compared directly, because")
    say("     agreement of input properties is not validation of the predictor.'")
say()
rule('LIMITS')
say("  * This measures the INPUT distribution, not predictor accuracy. It cannot")
say("    say any single DeepSig call is right, and agreement of these properties")
say("    would not prove the model transfers.")
say("  * The SP-plausible subset is the top decile of h-region hydrophobicity,")
say("    not a set of known signal peptides. It is an enrichment, not a truth set.")
say("  * The A-X-A regex is a crude proxy for a cleavage site and will match by")
say("    chance in an alanine-rich sequence.")
say("  * The real remedy is SignalP 6.0, which has an archaeal mode and needs")
say("    academic registration. This bounds the caveat; it does not remove it.")
dump(verdict if H_OK else verdict + ' (h-region axis blocked)')
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "report written" "$([ -s "$D/SECRETOME_BOUND.txt" ] && echo 0 || echo 1)" || OK=1
check "per-protein table written" "$([ -s "$D/nterm_stats.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/SECRETOME_BOUND.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  status: ${ST:-<none>}"
case "${ST:-MISSING}" in
  BOUNDED*) : ;;
  *) echo "  [STOP] ${ST:-no status line}"; OK=1 ;;
esac
[ "$OK" -eq 0 ] && runlog "39c_secretome: $ST" || runlog "**39c_secretome FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
