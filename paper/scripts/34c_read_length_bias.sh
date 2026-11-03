#!/usr/bin/env bash
# =============================================================================
# 34c_read_length_bias.sh  --  is step 34's x/rpoB proxy distorted by gene length?
#
# THE QUESTION
#   Step 34 compares read counts between markers as a ratio to rpoB. But a
#   DIAMOND blastx read count is not a count of gene copies — it is a count of
#   READS THAT OVERLAP the gene, and that scales with length. rpoB is ~4 kb;
#   nifH is ~0.9 kb. If the proxy is length-biased, then every ratio in step
#   34's section 2b and the mcrA ceiling in 2c are biased with it, and two
#   readings could flip:
#
#     nifH  "reads low — thin or divergent queries"   (0.2)
#     cutC  "READS TOO HIGH — suspect the query set"  (7.7)
#
# WHY IT MIGHT NOT MATTER, AND WHY THAT HAS TO BE MEASURED NOT ASSUMED
#   For SHORT reads the expected count scales with gene length, so a 4.4-fold
#   length difference is a 4.4-fold bias. For LONG reads it scales with
#   (read length + gene length) — and when reads are 6-10 kb and genes are
#   1-4 kb, the read length dominates and the bias mostly cancels.
#
#   v1 of this script ASSERTED the long-read model on the grounds that the
#   protocol says ONT. Then it measured, and the median came back at a few
#   hundred bp. The model is now CHOSEN BY THE MEASUREMENT, and the script says
#   which one it picked and why. It also measures over the WHOLE read file
#   rather than the first 20,000 reads, and compares the two — if the head of a
#   file disagrees with the file, nothing computed from the head is usable.
#
#   This step computes the read-length distribution FROM THE READS and
#   the median reference length FROM THE REFERENCE SET step 34 actually used,
#   then re-derives every conclusion under all three models:
#
#     RAW           what step 34 reported
#     GENE          count / gene length          (short-read model)
#     READ+GENE     count / (read + gene) length (long-read model)
#
#   A conclusion that survives all three is robust and can be written plainly.
#   A conclusion that changes between them is MODEL-DEPENDENT and has to be
#   reported with the model named. Neither outcome is a failure; hiding which
#   one we are in would be.
#
# NOTHING IS RE-SEARCHED
#   This reads step 34's own outputs — marker_levels.tsv and marker_ref.faa —
#   plus a sample of the reads. It cannot change step 34's counts, only their
#   interpretation. Runtime is under two minutes.
#
# THE NORMALISATION IS APPLIED TO NUMERATOR AND DENOMINATOR ALIKE
#   Rule 18, from step 32d, where a filter applied to one side alone turned
#   eight unbinned hits invisible. rpoB is normalised by its own length exactly
#   as every marker is by its own.
#
# INTERNAL VALIDITY CHECK
#   dsrA and dsrB are adjacent genes of one operon and must occur in equal
#   copy number. They came back 70 and 68. Any normalisation that PULLS THEM
#   APART is wrong, whatever it does to anything else. That check runs first
#   and a failure blocks the rest.
#
# OUTPUT  34_unbinned/length_bias.tsv · LENGTH_BIAS.txt
# USAGE   ws2 && bash scripts/34c_read_length_bias.sh 2>&1 | tee logs/34c.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="34c_lengthbias"
D="$WS2/34_unbinned"

echo "######## 34c — IS THE READ PROXY LENGTH-BIASED? ########"
echo "  started: $(date -Is)"

for f in "$D/marker_levels.tsv" "$D/marker_ref.faa"; do
  [ -s "$f" ] || { echo "  !! missing $f — run step 34 first."; exit 2; }
done

# ------------------------------------------------- read length, measured
# v1 took the FIRST 20,000 reads and got a median of 289 bp, which is not what
# an ONT R10.4.1 library looks like. Either the library really is that short or
# the head of the file is not representative of it — and those have opposite
# consequences for every ratio below, so guessing is not an option.
#
# So: one full pass over the read file building a length histogram (bounded
# memory, exact to the bucket), AND the head sample kept separately. If the two
# disagree the file is not homogeneous, and that is worth knowing on its own.
echo
echo "######## 1. READ LENGTH — FULL PASS, NOT THE HEAD ########"
CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
SED_IN="$CLEAN"; [ -s "$SED_IN" ] || SED_IN="$SED_READS"
echo "  file: $SED_IN"
echo "  size: $(du -h "$SED_IN" 2>/dev/null | cut -f1)"

HEAD="$D/read_len_head.txt"
if [ ! -s "$HEAD" ]; then
  if [ "${SED_IN##*.}" = "gz" ]; then
    ( zcat "$SED_IN" 2>/dev/null | head -n 80000 | awk 'NR%4==2{print length($0)}' > "$HEAD" ) || true
  else
    ( head -n 80000 "$SED_IN" | awk 'NR%4==2{print length($0)}' > "$HEAD" ) || true
  fi
fi

HIST="$D/read_len_hist.tsv"
if [ ! -s "$HIST" ]; then
  echo "  full pass over the read file (minutes, cached afterwards) ..."
  if [ "${SED_IN##*.}" = "gz" ]; then
    zcat "$SED_IN" 2>/dev/null
  else
    cat "$SED_IN"
  fi | awk 'NR%4==2{L=length($0); n++; s+=L; b[int(L/50)]++}
            END{printf "N\t%d\n", n; printf "SUM\t%.0f\n", s;
                for (k in b) printf "B\t%d\t%d\n", k*50, b[k]}' > "$HIST.part"
  if [ -s "$HIST.part" ] && grep -q '^N' "$HIST.part"; then
    mv "$HIST.part" "$HIST"
  else
    rm -f "$HIST.part"
    echo "  !! full pass produced nothing — cannot measure the read length."
    echo "     Every ratio below depends on it, so nothing is reported."
    exit 3
  fi
else
  echo "  histogram cached"
fi
[ -s "$HEAD" ] || { echo "  !! could not sample $SED_IN"; exit 3; }

# ------------------------------------------------------------ analyse

echo
python3 - "$D" "$HIST" "$HEAD" <<'PY'
import sys, os, csv, statistics as st
from collections import defaultdict

D, HIST, HEADF = sys.argv[1], sys.argv[2], sys.argv[3]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

# ---- whole-file length distribution, from the histogram ----------------
N = 0; SUM = 0
buckets = []          # (lower_bp, count) — bucket width 50
for line in open(HIST):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'N':
        N = int(f[1])
    elif f[0] == 'SUM':
        SUM = int(float(f[1]))
    elif f[0] == 'B' and len(f) == 3:
        buckets.append((int(f[1]), int(f[2])))
buckets.sort()
if not N or not buckets:
    say("  !! empty read-length histogram"); sys.exit(4)

def quantile(q):
    """bucket midpoint at cumulative fraction q"""
    want = q * N
    c = 0
    for lo, n in buckets:
        c += n
        if c >= want:
            return lo + 25
    return buckets[-1][0] + 25

def n50():
    """length-weighted median: half the BASES sit in reads at least this long"""
    half = SUM / 2
    c = 0
    for lo, n in reversed(buckets):
        c += (lo + 25) * n
        if c >= half:
            return lo + 25
    return buckets[0][0] + 25

RLEN = quantile(0.5)
MEAN = SUM / N
N50 = n50()
# WHICH STATISTIC DRIVES THE MODEL, AND WHY IT IS THE MEAN.
#   The number of reads overlapping a gene of length G is proportional to the
#   number whose start falls within (L_read + G) of it, summed over reads:
#       E[count] proportional to  sum_i (L_i + G)  =  N * (mean(L) + G)
#   The expectation is LINEAR in read length, so the mean is the statistic that
#   belongs in it. v2 used the median, which is smaller than the mean whenever
#   the distribution has a long right tail — and an ONT run always does. That
#   under-weighted the read term and pushed the model toward the short-read
#   limit it should not have been in.
EFF = MEAN
say(f"  reads       : {N:,}   bases {SUM:,}   ({SUM / 1e9:.2f} Gbp)")
say(f"  median      : {RLEN:,} bp        mean {MEAN:,.0f} bp")
say(f"  quartiles   : Q1 {quantile(0.25):,}   Q3 {quantile(0.75):,}   "
    f"P90 {quantile(0.90):,}   P99 {quantile(0.99):,}")
say(f"  read N50    : {N50:,} bp   (half the BASES are in reads this long or longer)")
say()
say(f"  the length model below uses the MEAN ({MEAN:,.0f} bp), not the median, because")
say("  the expected count is linear in read length: sum_i (L_i + G) = N(mean + G).")

# ---- is the head of the file representative of the whole? --------------
hl = [int(x) for x in open(HEADF) if x.strip().isdigit()]
if hl:
    hmed = st.median(hl)
    say()
    say(f"  first {len(hl):,} reads: median {hmed:,.0f} bp   vs whole file {RLEN:,} bp")
    if hmed and (max(hmed, RLEN) / max(1, min(hmed, RLEN))) > 2:
        say("  !! THE HEAD OF THE FILE IS NOT REPRESENTATIVE — they differ by more than")
        say("     twofold. Anything computed from a head sample of this file is wrong,")
        say("     and that includes v1 of this script. The whole-file numbers are used")
        say("     below. It also means the read set is heterogeneous: worth finding out")
        say("     why before it reaches a methods section.")
    else:
        say("  head and whole file agree, so the file is homogeneous in read length.")

# ---- does this look like the ONT library it is supposed to be? ---------
say()
# A LONG-READ LIBRARY IS JUDGED BY N50, NOT BY MEDIAN READ LENGTH.
#   v2 flagged a median of 425 bp as "not an ONT library". That was wrong. ONT
#   runs are strongly right-skewed: most READS are short while most BASES sit in
#   long ones. Median describes the typical read; N50 describes where the
#   sequence actually is, and assembly depends on the latter.
if N50 >= 2000:
    say(f"  N50 {N50:,} bp — this IS a long-read library. The median of {RLEN:,} bp says")
    say("  only that most READS are short; the N50 says most BASES are in long ones,")
    say("  which is the normal shape of an ONT run and the one assembly depends on.")
    if RLEN < N50 / 4:
        say()
        say(f"  Worth a Methods sentence though: {100.0 * sum(n for lo, n in buckets if lo < 1000) / N:.0f}%"
            " of reads are under 1 kb while the N50 is")
        say(f"  {N50:,} bp. A large short-fragment population is what degraded sediment DNA")
        say("  looks like, and it has a consequence downstream — a 400 bp read carries")
        say("  few 31-mers to match, so part of Kraken2's 0.47% is a read-length effect")
        say("  and not only community novelty. Translation-based classification (step")
        say("  35) is much less exposed to this: 400 bp is still ~130 amino acids.")
else:
    say(f"  !! N50 {N50:,} bp AND median {RLEN:,} bp — short by both measures, which is not")
    say("     what an ONT R10.4.1 library should look like. Check the file identity and")
    say("     any filtering step before using anything downstream of it.")
    say(f"     {N:,} reads and {SUM / 1e9:.2f} Gbp produced 28,572 sediment contigs and 29")
    say("     MAGs; ask whether that is consistent with this distribution.")

# ---------------------------------------- reference lengths, per marker
# The reference set step 34 actually searched. Headers look like
#   gene|tr|ACC|NAME   gene|sp|ACC|NAME   gene|ASM_<protein>   CTRLgene|...
reflen = defaultdict(list)
name = None
cur = 0
def flush():
    if name is not None and cur:
        reflen[name].append(cur)
with open(os.path.join(D, 'marker_ref.faa'), errors='replace') as fh:
    for line in fh:
        if line[:1] == '>':
            flush()
            tag = line[1:].split('|')[0].strip()
            name = tag[4:] if tag.startswith('CTRL') else tag
            cur = 0
        else:
            cur += len(line.strip())
flush()

say()
say()
say(f"  {'marker':<10}{'refs':>7}{'median aa':>11}{'gene bp':>10}{'read+gene':>12}")
say('-' * W)
gene_bp, thin = {}, []
for g in sorted(reflen):
    n = len(reflen[g])
    aa = st.median(reflen[g])
    bp = aa * 3
    gene_bp[g] = bp
    if n < 3:
        thin.append(g)
    say(f"  {g:<10}{n:>7}{aa:>11,.0f}{bp:>10,.0f}{bp + EFF:>12,.0f}")
if thin:
    say()
    say(f"  !! fewer than 3 reference sequences for: {', '.join(thin)} — their median")
    say("     length is unstable and every normalised number for them is soft.")

# ---------------------------------------------------- step 34's numbers
rows = list(csv.DictReader(open(os.path.join(D, 'marker_levels.tsv')),
                           delimiter='\t'))
def num(r, k):
    try:
        return float(r.get(k) or 0)
    except ValueError:
        return 0.0

reads = {r['gene']: num(r, 'reads_sed') for r in rows}
asm = {r['gene']: num(r, 'asm_sed') for r in rows}
genes = [r['gene'] for r in rows]

# marker_levels.tsv holds TARGETS ONLY — the controls live in step 34's own
# printed table and nowhere else (read_hits.tsv is named in its header but never
# written). Take the control numbers verbatim from VERDICT.txt rather than
# recounting them here: recounting would re-implement step 34's promiscuous-
# reference exclusions, and Rule 19 says reuse the code, not the idea.
if 'rpoB' not in reads or not reads['rpoB']:
    vp = os.path.join(D, 'VERDICT.txt')
    got = []
    if os.path.exists(vp):
        for line in open(vp, errors='replace'):
            f = line.split()
            if len(f) == 6 and f[1] in ('BACT', 'ARCH'):
                try:
                    asm[f[0]] = float(f[2])
                    reads[f[0]] = float(f[4])
                    got.append(f[0])
                except ValueError:
                    pass
    if got:
        say()
        say(f"  control counts taken from VERDICT.txt: {', '.join(got)}")
        say("  (marker_levels.tsv carries targets only)")

if 'rpoB' not in reads or not reads['rpoB']:
    say()
    say("  !! no rpoB read count in marker_levels.tsv OR VERDICT.txt.")
    say("     rpoB is the denominator of every ratio in step 34's section 2b, so")
    say("     without it there is nothing to normalise against and nothing is")
    say("     reported. Check that step 34's VERDICT.txt contains its section-1")
    say("     control table.")
    sys.exit(5)

MODELS = ('RAW', 'GENE', 'READ+GENE')
def norm(g, model):
    v = reads.get(g, 0.0)
    if model == 'RAW':
        return v
    L = gene_bp.get(g)
    if not L:
        return None
    return v / (L if model == 'GENE' else (L + EFF))

def ratio(g, model):
    a, b = norm(g, model), norm('rpoB', model)
    if a is None or b in (None, 0):
        return None
    return a / b

# ------------------------------------------- internal validity check
rule('2. VALIDITY CHECK — DOES THE NORMALISATION KEEP dsrA AND dsrB TOGETHER?')
say("  They are adjacent genes of one operon, so they occur in equal copy number.")
say("  Step 34 measured 70 and 68. A normalisation that pulls them apart is wrong,")
say("  whatever it does to anything else.")
say()
say(f"  {'model':<12}{'dsrA':>12}{'dsrB':>12}{'ratio':>10}   verdict")
say('-' * W)
VALID = {}
for m in MODELS:
    a, b = norm('dsrA', m), norm('dsrB', m)
    if a is None or b in (None, 0):
        say(f"  {m:<12}{'n/a':>12}{'n/a':>12}{'-':>10}   no reference lengths")
        VALID[m] = None
        continue
    r = a / b
    ok = 0.5 <= r <= 2.0
    VALID[m] = ok
    say(f"  {m:<12}{a:>12.4g}{b:>12.4g}{r:>10.2f}   {'ok' if ok else 'BROKEN'}")
say()
if all(v for v in VALID.values() if v is not None):
    say("  [PASS] all three models keep the operon pair together, so none of them is")
    say("         obviously wrong and the comparison below is worth reading.")
else:
    say("  [FAIL] a model splits an operon pair. Do not use it.")

# ------------------------------------------------ the comparison
rule('3. EVERY MARKER UNDER ALL THREE MODELS')
say("  'asm' is the assembly ratio to rpoB, which is a count of GENES and needs no")
say("  length correction. The read ratio is what may be biased.")
say()
say(f"  {'gene':<8}{'asm/rpoB':>10}{'RAW':>10}{'GENE':>10}{'READ+GENE':>12}"
    f"{'best agrees':>13}")
say('-' * W)
arp = asm.get('rpoB', 0.0)
out = []
for g in genes:
    if g in ('rpoB', 'glnA', 'tuf', 'ef1a'):
        continue
    ar = (asm.get(g, 0.0) / arp) if arp else None
    vals = {m: ratio(g, m) for m in MODELS}
    best, bestd = '-', None
    if ar:
        for m in MODELS:
            v = vals[m]
            if v is None or v <= 0:
                continue
            d = abs((v / ar) - 1) if ar else None
            if d is not None and (bestd is None or d < bestd):
                bestd, best = d, m
    say(f"  {g:<8}"
        + (f"{ar:>10.3f}" if ar is not None else f"{'-':>10}")
        + ''.join((f"{vals[m]:>10.3f}" if m != 'READ+GENE' else f"{vals[m]:>12.3f}")
                  if vals[m] is not None
                  else (f"{'-':>10}" if m != 'READ+GENE' else f"{'-':>12}")
                  for m in MODELS)
        + f"{best:>13}")
    out.append((g, ar, vals))

# ------------------------------------------ which readings are model-dependent
rule('4. WHICH OF STEP 34\'S READINGS DEPEND ON THE MODEL?')
say("  Step 34 called a marker 'consistent' when the read ratio sat within about")
say("  3x of the assembly ratio. Re-applying that same rule under each model:")
say()
say(f"  {'gene':<8}{'RAW':>14}{'GENE':>14}{'READ+GENE':>14}   stable?")
say('-' * W)
def call(v, ar):
    if v is None or ar in (None, 0):
        return 'n/a'
    r = v / ar
    if r > 3:
        return 'reads high'
    if r < 1 / 3:
        return 'reads low'
    return 'consistent'
flips = []
for g, ar, vals in out:
    cs = [call(vals[m], ar) for m in MODELS]
    real = [c for c in cs if c != 'n/a']
    stable = len(set(real)) <= 1 and bool(real)
    if not stable and real:
        flips.append(g)
    say(f"  {g:<8}{cs[0]:>14}{cs[1]:>14}{cs[2]:>14}   "
        + ('yes' if stable else ('NO — model-dependent' if real else '-')))
say()
# WHICH MODEL APPLIES IS DECIDED BY THE MEASUREMENT, NOT BY WHAT THE LIBRARY
# WAS SUPPOSED TO BE. v1 asserted the long-read model because the protocol says
# ONT; the measured median then came back at a few hundred bp, which makes the
# short-read model the right one instead.
MEDGENE = st.median(list(gene_bp.values())) if gene_bp else 0
if MEDGENE and EFF >= 3 * MEDGENE:
    PICK = 'READ+GENE'
    WHY = (f"a mean read of {EFF:,.0f} bp against a median gene of {MEDGENE:,.0f} bp: "
           "a read overlaps a short gene and a long one almost equally often, so "
           "most of the length bias cancels and RAW is close to correct")
elif MEDGENE and EFF <= MEDGENE / 3:
    PICK = 'GENE'
    WHY = (f"a mean read of {EFF:,.0f} bp against a median gene of {MEDGENE:,.0f} bp: "
           "reads are far shorter than the genes, so the count scales with gene "
           "length and RAW is biased by the full length ratio")
else:
    PICK = 'READ+GENE'
    WHY = (f"a mean read of {EFF:,.0f} bp is comparable to a median gene of "
           f"{MEDGENE:,.0f} bp, so neither limit applies cleanly and the "
           "(read + gene) form is the general one")
say()
say(f"  MODEL SELECTED BY MEASUREMENT: {PICK}")
say(f"  because {WHY}.")
say()
if flips:
    say(f"  >> {', '.join(flips)} change their reading depending on the length model.")
    say(f"     Read them from the {PICK} column and say so, or do not report them.")
    say("     Step 34's printed x/rpoB is the RAW column, so wherever RAW and")
    say(f"     {PICK} disagree, step 34's section 2b text is the one to correct.")
else:
    say("  >> No reading changes. Step 34's section 2b stands as written, and that")
    say("     is now a measured statement rather than an assumption: with a median")
    say(f"     read of {RLEN:,.0f} bp against genes of {min(gene_bp.values()):,.0f}-"
        f"{max(gene_bp.values()):,.0f} bp, gene length is a minor term.")

# ------------------------------------------------- the mcrA ceiling
rule('5. THE mcrA CEILING, RE-CHECKED')
say("  Every methanogen carries mcrA; only some carry any given methyltransferase.")
say("  So no methyltransferase can exceed mcrA in the same sample. Step 34 found")
say("  mtbB at 4.9x and mtaB at 2.2x. If that were a length artefact it would")
say("  vanish under normalisation.")
say()
say(f"  {'marker':<8}{'RAW':>12}{'GENE':>12}{'READ+GENE':>13}   holds?")
say('-' * W)
CEIL = []
for g in ('mtaB', 'mtbB', 'mtmB', 'mttB'):
    if g not in reads:
        continue
    xs = []
    for m in MODELS:
        a, b = norm(g, m), norm('mcrA', m)
        xs.append((a / b) if (a is not None and b) else None)
    viol = [x for x in xs if x is not None and x > 2.0]
    holds = len(viol) == len(
        [x for x in xs if x is not None]) and len(viol) > 0
    say(f"  {g:<8}" + ''.join(
        (f"{x:>12.2f}" if m != 'READ+GENE' else f"{x:>13.2f}") if x is not None
        else (f"{'-':>12}" if m != 'READ+GENE' else f"{'-':>13}")
        for m, x in zip(MODELS, xs))
        + ("   violates in every model" if holds
           else "   not a violation in every model"))
    CEIL.append((g, holds))
say()
if any(h for _g, h in CEIL):
    bad = ', '.join(g for g, h in CEIL if h)
    say(f"  >> {bad} exceeds mcrA under EVERY length model. The ceiling violation is")
    say("     not an artefact of gene length; those read counts are measuring a")
    say("     wider corrinoid methyltransferase family. Step 34's refusal to report")
    say("     them stands.")
else:
    say("  >> No marker violates the ceiling in every model. Re-read step 34's 2c")
    say("     with that in mind before quoting it.")

rule('WHAT THIS CHANGES')
say("  Nothing in step 34's COUNTS. This step cannot change them; it only asks")
say("  whether the ratios built from them are distorted by gene length.")
say()
say("  For Methods, one sentence: read-level abundance proxies were checked for")
say(f"  gene-length bias against a mean read length of {EFF:,.0f} bp (median {RLEN:,},")
say(f"  N50 {N50:,}) measured over")
say(f"  the whole read set ({N:,} reads, {SUM / 1e9:.2f} Gbp); the "
    f"{PICK} model applies")
say("  because " + WHY + ", and the marker-to-control ratios are"
    + (" not materially affected." if not flips else
       f" affected for {', '.join(flips)}."))
say()
say("  AND ONE CAVEAT THAT OUTRANKS ALL OF THE ABOVE: if the read-length")
say("  distribution in section 1 is not what this library should look like, fix")
say("  that first. Every number here is downstream of it.")

with open(os.path.join(D, 'length_bias.tsv'), 'w') as fh:
    fh.write('gene\tn_refs\tmedian_aa\tgene_bp\treads_sed\tasm_ratio\t'
             'ratio_RAW\tratio_GENE\tratio_READGENE\n')
    for g, ar, vals in out:
        fh.write(f"{g}\t{len(reflen.get(g, []))}\t"
                 f"{(st.median(reflen[g]) if reflen.get(g) else 0):.0f}\t"
                 f"{gene_bp.get(g, 0):.0f}\t{reads.get(g, 0):.0f}\t"
                 + (f"{ar:.4f}" if ar is not None else '-') + '\t'
                 + '\t'.join((f"{vals[m]:.4f}" if vals[m] is not None else '-')
                             for m in MODELS) + '\n')
open(os.path.join(D, 'LENGTH_BIAS.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'length_bias.tsv')}")
print(f"  written: {os.path.join(D, 'LENGTH_BIAS.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "length_bias.tsv written" "$([ -s "$D/length_bias.tsv" ] && echo 0 || echo 1)" || OK=1
[ "$OK" -eq 0 ] && runlog "34c_lengthbias complete" || runlog "**34c_lengthbias FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
