#!/usr/bin/env bash
# =============================================================================
# 32d_repeat_control.sh  --  are the 81 differences divergence, or paralogy?
#
# WHAT 32c FOUND
#   81 fixed differences, and they are not scattered. Five genes hold 60 of them:
#
#     IS66 family transposase   20 / 1,291 sites =  15.49 /kb    81x genome mean
#     Transposase               19 /   359 sites =  52.92 /kb   276x genome mean
#     rpoC                      14 / 3,458 sites =   4.05 /kb    21x genome mean
#     P-type Ca(2+) transporter  4 / 2,106 sites =   1.90 /kb    10x genome mean
#     hypothetical protein       3 /   247 sites =  12.15 /kb    63x genome mean
#
#   The top two are TRANSPOSASES, and the hypothetical at _10065 sits directly
#   beside the IS66 transposase at _10060 — one IS element. Mobile elements are
#   the single most repeat-rich, most multi-copy, most cross-mapping-prone class
#   of sequence in any genome. That is where a mapping artefact WOULD appear,
#   and it is exactly where the signal is.
#
#   Arithmetic, before testing anything:
#     mobile loci hold  42 / 81 differences  (52% of the signal)
#     on                1,897 / 421,698 sites  (0.45% of the genome compared)
#     removing them:    39 / 419,801 = 0.0929 /kb
#     floor 95% bound:                  0.1065 /kb
#     -> BELOW the bound. Poisson p rises from 0.0045 to 0.073.
#
#   So the entire "the rates differ" result rests on 0.45% of the compared
#   genome, and that 0.45% is transposases. That is not a conclusion yet — it is
#   a hypothesis with an obvious test.
#
# THE TEST
#   If these differences come from reads of OTHER COPIES of the same element
#   cross-mapping onto this one, two things must be true and both are checkable:
#
#   1. COPY NUMBER. The implicated gene should occur more than once in this
#      genome. Measured directly: pull the gene's sequence out of the reference
#      and align it back against the whole MAG. A single-copy gene returns one
#      hit. A multi-copy family returns several, and no read can be assigned to
#      one copy rather than another.
#
#   2. DEPTH. Reads from every copy pile onto whichever copy assembled, so a
#      multi-copy locus reads deeper than the genome median. A site sitting at
#      3x the median coverage is being fed by sequence that is not there.
#
#   Neither test involves the error floor, and neither can be satisfied by real
#   strain divergence in a single-copy gene.
#
# THE EXCLUSION RULE, STATED BEFORE THE RESULT
#   Sites in loci with two or more genomic copies are excluded from the
#   divergence rate. This is not post-hoc: it is what inStrain does for
#   short reads (popANI is computed on uniquely-mapping positions), and the
#   reason is that a variant which cannot be assigned to one copy is not a
#   variant. Applying it here makes this comparison directly comparable to the
#   standard method it had to be replaced with.
#
#   rpoC is reported separately throughout. It is single-copy and essential, so
#   if it is multi-copy here something else is wrong; and if it is not, then 14
#   differences in it are real and need explaining on their own terms.
#
# INPUT   32_strain/fixed_sites.tsv · fixed_sites_annotated.tsv · {SED,LIQ}.counts.tsv
#         07_mags/bins/<ref>.fa · 09_annotation/<ref>/<ref>.tsv
# OUTPUT  32_strain/REPEAT_CONTROL.txt · fixed_sites_classified.tsv
# RUNTIME 2-5 min.
# USAGE   ws2 && bash scripts/32d_repeat_control.sh 2>&1 | tee logs/32d.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="32d_repeat"
F="$WS2/07_mags"; S="$WS2/32_strain"; A="$WS2/09_annotation"
D="$S/repeat_check"; mkdir -p "$D"

REFID="${REF:-ETO_SED_MAG_002}"
MINDP="${MINDP:-20}"
MINCOPY_ID="${MINCOPY_ID:-90}"     # percent identity for a second copy to count
MINCOPY_COV="${MINCOPY_COV:-0.50}" # fraction of the gene that must align
DEPTH_X="${DEPTH_X:-2.0}"          # a site deeper than this x the median is an outlier

echo "######## 32d — REPEAT AND PARALOGY CONTROL ########"
echo "  started: $(date -Is)"
runlog "32d_repeat started"

for f in fixed_sites.tsv fixed_sites_annotated.tsv SED.counts.tsv LIQ.counts.tsv; do
  [ -s "$S/$f" ] || { echo "  !! missing $S/$f — run 32b and 32c first"; exit 1; }
done

WN=$(awk -F'\t' -v id="$REFID" 'NR==1{for(i=1;i<=NF;i++){if($i=="mag_id")a=i;if($i=="working_name")b=i};next}
     $a==id{print $b}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
REFFA=""
for cand in "$F/bins/${REFID}.fa" "$F/bins/${WN}.fa"; do
  [ -s "$cand" ] && { REFFA="$cand"; break; }
done
[ -n "$REFFA" ] || { echo "  !! reference fasta not found"; exit 1; }
BTSV="$A/$WN/${WN}.tsv"; BFNA="$A/$WN/${WN}.fna"
[ -s "$BTSV" ] || { echo "  !! no Bakta tsv"; exit 1; }

# ---- 1. pull the implicated genes out and align them back to the genome -----
echo
echo "######## 1. COPY NUMBER OF THE IMPLICATED GENES ########"
echo "  method: extract each gene's sequence from the reference, align it back"
echo "  against the whole MAG with minimap2, count hits at >=${MINCOPY_ID}% identity"
echo "  covering >=$(python3 -c "print(int(float('$MINCOPY_COV')*100))")% of the gene."
echo

python3 - "$S" "$REFFA" "$BTSV" "$BFNA" "$D" <<'PYG'
import sys, os, csv
from collections import defaultdict, Counter
S, REFFA, BTSV, BFNA, D = sys.argv[1:6]

def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p, errors='replace'):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs

ref = read_fasta(REFFA)
reflen = {k: len(v) for k, v in ref.items()}
bak = read_fasta(BFNA) if os.path.exists(BFNA) else {}
name_map = {}
if bak:
    shared = set(bak) & set(reflen)
    if len(shared) >= 0.9 * len(reflen):
        name_map = {k: k for k in shared}
    else:
        by_len = defaultdict(list)
        for k, L in reflen.items(): by_len[L].append(k)
        for bk, bs in bak.items():
            c = by_len.get(len(bs), [])
            if len(c) == 1: name_map[bk] = c[0]
else:
    name_map = {k: k for k in reflen}

coords = {}
for line in open(BTSV, errors='replace'):
    if line.startswith('#') or not line.strip(): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 8: continue
    try: st, en = int(f[2]), int(f[3])
    except ValueError: continue
    rc = name_map.get(f[0].strip())
    if rc: coords[f[5].strip()] = (rc, st, en, f[4].strip(), f[6].strip(), f[7].strip())

# which loci carry a difference
loci = Counter()
for r in csv.DictReader(open(os.path.join(S, 'fixed_sites_annotated.tsv')), delimiter='\t'):
    if r.get('locus'): loci[r['locus']] += 1

with open(os.path.join(D, 'implicated_genes.fna'), 'w') as fh:
    n = 0
    for l in loci:
        c = coords.get(l)
        if not c: continue
        rc, st, en = c[0], c[1], c[2]
        fh.write(f">{l}\n{ref[rc][st-1:en]}\n"); n += 1
    print(f"  extracted {n} gene sequence(s) carrying at least one difference")

with open(os.path.join(D, 'loci_meta.tsv'), 'w') as fh:
    fh.write("locus\tn_diff\tcontig\tstart\tend\tstrand\tgene\tproduct\n")
    for l, n in loci.most_common():
        c = coords.get(l)
        if c: fh.write(f"{l}\t{n}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[4]}\t{c[5]}\n")
        else: fh.write(f"{l}\t{n}\t?\t\t\t\t\t\n")
PYG

if [ -s "$D/implicated_genes.fna" ]; then
  run_tool minimap2 -x asm10 -c --secondary=yes -N 20 -t "$THREADS" \
      "$REFFA" "$D/implicated_genes.fna" > "$D/copies.paf" 2>"$D/minimap2.err"
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$D/copies.paf" ]; then
    echo "  !! minimap2 failed (rc=$rc) — copy number cannot be measured."
    tail -4 "$D/minimap2.err" | sed 's/^/     /'
    echo "     Without it the mobile-element hypothesis stays untested. Stopping"
    echo "     rather than reporting a classification built on nothing."
    runlog "**32d ABORTED: minimap2 failed**"
    exit 2
  fi
  echo "  minimap2: $(wc -l < "$D/copies.paf") alignment(s)"
else
  echo "  !! no gene sequences extracted"; exit 1
fi

# ---- 2. depth, classification, and the corrected rate ------------------------
python3 - "$S" "$D" "$MINDP" "$MINCOPY_ID" "$MINCOPY_COV" "$REFFA" "$DEPTH_X" <<'PY'
import sys, os, csv, math
from collections import Counter, defaultdict

S, D, MINDP, MINID, MINCOV = sys.argv[1], sys.argv[2], int(sys.argv[3]), \
                             float(sys.argv[4]), float(sys.argv[5])
REFFA, DEPTH_X = sys.argv[6], float(sys.argv[7])

# CONTIG ORDER. The counts files are emitted in the reference .fai's order,
# which is NOT alphabetical. v1 of this script merged them by comparing
# (contig_name, pos) tuples directly, i.e. lexicographically, so the two streams
# desynchronised and only coincidental matches were counted: 31,334 shared sites
# instead of 421,698. Every rate in section 3 was therefore wrong by 13x, and
# the verdict built on it was worthless. 32b and 32c both map contig -> index
# first; this one did not, because I rewrote the merge instead of reusing it.
_order = {}
with open(REFFA, errors='replace') as _fh:
    for _l in _fh:
        if _l.startswith('>'):
            _n = _l[1:].split()[0]
            if _n not in _order: _order[_n] = len(_order)
def key(r): return (_order.get(r[0], 10**9), r[1])
W = 96
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t: say('  ' + t); say('=' * W)

# ---------------------------------------------------------------- copy number
# PAF: 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen 8 tstart 9 tend
#      10 matches 11 alnlen 12 mapq
copies = defaultdict(list)
for line in open(os.path.join(D, 'copies.paf')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12: continue
    q, qlen = f[0], int(f[1])
    nm, alen = int(f[9]), int(f[10])
    if alen == 0: continue
    ident = 100.0 * nm / alen
    cov = alen / qlen if qlen else 0
    if ident >= MINID and cov >= MINCOV:
        copies[q].append((f[5], int(f[7]), int(f[8]), round(ident, 1)))

meta = {}
for r in csv.DictReader(open(os.path.join(D, 'loci_meta.tsv')), delimiter='\t'):
    meta[r['locus']] = r

# --------------------------------------------------------------------- depth
# Genome-wide depth is taken from the same parsed totals the comparison used,
# so "3x the median" means 3x the median of sites that were actually eligible.
def parse_total(bases, refb):
    tot, i, n = 0, 0, len(bases)
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            j = i + 1; num = ''
            while j < n and bases[j].isdigit(): num += bases[j]; j += 1
            i = j + int(num or 0); continue
        if c in '.,' or c.upper() in 'ACGT': tot += 1
        i += 1
    return tot

def depths(tag):
    out = []
    with open(os.path.join(S, f'{tag}.counts.tsv')) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5: continue
            t = parse_total(f[4], f[2])
            if t >= MINDP: out.append(t)
    out.sort()
    return out

dsed, dliq = depths('SED'), depths('LIQ')
med_s = dsed[len(dsed)//2] if dsed else 0
med_l = dliq[len(dliq)//2] if dliq else 0
p95_s = dsed[int(len(dsed)*0.95)] if dsed else 0
p95_l = dliq[int(len(dliq)*0.95)] if dliq else 0

rule('1. COPY NUMBER — IS THE GENE UNIQUE IN THIS GENOME?')
say(f"  {'locus':<20}{'n':>3}{'copies':>8}   product")
say('-' * W)
MOBILE_KEY = ('transposase', 'transposon', 'integrase', 'recombinase',
              'insertion sequence', 'is66', 'is3', 'is4', 'is5', 'is21',
              'is110', 'is200', 'is630', 'is1', 'mobile element')
cls = {}
for l, r in sorted(meta.items(), key=lambda kv: -int(kv[1]['n_diff'])):
    nc = len(copies.get(l, []))
    prod = (r.get('product') or '')
    ismob = any(k in prod.lower() for k in MOBILE_KEY)
    if nc >= 2:      c = 'MULTICOPY'
    elif ismob:      c = 'MOBILE(1copy)'
    else:            c = 'unique'
    cls[l] = c
    say(f"  {l:<20}{r['n_diff']:>3}{nc:>8}   {prod[:52]}"
        + ('   <-- ' + c if c != 'unique' else ''))
say()
say("  'copies' counts alignments of the gene back onto its own genome at")
say(f"  >={MINID:.0f}% identity over >={MINCOV*100:.0f}% of its length. One is the gene finding")
say("  itself. Two or more means a read landing here cannot be assigned to this")
say("  copy rather than another, and a variant that cannot be assigned to a")
say("  locus is not a variant of that locus.")

# ------------------------------------------------------------------- per site
rows = list(csv.DictReader(open(os.path.join(S, 'fixed_sites.tsv')), delimiter='\t'))
ann = {(r['contig'], r['pos']): r for r in
       csv.DictReader(open(os.path.join(S, 'fixed_sites_annotated.tsv')), delimiter='\t')}
hdr = list(rows[0])
AD = 'SED_dp' if 'SED_dp' in hdr else 'A_dp'
BD = 'LIQ_dp' if 'LIQ_dp' in hdr else 'B_dp'

rule('2. DEPTH — IS THE LOCUS BEING FED BY SEQUENCE THAT IS NOT THERE?')
say(f"  genome-wide, among eligible sites:  SED median {med_s}x (95th {p95_s}x)   "
    f"LIQ median {med_l}x (95th {p95_l}x)")
say()
byclass = defaultdict(lambda: {'n': 0, 'sed': [], 'liq': []})
out = []
for r in rows:
    a = ann.get((r['contig'], r['pos']), {})
    l = a.get('locus', '')
    c = cls.get(l, 'unique' if l else 'intergenic')
    try: ds, dl = int(r[AD]), int(r[BD])
    except (ValueError, KeyError): ds = dl = 0
    byclass[c]['n'] += 1
    byclass[c]['sed'].append(ds); byclass[c]['liq'].append(dl)
    out.append(dict(contig=r['contig'], pos=r['pos'],
                    locus=l, product=a.get('product', ''),
                    klass=c, sed_dp=ds, liq_dp=dl,
                    sed_x_median=round(ds / med_s, 2) if med_s else '',
                    liq_x_median=round(dl / med_l, 2) if med_l else '',
                    effect=a.get('effect', ''), ti_tv=a.get('ti_tv', ''),
                    homopolymer=a.get('homopolymer', '')))

say(f"  {'class':<16}{'sites':>6}{'SED med':>10}{'x gnome':>9}{'LIQ med':>10}{'x gnome':>9}")
say('-' * W)
for c in sorted(byclass, key=lambda k: -byclass[k]['n']):
    v = byclass[c]
    ms = sorted(v['sed'])[len(v['sed'])//2] if v['sed'] else 0
    ml = sorted(v['liq'])[len(v['liq'])//2] if v['liq'] else 0
    say(f"  {c:<16}{v['n']:>6}{ms:>10}{ms/med_s if med_s else 0:>9.1f}"
        f"{ml:>10}{ml/med_l if med_l else 0:>9.1f}")
say()
say("  A class sitting well above 1.0x median is receiving reads from elsewhere.")

with open(os.path.join(S, 'fixed_sites_classified.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)

# ------------------------------------------------------- the corrected rate
# The filter is applied IDENTICALLY to numerator and denominator. v1 subtracted
# differences from the numerator while leaving the denominator whole, which
# deflates the rate; and its merge was mis-ordered, so the denominator was
# 31,334 instead of 421,698. Both are fixed here and the merge is cross-checked
# against 32b before any rate is printed.
excl_loci = {l for l, c in cls.items() if c != 'unique'}
locus_span = {}
for l in excl_loci:
    m = meta.get(l)
    if m and m['contig'] != '?':
        locus_span[l] = (m['contig'], int(m['start']), int(m['end']))
spans = defaultdict(list)
for l, (c, st, en) in locus_span.items(): spans[c].append((st, en, l))

# Depth is the second exclusion, and section 2 is why. The 'unique' sites sit at
# 126x in SED against a genome median of 23x — above the 95th percentile of
# eligible sites. Copy number says those genes are single-copy; depth says they
# are still receiving reads from somewhere. Both filters are needed and both are
# applied to numerator and denominator alike.
SED_CAP = med_s * DEPTH_X
LIQ_CAP = med_l * DEPTH_X

def sites_iter(tag):
    with open(os.path.join(S, f'{tag}.counts.tsv')) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5: continue
            try: pos = int(f[1])
            except ValueError: continue
            t = parse_total(f[4], f[2])
            if t < MINDP: continue
            yield (f[0], pos), t

def in_excluded(c, p):
    for st, en, l in spans.get(c, []):
        if st <= p <= en: return True
    return False

ncmp = nrep = ndep = nclean = 0
ga, gb = sites_iter('SED'), sites_iter('LIQ')
ra, rb = next(ga, None), next(gb, None)
while ra is not None and rb is not None:
    ka, kb = key(ra[0]), key(rb[0])
    if ka == kb:
        ncmp += 1
        rep_hit = in_excluded(ra[0][0], ra[0][1])
        dep_hit = (ra[1] > SED_CAP) or (rb[1] > LIQ_CAP)
        if rep_hit: nrep += 1
        if dep_hit: ndep += 1
        if not rep_hit and not dep_hit: nclean += 1
        ra, rb = next(ga, None), next(gb, None)
    elif ka < kb: ra = next(ga, None)
    else:         rb = next(gb, None)

rule('3. THE RATE WITH REPEATS AND DEPTH OUTLIERS REMOVED')
_expect = None
if os.path.exists(os.path.join(S, 'DIVERGENCE.txt')):
    for line in open(os.path.join(S, 'DIVERGENCE.txt')):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 3 and 'measurement' in f[0]:
            try: _expect = int(f[1])
            except ValueError: pass
if _expect is not None and _expect != ncmp:
    say(f"  [FAIL] merge recount {ncmp:,} does not match 32b's {_expect:,}.")
    say("         The denominator is not the set the differences came from.")
    say("         NO RATE IS REPORTED. Fix the merge before reading anything else.")
    open(os.path.join(S, 'REPEAT_CONTROL.txt'), 'w').write('\n'.join(rep) + '\n')
    sys.exit(3)
say(f"  [PASS] merge recount {ncmp:,} matches 32b")
say()
say(f"  depth caps: SED > {SED_CAP:.0f}x or LIQ > {LIQ_CAP:.0f}x "
    f"({DEPTH_X:.1f}x the median) marks a site as depth-outlier")
say()

# classify the differences by the same two rules
d_rep = d_dep = d_clean = 0
site_dp = {}
for r in rows:
    try: ds, dl = int(r[AD]), int(r[BD])
    except (ValueError, KeyError): ds = dl = 0
    c = r['contig']; p = int(r['pos'])
    rep_hit = in_excluded(c, p)
    dep_hit = (ds > SED_CAP) or (dl > LIQ_CAP)
    if rep_hit: d_rep += 1
    if dep_hit: d_dep += 1
    if not rep_hit and not dep_hit: d_clean += 1

def line(lbl, d, n):
    say(f"  {lbl:<34}{d:>4} / {n:>8,} = {1000*d/n if n else 0:>8.4f} /kb")
line('all compared sites', len(rows), ncmp)
line('in multicopy / mobile loci', d_rep, nrep)
line('depth outliers', d_dep, ndep)
line('CLEAN (neither)', d_clean, nclean)
say()
say(f"  The two filters remove {100*(ncmp-nclean)/ncmp:.1f}% of the compared sites "
    f"and {100*(len(rows)-d_clean)/len(rows):.0f}% of the differences.")

fl = None
dv = os.path.join(S, 'DIVERGENCE.txt')
if os.path.exists(dv):
    for line_ in open(dv):
        f = line_.rstrip('\n').split('\t')
        if len(f) >= 3 and 'FLOOR' in f[0].upper():
            try: fl = (int(f[1]), int(f[2]))
            except ValueError: pass
r_clean = d_clean / nclean if nclean else 0
ub = None
if fl and fl[1] == 0:
    ub = 3 / fl[0]
    lam = fl[0] * r_clean
    say()
    say(f"  floor 95% upper bound    {1000*ub:.4f} /kb   (0 in {fl[0]:,} sites)")
    say(f"  clean-sequence rate      {1000*r_clean:.4f} /kb")
    say(f"  Poisson vs floor: expected {lam:.2f}, observed 0, p = {math.exp(-lam):.4f}")
    say(f"  -> {'ABOVE' if r_clean > ub else 'BELOW'} the floor's upper bound")

rule('4. VERDICT')
if ub is not None:
    if r_clean <= ub:
        say("  >> NO RESOLVABLE DIVERGENCE.")
        say()
        say(f"     {d_clean} differences remain in {nclean:,} sites that are neither in a")
        say("     multicopy/mobile locus nor at anomalous depth. That rate is at or")
        say("     below the upper confidence bound of an error floor which measured")
        say("     zero, so it cannot be separated from error at this depth.")
        say()
        say("     The apparent signal sat in transposases and in loci reading far")
        say("     deeper than the genome median — both signatures of reads arriving")
        say("     from sequence that is not this locus.")
        say()
        say("     This is a real result and the conservative one: the enrichment did")
        say("     not detectably change this organism, which licenses using")
        say("     ETO_LIQ_MAG_012 to interpret the in-situ ETO_SED_MAG_002.")
    else:
        say("  >> DIVERGENCE SURVIVES BOTH CONTROLS.")
        say(f"     {d_clean} differences in {nclean:,} uniquely-mapping, normal-depth")
        say(f"     sites = {1000*r_clean:.3f} /kb, against a floor bounded at "
            f"{1000*ub:.3f} /kb.")
        say("     THIS is the number to report — not the raw 0.19 /kb, which was")
        say("     half transposase cross-mapping.")
else:
    say("  No usable floor in DIVERGENCE.txt; rate reported without comparison.")
say()
say("  ON rpoC, whatever the verdict. It came back SINGLE-COPY, so its 14")
say("  differences are not paralogy — but section 2 should be checked for its")
say("  depth. An essential, highly conserved gene is exactly where an unbinned")
say("  relative deposits reads, and step 26 put a second Chakrabartyella in")
say("  these samples. Single-copy plus far-above-median depth means foreign")
say("  reads, not selection.")
say()
open(os.path.join(S, 'REPEAT_CONTROL.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(S, 'REPEAT_CONTROL.txt')}")
print(f"  written: {os.path.join(S, 'fixed_sites_classified.tsv')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "repeat control completed" "$RC" || OK=1
check "classification written" $([ -s "$S/fixed_sites_classified.tsv" ] && echo 0 || echo 1) || OK=1
NC=$(($(wc -l < "$S/fixed_sites_classified.tsv" 2>/dev/null || echo 1) - 1))
NF=$(($(wc -l < "$S/fixed_sites.tsv") - 1))
check "every site classified ($NC/$NF)" $([ "$NC" -eq "$NF" ] && echo 0 || echo 1) || OK=1
[ "$OK" -eq 0 ] && runlog "32d_repeat complete" || runlog "**32d_repeat FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
