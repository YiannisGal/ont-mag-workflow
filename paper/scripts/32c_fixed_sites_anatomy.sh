#!/usr/bin/env bash
# =============================================================================
# 32c_fixed_sites_anatomy.sh  --  what ARE the 81 fixed differences?
#
# WHY THIS IS NEEDED, AND WHY THE VERDICT IN 32b IS NOT ENOUGH
#   32b returned:
#       sedA vs sedB  (ERROR FLOOR)    28,163 sites      0 diffs   0.0000 /kb
#       SED  vs LIQ   (measurement)   421,698 sites     81 diffs   0.1921 /kb
#
#   Read carelessly that is "81 versus zero" and the answer is obvious. It is
#   not obvious, because the two comparisons DO NOT HAVE THE SAME POWER. The
#   split halves each carry half the depth, so far fewer of their sites reach
#   20x in BOTH members of the pair: 28,163 against 421,698, a 15-fold
#   difference. AN OBSERVED ZERO IS NOT A RATE OF ZERO.
#
#   Two defensible tests disagree:
#
#   1. Poisson. If the floor's true rate equalled the measurement's, we would
#      expect 5.41 differences in 28,163 sites. Observing 0 has p = 0.0045.
#      -> the measurement is above the floor.
#
#   2. This project's own 3x criterion, which 32b applies whenever the floor is
#      non-zero. The rule of three puts the 95% upper bound on a zero-count
#      floor at 3/28,163 = 0.1065 /kb. Three times that is 0.3196 /kb. The
#      measurement is 0.1921 /kb. -> NOT DISTINGUISHABLE.
#
#   The zero-floor branch in 32b sidestepped criterion 2 by not having a number
#   to multiply. That was an accident of coding, not a decision, and it must not
#   become the reason a borderline result reads as a clean one.
#
#   So the rate alone cannot settle this. What can is WHERE the 81 sites are.
#   Sequencing error is positionally indifferent. Selection is not. If the sites
#   concentrate in a few genes, that is a signal no error process produces, and
#   it survives the power problem entirely. If they are spread evenly and sit in
#   homopolymers, the conservative reading wins.
#
# WHAT THIS MEASURES
#   1. Per-contig and per-gene distribution, against the RIGHT denominator —
#      the number of sites actually compared in each gene, not gene length.
#      A gene covered at 5x contributes no compared sites and cannot host a
#      difference; scoring it against its length would invent enrichment.
#   2. Homopolymer context. ONT's characteristic error is the homopolymer, and
#      although this analysis counts only substitutions (indels are stripped by
#      the pileup parser), a substitution called at the edge of a 5-mer run is
#      still the most suspect class of call in the set.
#   3. Transition/transversion ratio. Real divergence runs Ti-heavy (~2:1 or
#      more). A ratio near 0.5 is what random error looks like.
#   4. Synonymous vs non-synonymous, where the site falls in a CDS.
#
# THE CONTIG-NAME HAZARD
#   Bakta renames contigs unless told not to, and step 10 did not tell it not
#   to. The fixed sites are named in the MAG fasta's coordinates. If the two
#   naming schemes do not match, this script maps them BY LENGTH and refuses
#   to proceed on any ambiguity. Joining on names that silently do not
#   correspond is how the all_mags_contigs.fna problem happened.
#
# INPUT   32_strain/fixed_sites.tsv · {SED,LIQ}.counts.tsv
#         07_mags/bins/<ref>.fa     · 09_annotation/<ref>/<ref>.tsv
# OUTPUT  32_strain/fixed_sites_annotated.tsv · ANATOMY.txt
# RUNTIME 1-2 min.
# USAGE   ws2 && bash scripts/32c_fixed_sites_anatomy.sh 2>&1 | tee logs/32c.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="32c_anatomy"
F="$WS2/07_mags"; S="$WS2/32_strain"; A="$WS2/09_annotation"

REFID="${REF:-ETO_SED_MAG_002}"
MINDP="${MINDP:-20}"

echo "######## 32c — ANATOMY OF THE FIXED DIFFERENCES ########"
echo "  started: $(date -Is)"
runlog "32c_anatomy started"

[ -s "$S/fixed_sites.tsv" ] || { echo "  !! no fixed_sites.tsv — run 32b first"; exit 1; }

WN=$(awk -F'\t' -v id="$REFID" 'NR==1{for(i=1;i<=NF;i++){if($i=="mag_id")a=i;if($i=="working_name")b=i};next}
     $a==id{print $b}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
REFFA=""
for cand in "$F/bins/${REFID}.fa" "$F/bins/${WN}.fa"; do
  [ -s "$cand" ] && { REFFA="$cand"; break; }
done
[ -n "$REFFA" ] || { echo "  !! reference fasta not found"; exit 1; }
BAK="$A/$WN"
[ -s "$BAK/${WN}.tsv" ] || { echo "  !! no Bakta annotation at $BAK"; exit 1; }
echo "  reference: ${REFFA#$WS2/}"
echo "  annotation: ${BAK#$WS2/}/${WN}.tsv"

python3 - "$S" "$REFFA" "$BAK/${WN}.tsv" "$BAK/${WN}.fna" "$MINDP" <<'PY'
import sys, os, csv, re
from collections import defaultdict, Counter

S, REFFA, BTSV, BFNA, MINDP = sys.argv[1:5] + [int(sys.argv[5])]
W = 96
def rule(t=''):
    print('=' * W)
    if t: print('  ' + t); print('=' * W)

# ------------------------------------------------------------ reference seqs
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
print(f"  reference: {len(ref)} contigs, {sum(len(v) for v in ref.values())} bp")

# ---------------------------------------------- Bakta contig name reconciliation
bakta_names = {}
if os.path.exists(BFNA):
    bs = read_fasta(BFNA)
    bakta_names = {k: len(v) for k, v in bs.items()}

reflen = {k: len(v) for k, v in ref.items()}
name_map = {}                      # bakta contig -> reference contig
if bakta_names:
    shared = set(bakta_names) & set(reflen)
    if len(shared) >= 0.9 * len(reflen):
        name_map = {k: k for k in shared}
        print(f"  contig names match directly ({len(shared)}/{len(reflen)})")
    else:
        by_len = defaultdict(list)
        for k, L in reflen.items(): by_len[L].append(k)
        amb = 0
        for bk, L in bakta_names.items():
            cands = by_len.get(L, [])
            if len(cands) == 1: name_map[bk] = cands[0]
            else: amb += 1
        print(f"  contig names differ — mapped {len(name_map)}/{len(bakta_names)} by length"
              f"{f', {amb} ambiguous' if amb else ''}")
        if amb > 0.1 * len(bakta_names):
            print("  !! too many contigs share a length to map safely. Gene-level")
            print("     results would be joins onto the wrong sequence. Re-run Bakta")
            print("     with --keep-contig-headers, or accept the positional results")
            print("     below and ignore every gene column.")
else:
    name_map = {k: k for k in reflen}
    print("  no Bakta .fna found — assuming contig names already agree")

# ---------------------------------------------------------------- annotation
genes = defaultdict(list)          # ref contig -> [(start, stop, strand, locus, gene, product, type)]
nfeat = 0
for line in open(BTSV, errors='replace'):
    if line.startswith('#') or not line.strip(): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 8: continue
    try: st, en = int(f[2]), int(f[3])
    except ValueError: continue
    rc = name_map.get(f[0].strip())
    if rc is None: continue
    genes[rc].append((st, en, f[4].strip(), f[5].strip(), f[6].strip(), f[7].strip(), f[1].strip()))
    nfeat += 1
for k in genes: genes[k].sort()
print(f"  annotation: {nfeat} features on {len(genes)} contigs")

def find_gene(c, p):
    for st, en, strand, locus, gname, prod, ftype in genes.get(c, []):
        if st <= p <= en: return (st, en, strand, locus, gname, prod, ftype)
        if st > p: break
    return None

# ------------------------------------- the denominator: sites actually compared
# Per contig AND per gene. Recomputed by streaming the two counts files in
# lockstep — the same merge 32b used, so the denominator is the identical set of
# sites the 81 came out of.
# This MUST be 32b's parser, applied with 32b's threshold. The counts files were
# pre-filtered in awk on the pileup's own depth column, but that column counts
# read-start markers and indel runs which the parser then strips — so the
# post-parse total is often lower and a site that passed the awk filter can fail
# the real one. Yielding every line inflated the denominator by 5,353 sites
# (427,051 against 32b's 421,698). A denominator that is not the same set the
# numerator came from makes every per-gene rate below slightly wrong.
def parse(bases, refb):
    out, i, n = Counter(), 0, len(bases)
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            j = i + 1; num = ''
            while j < n and bases[j].isdigit(): num += bases[j]; j += 1
            i = j + int(num or 0); continue
        if c in '.,': out[refb.upper()] += 1
        elif c.upper() in 'ACGT': out[c.upper()] += 1
        i += 1
    return out

def sites(tag):
    p = os.path.join(S, f'{tag}.counts.tsv')
    with open(p) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5: continue
            try: pos = int(f[1])
            except ValueError: continue
            if sum(parse(f[4], f[2]).values()) < MINDP: continue
            yield f[0], pos

order = {c: i for i, c in enumerate(ref)}
def key(r): return (order.get(r[0], 10**9), r[1])

cmp_contig = Counter(); cmp_gene = Counter(); ncmp = 0
ga, gb = sites('SED'), sites('LIQ')
ra, rb = next(ga, None), next(gb, None)
while ra is not None and rb is not None:
    ka, kb = key(ra), key(rb)
    if ka == kb:
        ncmp += 1
        cmp_contig[ra[0]] += 1
        g = find_gene(ra[0], ra[1])
        if g: cmp_gene[g[3]] += 1
        ra, rb = next(ga, None), next(gb, None)
    elif ka < kb: ra = next(ga, None)
    else:         rb = next(gb, None)
print(f"  sites compared (recomputed): {ncmp}")

# Cross-check against 32b's own count. If these disagree, the per-gene rates
# below are computed over a different set of sites than the differences were,
# and every one of them is wrong by that margin.
_dv = os.path.join(S, 'DIVERGENCE.txt')
_expect = None
if os.path.exists(_dv):
    for line in open(_dv):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 3 and 'measurement' in f[0]:
            try: _expect = int(f[1])
            except ValueError: pass
if _expect is not None:
    if _expect == ncmp:
        print(f"  [PASS] matches 32b's {_expect} compared sites")
    else:
        print(f"  [WARN] 32b reported {_expect} compared sites, this recount says "
              f"{ncmp} (difference {abs(_expect - ncmp)})")
        print("         The two are meant to be the identical set. Per-gene rates")
        print("         below are only as good as this agreement.")

# -------------------------------------------------------------- the 81 sites
# 32b's dump names its columns A_/B_ because compare() is generic over any pair;
# the earlier version of 32 wrote SED_/LIQ_. Read whichever is present rather
# than assuming, and say which mapping was used.
_raw = list(csv.DictReader(open(os.path.join(S, 'fixed_sites.tsv')), delimiter='\t'))
if not _raw:
    print("  !! fixed_sites.tsv is empty"); sys.exit(1)
_hdr = list(_raw[0])
if 'SED_allele' in _hdr:   A_COL, B_COL = 'SED_allele', 'LIQ_allele'
elif 'A_allele' in _hdr:   A_COL, B_COL = 'A_allele', 'B_allele'
else:
    print(f"  !! cannot find allele columns in fixed_sites.tsv: {_hdr}")
    sys.exit(1)
rows = _raw
print(f"  fixed differences: {len(rows)}   (allele columns: {A_COL} / {B_COL},")
print(f"                                    which for this dump are SED / LIQ)")
print()

TS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
def homopolymer(c, p, w=2, minrun=4):
    """Longest homopolymer run touching a window around the site."""
    s = ref.get(c, '')
    if not s: return 0
    lo, hi = max(0, p - 1 - w - 8), min(len(s), p + w + 8)
    seg = s[lo:hi].upper()
    best, run, prev = 0, 0, None
    idx = p - 1 - lo
    for i, ch in enumerate(seg):
        run = run + 1 if ch == prev else 1
        prev = ch
        if abs(i - idx) <= w and run >= minrun: best = max(best, run)
        # also catch a run that ends after the site
    # second pass: run lengths spanning the site
    j = idx
    if 0 <= j < len(seg):
        ch = seg[j]; L = 1
        k = j - 1
        while k >= 0 and seg[k] == ch: L += 1; k -= 1
        k = j + 1
        while k < len(seg) and seg[k] == ch: L += 1; k += 1
        best = max(best, L if L >= minrun else 0)
    return best

CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[i]; i += 1
COMP = str.maketrans('ACGT', 'TGCA')

out = []
tv = Counter(); hp = 0; syn = nonsyn = noncoding = nocds = 0
per_contig = Counter(); per_gene = Counter(); gene_info = {}

for r in rows:
    c, p = r['contig'], int(r['pos'])
    a, b = r[A_COL], r[B_COL]
    tv['Ti' if (a, b) in TS else 'Tv'] += 1
    h = homopolymer(c, p)
    if h: hp += 1
    per_contig[c] += 1
    g = find_gene(c, p)
    eff = ''
    if g:
        st, en, strand, locus, gname, prod, ftype = g
        per_gene[locus] += 1; gene_info[locus] = (gname, prod, ftype, en - st + 1)
        if ftype.lower() == 'cds':
            seq = ref[c][st - 1:en]
            off = p - st
            if strand == '-':
                seq = seq.translate(COMP)[::-1]
                off = (en - st) - off
                a2 = a.translate(COMP); b2 = b.translate(COMP)
            else:
                a2, b2 = a, b
            ci = off // 3; cp = off % 3
            cod = seq[ci * 3:ci * 3 + 3].upper()
            if len(cod) == 3 and cod in CODON:
                mut = cod[:cp] + b2 + cod[cp + 1:]
                if mut in CODON:
                    eff = 'syn' if CODON[cod] == CODON[mut] else \
                          f'nonsyn {CODON[cod]}->{CODON[mut]}'
                    if CODON[cod] == CODON[mut]: syn += 1
                    else: nonsyn += 1
                else: nocds += 1
            else: nocds += 1
        else:
            noncoding += 1
    else:
        noncoding += 1
    out.append(dict(contig=c, pos=p, SED=a, LIQ=b,
                    ti_tv='Ti' if (a, b) in TS else 'Tv',
                    homopolymer=h,
                    locus=g[3] if g else '', gene=g[4] if g else '',
                    product=g[5] if g else 'intergenic',
                    effect=eff))

with open(os.path.join(S, 'fixed_sites_annotated.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)

# ------------------------------------------------------------------- report
rep = []
def say(s=''):
    print(s); rep.append(s)

rule('1. IS IT ERROR-SHAPED OR BIOLOGY-SHAPED?')
say(f"  transitions {tv['Ti']}   transversions {tv['Tv']}   "
    f"Ti/Tv = {tv['Ti']/tv['Tv']:.2f}" if tv['Tv'] else
    f"  transitions {tv['Ti']}   transversions 0")
say("    Real divergence is transition-biased, typically 2:1 or higher, because")
say("    transitions are chemically easier and more often synonymous. Random")
say("    sequencing error has no such bias and lands near 0.5 (there are twice")
say("    as many possible transversions as transitions).")
say()
say(f"  sites in or beside a homopolymer run of >=4: {hp}/{len(rows)} "
    f"({100*hp/len(rows):.0f}%)")
say("    This analysis counts substitutions only — the pileup parser strips")
say("    indels — so homopolymers are not the usual ONT indel problem here.")
say("    They remain the least trustworthy substitution calls in the set.")

rule('2. DISTRIBUTION — THE TEST THAT SURVIVES THE POWER PROBLEM')
say("  Sequencing error is positionally indifferent. Selection is not. Genes are")
say("  scored against the number of sites ACTUALLY COMPARED in them, not their")
say("  length: a gene with no high-depth sites cannot host a difference, and")
say("  scoring it by length would invent enrichment where there was no coverage.")
say()
say(f"  {len(rows)} differences over {ncmp} compared sites "
    f"= {1000*len(rows)/ncmp:.4f} /kb")
say(f"  spread over {len(per_contig)} contig(s) and {len(per_gene)} gene(s)")
say()
hot = [(l, n) for l, n in per_gene.most_common() if n >= 3]
if hot:
    say("  Genes carrying 3 or more differences:")
    say(f"    {'locus':<16}{'n':>3}{'compared':>10}{'per kb':>9}  gene / product")
    for l, n in hot:
        cs = cmp_gene.get(l, 0)
        gname, prod, ftype, glen = gene_info[l]
        rate = 1000 * n / cs if cs else 0
        say(f"    {l:<16}{n:>3}{cs:>10}{rate:>9.2f}  "
            f"{(gname or '-')[:12]:<13}{prod[:38]}")
    say()
    say(f"  A gene whose per-kb rate is many times the genome-wide {1000*len(rows)/ncmp:.2f} /kb")
    say("  is the claim worth making. It does not depend on the error floor at all.")
else:
    say("  NO GENE CARRIES 3 OR MORE DIFFERENCES.")
    say()
    say("  That is the result, and it points the conservative way. The 81 sites")
    say("  are scattered, which is what a low-level uniform process looks like —")
    say("  residual error, or genuine but undirected drift. There is no gene here")
    say("  to name as a target of selection.")

rule('3. CODING EFFECT')
tot_cds = syn + nonsyn
if tot_cds:
    say(f"  synonymous {syn}   non-synonymous {nonsyn}   "
        f"dN/dS-like ratio {nonsyn/syn:.2f}" if syn else
        f"  synonymous 0   non-synonymous {nonsyn}")
    say(f"  non-coding or unresolvable: {noncoding + nocds}")
    say()
    say("  Under neutrality roughly 3 in 4 random substitutions in a CDS are")
    say("  non-synonymous, so a ratio near 3 means nothing is acting on these")
    say("  sites. Well below 3 is purifying selection; well above, with the")
    say("  clustering from section 2, would be positive selection.")
else:
    say("  No site resolved to a codon. Check the contig-name mapping above.")

rule('4. WHAT TO REPORT')
# Read the two comparisons back out of 32b's own output rather than writing
# their numbers into this script. A figure typed into a report by hand is a
# figure that goes stale the first time a threshold changes.
import math
fl = me = None
dvp = os.path.join(S, 'DIVERGENCE.txt')
if os.path.exists(dvp):
    for line in open(dvp):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4: continue
        try: n, k = int(f[1]), int(f[2])
        except ValueError: continue
        if 'FLOOR' in f[0].upper():        fl = (f[0].strip(), n, k)
        elif 'measurement' in f[0]:        me = (f[0].strip(), n, k)

if fl and me:
    r_me = me[2] / me[1] if me[1] else 0
    say("  The rate comparison is BORDERLINE. It must not be written as "
        f"'{me[2]} vs {fl[2]}'.")
    say(f"    measurement            {1000*r_me:.4f} /kb  ({me[2]} / {me[1]:,})")
    say(f"    floor, observed        {1000*fl[2]/fl[1] if fl[1] else 0:.4f} /kb  "
        f"({fl[2]} / {fl[1]:,})")
    if fl[2] == 0 and fl[1]:
        ub = 3 / fl[1]
        lam = fl[1] * r_me
        p = math.exp(-lam)
        say(f"    floor, 95% upper bound {1000*ub:.4f} /kb  "
            f"(rule of three: 3/{fl[1]:,})")
        say(f"    Poisson test of equal rates:  expected {lam:.2f}, observed 0, "
            f"p = {p:.4f}")
        verdict3 = "clears" if r_me > 3 * ub else "does NOT clear"
        say(f"    this project's 3x criterion:  {1000*r_me:.4f} vs "
            f"3 x {1000*ub:.4f} = {3000*ub:.4f}  -> {verdict3}")
        say()
        say("  AN OBSERVED ZERO IS NOT A RATE OF ZERO. The two tests disagree")
        say(f"  because the split-half control has {me[1]/fl[1]:.0f}x less power:")
        say("  halving the depth halves the sites reaching the threshold in BOTH")
        say("  members of the pair.")
    else:
        say("  (floor is non-zero — 32b's own 3x rule applies directly)")
    say()
else:
    say("  DIVERGENCE.txt not readable — rate comparison not restated here.")
    say()
if hot:
    say("  BUT section 2 found clustering, and clustering is not a rate argument.")
    say("  Lead with the genes; report the rate with both tests stated.")
else:
    say("  Section 2 found no clustering either. The defensible sentence is:")
    say()
    say("    'Fixed differences between the sediment and enrichment populations")
    say("     of ETO_SED_MAG_002 occur at 0.19 per kb across 421,698 sites")
    say("     compared at >=20x. A split-half control on the sediment library")
    say("     alone yielded no fixed differences in 28,163 sites, bounding the")
    say("     error floor below 0.11 per kb. The measurement exceeds that bound")
    say("     but by under a factor of two, and no gene carries more than two")
    say("     differences. We therefore report the enrichment as having produced")
    say("     no resolvable divergence in this organism.'")
    say()
    say("  That is the conservative reading and it is still a real result: it")
    say("  licenses using ETO_LIQ_MAG_012 to interpret the in-situ population.")
say()
say("  To strengthen the floor rather than argue about it, raise its power:")
say("  MINDP=10 on the split halves alone would bring their site count near the")
say("  measurement's. That changes the error regime, so it is a second control,")
say("  not a replacement — but two controls bracketing the answer beats one.")

open(os.path.join(S, 'ANATOMY.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(S, 'fixed_sites_annotated.tsv')}")
print(f"  written: {os.path.join(S, 'ANATOMY.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "anatomy completed" "$RC" || OK=1
check "annotated sites written" $([ -s "$S/fixed_sites_annotated.tsv" ] && echo 0 || echo 1) || OK=1
NANN=$(($(wc -l < "$S/fixed_sites_annotated.tsv" 2>/dev/null || echo 1) - 1))
NFIX=$(($(wc -l < "$S/fixed_sites.tsv") - 1))
check "every fixed site annotated ($NANN/$NFIX)" \
      $([ "$NANN" -eq "$NFIX" ] && echo 0 || echo 1) || OK=1
[ "$OK" -eq 0 ] && runlog "32c_anatomy complete" || runlog "**32c_anatomy FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
