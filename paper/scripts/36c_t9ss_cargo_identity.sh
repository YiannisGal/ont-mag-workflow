#!/usr/bin/env bash
# =============================================================================
# 36c_t9ss_cargo_identity.sh  --  WHAT are the 254 T9SS cargo proteins?
#
# WHERE THIS COMES FROM
#   Step 36 gap B put a floor under T9SS cargo: 254 proteins carrying TIGR04183
#   or TIGR04131, all 254 in the three Bacteroidota genomes that carry T9SS core
#   genes, none in the other 26. That is a COUNT. It says nothing about what the
#   cargo does, and the cargo identity is the part that bears on the chapter —
#   a secretion system exporting glycoside hydrolases into anoxic sediment is a
#   polysaccharide-degradation argument; one exporting adhesins is not.
#
#   Step 36's attempt to cross the CTD set with step 23d's CAZyme call is dead:
#   cazyme_endo_exo.tsv is a PER-GENOME summary (mag_id, genus, phylum,
#   completeness, n_endo, n_exo, exo_per_endo, role). There are no protein ids in
#   it to intersect, which is why the intersection was 0 and why v4 reports that
#   as NOT TESTED rather than as a finding.
#
#   So this step gets the cargo identity from the annotation the proteins
#   themselves carry, which is guaranteed to be in the same id space as the
#   hmmsearch targets because it is the same FASTA.
#
# WHAT IT DOES
#   1. Reads the CTD protein ids out of 36_gapclose/ctd/*.tbl.
#   2. Reads every protein's description from 09_annotation/<bin>/<bin>.faa —
#      the identical file hmmsearch ran on, so the join cannot fail. If a Bakta
#      .tsv is present it is used as well, for gene names and dbxrefs.
#   3. Classifies products into functional classes by keyword.
#   4. Compares the CTD set against the genome background OF THE SAME THREE
#      GENOMES.
#   5. Optionally crosses with any per-protein CAZyme table it can find.
#
# THE CONTROL, AND IT CAN FAIL
#   A secretion system's cargo must not look like a random draw from its own
#   genome. If the CTD set's functional composition is indistinguishable from the
#   background of the same three genomes, then either the domain is not marking
#   cargo or the annotation is too coarse to see it — and in both cases no
#   statement about WHAT is exported is supportable.
#
#   The background is the SAME genomes, not the catalogue: Bacteroidota are
#   CAZyme-rich and Clostridia are not, so a catalogue-wide background would
#   manufacture enrichment out of phylum composition. That is the step-19b error
#   and it is designed out here.
#
# 🔴 TWO DEFECTS IN v1, BOTH FOUND IN ITS OWN OUTPUT AND BOTH MINE
#   1. v1 asked ONE pooled question — are 'degradative enzymes' over-represented?
#      It merged peptidases (enriched 1.9x) with lipases and sulfatases (absent
#      from the cargo entirely) and they cancelled to fold 1.00, p = 0.54. It
#      printed NO DIFFERENCE for a set whose composition is plainly different.
#      A pooled binary can cancel; a composition test cannot. v2 runs a G-test
#      over the whole class distribution, plus two-sided Fisher per class with
#      Benjamini-Hochberg. Two-sided because v1's one-sided test could not see a
#      DEPLETION, and depletion to zero (SusC/SusD 0/254, lipase 0/254) is the
#      strongest signal in this cargo.
#   2. v1 filed "T9SS type A sorting domain-containing protein" — Bakta naming a
#      protein after the very domain being searched — under 'other, annotated',
#      i.e. counted a circular name as annotated evidence AGAINST the hypothesis.
#      v2 gives those their own class, reports it descriptively, and excludes it
#      from every statistic: that class is in the CTD set by construction and
#      testing it proves only the tautology.
#
# WHAT IT WILL NOT DO
#   It will not call a hypothetical protein cargo of known function. Bacteroidota
#   polysaccharide utilisation loci are full of uncharacterised proteins, and a
#   high 'hypothetical' fraction is a real property of the cargo, not a defect to
#   be smoothed over. It is reported as its own class.
#
# OUTPUT  36_gapclose/ctd_identity.tsv · CARGO.txt
# RUNTIME < 1 min
# USAGE   ws2 && bash scripts/36c_t9ss_cargo_identity.sh 2>&1 | tee logs/36c.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="36c_cargo"
A="$WS2/09_annotation"; F="$WS2/07_mags"; T9="$WS2/23_extras/t9ss"
D="$WS2/36_gapclose"

echo "######## 36c — WHAT IS THE T9SS EXPORTING? ########"
echo "  started: $(date -Is)"

if ! ls -1 "$D/ctd"/*.tbl >/dev/null 2>&1; then
  echo "  !! no CTD tables in $D/ctd — run step 36 (gap B) first."
  exit 2
fi
echo "  CTD tables : $(ls -1 "$D/ctd"/*.tbl | wc -l)"
echo "  proteomes  : $(ls -1 "$A" 2>/dev/null | wc -l)"

echo
rm -f "$D/CARGO.txt" "$D/ctd_identity.tsv"
python3 - "$D" "$A" "$F" "$T9" "$WS2" <<'PY'
import sys, os, re, csv, glob, math
from collections import defaultdict, Counter

# 🔴 THE TOKEN FILTER, TIGHTENED -- see Rule 46.
# A filter accepting any string of >=4 characters made 'Bacteria',
# 'Methanosarcina' and '99.85' into alias keys; a global alias merge then linked
# unrelated genomes through them, and reading the resulting set gave a DIFFERENT
# answer per process because set order follows hash randomisation. Two identical
# runs of 16b returned different genera for the same genome.
# A genome id in this project always contains an underscore AND a digit
# (ETO_LIQ_MAG_002, LIQ_das_metabat_075); a genus never does in a way that
# collides (Lentimicrobium, GCA-2746605, Methanobacterium_C, JBAIYS01).
IDTOK = re.compile(r'^(?=.*_)(?=.*\d)[A-Za-z0-9._-]{6,}$')

D, A, F, T9, WS2 = sys.argv[1:6]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

STATUS = ['INCOMPLETE']

def dump_and_exit(code=0, why=None):
    """A run that stops early must still leave its reason on disk. Step 33's
    29 genomes exited 0 having written nothing, and the absence read as a
    result; an early exit that writes no report repeats that.
    The status line is written HERE and nowhere else, so a crash cannot leave
    a stale OK for the validation gate."""
    if why:
        STATUS[0] = why
    with open(os.path.join(D, 'CARGO.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'CARGO.txt')}")
    sys.exit(code)

# ------------------------------------------------------------- fisher
def _lf(n):
    return math.lgamma(n + 1)

def fisher_greater(a, b, c, d):
    """one-sided P(X >= a) for the 2x2 [[a,b],[c,d]]"""
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, r2 = a + b, c + d
    c1 = a + c
    lo, hi = max(0, c1 - r2), min(r1, c1)
    def lp(x):
        return (_lf(r1) + _lf(r2) + _lf(c1) + _lf(n - c1)
                - _lf(n) - _lf(x) - _lf(r1 - x) - _lf(c1 - x) - _lf(r2 - c1 + x))
    tot = sum(math.exp(lp(x)) for x in range(lo, hi + 1))
    num = sum(math.exp(lp(x)) for x in range(a, hi + 1))
    return min(1.0, num / tot) if tot else 1.0

def fisher_two_sided(a, b, c, d):
    """sum of every table at most as probable as the observed one.
    v1 used a ONE-SIDED test for every class, which cannot see a depletion —
    and three of this cargo's classes are depleted to zero."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, r2, c1 = a + b, c + d, a + c
    lo, hi = max(0, c1 - r2), min(r1, c1)
    def lp(x):
        return (_lf(r1) + _lf(r2) + _lf(c1) + _lf(n - c1)
                - _lf(n) - _lf(x) - _lf(r1 - x) - _lf(c1 - x) - _lf(r2 - c1 + x))
    ps = {x: math.exp(lp(x)) for x in range(lo, hi + 1)}
    tot = sum(ps.values())
    obs = ps.get(a, 0.0)
    if not tot:
        return 1.0
    return min(1.0, sum(v for v in ps.values() if v <= obs * (1 + 1e-9)) / tot)

# ---- chi-square upper tail, for the G-test (no scipy on this machine) -------
def _gser(a, x):
    ap, s, dl = a, 1.0 / a, 1.0 / a
    for _ in range(500):
        ap += 1; dl *= x / ap; s += dl
        if abs(dl) < abs(s) * 1e-14:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a))

def _gcf(a, x):
    tiny = 1e-300
    b, c, d = x + 1 - a, 1 / tiny, 1 / (x + 1 - a)
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < tiny: d = tiny
        c = b + an / c
        if abs(c) < tiny: c = tiny
        d = 1 / d
        dl = d * c
        h *= dl
        if abs(dl - 1) < 1e-14:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h

def chi2_sf(x, k):
    """P(chi2_k > x)"""
    if x <= 0 or k <= 0:
        return 1.0
    a = k / 2.0
    return 1.0 - _gser(a, x / 2.0) if (x / 2.0) < a + 1 else _gcf(a, x / 2.0)

def bh(pvals):
    """Benjamini-Hochberg q-values, order preserved"""
    idx = sorted(range(len(pvals)), key=lambda i: pvals[i])
    m = len(pvals)
    q = [1.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(idx), start=1):
        r = m - rank + 1
        val = min(prev, pvals[i] * m / r)
        q[i] = val
        prev = val
    return q

# ------------------------------------------------------- CTD protein ids
ctd = defaultdict(set)
for fn in sorted(os.listdir(os.path.join(D, 'ctd'))):
    if not fn.endswith('.tbl'):
        continue
    mag = fn[:-4]
    for line in open(os.path.join(D, 'ctd', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) >= 3:
            ctd[mag].add(f[0])
ctd = {m: s for m, s in ctd.items() if s}
NCTD = sum(len(s) for s in ctd.values())
say(f"  CTD-bearing genomes : {len(ctd)}")
say(f"  CTD proteins        : {NCTD:,}")
if not NCTD:
    say("  !! no CTD hits. Nothing to identify.")
    dump_and_exit(0, 'NO_CTD_HITS')

# ------------------------------------------------------ product strings
# The .faa hmmsearch ran on is the only source guaranteed to be in the same id
# space as its own output. Bakta's .tsv is read too when present, for gene names.
prod, gene = {}, {}
missing = []
for mag in ctd:
    faa = os.path.join(A, mag, mag + '.faa')
    if not os.path.exists(faa):
        cand = glob.glob(os.path.join(A, mag, '*.faa'))
        faa = cand[0] if cand else None
    if not faa:
        missing.append(mag); continue
    for line in open(faa, errors='replace'):
        if line.startswith('>'):
            h = line[1:].rstrip('\n')
            pid, _, desc = h.partition(' ')
            prod[pid] = desc.strip()
    tsv = os.path.join(A, mag, mag + '.tsv')
    if os.path.exists(tsv):
        hdr = None
        for line in open(tsv, errors='replace'):
            f = line.rstrip('\n').split('\t')
            if hdr is None:
                if line.startswith('#') and any('Locus' in c for c in f):
                    hdr = [c.lstrip('#').strip() for c in f]
                continue
            if line.startswith('#') or len(f) < len(hdr):
                continue
            r = dict(zip(hdr, f))
            lt = (r.get('Locus Tag') or '').strip()
            if not lt:
                continue
            g = (r.get('Gene') or '').strip()
            p = (r.get('Product') or '').strip()
            if g and g != 'hypothetical protein':
                gene[lt] = g
            # only fill a MISSING product for a protein the FASTA has; never
            # introduce an id. Bakta's .tsv also lists tRNA/rRNA/ncRNA and
            # pseudogenes, which have a Product but no entry in the .faa, and
            # adding them inflated this step's background to 10,701 against
            # 36d's 10,542 for the identical set of genomes.
            if p and lt in prod and not prod[lt]:
                prod[lt] = p
if missing:
    say(f"  !! no proteome FASTA for: {', '.join(missing)}")

allids = set(prod)
hit = sum(1 for m in ctd for p in ctd[m] if p in prod)
say(f"  proteins with a description : {len(prod):,} in the CTD-bearing genomes")
say(f"  CTD proteins matched to one : {hit:,}/{NCTD:,}")
if hit < 0.9 * NCTD:
    say("  !! most CTD ids did not match a FASTA header. That is a join failure —")
    say("     the hmmsearch target names and the FASTA ids must be identical.")
    say("     Nothing below is usable.")
    dump_and_exit(0, 'JOIN_FAILED')
say("  [PASS] the CTD ids join onto the proteome they were searched in.")
say()

# ---------------------------------------------------------- classify
CLASSES = [
    ('carbohydrate-active', re.compile(
        r'glycosid|glycosyl hydrolase|glucanase|xylanase|xylosidase|amylase|'
        r'chitinase|cellulase|cellobio|pectat|pectin|mannanase|mannosidase|'
        r'galactosidase|glucosidase|arabino|fucosidase|sialidase|neuraminidase|'
        r'hexosaminidase|agarase|alginate|ulvan|laminarin|levanase|inulinase|'
        r'polysaccharide lyase|carbohydrate esterase|glucuronyl|rhamnosidase|'
        r'chondroitin|heparin|hyaluron|dextranase|pullulanase|glycoside', re.I)),
    ('sulfatase', re.compile(r'sulfatase|sulphatase', re.I)),
    ('peptidase / protease', re.compile(
        r'peptidase|protease|proteinase|collagenase|gingipain', re.I)),
    ('lipase / esterase', re.compile(
        r'lipase|esterase|phospholipase|cutinase|acylhydrolase', re.I)),
    ('nuclease', re.compile(r'nuclease|dnase|rnase|endonuclease', re.I)),
    ('SusC / SusD / TonB', re.compile(
        r'\bsusc\b|\bsusd\b|tonb-dependent|ton_b|starch-binding', re.I)),
    # Ig-like, PKD, cadherin, LamG/jellyroll, fibronectin-III and CARDB are
    # cell-surface adhesion MODULES. v1 listed only four of them, so most of the
    # surface proteins in the cargo fell into 'other, annotated'.
    ('adhesion / surface module', re.compile(
        r'adhesin|fibronectin|\bfn3\b|cadherin|hemagglutinin|s-layer|'
        r'surface layer|internalin|leucine-rich repeat|\bIg-like\b|'
        r'immunoglobulin|\bPKD\b|LamG|laminin_g|jellyroll|jelly roll|'
        r'thrombospondin|\bCARDB\b|\bbig_?[0-9]\b|invasin|intimin', re.I)),
    ('gliding motility', re.compile(r'gliding|\bspr[A-Z]\b|\bgld[A-Z]\b', re.I)),
    ('hypothetical / uncharacterised', re.compile(
        r'hypothetical|uncharacteri|\bDUF\d+|putative uncharacteri|unknown function', re.I)),
    # ---- CIRCULAR ANNOTATION, and it must be its own class ------------------
    # Bakta names a protein whose only recognisable feature is the C-terminal
    # sorting domain after that domain: "T9SS type A sorting domain-containing
    # protein", "Protein containing Por secretion system C-terminal sorting
    # domain", "Putative secreted protein (Por secretion system target)".
    # v1 filed all of these under 'other, annotated' — i.e. counted a protein
    # NAMED AFTER THE SEARCH as annotated evidence against the hypothesis, which
    # inflated the CTD denominator with proteins carrying no functional
    # information at all. They are functionally unannotated and are reported as
    # such. Placed after every functional class so a real function still wins.
    ('CTD-only (named for the domain)', re.compile(
        r't9ss[^a-z]*type|sorting domain[- ]containing|'
        r'por secretion system|secretion system c-terminal|'
        r'secretion system target', re.I)),
]

def classify(p):
    p = p or ''
    for name, rx in CLASSES:
        if rx.search(p):
            return name
    return 'other, annotated'

# a protein in one of these carries no functional information
UNINFORMATIVE = ('hypothetical / uncharacterised', 'CTD-only (named for the domain)')

ctdids = set().union(*ctd.values())
bg = allids - ctdids

cc = Counter(classify(prod.get(p, '')) for p in sorted(ctdids))
bc = Counter(classify(prod.get(p, '')) for p in sorted(bg))

# ---- what may NOT enter the composition test --------------------------------
# A protein named "T9SS type A sorting domain-containing protein" is in the CTD
# set BY CONSTRUCTION: the domain this step searched for is the only reason the
# annotator could name it. Testing that class against the background would
# manufacture a large p-value out of the tautology and prove nothing about
# cargo. It is displayed as its own descriptive row and excluded from every
# statistic. (Hypotheticals are NOT excluded — 'the cargo is less nameable than
# the genome' is a real difference, not a circular one.)
CIRCULAR = 'CTD-only (named for the domain)'
ctd_t = {p for p in ctdids if classify(prod.get(p, '')) != CIRCULAR}
bg_t = {p for p in bg if classify(prod.get(p, '')) != CIRCULAR}
cc_t = Counter(classify(prod.get(p, '')) for p in sorted(ctd_t))
bc_t = Counter(classify(prod.get(p, '')) for p in sorted(bg_t))
NC, NB = len(ctd_t), len(bg_t)

rule('1. WHAT THE CTD PROTEINS ARE')
say(f"  background = all {len(bg):,} other proteins in the SAME {len(ctd)} "
    f"CTD-bearing genome(s).")
say("  Protein ids come from the .faa hmmsearch ran on, and ONLY from it. Bakta's")
say("  .tsv is read for missing product names but may not introduce ids: it also")
say("  lists tRNA, rRNA and pseudogenes, which have a product and no protein.")
say("  Comparing against the whole catalogue would turn Bacteroidota's CAZyme")
say("  richness into fake enrichment; that is the step-19b error.")
say()
say(f"  {'functional class':<34}{'CTD':>7}{'% CTD':>8}"
    f"{'backgr':>9}{'% bg':>8}{'fold':>8}{'dir':>5}{'q (BH)':>11}")
say('-' * W)
order = sorted(set(list(cc) + list(bc)), key=lambda k: (-cc.get(k, 0), k))
tested = [k for k in order if k != CIRCULAR]
raw = [fisher_two_sided(cc_t.get(k, 0), NC - cc_t.get(k, 0),
                        bc_t.get(k, 0), NB - bc_t.get(k, 0)) for k in tested]
qs = dict(zip(tested, bh(raw)))
PVAL = {}
for k in order:
    a = cc.get(k, 0); c = bc.get(k, 0)
    if k == CIRCULAR:
        fa = a / max(1, len(ctdids)); fb = c / max(1, len(bg))
        # this row's percentages are of the FULL sets, not the tested ones —
        # marked, because a column with two denominators is a misreading waiting
        # to happen
        say(f"  {k:<34}{a:>7}{100 * fa:>6.1f}%†{c:>9,}{100 * fb:>6.1f}%†"
            f"{'  —':>8}{'  —':>5}{'circular':>11}")
        continue
    at, ct = cc_t.get(k, 0), bc_t.get(k, 0)
    fa = at / max(1, NC); fb = ct / max(1, NB)
    fold = (fa / fb) if fb else float('inf')
    q = qs[k]
    PVAL[k] = (at, fa, ct, fb, fold, q)
    fs = '    inf' if fold == float('inf') else f"{fold:>8.2f}"
    dr = '  up' if fold > 1 else ('down' if fold < 1 else '   =')
    star = '*' if q < 0.05 else ' '
    say(f"  {k:<34}{at:>7}{100 * fa:>7.1f}%{ct:>9,}{100 * fb:>7.1f}%{fs}{dr:>5}"
        f"{q:>10.2e}{star}")
say()
say(f"  † that row alone is a percentage of the FULL {len(ctdids)} / {len(bg):,};")
say(f"  every other row is a percentage of the {NC} CTD and {NB:,} background")
say("  proteins that are")
say("  NOT named for the sorting domain itself. That class is shown for")
say("  completeness and excluded from every test: a protein the annotator could")
say("  only name because it carries the domain this step searched for is in the")
say("  CTD set by construction, and testing it would prove the tautology.")
say()
say("  Two-sided Fisher per class, Benjamini-Hochberg across the classes. v1 used")
say("  a one-sided test, which cannot see a DEPLETION — and depletion to zero is")
say("  the strongest signal in this cargo set.")
say()
uninf_c = sum(cc.get(k, 0) for k in UNINFORMATIVE)
uninf_b = sum(bc.get(k, 0) for k in UNINFORMATIVE)
say(f"  functionally uninformative names (hypothetical + CTD-only):")
say(f"    CTD set   : {uninf_c}/{len(ctdids)}  ({100 * uninf_c / max(1, len(ctdids)):.1f}%)")
say(f"    background: {uninf_b:,}/{len(bg):,}  ({100 * uninf_b / max(1, len(bg)):.1f}%)")
say("  A cargo set that is mostly unnameable is a real property of Bacteroidota")
say("  secretomes, not a defect to smooth over — but it caps how much of the")
say("  cargo any annotation-based statement can cover.")
say()

# ---------------------------------------------------------- the control
rule('2. CONTROL — IS THE CARGO SET DIFFERENT FROM ITS OWN GENOMES?')
say("  A secretion system's cargo must not look like a random draw from the")
say("  genome that encodes it. If the composition above is indistinguishable")
say("  from the background, nothing about WHAT is exported is supportable —")
say("  whatever the counts.")
say()
# ---- the composition test proper -------------------------------------------
# v1 asked ONE pooled question: are 'degradative enzymes' over-represented? It
# merged peptidases (enriched 1.9x) with lipases and sulfatases (absent from the
# cargo entirely) and they cancelled to fold 1.00, p = 0.54 — reporting NO
# DIFFERENCE for a set whose composition is plainly different. A pooled binary
# can cancel; a composition test cannot. G-test (likelihood-ratio chi-square)
# over the whole class distribution, CTD vs background.
keys = [k for k in sorted(set(list(cc_t) + list(bc_t)))
        if k != CIRCULAR and (cc_t.get(k, 0) + bc_t.get(k, 0)) > 0]
G = 0.0
NT = NC + NB
low_exp = 0
for k in keys:
    o1, o2 = cc_t.get(k, 0), bc_t.get(k, 0)
    tot = o1 + o2
    e1, e2 = tot * NC / NT, tot * NB / NT
    if min(e1, e2) < 5:
        low_exp += 1
    if o1:
        G += 2 * o1 * math.log(o1 / e1)
    if o2:
        G += 2 * o2 * math.log(o2 / e2)
DF = len(keys) - 1
PG = chi2_sf(G, DF)
say(f"  classes compared          : {len(keys)}")
say(f"  G (likelihood-ratio chi2) : {G:.1f} on {DF} df")
say(f"  p                         : {PG:.2e}")
if low_exp:
    say(f"  note: {low_exp} class(es) have an expected count < 5; the per-class")
    say("        Fisher tests above are exact and are what should be quoted for")
    say("        those, the G-test is the overall statement.")
say()
sig = [(k, PVAL[k]) for k in sorted(keys) if PVAL.get(k) and PVAL[k][5] < 0.05]
up = [k for k, v in sig if v[4] > 1]
dn = [k for k, v in sig if v[4] < 1]
if PG < 1e-3 and sig:
    VERDICT = 'PASS'
    say("  [PASS] the cargo composition differs from the background of the same")
    say("         genomes. The CTD set is not a random draw from its own genomes,")
    say("         so statements about what is exported are supportable — for the")
    say("         classes that individually survive correction, and no others.")
    say()
    if up:
        say(f"    enriched : {', '.join(up)}")
    if dn:
        say(f"    depleted : {', '.join(dn)}")
elif PG < 1e-3:
    VERDICT = 'NOCLASS'
    say("  [PARTIAL] the overall composition differs, but no single class survives")
    say("         multiple-testing correction. Report that the cargo is")
    say("         compositionally distinct and name no class.")
else:
    VERDICT = 'FAIL'
    say("  [FAIL] the CTD set's composition is not distinguishable from its own")
    say("         genomes' background. Report the COUNT of cargo proteins (step 36")
    say("         gap B) and make NO claim about what is exported.")
say()
# ---- the pooled test v1 used, kept and labelled -----------------------------
DEGRAD = ('carbohydrate-active', 'sulfatase', 'peptidase / protease',
          'lipase / esterase', 'nuclease')
a = sum(cc_t.get(k, 0) for k in DEGRAD)
c = sum(bc_t.get(k, 0) for k in DEGRAD)
fa = a / max(1, NC); fb = c / max(1, NB)
say("  For the record, the pooled binary v1 used, which is why it said FAIL:")
say(f"    'degradative enzymes' CTD {a}/{NC} ({100 * fa:.1f}%) vs "
    f"background {c:,}/{NB:,} ({100 * fb:.1f}%), fold "
    f"{(fa / fb) if fb else float('inf'):.2f}")
say("    Pooling classes that move in OPPOSITE directions cancels them. That is")
say("    a defect of the statistic, not a property of the cargo.")
say()

# --------------------------------------------------- per genome
rule('3. PER GENOME')
gen = {}
tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
alias = {}
for path in (os.path.join(F, 'MAG_CROSSWALK.tsv'), tp):
    if not os.path.exists(path):
        continue
    with open(path, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            ids = {(v or '').strip() for v in r.values()
                   if v and ' ' not in str(v) and ';' not in str(v)
                   and IDTOK.match(str(v).strip())
                   and not re.fullmatch(r'[0-9._-]+', str(v).strip())
                   and not str(v).strip().startswith(('d__', 'p__', 'c__', 'o__',
                                                      'f__', 'g__', 's__'))}
            g = ''
            for c in r:
                if c and c.strip().lower() == 'genus':
                    g = (r.get(c) or '').strip()
            if not g:
                m = re.search(r'g__([^;\t]+)', ' '.join(str(v or '') for v in r.values()))
                g = m.group(1).strip() if m else ''
            for i in sorted(ids):
                alias[i] = set(ids)    # row-local; never merged
                if g:
                    gen.setdefault(i, g)
say(f"  {'genome':<26}{'genus':<20}{'CTD':>6}{'CAZy':>7}{'sulf':>6}"
    f"{'pept':>6}{'SusCD':>7}{'hypo':>6}")
say('-' * W)
for m in sorted(ctd, key=lambda x: (-len(ctd[x]), x)):
    cl = Counter(classify(prod.get(p, '')) for p in ctd[m])
    g = gen.get(m) or next((gen[a2] for a2 in sorted(alias.get(m, ())) if a2 in gen), '-')
    say(f"  {m:<26}{g[:19]:<20}{len(ctd[m]):>6}"
        f"{cl.get('carbohydrate-active', 0):>7}{cl.get('sulfatase', 0):>6}"
        f"{cl.get('peptidase / protease', 0):>6}"
        f"{cl.get('SusC / SusD / TonB', 0):>7}"
        f"{cl.get('hypothetical / uncharacterised', 0):>6}")
say()

# ------------------------------------------- optional per-protein CAZyme table
rule('4. CROSS WITH A PER-PROTEIN CAZyme TABLE, IF ONE EXISTS')
_d36 = os.path.join(D, 'cazy_ctd.tsv')
if os.path.exists(_d36):
    say("  >> STEP 36d HAS ALREADY DONE THIS CROSS, and its result SUPERSEDES the")
    say("     carbohydrate-active row in section 1 above, which came from Bakta")
    say("     product strings. 36d used dbCAN family HMMs and found 53 cargo")
    say("     CAZymes where the keyword method found 8 — because Bakta names many")
    say("     of them after the sorting domain rather than the catalytic one.")
    say("     It also ran two controls this step cannot: a co-location test")
    say("     (is the CAZy domain the sorting domain?) and a per-genome")
    say("     stratification. See 36_gapclose/CAZY_CROSS.txt.")
    say(f"     table: {os.path.relpath(_d36, WS2)}")
    say()
    say()
# 🔴 v2 widened the sniff and then SILENTLY PICKED THE FIRST of five "usable"
# tables -- which was 23_extras/pul/nonPUL_cazymes.tsv, a NON-PUL SUBSET. It
# reported 7/209 (3.3%) as if it were the answer, contradicting 36d's 53/254
# from the full set. A search that finds five candidates and quietly takes one
# is worse than a search that finds none, because the output looks like a
# result. Two rules now: prefer 36d's table outright, and if several others
# qualify, REFUSE to choose.
_d36 = os.path.join(D, 'cazy_ctd.tsv')
if os.path.exists(_d36):
    czp, fam_of, ctd_cz, bg_cz = set(), {}, set(), set()
    with open(_d36, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            pid = (r.get('protein') or '').strip()
            fam = (r.get('cazy_family') or '').strip()
            if not pid or not fam:
                continue
            czp.add(pid); fam_of[pid] = fam
            (ctd_cz if (r.get('ctd') or '').strip() == 'yes' else bg_cz).add(pid)
    say(f"  source: {os.path.relpath(_d36, WS2)}  (written by step 36d, dbCAN family HMMs)")
    say()
    say(f"  cargo proteins carrying a CAZy family     : {len(ctd_cz):,}/{len(ctdids):,}"
        f"   ({100.0 * len(ctd_cz) / max(1, len(ctdids)):.1f}%)")
    say(f"  same genomes, non-cargo, carrying one     : {len(bg_cz):,}/{len(bg):,}"
        f"   ({100.0 * len(bg_cz) / max(1, len(bg)):.1f}%)")
    say()
    say("  Compare with section 1's carbohydrate-active row, which is the KEYWORD")
    say(f"  count: {cc.get('carbohydrate-active', 0)} of the same cargo proteins. The"
        f" difference is the")
    say("  measurement 36d exists to make, and 36d's number is the one to quote —")
    say("  it also carries the co-location control and the per-genome stratification.")
    if ctd_cz:
        say()
        say("  families among the cargo CAZymes:")
        for f2, n in sorted(Counter(fam_of[p2] for p2 in sorted(ctd_cz)).items(), key=lambda kv: (-kv[1], kv[0]))[:15]:
            say(f"    {f2:<12}{n:>5}")
else:
    say("  36d has not been run (no 36_gapclose/cazy_ctd.tsv).")
    say()
    say("  23d's cazyme_endo_exo.tsv is a PER-GENOME summary — mag_id, genus,")
    say("  phylum, completeness, n_endo, n_exo, exo_per_endo, role. No protein ids,")
    say("  so it cannot be intersected with a protein set; step 36 v3 printed")
    say("  '0 CAZymes carry a CTD' from that empty intersection and that number")
    say("  never meant anything. Looking for a table that carries BOTH a CAZy")
    say("  family and a protein id:")
    say()
    cands = []
    for pat in ('*cazy*', '*dbcan*', '*CAZy*', '*overview*', '*hmmer*'):
        for root in (os.path.join(WS2, '30_cazy'), os.path.join(WS2, '23_extras'),
                     T9, os.path.join(WS2, '30_polymer'), WS2):
            cands += glob.glob(os.path.join(root, '**', pat + '.tsv'), recursive=True)
            cands += glob.glob(os.path.join(root, '**', pat + '.txt'), recursive=True)
    cands = sorted(set(cands))
    FAM = re.compile(r'^(GH|GT|PL|CE|AA|CBM)\d+')
    usable = []
    for p2 in cands[:60]:
        try:
            with open(p2, errors='replace') as fh:
                head = [l for _, l in zip(range(2000), fh)]
        except Exception:
            continue
        ids, fams = set(), 0
        for line in head[1:]:
            for v in line.rstrip('\n').split('\t'):
                v = v.strip()
                if FAM.match(v):
                    fams += 1
                elif v and re.fullmatch(r'[A-Za-z0-9_.|-]{6,}', v) and '_' in v:
                    ids.add(v)
        shared = {re.split(r'[_.|]', x)[0] for x in ids} & \
                 {re.split(r'[_.|]', x)[0] for x in ctdids}
        ok = bool(fams and shared)
        if ok:
            usable.append(p2)
        say(f"  {os.path.relpath(p2, WS2)[:70]:<72}"
            f"{'USABLE' if ok else ('no protein ids' if fams else 'no CAZy families')}")
    say()
    if not usable:
        say("  >> NOT TESTED. No table carries both. Recorded as not tested, never")
        say("     as 'no CAZyme carries a CTD'. Run step 36d, which derives the")
        say("     per-protein calls directly from dbCAN.")
    elif len(usable) > 1:
        say(f"  >> NOT TESTED — {len(usable)} candidate tables qualify and they are NOT")
        say("     interchangeable. Some are SUBSETS (nonPUL_cazymes.tsv holds only")
        say("     non-PUL CAZymes; secreted_cazymes.tsv only predicted-secreted ones),")
        say("     so picking one silently would produce a fraction with a denominator")
        say("     nobody chose. Run step 36d instead; it derives the full set.")
    else:
        say(f"  >> One candidate: {os.path.relpath(usable[0], WS2)}. Even so, prefer")
        say("     step 36d — it applies dbCAN's own parser rule and adds two controls.")
say()


rule('5. THE LARGEST CARGO PROTEINS, BY GENOME')
for m in sorted(ctd, key=lambda x: (-len(ctd[x]), x)):
    named = [(p, prod.get(p, '')) for p in sorted(ctd[m])
             if prod.get(p) and 'hypothetical' not in prod.get(p, '').lower()]
    say(f"  {m}  —  {len(named)} of {len(ctd[m])} CTD proteins have a product name")
    for p, d in named[:12]:
        say(f"    {p:<26}{d[:66]}")
    say()

with open(os.path.join(D, 'ctd_identity.tsv'), 'w') as oh:
    oh.write('protein\tmag\tgenus\tclass\tproduct\n')
    for m in sorted(ctd):
        g = gen.get(m) or next((gen[a2] for a2 in sorted(alias.get(m, ())) if a2 in gen), '')
        for p in sorted(ctd[m]):
            d = prod.get(p, '')
            oh.write(f"{p}\t{m}\t{g}\t{classify(d)}\t{d}\n")

STATUS[0] = {'PASS': 'OK', 'NOCLASS': 'PARTIAL'}.get(VERDICT, VERDICT)
rule('WHAT TO REPORT')
if VERDICT == 'PASS':
    say("  1. The T9SS cargo of these Bacteroidota is COMPOSITIONALLY DISTINCT from")
    say(f"     the background of the same genomes (G = {G:.0f}, {DF} df, p = {PG:.1e}).")
    say("     Name only the classes with q < 0.05 in section 1, and give the")
    say("     direction — a class depleted to zero is as much a result as an")
    say("     enriched one, and says the domain is not simply tagging abundant")
    say("     protein families.")
    say("  2. Quote the uninformative fraction. A cargo set that is mostly")
    say("     hypothetical or named only for its own sorting domain caps what any")
    say("     annotation-based statement can cover, and hiding it inflates the")
    say("     apparent reach of the claim.")
elif VERDICT == 'NOCLASS':
    say("  1. The cargo is compositionally distinct from the background of the same")
    say(f"     genomes (G = {G:.0f}, {DF} df, p = {PG:.1e}) but no individual class")
    say("     survives correction. Report the difference; name no class.")
elif VERDICT == 'FAIL':
    say(f"  1. Report step 36 gap B's COUNT ({NCTD:,} CTD proteins in {len(ctd)}")
    say("     machinery-bearing genomes) and make no claim about cargo identity.")
    say("     The composition test failed and that is the honest stopping point.")
else:
    say("  1. Too few annotated degradative enzymes to test enrichment. Report the")
    say("     composition table descriptively, without a claim.")
say("  3. Every count here is a FLOOR twice over: hmmsearch finds what TIGR04183")
say("     and TIGR04131 match, and Bakta names what it can name.")
say("  4. All three genomes are culture-derived (LIQ). Quote their step-26")
say("     sediment breadth with any in-situ statement about degradation.")

print()
print(f"  written: {os.path.join(D, 'ctd_identity.tsv')}")
dump_and_exit(0)          # writes CARGO.txt AND the status line
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "cargo report written" "$([ -s "$D/CARGO.txt" ] && echo 0 || echo 1)" || OK=1
check "identity table written" "$([ -s "$D/ctd_identity.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/CARGO.txt" 2>/dev/null | awk '{print $2}')
echo "  status line: ${ST:-<none>}"
case "${ST:-MISSING}" in
  OK)      ;;
  PARTIAL) echo "  [WARN] composition differs but no class survived correction." ;;
  *)       echo "  [STOP] ${ST:-no status line} — nothing here is quotable."; OK=1 ;;
esac
[ "$OK" -eq 0 ] && runlog "36c_cargo complete" || runlog "**36c_cargo FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
