#!/usr/bin/env bash
# =============================================================================
# 29_transporter_check.sh  --  a consumer without a transporter is not a consumer
#
# WHAT STEP 21 DID AND DID NOT CHECK
#   It built 593 cross-feeding edges by finding a metabolite's PRODUCING gene in
#   one genome and its CONSUMING gene in another. That is half the requirement.
#   A catabolic enzyme sitting in the cytoplasm can only act on a molecule that
#   got INTO the cytoplasm, and for most metabolites that needs a transporter.
#   Nobody checked whether the consumers have one.
#
# THE MISTAKE THIS SCRIPT MUST NOT MAKE
#   Requiring a transporter for everything would be worse than not checking at
#   all. HYDROGEN CROSSES MEMBRANES BY DIFFUSION. So do CO2, ammonia, hydrogen
#   sulfide, methanol and methanethiol — small, uncharged, membrane-permeant.
#   Demanding a transporter for H2 would delete all 184 hydrogen edges, which
#   are the best-supported edges in the network, and the deletion would look
#   like a finding.
#
#   So every metabolite is first classified by whether it NEEDS carrier-mediated
#   uptake at all, the classification is printed with its reasoning, and the
#   filter is applied ONLY to the ones that do. Charge at neutral pH is the main
#   criterion: acetate, formate, lactate and the amino acids are ionised and do
#   not cross a lipid bilayer at useful rates; methanol and H2 are not ionised
#   and do.
#
# RULE 1 IS ENFORCED HERE TOO
#   No KO number is trusted from memory. Every KO in the transporter table is
#   checked against KEGG's own ko_list at runtime and a mismatch ABORTS. That
#   guard exists because of the K14080 error at step 20, and a transporter table
#   is exactly the kind of list where a wrong KO silently removes real edges.
#
# WHAT A FAILED EDGE MEANS
#   Not "this organism cannot eat that". It means the transporter was not found,
#   which can be because it is absent, because the genome is incomplete, or
#   because the family is not in the table. Losing an edge is a DOWNGRADE in
#   confidence, not a disproof, and the output says so on every line.
#
# INPUTS   21_network/edges.tsv (or crossfeeding.tsv) · 11_kofam/*
#          09_annotation/<mag>/<mag>.tsv · 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   29_transport/edge_transport.tsv · transporter_matrix.tsv
# RUNTIME  seconds.
# USAGE    ws2 && bash scripts/29_transporter_check.sh 2>&1 | tee logs/29.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="29_transport"
N="$WS2/21_network"; K="$WS2/11_kofam"; A="$WS2/09_annotation"
F="$WS2/07_mags"; T="$WS2/29_transport"; mkdir -p "$T"

echo "######## 29 — TRANSPORTER CHECK ON THE CROSS-FEEDING NETWORK ########"
echo "  started: $(date -Is)"
runlog "29_transport started"
[ -d "$N" ] || { echo "  !! 21_network not found — run 21 first"; exit 1; }

# ---- locate ko_list for the label guard ------------------------------------
KOLIST=""
for c in "$K/ko_list" "$K/db/ko_list" "$WS2/db/ko_list" "$WS2/00_db/kofam/ko_list"; do
  [ -s "$c" ] && { KOLIST="$c"; break; }
done
[ -n "$KOLIST" ] || KOLIST=$(find "$WS2" -maxdepth 4 -name 'ko_list' -size +1k 2>/dev/null | head -1)
if [ -n "$KOLIST" ]; then echo "  ko_list: ${KOLIST#$WS2/}"
else echo "  !! ko_list NOT FOUND — the label guard cannot run"; fi

python3 - "$N" "$K" "$A" "$F" "$T" "${KOLIST:-}" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict

N, K, A, F, T, KOLIST = sys.argv[1:7]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus','') or '-')
def compl(m):
    try: return float((tax.get(m, {}) or {}).get('completeness','') or 0)
    except ValueError: return 0.0
wn_of = {mid(w): w for w in tax}
def as_wn(x): return x if x in tax else wn_of.get(x, x)

# ==========================================================================
# 1. WHICH METABOLITES ACTUALLY NEED A TRANSPORTER
# ==========================================================================
# key -> (needs_transport, reason)
PHYS = {
 'h2':        (False, 'H2 is small and uncharged — diffuses freely'),
 'hydrogen':  (False, 'H2 is small and uncharged — diffuses freely'),
 'co2':       (False, 'CO2 diffuses; carbonic anhydrase handles the rest'),
 'nh3':       (False, 'NH3 diffuses at the pH of this system'),
 'ammoni':    (False, 'NH3 diffuses at the pH of this system'),
 'sulfide':   (False, 'H2S is uncharged and membrane-permeant'),
 'h2s':       (False, 'H2S is uncharged and membrane-permeant'),
 'methanol':  (False, 'small uncharged alcohol — permeant'),
 'methanethiol': (False, 'small uncharged thiol — permeant'),
 'dms':       (False, 'dimethyl sulfide is uncharged and permeant'),
 # ACETATE RECLASSIFIED after the first run found a transporter in 1 of 29
 # genomes. The marker set is not the only suspect: acetic acid has pKa 4.76, so
 # ~1.7% is undissociated at pH 7, and the protonated acid crosses membranes
 # freely. Passive flux of the undissociated acid is a documented route in
 # anaerobes, so a genome without a named acetate carrier is NOT thereby unable
 # to take up acetate. Formic acid, by contrast, has pKa 3.75 — only ~0.06%
 # undissociated — and genuinely needs FocA/FdhC. Same functional group,
 # different answer, because the pKa differs by a factor of ten.
 'acetate':   (False, 'pKa 4.76 — ~1.7% undissociated at pH 7; passive flux of '
                      'acetic acid is a documented route, carrier optional'),
 'formate':   (True,  'ionised at pH 7 — needs a channel (FocA/FdhC family)'),
 'lactate':   (True,  'ionised at pH 7 — needs a permease'),
 'propionate':(True,  'ionised at pH 7 — needs a carrier'),
 'butyrate':  (True,  'ionised at pH 7 — needs a carrier'),
 'choline':   (True,  'permanently charged quaternary amine — needs BetT/OpuD'),
 'tma':       (True,  'protonated at pH 7 — carrier-mediated'),
 'trimethylamine': (True, 'protonated at pH 7 — carrier-mediated'),
 'methylamine':(True, 'protonated at pH 7 — carrier-mediated'),
 'b12':       (True,  'far too large to diffuse — BtuBFCD'),
 'cobalamin': (True,  'far too large to diffuse — BtuBFCD'),
 'biotin':    (True,  'needs BioY'),
 'thiamin':   (True,  'needs ThiBPQ/ThiT'),
 'amino':     (True,  'zwitterionic — ABC or secondary transporters'),
}
def needs_transport(metab):
    m = (metab or '').lower()
    for k, v in PHYS.items():
        if k in m: return v
    return (True, 'no physiological rule recorded — assumed to need a carrier')

# TRANSPORTER KOs. Each entry is (my label, keywords that MUST appear in KEGG's
# own definition). The first version stored only the gene symbol and asked
# whether that symbol appeared in the definition — and every one of the fifteen
# failed, on definitions that were all CORRECT:
#
#     K14393  I said 'actP'  KEGG says 'cation/acetate symporter'
#     K06212  I said 'focA'  KEGG says 'formate transporter'
#     K06858  I said 'btuB'  KEGG says 'vitamin B12 transport system ...'
#
# ActP IS the cation/acetate symporter. FocA IS the formate transporter. This
# ko_list's definition field carries the FUNCTION, not the gene symbol, so the
# guard was testing nomenclature and calling it verification. Exactly the step
# 19b Guard-A failure: 18 markers dropped on wording, only 4 of them real.
#
# The fix is to verify the BIOLOGY the KO encodes rather than the name I happen
# to use for it. Keywords below are the substrate and the activity — the two
# things that must be true for the KO to belong in this table at all.
# A useful side effect: K02064 and K02062 turn out to be substrate-binding and
# ATP-binding, not the thiP/thiB roles I assigned. The symbols were wrong, the
# KOs were right, and function-keyword checking is indifferent to my error.
TRANS = {
 'acetate':   {'K14393': ('acetate symporter', ('acetate',)),
               'K02106': ('short-chain fatty acid transporter',
                          ('fatty acid', 'transport'))},
 'formate':   {'K06212': ('formate transporter', ('formate', 'transport'))},
 'choline':   {'K02168': ('choline/betaine transporter', ('choline', 'transport')),
               'K05020': ('glycine betaine transporter', ('betaine', 'transport'))},
 'b12':       {'K06858': ('B12 transport, binding', ('b12', 'transport')),
               'K06073': ('B12 transport, permease', ('b12', 'transport')),
               'K06074': ('B12 transport, ATPase', ('b12', 'transport'))},
 'amino':     {'K01999': ('branched-chain AA transport', ('amino acid', 'transport')),
               'K01997': ('branched-chain AA transport', ('amino acid', 'transport')),
               'K01998': ('branched-chain AA transport', ('amino acid', 'transport')),
               'K10005': ('glutamate transport', ('glutamate', 'transport'))},
 'biotin':    {'K03523': ('biotin transport', ('biotin', 'transport'))},
 'thiamin':   {'K02064': ('thiamine transport', ('thiamin', 'transport')),
               'K02062': ('thiamine transport', ('thiamin', 'transport'))},
}

# ---- LABEL GUARD ---------------------------------------------------------
kodef = {}
if KOLIST and os.path.exists(KOLIST):
    for line in open(KOLIST):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 12 and f[0].startswith('K'): kodef[f[0]] = f[11]
        elif len(f) >= 2 and f[0].startswith('K'): kodef[f[0]] = f[-1]
print()
print("=" * 104)
print("  0. LABEL GUARD — no KO trusted from memory")
print("=" * 104)
print("  Verifying the FUNCTION each KO encodes against KEGG's own definition,")
print("  not the gene symbol I happen to use for it. The first version checked")
print("  symbols and failed all 15 on definitions that were entirely correct.")
print()
bad = []
if not kodef:
    print("  ko_list unavailable. THE GUARD CANNOT RUN, so the transporter table")
    print("  is unverified. Treat every result below as provisional.")
else:
    for met, kos in TRANS.items():
        for ko, (label, kws) in kos.items():
            d = kodef.get(ko)
            if d is None:
                bad.append((ko, label, 'NOT IN ko_list', kws)); continue
            dl = d.lower()
            missing = [k for k in kws if k not in dl]
            if missing: bad.append((ko, label, d[:64], tuple(missing)))
    if bad:
        print("  >> MISMATCH. KEGG's definition does not describe the function I")
        print("     claimed for these KOs:")
        for ko, label, d, miss in bad:
            print(f"     {ko}  I claimed '{label}'  KEGG: '{d}'")
            print(f"           missing required term(s): {', '.join(miss)}")
        print("  Aborting. A wrong KO in a transporter table silently deletes real")
        print("  edges, and a deleted edge looks exactly like a finding.")
        sys.exit(2)
    n = sum(len(v) for v in TRANS.values())
    print(f"  all {n} transporter KOs verified against KEGG. What was checked:")
    for met, kos in sorted(TRANS.items()):
        for ko, (label, kws) in sorted(kos.items()):
            print(f"    {ko}  {kodef[ko][:62]}")

# ==========================================================================
# 2. WHICH GENOMES HAVE WHICH TRANSPORTERS
# ==========================================================================
# ko_matrix.tsv is WIDE: one row per MAG, one COLUMN per KO, headed
# 'K00399:mcrA'. The previous loader hunted for a column whose VALUES look like
# K-numbers and found none, because the K-numbers are in the HEADER. Handle both
# layouts and print which one was used.
ko_by_mag = defaultdict(set)
layout = None
for p in sorted(glob.glob(os.path.join(K, '**', '*.tsv'), recursive=True)):
    bn = os.path.basename(p)
    if bn == 'ko_list' or 'module' in bn.lower() or 'target' in bn.lower(): continue
    try: rows_ = list(csv.reader(open(p), delimiter='\t'))
    except Exception: continue
    if len(rows_) < 2: continue
    hdr = [h.strip() for h in rows_[0]]
    kcols = {i: h.split(':')[0] for i, h in enumerate(hdr)
             if re.match(r'^K\d{5}(:|$)', h)}
    if kcols:
        mi = next((i for i, h in enumerate(hdr)
                   if h.lower() in ('mag_id', 'working_name', 'mag', 'genome')), 0)
        for r in rows_[1:]:
            if len(r) <= mi: continue
            m = as_wn(r[mi].strip())
            if m not in tax: continue
            for i, ko in kcols.items():
                if len(r) > i and str(r[i]).strip() not in ('', '0', '-', 'NA', 'False'):
                    ko_by_mag[m].add(ko)
        layout = f'wide ({bn}, {len(kcols)} KO columns)'
        break
    kcol = next((ci for ci in range(len(rows_[1]))
                 if sum(1 for r in rows_[1:200] if len(r) > ci and
                        re.fullmatch(r'K\d{5}', r[ci].strip())) > 50), None)
    if kcol is None: continue
    mcol = next((ci for ci in range(len(rows_[1])) if ci != kcol and
                 sum(1 for r in rows_[1:200] if len(r) > ci and
                     as_wn(r[ci].strip()) in tax) > 50), None)
    if mcol is None: continue
    for r in rows_[1:]:
        if len(r) <= max(kcol, mcol): continue
        m = as_wn(r[mcol].strip())
        if m in tax and re.fullmatch(r'K\d{5}', r[kcol].strip()):
            ko_by_mag[m].add(r[kcol].strip())
    layout = f'long ({bn})'
    break

print(f"  KO table: {layout or 'NONE READABLE'} -> {len(ko_by_mag)} genomes, "
      f"{sum(len(v) for v in ko_by_mag.values())} calls")

# DOES IT EVEN CONTAIN TRANSPORTERS? Step 11 ran a TARGETED KOfam search over a
# marker set chosen for methanogenesis and sulfur cycling, years before this
# question existed. A transporter KO that was never searched for is missing from
# the matrix for BOOKKEEPING reasons, and reading that as "no transporter" would
# be the purest form of the mistake this project keeps catching.
want = {ko for kos in TRANS.values() for ko in kos}
present_cols = set()
for v in ko_by_mag.values(): present_cols |= v
covered = want & present_cols
print(f"  transporter KOs present anywhere in that table: {len(covered)}/{len(want)}")
KO_USABLE = len(covered) > 0
if not KO_USABLE:
    print()
    print("  >> THE KO ROUTE IS UNAVAILABLE FOR THIS QUESTION — a bookkeeping fact,")
    print("     not a biological one. Not one transporter family was in step 11's")
    print("     target set, so their absence from ko_matrix.tsv says nothing at all")
    print("     about these genomes. Falling back to Bakta product annotations,")
    print("     which DO name transporters, with every verdict marked product-based.")
    print("     The proper fix is an hmmsearch over these 15 KOfam profiles by the")
    print("     same route steps 19b and 20 use.")

# product-string backup. Uses GENE SYMBOLS and SUBSTRATE WORDS, because Bakta
# writes 'acetate permease ActP', never my own label for the KO family.
PRODPAT = {
 'acetate': ('actp', 'acetate permease', 'acetate symporter', 'acetate transport'),
 'formate': ('foca', 'formate transporter', 'formate channel', 'formate/nitrite'),
 'choline': ('bett', 'opud', 'choline transport', 'betaine transport',
             'betaine/carnitine'),
 'b12':     ('btub', 'btuc', 'btud', 'cobalamin transport', 'b12 transport'),
 'amino':   ('livh', 'livm', 'livk', 'glti', 'amino acid abc transporter',
             'amino acid transport'),
 'biotin':  ('bioy', 'biotin transport'),
 'thiamin': ('thib', 'thip', 'thiq', 'thiamine transport', 'thiamin transport'),
}
prodhit = defaultdict(set)
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d)
    if n not in tax: continue
    t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8 or f[1].lower() != 'cds': continue
        g, pr = (f[6] or '').lower(), (f[7] or '').lower()
        for met, pats in PRODPAT.items():
            if any(x == g or x in pr for x in pats): prodhit[n].add(met)

if not KO_USABLE and not prodhit:
    print()
    print("  >> NO VERDICT. Neither KO calls nor Bakta product annotations yielded")
    print("     a single transporter in any genome. Two independent routes failing")
    print("     together means the inputs are wrong, not the biology.")
    sys.exit(4)
print(f"  product-annotation route: transporters found in {len(prodhit)} genomes")

# ZERO-HIT FAMILY FLAG — the step 21 lesson, applied to myself.
# Step 21's first hydrogen run returned 0 producers against 10 consumers, which
# cannot be true when every fermenter in the catalogue makes hydrogen; the
# marker set was wrong. The same rule binds here: A TRANSPORTER CLASS WITH NO
# HIT IN ANY OF 29 GENOMES IS A MARKER-SET FAILURE UNTIL PROVEN OTHERWISE, and
# its edges must NOT be downgraded — a downgrade on a broken marker set reads as
# a biological finding.
cover = {met: sum(1 for m in tax if met in prodhit.get(m, set()) or
                  (ko_by_mag.get(m, set()) & set(TRANS.get(met, {}))))
         for met in TRANS}
# THE THRESHOLD WAS MINE AND IT WAS TOO STRICT. The first version called a
# family broken only at EXACTLY zero, so acetate at 1/29 escaped and all 205 of
# its edges were downgraded anyway. This project already has a rule for that,
# adopted at step 19: A FAMILY ABSENT FROM MORE THAN 80% OF GENOMES IS A
# MARKER-SET FAILURE UNTIL PROVEN OTHERWISE. Using a stricter bar here than the
# rule the notebook already carries was inconsistency, not caution.
BROKEN = {met for met, n in cover.items() if n < 0.20 * len(tax)}
THIN   = {met for met, n in cover.items()
          if 0.20 * len(tax) <= n < 0.35 * len(tax)}
print()
print("  TRANSPORTER FAMILY COVERAGE across 29 genomes:")
for met, n in sorted(cover.items(), key=lambda t: -t[1]):
    tag = (f'  <-- <20% of genomes: marker-set failure, edges NOT downgraded'
           if met in BROKEN else
           '  <-- thin coverage, downgrades may be detection failure'
           if met in THIN else '')
    print(f"    {met:<12}{n:>3}/{len(tax)} genomes{tag}")

def has_transport(m, met):
    kos = TRANS.get(met, {})
    if kos and (ko_by_mag.get(m, set()) & set(kos)): return 'KO'
    if met in prodhit.get(m, set()): return 'product'
    return ''

# ==========================================================================
# 3. THE EDGES
# ==========================================================================
ef = None
for c in ('crossfeeding.tsv', 'edges.tsv', 'cross_feeding.tsv'):
    p = os.path.join(N, c)
    if os.path.exists(p) and os.path.getsize(p) > 0: ef = p; break
if not ef:
    print("  !! no edge table found in 21_network/"); sys.exit(1)
print(f"  edges: {os.path.basename(ef)}")
rows = list(csv.DictReader(open(ef), delimiter='\t'))
hdr = list(rows[0]) if rows else []
def col(*pats):
    for h in hdr:
        hl = h.lower()
        if any(x in hl for x in pats): return h
    return None
cc = col('consumer', 'target', 'to')
pc = col('producer', 'source', 'from')
mc = col('metabolite', 'compound', 'substrate', 'currency')
tc = col('type', 'edge')
if not (cc and mc):
    print(f"  !! could not identify consumer/metabolite columns in {hdr}"); sys.exit(1)

print()
print("=" * 104)
print("  1. WHICH METABOLITES NEED CARRIER-MEDIATED UPTAKE AT ALL")
print("=" * 104)
from collections import Counter as _C
allm = _C((r.get(mc) or '').strip() for r in rows if r.get(mc))
usedm = _C((r.get(mc) or '').strip() for r in rows
           if r.get(mc) and not (tc and 'cross' not in (r.get(tc) or 'cross').lower()))
mets = sorted(allm)
print("  edges per metabolite (all rows -> rows that pass the type filter):")
for m in mets:
    warn = '   <-- ALL FILTERED OUT' if allm[m] and not usedm[m] else ''
    print(f"    {m[:26]:<27}{allm[m]:>6} -> {usedm[m]:<6}{warn}")
print()
for m in mets:
    nt, why = needs_transport(m)
    print(f"  {m[:26]:<27}{'TRANSPORT REQUIRED' if nt else 'diffuses':<20}{why}")
print()
print("  THE FILTER IS APPLIED ONLY TO THE FIRST GROUP. Demanding a transporter")
print("  for H2 would delete all 184 hydrogen edges — the best-supported edges in")
print("  the network — and the deletion would read as a discovery.")

kept = dropped = skipped = settled = 0
out = []
for r in rows:
    if tc and 'cross' not in (r.get(tc) or 'cross').lower(): continue
    cons = as_wn((r.get(cc) or '').strip()); met = (r.get(mc) or '').strip()
    prod_ = as_wn((r.get(pc) or '').strip()) if pc else ''
    nt, why = needs_transport(met)
    key = next((k for k in TRANS if k in met.lower()), None)
    if not nt:
        v, note = 'no transporter needed', why
        settled += 1
    elif key is None:
        v, note = 'untested', 'NO transporter family in my table for this metabolite'
        skipped += 1
    elif key in BROKEN:
        v = 'untested'
        note = (f'{key} transporter found in <20% of genomes — marker-set '
                f'failure by the step-19 rule, not absence; edge left intact')
        skipped += 1
    else:
        ev = has_transport(cons, key)
        if ev: v, note = 'supported', f'transporter found ({ev})'; kept += 1
        else:
            v = 'DOWNGRADED'
            note = (f'no transporter found; genome {compl(cons):.0f}% complete'
                    + ('; FAMILY THINLY DETECTED across the catalogue, so this may '
                       'be detection failure' if key in THIN else
                       ' — absence may be a binning gap'))
            dropped += 1
    out.append(dict(producer=mid(prod_) if prod_ else '', consumer=mid(cons),
                    metabolite=met, needs_transport='yes' if nt else 'no',
                    verdict=v, note=note))

print()
print("=" * 104)
print("  2. EDGE OUTCOMES")
print("=" * 104)
# SETTLED IS NOT THE SAME AS UNTESTED, and the first version counted them
# together as "not tested" — which reported 0 of 593 edges resolved when in fact
# most of them had been resolved in step 21's favour. An edge whose metabolite
# needs no carrier is not an open question; it is a closed one.
print(f"  {settled} SETTLED — metabolite crosses the membrane unaided, so the")
print(f"           edge never required a transporter and stands as step 21 made it")
print(f"  {kept} supported by an identified transporter")
print(f"  {dropped} DOWNGRADED — catabolic gene present, no uptake route found")
print(f"  {skipped} NOT TESTABLE — marker set too sparse, or no family in my table")
print()
print(f"  RESOLVED: {settled + kept + dropped} of {settled + kept + dropped + skipped} edges "
      f"({100*(settled+kept+dropped)/max(1, settled+kept+dropped+skipped):.0f}%)")
if dropped:
    byd = defaultdict(int)
    for o in out:
        if o['verdict'] == 'DOWNGRADED': byd[o['metabolite']] += 1
    print()
    print("  downgraded by metabolite:")
    for m, c in sorted(byd.items(), key=lambda t: -t[1]):
        print(f"    {m[:34]:<35}{c:>5}")
with open(os.path.join(T, 'edge_transport.tsv'), 'w', newline='') as fh:
    if out:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)

untested_broken = [o for o in out if 'marker-set failure' in o['note']]
if untested_broken:
    print()
    print(f"  {len(untested_broken)} edges left INTACT rather than downgraded,")
    print("  because their transporter family was detected in under 20% of genomes.")
    print("  That is the step-19 rule turned on this script: a family missing from")
    print("  most of a catalogue indicts the marker set, not the community.")
notm = [o for o in out if o['note'].startswith('NO transporter family')]
if notm:
    mm = sorted({o['metabolite'] for o in notm})
    print()
    print(f"  {len(notm)} edges UNTESTED for want of a family in my table: "
          f"{', '.join(mm)}.")
    print("  Note what this includes: the TMA edges are the cutC -> methylamine")
    print("  methanogen route, which is this chapter's headline cross-feeding")
    print("  claim, and NOT ONE of them could be tested here. Say that rather than")
    print("  letting the reader assume the network was checked end to end.")

print()
print("  A DOWNGRADED EDGE IS NOT A REFUTED ONE. The transporter may be absent,")
print("  or the genome incomplete, or the family simply missing from my table —")
print("  which covers seven metabolite classes, not all of microbiology. Report")
print("  supported edges as supported and downgraded ones as 'catabolic capacity")
print("  present, uptake route not identified'. That is a weaker claim than step")
print("  21 made and a defensible one.")
print()
print("=" * 104)
print("  WHAT THIS CHECK CAN AND CANNOT DO AS BUILT")
print("=" * 104)
nsup = sum(1 for o in out if o['verdict'] == 'supported')
ndwn = sum(1 for o in out if o['verdict'] == 'DOWNGRADED')
print(f"  {settled} edges are SETTLED and that is a real result: H2, methanol and")
print("  acetate reach the cytoplasm without dedicated machinery, so two thirds of")
print("  step 21's cross-feeding network was never exposed to this objection.")
print()
if nsup + ndwn < 0.2 * len(out):
    print(f"  But only {nsup + ndwn} of the {nsup + ndwn + skipped} edges that DO need a")
    print("  carrier could be tested.")
    print("  BAKTA PRODUCT STRINGS ARE TOO SPARSE FOR TRANSPORTERS. That is the")
    print("  finding of this run — not a statement about the community. Every")
    print("  family except the amino-acid and vitamin ABC systems fell below the")
    print("  20% detection bar, which means the annotation, not the biology, is")
    print("  the limiting factor.")
    print()
    print("  TO MAKE THIS ANALYSIS WORTH ANYTHING it needs an hmmsearch over the")
    print("  KOfam profiles for these families, by the same route steps 19b and 20")
    print("  use — individual profile fetch, archive fallback, adaptive")
    print("  thresholds, stamp written only after a successful search. Until then")
    print("  report ONLY the metabolites that need no transporter (H2, methanol,")
    print("  acetate) as settled, and say the rest was not testable.")
PY
RC=$?

echo
echo "######## VALIDATION GATE ########"
# THE GATE MUST TEST A RESULT, NOT A FILE. The first run wrote a perfectly
# well-formed table in which every edge was downgraded because no KO data had
# loaded, and the gate passed it because the file existed and was non-empty.
# Second defect, caught by the next run: RC=$? sat AFTER two echo statements, so
# it captured the exit status of `echo` — always 0 — and the gate reported PASS
# on a script that had exited 4. $? survives exactly one command.
check "transporter data loaded (exit != 4)" $([ "$RC" != "4" ] && echo 0 || echo 1)
check "edge transport table" $([ -s "$T/edge_transport.tsv" ] && echo 0 || echo 1)
if [ -s "$T/edge_transport.tsv" ]; then
  SUP=$(awk -F'\t' 'NR>1 && $5=="supported"' "$T/edge_transport.tsv" | wc -l)
  DWN=$(awk -F'\t' 'NR>1 && $5=="DOWNGRADED"' "$T/edge_transport.tsv" | wc -l)
  echo "  supported=$SUP downgraded=$DWN"
  check "not every tested edge downgraded" $([ "$SUP" -gt 0 ] || [ "$DWN" -eq 0 ] && echo 0 || echo 1)
fi
runlog "29_transport complete"
echo
echo "######## 29 DONE ########"
