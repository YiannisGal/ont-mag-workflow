#!/usr/bin/env bash
# =============================================================================
# 23d_t9ss_export.sh  --  does Gaoshiqia have a route OUT of the cell at all?
#
# WHERE 23c LEFT IT
#   161 signal-peptide CAZymes, 137 PUL-associated, 86 both, 75 with a signal
#   peptide but outside any detected PUL. Those 75 are the whole question. A
#   signal peptide gets a protein across the INNER membrane into the periplasm
#   and no further. To reach the medium a Bacteroidota protein needs a second,
#   independent machine.
#
#   That machine has a name: the TYPE IX SECRETION SYSTEM (T9SS, also called
#   the Por secretion system). It is a Bacteroidetes invention. Its core is
#   GldK/GldL/GldM/GldN plus SprA(SovA), PorU, PorV, and cargo proteins carry a
#   conserved C-TERMINAL DOMAIN that PorU cleaves off during export.
#
#   So the question becomes falsifiable:
#     NO T9SS      -> the 75 have no route across the outer membrane by this
#                     mechanism. They are periplasmic. 23b stands unqualified
#                     and the public-good reading is dead.
#     T9SS present -> export is POSSIBLE. Not demonstrated — see the warning in
#                     section 2 — but possible, and worth writing as such.
#
# SECTION 3 ATTACKS THE SAME QUESTION FROM THE OTHER END
#   Trafficking prediction is indirect. There is a direct ecological test that
#   needs no trafficking at all: if Gaoshiqia releases oligosaccharides, then
#   SOMEBODY has to eat them. A genome that can finish oligosaccharides but
#   cannot open polymers is a genome that depends on someone else's enzymes.
#   Counting endo-acting against exo-acting CAZymes per genome finds those
#   dependants. If the Bacillota fermenters carry exo-glycosidases and central
#   sugar catabolism but almost no endo-acting enzymes, the community has a
#   demand for a public good regardless of what DeepSig thinks.
#
# WHAT THIS SCRIPT DOES NOT DO
#   It does not run an HMM search for the T9SS C-terminal domain, because the
#   TIGRFAM models for it (TIGR04183 type A, TIGR04131 type B) are not part of
#   any database this pipeline already downloads, so it is not fetched
#   here. Cargo counts here come from Bakta's product
#   strings and are therefore a FLOOR, clearly labelled as one. The MACHINERY
#   result does not depend on that limitation and is the load-bearing one.
#
# CONTROL
#   T9SS is confined to Bacteroidetes. If a Bacillota or an archaeon scores a
#   full machinery complement, the detector is matching noise and the section
#   is void. That check runs automatically and will say so.
#
# INPUTS   09_annotation/<mag>/<mag>.tsv   (Bakta)
#          23_extras/deepsig/*.tsv          (23)
#          23_extras/pul/pul_loci.tsv       (23b)
#          07_mags/MAG_TAXONOMY.tsv
# OUTPUT   23_extras/t9ss/t9ss_machinery.tsv
#          23_extras/t9ss/cazyme_endo_exo.tsv
#          23_extras/t9ss/gaoshiqia_nonpul_cazymes.tsv
#
# RUNTIME  seconds. No new tools, no downloads.
# USAGE    ws2 && bash scripts/23d_t9ss_export.sh 2>&1 | tee logs/23d.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="23d_t9ss_export"
A="$WS2/09_annotation"; F="$WS2/07_mags"; E="$WS2/23_extras"; T="$E/t9ss"
mkdir -p "$T"

echo "######## 23d — T9SS AND THE DEMAND SIDE ########"
echo "  started: $(date -Is)"
runlog "23d_t9ss_export started"
[ -d "$A" ] || { echo "  !! 09_annotation not found — run Bakta first"; exit 1; }

python3 - "$A" "$F" "$E" "$T" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter

A, F, E, T = sys.argv[1:5]

# --------------------------------------------------------------------------
# taxonomy
# --------------------------------------------------------------------------
tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus', '') or '-')
def phy(m): return ((tax.get(m, {}) or {}).get('phylum', '') or '')
def dom(m): return ((tax.get(m, {}) or {}).get('domain', '') or '')
def is_bacteroidota(m): return 'bacteroidota' in phy(m).lower()

# --------------------------------------------------------------------------
# Bakta tables: locus -> (gene, product)
# --------------------------------------------------------------------------
gene, prod = defaultdict(dict), defaultdict(dict)
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].lower() == 'cds':
            gene[n][f[5]] = f[6]
            prod[n][f[5]] = f[7]

if not prod:
    print("  !! no Bakta tables parsed — nothing to do"); sys.exit(1)

# ==========================================================================
# 1. T9SS MACHINERY
# ==========================================================================
# Core translocon and the shared gliding-motility motor. GldK/L/M/N power both
# T9SS and gliding in Bacteroidetes and are the most reliably annotated of the
# set; SprA/SovA is the outer-membrane translocon itself; PorU is the sortase
# that cleaves the cargo C-terminal domain and PorV the shuttle that delivers
# cargo to it. PorU + PorV together are the attachment branch: cargo they
# process ends up ANCHORED TO THE SURFACE, not free in the medium.
CORE   = {'gldk','gldl','gldm','gldn','spra','sova'}
# 'sprt' REMOVED after the first run. It fired in Lysinibacillus and
# Lacrimispora — both Bacillota, neither of which can have a T9SS. The name
# collides: in most bacteria `sprT` is the SprT-like DNA-repair protein and has
# nothing to do with secretion. The control caught it only because sprT is
# accessory-only and never reached a core complement; had it been in CORE the
# whole section would have been quietly wrong.
ACCESS = {'poru','porv','porp','porq','porz','port','spre','sprf','sprb'}
PRODPAT = ('type ix secretion', 't9ss', 'por secretion',
           'gliding motility protein gld', 'gliding motility-associated')

mach = {}
for m in sorted(prod, key=mid):
    found_core, found_acc = set(), set()
    for lt, p in prod[m].items():
        g = (gene[m].get(lt) or '').strip().lower()
        pl = (p or '').lower()
        g2 = re.sub(r'[^a-z]', '', g)
        if g2 in CORE:   found_core.add(g2)
        if g2 in ACCESS: found_acc.add(g2)
        if not g2 and any(x in pl for x in PRODPAT):
            mm = re.search(r'\b(gld[klmn]|spr[aetbf]|por[uvpqzt]|sova)\b', pl)
            if mm:
                w = mm.group(1)
                (found_core if w in CORE else found_acc).add(w)
    mach[m] = (found_core, found_acc)

print()
print("=" * 104)
print("  1. TYPE IX SECRETION SYSTEM — IS THERE A ROUTE ACROSS THE OUTER MEMBRANE?")
print("=" * 104)
print("  T9SS is a Bacteroidetes-only machine. Core = GldK/L/M/N + SprA(SovA).")
print("  PorU/PorV are the attachment branch: what they process is ANCHORED to")
print("  the cell surface. Surface-attached is still not released into porewater.")
print()
print(f"  {'MAG ID':<17}{'genus':<20}{'phylum':<20}{'core':>6}{'access':>8}  genes found")
rows = []
for m in sorted(prod, key=mid):
    c, a = mach[m]
    if not c and not a: continue
    print(f"  {mid(m):<17}{gen(m)[:19]:<20}{phy(m)[:19]:<20}{len(c):>6}{len(a):>8}"
          f"  {','.join(sorted(c | a))}")
    rows.append(dict(mag_id=mid(m), genus=gen(m), phylum=phy(m),
                     n_core=len(c), n_accessory=len(a),
                     core_genes=','.join(sorted(c)),
                     accessory_genes=','.join(sorted(a))))
if not rows:
    print("  (no T9SS genes detected in any genome)")
with open(os.path.join(T, 't9ss_machinery.tsv'), 'w', newline='') as fh:
    if rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

# ---- CONTROL: T9SS outside Bacteroidota means the detector is matching noise
bad = [m for m in prod if not is_bacteroidota(m) and len(mach[m][0]) >= 3]
print()
if bad:
    print("  >> CONTROL FAILED. These non-Bacteroidota scored >=3 core genes:")
    for m in bad:
        print(f"     {mid(m)}  {gen(m)}  {phy(m)}  {','.join(sorted(mach[m][0]))}")
    print("     T9SS does not occur outside Bacteroidetes. Either the gene names")
    print("     are colliding with something else or Bakta is over-calling. DO")
    print("     NOT use section 1 until this is understood.")
    T9SS_OK = False
else:
    nb = sum(1 for m in prod if not is_bacteroidota(m))
    print(f"  >> CONTROL PASSED. {nb} non-Bacteroidota genomes, none with a core")
    print("     complement. The detector is behaving like a Bacteroidetes-specific")
    print("     marker, which is what T9SS is.")
    T9SS_OK = True

# ==========================================================================
# 2. CARGO — A FLOOR, NOT A COUNT
# ==========================================================================
CARGO_PAT = ('t9ss type a sorting domain', 't9ss type b', 't9ss c-terminal',
             'sorting domain-containing', 'por secretion system c-terminal',
             'secretion system c-terminal sorting domain')
cargo = defaultdict(set)
for m in prod:
    for lt, p in prod[m].items():
        if any(x in (p or '').lower() for x in CARGO_PAT): cargo[m].add(lt)

print()
print("=" * 104)
print("  2. T9SS CARGO — PROTEINS ANNOTATED WITH A C-TERMINAL SORTING DOMAIN")
print("=" * 104)
print("  READ THIS AS A FLOOR. These come from Bakta product strings, not from an")
print("  HMM search against TIGR04183/TIGR04131. Any cargo protein whose product")
print("  string does not happen to mention the domain is invisible here. A zero")
print("  in this section is NOT evidence of absence.")
print()
print(f"  {'MAG ID':<17}{'genus':<20}{'cargo (floor)':>14}")
for m in sorted(prod, key=mid):
    if cargo.get(m): print(f"  {mid(m):<17}{gen(m)[:19]:<20}{len(cargo[m]):>14}")
if not cargo: print("  (none annotated — see the floor warning above)")

# ==========================================================================
# 3. THE DEMAND SIDE — WHO COULD NOT FEED THEMSELVES?
# ==========================================================================
# Endo-acting enzymes attack a polymer in the middle and open it up. Exo-acting
# enzymes remove terminal sugars from something already opened. A genome rich in
# exo but poor in endo can finish a meal it cannot start.
ENDO = ('cellulase','xylanase','chitinase','pectate lyase','alginate lyase',
        'agarase','glucanase','polysaccharide lyase','alpha-amylase',
        'endoglucanase','endo-','laminarinase','pullulanase','mannanase',
        'arabinanase','galactanase','amylase','chondroitin')
EXO  = ('galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
        'arabinofuranosidase','glucuronidase','hexosaminidase','sialidase',
        'neuraminidase','exo-','glucoamylase','trehalase','levanase',
        'carbohydrate esterase')

print()
print("=" * 104)
print("  3. ENDO vs EXO CAZymes — WHO OPENS POLYMERS, WHO ONLY FINISHES THEM")
print("=" * 104)
print("  Coarse, from product strings. 'Amylase' is counted endo even though the")
print("  beta- and gluco- forms are exo, so the endo column is if anything")
print("  generous. The ratio is the signal, not the absolute numbers.")
print()
print("  COMPLETENESS GATES ABSENCE, NOT PRESENCE. The first version of this")
print("  table applied the >=90% rule to every role and thereby refused to call")
print("  the four polymer degraders, three of which sit at 75-88%. That was")
print("  wrong in one direction: if an 86%-complete genome CARRIES 15 endo-acting")
print("  enzymes, those 15 genes exist — incompleteness can only have made the")
print("  count too LOW. So:")
print("    'POLYMER DEGRADER' is a PRESENCE claim  -> called at any completeness")
print("    'dependant / limited' are ABSENCE claims -> only called at >=90%")
print("  'Carbohydrate esterase' is counted exo though it strips decorations")
print("  rather than terminal sugars — it inflates the exo column slightly.")
print()
print(f"  {'MAG ID':<17}{'genus':<20}{'phylum':<18}{'compl':>7}{'endo':>6}{'exo':>6}"
      f"{'exo/endo':>10}  role")
erows = []
def _compl(m):
    try: return float((tax.get(m, {}) or {}).get('completeness', '') or 0)
    except ValueError: return 0.0
for m in sorted(prod, key=mid):
    ps = [(p or '').lower() for p in prod[m].values()]
    ne = sum(1 for p in ps if any(x in p for x in ENDO))
    nx = sum(1 for p in ps if any(x in p for x in EXO) and
                              not any(y in p for y in ENDO))
    if ne == 0 and nx == 0: continue
    ratio = (nx / ne) if ne else float('inf')
    cp = _compl(m)
    if ne >= 10:                 role = 'POLYMER DEGRADER'        # presence
    elif cp < 90:                role = f'absence not callable ({cp:.1f}%)'
    elif ne <= 2 and nx >= 5:    role = 'dependant — finishes, cannot start'
    elif ne == 0 and nx > 0:     role = 'dependant — no endo capacity at all'
    else:                        role = 'limited'
    rs = '  inf' if ratio == float('inf') else f"{ratio:>6.1f}"
    print(f"  {mid(m):<17}{gen(m)[:19]:<20}{phy(m)[:17]:<18}{cp:>6.1f}%{ne:>6}{nx:>6}"
          f"{rs:>10}  {role}")
    erows.append(dict(mag_id=mid(m), genus=gen(m), phylum=phy(m),
                      completeness=f"{cp:.1f}", n_endo=ne, n_exo=nx,
                      exo_per_endo=('inf' if ne == 0 else f"{ratio:.2f}"),
                      role=role))
with open(os.path.join(T, 'cazyme_endo_exo.tsv'), 'w', newline='') as fh:
    if erows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(erows[0]))
        w.writeheader(); w.writerows(erows)

ndep = sum(1 for r in erows if r['role'].startswith('dependant'))
ndeg = sum(1 for r in erows if r['role'] == 'POLYMER DEGRADER')

# ==========================================================================
# 4. THE 75 — WHAT ARE THEY, ACTUALLY?
# ==========================================================================
CAZY = ('glycoside hydrolase','glycosyl hydrolase','polysaccharide lyase',
        'carbohydrate esterase','glucanase','xylanase','chitinase','amylase',
        'pectinesterase','pectate lyase','cellulase','agarase','alginate lyase',
        'galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
        'arabinofuranosidase','levanase','pullulanase','chondroitin')
pul_loci = defaultdict(set)
id2wn = {mid(w): w for w in tax}
pf = os.path.join(E, 'pul', 'pul_loci.tsv')
if os.path.exists(pf):
    for r in csv.DictReader(open(pf), delimiter='\t'):
        wn = id2wn.get(r.get('mag_id', ''), r.get('mag_id', ''))
        for lt in (r.get('cazymes') or '').split(';'):
            if lt.strip(): pul_loci[wn].add(lt.strip())

gao = next((m for m in prod if 'gaoshiqia' in gen(m).lower()), None)
n75 = 0
if gao:
    sec = set()
    dp = os.path.join(E, 'deepsig', f'{gao}.tsv')
    if os.path.exists(dp):
        for line in open(dp):
            if line.startswith('#') or not line.strip(): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3 or 'signal' not in ' '.join(f[1:]).lower(): continue
            if any(x in (prod[gao].get(f[0]) or '').lower() for x in CAZY):
                sec.add(f[0])
    out = sorted(sec - pul_loci.get(gao, set()))
    n75 = len(out)
    print()
    print("=" * 104)
    print(f"  4. THE {n75} SIGNAL-PEPTIDE CAZymes OF {mid(gao)} THAT SIT OUTSIDE ANY PUL")
    print("=" * 104)
    print("  Written to gaoshiqia_nonpul_cazymes.tsv. The families here decide how")
    print("  much the number is worth: a set dominated by exo-glycosidases is a")
    print("  periplasmic finishing crew and says nothing about export, while")
    print("  endo-acting enzymes on large polymers would have to work outside the")
    print("  cell because the substrate cannot cross the outer membrane.")
    print()
    fam = Counter()
    with open(os.path.join(T, 'gaoshiqia_nonpul_cazymes.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['mag_id', 'locus_tag', 'gene', 'product', 'endo_or_exo',
                    't9ss_cargo_annotated'])
        for lt in out:
            p = prod[gao].get(lt, '')
            pl = p.lower()
            ee = 'endo' if any(x in pl for x in ENDO) else (
                 'exo' if any(x in pl for x in EXO) else 'unclear')
            fam[ee] += 1
            w.writerow([mid(gao), lt, gene[gao].get(lt, ''), p, ee,
                        'yes' if lt in cargo.get(gao, set()) else 'no'])
    tot75 = sum(fam.values()) or 1
    for k in ('endo', 'exo', 'unclear'):
        print(f"    {k:<9}{fam.get(k,0):>4}  ({100*fam.get(k,0)/tot75:.0f}%)")
    n_endo75 = fam.get('endo', 0)
    endo_frac = n_endo75 / tot75
    n_cargo75 = len([lt for lt in out if lt in cargo.get(gao, set())])
    # DO NOT OVER-READ A ZERO HERE. If only a handful of proteins in the whole
    # genome carry an annotated sorting domain, then zero among 75 is what
    # chance predicts and says nothing either way. State the expectation.
    ncargo_tot, nprot = len(cargo.get(gao, set())), len(prod[gao]) or 1
    exp = tot75 * ncargo_tot / nprot
    print()
    print(f"    annotated T9SS cargo among these {tot75}: {n_cargo75}")
    print(f"    expected by chance: {exp:.2f}  ({ncargo_tot} annotated cargo in "
          f"{nprot} proteins)")
    if exp < 1:
        print("    -> the expectation is below one, so this count CANNOT")
        print("       discriminate. Ignore it until an HMM search replaces it.")
else:
    n_endo75 = n_cargo75 = 0; endo_frac = 0.0

# ==========================================================================
# VERDICT
# ==========================================================================
print()
print("=" * 104)
print("  VERDICT")
print("=" * 104)
if gao:
    c, a = mach[gao]
    print(f"  {mid(gao)} ({gen(gao)}): T9SS core {len(c)}/5 ({','.join(sorted(c)) or 'none'}), "
          f"accessory {len(a)} ({','.join(sorted(a)) or 'none'}).")
    print(f"  Of the {n75} non-PUL signal-peptide CAZymes, {n_endo75} "
          f"{'is' if n_endo75 == 1 else 'are'} endo-acting and {n_cargo75} "
          f"{'carries' if n_cargo75 == 1 else 'carry'} an annotated T9SS sorting domain.")
    print()
    if not T9SS_OK:
        print("  >> NO VERDICT. The T9SS control failed; section 1 cannot be used.")
    elif len(c) == 0:
        print("  >> NO EXPORT ROUTE FOUND. Without T9SS core genes there is no")
        print("     Bacteroidetes mechanism to move these proteins across the outer")
        print("     membrane. The 75 are periplasmic until shown otherwise, and the")
        print("     public-good reading of step 20 does not recover. Write Gaoshiqia")
        print("     as a private digester and say so plainly.")
    elif len(c) >= 3 and endo_frac >= 0.25:
        # The first version of this branch fired on n_endo75 > 0, which is not a
        # test — one endo-acting enzyme in a hundred would have passed it. A
        # quarter of the set is still an arbitrary line, but it is at least a
        # line about proportion rather than mere existence.
        print("  >> EXPORT IS POSSIBLE AND THERE IS SOMETHING WORTH EXPORTING.")
        print("     The machinery is present and some of the non-PUL enzymes act on")
        print("     polymers too large to cross the outer membrane, so they would")
        print("     have to work outside to work at all. This is the strongest")
        print("     version of the public-good case the genome can support — and it")
        print("     is still POSSIBILITY, not demonstration. PorU/PorV-processed")
        print("     cargo ends up anchored to the cell surface; surface-attached")
        print("     digestion feeds the cell that owns the surface first.")
    elif len(c) >= 3:
        print("  >> MACHINERY PRESENT, CARGO UNCONVINCING. T9SS exists, but only")
        print(f"     {100*endo_frac:.0f}% of the non-PUL enzymes act on polymers. The rest are")
        print("     exo-acting or unclassifiable — a periplasmic finishing crew")
        print("     working on oligosaccharides the PULs already imported. That is")
        print("     the private-digester picture with extra steps. 23b stands.")
        print("     Having the machinery is not the same as having freight for it.")
    else:
        print("  >> PARTIAL MACHINERY. Too few core genes to call a working T9SS,")
        print("     and a partial complement in a 99.9% complete genome is more")
        print("     likely annotation drift than a real half-system. Treat as")
        print("     no route until an HMM search says otherwise.")
print()
nlow = sum(1 for r in erows if r['role'].startswith('absence not callable'))
print(f"  DEMAND SIDE: {ndeg} genome{'' if ndeg == 1 else 's'} carr"
      f"{'ies' if ndeg == 1 else 'y'} enough endo-acting enzymes to open polymers; "
      f"{ndep} look{'s' if ndep == 1 else ''}")
print("  like a dependant that can finish oligosaccharides but not open them")
print(f"  ({nlow} genome{'' if nlow == 1 else 's'} too incomplete for an absence call).")
print("  This holds whatever the trafficking prediction says. If the dependants")
print("  are the Bacillota fermenters that feed the methanogens, the community")
print("  needs a polymer opener whether or not Gaoshiqia is a willing one —")
print("  and note that a private digester still leaks. Lysis and imperfect")
print("  capture release oligosaccharides without anything being secreted.")
print()
print("  WHAT WOULD SETTLE IT PROPERLY: hmmsearch with TIGR04183 and TIGR04131")
print("  against the three Bacteroidota proteomes. That needs those two models")
print("  fetched onto the server first.")
PY

echo
echo "######## VALIDATION GATE ########"
check "T9SS machinery table" $([ -s "$T/t9ss_machinery.tsv" ] && echo 0 || echo 1)
check "endo/exo table"       $([ -s "$T/cazyme_endo_exo.tsv" ] && echo 0 || echo 1)
runlog "23d_t9ss_export complete"
echo
echo "######## 23d DONE ########"
