#!/usr/bin/env bash
# =============================================================================
# 23e_secretion_systems.sh  --  the census 23d should have been
#
# WHY THIS EXISTS
#   23d looked at ONE secretion system because it was built to answer one
#   question: can a Gaoshiqia CAZyme cross the outer membrane. T9SS was the
#   right tool for that question and the wrong scope for the chapter. A
#   community's export machinery is a trait in its own right, and there are ten
#   other ways out of a cell that were never looked at.
#
# WHAT IS ADDED, AND WHY EACH ONE EARNS ITS PLACE HERE
#
#   Sec        Universal, essential. Present in all 29 or the detector is
#              broken — this is the POSITIVE CONTROL for the whole script.
#
#   Tat        Exports FOLDED proteins. This matters more here than anywhere
#              else in the analysis: hydrogenases, formate dehydrogenases and
#              the anaerobic reductases bind their cofactors in the cytoplasm
#              and CANNOT go through Sec. Step 21 built 184 hydrogen-exchange
#              pairs without ever asking whether the hydrogenases are
#              periplasmic. Tat is how you ask.
#
#   T6SS       An antibacterial weapon — a needle that injects toxins into a
#              neighbouring cell. Step 21's competition layer is inferred from
#              resource overlap, which is potential competition. A T6SS is
#              actual antagonism. If the competing pairs also carry T6SS the
#              competition edges stop being an inference.
#
#   T4SS       Conjugation. With 68 proviruses across 18 MAGs this community is
#              already established as mobile-element rich; T4SS is the other
#              half of horizontal transfer.
#
#   T1/T2/T5   The remaining Gram-negative routes across the outer membrane.
#              T2SS is the one that would give Gaoshiqia's 75 CAZymes a second
#              exit if T9SS did not — so its ABSENCE is part of 23d's argument
#              and was never checked. Checking it now.
#
#   T3SS       Expected ABSENT. An injectisome in an anoxic lagoon sediment
#              would mean either a eukaryote-associated organism or a false
#              positive, and either is worth knowing. NEGATIVE CONTROL.
#
#   T7SS       The Gram-positive equivalent, for the Bacillota that make up
#              most of this catalogue.
#
#   Sortase    Gram-positive surface anchoring. Not secretion but the
#              destination question again: an LPXTG protein is covalently
#              stitched to peptidoglycan, so it is surface-displayed, not
#              released. Same distinction that decided the Gaoshiqia question.
#
# THE MOTIF SCANS, AND WHY THEY ARE THE HONEST PART
#   Machinery tells you the cell HAS a route. It does not tell you how much
#   traffic uses it. Two routes leave a signature in the protein sequence:
#     - Tat substrates carry a TWIN-ARGININE motif (S/T-R-R-x-Φ-Φ) in the
#       signal peptide.
#     - Sortase substrates carry LPXTG near the C-terminus.
#   Both are scanned directly from the Bakta .faa files, and each is
#   CROSS-VALIDATED against its machinery: genomes with tatC should carry many
#   twin-arginine proteins and genomes without it should carry few. If that
#   correlation fails, the motif scan is matching noise and says so.
#
# LIMITS, STATED UP FRONT
#   Gene-name and product-string detection, not HMM profiles. It will miss
#   diverged systems and it can be fooled by name collisions — 23d found
#   exactly that with `sprT`, which is a DNA-repair protein in Bacillota and
#   has nothing to do with secretion. Every system here therefore carries an
#   EXPECTATION, and a system found where it should not be is reported as a
#   detector failure rather than a discovery.
#
# INPUTS   09_annotation/<mag>/<mag>.tsv and <mag>.faa   ·   07_mags/MAG_TAXONOMY.tsv
# OUTPUT   23_extras/secretion/secretion_systems.tsv
#          23_extras/secretion/motif_substrates.tsv
#          23_extras/secretion/t6ss_detail.tsv
# RUNTIME  under a minute. No new tools, no downloads.
# USAGE    ws2 && bash scripts/23e_secretion_systems.sh 2>&1 | tee logs/23e.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="23e_secretion"
A="$WS2/09_annotation"; F="$WS2/07_mags"; S="$WS2/23_extras/secretion"
mkdir -p "$S"

echo "######## 23e — SECRETION SYSTEM CENSUS ########"
echo "  started: $(date -Is)"
runlog "23e_secretion started"
[ -d "$A" ] || { echo "  !! 09_annotation not found"; exit 1; }

python3 - "$A" "$F" "$S" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict

A, F, S = sys.argv[1:4]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus', '') or '-')
def phy(m): return ((tax.get(m, {}) or {}).get('phylum', '') or '')
def dom(m): return ((tax.get(m, {}) or {}).get('domain', '') or '')
def is_arc(m): return 'archae' in dom(m).lower()
def gramneg(m):
    return any(x in phy(m).lower() for x in
               ('bacteroidota','pseudomonadota','desulfobacterota',
                'verrucomicrobiota','planctomycetota','spirochaetota'))

# --------------------------------------------------------------------------
# Bakta tables. NOTE ON GENE-NAME NORMALISATION, which 23d got away with:
# stripping every non-letter turns virB4 and virB11 into the same token 'virb'.
# 23d's markers had no digits so it did not matter there. Here they do, so the
# trailing paralog suffix Bakta appends (_1, _2) is removed FIRST and digits
# are then kept.
# --------------------------------------------------------------------------
def norm(g):
    g = (g or '').strip().lower()
    g = re.sub(r'_\d+$', '', g)
    return re.sub(r'[^a-z0-9]', '', g)

gene, prod = defaultdict(dict), defaultdict(dict)
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].lower() == 'cds':
            gene[n][f[5]] = f[6]; prod[n][f[5]] = f[7]
if not prod:
    print("  !! no Bakta tables parsed"); sys.exit(1)
MAGS = sorted(prod, key=mid)

# --------------------------------------------------------------------------
# systems: (label, gene tokens, product substrings, min distinct hits, expectation)
# --------------------------------------------------------------------------
# THE RULE THAT FIXES THE FIRST RUN'S FALSE POSITIVES:
#   A SYSTEM THAT CROSSES THE OUTER MEMBRANE MUST HAVE ITS OUTER-MEMBRANE
#   COMPONENT. A Gram-positive has no outer membrane, so it cannot have a
#   secretin or a TolC channel — and if the cytoplasmic ATPase is found without
#   one, what has been found is a HOMOLOGUE OF SOMETHING ELSE.
# The first run called T2SS in ten Bacillota (competence and type IV pilus genes
# are genuinely homologous to T2SS and Bakta names them gspE/gspF), and T3SS in
# eight (the FLAGELLAR export apparatus is the ancestor of the injectisome and
# shares the Sct nomenclature — SctN=FliI, SctV=FlhA, SctU=FlhB). Neither was a
# bug in the annotation. Both were my marker sets ignoring the architecture.
# The flagellar-homologous tokens are removed outright and each Gram-negative
# system now carries a REQUIRED component in the sixth field.
SYS = [
 ('Sec',      {'seca','secy','sece','secg','secd','secf','yidc','yajc','ffh','ftsy',
               'sec61','srp54','srpr'},
              ('preprotein translocase','protein translocase subunit',
               'protein-export membrane protein'), 3, 'BACTERIA', set()),
 ('Tat',      {'tata','tatb','tatc','tate','tatd'},
              ('twin-arginine translocase','sec-independent protein translocase'),
              1, '-', set()),
 ('T1SS',     {'tolc','hlyb','hlyd','prtd','prte','prtf','cyab','lapb'},
              ('type i secretion',), 2, 'GRAMNEG', {'tolc'}),
 ('T2SS',     {'gspd','gspe','gspf','gspg','gsph','gspi','gspj','gspk','gspl','gspm',
               'puld','pule','xcpq','epsd','epse'},
              ('type ii secretion','general secretion pathway'), 3, 'GRAMNEG',
              {'gspd','puld','xcpq','epsd'}),                      # the secretin
 # sctN/sctU/sctV/hrcN/hrcV/yscN/yscV/invA/spaO REMOVED — all flagellar homologues.
 ('T3SS',     {'sctc','sctj','sctd','sctw','yscc','hrcc','invg','prgh'},
              ('type iii secretion',), 2, 'NONE', {'sctc','yscc','hrcc','invg'}),
 ('T4SS',     {'virb2','virb3','virb4','virb5','virb6','virb8','virb9','virb10','virb11',
               'vird4','trab','trag','traj','tram','trau','trbl'},
              ('type iv secretion','conjugal transfer','conjugative transfer'),
              3, '-', {'virb4','virb10','virb11','vird4','trab'}),
 ('T5SS',     set(), ('autotransporter','type v secretion'), 1, 'GRAMNEG', set()),
 # hcp and vgrG are homologues of the phage tail tube and spike. In a catalogue
 # where 18 of 29 genomes carry a provirus they are the expected false positive,
 # so the membrane complex is now mandatory.
 ('T6SS',     {'tssa','tssb','tssc','tssf','tssg','tssh','tssj','tssk','tssl','tssm',
               'hcp','vgrg','clpv','paar','impa','impb','impc','evpa','evpb','evpc'},
              ('type vi secretion',), 3, 'GRAMNEG', {'tssm','tssj','tssl','tssb','tssc'}),
 ('T7SS/ESX', {'esxa','esxb','eccb','eccc','eccd','ecce','essa','essb','essc','esaa','eslb'},
              ('type vii secretion','esx-1 secretion'), 2, 'GRAMPOS', set()),
 ('T9SS',     {'gldk','gldl','gldm','gldn','spra','sova','poru','porv','porp','porq',
               'porz','port'},
              ('type ix secretion','por secretion'), 3, 'BACTEROIDOTA',
              {'gldl','gldm','gldn','spra','sova'}),
 ('Sortase',  {'srta','srtb','srtc','srtd','srte'}, ('sortase',), 1, 'GRAMPOS', set()),
]

# Product strings that must NOT count towards a system even though they match
# its keyword. Found by the smoke test: PorU, the T9SS cargo protease, is
# annotated "type IX secretion system sortase PorU" — so every T9SS-bearing
# Bacteroidota was being called Gram-positive-sortase-positive off one protein.
# A keyword is a substring, not a concept, and "sortase" appears inside the
# name of a machine that has nothing to do with LPXTG anchoring.
EXCLUDE = {
 'Sortase': ('type ix', 'poru', 'por secretion'),
 'T6SS':    ('phage tail', 'tail sheath', 'baseplate', 'tail tube'),
}

# NAMED GENES AND PRODUCT-STRING HITS ARE COUNTED SEPARATELY.
# The first run reported "37 T9SS markers" in Lentimicrobium. There are only
# twelve markers in the set — the other 25 were distinct PROTEINS whose product
# string mentioned T9SS, each generating its own token. Product hits therefore
# satisfied a "3 distinct markers" minimum trivially, and the column stopped
# meaning what its header said. Presence is now decided on named genes only;
# product hits are shown in brackets as supporting evidence.
gh = defaultdict(lambda: defaultdict(set))    # mag -> system -> {gene tokens}
ph = defaultdict(lambda: defaultdict(int))    # mag -> system -> product-only count
for m in MAGS:
    for lt in prod[m]:
        g = norm(gene[m].get(lt)); p = (prod[m].get(lt) or '').lower()
        for label, genes, pats, _mn, _exp, _rq in SYS:
            if any(x in p for x in EXCLUDE.get(label, ())): continue
            if g and g in genes: gh[m][label].add(g)
            elif pats and any(x in p for x in pats): ph[m][label] += 1

REQ = {s[0]: s[5] for s in SYS}
EXP = {s[0]: s[4] for s in SYS}

def forbidden(m, label):
    """Architecture the markers cannot override. A Gram-positive has no outer
    membrane, so it cannot have a TolC channel, a secretin, or an autotransporter
    barrel — whatever the gene names say. The first run reported six Bacillota as
    T1SS-positive on a `tolC` call; a Firmicute has nowhere to put a TolC, so
    those are efflux-pump homologues. Suppress the CALL, count the suppression,
    and print it. Suppressing quietly would be worse than the false positive."""
    e = EXP.get(label)
    if e == 'GRAMNEG':      return not gramneg(m)
    if e == 'GRAMPOS':      return gramneg(m)
    if e == 'BACTEROIDOTA': return 'bacteroidota' not in phy(m).lower()
    return False

suppressed = defaultdict(list)
def present(m, label, mn):
    h = gh[m].get(label, set())
    if len(h) < mn: return False
    r = REQ.get(label) or set()
    if r and not (h & r): return False
    if forbidden(m, label):
        if m not in suppressed[label]: suppressed[label].append(m)
        return False
    return True
hits = gh                                     # T6SS section reads named genes only

print()
print("=" * 112)
print("  1. SECRETION AND EXPORT MACHINERY ACROSS 29 GENOMES")
print("=" * 112)
print("  Detection is by gene name and product string, NOT by HMM — it misses")
print("  diverged systems and can be fooled by name collisions, which is why every")
print("  system carries an expectation and every outer-membrane system now")
print("  requires its outer-membrane component.")
print()
print("  Cell shows  namedGenes(productHits).  * = system called present:")
print("  enough named genes AND, for the outer-membrane systems, the required")
print("  component. ARCHAEA ARE BLANK, not zero — see the control below.")
print()
hdr = f"  {'MAG ID':<17}{'genus':<18}{'gram':<5}"
for label, _g, _p, mn, _e, _r in SYS: hdr += f"{label[:9]:>11}"
print(hdr)
rows = []
for m in MAGS:
    gr = 'arc' if is_arc(m) else ('neg' if gramneg(m) else 'pos')
    line = f"  {mid(m):<17}{gen(m)[:17]:<18}{gr:<5}"
    row = dict(mag_id=mid(m), genus=gen(m), phylum=phy(m), gram=gr)
    for label, _g, _p, mn, _e, _r in SYS:
        n = len(gh[m].get(label, ())); q = ph[m].get(label, 0)
        col = label.replace('/', '_')
        if gr == 'arc':
            # A zero that means "not measurable" is indistinguishable from a
            # zero that means "absent". Print neither.
            line += f"{'n/a':>11}"; row[col] = ''; row[col + '_present'] = 'n/a'
        else:
            cell = f"{n}({q})" if q else str(n)
            if present(m, label, mn): cell += '*'
            line += f"{cell:>11}"
            row[col] = n; row[col + '_prod'] = q
            row[col + '_present'] = 'yes' if present(m, label, mn) else 'no'
    print(line); rows.append(row)
with open(os.path.join(S, 'secretion_systems.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# --------------------------------------------------------------------------
# 2. CONTROLS
# --------------------------------------------------------------------------
print()
print("=" * 112)
print("  2. CONTROLS — DOES THE DETECTOR BEHAVE?")
print("=" * 112)
BACT = [m for m in MAGS if not is_arc(m)]
ARC  = [m for m in MAGS if is_arc(m)]
arc_sec = sum(len(gh[m].get('Sec', ())) for m in ARC)
print(f"  SCOPE: {len(BACT)} bacterial genomes assessed, {len(ARC)} archaeal EXCLUDED.")
print(f"  The first run applied the universal-Sec control to all 29 and failed on")
print(f"  the 8 archaea, which scored 0 for Sec and near-0 for everything else.")
print(f"  Archaea DO translocate protein — they have SecYE and SRP, but no SecA,")
print(f"  no SecG, and Bakta's databases name their components differently or not")
print(f"  at all (total archaeal Sec gene hits here: {arc_sec}).")
print(f"  APPLYING A BACTERIAL MARKER SET TO ARCHAEA AND REPORTING THE ZEROS IS")
print(f"  EXACTLY THE 19b rpoB/tuf ERROR. Archaeal rows are blank, not zero, and")
print(f"  archaeal secretion is UNASSESSED by this script rather than absent.")
print()
ok = True
for label, _g, _p, mn, exp, _r in SYS:
    pres = [m for m in BACT if present(m, label, mn)]
    if exp == 'BACTERIA':
        miss = [m for m in BACT if m not in pres]
        if miss:
            ok = False
            print(f"  [FAIL] {label}: essential and universal in Bacteria, but ABSENT in "
                  f"{len(miss)}: {', '.join(mid(x) for x in miss[:6])}")
            print(f"         Sec cannot be missing from a living bacterium. The gene-name")
            print(f"         detection is failing and NO ROW of section 1 is safe.")
        else:
            print(f"  [PASS] {label}: present in all {len(BACT)} bacterial genomes, "
                  f"as it must be.")
    elif exp == 'NONE':
        if pres:
            ok = False
            print(f"  [FLAG] {label}: expected absent, found in "
                  f"{', '.join(mid(x) for x in pres)} — inspect before believing.")
        else:
            print(f"  [PASS] {label}: absent everywhere, as expected for an anoxic "
                  f"sediment with no eukaryote-associated members.")
    elif exp in ('GRAMNEG', 'GRAMPOS', 'BACTEROIDOTA'):
        def wrong(m, _e=exp):
            if _e == 'GRAMNEG':      return not gramneg(m)
            if _e == 'GRAMPOS':      return gramneg(m)
            return 'bacteroidota' not in phy(m).lower()
        sup = suppressed.get(label, [])
        print(f"  [PASS] {label}: {len(pres)} genome(s), all within {exp}.")
        if sup:
            gt = 'Gram-positive' if exp == 'GRAMNEG' else (
                 'Gram-negative' if exp == 'GRAMPOS' else 'non-Bacteroidota')
            need = ('an outer membrane' if exp == 'GRAMNEG' else
                    'the right envelope' if exp == 'GRAMPOS' else 'the right lineage')
            print(f"         {len(sup)} SUPPRESSED — markers met the minimum in "
                  f"{', '.join(mid(x) for x in sup[:8])},")
            print(f"         but those genomes are {gt} and have no {need}, so the")
            print(f"         system cannot be there. The hits are homologues of")
            print(f"         something else — efflux pumps, competence machinery,")
            print(f"         flagellar export, phage tails.")

# --------------------------------------------------------------------------
# 3. MOTIF SCANS
# --------------------------------------------------------------------------
# Consensus is (S/T)-R-R-x-PHI-PHI. The first version excluded D/E at x, which
# is wrong — the variable position is genuinely variable and real Tat signals
# carry acidic residues there (the smoke-test sequence SRRDFLK, a perfectly
# ordinary motif, scored zero). Constrain the RR pair and the two hydrophobics,
# leave x alone.
TAT = re.compile(r'[ST]RR.[FGAVLIM][LVIMAFK]')
LPXTG = re.compile(r'LP[A-Z]TG')

def read_faa(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith('>'):
            if name: yield name, ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: yield name, ''.join(buf)

print()
print("=" * 112)
print("  3. SUBSTRATE MOTIFS — HOW MUCH TRAFFIC, NOT JUST WHICH ROUTES")
print("=" * 112)
print("  Twin-arginine (S/T-RR-x-PHI-PHI in the first 40 residues) marks Tat cargo:")
print("  proteins exported ALREADY FOLDED, which is the only way a cofactor-bound")
print("  enzyme can reach the periplasm. LPXTG near the C-terminus marks sortase")
print("  cargo: covalently stitched to peptidoglycan, so SURFACE-DISPLAYED rather")
print("  than released — the same destination distinction that decided 23d.")
print()
print(f"  {'MAG ID':<17}{'genus':<18}{'proteins':>9}{'twin-Arg':>10}{'tatC?':>7}"
      f"{'LPXTG':>8}{'srt?':>6}")
mrows = []
for m in MAGS:
    fa = os.path.join(A, m, f'{m}.faa')
    if not os.path.exists(fa): continue
    ntot = ntat = nlpx = 0
    for nm, sq in read_faa(fa):
        ntot += 1
        if TAT.search(sq[:40]): ntat += 1
        if len(sq) > 40 and LPXTG.search(sq[-60:]): nlpx += 1
    has_tat = 'yes' if 'tatc' in hits[m].get('Tat', ()) else 'no'
    has_srt = 'yes' if present(m, 'Sortase', 1) else 'no'
    print(f"  {mid(m):<17}{gen(m)[:17]:<18}{ntot:>9}{ntat:>10}{has_tat:>7}"
          f"{nlpx:>8}{has_srt:>6}")
    mrows.append(dict(mag_id=mid(m), genus=gen(m), n_proteins=ntot,
                      twin_arginine=ntat, has_tatC=has_tat,
                      lpxtg=nlpx, has_sortase=has_srt))
if mrows:
    with open(os.path.join(S, 'motif_substrates.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(mrows[0]))
        w.writeheader(); w.writerows(mrows)

# ---- cross-validation of the motif scan against its own machinery
def _rate(sel, key):
    g = [r for r in mrows if r[key] == 'yes'] if sel else [r for r in mrows if r[key] == 'no']
    if not g: return None
    return sum(r['twin_arginine' if key == 'has_tatC' else 'lpxtg'] for r in g) / len(g)

print()
for key, motif in (('has_tatC', 'twin-arginine'), ('has_sortase', 'LPXTG')):
    wi, wo = _rate(True, key), _rate(False, key)
    if wi is None or wo is None:
        print(f"  {motif}: only one group present — correlation not testable.")
    elif wi > wo * 1.5:
        print(f"  [PASS] {motif}: mean {wi:.1f} per genome WITH the machinery vs "
              f"{wo:.1f} without.")
        print(f"         The motif tracks the transporter, so it is finding real cargo.")
    else:
        print(f"  [FLAG] {motif}: mean {wi:.1f} with the machinery vs {wo:.1f} "
              f"without — no enrichment.")
        print(f"         A motif that appears just as often without its transporter is")
        print(f"         matching background sequence. Do not report these counts.")

# --------------------------------------------------------------------------
# 4. T6SS DETAIL — competition as mechanism rather than inference
# --------------------------------------------------------------------------
print()
print("=" * 112)
print("  4. T6SS — DOES ANYTHING HERE CARRY A WEAPON?")
print("=" * 112)
t6 = [m for m in MAGS if present(m, 'T6SS', 3)]
part = [m for m in MAGS if 0 < len(hits[m].get('T6SS', ())) < 3]
if t6:
    for m in t6:
        print(f"  {mid(m):<17}{gen(m)[:19]:<20}{','.join(sorted(hits[m]['T6SS']))}")
    print()
    print("  Step 21 inferred competition from resource-niche overlap, which is")
    print("  POTENTIAL competition. A T6SS is an actual apparatus for killing a")
    print("  neighbour. Cross these genomes against the competition edge list: a")
    print("  high-overlap pair in which one member is armed is a different claim")
    print("  from a high-overlap pair where neither is.")
else:
    print("  No genome carries a callable T6SS.")
    if part:
        print(f"  Partial hits in {len(part)} genome(s): "
              f"{', '.join(mid(x) for x in part)} — most likely Hcp/VgrG-like")
        print("  phage-tail homologues, which are structurally related to T6SS and")
        print("  are exactly the false positive to expect in a genome full of")
        print("  proviruses. Not evidence of a secretion system.")
    print()
    print("  THIS IS A RESULT. Interbacterial antagonism by T6SS appears absent")
    print("  from this catalogue, so step 21's competition edges stay what they")
    print("  are: overlapping resource use, not documented aggression. Say that")
    print("  explicitly rather than leaving the reader to assume conflict.")
with open(os.path.join(S, 't6ss_detail.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id', 'genus', 'n_markers', 'markers', 'call'])
    for m in MAGS:
        h = hits[m].get('T6SS', ())
        if h: w.writerow([mid(m), gen(m), len(h), ','.join(sorted(h)),
                          'present' if len(h) >= 3 else 'partial'])

print()
print("=" * 112)
print("  WHAT TO TAKE FROM THIS")
print("=" * 112)
if not ok:
    print("  >> AT LEAST ONE CONTROL FAILED. Fix that before reading section 1.")
else:
    print("  >> Controls passed. Section 1 is usable as a comparative table.")
print("  The chapter sentence this supports is about DESTINATION, not quantity:")
print("  which compartments this community can reach, and by what route. Read it")
print("  beside 23c and 23d — a signal peptide says 'left the cytoplasm', Tat says")
print("  'left folded', sortase and T9SS/PorU say 'stopped at the surface', and")
print("  only T1SS, T2SS and T5SS say 'reached the medium' in a Gram-negative.")
print()
print("  STILL NOT LOCALISATION. This is gene content. Every statement above is")
print("  about what the genome could do, not what the cell did.")
PY

echo
echo "######## VALIDATION GATE ########"
check "secretion table" $([ -s "$S/secretion_systems.tsv" ] && echo 0 || echo 1)
check "motif table"     $([ -s "$S/motif_substrates.tsv" ] && echo 0 || echo 1)
runlog "23e_secretion complete"
echo
echo "######## 23e DONE ########"
