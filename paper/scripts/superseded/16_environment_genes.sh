#!/usr/bin/env bash
# =============================================================================
# 16_environment_genes.sh  --  what makes these genomes fit ETOLIKO specifically
#
# THE GAP THIS FILLS
#   Steps 11/11b covered metabolism — who eats what, who makes methane. Nothing
#   so far addresses the physical environment. Etoliko is a permanently
#   stratified coastal lagoon: brackish-to-saline, anoxic below the chemocline,
#   PERMANENTLY SULFIDIC, sulfate-replete, and known for metal accumulation in
#   the sediment. Organisms living there need traits that a generic metabolic
#   reconstruction never asks about.
#
#   Screens applied to all 29 MAGs:
#     OSMO   compatible solutes — ectoine, glycine betaine, trehalose, proline,
#            and their transporters. A meromictic lagoon has a salinity
#            gradient; osmolyte strategy is how organisms live across it.
#     SULF   sulfide tolerance and detoxification — sulfide:quinone
#            oxidoreductase, persulfide dioxygenase, rhodanese, sulfurtransferase.
#            Sulfide is toxic; permanent exposure demands machinery.
#     METAL  arsenic, mercury, cadmium/zinc/cobalt, copper, chromate resistance.
#     OXID   oxidative stress — but the ANAEROBIC set (rubrerythrin, superoxide
#            reductase, rubredoxin) as well as the aerobic one (SOD, catalase).
#            Which set an organism carries says whether it tolerates transient
#            oxygen or is a strict anaerobe.
#     PHAGE  prophage and mobile-element signal from the annotations.
#     CRISPR arrays already called by Bakta — a record of past phage exposure.
#
# METHOD AND ITS LIMIT
#   Screens run on the Bakta annotations from step 10, matched on gene symbol
#   and product keywords. This is a SCREEN. Unlike the step-15 markers, nothing
#   here is phylogenetically verified, so results are reported as "annotation
#   evidence for", not as demonstrated capability. Counts are per genome and
#   should be compared BETWEEN genomes rather than read as absolutes.
#
# INPUT    09_annotation/<MAG>/<MAG>.tsv, 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   16_environment/environment_matrix.tsv + printed tables
#
# RUNTIME  under a minute.
# USAGE    ws2 && bash scripts/16_environment_genes.sh 2>&1 | tee logs/16_env.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="16_environment"
A="$WS2/09_annotation"; F="$WS2/07_mags"; E="$WS2/16_environment"; mkdir -p "$E"

echo "######## 16 — ENVIRONMENT-SPECIFIC GENE SCREENS ########"
echo "  started: $(date -Is)"
runlog "16_environment_genes started"
[ -d "$A" ] || { echo "  !! $A missing — run 10_annotate.sh first"; exit 1; }

python3 - "$A" "$F" "$E" <<'PY'
import sys, os, csv, glob, re
from collections import Counter

A, F, E = sys.argv[1], sys.argv[2], sys.argv[3]

def M(cat, key, label, genes='', prods=()):
    return dict(cat=cat, key=key, label=label,
                genes=re.compile(rf'^({genes})$', re.I) if genes else None,
                prods=[p.lower() for p in prods])

MARKERS = [
 # ---- compatible solutes / osmotic adaptation ------------------------------
 M('OSMO','ectABC','ectoine biosynthesis','ectA|ectB|ectC',
   ('ectoine synthase','diaminobutyrate','ectoine')),
 M('OSMO','betAB','glycine betaine synthesis','betA|betB|gbsA|gbsB',
   ('choline dehydrogenase','betaine-aldehyde dehydrogenase')),
 M('OSMO','opuABC','betaine/carnitine/choline transporter','opuA|opuB|opuC|opuD|betT|betP|proP',
   ('glycine betaine','choline transport','betaine/carnitine','osmoprotectant')),
 M('OSMO','proVWX','proline/glycine betaine ABC transport','proV|proW|proX',
   ('glycine betaine/l-proline','proline abc')),
 M('OSMO','otsAB','trehalose biosynthesis','otsA|otsB|treS|treY|treZ',
   ('trehalose-6-phosphate','trehalose synthase')),
 M('OSMO','mscLS','mechanosensitive channel','mscL|mscS|mscM',
   ('mechanosensitive ion channel',)),
 M('OSMO','kdpA','K+ transport (turgor)','kdpA|kdpB|kdpC|trkA|trkH|kup',
   ('potassium-transporting','potassium uptake')),
 # ---- sulfide tolerance / detoxification -----------------------------------
 M('SULF','sqr','sulfide:quinone oxidoreductase','sqr|sqrA',
   ('sulfide-quinone','sulfide:quinone')),
 M('SULF','sdo','persulfide dioxygenase (ETHE1)','ethe1|sdo|blh',
   ('persulfide dioxygenase','sulfur dioxygenase','hydroxyacylglutathione')),
 M('SULF','tst','rhodanese / thiosulfate sulfurtransferase','tst|sseA|glpE|pspE',
   ('rhodanese','thiosulfate sulfurtransferase','sulfurtransferase')),
 M('SULF','fccB','flavocytochrome c sulfide dehydrogenase','fccA|fccB',
   ('flavocytochrome c','sulfide dehydrogenase')),
 M('SULF','cysK','cysteine synthase (sulfide assimilation)','cysK|cysM|cysE',
   ('cysteine synthase','o-acetylserine')),
 M('SULF','sudA','sulfide dehydrogenase / NADPH-sulfur','sudA|sudB|npsr',
   ('sulfide dehydrogenase','sulfur reductase')),
 # ---- metal resistance ------------------------------------------------------
 M('METAL','ars','arsenic resistance','arsA|arsB|arsC|arsR|arsD|acr3',
   ('arsenical','arsenate reductase','arsenite')),
 M('METAL','mer','mercury resistance','merA|merB|merR|merT|merP',
   ('mercuric','organomercurial')),
 M('METAL','czc','cobalt/zinc/cadmium efflux','czcA|czcB|czcC|czcD|cadA|zntA',
   ('cobalt-zinc-cadmium','cadmium-translocating','zinc/cadmium')),
 M('METAL','cop','copper resistance/efflux','copA|copB|copC|copD|cusA|cusB|cueO',
   ('copper-translocating','copper resistance','multicopper oxidase')),
 M('METAL','chrA','chromate efflux','chrA|chrR',
   ('chromate transport','chromate reductase')),
 M('METAL','terC','tellurite resistance','terC|terD|terZ|tehB',
   ('tellurite resistance','tellurium')),
 M('METAL','feoB','ferrous iron transport','feoA|feoB|efeU',
   ('ferrous iron transport',)),
 # ---- oxidative stress: aerobic vs anaerobic strategies ---------------------
 M('OXID','sod','superoxide dismutase (aerobic-type)','sodA|sodB|sodC|sodN',
   ('superoxide dismutase',)),
 M('OXID','kat','catalase','katA|katB|katE|katG|srpA',
   ('catalase',)),
 M('OXID','ahpC','peroxiredoxin / alkyl hydroperoxide reductase','ahpC|ahpF|tpx|bcp|prx',
   ('peroxiredoxin','alkyl hydroperoxide')),
 M('OXID','rbr','rubrerythrin (ANAEROBIC peroxide defence)','rbr|rbr2|rbr3',
   ('rubrerythrin',)),
 M('OXID','sor','superoxide reductase / desulfoferrodoxin (ANAEROBIC)','sor|dfx|rbo|neelaredoxin',
   ('superoxide reductase','desulfoferrodoxin','neelaredoxin')),
 M('OXID','rd','rubredoxin / rubredoxin oxidoreductase','rub|rubA|roo|norV',
   ('rubredoxin','flavodiiron','flavorubredoxin')),
 M('OXID','msrA','methionine sulfoxide reductase','msrA|msrB|msrP|msrQ',
   ('methionine sulfoxide reductase',)),
 # ---- mobile elements / phage ----------------------------------------------
 M('PHAGE','integrase','integrase / recombinase','int|xerC|xerD|intA',
   ('integrase','site-specific recombinase','tyrosine recombinase')),
 M('PHAGE','capsid','phage structural protein','',
   ('phage','capsid','tail fiber','portal protein','terminase','baseplate','holin')),
 M('PHAGE','transposase','transposase / IS element','tnpA|tniA|istA',
   ('transposase','is element','insertion sequence')),
 M('PHAGE','toxin_antitoxin','toxin-antitoxin system','relE|relB|mazE|mazF|vapB|vapC|higA|higB',
   ('toxin-antitoxin','addiction module')),
 M('PHAGE','restriction','restriction-modification','hsdR|hsdM|hsdS|mrr',
   ('restriction endonuclease','type i restriction','dna methyltransferase')),
 M('PHAGE','cas','CRISPR-associated protein','cas1|cas2|cas3|cas9|csn1|cmr|cse',
   ('crispr-associated',)),
]

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

def hit(gene, prod, m):
    if m['genes'] and gene and m['genes'].match(gene): return True
    if m['prods'] and prod and any(x in prod for x in m['prods']): return True
    return False

mags, mat, crispr = [], {}, {}
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)
    tsv = os.path.join(d, f'{name}.tsv')
    if not os.path.exists(tsv): continue
    mags.append(name)
    c = {m['key']: 0 for m in MARKERS}
    ncr = 0
    for line in open(tsv):
        if line.startswith('#') or not line.strip(): continue
        p_ = line.rstrip('\n').split('\t')
        if len(p_) < 8: continue
        ftype = p_[1]
        if ftype.lower() == 'crispr':
            ncr += 1; continue
        if ftype not in ('cds', 'CDS'): continue
        gene = p_[6]; prod = p_[7].lower()
        for m in MARKERS:
            if hit(gene, prod, m): c[m['key']] += 1
    mat[name] = c; crispr[name] = ncr

if not mags:
    print("  !! no annotations parsed"); raise SystemExit(1)
print(f"  screened {len(mags)} MAGs across {len(MARKERS)} gene families")

keys = [m['key'] for m in MARKERS]
with open(os.path.join(E, 'environment_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','genus','crispr_arrays'] + keys)
    for n in mags:
        t = tax.get(n, {})
        w.writerow([t.get('mag_id',''), n, t.get('sample',''), t.get('domain',''),
                    t.get('genus',''), crispr[n]] + [mat[n][k] for k in keys])

def show(title, cat, extra=None):
    ks = [m['key'] for m in MARKERS if m['cat'] == cat]
    print()
    print("=" * 104)
    print(f"  {title}")
    print("=" * 104)
    hdr = ''.join(f"{k[:8]:>10}" for k in ks) + (f"{extra[0]:>10}" if extra else '')
    print(f"  {'MAG ID':<17}{'dom':<4}{'genus':<18}{hdr}")
    for n in sorted(mags, key=lambda n: (tax.get(n,{}).get('sample','')!='SED',
                                         tax.get(n,{}).get('mag_id',''))):
        t = tax.get(n, {})
        cells = ''.join(f"{(mat[n][k] or '·'):>10}" for k in ks)
        if extra: cells += f"{(crispr[n] or '·'):>10}"
        print(f"  {t.get('mag_id',n)[:17]:<17}{(t.get('domain','?') or '?')[:3]:<4}"
              f"{(t.get('genus','-') or '-')[:18]:<18}{cells}")

show("OSMOTIC ADAPTATION — compatible solutes and transport", 'OSMO')
show("SULFIDE TOLERANCE AND DETOXIFICATION", 'SULF')
show("METAL RESISTANCE", 'METAL')
show("OXIDATIVE STRESS — aerobic (sod/kat/ahpC) vs anaerobic (rbr/sor/rd)", 'OXID')
show("MOBILE ELEMENTS AND PHAGE DEFENCE", 'PHAGE', extra=('CRISPR',))

# ------------------------------------------------- interpretation helpers
print()
print("=" * 104)
print("  OXYGEN STRATEGY — which defence system does each genome carry?")
print("=" * 104)
print(f"  {'MAG ID':<17}{'genus':<20}{'aerobic':>9}{'anaerobic':>11}   inference")
for n in sorted(mags, key=lambda n: tax.get(n,{}).get('mag_id','')):
    t = tax.get(n, {})
    aer = mat[n]['sod'] + mat[n]['kat'] + mat[n]['ahpC']
    ana = mat[n]['rbr'] + mat[n]['sor'] + mat[n]['rd']
    if aer and not ana:   inf = 'aerotolerant / facultative'
    elif ana and not aer: inf = 'strict anaerobe (anaerobic ROS defence only)'
    elif aer and ana:     inf = 'both systems — transient oxygen exposure'
    else:                 inf = 'neither detected'
    print(f"  {t.get('mag_id',n):<17}{(t.get('genus','-') or '-')[:20]:<20}"
          f"{aer:>9}{ana:>11}   {inf}")

print()
print("=" * 104)
print("  SUMMARY BY CATEGORY (MAGs carrying at least one gene in the family)")
print("=" * 104)
for cat, lab in (('OSMO','osmotic adaptation'), ('SULF','sulfide handling'),
                 ('METAL','metal resistance'), ('OXID','oxidative stress'),
                 ('PHAGE','mobile elements')):
    ks = [m['key'] for m in MARKERS if m['cat'] == cat]
    print(f"\n  {lab.upper()}")
    for k in ks:
        lab2 = next(m['label'] for m in MARKERS if m['key'] == k)
        n_mags = sum(1 for n in mags if mat[n][k] > 0)
        tot = sum(mat[n][k] for n in mags)
        bar = '#' * int(round(20 * n_mags / len(mags)))
        print(f"    {k:<16}{n_mags:>3}/{len(mags)} MAGs  {tot:>5} genes  {bar:<20} {lab2}")

ncr = sum(1 for n in mags if crispr[n] > 0)
print()
print(f"  CRISPR arrays: present in {ncr}/{len(mags)} MAGs "
      f"({sum(crispr.values())} arrays total) — a record of past phage exposure "
      f"in a system where no viral analysis has been done.")
print()
print("  CAVEAT: these are ANNOTATION screens, not phylogenetically verified")
print("  markers like the step-15 mcrA/mttB calls. Report them as 'genomic")
print("  evidence for', compare counts BETWEEN genomes rather than in absolute")
print("  terms, and verify any family that becomes central to an argument.")
PY

echo
echo "######## VALIDATION GATE ########"
check "environment matrix written" $([ -s "$E/environment_matrix.tsv" ] && echo 0 || echo 1)
n=$(( $(wc -l < "$E/environment_matrix.tsv") - 1 ))
nm=$(ls -1 "$WS2/07_mags/bins"/*.fa | wc -l)
check "every MAG screened ($n/$nm)" $([ "$n" -eq "$nm" ] && echo 0 || echo 1)

runlog "16_environment_genes complete"
echo
echo "######## 16 DONE ########"
