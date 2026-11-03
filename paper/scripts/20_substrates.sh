#!/usr/bin/env bash
# =============================================================================
# 20_substrates.sh  --  what this community actually degrades
#
# THE HYPOTHESIS THIS SCRIPT EXISTS TO TEST
#   Step 11b found K14080 (methylated-thiol:CoM methyltransferase, mtaA/mtsA
#   family) at 8 copies in MTP4, 6 in ETO_LIQ_MAG_010 and 5 in _018 — the three
#   methylotrophic Methanosarcinaceae, and nowhere else. Multiple copies of that
#   family indicate DIVERSE methylated substrates, not methanol alone:
#   specifically dimethylsulfide (DMS) and methanethiol (MeSH).
#
#   Etoliko is a coastal saline lagoon with algal input, so DMSP
#   (dimethylsulfoniopropionate) should be abundant. Bacteria cleave DMSP to DMS
#   (ddd genes) or demethylate it to MeSH (dmdA). Both products are
#   NON-COMPETITIVE methanogenesis substrates — sulfate reducers cannot use
#   them. In a permanently sulfate-replete lagoon with no SMTZ, that is exactly
#   the pathway audit finding F20 predicts.
#
#   PREDICTION: bacterial MAGs carry DMSP degradation genes AND the three
#   Methanosarcinaceae carry methylated-thiol methyltransferases. If both hold,
#   the community encodes a complete algal-DMSP -> DMS/MeSH -> CH4 route that
#   bypasses sulfate reduction entirely. If the bacterial half is absent, the
#   archaeal copies are for methanol/methylamines only and the DMSP idea dies.
#
# EVERYTHING ELSE IT COVERS — degradation beyond the obvious
#   CAZY    glycoside hydrolases, polysaccharide lyases, carbohydrate esterases,
#           auxiliary activities, plus named targets (cellulose, xylan, chitin,
#           pectin, alginate, agar, starch, laminarin, fucoidan)
#   PEPT    peptidases by catalytic class, extracellular vs cytoplasmic
#   LIPID   lipases, esterases, phospholipases, fatty acid beta-oxidation
#   AROM    ANAEROBIC aromatic degradation — benzoyl-CoA reductase and the
#           central benzoyl-CoA pathway, not the aerobic ring dioxygenases that
#           dominate database keyword searches and cannot operate here
#   METHYL  DMSP/DMS/MeSH, TMA/TMAO, glycine betaine, choline, methanol
#   ONE_C   formate, CO, methylamines
#   P_N     phosphonates (phn), sulfonates (ssu/tau), urea, cyanate
#   NUCL    nucleoside and purine/pyrimidine degradation
#
# METHOD
#   KOfam HMMs where a KO exists (profiles cached from step 19's archive), plus
#   Bakta product/gene matching for families KEGG does not resolve. Each call
#   records WHICH evidence supported it, so KO-backed and keyword-backed hits
#   are never conflated.
#
# INPUT    09_annotation/, 11_kofam/profiles.tar.gz, 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   20_substrates/
#
# RUNTIME  20-40 min.
# USAGE    ws2 && bash scripts/20_substrates.sh 2>&1 | tee logs/20_substrates.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="20_substrates"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"; U="$WS2/20_substrates"
mkdir -p "$U/profiles" "$U/hits"

echo "######## 20 — SUBSTRATE DEGRADATION AND METHYLATED COMPOUNDS ########"
echo "  started: $(date -Is)"
runlog "20_substrates started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }

# ---- KO families ------------------------------------------------------------
cat > "$U/ko_targets.tsv" <<'KOS'
K17486	dmdA	METHYL_S	DMSP demethylase — DMSP to MeSH
K16953	dddD	METHYL_S	DMSP CoA-transferase/lyase
K17487	dddP	METHYL_S	DMSP lyase DddP — DMSP to DMS
K16954	dddL	METHYL_S	DMSP lyase DddL
K17950	dddQ	METHYL_S	DMSP lyase DddQ
K17951	dddW	METHYL_S	DMSP lyase DddW
K14080	mtsA	METHYL_S	methylated-thiol:CoM methyltransferase (DMS/MeSH)
K16179	mtsB	METHYL_S	methylated-thiol methyltransferase corrinoid
K17218	sqr	METHYL_S	sulfide:quinone oxidoreductase
K00399	mcrA	METHYL_S	methyl-CoM reductase (terminal step, context)
K00551	betA	METHYL_N	choline dehydrogenase
K00130	betB	METHYL_N	betaine-aldehyde dehydrogenase
K14083	mttB	METHYL_N	MttB superfamily (TMA MT in Methanosarcinaceae; betaine MT in bacteria)
K16177	mtbB	METHYL_N	dimethylamine methyltransferase
K16157	mtmB	METHYL_N	monomethylamine methyltransferase
K07811	torA	METHYL_N	TMAO reductase
K03532	tdm	METHYL_N	trimethylamine dehydrogenase
K04480	mtaB	METHYL_N	methanol:corrinoid methyltransferase
K17067	mdh2	METHYL_N	methanol dehydrogenase
K04112	bcrA	AROM	benzoyl-CoA reductase alpha (ANAEROBIC ring reduction)
K04113	bcrB	AROM	benzoyl-CoA reductase beta
K04114	bcrC	AROM	benzoyl-CoA reductase gamma
K04115	bcrD	AROM	benzoyl-CoA reductase delta
K04116	bamB	AROM	class II benzoyl-CoA reductase BamB
K07543	badH	AROM	2-hydroxycyclohexanecarboxyl-CoA dehydrogenase
K01617	hbd	AROM	4-hydroxybenzoyl-CoA reductase
K01912	bclA	AROM	benzoate-CoA ligase
K05939	phnJ	P_N	C-P lyase PhnJ (phosphonate)
K06163	phnG	P_N	C-P lyase PhnG
K02041	phnC	P_N	phosphonate ABC transporter
K00299	ssuE	P_N	FMN reductase, sulfonate desulfonation
K04091	ssuD	P_N	alkanesulfonate monooxygenase
K03119	tauD	P_N	taurine dioxygenase
K01427	ureC	P_N	urease alpha
K01725	cynS	P_N	cyanase
K01179	celA	CAZY	endoglucanase (cellulose)
K01181	xynA	CAZY	endo-1,4-beta-xylanase
K01183	chiA	CAZY	chitinase
K01213	pelA	CAZY	pectate lyase / pectinesterase
K01176	amyA	CAZY	alpha-amylase
K01207	nagZ	CAZY	beta-N-acetylhexosaminidase
K01190	lacZ	CAZY	beta-galactosidase
K05349	bglX	CAZY	beta-glucosidase
K01238	aglA	CAZY	alpha-glucosidase
K01811	melA	CAZY	alpha-galactosidase
K01192	manB	CAZY	beta-mannosidase
K01187	malZ	CAZY	alpha-glucosidase (maltose)
K01046	lip	LIPID	triacylglycerol lipase
K01175	pldB	LIPID	lysophospholipase
K01058	plaA	LIPID	phospholipase A
K00232	fadE	LIPID	acyl-CoA oxidase / dehydrogenase
K01782	fadB	LIPID	fatty acid beta-oxidation multifunctional
K00626	atoB	LIPID	acetyl-CoA acetyltransferase
K01442	cbh	PEPT	choloylglycine hydrolase
K01417	pepD	PEPT	dipeptidase D
K01262	pepP	PEPT	Xaa-Pro aminopeptidase
K01179	x1	PEPT	placeholder
K00762	pyrE	NUCL	orotate phosphoribosyltransferase
K01464	dht	NUCL	dihydropyrimidinase
K01487	guaD	NUCL	guanine deaminase
K01477	allB	NUCL	allantoinase
K00365	uox	NUCL	urate oxidase
KOS
NKO=$(wc -l < "$U/ko_targets.tsv")
echo "  KO families: $NKO"

# ---- profiles from the cached archive ---------------------------------------
BASE="https://www.genome.jp/ftp/db/kofam"
[ -s "$U/ko_list" ] || cp -n "$K/ko_list" "$U/ko_list" 2>/dev/null || true
got=0
while IFS=$'\t' read -r ko rest; do
  [ -z "$ko" ] && continue
  [ -s "$U/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  if [ -s "$K/profiles/${ko}.hmm" ]; then cp "$K/profiles/${ko}.hmm" "$U/profiles/"; got=$((got+1)); fi
done < "$U/ko_targets.tsv"
if [ "$got" -lt "$NKO" ]; then
  TARB="$K/profiles.tar.gz"
  if [ -s "$TARB" ]; then
    echo "  extracting missing families from the cached KOfam archive ..."
    awk -F'\t' '{print "profiles/"$1".hmm"}' "$U/ko_targets.tsv" > "$U/wanted.txt"
    tar -xzf "$TARB" -C "$U" -T "$U/wanted.txt" 2>/dev/null || true
    got=$(ls -1 "$U/profiles"/*.hmm 2>/dev/null | wc -l)
  else
    echo "  !! $TARB missing — run 19_auxotrophy.sh first (it caches the archive)"
  fi
fi
echo "  profiles available: $got/$NKO"
[ "$got" -gt 10 ] || { echo "  !! too few profiles"; exit 1; }
cat "$U/profiles"/*.hmm > "$U/targets.hmm"

# ---- search ------------------------------------------------------------------
echo
echo "######## hmmsearch ($(date +%H:%M)) ########"
n=0; tot=$(ls -1d "$A"/*/ | wc -l)
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1)); out="$U/hits/${name}.tbl"
  [ -s "$out" ] && { printf '\r  [%2d/%2d] %-28s cached' "$n" "$tot" "$name"; continue; }
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 "$U/targets.hmm" "$faa" \
    >/dev/null 2>&1 || echo "  !! failed on $name"
done
echo

python3 - "$U" "$F" "$A" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter

U, F, A = sys.argv[1], sys.argv[2], sys.argv[3]

meta = {}
for line in open(os.path.join(U, 'ko_targets.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 4: meta[p[0]] = dict(gene=p[1], cat=p[2], desc=p[3])

thr = {}
kl = os.path.join(U, 'ko_list')
if os.path.exists(kl):
    for r in csv.DictReader(open(kl), delimiter='\t'):
        ko = r.get('knum') or r.get('#knum') or ''
        t = r.get('threshold', '')
        if ko in meta and t not in ('', '-'):
            try: thr[ko] = float(t)
            except ValueError: pass

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

# ---- KO hits
ko_hits = {}
for tbl in sorted(glob.glob(os.path.join(U, 'hits', '*.tbl'))):
    name = os.path.basename(tbl)[:-4]
    c = Counter()
    for line in open(tbl):
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) < 6: continue
        ko = f[2].split('.')[0]
        if ko not in meta: continue
        try: score, ev = float(f[5]), float(f[4])
        except ValueError: continue
        if (score >= thr[ko]) if ko in thr else (ev <= 1e-5):
            c[ko] += 1
    ko_hits[name] = c

# ---- Bakta keyword evidence for families KEGG does not resolve
KEYWORDS = {
 'CAZY_GH':   ('glycoside hydrolase','glycosyl hydrolase','glucanase','xylanase','chitinase',
               'amylase','glucosidase','galactosidase','mannanase','fucosidase','agarase',
               'laminarinase','cellulase','pullulanase'),
 'CAZY_PL':   ('polysaccharide lyase','pectate lyase','alginate lyase','heparin lyase'),
 'CAZY_CE':   ('carbohydrate esterase','acetyl xylan esterase','pectin methylesterase'),
 'CAZY_CBM':  ('carbohydrate-binding','cellulosome','cohesin','dockerin'),
 'SUS':       ('susc','susd','ton-b dependent receptor','tonb-dependent'),
 'PEPT_ser':  ('serine protease','subtilisin','s8 family','trypsin-like'),
 'PEPT_met':  ('metallopeptidase','metalloprotease','m20 family','m24 family','zinc protease'),
 'PEPT_cys':  ('cysteine protease','papain','c1 family'),
 'PEPT_asp':  ('aspartic protease','pepsin'),
 'PEPT_exo':  ('aminopeptidase','carboxypeptidase','dipeptidase','tripeptidase','oligopeptidase'),
 'LIPASE':    ('lipase','esterase','phospholipase','cutinase'),
 'AROM_anae': ('benzoyl-coa','benzoate-coa','4-hydroxybenzoyl','phenylacetyl-coa',
               'cyclohexanecarboxyl'),
 'DMSP':      ('dmsp','dimethylsulfoniopropionate','dimethylsulfide','methanethiol',
               'methylmercaptan','dimethyl sulfide'),
 'TMA':       ('trimethylamine','tmao','dimethylamine','methylamine'),
 'BETAINE':   ('glycine betaine','betaine'),
 'PHOSPHONATE': ('phosphonate','c-p lyase','phosphonatase'),
 'SULFONATE': ('sulfonate','taurine','sulfoacetaldehyde'),
 'UREA':      ('urease','urea '),
}
kw_hits = {}
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    name = os.path.basename(d); tsv = os.path.join(d, f'{name}.tsv')
    if not os.path.exists(tsv): continue
    c = Counter()
    for line in open(tsv):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8 or f[1] not in ('cds','CDS'): continue
        prod = f[7].lower()
        for k, pats in KEYWORDS.items():
            if any(x in prod for x in pats): c[k] += 1
    kw_hits[name] = c

mags = sorted(set(ko_hits) | set(kw_hits))
print(f"  MAGs profiled: {len(mags)}")

# =====================================================================
# THE DMSP TEST
# =====================================================================
print()
print("=" * 104)
print("  THE DMSP -> DMS/MeSH -> CH4 HYPOTHESIS")
print("=" * 104)
BACT_DMSP = ['K17486','K16953','K17487','K16954','K17950','K17951']
ARCH_MT   = ['K14080','K16179']
print("  BACTERIAL HALF — DMSP cleavage / demethylation")
found_b = []
for m in mags:
    n = sum(ko_hits.get(m, Counter())[k] for k in BACT_DMSP)
    kw = kw_hits.get(m, Counter())['DMSP']
    if n or kw:
        t = tax.get(m, {})
        genes = ','.join(f"{meta[k]['gene']}x{ko_hits[m][k]}" for k in BACT_DMSP if ko_hits.get(m, Counter())[k])
        print(f"    {t.get('mag_id',m):<17} {(t.get('domain','?') or '?')[:3]:<4} "
              f"{(t.get('genus','-') or '-')[:20]:<20} KO:{n} kw:{kw}  {genes}")
        found_b.append(m)
if not found_b:
    print("    NONE — no bacterial DMSP degradation detected")

print()
print("  ARCHAEAL HALF — methylated-thiol:CoM methyltransferase (DMS/MeSH -> CH4)")
found_a = []
for m in mags:
    n = sum(ko_hits.get(m, Counter())[k] for k in ARCH_MT)
    if n:
        t = tax.get(m, {})
        genes = ','.join(f"{meta[k]['gene']}x{ko_hits[m][k]}" for k in ARCH_MT if ko_hits[m][k])
        print(f"    {t.get('mag_id',m):<17} {(t.get('domain','?') or '?')[:3]:<4} "
              f"{(t.get('genus','-') or '-')[:20]:<20} {genes}")
        found_a.append(m)

print()
if found_b and found_a:
    print("  >> BOTH HALVES PRESENT. The community encodes a complete route from")
    print("     DMSP to methane that bypasses sulfate reduction entirely — the")
    print("     non-competitive pathway predicted for a sulfate-replete lagoon.")
    print("     NEXT: this is a genomic potential, not a measured flux. Confirm by")
    print("     measuring DMSP/DMS in porewater, or by amending an enrichment with")
    print("     DMSP and following methane. That is a small, decisive experiment.")
elif found_a and not found_b:
    print("  >> ARCHAEAL HALF ONLY. The methyltransferase copies are most likely")
    print("     for methanol and methylamines. The DMSP route is NOT supported;")
    print("     do not claim it.")
else:
    print("  >> Hypothesis not supported by the gene content.")

# =====================================================================
# SUBSTRATE TABLES
# =====================================================================
def ko_table(cat, title):
    ks = [k for k, v in meta.items() if v['cat'] == cat]
    ks = [k for k in ks if any(ko_hits.get(m, Counter())[k] for m in mags)]
    if not ks: return
    print()
    print("=" * 104)
    print(f"  {title}")
    print("=" * 104)
    print(f"  {'MAG ID':<17}{'genus':<20}" + ''.join(f"{meta[k]['gene'][:7]:>9}" for k in ks))
    for m in sorted(mags, key=lambda m: (tax.get(m,{}).get('sample','')!='SED',
                                         tax.get(m,{}).get('mag_id',''))):
        t = tax.get(m, {})
        cells = ''.join(f"{(ko_hits.get(m,Counter())[k] or '·'):>9}" for k in ks)
        print(f"  {t.get('mag_id',m)[:17]:<17}{(t.get('genus','-') or '-')[:20]:<20}{cells}")

ko_table('METHYL_S', 'METHYLATED SULFUR COMPOUNDS')
ko_table('METHYL_N', 'METHYLATED NITROGEN COMPOUNDS + METHANOL')
ko_table('AROM',     'ANAEROBIC AROMATIC DEGRADATION (benzoyl-CoA pathway)')
ko_table('P_N',      'PHOSPHONATES, SULFONATES, UREA, CYANATE')
ko_table('CAZY',     'NAMED CARBOHYDRATE-ACTIVE ENZYMES (KO-backed)')
ko_table('LIPID',    'LIPID AND FATTY ACID DEGRADATION')
ko_table('NUCL',     'NUCLEOTIDE / PURINE DEGRADATION')

print()
print("=" * 104)
print("  DEGRADATION CAPACITY BY CLASS (annotation keyword evidence, gene counts)")
print("=" * 104)
cols = ['CAZY_GH','CAZY_PL','CAZY_CE','CAZY_CBM','SUS','PEPT_ser','PEPT_met',
        'PEPT_cys','PEPT_exo','LIPASE','AROM_anae']
print(f"  {'MAG ID':<17}{'genus':<20}" + ''.join(f"{c.replace('CAZY_','').replace('PEPT_','p')[:6]:>8}" for c in cols))
for m in sorted(mags, key=lambda m: -sum(kw_hits.get(m, Counter())[c] for c in cols)):
    t = tax.get(m, {})
    cells = ''.join(f"{(kw_hits.get(m,Counter())[c] or '·'):>8}" for c in cols)
    print(f"  {t.get('mag_id',m)[:17]:<17}{(t.get('genus','-') or '-')[:20]:<20}{cells}")

# ---- write everything
allk = sorted(meta)
with open(os.path.join(U, 'substrate_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','genus']
               + [f"KO_{meta[k]['gene']}" for k in allk]
               + [f"KW_{c}" for c in sorted(KEYWORDS)])
    for m in mags:
        t = tax.get(m, {})
        w.writerow([t.get('mag_id',''), m, t.get('sample',''), t.get('domain',''),
                    t.get('genus','')]
                   + [ko_hits.get(m, Counter())[k] for k in allk]
                   + [kw_hits.get(m, Counter())[c] for c in sorted(KEYWORDS)])

print()
print("  EVIDENCE CLASSES ARE KEPT SEPARATE in substrate_matrix.tsv: KO_* columns")
print("  are HMM hits against KOfam adaptive thresholds; KW_* columns are Bakta")
print("  product-string matches. KO evidence is much stronger. Never merge them")
print("  into a single 'this MAG degrades X' count.")
PY

echo
echo "######## VALIDATION GATE ########"
check "substrate matrix" $([ -s "$U/substrate_matrix.tsv" ] && echo 0 || echo 1)
n=$(( $(wc -l < "$U/substrate_matrix.tsv") - 1 ))
nm=$(ls -1 "$WS2/07_mags/bins"/*.fa | wc -l)
check "every MAG profiled ($n/$nm)" $([ "$n" -eq "$nm" ] && echo 0 || echo 1)

runlog "20_substrates complete"
echo
echo "######## 20 DONE ########"
echo "The DMSP block is the one that decides a new hypothesis."
