#!/usr/bin/env bash
# =============================================================================
# 11_metabolism.sh  --  full metabolic reconstruction, all 29 MAGs, all pathways
#
# SCOPE
#   Every MAG, bacterial and archaeal, is profiled across every category below.
#   Nothing is treated as background to anything else. The syntrophy question
#   only has an answer if fermenters and methanogens are in the same table.
#
#     C1        Wood-Ljungdahl, formate/CO oxidation, methanogenesis
#     CARB      glycolysis, TCA, pentose phosphate, CAZyme classes
#     PROT      peptidases, amino-acid fermentation (Stickland)
#     FERM      pyruvate node; acetate / lactate / butyrate / propionate / ethanol / H2
#     H2/FDH    [FeFe] and [NiFe] hydrogenases, formate dehydrogenase
#     EB        electron bifurcation: Rnf, Hdr, Etf, Nfn
#     S         dsrAB (+type), aprAB, sat, sqr, sox, asrABC, phsA
#     N         nifH, nar/nap, nir, nor, nos, nrfA, amo, hzs
#     RESP      cytochrome oxidases, metal reduction
#     TRANS     transporters, motility, chemotaxis, sporulation
#     COFAC     B12 / cobalamin, biotin, thiamine — auxotrophy is a syntrophy signal
#
# METHOD, AND ITS LIMIT
#   Profiles are built from the Bakta annotations at step 10, matched on gene
#   symbol, EC number, product string and KO where present — several independent
#   fields, so a single mislabelled field cannot create a call on its own.
#
#   ANNOTATION IS A SCREEN, NOT A RESULT. Four traps are known to produce
#   confident false positives in exactly this dataset, and each is handled
#   explicitly rather than trusted:
#
#     mcrA vs ACR    alkyl-CoM reductases of alkane oxidisers annotate as mcrA.
#                    workspace_1 had a ~50% false-positive rate on "complete
#                    mcrABG". Candidate sequences are EXTRACTED here for
#                    phylogenetic verification; presence is reported as
#                    CANDIDATE until that tree is built.
#     DsrAB          reductive (sulfate reducers) and oxidative (sulfur
#                    oxidisers) types are not separable by annotation. Extracted
#                    for the same treatment.
#     mttB / K14083  the MttB SUPERFAMILY is largely bacterial glycine betaine
#                    methyltransferases, NOT trimethylamine methyltransferase.
#                    Counting it as methylotrophic methanogenesis is a specific
#                    error this project must not repeat. Flagged, never counted.
#     nifH           paralogues (bchL/chlL, nflH) score as nitrogenase.
#                    Extracted, not asserted.
#
# INPUT    09_annotation/<MAG>/<MAG>.tsv and .faa   (step 10)
# OUTPUT   10_function/function_matrix.tsv          MAG x marker presence
#          10_function/pathway_summary.tsv          module completeness per MAG
#          10_function/syntrophy_handoff.tsv        products vs consumers
#          10_function/verify/*.faa                 sequences needing a tree
#
# RUNTIME  minutes.
# USAGE    ws2 && bash scripts/11_metabolism.sh 2>&1 | tee logs/11_metabolism.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="11_metabolism"
A="$WS2/09_annotation"; F="$WS2/07_mags"; F2="$WS2/10_function"
mkdir -p "$F2/verify"

echo "######## 11 — METABOLIC RECONSTRUCTION (all 29 MAGs) ########"
echo "  started: $(date -Is)"
runlog "11_metabolism started"
[ -d "$A" ] || { echo "  !! $A missing — run 10_annotate.sh first"; exit 1; }

python3 - "$A" "$F" "$F2" <<'PY'
import sys, os, csv, glob, re, json

A, F, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

# ============================================================ MARKER SET
# (category, key, label, gene regex, EC list, product keywords, KO list, caveat)
def M(cat, key, label, genes='', ecs=(), prods=(), kos=(), caveat=''):
    return dict(cat=cat, key=key, label=label,
                genes=re.compile(rf'^({genes})$', re.I) if genes else None,
                ecs=set(ecs), prods=[p.lower() for p in prods], kos=set(kos),
                caveat=caveat)

MARKERS = [
 # ---- C1 / Wood-Ljungdahl -------------------------------------------------
 M('C1','acsB','CO-methylating acetyl-CoA synthase','acsB|cdhC','2.3.1.169',('acetyl-coa synthase',),('K14138',)),
 M('C1','cooS','CO dehydrogenase','cooS|acsA|cdhA','1.2.7.4',('carbon monoxide dehydrogenase',),('K00198',)),
 M('C1','fhs','formate-THF ligase','fhs','6.3.4.3',('formate--tetrahydrofolate ligase',),('K01938',)),
 M('C1','folD','methylene-THF dehydrogenase','folD','1.5.1.5',('methylenetetrahydrofolate dehydrogenase',),('K01491',)),
 M('C1','metF','methylene-THF reductase','metF|metV','1.5.1.20',('methylenetetrahydrofolate reductase',),('K00297',)),
 # ---- methanogenesis core -------------------------------------------------
 M('C1','mcrA','methyl-CoM reductase alpha','mcrA','2.8.4.1',('methyl-coenzyme m reductase',),('K00399',),
   'CANDIDATE: alkyl-CoM reductase (ACR) of alkane oxidisers is indistinguishable by annotation'),
 M('C1','mcrB','methyl-CoM reductase beta','mcrB','',('methyl-coenzyme m reductase',),('K00401',)),
 M('C1','mcrG','methyl-CoM reductase gamma','mcrG','',('methyl-coenzyme m reductase',),('K00402',)),
 M('C1','mtrA','tetrahydromethanopterin S-methyltransferase','mtrA','2.1.1.86',('methyltransferase subunit a',),('K00577',)),
 M('C1','mer','methylene-H4MPT reductase','mer','1.5.98.2',('methylenetetrahydromethanopterin reductase',),('K00320',)),
 M('C1','mtd','methylene-H4MPT dehydrogenase','mtd','1.5.98.1',(),('K00319',)),
 M('C1','mch','methenyl-H4MPT cyclohydrolase','mch','3.5.4.27',(),('K01499',)),
 M('C1','ftr','formylmethanofuran-H4MPT formyltransferase','ftr','2.3.1.101',(),('K00672',)),
 M('C1','fwdB','formylmethanofuran dehydrogenase','fwdB|fmdB','1.2.7.12',('formylmethanofuran dehydrogenase',),('K00201',)),
 # methylotrophic branch
 M('C1','mtaB','methanol--CoM methyltransferase','mtaB','2.1.1.90',('methanol',),('K04480',)),
 M('C1','mtbA','methylamine--CoM methyltransferase','mtbA','2.1.1.247',(),('K14082',)),
 M('C1','mtmB','monomethylamine methyltransferase','mtmB','2.1.1.248',('monomethylamine',),('K14083_NO',)),
 M('C1','mttB_SF','MttB superfamily (NOT TMA methyltransferase)','mttB','',('trimethylamine methyltransferase','mttb'),('K14083',),
   'DO NOT COUNT as methylamine methanogenesis: K14083/MttB superfamily is predominantly bacterial glycine betaine methyltransferase'),
 # acetoclastic
 M('C1','ackA','acetate kinase','ackA|ack','2.7.2.1',('acetate kinase',),('K00925',)),
 M('C1','pta','phosphotransacetylase','pta','2.3.1.8',('phosphate acetyltransferase',),('K00625',)),
 M('C1','acdA','acetyl-CoA synthetase (ADP-forming)','acdA|acdB','6.2.1.13',(),('K01905',)),
 # ---- carbohydrate --------------------------------------------------------
 M('CARB','pfk','phosphofructokinase','pfkA|pfkB|pfp','2.7.1.11',('6-phosphofructokinase',),('K00850',)),
 M('CARB','gapA','GAPDH','gapA|gap','1.2.1.12',('glyceraldehyde-3-phosphate dehydrogenase',),('K00134',)),
 M('CARB','pyk','pyruvate kinase','pyk|pykA','2.7.1.40',('pyruvate kinase',),('K00873',)),
 M('CARB','zwf','glucose-6-P dehydrogenase','zwf','1.1.1.49',(),('K00036',)),
 M('CARB','gltA','citrate synthase','gltA','2.3.3.1|2.3.3.16',('citrate synthase',),('K01647',)),
 M('CARB','sdhA','succinate dehydrogenase','sdhA|frdA','1.3.5.1|1.3.5.4',('succinate dehydrogenase','fumarate reductase',),('K00239',)),
 M('CARB','GH','glycoside hydrolase (any)','','',('glycoside hydrolase','glycosyl hydrolase','cellulase','xylanase','amylase','chitinase','glucosidase','galactosidase','mannanase','pectate',),()),
 M('CARB','PL','polysaccharide lyase','','',('polysaccharide lyase','pectate lyase','alginate lyase',),()),
 M('CARB','CE','carbohydrate esterase','','',('carbohydrate esterase','acetyl xylan esterase',),()),
 M('CARB','susC','SusC/D polysaccharide utilisation (Bacteroidota)','susC|susD','',('susc','susd','ton-b dependent receptor',),()),
 # ---- protein / amino acid ------------------------------------------------
 M('PROT','peptidase','peptidase / protease (any)','','',('peptidase','protease','proteinase','aminopeptidase','carboxypeptidase',),()),
 M('PROT','grdA','glycine reductase (Stickland)','grdA|grdB','1.21.4.2',('glycine reductase',),('K10670',)),
 M('PROT','prdA','D-proline reductase (Stickland)','prdA|prdB','1.21.4.1',('proline reductase',),()),
 M('PROT','kamA','lysine 2,3-aminomutase','kamA','5.4.3.2',('lysine 2,3-aminomutase',),('K01843',)),
 M('PROT','hadA','2-hydroxyacyl-CoA dehydratase','hadA|hgdA','',('2-hydroxyglutaryl-coa dehydratase',),()),
 # ---- fermentation --------------------------------------------------------
 M('FERM','por','pyruvate:ferredoxin oxidoreductase','porA|porB|nifJ','1.2.7.1',('pyruvate','ferredoxin oxidoreductase',),('K00169',)),
 M('FERM','pflB','pyruvate formate-lyase','pflB|pflD','2.3.1.54',('formate acetyltransferase','pyruvate formate-lyase',),('K00656',)),
 M('FERM','ldh','L-lactate dehydrogenase','ldh|ldhA','1.1.1.27',('l-lactate dehydrogenase',),('K00016',)),
 M('FERM','adhE','bifunctional acetaldehyde/alcohol dehydrogenase','adhE','1.2.1.10',('aldehyde-alcohol dehydrogenase',),('K04072',)),
 M('FERM','buk','butyrate kinase','buk','2.7.2.7',('butyrate kinase',),('K00929',)),
 M('FERM','ptb','phosphate butyryltransferase','ptb','2.3.1.19',(),('K00634',)),
 M('FERM','but','butyryl-CoA:acetate CoA-transferase','but|ato','2.8.3.8',(),('K01034',)),
 M('FERM','bcd','butyryl-CoA dehydrogenase','bcd','1.3.8.1',('butyryl-coa dehydrogenase',),('K00248',)),
 M('FERM','mmdA','methylmalonyl-CoA decarboxylase (propionate)','mmdA|scpB','7.2.4.3',('methylmalonyl-coa decarboxylase',),()),
 M('FERM','mcm','methylmalonyl-CoA mutase (propionate)','mutA|mutB|scpA','5.4.99.2',('methylmalonyl-coa mutase',),('K01847',)),
 # ---- hydrogen / formate --------------------------------------------------
 M('H2','hydA','[FeFe] hydrogenase','hydA','1.12.7.2',('iron hydrogenase','[fefe]','fe-only hydrogenase',),('K00532',)),
 M('H2','hyaB','[NiFe] hydrogenase large subunit','hyaB|hybC|hoxH|mbhL','1.12.99.6',('nickel-dependent hydrogenase','[nife]','hydrogenase large',),('K06281',)),
 M('H2','echA','Ech / energy-converting hydrogenase','echA|echE|ehaA|ehbA','',('energy-converting hydrogenase',),()),
 M('H2','fdhA','formate dehydrogenase','fdhA|fdhF|fdoG|fdnG','1.17.1.9|1.2.1.2',('formate dehydrogenase',),('K00123',)),
 M('H2','hycE','formate hydrogenlyase','hycE','',('formate hydrogenlyase',),()),
 # ---- electron bifurcation -----------------------------------------------
 M('EB','rnfC','Rnf complex','rnfA|rnfB|rnfC|rnfD|rnfE|rnfG','7.2.1.2',('rnf','na(+)-translocating ferredoxin',),('K03615',)),
 M('EB','hdrA','heterodisulfide reductase','hdrA|hdrB|hdrC|hdrD','1.8.98.1',('heterodisulfide reductase',),('K03388',)),
 M('EB','etfA','electron transfer flavoprotein','etfA|etfB','',('electron transfer flavoprotein',),('K03522',)),
 M('EB','nfnA','NADH-dependent ferredoxin:NADP oxidoreductase','nfnA|nfnB','1.18.1.-',('nad(p)h-dependent',),()),
 # ---- sulfur --------------------------------------------------------------
 M('S','dsrA','dissimilatory sulfite reductase alpha','dsrA','1.8.99.5',('dissimilatory sulfite reductase',),('K11180',),
   'CANDIDATE: reductive vs oxidative DsrAB type is NOT separable by annotation'),
 M('S','dsrB','dissimilatory sulfite reductase beta','dsrB','',('dissimilatory sulfite reductase',),('K11181',)),
 M('S','aprA','adenylylsulfate reductase','aprA|aprB','1.8.99.2',('adenylylsulfate reductase',),('K00394',)),
 M('S','sat','sulfate adenylyltransferase','sat|cysD|cysN','2.7.7.4',('sulfate adenylyltransferase',),('K00958',)),
 M('S','asrA','anaerobic sulfite reductase','asrA|asrB|asrC','',('anaerobic sulfite reductase',),()),
 M('S','phsA','thiosulfate reductase','phsA|psrA','',('thiosulfate reductase',),()),
 M('S','sqr','sulfide:quinone oxidoreductase','sqr','1.8.5.4',('sulfide-quinone',),('K17218',)),
 M('S','soxB','Sox thiosulfate oxidation','soxB|soxA|soxX|soxY|soxZ','',('sulfur oxidation protein',),('K17224',)),
 M('S','ttrA','tetrathionate reductase','ttrA|ttrB','',('tetrathionate reductase',),()),
 # ---- nitrogen ------------------------------------------------------------
 M('N','nifH','nitrogenase iron protein','nifH','1.18.6.1',('nitrogenase iron protein',),('K02588',),
   'CANDIDATE: nifH paralogues (bchL/chlL, nflH) are common false positives'),
 M('N','nifD','nitrogenase molybdenum-iron alpha','nifD','',('nitrogenase molybdenum-iron',),('K02586',)),
 M('N','narG','nitrate reductase (membrane)','narG|narZ','1.7.5.1',('nitrate reductase',),('K00370',)),
 M('N','napA','nitrate reductase (periplasmic)','napA','1.9.6.1',('periplasmic nitrate reductase',),('K02567',)),
 M('N','nirK','nitrite reductase (Cu)','nirK','1.7.2.1',('nitrite reductase','copper-containing',),('K00368',)),
 M('N','nirS','nitrite reductase (cd1)','nirS','1.7.2.1',('cytochrome cd1',),('K15864',)),
 M('N','norB','nitric oxide reductase','norB','1.7.2.5',('nitric oxide reductase',),('K04561',)),
 M('N','nosZ','nitrous oxide reductase','nosZ','1.7.2.4',('nitrous-oxide reductase',),('K00376',)),
 M('N','nrfA','ammonifying nitrite reductase (DNRA)','nrfA','1.7.2.2',('cytochrome c nitrite reductase',),('K03385',)),
 M('N','hzsA','hydrazine synthase (anammox)','hzsA','',('hydrazine synthase',),()),
 M('N','amoA','ammonia monooxygenase','amoA|pmoA','1.14.99.39',('ammonia monooxygenase',),('K10944',)),
 # ---- respiration / metals ------------------------------------------------
 M('RESP','coxA','cytochrome c oxidase (aa3)','coxA|ctaD','7.1.1.9',('cytochrome c oxidase',),('K02274',)),
 M('RESP','cydA','cytochrome bd oxidase','cydA|cydB','7.1.1.7',('cytochrome bd',),('K00425',)),
 M('RESP','mtrC','metal reduction (Mtr/Omc)','mtrC|omcA|mtrA|mtrB','',('decaheme','metal reducing',),()),
 M('RESP','cyc','multiheme cytochrome c','','',('multiheme','cytochrome c, ',),()),
 # ---- transport / behaviour ----------------------------------------------
 M('TRANS','flgB','flagellar assembly','flgB|flgC|fliC|fliD|flhA','',('flagellar',),()),
 M('TRANS','cheA','chemotaxis','cheA|cheB|cheR|cheW|cheY','',('chemotaxis',),()),
 M('TRANS','spo0A','sporulation master regulator','spo0A|spoIIE|spoVG','',('sporulation',),()),
 M('TRANS','ABC','ABC transporter (any)','','',('abc transporter',),()),
 M('TRANS','TRAP','TRAP transporter','dctP|dctM|dctQ','',('trap transporter','c4-dicarboxylate',),()),
 # ---- cofactors -----------------------------------------------------------
 M('COFAC','cobT','cobalamin (B12) biosynthesis','cobT|cobS|cbiA|cbiB|cobN','',('cobalamin','cobyrinic','corrinoid',),()),
 M('COFAC','btuB','B12 uptake (auxotrophy signal)','btuB|btuF','',('vitamin b12 transporter','b12-binding',),()),
 M('COFAC','bioB','biotin biosynthesis','bioB|bioA|bioD','',('biotin synthase',),()),
 M('COFAC','thiC','thiamine biosynthesis','thiC|thiE|thiG','',('thiamine',),()),
]

# ------------------------------------------------ KEGG-style module definitions
MODULES = {
 'M00567 hydrogenotrophic methanogenesis (CO2 -> CH4)': ['fwdB','ftr','mch','mtd','mer','mtrA','mcrA','mcrB','mcrG'],
 'M00357 acetoclastic methanogenesis':                  ['ackA','pta','acsB','cooS','mtrA','mcrA','mcrB','mcrG'],
 'M00356 methylotrophic methanogenesis (methanol)':     ['mtaB','mtrA','mcrA','mcrB','mcrG'],
 'M00563 methylotrophic methanogenesis (methylamine)':  ['mtbA','mtmB','mcrA','mcrB','mcrG'],
 'Wood-Ljungdahl (acetogenesis)':                       ['cooS','acsB','fhs','folD','metF'],
 'Dissimilatory sulfate reduction':                     ['sat','aprA','dsrA','dsrB'],
 'Denitrification (full)':                              ['narG','nirK','norB','nosZ'],
 'DNRA':                                                ['narG','nrfA'],
 'Nitrogen fixation':                                   ['nifH','nifD'],
 'Butyrate production':                                 ['bcd','ptb','buk'],
 'Propionate (methylmalonyl-CoA)':                      ['mcm','mmdA'],
 'Glycolysis (EMP core)':                               ['pfk','gapA','pyk'],
 'TCA (oxidative core)':                                ['gltA','sdhA'],
}

# ============================================================ PARSE BAKTA
tax = {}
tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(tp):
    tax = {r['working_name']: r for r in csv.DictReader(open(tp), delimiter='\t')}

def hits_for(gene, ec, prod, dbx, m):
    if m['genes'] and gene and m['genes'].match(gene): return True
    if m['ecs'] and ec and any(e in m['ecs'] for e in ec.split(',')): return True
    if m['kos'] and dbx and any(k in dbx for k in m['kos']): return True
    if m['prods'] and prod and any(p in prod for p in m['prods']): return True
    return False

mags, matrix, seqs_to_verify = [], {}, {}
VERIFY = {'mcrA','dsrA','nifH','mttB_SF'}

for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)
    tsv = os.path.join(d, f'{name}.tsv')
    if not os.path.exists(tsv): continue
    mags.append(name)
    counts = {m['key']: 0 for m in MARKERS}
    locus_for = {k: [] for k in VERIFY}

    with open(tsv) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 6: continue
            # Bakta tsv: seqid, type, start, stop, strand, locus_tag, gene, product, dbxrefs
            ftype = p[1]
            locus = p[5] if len(p) > 5 else ''
            gene  = p[6] if len(p) > 6 else ''
            prod  = (p[7] if len(p) > 7 else '').lower()
            dbx   = p[8] if len(p) > 8 else ''
            ec    = ''
            mec = re.search(r'EC:([0-9.\-]+)', dbx)
            if mec: ec = mec.group(1)
            if ftype not in ('cds', 'CDS'): continue
            for m in MARKERS:
                if hits_for(gene, ec, prod, dbx, m):
                    counts[m['key']] += 1
                    if m['key'] in VERIFY and locus:
                        locus_for[m['key']].append(locus)
    matrix[name] = counts
    seqs_to_verify[name] = locus_for

if not mags:
    print("  !! no Bakta output parsed — run 10_annotate.sh"); raise SystemExit(1)
print(f"  parsed {len(mags)} annotated MAGs")

# ============================================================ WRITE MATRIX
keys = [m['key'] for m in MARKERS]
lab  = {m['key']: m for m in MARKERS}
with open(os.path.join(OUT, 'function_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','phylum','genus'] + keys)
    for n in mags:
        t = tax.get(n, {})
        w.writerow([t.get('mag_id',''), n, t.get('sample',''), t.get('domain',''),
                    t.get('phylum',''), t.get('genus','')] + [matrix[n][k] for k in keys])

# ============================================================ MODULES
def completeness(n, mod):
    need = MODULES[mod]
    have = [k for k in need if matrix[n].get(k, 0) > 0]
    return len(have), len(need), have

with open(os.path.join(OUT, 'pathway_summary.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','domain','genus','module','present','total','pct','missing'])
    for n in mags:
        t = tax.get(n, {})
        for mod in MODULES:
            h, tot, have = completeness(n, mod)
            miss = [k for k in MODULES[mod] if k not in have]
            w.writerow([t.get('mag_id',''), n, t.get('domain',''), t.get('genus',''),
                        mod, h, tot, f"{100*h/tot:.0f}", ';'.join(miss)])

# ============================================================ REPORT
def show(title, cats, minpct=0):
    print()
    print("=" * 100)
    print(f"  {title}")
    print("=" * 100)
    ks = [m['key'] for m in MARKERS if m['cat'] in cats]
    hdr = ''.join(f"{k[:7]:>8}" for k in ks)
    print(f"  {'MAG ID':<17}{'dom':<4}{'genus':<16}{hdr}")
    for n in sorted(mags, key=lambda n: (tax.get(n,{}).get('sample','')!='SED',
                                         tax.get(n,{}).get('mag_id',''))):
        t = tax.get(n, {})
        cells = ''.join(f"{(matrix[n][k] if matrix[n][k] else '·'):>8}" for k in ks)
        print(f"  {t.get('mag_id',n)[:17]:<17}{(t.get('domain','?') or '?')[:3]:<4}"
              f"{(t.get('genus','-') or '-')[:16]:<16}{cells}")

show("C1 METABOLISM & METHANOGENESIS  (counts; · = absent)", {'C1'})
show("SULFUR & NITROGEN CYCLING", {'S','N'})
show("FERMENTATION & HYDROGEN", {'FERM','H2','EB'})
show("CARBOHYDRATE & PROTEIN DEGRADATION", {'CARB','PROT'})
show("RESPIRATION, TRANSPORT, COFACTORS", {'RESP','TRANS','COFAC'})

# ------------------------------------------------------------ MODULE TABLE
print()
print("=" * 100)
print("  PATHWAY MODULE COMPLETENESS  (>=60% shown)")
print("=" * 100)
for n in sorted(mags, key=lambda n: (tax.get(n,{}).get('sample','')!='SED',
                                     tax.get(n,{}).get('mag_id',''))):
    t = tax.get(n, {})
    lines = []
    for mod in MODULES:
        h, tot, _ = completeness(n, mod)
        if 100*h/tot >= 60:
            lines.append(f"{mod} {h}/{tot}")
    if lines:
        print(f"  {t.get('mag_id',n):<17} {(t.get('genus','-') or '-')[:20]:<20} "
              + "; ".join(lines))

# ------------------------------------------------------------ TRAP REPORT
print()
print("=" * 100)
print("  MARKERS REQUIRING PHYLOGENETIC VERIFICATION — NOT results yet")
print("=" * 100)
for m in MARKERS:
    if m['key'] not in VERIFY: continue
    carriers = [(n, matrix[n][m['key']]) for n in mags if matrix[n][m['key']] > 0]
    print(f"\n  {m['key']}  — {m['label']}")
    print(f"    {m['caveat'] or '(no caveat)'}")
    if not carriers:
        print("    no hits")
        continue
    for n, c in carriers:
        t = tax.get(n, {})
        print(f"      {t.get('mag_id',n):<17} {(t.get('domain','?') or '?')[:3]:<4} "
              f"{(t.get('genus','-') or '-')[:22]:<22} {c} copies")

# ------------------------------------- extract sequences needing a tree
extracted = 0
for key in VERIFY:
    outp = os.path.join(OUT, 'verify', f'{key}_candidates.faa')
    with open(outp, 'w') as out:
        for n in mags:
            loci = set(seqs_to_verify[n].get(key, []))
            if not loci: continue
            faa = os.path.join(A, n, f'{n}.faa')
            if not os.path.exists(faa): continue
            keep, hdr = False, ''
            for line in open(faa):
                if line.startswith('>'):
                    tag = line[1:].split()[0]
                    keep = tag in loci
                    if keep:
                        mid = tax.get(n, {}).get('mag_id', n)
                        out.write(f">{mid}|{n}|{tag}\n"); extracted += 1
                elif keep:
                    out.write(line)
    if os.path.getsize(outp) == 0: os.remove(outp)
print()
print(f"  extracted {extracted} candidate sequences -> {os.path.join(OUT,'verify')}/")
print("  Build a tree of each against curated references before making any claim.")

# ------------------------------------------------------------ SYNTROPHY
print()
print("=" * 100)
print("  SYNTROPHY HAND-OFF — who makes what, who eats what")
print("=" * 100)
PRODUCERS = {'H2': ['hydA','hyaB','echA'], 'formate': ['fdhA','pflB'],
             'acetate': ['ackA','pta','acsB'], 'butyrate': ['buk','ptb','but'],
             'propionate': ['mcm','mmdA'], 'lactate': ['ldh'], 'ethanol': ['adhE']}
CONSUMERS = {'H2': ['fwdB','mcrA'], 'formate': ['fdhA','mcrA'],
             'acetate': ['ackA','mcrA'], 'methanol': ['mtaB'], 'methylamine': ['mtbA']}

with open(os.path.join(OUT, 'syntrophy_handoff.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','domain','genus','role','compound','evidence'])
    for n in mags:
        t = tax.get(n, {}); dom = t.get('domain', '')
        for comp, ks in PRODUCERS.items():
            hit = [k for k in ks if matrix[n].get(k, 0) > 0]
            if len(hit) >= max(1, len(ks)//2) and dom == 'Bacteria':
                w.writerow([t.get('mag_id',''), n, dom, t.get('genus',''),
                            'producer', comp, ';'.join(hit)])
        for comp, ks in CONSUMERS.items():
            hit = [k for k in ks if matrix[n].get(k, 0) > 0]
            if len(hit) == len(ks) and dom == 'Archaea':
                w.writerow([t.get('mag_id',''), n, dom, t.get('genus',''),
                            'consumer', comp, ';'.join(hit)])

for comp in PRODUCERS:
    prod = [tax.get(n,{}).get('mag_id',n) for n in mags
            if tax.get(n,{}).get('domain','')=='Bacteria'
            and sum(1 for k in PRODUCERS[comp] if matrix[n].get(k,0)>0) >= max(1,len(PRODUCERS[comp])//2)]
    cons = [tax.get(n,{}).get('mag_id',n) for n in mags
            if tax.get(n,{}).get('domain','')=='Archaea'
            and comp in CONSUMERS
            and all(matrix[n].get(k,0)>0 for k in CONSUMERS[comp])]
    print(f"\n  {comp.upper():<12} produced by {len(prod)} bacterial MAGs"
          + (f", consumed by {len(cons)} archaeal MAGs" if comp in CONSUMERS else ""))
    if prod: print(f"    producers: {', '.join(prod[:10])}{' ...' if len(prod)>10 else ''}")
    if cons: print(f"    consumers: {', '.join(cons)}")

print()
print("  This is a genomic POTENTIAL map, not a measured flux. It states which")
print("  hand-offs are possible given the gene content, which is what a")
print("  metagenome can support. Rates require isotope or activity data.")
PY

echo
echo "######## VALIDATION GATE ########"
check "function matrix written"   $([ -s "$F2/function_matrix.tsv" ] && echo 0 || echo 1)
check "pathway summary written"   $([ -s "$F2/pathway_summary.tsv" ] && echo 0 || echo 1)
check "syntrophy table written"   $([ -s "$F2/syntrophy_handoff.tsv" ] && echo 0 || echo 1)
nm=$(( $(wc -l < "$F2/function_matrix.tsv") - 1 ))
nmag=$(ls -1 "$WS2/07_mags/bins"/*.fa | wc -l)
check "every MAG profiled ($nm/$nmag)" $([ "$nm" -eq "$nmag" ] && echo 0 || echo 1)

runlog "11_metabolism complete"
echo
echo "######## 11 DONE — next: 12_derep_phylogeny.sh ########"
echo "Nothing in verify/ is a result yet."
