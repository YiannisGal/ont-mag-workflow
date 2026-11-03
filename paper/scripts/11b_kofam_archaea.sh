#!/usr/bin/env bash
# =============================================================================
# 11b_kofam_archaea.sh  --  diagnose annotation bias, then fix it with KOfam
#
# WHY THIS EXISTS
#   Bakta's reference clusters are heavily bacteria-weighted. In this dataset:
#       Bacteria  9.3-48.4% hypothetical proteins (median ~24%)
#       Archaea   84.0-89.9% hypothetical
#   Gene CALLING was fine for the archaea (mean CDS 233-294 aa is normal); gene
#   NAMING largely failed. So a step-11 table showing few methanogenesis genes
#   in the eight methanogens would be an ANNOTATION ARTEFACT, not biology — and
#   is exactly the kind of empty result that reads as a finding.
#
#   This script (1) measures that bias objectively, and (2) if it is real,
#   re-annotates with KOfam HMMs, which are domain-agnostic.
#
# WHY KOfam AND NOT DRAM / eggNOG
#   KEGG serves INDIVIDUAL profiles at .../kofam/profiles/<KO>.hmm, so only the
#   ~70 families needed here are downloaded (~20 MB, not multi-GB). KOfam also
#   ships a per-family adaptive score THRESHOLD, so calls use the family's own
#   cutoff rather than one arbitrary e-value. Runs on plain HMMER 3.4.
#   No DRAM, no Nextflow, no 50 GB database.
#
# ALL 29 MAGs ARE SEARCHED, not only the archaea — the bacteria get the same
# KO-level evidence so the two domains stay directly comparable in one table.
#
# INPUT    09_annotation/<MAG>/<MAG>.faa
#          10_function/function_matrix.tsv, 09_annotation/annotation_summary.tsv
# OUTPUT   11_kofam/ko_matrix.tsv, ko_modules.tsv, ko_vs_bakta.tsv
#
# RUNTIME  download ~2 min, hmmsearch ~10-20 min.
# USAGE    ws2 && bash scripts/11b_kofam_archaea.sh 2>&1 | tee logs/11b_kofam.log
#          FORCE_KOFAM=1 bash scripts/11b_kofam_archaea.sh   # skip the gate
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="11b_kofam"
A="$WS2/09_annotation"; F="$WS2/07_mags"; FN="$WS2/10_function"; K="$WS2/11_kofam"
mkdir -p "$K/profiles" "$K/hits"

echo "######## 11b — ANNOTATION BIAS CHECK + KOfam RE-ANNOTATION ########"
echo "  started: $(date -Is)"
runlog "11b_kofam started"

# =====================================================================
# PHASE 1 — DIAGNOSE. Is the archaeal annotation actually deficient?
# =====================================================================
echo
echo "######## PHASE 1 — ANNOTATION COVERAGE BY DOMAIN ########"
NEED_KOFAM=$(python3 - "$FN" "$A" <<'PY'
import sys, os, csv, statistics as st
FN, A = sys.argv[1], sys.argv[2]

summ = {}
p = os.path.join(A, 'annotation_summary.tsv')
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        summ[r['working_name']] = r

mat = []
p = os.path.join(FN, 'function_matrix.tsv')
if os.path.exists(p):
    mat = list(csv.DictReader(open(p), delimiter='\t'))

C1 = ['mcrA','mcrB','mcrG','mtrA','mer','mtd','mch','ftr','fwdB',
      'mtaB','mtbA','mtmB','acsB','cooS']

print(f"  {'domain':<10} {'MAGs':>5} {'median hyp%':>12} {'median C1 markers':>19} "
      f"{'MAGs w/ >=1 C1':>15}")
verdict = 0
for dom in ('Bacteria', 'Archaea'):
    sub = [r for r in mat if r.get('domain') == dom]
    if not sub: continue
    hyps = []
    for r in sub:
        s = summ.get(r['working_name'], {})
        if s.get('pct_hypothetical'):
            hyps.append(float(s['pct_hypothetical']))
    c1n = [sum(1 for k in C1 if int(r.get(k, 0) or 0) > 0) for r in sub]
    withc1 = sum(1 for n in c1n if n > 0)
    print(f"  {dom:<10} {len(sub):>5} {(st.median(hyps) if hyps else float('nan')):>12.1f} "
          f"{(st.median(c1n) if c1n else 0):>19.1f} {withc1:>10}/{len(sub)}")
    if dom == 'Archaea':
        med_h = st.median(hyps) if hyps else 0
        med_c = st.median(c1n) if c1n else 0
        # a methanogen must carry mcrABG + several C1 steps. Fewer than 6 of 14
        # C1 markers, or >60% hypothetical, means the annotation failed.
        if med_h > 60 or med_c < 6:
            verdict = 1

print()
if verdict:
    print("  >> ARCHAEAL ANNOTATION IS DEFICIENT.")
    print("     The methanogenesis genes are present in these genomes (GTDB placed")
    print("     them as methanogens from ribosomal markers); Bakta simply could not")
    print("     name them. Proceeding to KOfam.")
else:
    print("  >> Archaeal annotation looks adequate. KOfam would add little.")
    print("     Set FORCE_KOFAM=1 to run anyway (still useful: KO numbers are")
    print("     needed for KEGG module completeness and for the manuscript).")
sys.exit(verdict)
PY
echo $?)

if [ "${NEED_KOFAM:-1}" -eq 0 ] && [ -z "${FORCE_KOFAM:-}" ]; then
  echo
  echo "  Gate not tripped and FORCE_KOFAM unset — stopping here."
  echo "  Re-run with FORCE_KOFAM=1 to annotate anyway."
  runlog "11b_kofam: gate not tripped, skipped"
  exit 0
fi

# =====================================================================
# PHASE 2 — FETCH ONLY THE PROFILES WE NEED
# =====================================================================
echo
echo "######## PHASE 2 — KOfam PROFILES ########"

# ---- the KO set. Curated; grouped by pathway. -----------------------------
cat > "$K/ko_targets.tsv" <<'KOS'
K00399	mcrA	methanogenesis	methyl-CoM reductase alpha (VERIFY: ACR look-alike)
K00401	mcrB	methanogenesis	methyl-CoM reductase beta
K00402	mcrG	methanogenesis	methyl-CoM reductase gamma
K00577	mtrA	methanogenesis	tetrahydromethanopterin S-methyltransferase A
K00578	mtrB	methanogenesis	tetrahydromethanopterin S-methyltransferase B
K00579	mtrC	methanogenesis	tetrahydromethanopterin S-methyltransferase C
K00580	mtrD	methanogenesis	tetrahydromethanopterin S-methyltransferase D
K00581	mtrE	methanogenesis	tetrahydromethanopterin S-methyltransferase E
K00584	mtrH	methanogenesis	tetrahydromethanopterin S-methyltransferase H
K00320	mer	hydrogenotrophic	methylene-H4MPT reductase
K00319	mtd	hydrogenotrophic	methylene-H4MPT dehydrogenase
K01499	mch	hydrogenotrophic	methenyl-H4MPT cyclohydrolase
K00672	ftr	hydrogenotrophic	formylmethanofuran-H4MPT formyltransferase
K00200	fwdA	hydrogenotrophic	formylmethanofuran dehydrogenase A
K00201	fwdB	hydrogenotrophic	formylmethanofuran dehydrogenase B
K00202	fwdC	hydrogenotrophic	formylmethanofuran dehydrogenase C
K00203	fwdD	hydrogenotrophic	formylmethanofuran dehydrogenase D
K11261	fwdF	hydrogenotrophic	formylmethanofuran dehydrogenase F
K00440	frhA	hydrogenotrophic	F420-reducing hydrogenase alpha
K00441	frhB	hydrogenotrophic	F420-reducing hydrogenase beta
K00443	frhG	hydrogenotrophic	F420-reducing hydrogenase gamma
K14126	mvhA	hydrogenotrophic	F420-non-reducing hydrogenase A
K14127	mvhD	hydrogenotrophic	F420-non-reducing hydrogenase D
K03388	hdrA2	energy	heterodisulfide reductase A2
K03389	hdrB2	energy	heterodisulfide reductase B2
K03390	hdrC2	energy	heterodisulfide reductase C2
K08264	hdrD	energy	heterodisulfide reductase D
K08265	hdrE	energy	heterodisulfide reductase E
K04480	mtaB	methylotrophic	methanol--corrinoid methyltransferase
K14080	mtaA	methylotrophic	methylated-thiol--CoM methyltransferase
K16176	mtaC	methylotrophic	methanol corrinoid protein
K14082	mtbA	methylotrophic	methylamine--CoM methyltransferase
K16157	mtmB	methylotrophic	monomethylamine methyltransferase
K16177	mtbB	methylotrophic	dimethylamine methyltransferase
K14083	mttB	TRAP	MttB superfamily — NOT trimethylamine MT; mostly bacterial betaine MT
K00925	ackA	acetoclastic	acetate kinase
K00625	pta	acetoclastic	phosphotransacetylase
K01895	acs	acetoclastic	acetyl-CoA synthetase
K00193	cdhC	acetoclastic	CO dehydrogenase/acetyl-CoA synthase beta
K00194	cdhD	acetoclastic	CO dehydrogenase/acetyl-CoA synthase delta
K00197	cdhE	acetoclastic	CO dehydrogenase/acetyl-CoA synthase gamma
K00198	cooS	wood-ljungdahl	anaerobic carbon monoxide dehydrogenase
K14138	acsB	wood-ljungdahl	acetyl-CoA synthase
K01938	fhs	wood-ljungdahl	formate--tetrahydrofolate ligase
K01491	folD	wood-ljungdahl	methylene-THF dehydrogenase/cyclohydrolase
K00297	metF	wood-ljungdahl	methylene-THF reductase
K11180	dsrA	sulfur	dissimilatory sulfite reductase alpha (VERIFY: reductive vs oxidative)
K11181	dsrB	sulfur	dissimilatory sulfite reductase beta
K00394	aprA	sulfur	adenylylsulfate reductase A
K00395	aprB	sulfur	adenylylsulfate reductase B
K00958	sat	sulfur	sulfate adenylyltransferase
K17218	sqr	sulfur	sulfide:quinone oxidoreductase
K17224	soxB	sulfur	sulfur oxidation protein SoxB
K02588	nifH	nitrogen	nitrogenase iron protein (VERIFY: bchL/chlL paralogues)
K02586	nifD	nitrogen	nitrogenase molybdenum-iron alpha
K02591	nifK	nitrogen	nitrogenase molybdenum-iron beta
K00370	narG	nitrogen	nitrate reductase alpha
K02567	napA	nitrogen	periplasmic nitrate reductase
K00368	nirK	nitrogen	nitrite reductase (NO-forming, Cu)
K15864	nirS	nitrogen	nitrite reductase (NO-forming, cd1)
K04561	norB	nitrogen	nitric oxide reductase B
K00376	nosZ	nitrogen	nitrous-oxide reductase
K03385	nrfA	nitrogen	nitrite reductase (ammonia-forming)
K00169	porA	fermentation	pyruvate:ferredoxin oxidoreductase alpha
K00656	pflB	fermentation	formate C-acetyltransferase
K00016	ldh	fermentation	L-lactate dehydrogenase
K04072	adhE	fermentation	acetaldehyde/alcohol dehydrogenase
K00929	buk	fermentation	butyrate kinase
K00634	ptb	fermentation	phosphate butyryltransferase
K00248	bcd	fermentation	butyryl-CoA dehydrogenase
K01034	atoD	fermentation	acetate CoA-transferase alpha
K01847	mutA	fermentation	methylmalonyl-CoA mutase
K00123	fdhF	hydrogen	formate dehydrogenase major subunit
K00532	hydA	hydrogen	ferredoxin hydrogenase
K06281	hyaB	hydrogen	[NiFe] hydrogenase large subunit
K03615	rnfC	energy	Rnf electron transport complex C
K03612	rnfA	energy	Rnf electron transport complex A
K03522	etfA	energy	electron transfer flavoprotein alpha
KOS
NKO=$(wc -l < "$K/ko_targets.tsv")
echo "  target KO families: $NKO"

BASE="https://www.genome.jp/ftp/db/kofam"

# ---- ko_list gives each family its own adaptive score threshold ----------
if [ ! -s "$K/ko_list" ]; then
  echo "  fetching ko_list (adaptive thresholds) ..."
  wget -q -O "$K/ko_list.gz" "$BASE/ko_list.gz" \
    && gunzip -f "$K/ko_list.gz" \
    || echo "  !! could not fetch ko_list — will fall back to E-value 1e-10"
fi
[ -s "$K/ko_list" ] && echo "  ko_list: $(wc -l < "$K/ko_list") families" \
                    || echo "  ko_list: ABSENT (E-value fallback)"

# ---- ROUTE A: individual profiles (cheap; ~20 MB) ------------------------
miss=0; got=0; first_err=""
while IFS=$'\t' read -r ko sym cat desc; do
  [ -z "$ko" ] && continue
  if [ -s "$K/profiles/${ko}.hmm" ]; then got=$((got+1)); continue; fi
  err=$(wget -O "$K/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" 2>&1)
  if [ -s "$K/profiles/${ko}.hmm" ] && head -1 "$K/profiles/${ko}.hmm" | grep -q HMMER; then
    got=$((got+1))
  else
    rm -f "$K/profiles/${ko}.hmm"; miss=$((miss+1))
    [ -z "$first_err" ] && first_err="$err"
  fi
done < "$K/ko_targets.tsv"
echo "  route A (individual profiles): $got present, $miss missing"

# ---- ROUTE B: whole archive, extract only what we need --------------------
# KEGG does not always expose profiles/ as a browsable directory. When route A
# fails, fetch profiles.tar.gz once (~400 MB) and pull out just our families.
if [ "$got" -eq 0 ]; then
  echo
  echo "  route A failed. First wget error was:"
  echo "$first_err" | tail -4 | sed 's/^/      /'
  echo
  echo "  route B: downloading the full KOfam profile archive (once, ~400 MB)"
  TARB="$K/profiles.tar.gz"
  if [ ! -s "$TARB" ]; then
    wget -c -O "$TARB" "$BASE/profiles.tar.gz" \
      || { echo "  !! archive download failed — no network route to genome.jp."
           echo "     Fetch $BASE/profiles.tar.gz manually and place it at:"
           echo "     $TARB   then re-run this script."; exit 1; }
  fi
  echo "  extracting the $NKO target families ..."
  awk -F'\t' '{print "profiles/"$1".hmm"}' "$K/ko_targets.tsv" > "$K/wanted.txt"
  tar -xzf "$TARB" -C "$K" --wildcards --no-anchored \
      -T "$K/wanted.txt" 2>/dev/null || \
  tar -xzf "$TARB" -C "$K" $(tr '\n' ' ' < "$K/wanted.txt") 2>/dev/null || true
  got=$(ls -1 "$K/profiles"/*.hmm 2>/dev/null | wc -l)
  echo "  route B recovered: $got profiles"
fi

[ "$got" -gt 0 ] || { echo "  !! no profiles obtained by either route"; exit 1; }

cat "$K/profiles"/*.hmm > "$K/kofam_targets.hmm"
echo "  combined HMM: $(du -h "$K/kofam_targets.hmm" | cut -f1)"
record_version "hmmsearch" "${ENV_hmmsearch:--}" hmmsearch -h

# =====================================================================
# PHASE 3 — SEARCH EVERY PROTEOME (both domains)
# =====================================================================
echo
echo "######## PHASE 3 — hmmsearch, all 29 MAGs ($(date +%H:%M)) ########"
n=0; tot=$(ls -1d "$A"/*/ 2>/dev/null | wc -l)
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1))
  out="$K/hits/${name}.tbl"
  if [ -s "$out" ]; then printf '\r  [%2d/%2d] %-28s cached' "$n" "$tot" "$name"; continue; fi
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 \
      "$K/kofam_targets.hmm" "$faa" > /dev/null 2>&1 \
    || echo "  !! hmmsearch failed on $name"
done
echo

# =====================================================================
# PHASE 4 — KO MATRIX, MODULES, AND COMPARISON WITH BAKTA
# =====================================================================
python3 - "$K" "$F" "$FN" <<'PY'
import sys, os, csv, glob, statistics as st

K, F, FN = sys.argv[1], sys.argv[2], sys.argv[3]

targets = {}
for line in open(os.path.join(K, 'ko_targets.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 4: targets[p[0]] = dict(sym=p[1], cat=p[2], desc=p[3])

# KOfam adaptive thresholds
thr = {}
kl = os.path.join(K, 'ko_list')
if os.path.exists(kl):
    for r in csv.DictReader(open(kl), delimiter='\t'):
        ko = r.get('knum') or r.get('#knum') or ''
        t = r.get('threshold', '')
        if ko in targets and t not in ('', '-'):
            try: thr[ko] = (float(t), r.get('score_type', 'full'))
            except ValueError: pass
print(f"  adaptive thresholds available for {len(thr)}/{len(targets)} families")
print("  (families without one fall back to E-value 1e-5 — flagged in the output)")

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

hits, loci = {}, {}
for tbl in sorted(glob.glob(os.path.join(K, 'hits', '*.tbl'))):
    name = os.path.basename(tbl)[:-4]
    counts, lc = {}, {}
    for line in open(tbl):
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) < 8: continue
        target, query = f[0], f[2]
        ko = query.split('.')[0]
        if ko not in targets: continue
        try: fullscore, evalue = float(f[5]), float(f[4])
        except ValueError: continue
        ok = (fullscore >= thr[ko][0]) if ko in thr else (evalue <= 1e-5)
        if ok:
            counts[ko] = counts.get(ko, 0) + 1
            lc.setdefault(ko, []).append(target)
    hits[name] = counts; loci[name] = lc

mags = sorted(hits)
kos  = list(targets)

with open(os.path.join(K, 'ko_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','genus']
               + [f"{k}:{targets[k]['sym']}" for k in kos])
    for m in mags:
        t = tax.get(m, {})
        w.writerow([t.get('mag_id',''), m, t.get('sample',''), t.get('domain',''),
                    t.get('genus','')] + [hits[m].get(k, 0) for k in kos])

# ------------------------------------------------------------- MODULES
MODULES = {
 'M00567 hydrogenotrophic (CO2->CH4)': ['K00200','K00201','K00672','K01499','K00319','K00320','K00577','K00399','K00401','K00402'],
 'M00357 acetoclastic':                ['K00925','K00625','K00193','K00197','K00577','K00399','K00401','K00402'],
 'M00356 methylotrophic (methanol)':   ['K04480','K16176','K14080','K00399','K00401','K00402'],
 'M00563 methylotrophic (methylamine)':['K16157','K14082','K00399','K00401','K00402'],
 'F420 / Mvh-Hdr energy conservation': ['K00440','K00441','K00443','K14126','K14127','K03388','K03389','K03390'],
 'Wood-Ljungdahl':                     ['K00198','K14138','K01938','K01491','K00297'],
 'Dissimilatory sulfate reduction':    ['K00958','K00394','K00395','K11180','K11181'],
 'Denitrification (full)':             ['K00370','K00368','K04561','K00376'],
 'Nitrogen fixation':                  ['K02588','K02586','K02591'],
 'Butyrate production':                ['K00248','K00634','K00929'],
}
with open(os.path.join(K, 'ko_modules.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','domain','genus','module','present','total','pct','missing'])
    for m in mags:
        t = tax.get(m, {})
        for mod, ks in MODULES.items():
            have = [k for k in ks if hits[m].get(k, 0) > 0]
            miss = [f"{k}:{targets[k]['sym']}" for k in ks if k not in have]
            w.writerow([t.get('mag_id',''), m, t.get('domain',''), t.get('genus',''),
                        mod, len(have), len(ks),
                        f"{100*len(have)/len(ks):.0f}", ';'.join(miss)])

# ------------------------------------------------------------- REPORT
def block(title, cats):
    ks = [k for k in kos if targets[k]['cat'] in cats]
    if not ks: return
    print()
    print("=" * 108)
    print(f"  {title}")
    print("=" * 108)
    print(f"  {'MAG ID':<17}{'dom':<4}{'genus':<18}" + ''.join(f"{targets[k]['sym'][:6]:>7}" for k in ks))
    for m in sorted(mags, key=lambda m: (tax.get(m,{}).get('domain','')!='Archaea',
                                         tax.get(m,{}).get('mag_id',''))):
        t = tax.get(m, {})
        cells = ''.join(f"{(hits[m].get(k,0) or '·'):>7}" for k in ks)
        print(f"  {t.get('mag_id',m)[:17]:<17}{(t.get('domain','?') or '?')[:3]:<4}"
              f"{(t.get('genus','-') or '-')[:18]:<18}{cells}")

block("METHANOGENESIS CORE (mcr / mtr)", {'methanogenesis'})
block("HYDROGENOTROPHIC PATHWAY + F420/Mvh HYDROGENASES", {'hydrogenotrophic'})
block("METHYLOTROPHIC + ACETOCLASTIC BRANCHES", {'methylotrophic','acetoclastic','TRAP'})
block("ENERGY CONSERVATION (Hdr / Rnf / Etf)", {'energy'})
block("SULFUR + NITROGEN", {'sulfur','nitrogen'})
block("FERMENTATION + HYDROGEN + WOOD-LJUNGDAHL", {'fermentation','hydrogen','wood-ljungdahl'})

print()
print("=" * 108)
print("  MODULE COMPLETENESS (>=50%)")
print("=" * 108)
for m in sorted(mags, key=lambda m: (tax.get(m,{}).get('domain','')!='Archaea',
                                     tax.get(m,{}).get('mag_id',''))):
    t = tax.get(m, {})
    out = []
    for mod, ks in MODULES.items():
        have = sum(1 for k in ks if hits[m].get(k, 0) > 0)
        if 100*have/len(ks) >= 50: out.append(f"{mod} {have}/{len(ks)}")
    if out:
        print(f"  {t.get('mag_id',m):<17} {(t.get('genus','-') or '-')[:20]:<20} " + "; ".join(out))

# ------------------------------------------------------------- BAKTA vs KOfam
print()
print("=" * 108)
print("  WHAT KOfam RECOVERED THAT BAKTA MISSED")
print("=" * 108)
bak = {}
p = os.path.join(FN, 'function_matrix.tsv')
if os.path.exists(p):
    bak = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}
PAIRS = [('K00399','mcrA'),('K00401','mcrB'),('K00402','mcrG'),('K00577','mtrA'),
         ('K00320','mer'),('K00319','mtd'),('K01499','mch'),('K00672','ftr'),
         ('K00201','fwdB'),('K04480','mtaB'),('K14082','mtbA'),('K11180','dsrA'),
         ('K02588','nifH')]
print(f"  {'marker':<8} {'Bakta MAGs':>11} {'KOfam MAGs':>11} {'recovered':>10}")
for ko, sym in PAIRS:
    nb = sum(1 for m in mags if int((bak.get(m, {}) or {}).get(sym, 0) or 0) > 0)
    nk = sum(1 for m in mags if hits[m].get(ko, 0) > 0)
    print(f"  {sym:<8} {nb:>11} {nk:>11} {nk-nb:>+10}")

with open(os.path.join(K, 'ko_vs_bakta.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(['marker','ko','n_bakta','n_kofam','delta'])
    for ko, sym in PAIRS:
        nb = sum(1 for m in mags if int((bak.get(m, {}) or {}).get(sym, 0) or 0) > 0)
        nk = sum(1 for m in mags if hits[m].get(ko, 0) > 0)
        w.writerow([sym, ko, nb, nk, nk-nb])

# ------------------------------------------------------------- TRAPS
print()
print("=" * 108)
print("  STILL REQUIRES PHYLOGENETIC VERIFICATION")
print("=" * 108)
os.makedirs(os.path.join(K, 'verify'), exist_ok=True)
for ko, note in (('K00399','mcrA — alkyl-CoM reductase (ACR) of alkane oxidisers is identical by HMM'),
                 ('K11180','dsrA — reductive (SRB) vs oxidative type not separable by HMM'),
                 ('K02588','nifH — bchL/chlL/nflH paralogues score as nitrogenase'),
                 ('K14083','mttB — MttB SUPERFAMILY, predominantly bacterial glycine betaine MT; NOT TMA methyltransferase, NEVER count as methylamine methanogenesis')):
    carriers = [(m, hits[m].get(ko, 0)) for m in mags if hits[m].get(ko, 0) > 0]
    sym = targets.get(ko, {}).get('sym', ko)
    print(f"\n  {ko} {sym}")
    print(f"    {note}")
    if not carriers: print("    no hits"); continue
    for m, c in carriers:
        t = tax.get(m, {})
        print(f"      {t.get('mag_id',m):<17} {(t.get('domain','?') or '?')[:3]:<4} "
              f"{(t.get('genus','-') or '-')[:22]:<22} {c} copies")
    # extract sequences for tree-building
    outp = os.path.join(K, 'verify', f'{sym}_{ko}_candidates.faa')
    with open(outp, 'w') as out:
        for m, _ in carriers:
            want = set(loci[m].get(ko, []))
            faa = os.path.join(os.path.dirname(F), '09_annotation', m, f'{m}.faa')
            if not os.path.exists(faa): continue
            keep = False
            for line in open(faa):
                if line.startswith('>'):
                    tag = line[1:].split()[0]; keep = tag in want
                    if keep:
                        out.write(f">{tax.get(m,{}).get('mag_id',m)}|{m}|{tag}\n")
                elif keep: out.write(line)
    if os.path.getsize(outp) == 0: os.remove(outp)
    else: print(f"    -> sequences: {outp}")

print()
print("  A complete mcrABG set here is a CANDIDATE methanogen, not a confirmed one,")
print("  until the mcrA tree separates true MCR from ACR. workspace_1 had ~50%")
print("  false positives on exactly this call.")
PY

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "KO matrix written"     $([ -s "$K/ko_matrix.tsv" ] && echo 0 || echo 1)
check "KO modules written"    $([ -s "$K/ko_modules.tsv" ] && echo 0 || echo 1)
nh=$(ls -1 "$K/hits"/*.tbl 2>/dev/null | wc -l)
nm=$(ls -1 "$WS2/07_mags/bins"/*.fa | wc -l)
check "every MAG searched ($nh/$nm)" $([ "$nh" -eq "$nm" ] && echo 0 || echo 1)
# the whole point: methanogens must now carry mcrA
narc=$(awk -F'\t' 'NR>1 && $4=="Archaea"' "$K/ko_matrix.tsv" | wc -l)
mcr=$(python3 - "$K" <<'PY'
import csv,sys,os
p=os.path.join(sys.argv[1],'ko_matrix.tsv')
rows=list(csv.DictReader(open(p),delimiter='\t'))
col=[c for c in rows[0] if c.startswith('K00399')][0] if rows else None
print(sum(1 for r in rows if r.get('domain')=='Archaea' and int(r.get(col,0) or 0)>0) if col else 0)
PY
)
check "archaeal MAGs with mcrA ($mcr/$narc)" $([ "${mcr:-0}" -gt 0 ] && echo 0 || echo 1)

runlog "11b_kofam complete"
echo
echo "######## 11b DONE ########"
