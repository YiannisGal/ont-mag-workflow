#!/usr/bin/env bash
# =============================================================================
# 19_auxotrophy.sh  --  who depends on whom
#
# THE ARGUMENT THIS ENABLES
#   Everything so far shows organisms CO-OCCUR. Co-occurrence is weak evidence
#   of interaction — organisms share a habitat for many reasons. Auxotrophy is
#   different: a genome that cannot synthesise an essential amino acid or
#   cofactor MUST obtain it from outside. In a closed enrichment culture, the
#   only outside is the other members. That converts "these organisms were
#   found together" into "these organisms could not grow apart", which is the
#   strongest form of the syntrophy argument and needs no new sequencing.
#
#   It also bears directly on the culturomics chapter: an isolate that grew in
#   pure culture cannot be auxotrophic for anything the medium lacked, so
#   predicted auxotrophies explain WHICH organisms resisted isolation.
#
# WHY KOfam AGAIN, NOT BAKTA
#   Step 11b established that Bakta names ~86% of archaeal proteins
#   "hypothetical" and invented 13 false mtaB calls. Biosynthesis pathways
#   scored on Bakta annotations would be wrong in exactly the same direction.
#   KOfam profiles with per-family adaptive thresholds are used instead, the
#   same method that recovered mcrA in 7 genomes where Bakta found none.
#
# THE HARD PART: ABSENCE vs INCOMPLETENESS
#   A missing gene in an 80%-complete MAG may be missing from the BIN, not the
#   organism. Every auxotrophy call here is therefore reported WITH the genome's
#   completeness, and a "confidence" column marks calls in genomes below 90%
#   complete as tentative. Auxotrophy inferred from an incomplete genome is a
#   hypothesis; inferred from a >95% complete genome it is close to a result.
#
# INPUT    09_annotation/<MAG>/<MAG>.faa, 07_mags/MAG_TAXONOMY.tsv, 11_kofam/
# OUTPUT   19_auxotrophy/
#
# RUNTIME  20-40 min (profile download + hmmsearch).
# USAGE    ws2 && bash scripts/19_auxotrophy.sh 2>&1 | tee logs/19_auxotrophy.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="19_auxotrophy"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"; X="$WS2/19_auxotrophy"
mkdir -p "$X/profiles" "$X/hits"

echo "######## 19 — AUXOTROPHY AND METABOLIC COMPLEMENTARITY ########"
echo "  started: $(date -Is)"
runlog "19_auxotrophy started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }

BASE="https://www.genome.jp/ftp/db/kofam"

# ---- pathway definitions: KO<TAB>gene<TAB>pathway --------------------------
cat > "$X/ko_pathways.tsv" <<'KOS'
K00765	hisG	His
K00013	hisD	His
K01693	hisB	His
K01657	trpE	Trp
K01696	trpB	Trp
K01609	trpC	Trp
K04518	pheA	Phe
K14170	pheA2	Phe
K00220	tyrC	Tyr
K04517	tyrA	Tyr
K00800	aroA	AAA_common
K01626	aroG	AAA_common
K01649	leuA	Leu
K00052	leuB	Leu
K01703	leuC	Leu
K01687	ilvD	Ile_Val
K00826	ilvE	Ile_Val
K00053	ilvC	Ile_Val
K01754	ilvA	Ile
K01714	dapA	Lys
K01586	lysA	Lys
K00674	dapD	Lys
K00872	thrB	Thr
K01733	thrC	Thr
K00003	hom	Thr_Met
K01739	metB	Met
K01760	metC	Met
K00549	metE	Met
K00548	metH	Met
K01738	cysK	Cys
K00640	cysE	Cys
K01940	argG	Arg
K01755	argH	Arg
K00930	argB	Arg
K00147	proA	Pro
K00286	proC	Pro
K00931	proB	Pro
K00058	serA	Ser
K01079	serB	Ser
K00831	serC	Ser
K00600	glyA	Gly
K00825	avtA	Ala
K14260	alaC	Ala
K00812	aspB	Asp
K01953	asnB	Asn
K00262	gdhA	Glu
K00265	gltB	Glu
K01915	glnA	Gln
K03147	thiC	B1_thiamine
K00788	thiE	B1_thiamine
K00941	thiD	B1_thiamine
K02858	ribB	B2_riboflavin
K00794	ribE	B2_riboflavin
K00793	ribH	B2_riboflavin
K00606	panB	B5_pantothenate
K01918	panC	B5_pantothenate
K00097	pdxA	B6_pyridoxal
K03474	pdxJ	B6_pyridoxal
K06215	pdxS	B6_pyridoxal
K01012	bioB	B7_biotin
K00652	bioF	B7_biotin
K00833	bioA	B7_biotin
K00796	folP	B9_folate
K01495	folE	B9_folate
K11754	folK	B9_folate
K02224	cbiA	B12_cobalamin
K00768	cobT	B12_cobalamin
K02230	cobN	B12_cobalamin
K19712	cbiB	B12_cobalamin
KOS
NKO=$(wc -l < "$X/ko_pathways.tsv")
echo "  KO families: $NKO across $(cut -f3 "$X/ko_pathways.tsv" | sort -u | wc -l) pathways"

# ---- thresholds (reuse the ko_list already downloaded at 11b) ---------------
if [ -s "$K/ko_list" ]; then
  cp -n "$K/ko_list" "$X/ko_list" 2>/dev/null || true
  echo "  reusing ko_list from step 11b"
else
  wget -q -O "$X/ko_list.gz" "$BASE/ko_list.gz" && gunzip -f "$X/ko_list.gz" \
    || echo "  !! no ko_list — falling back to E-value 1e-5"
fi

# ---- profiles ---------------------------------------------------------------
# ROUTE A: individual profiles. KEGG serves these intermittently — they worked
# at step 11b (78/78) and returned 404 for all 69 families here on 2026-08-23.
# Never rely on this endpoint alone.
got=0; miss=0
while IFS=$'\t' read -r ko gene path; do
  [ -z "$ko" ] && continue
  if [ -s "$X/profiles/${ko}.hmm" ]; then got=$((got+1)); continue; fi
  if [ -s "$K/profiles/${ko}.hmm" ]; then          # already fetched at 11b
    cp "$K/profiles/${ko}.hmm" "$X/profiles/"; got=$((got+1)); continue
  fi
  if wget -q -O "$X/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" \
     && head -1 "$X/profiles/${ko}.hmm" | grep -q HMMER; then
    got=$((got+1))
  else
    rm -f "$X/profiles/${ko}.hmm"; miss=$((miss+1))
  fi
done < "$X/ko_pathways.tsv"
echo "  route A (individual profiles): $got present, $miss missing"

# ROUTE B: the full archive, fetched once and shared with step 11b. Extract only
# the families needed. ~400 MB, resumable, reusable by any later KO analysis.
if [ "$got" -lt "$NKO" ]; then
  TARB="$K/profiles.tar.gz"
  mkdir -p "$K"
  if [ ! -s "$TARB" ]; then
    echo "  route B: fetching the KOfam profile archive once (~400 MB, resumable)"
    wget -c -O "$TARB" "$BASE/profiles.tar.gz" \
      || { echo "  !! archive download failed. Fetch $BASE/profiles.tar.gz by hand,"
           echo "     place it at $TARB, and re-run."; exit 1; }
  else
    echo "  route B: using cached archive $TARB"
  fi
  echo "  extracting the $NKO target families ..."
  awk -F'\t' '{print "profiles/"$1".hmm"}' "$X/ko_pathways.tsv" > "$X/wanted.txt"
  tar -xzf "$TARB" -C "$X" -T "$X/wanted.txt" 2>/dev/null \
    || tar -xzf "$TARB" -C "$X" $(tr '\n' ' ' < "$X/wanted.txt") 2>/dev/null || true
  # tar writes into $X/profiles/ — the same directory route A uses
  got=$(ls -1 "$X/profiles"/*.hmm 2>/dev/null | wc -l)
  echo "  after route B: $got profiles present"
  # keep them for step 11b and anything later
  cp -n "$X/profiles"/*.hmm "$K/profiles/" 2>/dev/null || true
fi

[ "$got" -gt 0 ] || { echo "  !! no profiles obtained by either route"; exit 1; }
if [ "$got" -lt $((NKO / 2)) ]; then
  echo "  !! fewer than half the families are available ($got/$NKO)."
  echo "     Pathway completeness would be systematically underestimated and every"
  echo "     genome would look auxotrophic. Refusing to continue."
  exit 1
fi
cat "$X/profiles"/*.hmm > "$X/pathways.hmm"

# ---- search -----------------------------------------------------------------
echo
echo "######## hmmsearch across all MAGs ($(date +%H:%M)) ########"
n=0; tot=$(ls -1d "$A"/*/ | wc -l)
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1)); out="$X/hits/${name}.tbl"
  if [ -s "$out" ]; then printf '\r  [%2d/%2d] %-28s cached' "$n" "$tot" "$name"; continue; fi
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 \
      "$X/pathways.hmm" "$faa" >/dev/null 2>&1 || echo "  !! failed on $name"
done
echo

# ---- analyse ----------------------------------------------------------------
python3 - "$X" "$F" "$A" <<'PY'
import sys, os, csv, glob
from collections import defaultdict, Counter

X, F, A = sys.argv[1], sys.argv[2], sys.argv[3]
COMPLETE_STRICT = 90.0     # above this, an absence is fairly trustworthy
PATH_THRESHOLD  = 0.50     # pathway present if >= this fraction of markers found

ko2path, ko2gene = {}, {}
for line in open(os.path.join(X, 'ko_pathways.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 3: ko2path[p[0]] = p[2]; ko2gene[p[0]] = p[1]
paths = sorted(set(ko2path.values()))
path_kos = defaultdict(list)
for ko, pth in ko2path.items(): path_kos[pth].append(ko)

thr = {}
kl = os.path.join(X, 'ko_list')
if os.path.exists(kl):
    for r in csv.DictReader(open(kl), delimiter='\t'):
        ko = r.get('knum') or r.get('#knum') or ''
        t = r.get('threshold', '')
        if ko in ko2path and t not in ('', '-'):
            try: thr[ko] = float(t)
            except ValueError: pass
print(f"  adaptive thresholds for {len(thr)}/{len(ko2path)} families")

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

hits = {}
for tbl in sorted(glob.glob(os.path.join(X, 'hits', '*.tbl'))):
    name = os.path.basename(tbl)[:-4]
    found = set()
    for line in open(tbl):
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) < 6: continue
        ko = f[2].split('.')[0]
        if ko not in ko2path: continue
        try: score, ev = float(f[5]), float(f[4])
        except ValueError: continue
        if (score >= thr[ko]) if ko in thr else (ev <= 1e-5):
            found.add(ko)
    hits[name] = found
mags = sorted(hits)
print(f"  MAGs analysed: {len(mags)}")

def frac(m, pth):
    ks = path_kos[pth]
    return sum(1 for k in ks if k in hits[m]) / len(ks)

# ------------------------------------------------------------- matrix
with open(os.path.join(X, 'pathway_completeness.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','genus','completeness'] + paths)
    for m in mags:
        t = tax.get(m, {})
        w.writerow([t.get('mag_id',''), m, t.get('sample',''), t.get('domain',''),
                    t.get('genus',''), t.get('completeness','')]
                   + [f"{frac(m,p_):.2f}" for p_ in paths])

AA = [p_ for p_ in paths if not p_.startswith('B')]
VIT = [p_ for p_ in paths if p_.startswith('B')]

def block(title, plist):
    print()
    print("=" * 116)
    print(f"  {title}   (+ = pathway present, · = ABSENT -> candidate auxotrophy)")
    print("=" * 116)
    print(f"  {'MAG ID':<17}{'comp%':>6}{'genus':<18}" + ''.join(f"{p_[:6]:>8}" for p_ in plist))
    for m in sorted(mags, key=lambda m: (tax.get(m,{}).get('sample','')!='SED',
                                         tax.get(m,{}).get('mag_id',''))):
        t = tax.get(m, {})
        cells = ''.join(f"{'+' if frac(m,p_) >= PATH_THRESHOLD else '·':>8}" for p_ in plist)
        print(f"  {t.get('mag_id',m)[:17]:<17}{(t.get('completeness','') or '')[:5]:>6}"
              f"{(t.get('genus','-') or '-')[:18]:<18}{cells}")

block("AMINO ACID BIOSYNTHESIS", AA)
block("VITAMIN / COFACTOR BIOSYNTHESIS", VIT)

# ------------------------------------------------------------- auxotrophies
print()
print("=" * 116)
print("  CANDIDATE AUXOTROPHIES — ranked by how much the genome can be trusted")
print("=" * 116)
rows = []
for m in mags:
    t = tax.get(m, {})
    try: comp = float(t.get('completeness', 0) or 0)
    except ValueError: comp = 0.0
    missing = [p_ for p_ in paths if frac(m, p_) < PATH_THRESHOLD]
    conf = 'firm' if comp >= COMPLETE_STRICT else 'tentative'
    rows.append((comp, m, t, missing, conf))
    for p_ in missing:
        pass
rows.sort(reverse=True)
print(f"  {'MAG ID':<17}{'comp%':>6}{'conf':<11}{'n missing':>10}   missing pathways")
for comp, m, t, missing, conf in rows:
    print(f"  {t.get('mag_id',m):<17}{comp:>6.1f}{conf:<11}{len(missing):>10}   "
          f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}")

with open(os.path.join(X, 'auxotrophies.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','completeness','confidence','pathway','fraction_present'])
    for comp, m, t, missing, conf in rows:
        for p_ in missing:
            w.writerow([t.get('mag_id',''), m, f"{comp:.1f}", conf, p_, f"{frac(m,p_):.2f}"])

# ------------------------------------------------------------- complementarity
print()
print("=" * 116)
print("  COMPLEMENTARITY — for each compound, who cannot make it and who can")
print("=" * 116)
net = []
for p_ in paths:
    aux = [m for m in mags if frac(m, p_) < PATH_THRESHOLD]
    pro = [m for m in mags if frac(m, p_) >= PATH_THRESHOLD]
    if not aux: continue
    print(f"\n  {p_}: {len(aux)} cannot synthesise | {len(pro)} can")
    firm = [m for m in aux
            if (float(tax.get(m, {}).get('completeness', 0) or 0) >= COMPLETE_STRICT)]
    if firm:
        print(f"    firm auxotrophs (>= {COMPLETE_STRICT:.0f}% complete): "
              + ', '.join(tax.get(m, {}).get('mag_id', m) for m in firm))
    net.append(dict(pathway=p_, n_auxotroph=len(aux), n_prototroph=len(pro),
                    firm_auxotrophs=';'.join(tax.get(m,{}).get('mag_id',m) for m in firm)))
with open(os.path.join(X, 'complementarity.tsv'), 'w', newline='') as fh:
    if net:
        w = csv.DictWriter(fh, fieldnames=list(net[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(net)

# ------------------------------------------------------------- dependency index
print()
print("=" * 116)
print("  DEPENDENCY INDEX — fraction of pathways each genome must import")
print("=" * 116)
idx = []
for m in mags:
    t = tax.get(m, {})
    try: comp = float(t.get('completeness', 0) or 0)
    except ValueError: comp = 0.0
    miss = sum(1 for p_ in paths if frac(m, p_) < PATH_THRESHOLD)
    idx.append((miss/len(paths), comp, t.get('mag_id', m), t.get('sample',''),
                t.get('genus',''), t.get('domain','')))
idx.sort(reverse=True)
print(f"  {'MAG ID':<17}{'smp':<5}{'dom':<5}{'genus':<20}{'dependency':>11}{'comp%':>7}")
for dep, comp, mag, smp, gen, dom in idx:
    print(f"  {mag:<17}{smp:<5}{(dom or '?')[:3]:<5}{(gen or '-')[:20]:<20}"
          f"{dep:>10.0%}{comp:>7.1f}")

bac = [i[0] for i in idx if i[5] == 'Bacteria']
arch = [i[0] for i in idx if i[5] == 'Archaea']
if bac and arch:
    print()
    print(f"  mean dependency — Bacteria {sum(bac)/len(bac):.0%} (n={len(bac)}) | "
          f"Archaea {sum(arch)/len(arch):.0%} (n={len(arch)})")
    print("  Methanogens are expected to be highly dependent: they are metabolic")
    print("  specialists that import most building blocks. A high value there is")
    print("  biology. A high value in a fermenter is more interesting.")

print()
print("  HOW TO USE THIS")
print("   - 'firm' calls (>=90% complete genomes) are the ones to build on.")
print("   - Auxotrophy predicts which organisms CANNOT be isolated in minimal")
print("     medium — testable against the culturomics chapter directly.")
print("   - Absence of a pathway is not proof: a divergent enzyme can escape an")
print("     HMM. Verify any single claim by inspecting the genes around the gap.")
PY

# =====================================================================
# BONUS: AMR genes — already computed by AMRFinderPlus inside Bakta
# =====================================================================
echo
echo "######## ANTIMICROBIAL RESISTANCE (from the Bakta expert-system pass) ########"
python3 - "$A" "$F" "$X" <<'PY'
import sys, os, csv, glob
from collections import Counter
A, F, X = sys.argv[1], sys.argv[2], sys.argv[3]
tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}
rows, per = [], Counter()
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); tsv = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(tsv): continue
    for line in open(tsv):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[1] not in ('cds','CDS'): continue
        dbx = f[8]
        if 'NCBIFAM' in dbx.upper() or 'AMR' in dbx.upper() or 'CARD' in dbx.upper():
            prod = f[7]
            if any(k in prod.lower() for k in
                   ('resistance','beta-lactamase','efflux','antibiotic','macrolide',
                    'tetracycline','aminoglycoside','vancomycin','chloramphenicol')):
                rows.append((tax.get(n,{}).get('mag_id',n), n, f[6], prod))
                per[tax.get(n,{}).get('mag_id',n)] += 1
if rows:
    print(f"  {len(rows)} AMR-associated genes across {len(per)} MAGs")
    for mag, c in per.most_common(12):
        print(f"    {mag:<17} {c}")
    with open(os.path.join(X, 'amr_genes.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t'); w.writerow(['mag_id','working_name','gene','product'])
        w.writerows(rows)
    print(f"  written: {os.path.join(X,'amr_genes.tsv')}")
    print("  NOTE: environmental AMR genes are usually intrinsic efflux pumps and")
    print("  housekeeping enzymes, not acquired resistance. Do not present this as")
    print("  a public-health finding without checking each hit.")
else:
    print("  no AMR-associated annotations found")
PY

echo
echo "######## VALIDATION GATE ########"
check "pathway completeness table" $([ -s "$X/pathway_completeness.tsv" ] && echo 0 || echo 1)
check "auxotrophy table"           $([ -s "$X/auxotrophies.tsv" ] && echo 0 || echo 1)
check "complementarity table"      $([ -s "$X/complementarity.tsv" ] && echo 0 || echo 1)

runlog "19_auxotrophy complete"
echo
echo "######## 19 DONE ########"
