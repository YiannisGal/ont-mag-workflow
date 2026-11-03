#!/usr/bin/env bash
# =============================================================================
# 19b_auxotrophy_fixed.sh  --  auxotrophy, done properly
#
# WHAT WENT WRONG IN STEP 19
#   29/29 genomes "auxotrophic" for vitamin B6. 26/29 for alanine. Neither is
#   possible: a community in which no member can make pyridoxal does not exist,
#   and alanine is made by transamination through several redundant routes.
#   Those were marker-set failures reported as biology.
#
#   THE ROOT CAUSE WAS STRUCTURAL, not a matter of a few missing KOs:
#
#   1. ONE ROUTE ASSUMED PER COMPOUND. B6 was scored on pdxA/pdxJ/pdxD — the
#      DXP-DEPENDENT route, which is essentially confined to gammaproteobacteria.
#      Most bacteria, and all archaea, use the DXP-INDEPENDENT route (pdxS/pdxT),
#      which was not in the marker set at all. Every methanogen in this dataset
#      was therefore guaranteed to score zero for a pathway it certainly has.
#
#   2. FRACTION-OF-UNION SCORING. Step 19 called a pathway present when >=50% of
#      ALL its markers were found. When markers belong to ALTERNATIVE routes that
#      is mathematically wrong: an organism with a complete route B scores 0% if
#      the set is dominated by route A. Completeness must be computed PER ROUTE
#      and the best route taken.
#
#   3. TOO FEW MARKERS. With one marker per pathway, a 90%-complete genome misses
#      a present pathway 10% of the time; with three independent markers, 0.1%.
#      The script now prints that table and refuses to score any pathway with
#      fewer than 3 surviving markers, calling it UNSCORABLE instead.
#
# THE THREE GUARDS
#   A. LABEL GUARD    every KO is checked against KEGG's own definition from
#                     ko_list. A mismatch DROPS that marker with a warning rather
#                     than silently mis-scoring (step 20's K14080 error).
#   B. POSITIVE CONTROL  rpoB and tuf must be found in nearly every genome. If a
#                     universal gene is missing, the search is broken and no
#                     auxotrophy is reported at all.
#   C. PREVALENCE FLAG  any pathway absent in >80% of genomes is reported as a
#                     SUSPECTED MARKER-SET FAILURE, not as an auxotrophy. This is
#                     the rule that would have caught B6 automatically.
#
# HONEST LIMIT
#   Some compounds are not scorable from marker presence at all. Alanine is the
#   clearest case: it is produced by broad-specificity transaminases that also do
#   other jobs, so both presence and absence are uninformative. If Ala still
#   trips the prevalence flag after the expanded route set, it is reported as
#   UNSCORABLE and excluded. Reporting "we could not score this" is a result;
#   reporting a false auxotrophy is not.
#
# INPUT    09_annotation/<MAG>/<MAG>.faa, 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   19b_auxotrophy/
#     pathway_completeness.tsv   best-route completeness per genome per pathway
#     auxotrophies.tsv           firm auxotrophies only
#     marker_audit.tsv           every KO: kept, dropped, why
#     UNSCORABLE.txt             pathways this method cannot decide, and why
#
# RUNTIME  30-60 min.
# USAGE    ws2 && bash scripts/19b_auxotrophy_fixed.sh 2>&1 | tee logs/19b_aux.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="19b_auxotrophy"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"
X="$WS2/19_auxotrophy"; Y="$WS2/19b_auxotrophy"
mkdir -p "$Y/profiles" "$Y/hits"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 19b — AUXOTROPHY, ROUTE-AWARE ########"
echo "  started: $(date -Is)"
runlog "19b_auxotrophy started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }

# ---------------------------------------------------------------------------
# KO <tab> gene <tab> PATHWAY <tab> ROUTE <tab> KEGG-definition keyword
# ROUTE groups markers that must be present TOGETHER. Completeness is computed
# per route and the BEST route wins — the fix for failure mode 2 above.
# CONTROL rows are universal genes that prove the search works.
# ---------------------------------------------------------------------------
cat > "$Y/routes.tsv" <<'KOS'
K03043	rpoB	CONTROL	bacteria	DNA-directed RNA polymerase subunit beta
K02358	tuf	CONTROL	bacteria	elongation factor Tu
K03231	ef1a	CONTROL	archaea	elongation factor 1-alpha
K01915	glnA	CONTROL	both	glutamine synthetase
K02886	rplB	CONTROL	both	ribosomal protein L2
K02982	rpsC	CONTROL	both	ribosomal protein S3
K02931	rplE	CONTROL	both	ribosomal protein L5
K06215	pdxS	B6_pyridoxal	DXP_independent	pyridoxal 5'-phosphate synthase
K08681	pdxT	B6_pyridoxal	DXP_independent	phosphate synthase pdxt
K00097	pdxA	B6_pyridoxal	DXP_dependent	4-hydroxythreonine-4-phosphate dehydrogenase
K03474	pdxJ	B6_pyridoxal	DXP_dependent	pyridoxine 5-phosphate synthase
K00275	pdxH	B6_pyridoxal	DXP_dependent	pyridoxamine 5'-phosphate oxidase
K01662	dxs	B6_pyridoxal	DXP_dependent	1-deoxy-D-xylulose-5-phosphate synthase
K00868	pdxK	B6_pyridoxal	salvage	pyridoxine kinase
K14260	alaC	Ala	transamination	alanine-synthesizing
K14261	alaA	Ala	transamination	alanine-synthesizing
K00259	ald	Ala	dehydrogenase	alanine dehydrogenase
K01775	alr	Ala	racemase_context	alanine racemase
K04517	tyrA	Tyr	prephenate_DH	prephenate dehydrogenase
K00220	tyrC	Tyr	prephenate_DH	prephenate dehydrogenase
K00832	tyrB	Tyr	transamination	aromatic-amino-acid transaminase
K00812	aspB	Tyr	transamination	aspartate aminotransferase
K01850	aroQ	Tyr	chorismate_mutase	chorismate mutase
K04092	aroH	Tyr	chorismate_mutase	chorismate mutase
K00765	hisG	His	main	ATP phosphoribosyltransferase
K00013	hisD	His	main	histidinol dehydrogenase
K01693	hisB	His	main	imidazoleglycerol-phosphate dehydratase
K01523	hisIE	His	main	phosphoribosyl-atp
K01657	trpE	Trp	main	anthranilate synthase
K01696	trpB	Trp	main	tryptophan synthase
K01609	trpC	Trp	main	indole-3-glycerol phosphate synthase
K00766	trpD	Trp	main	anthranilate phosphoribosyltransferase
K01649	leuA	Leu	main	2-isopropylmalate synthase
K00052	leuB	Leu	main	3-isopropylmalate dehydrogenase
K01703	leuC	Leu	main	3-isopropylmalate
K00826	ilvE	Leu	main	branched-chain amino acid aminotransferase
K01687	ilvD	Ile_Val	main	dihydroxy-acid dehydratase
K00053	ilvC	Ile_Val	main	ketol-acid reductoisomerase
K01652	ilvB	Ile_Val	main	acetolactate synthase
K00826	ilvE2	Ile_Val	main	branched-chain amino acid aminotransferase
K01714	dapA	Lys	DAP	4-hydroxy-tetrahydrodipicolinate synthase
K01586	lysA	Lys	DAP	diaminopimelate decarboxylase
K00674	dapD	Lys	DAP	succinyltransferase
K01778	dapF	Lys	DAP	diaminopimelate epimerase
K00872	thrB	Thr	main	homoserine kinase
K01733	thrC	Thr	main	threonine synthase
K00003	hom	Thr	main	homoserine dehydrogenase
K01739	metB	Met	transsulfuration	cystathionine gamma-synthase
K01760	metC	Met	transsulfuration	beta-lyase
K00549	metE	Met	methyl_transfer	5-methyltetrahydropteroyltriglutamate
K00548	metH	Met	methyl_transfer	homocysteine methyltransferase
K01738	cysK	Cys	main	cysteine synthase
K00640	cysE	Cys	main	serine O-acetyltransferase
K10150	cysM	Cys	main	cysteine synthase
K01940	argG	Arg	main	argininosuccinate synthase
K01755	argH	Arg	main	argininosuccinate lyase
K00930	argB	Arg	main	acetylglutamate kinase
K00145	argC	Arg	main	N-acetyl-gamma-glutamyl-phosphate reductase
K00147	proA	Pro	main	glutamate-5-semialdehyde dehydrogenase
K00286	proC	Pro	main	pyrroline-5-carboxylate reductase
K00931	proB	Pro	main	glutamate 5-kinase
K00058	serA	Ser	phosphorylated	D-3-phosphoglycerate dehydrogenase
K01079	serB	Ser	phosphorylated	phosphoserine phosphatase
K00831	serC	Ser	phosphorylated	phosphoserine aminotransferase
K00600	glyA	Gly	main	glycine hydroxymethyltransferase
K00281	gcvP	Gly	main	glycine dehydrogenase
K00605	gcvT	Gly	main	aminomethyltransferase
K02437	gcvH	Gly	main	glycine cleavage system H
K03147	thiC	B1_thiamine	main	phosphomethylpyrimidine synthase
K00788	thiE	B1_thiamine	main	thiamine-phosphate
K00941	thiD	B1_thiamine	main	phosphomethylpyrimidine kinase
K03149	thiG	B1_thiamine	main	thiazole synthase
K02858	ribB	B2_riboflavin	main	butanone 4-phosphate synthase
K00794	ribE	B2_riboflavin	main	ribityllumazine synthase
K00793	ribH	B2_riboflavin	main	riboflavin synthase
K11752	ribD	B2_riboflavin	main	deaminase
K00606	panB	B5_pantothenate	main	hydroxymethyltransferase
K01918	panC	B5_pantothenate	main	pantoate
K00077	panE	B5_pantothenate	main	dehydropantoate 2-reductase
K01579	panD	B5_pantothenate	main	aspartate 1-decarboxylase
K00652	bioF	B7_biotin	main	8-amino-7-oxononanoate synthase
K00833	bioA	B7_biotin	main	8-amino-7-oxononanoate aminotransferase
K01935	bioD	B7_biotin	main	dethiobiotin synthetase
K01012	bioB	B7_biotin	main	biotin synthase
K01633	folB	B9_folate	main	dihydroneopterin aldolase
K00796	folP	B9_folate	main	dihydropteroate synthase
K13939	folK	B9_folate	main	hydroxymethyldihydropteridine
K00287	folA	B9_folate	main	dihydrofolate reductase
K00768	cobT	B12_cobalamin	main	dimethylbenzimidazole
K02230	cobN	B12_cobalamin	main	cobaltochelatase
K19712	cobU	B12_cobalamin	main	adenosylcobinamide-phosphate guanylyltransferase
K02224	cbiA	B12_cobalamin	main	cobyrinic acid
KOS
NKO=$(awk -F'\t' '{print $1}' "$Y/routes.tsv" | sort -u | wc -l)
NPATH=$(awk -F'\t' '$3!="CONTROL"{print $3}' "$Y/routes.tsv" | sort -u | wc -l)
echo "  KO families: $NKO   pathways: $NPATH   (+3 positive controls)"

# ------------------------------------------------------------- ko_list
for src in "$K/ko_list" "$X/ko_list" "$WS2/20_substrates/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$Y/ko_list" 2>/dev/null || true; break; }
done
[ -s "$Y/ko_list" ] || { wget -q -O "$Y/ko_list.gz" "$BASE/ko_list.gz" && gunzip -f "$Y/ko_list.gz"; }
[ -s "$Y/ko_list" ] || { echo "  !! ko_list unavailable — the label guard cannot run."
                         echo "     Refusing to repeat step 19 on unverified markers."; exit 1; }

# ------------------------------------------------------------- profiles
got=0
for ko in $(cut -f1 "$Y/routes.tsv" | sort -u); do
  [ -s "$Y/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$X/profiles/${ko}.hmm" "$WS2/20_substrates/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$Y/profiles/"; break; }
  done
  [ -s "$Y/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  if wget -q -O "$Y/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" \
     && head -1 "$Y/profiles/${ko}.hmm" | grep -q HMMER; then got=$((got+1))
  else rm -f "$Y/profiles/${ko}.hmm"; fi
done
echo "  profiles from cache/endpoint: $got/$NKO"

if [ "$got" -lt "$NKO" ]; then
  TARB="$K/profiles.tar.gz"
  [ -s "$TARB" ] || { echo "  fetching the KOfam archive once (~400 MB, resumable)"
                      mkdir -p "$K"; wget -c -O "$TARB" "$BASE/profiles.tar.gz" || true; }
  if [ -s "$TARB" ]; then
    cut -f1 "$Y/routes.tsv" | sort -u | sed 's|^|profiles/|; s|$|.hmm|' > "$Y/wanted.txt"
    tar -xzf "$TARB" -C "$Y" -T "$Y/wanted.txt" 2>/dev/null \
      || tar -xzf "$TARB" -C "$Y" $(tr '\n' ' ' < "$Y/wanted.txt") 2>/dev/null || true
    got=$(ls -1 "$Y/profiles"/*.hmm 2>/dev/null | wc -l)
    echo "  after archive extraction: $got/$NKO"
    cp -n "$Y/profiles"/*.hmm "$K/profiles/" 2>/dev/null || true
  fi
fi
# a partial profile set is exactly what produced 29/29 B6
if [ "$got" -lt $(( NKO * 9 / 10 )) ]; then
  echo "  !! only $got/$NKO profiles available (<90%). This is the condition that"
  echo "     produced the 29/29 B6 result. Refusing to continue."
  exit 1
fi
cat "$Y/profiles"/*.hmm > "$Y/routes.hmm"

# ------------------------------------------------------------- search
echo
echo "######## hmmsearch across all MAGs ($(date +%H:%M)) ########"
# CACHE KEY = HASH OF THE MARKER SET.
# The previous run cached hits under the MAG name alone. When four control
# families were added, every .tbl was reused unchanged, so rplB/rpsC/rplE/ef1a
# were never actually searched and reported 0/7 and 0/5 — which then failed the
# control gate and looked like broken proteomes. A cache whose key omits the
# query is not a cache, it is a silent wrong answer. Hash the HMM set.
HMMHASH=$(md5sum "$Y/routes.hmm" | cut -c1-12)
echo "  marker-set hash: $HMMHASH"
# A MISSING hash is not "no change", it is "provenance unknown". The first
# version only purged on a hash MISMATCH, so on the very first run after this
# code was added there was no .hmmhash to compare against, the stale .tbl files
# from the previous marker set were reused, and rplB/rpsC/rplE/ef1a came back
# 0/29 again. Unverifiable cache = no cache.
if [ ! -f "$Y/hits/.hmmhash" ]; then
  if ls "$Y/hits"/*.tbl >/dev/null 2>&1; then
    echo "  cached hits exist but carry NO marker-set stamp — provenance unknown,"
    echo "  discarding them rather than trusting them"
    rm -f "$Y/hits"/*.tbl
  fi
elif [ "$(cat "$Y/hits/.hmmhash")" != "$HMMHASH" ]; then
  echo "  marker set CHANGED since the cached run — discarding stale hits"
  rm -f "$Y/hits"/*.tbl
fi
n=0; tot=$(ls -1d "$A"/*/ 2>/dev/null | wc -l)
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1)); out="$Y/hits/${name}.tbl"
  if [ -s "$out" ]; then printf '\r  [%2d/%2d] %-28s cached  ' "$n" "$tot" "$name"; continue; fi
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 \
      "$Y/routes.hmm" "$faa" >/dev/null 2>&1 || echo "  !! failed on $name"
done
echo
# stamp AFTER the loop: certify work that happened, not work intended
echo "$HMMHASH" > "$Y/hits/.hmmhash"


# ===========================================================================
python3 - "$Y" "$F" <<'PY'
import sys, os, csv, glob
from collections import defaultdict, Counter

Y, F = sys.argv[1], sys.argv[2]
# ROUTE_T was 0.60. On a TWO-marker route that silently demands 2/2, because
# ceil(0.6 x 2) = 2 — so a route was only "present" if BOTH markers were found,
# while a four-marker route needed 3/4. The three pathways that came back
# >=93% absent (Ala 100%, Tyr 100%, Met 93%) are EXACTLY the three whose every
# scorable route has two markers. That is not biology, it is the threshold
# interacting with route size. At 0.50 a two-marker route needs 1 and a
# four-marker route needs 2, which is consistent across sizes.
ROUTE_T   = 0.50
MIN_ROUTE = 2        # a single-marker route may not decide a pathway on its own
SENS_DELTA = 0.10    # sensitivity sweep around ROUTE_T
MIN_MARK  = 3        # a pathway needs >=3 surviving markers to be scorable
STRICT    = 90.0     # only genomes this complete get a firm auxotrophy call
PREVALENCE_FAIL = 0.80

kegg_def, thr = {}, {}
for r in csv.DictReader(open(os.path.join(Y, 'ko_list')), delimiter='\t'):
    kegg_def[r['knum']] = (r.get('definition') or '').lower()
    try:
        if r.get('score_type','') in ('full','domain') and r.get('threshold','-') != '-':
            thr[r['knum']] = float(r['threshold'])
    except (ValueError, KeyError):
        pass

rows_in = []
for line in open(os.path.join(Y, 'routes.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 5:
        rows_in.append(dict(ko=p[0], gene=p[1], path=p[2], route=p[3], kw=p[4].lower()))

# =====================================================================
# GUARD A — label check. Drop, do not silently mis-score.
# =====================================================================
print("=" * 104)
print("  GUARD A — KO labels vs KEGG definitions")
print("=" * 104)
kept, dropped = [], []
audit = []
for e in rows_in:
    d = kegg_def.get(e['ko'])
    if d is None:
        dropped.append(e); audit.append({**e, 'status':'DROPPED', 'reason':'not in ko_list'})
        print(f"  DROP  {e['ko']} {e['gene']:<8} {e['path']:<16} not in ko_list")
    elif e['kw'] not in d:
        dropped.append(e); audit.append({**e, 'status':'DROPPED', 'reason':f'KEGG says: {d[:60]}'})
        print(f"  DROP  {e['ko']} {e['gene']:<8} {e['path']:<16} I said '{e['kw'][:34]}'")
        print(f"        KEGG says: {d[:82]}")
    else:
        kept.append(e); audit.append({**e, 'status':'kept', 'reason':''})
print(f"  kept {len(kept)} / dropped {len(dropped)} of {len(rows_in)} markers")
with open(os.path.join(Y, 'marker_audit.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t',
                       fieldnames=['ko','gene','path','route','kw','status','reason'])
    w.writeheader(); w.writerows(audit)

# =====================================================================
# LOAD HITS
# =====================================================================
tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
for r in csv.DictReader(open(p), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
assert tax, "MAG_TAXONOMY.tsv unreadable"

hits = defaultdict(set)
for tbl in glob.glob(os.path.join(Y, 'hits', '*.tbl')):
    mag = os.path.basename(tbl)[:-4]
    hits.setdefault(mag, set())
    for line in open(tbl):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 6: continue
        ko, score = f[2], float(f[5])
        if ko not in kegg_def: ko = f[3]
        if ko not in kegg_def: continue
        if ko in thr and score < thr[ko]: continue
        hits[mag].add(ko)
mags = sorted(set(tax) | set(hits))

def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def comp(m):
    try: return float((tax.get(m, {}) or {}).get('completeness', 0) or 0)
    except ValueError: return 0.0

# =====================================================================
# GUARD B — positive control, DOMAIN-AWARE
#
# The first version of this guard FAILED, and it was right to: rpoB found in
# 14/29, tuf in 17/29. But the failure was mine, not the proteomes'. K03043
# (rpoB) and K02358 (tuf) are BACTERIAL. Archaea encode a split RNA polymerase
# B (rpoB'/rpoB'') and EF-1-alpha, not EF-Tu, so all 8 archaeal MAGs were
# guaranteed to miss both. Calling a bacterial gene "universal" and then
# reporting its absence is EXACTLY the error that produced 29/29 B6 at step 19,
# reproduced inside the guard written to prevent it.
#
# Fixed three ways: controls declare which domain they apply to; controls are a
# POOL and at least two must pass per domain, so one wrong KO cannot veto the
# run; and the hit rate is measured only over genomes >=90% complete, because a
# 55%-complete MAG can legitimately be missing anything.
# =====================================================================
print()
print("=" * 104)
print("  GUARD B — positive controls, per domain, over >=90%-complete genomes")
print("=" * 104)

def dom_of(m):
    d = ((tax.get(m, {}) or {}).get('domain', '') or '').lower()
    return 'archaea' if 'arch' in d else ('bacteria' if 'bact' in d else '?')

ctrl = [e for e in kept if e['path'] == 'CONTROL']
ctrl_ok = True
for dm in ('bacteria', 'archaea'):
    pool = [m for m in mags if dom_of(m) == dm and comp(m) >= STRICT]
    if not pool:
        print(f"\n  {dm}: no genome >=90% complete — control not evaluable, skipped")
        continue
    print(f"\n  {dm}: {len(pool)} genome(s) >=90% complete")
    npass = 0
    for e in ctrl:
        if e['route'] not in (dm, 'both'): continue
        n = sum(1 for m in pool if e['ko'] in hits.get(m, set()))
        # allow one miss, or 10%, whichever is larger. On a 7-genome pool a flat
        # 90% rule demands 7/7, so a single truncated gene vetoes the domain.
        allowed = max(1, int(0.1 * len(pool)))
        ok = n >= len(pool) - allowed
        npass += ok
        print(f"    {e['gene']:<8}{e['ko']:<9}{n:>3}/{len(pool)}   {'PASS' if ok else 'miss'}")
    if npass < 2:
        ctrl_ok = False
        print(f"    !! only {npass} control(s) passed for {dm} (need 2)")
    else:
        print(f"    {npass} controls pass — {dm} proteomes are readable")

if not ctrl_ok:
    print()
    print("  !! Fewer than two controls passed in a domain. Either the proteomes or")
    print("     the control KOs are wrong, and auxotrophies computed now would be")
    print("     uninterpretable. Check marker_audit.tsv, then the .faa files.")
    print("     REFUSING TO CONTINUE.")
    raise SystemExit(1)
print("\n  controls pass in every evaluable domain — absences below belong to the")
print("  genomes, not to the search.")

# =====================================================================
# WHY >=3 MARKERS — the arithmetic, printed
# =====================================================================
print()
print("=" * 104)
print("  WHY A PATHWAY NEEDS >=3 MARKERS")
print("=" * 104)
print("  For a genome C% complete, the chance of missing ALL n markers of a pathway")
print("  the organism actually HAS is (1-C)^n. That is the false-auxotrophy rate:")
print(f"  {'completeness':<14}" + ''.join(f"{f'n={n}':>10}" for n in (1,2,3,4)))
for c in (0.70, 0.80, 0.90, 0.95):
    print(f"  {f'{c:.0%}':<14}" + ''.join(f"{(1-c)**n:>9.2%} " for n in (1,2,3,4)))
print("  At 90% complete, one marker gives a 10% false-auxotrophy rate — across 29")
print(f"  genomes and {len(set(e['path'] for e in kept if e['path']!='CONTROL'))} pathways that is dozens of invented auxotrophies.")

# =====================================================================
# SCORING — per route, best route wins
# =====================================================================
by_path = defaultdict(lambda: defaultdict(list))
for e in kept:
    if e['path'] == 'CONTROL': continue
    by_path[e['path']][e['route']].append(e['ko'])

paths = sorted(by_path)
scorable, unscorable = [], []
for p_ in paths:
    nmark = sum(len(v) for v in by_path[p_].values())
    if nmark < MIN_MARK:
        unscorable.append((p_, f"only {nmark} marker(s) survived the label guard "
                               f"(minimum {MIN_MARK})"))
    else:
        scorable.append(p_)

def route_present(m, p_, thr):
    for rt, kos in by_path[p_].items():
        if rt in ('context', 'racemase_context', 'salvage'): continue
        if len(kos) < MIN_ROUTE: continue
        if sum(1 for k in kos if k in hits.get(m, set())) / len(kos) >= thr:
            return True
    return False

def best_route(m, p_):
    """Completeness of the BEST route, not the fraction of the pooled union."""
    best, name = 0.0, ''
    for rt, kos in by_path[p_].items():
        if rt in ('context', 'racemase_context', 'salvage'): continue
        # a one-marker route would call the pathway present on a single hit,
        # which is the mirror image of the step-19 error: a false PROTOtrophy.
        if len(kos) < MIN_ROUTE: continue
        f = sum(1 for k in kos if k in hits.get(m, set())) / len(kos)
        if f > best: best, name = f, rt
    return best, name

# prevalence flag (GUARD C)
print()
print("=" * 104)
print("  GUARD C — prevalence flag: a pathway absent in >80% of genomes is a")
print("  marker-set failure until proven otherwise")
print("=" * 104)
print(f"  {'pathway':<18}{'routes':<8}{'mk':>4}{'minRt':>6}"
      f"{'%abs@.4':>9}{'%abs@.5':>9}{'%abs@.6':>9}   verdict")
print("  minRt = markers in the SMALLEST scorable route. A pathway whose routes")
print("  are all 2 markers is the one a high threshold punishes; the sweep makes")
print("  that visible instead of letting the threshold decide the biology.")
final = []
for p_ in scorable:
    nabs = sum(1 for m in mags if not route_present(m, p_, ROUTE_T))
    frac = nabs / max(1, len(mags))
    nmark = sum(len(v) for v in by_path[p_].values())
    sweep = {}
    for t in (ROUTE_T - SENS_DELTA, ROUTE_T, ROUTE_T + SENS_DELTA):
        sweep[t] = sum(1 for m in mags if not route_present(m, p_, t)) / max(1, len(mags))
    minrt = min((len(v) for rt, v in by_path[p_].items()
                 if rt not in ('context','racemase_context','salvage') and len(v) >= MIN_ROUTE),
                default=0)
    spread = max(sweep.values()) - min(sweep.values())
    if frac > PREVALENCE_FAIL:
        unscorable.append((p_, f"absent in {100*frac:.0f}% of genomes — "
                               "suspected marker-set failure, not auxotrophy"))
        v = 'SUSPECT -> unscorable'
    else:
        final.append(p_); v = 'ok'
    if spread > 0.20 and v == 'ok':
        v = f'ok (THRESHOLD-SENSITIVE, {100*spread:.0f}pt spread)'
    print(f"  {p_:<18}{len(by_path[p_]):<8}{nmark:>4}{minrt:>6}"
          + ''.join(f"{100*sweep[t]:>8.0f}%" for t in sorted(sweep))
          + f"   {v}")

print()
print(f"  scorable pathways: {len(final)}   unscorable: {len(unscorable)}")
with open(os.path.join(Y, 'UNSCORABLE.txt'), 'w') as fh:
    fh.write("Pathways this method cannot decide, and why.\n")
    fh.write("Reporting 'not scorable' is a result. Reporting a false auxotrophy is not.\n\n")
    for p_, why in unscorable:
        fh.write(f"{p_}\t{why}\n")
        print(f"    UNSCORABLE  {p_}: {why}")

# =====================================================================
# OUTPUT
# =====================================================================
with open(os.path.join(Y, 'pathway_completeness.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','sample','domain','genus','completeness'] + final)
    for m in mags:
        t = tax.get(m, {}) or {}
        w.writerow([t.get('mag_id',''), m, t.get('sample',''), t.get('domain',''),
                    t.get('genus',''), t.get('completeness','')]
                   + [f"{best_route(m,p_)[0]:.2f}" for p_ in final])

print()
print("=" * 104)
print("  FIRM AUXOTROPHIES  (genomes >=90% complete only)")
print("=" * 104)
aux_rows = []
firm_mags = [m for m in mags if comp(m) >= STRICT]
print(f"  genomes complete enough to trust an absence: {len(firm_mags)}/{len(mags)}")
print(f"  {'MAG ID':<17}{'genus':<20}{'compl':>7}   missing")
for m in sorted(firm_mags, key=mid):
    miss = [p_ for p_ in final if not route_present(m, p_, ROUTE_T)]
    t = tax.get(m, {}) or {}
    print(f"  {mid(m):<17}{(t.get('genus','-') or '-')[:19]:<20}{comp(m):>6.1f}%   "
          + (', '.join(miss) if miss else 'none — prototrophic for every scorable pathway'))
    for p_ in miss:
        f_, rt = best_route(m, p_)
        aux_rows.append(dict(mag_id=mid(m), working_name=m, completeness=f"{comp(m):.1f}",
                             pathway=p_, best_route=rt, route_completeness=f"{f_:.2f}"))
with open(os.path.join(Y, 'auxotrophies.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=[
        'mag_id','working_name','completeness','pathway','best_route','route_completeness'])
    w.writeheader(); w.writerows(aux_rows)

# =====================================================================
# BEFORE / AFTER
# =====================================================================
print()
print("=" * 104)
print("  WHAT CHANGED vs STEP 19")
print("=" * 104)
for p_ in ('B6_pyridoxal', 'Ala', 'Tyr'):
    if p_ in final:
        nabs = sum(1 for m in mags if best_route(m, p_)[0] < ROUTE_T)
        routes = ', '.join(sorted(by_path[p_]))
        print(f"  {p_:<16} step 19: near-universal 'auxotrophy'")
        print(f"  {'':<16} step 19b: {nabs}/{len(mags)} absent, routes scored = {routes}")
    else:
        why = next((w for q, w in unscorable if q == p_), 'not in the marker set')
        print(f"  {p_:<16} step 19b: UNSCORABLE — {why}")
print()
print("  The B6 case is the diagnostic one: step 19 scored only the DXP-DEPENDENT")
print("  route, which most bacteria and ALL archaea do not use. Adding pdxS/pdxT")
print("  (DXP-independent) and scoring per route is the entire fix.")
print()
print("  REPORT AS: 'Biosynthetic capacity was scored per alternative route rather")
print("  than as a fraction of pooled markers, requiring at least three markers per")
print("  pathway; pathways absent from more than 80% of genomes were treated as")
print("  marker-set failures and excluded. Auxotrophies were called only in genomes")
print("  >=90% complete.'")
PY

echo
echo "######## VALIDATION GATE ########"
check "pathway completeness" $([ -s "$Y/pathway_completeness.tsv" ] && echo 0 || echo 1)
check "marker audit"         $([ -s "$Y/marker_audit.tsv" ] && echo 0 || echo 1)
check "unscorable list"      $([ -s "$Y/UNSCORABLE.txt" ] && echo 0 || echo 1)
NP=$(( $(head -1 "$Y/pathway_completeness.tsv" 2>/dev/null | tr '\t' '\n' | wc -l) - 6 ))
check "at least 8 scorable pathways" $([ "$NP" -ge 8 ] && echo 0 || echo 1)

runlog "19b_auxotrophy complete"
echo
echo "######## 19b DONE ########"
echo "Step 21 picks up 19b automatically if this directory exists."
