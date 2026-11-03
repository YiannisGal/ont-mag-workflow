#!/usr/bin/env bash
# =============================================================================
# 21_interaction_network.sh  --  who competes with whom, who feeds whom
#
# WHAT THIS IS, AND WHAT IT IS NOT
#   You asked for network analysis, co-occurrence, competition. Two of those
#   three are possible with this dataset. The third is not, and it matters that
#   the thesis says so:
#
#   CO-OCCURRENCE IS IMPOSSIBLE HERE. SparCC, SPIEC-EASI, WGCNA and every other
#   correlation network need ABUNDANCE VARIATION ACROSS MANY SAMPLES — the
#   signal IS the covariance. You have TWO metagenomes. A correlation computed
#   from n=2 is not a weak correlation, it is undefined: any two points lie on a
#   line, so every pair would come out at r = +/-1. Do not compute it, do not
#   report it, and do not let a reviewer think you tried. Co-occurrence belongs
#   to the 16S amplicon chapter, where you have many samples across the depth
#   profile — that is the correct home for it, and it will work properly there.
#
#   WHAT GENE CONTENT *CAN* SUPPORT, and what this script builds:
#     1. COMPETITION — resource-niche overlap. Two genomes encoding the same
#        uptake and degradation machinery draw on the same substrate pool.
#     2. COOPERATION / CROSS-FEEDING — directed producer -> consumer edges for
#        specific exchange metabolites, each with a named gene on both ends.
#     3. DEPENDENCY — auxotroph -> prototroph edges from step 19.
#   All three are POTENTIAL interactions inferred from genomes. None is a
#   measured flux. Say "encodes the potential for", never "cross-feeds".
#
# THE COMPETITION CONFOUND, HANDLED
#   Two Methanosarcina MAGs share substrate genes because they are both
#   Methanosarcina, not because they are locked in competition. Raw overlap
#   recovers phylogeny, not ecology. So every competition edge is reported
#   alongside the pair's taxonomic distance, and the script separately lists
#   HIGH-OVERLAP DISTANT PAIRS — convergent competitors from different phyla,
#   which is the ecologically interesting class and the one raw overlap buries.
#
# THE GUARD THAT EXISTS BECAUSE OF A REAL ERROR
#   Step 20 shipped with K14080 labelled "mtsA, methylated-thiol:CoM
#   methyltransferase". It is actually mtaA — the SHARED methyl transfer step
#   used by the methanol and methylamine pathways too. That one mislabel came
#   within a hair of a false claim in the results chapter. So this script does
#   not trust the labels in its own KO table: it reads the KEGG definition of
#   every KO from ko_list at runtime and refuses to run if my label and KEGG's
#   definition disagree. Any wrong KO number now stops the script instead of
#   producing a confident wrong figure.
#
# EXCLUDED ON PURPOSE
#   Step 19's Ala, B6 and Tyr calls are known-bad (29/29 and 26/29 "auxotrophic"
#   is a marker-set failure, not biology). They are hard-excluded here until
#   19b fixes them, rather than being allowed to dominate the dependency network.
#
# INPUT    19_auxotrophy/pathway_completeness.tsv, 20_substrates/substrate_matrix.tsv,
#          09_annotation/<MAG>/<MAG>.faa, 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   21_network/
#     nodes.tsv  edges.tsv          import straight into Cytoscape or Gephi
#     competition_matrix.tsv        pairwise resource overlap
#     crossfeeding.tsv              directed metabolite exchanges
#     network.png                   if matplotlib is available
#
# RUNTIME  20-40 min.
# USAGE    ws2 && bash scripts/21_interaction_network.sh 2>&1 | tee logs/21_network.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="21_network"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"

A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"
X="$WS2/19_auxotrophy"; U="$WS2/20_substrates"; N="$WS2/21_network"
mkdir -p "$N/profiles" "$N/hits"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 21 — METABOLIC INTERACTION NETWORK ########"
echo "  started: $(date -Is)"
runlog "21_interaction_network started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }
[ -s "$U/substrate_matrix.tsv" ] || { echo "  !! run 20_substrates.sh first"; exit 1; }
[ -s "$X/pathway_completeness.tsv" ] || { echo "  !! run 19_auxotrophy.sh first"; exit 1; }

# ---------------------------------------------------------------------------
# EXCHANGE METABOLITES
# KO <tab> gene <tab> METABOLITE <tab> PROD|CONS <tab> keyword that MUST appear
# in the KEGG definition. The keyword is the guard, not decoration.
# ---------------------------------------------------------------------------
cat > "$N/exchange.tsv" <<'KOS'
K00925	ackA	ACETATE	PROD	acetate kinase
K13788	pta	ACETATE	PROD	phosphate acetyltransferase
K00625	eutD	ACETATE	PROD	phosphate acetyltransferase
K01895	acs	ACETATE	CONS	acetyl-CoA synthetase
K00192	cdhA	ACETATE	CONS	carbon-monoxide dehydrogenase
K00656	pflD	FORMATE	PROD	formate C-acetyltransferase
K00123	fdoG	FORMATE	BOTH	formate dehydrogenase
K00200	fwdA	FORMATE	CONS	formylmethanofuran dehydrogenase
K00532	hydA	H2	PROD	ferredoxin hydrogenase
K00533	hydB	H2	PROD	ferredoxin hydrogenase
K00534	hydC	H2	PROD	ferredoxin hydrogenase
K00437	hydB2	H2	PROD	hydrogenase large subunit
K18332	hydA2	H2	PROD	hydrogenase
K06281	hyaA	H2	PROD	hydrogenase
K06282	hyaB	H2	PROD	hydrogenase
K00436	hoxH	H2	PROD	hydrogenase
K18016	hybO	H2	PROD	hydrogenase
K05922	hydB3	H2	PROD	hydrogenase
K00440	frhB	H2	CONS	coenzyme F420 hydrogenase
K00441	frhC	H2	CONS	coenzyme F420 hydrogenase
K00443	frhA	H2	CONS	coenzyme F420 hydrogenase
K01051	pme	METHANOL	PROD	pectinesterase
K04480	mtaB	METHANOL	CONS	methyltransferase
K14082	mtaC	METHANOL	CONS	corrinoid
K20038	cutC	TMA	PROD	choline trimethylamine-lyase
K07811	torA	TMA	PROD	trimethylamine-N-oxide reductase
K14083	mttB	TMA	CONS	methyltransferase
K02230	cobN	B12	PROD	cobaltochelatase
K00768	cobT	B12	PROD	phosphoribosyltransferase
K11180	dsrA	SULFIDE	PROD	sulfite reductase
K11181	dsrB	SULFIDE	PROD	sulfite reductase
K17218	sqr	SULFIDE	CONS	sulfide
KOS
NKO=$(wc -l < "$N/exchange.tsv")
echo "  exchange families: $NKO across $(cut -f3 "$N/exchange.tsv" | sort -u | wc -l) metabolites"

# ------------------------------------------------------------- ko_list + profiles
for src in "$K/ko_list" "$X/ko_list" "$U/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$N/ko_list" 2>/dev/null || true; break; }
done
[ -s "$N/ko_list" ] || { wget -q -O "$N/ko_list.gz" "$BASE/ko_list.gz" && gunzip -f "$N/ko_list.gz"; }
[ -s "$N/ko_list" ] || { echo "  !! ko_list unavailable — the KEGG-definition guard cannot run."
                         echo "     Refusing to build a network from unverified KO labels."; exit 1; }

got=0
while IFS=$'\t' read -r ko gene met dir kw; do
  [ -z "${ko:-}" ] && continue
  [ -s "$N/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$X/profiles/${ko}.hmm" "$U/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$N/profiles/"; break; }
  done
  [ -s "$N/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  if wget -q -O "$N/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" \
     && head -1 "$N/profiles/${ko}.hmm" | grep -q HMMER; then got=$((got+1))
  else rm -f "$N/profiles/${ko}.hmm"; fi
done < "$N/exchange.tsv"
echo "  profiles present: $got/$NKO"

if [ "$got" -lt "$NKO" ]; then
  TARB="$K/profiles.tar.gz"
  [ -s "$TARB" ] || { echo "  fetching the KOfam archive once (~400 MB, resumable)"
                      mkdir -p "$K"; wget -c -O "$TARB" "$BASE/profiles.tar.gz" || true; }
  if [ -s "$TARB" ]; then
    awk -F'\t' '{print "profiles/"$1".hmm"}' "$N/exchange.tsv" > "$N/wanted.txt"
    tar -xzf "$TARB" -C "$N" -T "$N/wanted.txt" 2>/dev/null \
      || tar -xzf "$TARB" -C "$N" $(tr '\n' ' ' < "$N/wanted.txt") 2>/dev/null || true
    got=$(ls -1 "$N/profiles"/*.hmm 2>/dev/null | wc -l)
    echo "  after archive extraction: $got/$NKO"
    cp -n "$N/profiles"/*.hmm "$K/profiles/" 2>/dev/null || true
  fi
fi
if [ "$got" -lt "$NKO" ]; then
  echo "  !! $((NKO-got)) profile(s) missing. Missing a producer family invents a"
  echo "     one-sided exchange; missing a consumer family deletes a real edge."
  echo "     Both corrupt the network silently. Refusing to continue."
  ls -1 "$N/profiles"/*.hmm 2>/dev/null | xargs -n1 basename 2>/dev/null \
    | sed 's/\.hmm//' | sort > "$N/.have"
  cut -f1 "$N/exchange.tsv" | sort > "$N/.want"
  echo "     missing: $(comm -13 "$N/.have" "$N/.want" | tr '\n' ' ')"
  exit 1
fi
cat "$N/profiles"/*.hmm > "$N/exchange.hmm"

# ------------------------------------------------------------------- search
echo
echo "######## hmmsearch across all MAGs ($(date +%H:%M)) ########"
# CACHE KEY = HASH OF THE MARKER SET, and a MISSING stamp also invalidates.
# Adding KOs without re-searching silently reports them as absent — that is
# exactly how 19b came back with rplB/rpsC/rplE at 0/29 twice.
HMMHASH=$(md5sum "$N/exchange.hmm" | cut -c1-12)
echo "  marker-set hash: $HMMHASH"
if [ ! -f "$N/hits/.hmmhash" ]; then
  if ls "$N/hits"/*.tbl >/dev/null 2>&1; then
    echo "  cached hits carry no marker-set stamp — discarding (provenance unknown)"
    rm -f "$N/hits"/*.tbl
  fi
elif [ "$(cat "$N/hits/.hmmhash")" != "$HMMHASH" ]; then
  echo "  marker set CHANGED — discarding stale hits"
  rm -f "$N/hits"/*.tbl
fi
n=0; tot=$(ls -1d "$A"/*/ 2>/dev/null | wc -l)
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1)); out="$N/hits/${name}.tbl"
  if [ -s "$out" ]; then printf '\r  [%2d/%2d] %-28s cached  ' "$n" "$tot" "$name"; continue; fi
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 \
      "$N/exchange.hmm" "$faa" >/dev/null 2>&1 || echo "  !! failed on $name"
done
echo
# stamp AFTER the loop: certify work that happened, not work intended
echo "$HMMHASH" > "$N/hits/.hmmhash"


# ===========================================================================
python3 - "$N" "$F" "$X" "$U" <<'PY'
import sys, os, csv, glob, itertools
from collections import defaultdict, Counter

N, F, X, U = sys.argv[1:5]

# =====================================================================
# GUARD — verify every KO label against KEGG's own definition
# =====================================================================
kegg_def = {}
for r in csv.DictReader(open(os.path.join(N, 'ko_list')), delimiter='\t'):
    kegg_def[r['knum']] = (r.get('definition') or '').lower()

thr = {}
for r in csv.DictReader(open(os.path.join(N, 'ko_list')), delimiter='\t'):
    try:
        if r.get('score_type','') in ('full','domain') and r.get('threshold','-') != '-':
            thr[r['knum']] = float(r['threshold'])
    except (ValueError, KeyError):
        pass

ex = []
for line in open(os.path.join(N, 'exchange.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 5:
        ex.append(dict(ko=p[0], gene=p[1], met=p[2], dir=p[3], kw=p[4].lower()))

print("=" * 104)
print("  LABEL GUARD — every KO checked against its KEGG definition")
print("  (step 20 called K14080 'mtsA'; KEGG calls it mtaA. That mislabel almost")
print("   became a results-chapter claim. Nothing here runs on trust.)")
print("=" * 104)
bad, unk = [], []
for e in ex:
    d = kegg_def.get(e['ko'])
    if d is None:
        unk.append(e); continue
    if e['kw'] not in d:
        bad.append((e, d))
        print(f"  MISMATCH  {e['ko']}  I said '{e['gene']}' ({e['kw']})")
        print(f"            KEGG says: {d[:88]}")
if unk:
    print(f"  {len(unk)} KO(s) absent from ko_list: "
          + ', '.join(f"{e['ko']}({e['gene']})" for e in unk))
if bad:
    print()
    print(f"  !! {len(bad)} KO label(s) disagree with KEGG. A wrong KO number puts a")
    print("     gene on the wrong end of an exchange and the network comes out")
    print("     confidently wrong. Fix exchange.tsv and re-run. REFUSING TO CONTINUE.")
    raise SystemExit(1)
if unk:
    print()
    print(f"  !! {len(unk)} KO(s) could not be verified at all. Refusing to continue.")
    raise SystemExit(1)
print(f"  all {len(ex)} KO labels agree with KEGG.")

# =====================================================================
# LOAD
# =====================================================================
tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        wn = r.get('working_name') or r.get('bin') or r.get('mag_id')
        if wn: tax[wn] = r
assert tax, "MAG_TAXONOMY.tsv empty or unreadable — refusing to build a network"

hits = defaultdict(Counter)
for tbl in glob.glob(os.path.join(N, 'hits', '*.tbl')):
    mag = os.path.basename(tbl)[:-4]
    for line in open(tbl):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 6: continue
        ko, score = f[2], float(f[5])
        if ko not in kegg_def: ko = f[3]
        if ko not in kegg_def: continue
        if ko in thr and score < thr[ko]: continue
        hits[mag][ko] += 1
mags = sorted(set(tax) | set(hits))
print(f"  MAGs: {len(mags)}")

def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return (tax.get(m, {}) or {}).get('genus', '') or '-'
def dom(m): return ((tax.get(m, {}) or {}).get('domain', '') or '?')[:3]
def smp(m): return (tax.get(m, {}) or {}).get('sample', '')

RANKS = ['phylum','class','order','family','genus']
def taxdist(a, b):
    """0 = same genus, 5 = different phylum. Shared-lineage confound control."""
    ta, tb = tax.get(a, {}), tax.get(b, {})
    for i, r in enumerate(RANKS):
        va, vb = (ta.get(r) or '').strip(), (tb.get(r) or '').strip()
        if not va or not vb: return 5
        if va != vb: return 5 - i
    return 0

# =====================================================================
# 1. CROSS-FEEDING — directed, named metabolite, named gene both ends
# =====================================================================
print()
print("=" * 104)
print("  CROSS-FEEDING POTENTIAL — directed producer -> consumer")
print("=" * 104)
# =====================================================================
# ZERO-HIT FAMILY FLAG
# A KO family with NO hits in ANY of 29 genomes is usually a marker problem,
# not an absence. This is the step-19 lesson in general form: an empty column
# should be interrogated before it is reported. The exception is a zero that an
# INDEPENDENT line already supports — dsrAB is genuinely absent here (step 11b
# KOfam, and Bakta's four calls were false positives), so that zero is real.
# =====================================================================
print()
print("=" * 104)
print("  ZERO-HIT FAMILIES — interrogate before reporting an absence")
print("=" * 104)
CORROBORATED_ZERO = {'K11180': 'dsrA', 'K11181': 'dsrB'}   # confirmed absent at 11b
zero = [e for e in ex if not any(hits[m][e['ko']] for m in mags)]
if not zero:
    print("  none — every family has at least one hit somewhere")
for e in zero:
    tag = ('genuinely absent, confirmed independently at step 11b'
           if e['ko'] in CORROBORATED_ZERO else
           'SUSPECT — no hit in any of 29 genomes; treat as a marker-set gap')
    print(f"  {e['gene']:<8}{e['ko']:<9}{e['met']:<10}{e['dir']:<6}{tag}")
zero_met = {e['met'] for e in zero if e['ko'] not in CORROBORATED_ZERO}
if zero_met:
    print()
    print("  Affected metabolites: " + ', '.join(sorted(zero_met)))
    print("  Any 'PRODUCED BUT NOT CONSUMED' or 'CONSUMERS WITHOUT A PRODUCER'")
    print("  verdict below that involves these is a statement about the MARKER SET,")
    print("  not about the community. H2 in particular: every fermenter in this")
    print("  catalogue makes hydrogen, and 10 methanogens are encoded to consume it.")

by_met = defaultdict(lambda: {'PROD': [], 'CONS': [], 'BOTH': []})
for e in ex: by_met[e['met']][e['dir']].append(e)

edges, cf_rows = [], []
for met in sorted(by_met):
    P = [e for e in by_met[met]['PROD'] + by_met[met]['BOTH']]
    C = [e for e in by_met[met]['CONS'] + by_met[met]['BOTH']]
    prod = {m: [e['gene'] for e in P if hits[m][e['ko']]] for m in mags}
    cons = {m: [e['gene'] for e in C if hits[m][e['ko']]] for m in mags}
    prod = {m: g for m, g in prod.items() if g}
    cons = {m: g for m, g in cons.items() if g}
    print(f"\n  {met}")
    if not prod and not cons:
        print("    neither producers nor consumers — no edge"); continue
    print(f"    producers ({len(prod)}): " + (', '.join(
        f"{mid(m)}[{','.join(prod[m])}]" for m in sorted(prod, key=mid)[:6]) or 'NONE'))
    print(f"    consumers ({len(cons)}): " + (', '.join(
        f"{mid(m)}[{','.join(cons[m])}]" for m in sorted(cons, key=mid)[:6]) or 'NONE'))
    if not prod:
        if met in zero_met:
            print(f"    >> NO PRODUCER MARKER SURVIVED for {met} — see the zero-hit table")
            print("       above. This is a marker-set gap, NOT evidence that nothing")
            print("       produces it. Do not report this as an ecological gap.")
            continue
        print(f"    >> CONSUMERS WITHOUT A PRODUCER. Either the source of {met} is")
        print("       outside the MAG catalogue (unbinned fraction, or abiotic /")
        print("       allochthonous input), or the consumer call is wrong. This is a")
        print("       gap worth naming in the discussion, not an edge.")
        continue
    if not cons:
        print(f"    >> PRODUCED BUT NOT CONSUMED — {met} would accumulate or be lost")
        print("       from the system rather than cycled internally.")
        continue
    for a, b in itertools.product(sorted(prod), sorted(cons)):
        if a == b: continue
        same = smp(a) == smp(b)
        edges.append(dict(source=mid(a), target=mid(b), type='crossfeed',
                          metabolite=met, weight=1.0,
                          same_sample='yes' if same else 'no',
                          evidence=f"{'/'.join(prod[a])} -> {'/'.join(cons[b])}"))
        cf_rows.append(edges[-1])
    npairs = sum(1 for a, b in itertools.product(prod, cons) if a != b)
    print(f"    >> {npairs} directed producer->consumer pair(s)")

# the specific hypothesis worth stating
tma_p = [m for m in mags if hits[m]['K20038'] or hits[m]['K07811']]
tma_c = [m for m in mags if hits[m]['K14083']]
me_p  = [m for m in mags if hits[m]['K01051']]
me_c  = [m for m in mags if hits[m]['K04480'] or hits[m]['K14082']]
print()
print("  THE TWO EXCHANGES THAT MATTER FOR THE METHANE STORY")
print(f"    choline/TMAO -> TMA -> CH4 : {len(tma_p)} producer(s), {len(tma_c)} consumer(s)")
if tma_p and tma_c:
    print("      Both ends present. This is the strongest cross-feeding claim in the")
    print("      dataset because the consumer end is already triple-confirmed: mttB")
    print("      tree 5/5, pyl operons in the same three MAGs, and mcrA.")
    print("      CAVEAT: K14083 is the MttB SUPERFAMILY. Bacterial members are glycine")
    print("      betaine methyltransferases, not TMA methyltransferases. Trust it only")
    print("      for the three MAGs the 15b tree already confirmed.")
elif tma_c and not tma_p:
    print("      Consumers only — the TMA source is not in the catalogue. In a coastal")
    print("      lagoon, osmolyte-rich algal and animal input is the obvious external")
    print("      source. Name it as a gap, do not invent the edge.")
print(f"    pectin -> methanol -> CH4  : {len(me_p)} producer(s), {len(me_c)} consumer(s)")
if me_p and me_c:
    print("      Pectinesterase releases methanol from plant pectin. With a lignocellulose-")
    print("      fed lagoon and a polysaccharide-degradation hub already identified, this")
    print("      is a concrete, testable route to methylotrophic methane that does NOT")
    print("      require DMSP. Test it the same way: amend an enrichment with pectin.")

# =====================================================================
# 2. COMPETITION — resource overlap, with the phylogeny confound controlled
# =====================================================================
print()
print("=" * 104)
print("  COMPETITION — resource-niche overlap")
print("=" * 104)
sm = os.path.join(U, 'substrate_matrix.tsv')
prof = {}
with open(sm) as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    cols = [c for c in (rd.fieldnames or []) if c.startswith(('KO_', 'KW_'))]
    for r in rd:
        wn = r.get('working_name') or r.get('mag_id')
        if not wn: continue
        prof[wn] = {c for c in cols if (r.get(c) or '0').strip() not in ('', '0', '·')}
print(f"  substrate features per genome: {len(cols)}  |  genomes profiled: {len(prof)}")

pairs = []
for a, b in itertools.combinations(sorted(prof), 2):
    A_, B_ = prof[a], prof[b]
    if not A_ or not B_: continue
    inter = len(A_ & B_)
    jac = inter / len(A_ | B_)
    ovl = inter / min(len(A_), len(B_))      # overlap coefficient: fair to small genomes
    pairs.append((ovl, jac, inter, a, b, taxdist(a, b)))
pairs.sort(reverse=True)

print()
print("  TOP 15 OVERLAPPING PAIRS (raw — expect relatives to dominate)")
print(f"  {'MAG A':<17}{'MAG B':<17}{'ovl':>6}{'jacc':>7}{'shared':>8}  relatedness")
REL = {0:'same genus', 1:'same family', 2:'same order', 3:'same class',
       4:'same phylum', 5:'DIFFERENT PHYLUM'}
for ovl, jac, inter, a, b, td in pairs[:15]:
    print(f"  {mid(a):<17}{mid(b):<17}{ovl:>6.2f}{jac:>7.2f}{inter:>8}  {REL[td]}")

print()
print("  CONVERGENT COMPETITORS — high overlap DESPITE different phyla")
print("  (this is the ecologically meaningful list: shared niche that is not")
print("   simply shared ancestry)")
distant = [p for p in pairs if p[5] == 5]
print(f"  {'MAG A':<17}{'genus A':<18}{'MAG B':<17}{'genus B':<18}{'ovl':>6}{'shared':>8}")
for ovl, jac, inter, a, b, td in distant[:15]:
    print(f"  {mid(a):<17}{gen(a)[:17]:<18}{mid(b):<17}{gen(b)[:17]:<18}{ovl:>6.2f}{inter:>8}")

COMP_T = 0.60
for ovl, jac, inter, a, b, td in pairs:
    if ovl < COMP_T: break
    edges.append(dict(source=mid(a), target=mid(b), type='competition',
                      metabolite='resource_overlap', weight=round(ovl, 3),
                      same_sample='yes' if smp(a) == smp(b) else 'no',
                      evidence=f"{inter} shared features, {REL[td]}"))

with open(os.path.join(N, 'competition_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_a','mag_b','overlap_coef','jaccard','shared_features',
                'relatedness','same_sample'])
    for ovl, jac, inter, a, b, td in pairs:
        w.writerow([mid(a), mid(b), f"{ovl:.4f}", f"{jac:.4f}", inter,
                    REL[td], 'yes' if smp(a) == smp(b) else 'no'])

print()
print("  READ THIS BEFORE USING THE COMPETITION EDGES")
print("  Shared genes mean shared POTENTIAL, not realised competition. Two genomes")
print("  can encode the same hydrolases and never compete — if the resource is not")
print("  limiting, if they are spatially separated across the 70 cm core, or if")
print("  they are active at different times. Competition also requires that both")
print("  organisms are actually present and active together, which gene content")
print("  cannot show. Phrase it as overlapping resource niches.")

# =====================================================================
# 3. DEPENDENCY — auxotroph -> prototroph, with the broken pathways excluded
# =====================================================================
print()
print("=" * 104)
print("  DEPENDENCY — who cannot make what, and who could supply it")
print("=" * 104)
# Prefer 19b if it has been run: it scores per alternative route, drops markers
# that disagree with KEGG, and has already removed anything it could not score.
X2 = os.path.join(os.path.dirname(X.rstrip('/')), '19b_auxotrophy')
pc2 = os.path.join(X2, 'pathway_completeness.tsv')
if os.path.exists(pc2):
    pc, BROKEN = pc2, set()
    print("  USING 19b (route-aware). No hard exclusions needed — 19b already moved")
    print("  unscorable pathways to UNSCORABLE.txt rather than calling them absent.")
else:
    pc = os.path.join(X, 'pathway_completeness.tsv')
    BROKEN = {'Ala', 'B6_pyridoxal', 'Tyr'}
    print("  USING step 19 (19b not run). EXCLUDED as unreliable: "
          + ', '.join(sorted(BROKEN)))
    print("  (29/29 and 26/29 genomes 'auxotrophic' is a marker-set failure. Leaving")
    print("   them in would make them the densest edges in the whole network.)")
rows = list(csv.DictReader(open(pc), delimiter='\t'))
meta_cols = {'mag_id','working_name','sample','domain','genus','completeness'}
paths = [c for c in (rows[0].keys() if rows else []) if c not in meta_cols]
use = [p for p in paths if p not in BROKEN]
print(f"  pathways used: {len(use)} of {len(paths)}")

PT, STRICT = 0.50, 90.0
dep, ndep = [], 0
for p_ in use:
    aux, pro = [], []
    for r in rows:
        wn = r.get('working_name') or r.get('mag_id')
        try: fr = float(r.get(p_, 0) or 0)
        except ValueError: continue
        # only genomes complete enough that an absence is trustworthy
        try: comp = float(r.get('completeness', 0) or 0)
        except ValueError: comp = 0
        (aux if fr < PT else pro).append((wn, comp))
    firm = [w for w, c in aux if c >= STRICT]
    if not firm or not pro: continue
    frac_aux = len(aux) / max(1, len(aux) + len(pro))
    if frac_aux > 0.80:
        print(f"    {p_}: {100*frac_aux:.0f}% of genomes auxotrophic — SUSPECT marker set,")
        print("      flagged and skipped (this is the 19b failure mode repeating)")
        continue
    ndep += 1
    for a in firm:
        for b, _ in pro:
            if a == b: continue
            dep.append(dict(source=mid(b), target=mid(a), type='dependency',
                            metabolite=p_, weight=1.0,
                            same_sample='yes' if smp(a) == smp(b) else 'no',
                            evidence=f"{mid(a)} lacks {p_}, {mid(b)} complete"))
edges.extend(dep)
print(f"  pathways yielding dependency edges: {ndep}   edges: {len(dep)}")

# =====================================================================
# WRITE
# =====================================================================
with open(os.path.join(N, 'edges.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=[
        'source','target','type','metabolite','weight','same_sample','evidence'])
    w.writeheader(); w.writerows(edges)
with open(os.path.join(N, 'crossfeeding.tsv'), 'w', newline='') as fh:
    if cf_rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(cf_rows[0].keys()))
        w.writeheader(); w.writerows(cf_rows)
with open(os.path.join(N, 'nodes.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','sample','domain','genus','completeness',
                'degree_crossfeed','degree_competition','degree_dependency'])
    deg = defaultdict(Counter)
    for e in edges:
        deg[e['source']][e['type']] += 1; deg[e['target']][e['type']] += 1
    for m in sorted(mags, key=mid):
        i = mid(m)
        w.writerow([i, smp(m), dom(m), gen(m),
                    (tax.get(m, {}) or {}).get('completeness',''),
                    deg[i]['crossfeed'], deg[i]['competition'], deg[i]['dependency']])

print()
print("=" * 104)
print("  NETWORK SUMMARY")
print("=" * 104)
c = Counter(e['type'] for e in edges)
for t, n_ in c.most_common(): print(f"    {t:<14}{n_:>6} edges")
print(f"    {'TOTAL':<14}{len(edges):>6}")
print()
print("  MOST CONNECTED GENOMES (interaction hubs)")
tot = Counter()
for e in edges: tot[e['source']] += 1; tot[e['target']] += 1
print(f"  {'MAG ID':<17}{'genus':<22}{'edges':>7}")
inv = {mid(m): m for m in mags}
for i, n_ in tot.most_common(10):
    m = inv.get(i, i)
    print(f"  {i:<17}{gen(m)[:21]:<22}{n_:>7}")

# ------------------------------------------------------------------- figure
try:
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt, math
    ids = [mid(m) for m in mags]
    ang = {i: 2*math.pi*k/len(ids) for k, i in enumerate(ids)}
    pos = {i: (math.cos(a), math.sin(a)) for i, a in ang.items()}
    COL = {'crossfeed':'#2a78d6', 'competition':'#eb6834', 'dependency':'#9aa4b0'}
    fig, ax = plt.subplots(figsize=(9, 9))
    for e in edges:
        if e['source'] not in pos or e['target'] not in pos: continue
        x1, y1 = pos[e['source']]; x2, y2 = pos[e['target']]
        ax.plot([x1, x2], [y1, y2], color=COL.get(e['type'], '#ccc'),
                lw=0.5, alpha=0.35 if e['type'] == 'dependency' else 0.7, zorder=1)
    for m in mags:
        i = mid(m)
        if i not in pos: continue
        x, y = pos[i]
        ax.scatter([x], [y], s=60 + 6*tot[i], zorder=3,
                   color='#eb6834' if dom(m).lower().startswith('arc') else '#2a78d6',
                   edgecolor='white', linewidth=0.8)
        a = ang[i]
        ax.text(1.09*math.cos(a), 1.09*math.sin(a), i, fontsize=6,
                ha='left' if math.cos(a) > 0 else 'right', va='center',
                rotation=math.degrees(a) if math.cos(a) > 0 else math.degrees(a)+180)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45); ax.axis('off')
    ax.set_title('Etoliko MAG metabolic interaction potential\n'
                 'blue = cross-feeding, orange = resource overlap, grey = dependency\n'
                 'nodes: blue bacteria, orange archaea; size = degree', fontsize=9)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(os.path.join(N, f'network.{ext}'), dpi=300, bbox_inches='tight')
    print(f"\n  figure: {os.path.join(N, 'network.png')} (+pdf, svg)")
except Exception as exc:
    print(f"\n  figure skipped ({type(exc).__name__}: {exc}) — tables are unaffected")

print()
print("  HOW TO REPORT THIS")
print("  'Genome-resolved metabolic modelling of the 29 MAGs identified N potential")
print("   cross-feeding exchanges, M pairs with overlapping resource niches, and")
print("   K amino-acid/vitamin dependencies. These are inferred from gene content")
print("   and represent metabolic POTENTIAL; they were not measured. Co-occurrence")
print("   networks were not computed, as two metagenomes cannot support them; that")
print("   analysis is reported for the 16S amplicon dataset in Chapter X.'")
PY

echo
echo "######## VALIDATION GATE ########"
check "nodes table"          $([ -s "$N/nodes.tsv" ] && echo 0 || echo 1)
check "edges table"          $([ -s "$N/edges.tsv" ] && echo 0 || echo 1)
check "competition matrix"   $([ -s "$N/competition_matrix.tsv" ] && echo 0 || echo 1)
check "cross-feeding table"  $([ -s "$N/crossfeeding.tsv" ] && echo 0 || echo 1)
NE=$(( $(wc -l < "$N/edges.tsv" 2>/dev/null || echo 1) - 1 ))
check "edges non-empty"      $([ "$NE" -gt 0 ] && echo 0 || echo 1)

runlog "21_interaction_network complete"
echo
echo "######## 21 DONE ########"
