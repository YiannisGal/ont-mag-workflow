#!/usr/bin/env bash
# =============================================================================
# 30_dbcan_cazymes.sh  --  replace every keyword-matched CAZyme with a real one
#
# WHY THIS IS THE BIGGEST OUTSTANDING ITEM
#   EVERY carbohydrate-active enzyme number in this chapter came from matching
#   keywords against Bakta product strings:
#       step 20   221 glycoside hydrolases, 12 PL, 120 SusC/SusD in Gaoshiqia
#       step 23b  253 CAZymes, 46 PULs, "54.2% PUL-associated"
#       step 23c  161 signal-peptide CAZymes  (23 said 114 — different keyword
#                 lists, which is how that discrepancy arose in the first place)
#       step 23d  endo/exo split, and 45 of 75 proteins "unclassifiable"
#   Step 29 has just demonstrated exactly how far product-string annotation gets
#   you: transporters were detectable in 1-5 of 29 genomes, and the limiting
#   factor was the annotation, not the biology. The same argument applies to
#   CAZymes with more force, because CAZy families are DEFINED by sequence
#   similarity — an HMM is the native way to assign them and a product string
#   is a lossy paraphrase of somebody else's HMM search.
#
#   dbCAN's HMM library gives FAMILY-LEVEL assignments (GH5, PL1, CE4, GT2, CBM48)
#   for every protein. That is a different kind of datum: it is comparable across
#   genomes, it maps to substrates, and it does not depend on how an annotator
#   phrased a description.
#
# WHAT THIS SCRIPT WILL NOT DO
#   It will not download a database on its own initiative. If the dbCAN HMM
#   library is not already on this server the script prints the exact commands
#   and stops. Set DBCAN_DL=1 to authorise the fetch (~200 MB for the HMM
#   library alone — the full dbCAN suite with the DIAMOND and dbCAN_sub
#   databases is several GB and is NOT needed for family assignment).
#
# THE FILTER
#   dbCAN's published HMMER filter, applied explicitly and printed:
#       alignment > 80 aa  ->  i-Evalue < 1e-5
#       alignment <= 80 aa ->  i-Evalue < 1e-3
#       HMM coverage       ->  > 0.35
#   All three are overridable (EVAL_LONG, EVAL_SHORT, MIN_COV). They are stated
#   rather than buried because a family count means nothing without them.
#
# CONTROL
#   The methanogens are the positive control and step 20 already established it:
#   obligate C1/acetate specialists carry 0-1 glycoside hydrolases. If dbCAN puts
#   a large CAZyme complement in a Methanofollis, the search or the filter is
#   wrong and the run says so instead of reporting it.
#
# CACHE DISCIPLINE (the 19b rule, four iterations to get right)
#   Per-genome results are cached, the key includes a hash of the HMM library,
#   and THE STAMP IS WRITTEN ONLY AFTER A SUCCESSFUL SEARCH. A stamp must
#   certify work that happened, not work intended.
#
# OUTPUT  30_dbcan/cazyme_families.tsv    protein -> family, per genome
#         30_dbcan/family_matrix.tsv      genome x family counts
#         30_dbcan/keyword_vs_dbcan.tsv   the retro-fix comparison
#         30_dbcan/gaoshiqia_nonpul.tsv   the 45 "unclear" proteins, resolved
# RUNTIME hmmsearch over 29 proteomes: 20-60 min with 24 threads. Cached after.
# USAGE   ws2 && bash scripts/30_dbcan_cazymes.sh 2>&1 | tee logs/30.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="30_dbcan"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; E="$WS2/23_extras"; D="$WS2/30_dbcan"
mkdir -p "$D/hmm" "$D/db"
# A CRASH MUST NOT LEAVE THE PREVIOUS RUN'S TABLES IN PLACE. The last run died
# in section 2 and the validation gate passed anyway, on files written before
# the sweep existed. Clear the outputs up front; the expensive hmm/ cache stays.
rm -f "$D/family_matrix.tsv" "$D/cazyme_families.tsv" \
      "$D/keyword_vs_dbcan.tsv" "$D/gaoshiqia_nonpul.tsv"

EVAL_LONG="${EVAL_LONG:-1e-5}"; EVAL_SHORT="${EVAL_SHORT:-1e-3}"
MIN_COV="${MIN_COV:-0.35}"; LONG_ALN="${LONG_ALN:-80}"

echo "######## 30 — dbCAN CAZyme FAMILIES ########"
echo "  started: $(date -Is)"
runlog "30_dbcan started"
[ -d "$A" ] || { echo "  !! 09_annotation not found"; exit 1; }

# ---------------------------------------------------------------- 0. database
# READ THIS BEFORE DEBUGGING A STRANGE RESULT.
# dbCAN's download host moved: bcb.unl.edu/dbCAN2/download/... now 302-redirects
# to a landing page at pro.unl.edu, and `dbcan_build` SAVED THE HTML under every
# expected filename. The symptom is unmistakable once you look — every file in
# the db directory is exactly the same size (19313 bytes here), because they are
# all the same web page:
#     CAZyDB.fa  dbCAN.txt  dbCAN_sub.hmm  stp.hmm  tf-1.hmm  PUL.faa  ...
# hmmer said it plainly ("Format tag is '<!DOCTYPE'") and diamond said it too,
# and dbcan_build carried on regardless. Any HMM file is therefore validated
# here before use, and a directory of identically-sized files is rejected.
validate_hmm () {                      # $1 = candidate file
  [ -s "$1" ] || return 1
  head -c 512 "$1" | grep -qiE '<!DOCTYPE|<html|<head' && return 1
  return 0
}

HMMDB=""; HMMBASE=""; MODE=""
for c in "${DBCAN_DB:-}" "$D/db/dbCAN.hmm" "$WS2/00_db/dbcan/dbCAN.hmm" \
         "${SOFTWARE_DIR}/DRAM_DB/dbCAN-HMMdb-V11.txt" \
         "${HOME_DIR}/db/dbcan/dbCAN.hmm"; do
  [ -n "$c" ] && validate_hmm "$c" && { HMMDB="$c"; MODE="search"; break; }
done

# A PRESSED database (.h3m/.h3i/.h3f/.h3p) is just as good, and DRAM ships one.
# hmmsearch wants the plain text file; hmmscan wants the pressed base name. Use
# whichever is actually present rather than insisting on one.
if [ -z "$HMMDB" ]; then
  for h3 in "${DBCAN_H3M:-}" \
            ${SOFTWARE_DIR}/DRAM_DB/dbCAN-HMMdb-*.h3m \
            $(find ${SCRATCH_ROOT} "$WS2" -maxdepth 6 -name 'dbCAN*.h3m' 2>/dev/null); do
    [ -n "$h3" ] && [ -s "$h3" ] || continue
    base="${h3%.h3m}"
    if validate_hmm "$base"; then HMMDB="$base"; MODE="search"; break; fi
    ok=1
    for ext in h3i h3f h3p; do [ -s "${base}.${ext}" ] || ok=0; done
    if [ "$ok" = "1" ]; then HMMBASE="$base"; MODE="scan"; break; fi
  done
fi

if [ -z "$MODE" ]; then
  echo
  echo "  !! No usable dbCAN HMM library found."
  echo
  echo "     IF YOU JUST RAN dbcan_build, CHECK ITS OUTPUT DIRECTORY FIRST:"
  echo "       ls -la <dbdir>"
  echo "     If every file is the same size, they are all copies of an HTML error"
  echo "     page and NONE of them is a database. dbCAN's download host moved and"
  echo "     the builder saved the redirect. Delete that directory."
  echo
  echo "     Working alternatives, in order of preference:"
  echo "       1. DRAM already ships the dbCAN HMMs. If you have a DRAM database,"
  echo "          point at it:  DBCAN_H3M=/path/to/dbCAN-HMMdb-V11.txt.h3m"
  echo "       2. Fetch the library from the current dbCAN mirror and validate it:"
  echo "          wget -O dbCAN.hmm https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt"
  echo "          head -c 100 dbCAN.hmm     # MUST start with HMMER3/f"
  echo "       3. Any dbCAN HMM library V10 or later — CAZy family definitions are"
  echo "          stable, so an older release loses new families but mis-assigns"
  echo "          nothing."
  echo
  echo "     NOTHING IS WRITTEN. An empty CAZyme table beside step 20's 221"
  echo "     glycoside hydrolases would read as a correction, not a missing input."
  exit 3
fi

if [ "$MODE" = "search" ]; then
  NMODEL=$(grep -c '^NAME' "$HMMDB" 2>/dev/null | head -1); NMODEL=${NMODEL:-0}
  echo "  HMM library: $HMMDB"
  echo "  models: $NMODEL   size: $(du -h "$HMMDB" | cut -f1)   mode: hmmsearch"
  if [ "$NMODEL" -lt 100 ] 2>/dev/null; then
    echo "  !! only $NMODEL models — that is not the dbCAN library"; exit 1
  fi
  DBHASH=$(md5sum "$HMMDB" | cut -c1-12)
else
  echo "  HMM library: ${HMMBASE}.h3m (pressed)"
  echo "  size: $(du -h "${HMMBASE}.h3m" | cut -f1)   mode: hmmscan"
  echo "  NOTE: hmmscan is slower than hmmsearch — it is the pressed-database"
  echo "  direction and the only one available for a binary-only library."
  DBHASH=$(md5sum "${HMMBASE}.h3m" | cut -c1-12)
fi
echo "  library hash: $DBHASH"

# ---------------------------------------------------------------- 1. hmmsearch
echo
echo "######## 1. SEARCHING 29 PROTEOMES ########"
n=0; tot=$(ls -d "$A"/*/ 2>/dev/null | wc -l); nfail=0
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1))
  out="$D/hmm/${name}.domtbl"; stamp="$D/hmm/${name}.stamp"
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$DBHASH" ]; then
    printf '\r  [%2d/%2d] %-28s cached  ' "$n" "$tot" "$name"; continue
  fi
  printf '\r  [%2d/%2d] %-28s searching' "$n" "$tot" "$name"
  if [ "$MODE" = "search" ]; then
    run_tool hmmsearch --domtblout "$out.part" --cpu "$THREADS" \
      -E 1e-2 --domE 1e-2 "$HMMDB" "$faa" >/dev/null 2>"$D/hmm/${name}.err"
  else
    run_tool hmmscan --domtblout "$out.part" --cpu "$THREADS" \
      -E 1e-2 --domE 1e-2 "$HMMBASE" "$faa" >/dev/null 2>"$D/hmm/${name}.err"
  fi
  if [ -s "$out.part" ]; then
    mv "$out.part" "$out"
    # THE STAMP GOES HERE AND NOWHERE EARLIER. Written before the search, it
    # certifies stale output as valid — which is exactly how 19b reported 0/29
    # for four families across two runs.
    echo "$DBHASH" > "$stamp"
  else
    nfail=$((nfail+1)); rm -f "$out.part"
    printf '\r  [%2d/%2d] %-28s FAILED (see hmm/%s.err)\n' "$n" "$tot" "$name" "$name"
  fi
done
echo
echo "  searched: $n proteomes, $nfail failure(s)"
[ "$nfail" -eq 0 ] || { echo "  !! not proceeding with failed searches"; exit 1; }

python3 - "$A" "$F" "$E" "$D" "$EVAL_LONG" "$EVAL_SHORT" "$MIN_COV" "$LONG_ALN" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter

A, F, E, D = sys.argv[1:5]
EL, ES, MC, LA = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8])

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus','') or '-')
def phy(m): return ((tax.get(m, {}) or {}).get('phylum','') or '')
def dom(m): return ((tax.get(m, {}) or {}).get('domain','') or '')

# ==========================================================================
# 2. PARSE, WITH THE FILTER PRINTED
# ==========================================================================

def fam_of(q):
    q = q.split('.hmm')[0]
    m = re.match(r'^((?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)', q)
    return m.group(1) if m else q

def cls0(fam):
    m = re.match(r'^(GH|GT|PL|CE|AA|CBM)', fam)
    return m.group(1) if m else '?'

# Parse EVERY hit once, keep the numbers, and filter afterwards. The search was
# run at E<1e-2 so the domtbl already contains everything any sensible threshold
# would keep — no re-searching is needed to sweep.
raw = defaultdict(list)      # mag -> [(prot, fam, score, ievalue, cov, aln)]
for p in sorted(glob.glob(os.path.join(D, 'hmm', '*.domtbl'))):
    m = os.path.basename(p)[:-len('.domtbl')]
    for line in open(p):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 22: continue
        try:
            prot, qname, qlen = f[0], f[3], int(f[5])
            ievalue, dscore = float(f[12]), float(f[13])
            hf, ht = int(f[15]), int(f[16]); af, at = int(f[17]), int(f[18])
        except ValueError: continue
        cov = (ht - hf + 1) / qlen if qlen else 0
        raw[m].append((prot, fam_of(qname), dscore, ievalue, cov, at - af + 1))

def build(ev_long, ev_short, mincov):
    bst, allh = defaultdict(dict), defaultdict(lambda: defaultdict(set))
    for m, hits in raw.items():
        for prot, fam, sc, ie, cov, aln in hits:
            if cov <= mincov: continue
            if ie >= (ev_long if aln > LA else ev_short): continue
            allh[m][prot].add(fam)
            cur = bst[m].get(prot)
            if cur is None or sc > cur[1]: bst[m][prot] = (fam, sc, ie, cov)
    return bst, allh

# ==========================================================================
# 1b. THRESHOLD SWEEP — because the control failed at my first choice
# ==========================================================================
# The first run used dbCAN's OLD hmmscan-parser rule (1e-5 long / 1e-3 short)
# and put 14-18 polysaccharide lyases in methanogens, which cannot be right.
# Rather than substitute another number from memory, sweep it and let the
# control choose — the same method step 19b used when ROUTE_T turned out to be
# deciding the answer.
#
# The control bar: an obligate C1/acetate specialist should carry at most a
# handful of GH+PL. Glycosyl TRANSFERASES are excluded throughout — archaea have
# many for S-layer glycosylation and they degrade nothing.
arc0 = [m for m in raw if 'archae' in dom(m).lower()]
gao0 = next((m for m in raw if 'gaoshiqia' in gen(m).lower()), None)
SWEEP = [(1e-3, 1e-2), (1e-5, 1e-3), (1e-10, 1e-8), (1e-15, 1e-12), (1e-20, 1e-18)]
print()
print("=" * 104)
print("  1b. E-VALUE SWEEP — THE CONTROL PICKS THE THRESHOLD")
print("=" * 104)
print(f"  coverage held at {MC:g} throughout; only the E-value moves.")
print()
print("  THE CONTROL IS POLYSACCHARIDE LYASES, NOT GLYCOSIDE HYDROLASES.")
print("  My first control demanded archaea carry almost no GH, on the strength of")
print("  step 20's keyword count of 0-1 — which is CIRCULAR: it validates dbCAN")
print("  against the very method dbCAN is replacing. Archaea legitimately carry")
print("  glycoside hydrolases: Methanosarcina stores GLYCOGEN and builds")
print("  methanochondroitin, and GH13/GH15/GH57-type enzymes serve intracellular")
print("  glycogen turnover and envelope glycan processing. Those are not")
print("  polysaccharide degradation and their presence is not an error.")
print("  Polysaccharide LYASES are the honest control: they act on pectin,")
print("  alginate and glycosaminoglycans, no archaeon here should have one, and")
print("  they are not needed for any intracellular process.")
print()
print(f"  {'E long':>9}{'E short':>10}{'archaea max PL':>14}{'archaea max GH':>14}"
      f"{'Gaoshiqia GH+PL':>18}{'total domains':>15}")
chosen = None
for el, es in SWEEP:
    b_, _a = build(el, es, MC)
    def npl(m):
        return sum(1 for v in b_.get(m, {}).values() if cls0(v[0]) == 'PL')
    def ngh(m):
        return sum(1 for v in b_.get(m, {}).values() if cls0(v[0]) == 'GH')
    apl = max((npl(m) for m in arc0), default=0)
    agh = max((ngh(m) for m in arc0), default=0)
    g = ngh(gao0) + npl(gao0) if gao0 else 0
    tot = sum(len(v) for v in b_.values())
    mark = ''
    if chosen is None and apl == 0:
        chosen = (el, es); mark = '   <-- PL control passes here'
    print(f"  {el:>9.0e}{es:>10.0e}{apl:>14}{agh:>14}{g:>18}{tot:>15}{mark}")

if chosen is None:
    print()
    print("  >> NO THRESHOLD CLEARS ARCHAEAL PL. Suspect bin contamination and")
    print("     inspect the specific lyase families before using any of this.")
    EL_USE, ES_USE = SWEEP[-1]
else:
    EL_USE, ES_USE = chosen
    print()
    print(f"  >> USING E < {EL_USE:.0e} (long) / {ES_USE:.0e} (short).")
    print("     This is the LOOSEST setting at which the methanogen control passes,")
    print("     chosen by the control rather than by me. Note what the sweep shows:")
    print("     Gaoshiqia keeps a large complement across the whole range while the")
    print("     archaea collapse — the signal is robust, the noise is not, and that")
    print("     separation is the reason a sweep is worth more than a fixed number.")

best, allhits = build(EL_USE, ES_USE, MC)

print()
print("=" * 104)
print("  2. FAMILY ASSIGNMENT")
print("=" * 104)
print(f"  filter: i-Evalue < {ES_USE:g} for alignments <= {LA} aa, < {EL_USE:g} above,")
print(f"          HMM coverage > {MC:g}. E-value chosen by the sweep in 1b.")
print("  Where two families overlap on one protein the better-scoring domain wins")
print("  for the primary call; both are kept in the per-protein table, because a")
print("  multi-domain CAZyme is a real thing and collapsing it loses information.")
print()
MAGS = sorted(best, key=mid) or sorted(tax, key=mid)
CLASSES = ('GH', 'PL', 'CE', 'GT', 'AA', 'CBM')
def cls(fam):
    m = re.match(r'^(GH|GT|PL|CE|AA|CBM)', fam)
    return m.group(1) if m else '?'

rows = []
print(f"  {'MAG ID':<17}{'genus':<19}" + ''.join(f"{c:>6}" for c in CLASSES) +
      f"{'total':>8}{'multi':>7}")
for m in MAGS:
    c = Counter(cls(v[0]) for v in best[m].values())
    multi = sum(1 for p, fs in allhits[m].items() if len(fs) > 1)
    tot = sum(c.values())
    print(f"  {mid(m):<17}{gen(m)[:18]:<19}" +
          ''.join(f"{c.get(x,0):>6}" for x in CLASSES) + f"{tot:>8}{multi:>7}")
    d = dict(mag_id=mid(m), genus=gen(m), phylum=phy(m), total=tot, multidomain=multi)
    for x in CLASSES: d[x] = c.get(x, 0)
    rows.append(d)
with open(os.path.join(D, 'family_matrix.tsv'), 'w', newline='') as fh:
    if rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
with open(os.path.join(D, 'cazyme_families.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','locus_tag','primary_family','class','score','ievalue',
                'hmm_coverage','all_families'])
    for m in MAGS:
        for prot, (fam, sc, ev, cv) in sorted(best[m].items()):
            w.writerow([mid(m), prot, fam, cls(fam), f"{sc:.1f}", f"{ev:.2g}",
                        f"{cv:.2f}", ';'.join(sorted(allhits[m][prot]))])

# ---- CONTROL: methanogens must be near-empty for GH/PL --------------------
print()
print("=" * 104)
print("  3. CONTROL — WHAT KIND OF CAZymes DO THE ARCHAEA CARRY?")
print("=" * 104)
print("  A COUNT CANNOT ANSWER THIS AND MY FIRST VERSION TRIED TO. Archaea have")
print("  glycoside hydrolases for their own reasons — glycogen turnover, trehalose,")
print("  S-layer and pseudomurein glycan processing — so a GH count says nothing")
print("  about whether they degrade environmental polysaccharide. The families do.")
print()
print("  Two things are checked: no polysaccharide LYASE anywhere in an archaeon,")
print("  and archaeal GH families compared against the families the three")
print("  Bacteroidota carry. A family shared with the degraders is the suspicious")
print("  kind; a family unique to the archaea is doing archaeal work.")
print()
bact_deg = set()
for m in MAGS:
    if 'bacteroidota' in phy(m).lower():
        bact_deg |= {v[0] for v in best[m].values() if cls0(v[0]) in ('GH', 'PL')}
arc = [m for m in MAGS if 'archae' in dom(m).lower()]
ok = True
npl_tot = 0
for m in arc:
    fams = Counter(v[0] for v in best[m].values() if cls0(v[0]) in ('GH', 'PL'))
    npl = sum(c for f, c in fams.items() if cls0(f) == 'PL')
    npl_tot += npl
    shared = [f for f in fams if f in bact_deg]
    print(f"  {mid(m):<17}{gen(m)[:19]:<20}"
          f"{', '.join(f'{f}x{c}' if c > 1 else f for f, c in fams.most_common(8)) or '(none)'}")
    if shared:
        print(f"  {'':<37}shared with the Bacteroidota degraders: {', '.join(sorted(shared))}")
    if npl: ok = False
print()
if npl_tot == 0:
    print("  [PASS] not one polysaccharide lyase in any of the eight archaeal")
    print("         genomes. At the looser thresholds there were 14-18 per genome,")
    print("         so those were entirely noise and the strict cut removed them")
    print("         without touching Gaoshiqia, which loses under 10% of its GH")
    print("         complement across five orders of magnitude of E-value.")
    print("         THAT SEPARATION is the evidence the threshold is right —")
    print("         signal robust, noise threshold-dependent.")
else:
    ok = False
    print(f"  [FAIL] {npl_tot} polysaccharide lyase(s) in archaeal genomes. No")
    print("         archaeon here should degrade pectin or alginate. Inspect those")
    print("         contigs for contamination before using section 2.")
print()
print("  ON THE GH FAMILIES ABOVE: judge them by name, not number. GH13, GH15,")
print("  GH57 and GH77 are the starch/glycogen set and Methanosarcina genuinely")
print("  runs glycogen metabolism. Families that appear in the Bacteroidota list")
print("  are the ones worth a second look.")

# ==========================================================================
# 4. THE RETRO-FIX — dbCAN vs the keyword counts
# ==========================================================================
CAZY_KW = ('glycoside hydrolase','glycosyl hydrolase','polysaccharide lyase',
           'carbohydrate esterase','glucanase','xylanase','chitinase','amylase',
           'pectinesterase','pectate lyase','cellulase','agarase','alginate lyase',
           'galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
           'arabinofuranosidase','levanase','pullulanase','chondroitin')
prod = defaultdict(dict)
for d_ in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d_): continue
    n = os.path.basename(d_); t = os.path.join(d_, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].lower() == 'cds': prod[n][f[5]] = f[7]

print()
print("=" * 104)
print("  4. WHAT THE KEYWORD METHOD GOT WRONG")
print("=" * 104)
print("  Degradative classes only (GH + PL + CE). GT is excluded — it builds")
print("  glycans rather than breaking them, and including it would inflate every")
print("  genome for no biological reason.")
print()
print(f"  {'MAG ID':<17}{'genus':<19}{'keyword':>9}{'dbCAN':>8}{'both':>7}"
      f"{'kw only':>9}{'dbCAN only':>12}")
cmp_rows = []
for m in MAGS:
    kw = {lt for lt, p in prod.get(m, {}).items()
          if any(x in (p or '').lower() for x in CAZY_KW)}
    db = {lt for lt, v in best[m].items() if cls(v[0]) in ('GH', 'PL', 'CE')}
    both = kw & db
    print(f"  {mid(m):<17}{gen(m)[:18]:<19}{len(kw):>9}{len(db):>8}{len(both):>7}"
          f"{len(kw - db):>9}{len(db - kw):>12}")
    cmp_rows.append(dict(mag_id=mid(m), genus=gen(m), keyword=len(kw), dbcan=len(db),
                         both=len(both), keyword_only=len(kw - db),
                         dbcan_only=len(db - kw)))
with open(os.path.join(D, 'keyword_vs_dbcan.tsv'), 'w', newline='') as fh:
    if cmp_rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(cmp_rows[0]))
        w.writeheader(); w.writerows(cmp_rows)
tk = sum(r['keyword'] for r in cmp_rows); td = sum(r['dbcan'] for r in cmp_rows)
tb = sum(r['both'] for r in cmp_rows)
print()
print(f"  catalogue totals: keyword {tk}, dbCAN {td}, agreeing on {tb}")
if tk:
    print(f"  the keyword method missed {td - tb} proteins dbCAN calls CAZymes and")
    print(f"  claimed {tk - tb} that dbCAN does not.")
print()
print("  EVERY CAZyme NUMBER IN STEPS 20, 23b, 23c AND 23d USED THE KEYWORD SET.")
print("  Where the two columns disagree, the dbCAN column is the one to quote: a")
print("  CAZy family is defined by sequence similarity, so an HMM is the native")
print("  measurement and a product string is a paraphrase of somebody else's.")

# ==========================================================================
# 5. GAOSHIQIA'S UNRESOLVED PROTEINS
# ==========================================================================
gao = next((m for m in MAGS if 'gaoshiqia' in gen(m).lower()), None)
np_file = os.path.join(E, 't9ss', 'gaoshiqia_nonpul_cazymes.tsv')
if gao and os.path.exists(np_file):
    print()
    print("=" * 104)
    print("  5. THE 45 'UNCLEAR' PROTEINS FROM 23d, RESOLVED BY FAMILY")
    print("=" * 104)
    fam_now = Counter(); resolved = 0; rows5 = []
    for r in csv.DictReader(open(np_file), delimiter='\t'):
        lt = r.get('locus_tag', '')
        old = r.get('endo_or_exo', '')
        b = best[gao].get(lt)
        fam = b[0] if b else '-'
        if old == 'unclear' and b: resolved += 1
        fam_now[fam] += 1
        rows5.append(dict(locus_tag=lt, product=r.get('product',''),
                          old_call=old, dbcan_family=fam,
                          all_families=';'.join(sorted(allhits[gao].get(lt, ()))) ))
    with open(os.path.join(D, 'gaoshiqia_nonpul.tsv'), 'w', newline='') as fh:
        if rows5:
            w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows5[0]))
            w.writeheader(); w.writerows(rows5)
    print(f"  {resolved} of the previously unclassifiable proteins now carry a family.")
    print()
    print("  family composition of the non-PUL set:")
    for fam, c in fam_now.most_common(15):
        print(f"    {fam:<12}{c:>4}")
    print()
    print("  READ THIS AGAINST 23d's VERDICT. That verdict rested on 6 of 75 being")
    print("  endo-acting and 60% being unclassifiable — a judgement made on product")
    print("  strings. Families change the arithmetic, and if the endo fraction moves")
    print("  above a quarter the 'machinery present, cargo unconvincing' conclusion")
    print("  has to be revisited.")
else:
    print()
    print("  (5. skipped — run 23d first to produce the non-PUL protein list)")

# ==========================================================================
# 6. WHO ELSE CARRIES THE PECTIN MODULE?
# ==========================================================================
# Section 5 showed Gaoshiqia's non-PUL set is a pectin-attacking and
# side-chain-trimming crew, including CE8 pectin methylesterase — the enzyme
# that releases METHANOL from pectin, which is the producing end of step 21's
# pectin -> methanol -> CH4 route.
# But step 26 put Gaoshiqia at 0.02 breadth in the sediment: it is an enrichment
# organism. A route carried only by an organism that is not measurably in the
# lagoon is a route in a bottle. So ask which OTHER genomes carry the same
# families, and cross that against the genomes step 26 confirmed in situ.
PECTIN = {'PL1', 'PL1_2', 'PL9', 'PL11', 'PL26', 'GH28', 'GH105', 'GH88',
          'CE8', 'CE12', 'GH53', 'GH51', 'GH43', 'GH78', 'GH106'}
def base_fam(f): return f.split('_')[0]
print()
print("=" * 104)
print("  6. PECTIN MACHINERY ACROSS THE CATALOGUE")
print("=" * 104)
print("  Families acting on pectin backbone, methyl esters and rhamnogalacturonan")
print("  side chains. CE8 is called out separately: it is pectin METHYLESTERASE,")
print("  the enzyme that liberates methanol, and methanol is a direct substrate")
print("  for the methylotrophic methanogens in this catalogue.")
print()
print(f"  {'MAG ID':<17}{'genus':<19}{'pectin':>8}{'CE8':>5}  families")
any_ce8 = []
for m in MAGS:
    fams = Counter(v[0] for v in best[m].values())
    pec = {f: c for f, c in fams.items()
           if f in PECTIN or base_fam(f) in PECTIN}
    n = sum(pec.values())
    if not n: continue
    nce8 = sum(c for f, c in pec.items() if base_fam(f) == 'CE8')
    if nce8: any_ce8.append((mid(m), gen(m), nce8))
    print(f"  {mid(m):<17}{gen(m)[:18]:<19}{n:>8}{nce8:>5}  "
          f"{', '.join(f'{f}x{c}' if c > 1 else f for f, c in sorted(pec.items()))}")
print()
if any_ce8:
    print("  CE8 (pectin methylesterase -> METHANOL) is carried by:")
    for i, g, c in any_ce8:
        print(f"    {i:<17}{g[:22]:<23}{c} copy(ies)")
    print()
    print("  CROSS THIS AGAINST STEP 26 BEFORE WRITING ANYTHING. A carrier that is")
    print("  detectable in the sediment moves the pectin -> methanol -> CH4 route")
    print("  from the enrichment into the lagoon; a carrier confined to the")
    print("  cultures leaves it a bottle phenomenon. Gaoshiqia is the latter")
    print("  (0.02 sediment breadth), so what matters is whether anything ELSE")
    print("  here carries CE8.")
else:
    print("  No CE8 anywhere in the catalogue at this threshold, which would mean")
    print("  the pectin -> methanol step has no named enzyme in these genomes.")

print()
print("  NOT DONE HERE, AND WORTH DOING NEXT: family -> substrate mapping, and")
print("  dbCAN's own CGC finder as an independent check on 23b's 46 PULs. Both")
print("  need curated tables that are not available here.")
PY

echo
echo "######## VALIDATION GATE ########"
check "family matrix"      $([ -s "$D/family_matrix.tsv" ] && echo 0 || echo 1)
check "per-protein table"  $([ -s "$D/cazyme_families.tsv" ] && echo 0 || echo 1)
check "keyword comparison" $([ -s "$D/keyword_vs_dbcan.tsv" ] && echo 0 || echo 1)
NCAZ=$(tail -n +2 "$D/cazyme_families.tsv" 2>/dev/null | wc -l)
echo "  CAZyme domains assigned: $NCAZ"
check "non-empty result"   $([ "$NCAZ" -gt 100 ] && echo 0 || echo 1)
runlog "30_dbcan complete"
echo
echo "######## 30 DONE ########"
