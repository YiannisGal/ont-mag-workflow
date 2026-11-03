#!/usr/bin/env bash
# =============================================================================
# 35_kaiju_protein_taxonomy.sh  --  classify the reads in PROTEIN space
#
# THE PROBLEM
#   Kraken2 against RefSeq standard-16 classified **0.47% of sediment reads**
#   (step 13, cleaned set). That is not a failure of Kraken2; it is what
#   nucleotide k-mer matching does to a community where 22 of 29 MAGs have no
#   GTDB species match and 81% of sediment genera have no cultured
#   representative (step 24). A 31-mer has to match almost exactly, and across
#   the divergence separating these lineages from anything sequenced, nothing
#   does.
#
#   **Kaiju translates the reads and searches PROTEIN space.** Third codon
#   positions are free to drift, so protein identity survives divergence that
#   destroys nucleotide identity — and it survives ONT indels better, because
#   greedy mode tolerates substitutions inside a match.
#
#   Step 24 costed this and never ran it.
#
# WHAT THIS CHANGES AND WHAT IT DOES NOT
#   It does NOT replace SingleM. SingleM counts marker windows and reports
#   GENOME EQUIVALENTS, the right currency for "how much of the community did we
#   recover" (2.40% sediment, step 22). Kaiju classifies READS — a different
#   quantity that inflates abundant organisms, the same distinction step 22 drew
#   between 38% of reads and 26% of genome equivalents.
#
#   What Kaiju adds is BREADTH OF NAMING: which lineages are here at all.
#
# THE HEADLINE COMPARISON IS CONFOUNDED, AND THIS SCRIPT SAYS SO
#   Kraken2/RefSeq vs Kaiju/progenomes changes TWO things at once: nucleotide vs
#   protein, and one database vs another. The fold-difference is therefore a
#   statement about the two pipelines as run, not a clean measurement of what
#   translation buys. Section 1 prints it flagged as confounded.
#   The uncontaminated measurement is section 4, DEPTH OF ASSIGNMENT, which is
#   internal to Kaiju: how far down the ranks its assignments actually reach.
#
# CONTROLS, BOTH FROM OUR OWN DATA
#   POSITIVE: we assembled genomes for several genera. If Kaiju cannot find
#   them, its naming of everything else is not to be read. Catalogue genera that
#   are GTDB placeholders (UBA*, alphanumeric codes) are EXCLUDED from this
#   control — they have no NCBI name to find, so their absence is nomenclature,
#   not failure. GTDB suffixed names (Clostridium_A) are tested on their base.
#   NEGATIVE: a euxinic sediment must not show a large eukaryote/aerobe
#   fraction, and no single taxon should dominate. Step 13d found 84.65%
#   *Erwinia* in this sample and it was a pooled 16S amplicon run.
#
# OUTPUT  35_kaiju/{SED,LIQ}.kaiju · *_summary_<rank>.tsv · VERDICT.txt
# RUNTIME 20-60 min per sample once the database exists. First build 1-3 h.
# USAGE   ws2 && bash scripts/35_kaiju_protein_taxonomy.sh 2>&1 | tee logs/35.log
#         KAIJU_DB=/path/to/kaiju_db.fmi   point at a specific index
#         SKIP_PROBE=1                     skip the 20k-read pre-flight
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="35_kaiju"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"

F="$WS2/07_mags"; D="$WS2/35_kaiju"; mkdir -p "$D"
DBDIR="${KAIJU_DBDIR:-${DB_DIR}/kaiju}"
MINLEN="${MINLEN:-11}"       # minimum match length, amino acids (Kaiju default)
MODE="${MODE:-greedy}"       # greedy tolerates substitutions; mem does not
MISMATCH="${MISMATCH:-5}"    # allowed substitutions in greedy mode

echo "######## 35 — KAIJU: TAXONOMY IN PROTEIN SPACE ########"
echo "  started: $(date -Is)"
runlog "35_kaiju started"

# ------------------------------------------------------------------ binaries
# Resolve by LOOKING FOR THE BINARY, never by assuming ENV_<tool> is set.
# This bug (run_tool: "no working env resolved for X") has now cost three steps
# — hmmscan in 33, prodigal in 34 — because ENV_<tool> only exists for tools the
# smoketest happened to probe. A newly installed tool is never in env_resolved.sh.
bin_env() {
  local tool="$1"
  local var="ENV_${tool//-/_}"
  local cands="${!var:-} ${ENV_kaiju:-} ${ENV_bakta:-} -"
  local e prefix root d
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then
      command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }
      continue
    fi
    if declare -f env_prefix >/dev/null 2>&1; then
      prefix=$(env_prefix "$e" 2>/dev/null) || continue
      [ -x "$prefix/bin/$tool" ] && { echo "$e"; return 0; }
    fi
  done
  # last resort: scan every conda root for an env that actually holds the binary.
  # This finds a hand-made env that env_resolved.sh has never heard of.
  for root in ${CONDA_ROOTS:-}; do
    [ -d "$root" ] || continue
    for d in "$root"/*; do
      [ -x "$d/bin/$tool" ] && { echo "$(basename "$d")"; return 0; }
    done
  done
  return 1
}
for t in kaiju kaiju2table; do
  e=$(bin_env "$t" || true)
  if [ -n "$e" ]; then export "ENV_${t}=$e"; echo "  $t: env '$e'"
  else
    echo "  !! $t not found in PATH or under: ${CONDA_ROOTS:-<unset>}"
    echo "     Install it into a writable prefix — NOT into a shared env:"
    echo "       conda create -y -p ${CONDA_ENVS_USER}/kaiju -c bioconda -c conda-forge kaiju"
    runlog "**35 ABORTED: kaiju not available**"
    exit 2
  fi
done
record_version "kaiju" "${ENV_kaiju:--}" kaiju -h 2>/dev/null || true

# ------------------------------------------------------------------ database
echo
echo "######## 0. DATABASE ########"
FMI=""; NODES=""; NAMES=""
for cand in "${KAIJU_DB:-}" "$DBDIR"/*.fmi "$DBDIR"/*/*.fmi; do
  [ -n "${cand:-}" ] && [ -s "$cand" ] && { FMI="$cand"; break; }
done
if [ -n "$FMI" ]; then
  base=$(dirname "$FMI")
  for n in "$base/nodes.dmp" "$DBDIR/nodes.dmp" "$base/../nodes.dmp"; do
    [ -s "$n" ] && { NODES="$n"; break; }; done
  for n in "$base/names.dmp" "$DBDIR/names.dmp" "$base/../names.dmp"; do
    [ -s "$n" ] && { NAMES="$n"; break; }; done
fi
if [ -z "$FMI" ] || [ -z "$NODES" ] || [ -z "$NAMES" ]; then
  echo "  !! no COMPLETE Kaiju database found under $DBDIR"
  echo "     found: index=${FMI:-none}  nodes=${NODES:-none}  names=${NAMES:-none}"
  DFP="$DBDIR"; [ -d "$DFP" ] || DFP=$(dirname "$DBDIR"); [ -d "$DFP" ] || DFP=${SCRATCH_ROOT}
  AV=$(df -BG --output=avail "$DFP" 2>/dev/null | tail -1 | tr -dc '0-9'); AV=${AV:-0}
  echo
  echo "     free space at $DFP: ${AV} GB    (progenomes: ~35 GB peak, ~13 GB kept)"
  echo
  echo "     mkdir -p $DBDIR && cd $DBDIR"
  echo "     kaiju-makedb -s progenomes -t $THREADS"
  echo
  echo "     progenomes is the right trade: representative-genome based, so it is"
  echo "     fast and its taxonomy is clean. nr_euk (~65 GB) sees more but is"
  echo "     dominated by well-studied organisms — the same bias that gave"
  echo "     Kraken2 0.47%."
  echo
  echo "     NOTHING IS WRITTEN. A run against a missing database classifies"
  echo "     nothing and would read as 'this community is unclassifiable', which"
  echo "     is precisely the claim this step exists to TEST."
  runlog "**35 ABORTED: no Kaiju database**"
  exit 3
fi
IDXGB=$(du -m "$FMI" | cut -f1); IDXGB=$(( IDXGB / 1024 ))
echo "  index : $FMI  (${IDXGB} GB)"
echo "  nodes : $NODES"
echo "  names : $NAMES"
echo "  mode  : $MODE (minimum match $MINLEN aa, $MISMATCH mismatches allowed)"

# Kaiju holds the FM-index in RAM. Step 32 lost three runs to an OOM kill that
# took four wrong diagnoses to even locate, so check BEFORE launching.
AVGB=$(free -g 2>/dev/null | awk '/^Mem:/{print $7}'); AVGB=${AVGB:-0}
echo "  memory: ${AVGB} GB available, index needs ~${IDXGB} GB resident"
if [ "$AVGB" -gt 0 ] && [ "$IDXGB" -gt 0 ] && [ "$AVGB" -lt $(( IDXGB + 4 )) ]; then
  echo "  !! not enough free memory to hold this index plus working space."
  echo "     Free memory or use a smaller database. Refusing to start a run that"
  echo "     will be OOM-killed after an hour."
  runlog "**35 ABORTED: insufficient memory for index**"
  exit 4
fi
DBTAG=$(basename "$FMI" .fmi)
printf '%s\t%s\t%s\t%s\t%s\n' "$STEP" "kaiju_db" "$DBTAG" "$FMI" "$(date -Is)" >> "$VERSIONS"

# ------------------------------------------------------------------- reads in
CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
SAMESET=1
if [ -s "$CLEAN" ]; then
  SED_IN="$CLEAN"
  echo "  sediment reads: cleaned set (RBK_test amplicon run excluded, step 13f)"
else
  SED_IN="$SED_READS"
  SAMESET=0
  echo "  !! using the UNCLEANED sediment reads — the amplicon run is still in."
  echo "     The Kraken2 comparison below is then across DIFFERENT read sets and"
  echo "     will be marked invalid. Run 13f with WRITE_CLEAN=1 first."
fi

# ------------------------------------------------------------------- pre-flight
# 20,000 reads through the real database before committing an hour to it.
# An empty or truncated index classifies nothing and looks exactly like a novel
# community — the failure mode this step is supposed to distinguish from.
if [ "${SKIP_PROBE:-0}" != "1" ] && [ ! -s "$D/SED.kaiju" ]; then
  echo
  echo "######## 0b. PRE-FLIGHT (20,000 reads) ########"
  PB="$TMPDIR/k35probe.fq"
  if [ "${SED_IN##*.}" = "gz" ]; then
    ( zcat "$SED_IN" 2>/dev/null | head -n 80000 > "$PB" ) || true
  else
    ( head -n 80000 "$SED_IN" > "$PB" ) || true
  fi
  if [ ! -s "$PB" ]; then
    echo "  !! could not read $SED_IN"; runlog "**35 ABORTED: unreadable reads**"; exit 5
  fi
  run_tool kaiju -t "$NODES" -f "$FMI" -i "$PB" -o "$PB.out" \
      -z "$THREADS" -a "$MODE" -m "$MINLEN" -e "$MISMATCH" \
      > "$D/probe.log" 2>&1
  prc=$?
  pn=$(wc -l < "$PB.out" 2>/dev/null || echo 0)
  pc=$(grep -c '^C' "$PB.out" 2>/dev/null || echo 0)
  echo "  exit $prc · $pn reads processed · $pc classified"
  if [ "$prc" -ne 0 ] || [ "$pn" -eq 0 ]; then
    echo "  !! probe failed. Last lines of $D/probe.log:"
    tail -8 "$D/probe.log" | sed 's/^/     /'
    echo
    echo "     If the complaint is about gzip: this Kaiju build cannot read .gz."
    echo "     Decompress once and set SED_READS/LIQ_READS to the plain files."
    runlog "**35 ABORTED: kaiju probe failed**"; exit 6
  fi
  if [ "$pc" -eq 0 ]; then
    echo "  !! 20,000 reads, ZERO classified. That is a broken or empty index, not"
    echo "     a novel community — even a euxinic sediment hits ribosomal proteins."
    echo "     Rebuild the database before reading anything from this step."
    runlog "**35 ABORTED: probe classified nothing**"; exit 7
  fi
  echo "  [PASS] index is live (probe rate $(awk -v c="$pc" -v n="$pn" 'BEGIN{printf "%.2f", 100*c/n}')%)"
  rm -f "$PB" "$PB.out"
fi

# ---------------------------------------------------------------------- run
echo
echo "######## 1. CLASSIFYING ########"
for pair in "SED:$SED_IN" "LIQ:$LIQ_READS"; do
  code="${pair%%:*}"; reads="${pair#*:}"
  out="$D/${code}.kaiju"
  [ -s "$reads" ] || { echo "  !! missing reads for $code"; continue; }
  if [ -s "$out" ]; then echo "  $code: cached ($(wc -l < "$out") reads)"; continue; fi
  echo "  $code: kaiju ($(date +%H:%M)) ..."
  run_tool kaiju -t "$NODES" -f "$FMI" -i "$reads" -o "$out.part" \
      -z "$THREADS" -a "$MODE" -m "$MINLEN" -e "$MISMATCH" \
      > "$D/${code}.kaiju.log" 2>&1
  rc=$?
  # exit 0 is not success — step 33 accepted 29 zero exits that wrote nothing.
  if [ "$rc" -ne 0 ] || [ ! -s "$out.part" ]; then
    echo "    !! kaiju exited $rc and/or wrote nothing. Last lines:"
    tail -6 "$D/${code}.kaiju.log" | sed 's/^/       /'
    rm -f "$out.part"; continue
  fi
  mv "$out.part" "$out"
  echo "    $(wc -l < "$out") reads processed"
done

# ------------------------------------------------------------------- tables
echo
echo "######## 2. SUMMARY TABLES ########"
for code in SED LIQ; do
  [ -s "$D/${code}.kaiju" ] || continue
  for rank in phylum class order family genus species; do
    t="$D/${code}_summary_${rank}.tsv"
    # a header-only leftover is not a cached table
    [ -s "$t" ] && [ "$(wc -l < "$t")" -gt 1 ] && continue
    run_tool kaiju2table -t "$NODES" -n "$NAMES" -r "$rank" \
        -o "$t" "$D/${code}.kaiju" >/dev/null 2>&1 \
      || echo "  !! kaiju2table failed for $code $rank"
  done
  echo "  $code: $(ls -1 "$D/${code}"_summary_*.tsv 2>/dev/null | wc -l) rank tables"
done

# --------------------------------------------------------------- interpret
echo
python3 - "$D" "$F" "$WS2" "$SAMESET" "$DBTAG" <<'PY'
import sys, os, re, csv

D, F, WS2, SAMESET, DBTAG = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == '1', sys.argv[5]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

PSEUDO = ('unclassified', 'cannot be assigned', 'viruses')

def load_table(code, rank):
    """kaiju2table columns: file, percent, reads, taxon_id, taxon_name.
    Returns {name: (percent, reads)} for REAL taxa only, plus the pseudo rows."""
    p = os.path.join(D, f'{code}_summary_{rank}.tsv')
    real, pseudo = {}, {}
    if not os.path.exists(p):
        return real, pseudo
    with open(p, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            name = (r.get('taxon_name') or '').strip().split(';')[-1].strip()
            if not name:
                continue
            try:
                pct = float(r.get('percent') or 0)
            except ValueError:
                pct = 0.0
            try:
                rds = int(float(r.get('reads') or 0))
            except ValueError:
                rds = 0
            low = name.lower()
            (pseudo if any(k in low for k in PSEUDO) else real)[name] = (pct, rds)
    return real, pseudo

# ------------------------------------------------- classification rate
rates = {}
for code in ('SED', 'LIQ'):
    p = os.path.join(D, f'{code}.kaiju')
    if not os.path.exists(p):
        continue
    tot = cls = 0
    with open(p, errors='replace') as fh:
        for line in fh:
            tot += 1
            if line[:1] == 'C':
                cls += 1
    rates[code] = (cls, tot)

rule('1. HOW MUCH GOT CLASSIFIED?   (confounded — read section 4 for the clean number)')
K2 = {'SED': 0.47, 'LIQ': 3.61}
say(f"  {'sample':<8}{'classified':>13}{'total':>13}{'rate':>9}{'Kraken2':>10}{'ratio':>9}")
say('-' * W)
for code, (c, t) in rates.items():
    r = 100.0 * c / t if t else 0.0
    k = K2.get(code)
    if k:
        say(f"  {code:<8}{c:>13,}{t:>13,}{r:>8.2f}%{k:>9.2f}%{r / k:>8.0f}x")
    else:
        say(f"  {code:<8}{c:>13,}{t:>13,}{r:>8.2f}%{'-':>10}{'-':>9}")
say()
if not SAMESET:
    say("  >> INVALID COMPARISON: this run used the UNCLEANED sediment reads while")
    say("     the Kraken2 figure is from the cleaned set. Different denominators.")
say("  CONFOUND, stated plainly: Kraken2 ran against RefSeq standard-16 and Kaiju")
say(f"  against {DBTAG}. Two things differ — nucleotide vs protein, and database vs")
say("  database — so this ratio describes the two pipelines as run. It is not a")
say("  clean estimate of what translation alone buys, and must not be written up")
say("  as one. The uncontaminated measurement is section 4.")

# -------------------------------------------------------- positive control
rule('2. POSITIVE CONTROL — DOES KAIJU FIND THE ORGANISMS WE ASSEMBLED?')

def genera_from_catalogue(path):
    """Find genus names however MAG_TAXONOMY.tsv happens to be shaped.
    32c read a column name that 32b never wrote; do not repeat that."""
    out, how = set(), 'none'
    if not os.path.exists(path):
        return out, 'missing'
    with open(path, errors='replace') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if not rows:
        return out, 'empty'
    cols = list(rows[0].keys())
    gcol = next((c for c in cols if c and c.strip().lower() == 'genus'), None)
    if gcol:
        for r in rows:
            v = (r.get(gcol) or '').strip()
            if v and v not in ('-', 'NA', 'None'):
                out.add(v)
        how = f"column '{gcol}'"
    if not out:  # fall back to pulling g__ out of any lineage string anywhere
        for r in rows:
            for v in r.values():
                for m in re.findall(r'g__([^;\t\r\n]+)', str(v or '')):
                    m = m.strip()
                    if m:
                        out.add(m)
        how = 'g__ parsed from lineage strings'
    return out, how

tax = os.path.join(F, 'MAG_TAXONOMY.tsv')
known, how = genera_from_catalogue(tax)
say(f"  source: {tax}")
say(f"  genus names read via: {how}   ({len(known)} distinct)")

def ncbi_keys(g):
    """GTDB name -> the names NCBI might use. Returns () if untestable."""
    g = g.strip().lstrip('g__')
    if not g or re.search(r'\d', g):        # UBA1234, JAAXXX01 — GTDB placeholders
        return ()
    if not re.fullmatch(r'[A-Z][a-z]+(_[A-Z])?', g):
        return ()
    base = re.sub(r'_[A-Z]$', '', g)        # Clostridium_A -> Clostridium
    return (g, base) if base != g else (g,)

testable = {g: ncbi_keys(g) for g in known}
untestable = sorted(g for g, k in testable.items() if not k)
testable = {g: k for g, k in testable.items() if k}

sed_gen, _ = load_table('SED', 'genus')
liq_gen, _ = load_table('LIQ', 'genus')

def lookup(keys, table):
    for k in keys:
        if k in table:
            return table[k][0]
    return None

hits = {g: lookup(k, sed_gen) for g, k in testable.items()}
hits_liq = {g: lookup(k, liq_gen) for g, k in testable.items()}
found = sorted((g for g, v in hits.items() if v is not None), key=lambda g: -hits[g])
missing = sorted(g for g, v in hits.items() if v is None
                 and hits_liq.get(g) is None)
liq_only = sorted(g for g, v in hits.items() if v is None
                  and hits_liq.get(g) is not None)

say(f"  testable against NCBI taxonomy: {len(testable)}"
    f"    GTDB placeholders excluded: {len(untestable)}")
if untestable:
    say(f"    excluded: {', '.join(untestable[:10])}"
        + (' ...' if len(untestable) > 10 else ''))
    say("    (these have no NCBI name for Kaiju to return; their absence is")
    say("     nomenclature, not a failure of the classifier.)")
say()
say(f"  {'catalogue genus':<28}{'SED %':>10}{'LIQ %':>10}")
say('-' * W)
for g in found[:20]:
    lv = hits_liq.get(g)
    say(f"  {g:<28}{hits[g]:>9.3f}%" + (f"{lv:>9.3f}%" if lv is not None else f"{'-':>10}"))
for g in liq_only[:10]:
    say(f"  {g:<28}{'-':>10}{hits_liq[g]:>9.3f}%")
say()
n_ok = len(found) + len(liq_only)
say(f"  recovered {n_ok}/{len(testable)} testable catalogue genera "
    f"({len(found)} in sediment, {len(liq_only)} only in the cultures)")

if not testable:
    CTRL = 'NOCTRL'
elif n_ok >= max(3, 0.5 * len(testable)):
    CTRL = 'PASS'
elif n_ok >= 2:
    CTRL = 'WEAK'
else:
    CTRL = 'FAIL'

if CTRL == 'PASS':
    say("  [PASS] Kaiju recovers organisms we independently assembled from these")
    say("         same samples, so its naming of the rest is worth reading.")
elif CTRL == 'WEAK':
    say("  [WEAK] Only a minority of testable catalogue genera came back. Read the")
    say("         classification as suggestive of BREADTH only; do not quote any")
    say("         individual genus percentage without a second line of evidence.")
elif CTRL == 'NOCTRL':
    say("  [NO CONTROL] No testable genus names could be read from MAG_TAXONOMY.tsv.")
    say("         This is a missing-input problem, not a Kaiju result. Check that")
    say("         file exists and carries either a 'genus' column or GTDB lineages.")
else:
    say("  [FAIL] Kaiju misses the genera whose genomes we hold. Either the")
    say("         database lacks these lineages or the run is broken. NO VERDICT —")
    say("         do not read the classification of anything else.")
if missing:
    say(f"  not found anywhere: {', '.join(missing[:12])}")
    say("  (a testable genus missing here is usually absent from the DATABASE:")
    say("   progenomes carries representative genomes only, and recently named")
    say("   genera may postdate the release.)")

# ------------------------------------------------- the negative control
rule('3. NEGATIVE CONTROL — DOES IT LOOK LIKE A EUXINIC SEDIMENT?')
sed_phy, _ = load_table('SED', 'phylum')
top = sorted(sed_phy.items(), key=lambda x: -x[1][0])[:12]
say(f"  {'phylum':<40}{'sediment %':>12}")
say('-' * W)
for n, (p, _r) in top:
    say(f"  {n[:39]:<40}{p:>11.3f}%")
say()
EUK = ('streptophyta', 'chordata', 'arthropoda', 'ascomycota', 'basidiomycota',
       'mollusca', 'nematoda', 'apicomplexa')
euk = sum(p for n, (p, _r) in sed_phy.items() if n.lower() in EUK)
sed_gen_all = sed_gen
dom_name, dom_pct = ('-', 0.0)
if sed_gen_all:
    dom_name, (dom_pct, _r) = max(sed_gen_all.items(), key=lambda x: x[1][0])
say(f"  obvious eukaryote phyla : {euk:.2f}%   (expect low single digits at most)")
say(f"  most abundant genus     : {dom_name} at {dom_pct:.2f}%")
NEG_OK = euk < 15.0 and dom_pct < 50.0
if NEG_OK:
    say("  [PASS] no eukaryote bloom, no single-taxon takeover.")
else:
    say("  [FLAG] this is the step-13d signature. That step found 84.65% *Erwinia*")
    say("         in this sample and the cause was a 16S amplicon run pooled into")
    say("         the shotgun set. Establish where this signal comes from before")
    say("         reporting anything else from this table.")

# ---------------------------------------------- depth of assignment
rule('4. DEPTH OF ASSIGNMENT — the measurement with no database confound')
say("  Kaiju assigns each read to the LCA of its best protein matches. A read")
say("  from an organism with a close relative in the database lands at genus or")
say("  species; a read from a lineage with nothing near it gets pulled up to")
say("  phylum or falls out entirely. The SHAPE of that decay is a property of")
say("  this community's novelty, measured entirely within one tool and one")
say("  database — so unlike section 1 it is not comparing two things at once.")
say()
say(f"  {'rank':<10}{'SED taxa':>10}{'SED % reads':>14}{'LIQ taxa':>10}{'LIQ % reads':>14}")
say('-' * W)
depth = {}
for rank in ('phylum', 'class', 'order', 'family', 'genus', 'species'):
    row = {}
    for code in ('SED', 'LIQ'):
        real, _ps = load_table(code, rank)
        row[code] = (len(real), sum(p for p, _r in real.values()))
    depth[rank] = row
    say(f"  {rank:<10}{row['SED'][0]:>10}{row['SED'][1]:>13.2f}%"
        f"{row['LIQ'][0]:>10}{row['LIQ'][1]:>13.2f}%")
say()
if depth.get('phylum') and depth['phylum']['SED'][1] > 0:
    ph = depth['phylum']['SED'][1]
    ge = depth['genus']['SED'][1]
    sp = depth['species']['SED'][1]
    say(f"  Sediment: {ph:.1f}% of reads get a phylum, {ge:.1f}% get a genus,"
        f" {sp:.1f}% a species.")
    say(f"  {ph - ge:.1f} percentage points of the community can be placed in a phylum")
    say("  and no further. That is the quantitative form of 'no cultured relative',")
    say("  and it is the number step 24 was reaching for when it reported 81% of")
    say("  sediment genera without a cultured representative.")

# ------------------------------------------------------------- breadth
rule('5. BREADTH vs SingleM (step 24) — DIFFERENT QUANTITIES, SHOWN TOGETHER')
say(f"  {'rank':<10}{'Kaiju SED':>11}{'Kaiju LIQ':>11}   SingleM SED / LIQ")
say('-' * W)
SM = {'phylum': (114, 65), 'class': (306, 167), 'order': (741, 347),
      'family': (1491, 622), 'genus': (3416, 1451), 'species': (769, 465)}
for rank in ('phylum', 'class', 'order', 'family', 'genus', 'species'):
    s = SM.get(rank, ('-', '-'))
    say(f"  {rank:<10}{depth[rank]['SED'][0]:>11}{depth[rank]['LIQ'][0]:>11}"
        f"   {s[0]} / {s[1]}")
say()
say("  These columns are NOT comparable as counts. SingleM counts marker windows")
say("  and estimates genome equivalents; Kaiju counts reads and reports every")
say("  taxon with at least one. Read-based counts over-represent abundant")
say("  organisms and produce a long tail of one-read taxa — the same distinction")
say("  step 22 drew between 38% of reads and 26% of genome equivalents. Show both")
say("  and label which is which; do not average them.")

rule('6. WHAT TO REPORT')
if CTRL == 'FAIL' or CTRL == 'NOCTRL':
    say("  >> NO VERDICT. The positive control did not pass; see section 2.")
elif not NEG_OK:
    say("  >> NO VERDICT. The negative control flagged; see section 3.")
else:
    if rates.get('SED'):
        c, t = rates['SED']
        r = 100.0 * c / t if t else 0.0
        say(f"  Kaiju classified {r:.2f}% of sediment reads in protein space"
            + (f" against {K2['SED']:.2f}% for Kraken2 in nucleotide space"
               if SAMESET else " (Kraken2 comparison invalid, see section 1)") + ".")
    say()
    say("  The defensible sentences:")
    say("   1. Protein-level classification places a large fraction of a community")
    say("      that nucleotide k-mer classification cannot. State the confound in")
    say("      the same sentence: the two pipelines differ in database as well as")
    say("      in search space.")
    say("   2. The rank-decay curve in section 4 is the citable novelty measure —")
    say("      it is internal to one tool and one database.")
    say(f"   3. Name the database in Methods: Kaiju {DBTAG}. Kaiju/progenomes and")
    say("      Kaiju/nr_euk are different measurements.")
    say("   4. Kaiju percentages are READ shares. For 'how much of the community")
    say("      did we recover', SingleM's genome equivalents remain the number.")
    if CTRL == 'WEAK':
        say("   5. The positive control was WEAK — breadth claims only, no per-genus")
        say("      percentages without corroboration.")
say()
say("  LIMITS:")
say("   1. LCA behaviour: a deep-rank name is evidence, a high-rank name is mostly")
say("      an admission that nothing close is in the database.")
say("   2. ONT reads are long. A 10 kb read spans many genes but Kaiju assigns it")
say("      from its single best-scoring match region, so the assignment rests on")
say("      one locus — which may be the most conserved one in the read. Long reads")
say("      raise sensitivity, not per-read specificity.")
say("   3. Read shares are not cell shares: genome size and rRNA content differ,")
say("      and step 27 showed how far rRNA-based shares can move.")
say("   4. progenomes covers representative genomes only. A genus absent from the")
say("      output may be absent from the database.")

open(os.path.join(D, 'VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'VERDICT.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "SED classified" "$([ -s "$D/SED.kaiju" ] && echo 0 || echo 1)" || OK=1
check "genus table written" "$([ -s "$D/SED_summary_genus.tsv" ] && echo 0 || echo 1)" || OK=1
check "verdict written" "$([ -s "$D/VERDICT.txt" ] && echo 0 || echo 1)" || OK=1
if grep -q 'NO VERDICT' "$D/VERDICT.txt" 2>/dev/null; then
  echo "  [STOP] a control did not pass — nothing in $D may be reported."
  OK=1
fi
[ "$OK" -eq 0 ] && runlog "35_kaiju complete" || runlog "**35_kaiju FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
