#!/usr/bin/env bash
# =============================================================================
# 22_singlem.sh  --  SingleM: how much of the community did we actually get?
#
# WHY SINGLEM IS DIFFERENT FROM KRAKEN2 AND SILVA
#   Kraken2 answers "which sequenced genomes are here" — blind to anything
#   without a genome. SILVA (13b) answers "which rRNA lineages are here" —
#   sensitive to uncultured taxa but in rRNA units, with copy-number bias.
#
#   SingleM does something neither can: it extracts a short conserved WINDOW
#   from ~60 universal single-copy ribosomal protein genes directly out of
#   unassembled reads, and clusters those windows into OTUs. Because the markers
#   are single-copy, an OTU count is a GENOME count — no rRNA copy-number bias —
#   and because it needs no reference genome for the organism, it sees lineages
#   that have never been sequenced. It is currently the best available estimate
#   of prokaryotic community composition from raw reads.
#
# THE ONE THING THIS IS REALLY FOR: `singlem appraise`
#   Step 08d put MAG recovery at ~2% (sediment) and ~38% (cultures) from read
#   mapping with CoverM. That estimate is circular in a specific way: reads only
#   map to genomes that were recovered, so the denominator is the assembly's
#   view of the community, and anything that neither assembled nor binned is
#   invisible to it.
#
#   `singlem appraise` computes the same quantity a completely different way. It
#   takes the OTU windows found in the READS, the OTU windows found in the MAGs,
#   and asks what fraction of read-derived windows have a match in the MAG set.
#   Nothing is assembled, nothing is mapped. An INDEPENDENT measurement of the
#   headline number in the chapter is worth having, and if the two agree the
#   ~2% figure becomes very hard to argue with.
#
#   It also names what is missing: `--output-unaccounted-for-otu-table` lists
#   the community members with no MAG, at genome-level resolution.
#
# ONT CAVEAT — STATE THIS
#   SingleM was designed for short reads. It still works on long reads (windows
#   are extracted per-ORF), but two things change: (1) a long read can carry
#   several markers, so read counts are not directly comparable to short-read
#   studies, and (2) ONT indels can corrupt a window and lose it. Both effects
#   reduce sensitivity rather than invent taxa, so the PRESENCE calls and the
#   appraise fraction are more trustworthy than the absolute abundances.
#   Report composition as relative, and cross-check against 13b SILVA.
#
# OUTPUT   22_singlem/
#     {SED,LIQ}_otu_table.tsv          community OTUs from reads
#     {SED,LIQ}_profile.tsv            taxonomic profile
#     mags_otu_table.tsv               OTUs from the 29 MAGs
#     appraise/                        THE KEY RESULT: % of community in MAGs
#     unaccounted_for.tsv              named lineages with no MAG
#
# RUNTIME  2-6 h on the read sets. RUN IN tmux.
# USAGE    ws2 && bash scripts/22_singlem.sh 2>&1 | tee logs/22_singlem.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="22_singlem"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
G="$WS2/22_singlem"; F="$WS2/07_mags"
mkdir -p "$G/appraise"

echo "######## 22 — SINGLEM COMMUNITY PROFILE ########"
echo "  started: $(date -Is)"
runlog "22_singlem started"

# ---------------------------------------------------------------- INSTALL
# SingleM is not in any environment on this machine. Build a dedicated one on
# ${SCRATCH_ROOT} — never into a shared env, and never onto / (133 GB free).
SM_ENV="${SINGLEM_ENV:-${ENVS_DIR}/singlem_env}"
SM_BIN="$SM_ENV/bin/singlem"
if [ ! -x "$SM_BIN" ]; then
  echo "  SingleM not found at $SM_BIN"
  # The first attempt failed with 'nothing provides __win needed by tqdm' and a
  # cascade of biopython/zlib conflicts. That is not a SingleM problem: the
  # `defaults` channel was in the solve and, without strict priority, the solver
  # mixed packages built against incompatible zlib/python and then reported the
  # symptom rather than the cause. --override-channels + --strict-channel-priority
  # with conda-forge FIRST is the fix, and pinning python keeps the solve small.
  echo "  installing (conda-forge first, strict priority, no defaults) ..."
  SOLVER=conda
  command -v mamba >/dev/null 2>&1 && SOLVER=mamba
  echo "  solver: $SOLVER"
  $SOLVER create -y -p "$SM_ENV" \
      --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda \
      "python=3.10" singlem \
    || {
      echo "  !! strict solve failed too. Trying without the python pin ..."
      $SOLVER create -y -p "$SM_ENV" \
          --override-channels --strict-channel-priority \
          -c conda-forge -c bioconda singlem \
        || {
          echo
          echo "  !! conda cannot solve SingleM on this machine."
          echo "     Do NOT retry the plain command — it will fail the same way."
          echo "     Options, in order of preference:"
          echo "       1. install mamba first, then re-run this script:"
          echo "            conda install -y -n base -c conda-forge mamba"
          echo "       2. pixi/micromamba, which use the same channels:"
          echo "            micromamba create -y -p $SM_ENV -c conda-forge -c bioconda singlem"
          echo "       3. skip SingleM. It is a CROSS-CHECK on the ~2%/~38% recovery"
          echo "          figures, not a dependency of any other step. Nothing else"
          echo "          in the pipeline is blocked by this."
          exit 1
        }
    }
fi
[ -x "$SM_BIN" ] || { echo "  !! still no singlem binary at $SM_BIN"; exit 1; }
export PATH="$SM_ENV/bin:$PATH"
echo "  singlem: $("$SM_BIN" --version 2>&1 | head -1)"
record_version "singlem" "$SM_ENV" "$SM_BIN" --version 2>/dev/null || true

# ---------------------------------------------------------------- METAPACKAGE
# The reference metapackage is a separate multi-GB download. Without it every
# command fails with an unhelpful error, so resolve it explicitly here.
MP="${SINGLEM_METAPACKAGE_PATH:-}"
if [ -z "$MP" ] || [ ! -e "$MP" ]; then
  MPDIR="${DB_DIR}/singlem"
  mkdir -p "$MPDIR"
  MP=$(ls -d "$MPDIR"/*.smpkg* 2>/dev/null | head -1 || true)
  if [ -z "$MP" ]; then
    echo "  metapackage not found — downloading into $MPDIR"
    echo "  (several GB; resumable. If this fails, see"
    echo "   https://wwood.github.io/singlem/ for the current download command.)"
    "$SM_BIN" data --output-directory "$MPDIR" \
      || { echo "  !! metapackage download failed."
           echo "     Fetch it manually into $MPDIR, then re-run with"
           echo "     SINGLEM_METAPACKAGE_PATH=/path/to/metapackage"; exit 1; }
    MP=$(ls -d "$MPDIR"/*.smpkg* 2>/dev/null | head -1 || true)
  fi
fi
[ -n "$MP" ] && [ -e "$MP" ] || { echo "  !! no metapackage resolved"; exit 1; }
export SINGLEM_METAPACKAGE_PATH="$MP"
echo "  metapackage: $MP"

# ---------------------------------------------------------------- READS
echo
echo "######## A. OTUs FROM THE UNASSEMBLED READS ########"
pipe_reads() {  # pipe_reads <CODE> <reads>
  local code="$1" reads="$2"
  local otu="$G/${code}_otu_table.tsv" prof="$G/${code}_profile.tsv"
  [ -s "$reads" ] || { echo "  !! missing $reads"; return 1; }
  if [ -s "$otu" ]; then echo "  $code: OTU table already present"; else
    echo "  $code: singlem pipe ($(date +%H:%M)) ..."
    "$SM_BIN" pipe --forward "$reads" --otu-table "$otu" \
        --threads "$THREADS" --metapackage "$MP" 2>&1 | sed 's/^/    /'
  fi
  if [ ! -s "$prof" ]; then
    # taxonomic profile is a separate output in current SingleM; older versions
    # only emit the OTU table, so a failure here is not fatal.
    "$SM_BIN" pipe --forward "$reads" --taxonomic-profile "$prof" \
        --threads "$THREADS" --metapackage "$MP" 2>/dev/null \
      || echo "    (no --taxonomic-profile in this version; OTU table is enough)"
  fi
  [ -s "$otu" ] && echo "    OTU rows: $(( $(wc -l < "$otu") - 1 ))"
}
pipe_reads SED "$SED_READS"
pipe_reads LIQ "$LIQ_READS"

# ---------------------------------------------------------------- MAGs
echo
echo "######## B. OTUs FROM THE 29 MAGs ########"
MAGOTU="$G/mags_otu_table.tsv"
if [ -s "$MAGOTU" ]; then
  echo "  already present"
else
  NB=$(ls -1 "$F"/bins/*.fa 2>/dev/null | wc -l)
  echo "  running singlem pipe over $NB MAG FASTAs ..."
  [ "$NB" -gt 0 ] || { echo "  !! no bins at $F/bins"; exit 1; }
  "$SM_BIN" pipe --genome-fasta-files "$F"/bins/*.fa \
      --otu-table "$MAGOTU" --threads "$THREADS" --metapackage "$MP" 2>&1 \
    | sed 's/^/    /'
fi
[ -s "$MAGOTU" ] && echo "  MAG OTU rows: $(( $(wc -l < "$MAGOTU") - 1 ))"

# ---------------------------------------------------------------- APPRAISE
echo
echo "######## C. APPRAISE — WHAT FRACTION OF THE COMMUNITY IS IN THE MAGs ########"
echo "  This is the independent check on the ~2% / ~38% figures from step 08d."
for code in SED LIQ; do
  otu="$G/${code}_otu_table.tsv"
  [ -s "$otu" ] || continue
  out="$G/appraise/${code}_appraise.txt"
  if [ -s "$out" ]; then echo "  $code: already appraised"; continue; fi
  echo "  $code:"
  "$SM_BIN" appraise \
      --metagenome-otu-tables "$otu" \
      --genome-otu-tables "$MAGOTU" \
      --output-unaccounted-for-otu-table "$G/appraise/${code}_unaccounted.tsv" \
      2>&1 | tee "$out" | sed 's/^/    /' \
    || echo "    !! appraise failed for $code"
  # --output-found-in is not a flag in this SingleM version; passing it made
  # appraise exit with "unrecognized arguments" and print its usage message.
  # The old validation gate then reported PASS because the .txt file existed.
  if grep -qE '^usage: singlem|unrecognized arguments|error:' "$out" 2>/dev/null; then
    echo "    !! appraise did NOT run — see $out"
  fi
done

# ===========================================================================
python3 - "$G" "$F" "$WS2" <<'PY'
import sys, os, csv, glob
from collections import Counter, defaultdict

G, F, WS2 = sys.argv[1], sys.argv[2], sys.argv[3]

def read_otu(p):
    if not os.path.exists(p): return []
    with open(p) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        return list(rd)

print()
print("=" * 100)
print("  SINGLEM SUMMARY")
print("=" * 100)

mag_rows = read_otu(os.path.join(G, 'mags_otu_table.tsv'))
mag_seqs = {r.get('sequence') for r in mag_rows if r.get('sequence')}
print(f"  distinct marker windows in the 29 MAGs: {len(mag_seqs)}")

for code in ('SED', 'LIQ'):
    rows = read_otu(os.path.join(G, f'{code}_otu_table.tsv'))
    if not rows:
        print(f"\n  {code}: no OTU table"); continue
    print("\n" + "-" * 100)
    print(f"  {code}")
    print("-" * 100)

    # coverage is SingleM's per-OTU estimate of genome equivalents
    def cov(r):
        for k in ('coverage', 'num_hits'):
            try: return float(r.get(k) or 0)
            except ValueError: pass
        return 0.0

    tot = sum(cov(r) for r in rows) or 1
    by_tax = Counter()
    for r in rows:
        t = (r.get('taxonomy') or '').strip()
        by_tax[t] += cov(r)

    dom = Counter()
    for t, c in by_tax.items():
        parts = [x.strip() for x in (t or '').split(';') if x.strip()]
        d = next((x[3:] for x in parts if x.startswith('d__')), '')
        if not d:
            d = next((x for x in parts if x.lower() not in ('root',)), '') or 'unassigned'
        dom[d] += c
    # SingleM taxonomy strings start with "Root; d__Bacteria; ...", so taking
    # field 0 returned "Root" for every row and printed "Root 100.0%".
    print("  domain composition (by marker coverage): " + ', '.join(
        f"{k} {100*v/tot:.1f}%" for k, v in dom.most_common(5)))

    print(f"\n  TOP 20 LINEAGES")
    print(f"  {'%':>7}   taxonomy")
    for t, c in by_tax.most_common(20):
        short = ' '.join(t.split(';')[-2:]).strip() or 'unassigned'
        print(f"  {100*c/tot:>6.2f}%   {short[:82]}")

    # ------------------------------------------------ appraise, parsed
    ap = os.path.join(G, 'appraise', f'{code}_appraise.txt')
    unacc = os.path.join(G, 'appraise', f'{code}_unaccounted.tsv')
    if os.path.exists(ap):
        print(f"\n  APPRAISE — raw output:")
        for line in open(ap):
            if line.strip(): print("    " + line.rstrip())
    # independent recovery estimate straight from the windows
    seqs_cov = defaultdict(float)
    for r in rows:
        s = r.get('sequence')
        if s: seqs_cov[s] += cov(r)
    inmag = sum(v for s, v in seqs_cov.items() if s in mag_seqs)
    allc = sum(seqs_cov.values()) or 1
    print(f"\n  WINDOW-LEVEL RECOVERY (exact-match, computed here independently)")
    print(f"    {100*inmag/allc:.2f}% of marker coverage in this sample corresponds to a")
    print(f"    window that is also present in the 29-MAG catalogue.")
    print("    NOTE: exact-window matching is STRICTER than singlem appraise, which")
    print("    allows near-matches. Treat this as a lower bound and the appraise")
    print("    percentage above as the figure to quote.")

    if os.path.exists(unacc):
        urows = read_otu(unacc)
        ub = Counter()
        for r in urows:
            t = (r.get('taxonomy') or '').strip()
            try: ub[t] += float(r.get('coverage') or 0)
            except ValueError: pass
        utot = sum(ub.values()) or 1
        print(f"\n  TOP 15 LINEAGES WITH NO MAG  ({len(urows)} unaccounted-for windows)")
        print(f"  {'%unacc':>8}   taxonomy")
        for t, c in ub.most_common(15):
            short = ' '.join(t.split(';')[-2:]).strip() or 'unassigned'
            print(f"  {100*c/utot:>7.2f}%   {short[:82]}")

print()
print("=" * 100)
print("  CROSS-CHECK AGAINST THE OTHER TWO METHODS")
print("=" * 100)
print("  Three independent estimates of MAG recovery now exist:")
print("    CoverM read mapping (08d)   sediment ~2%    cultures ~38%")
print("    SingleM appraise (here)     see above")
print("    SILVA genus overlap (13c)   see 13c output")
print()
print("  They measure different things and will NOT agree exactly — read mapping")
print("  counts reads, SingleM counts genome equivalents, SILVA counts rRNA. What")
print("  matters is whether they agree on the ORDER OF MAGNITUDE. If all three say")
print("  'a few percent of the sediment', that number is solid and no reviewer will")
print("  move it. If SingleM says 20%, something in the read-mapping estimate needs")
print("  explaining and it is much better to find that now.")
print()
print("  METHODS SENTENCE: 'Community composition was additionally profiled with")
print("  SingleM, which recovers conserved windows from universal single-copy")
print("  ribosomal marker genes directly from unassembled reads and is therefore")
print("  independent of both assembly and reference genome availability. The")
print("  fraction of the community represented by the MAG catalogue was estimated")
print("  with singlem appraise. Because SingleM was developed for short reads,")
print("  Nanopore-derived abundances are reported as relative values.'")
PY

echo
echo "######## VALIDATION GATE ########"
for c in SED LIQ; do
  check "$c OTU table" $([ -s "$G/${c}_otu_table.tsv" ] && echo 0 || echo 1)
done
check "MAG OTU table"  $([ -s "$G/mags_otu_table.tsv" ] && echo 0 || echo 1)
# check SUCCESS, not file existence — a usage/error message is also a file
APOK=1
for c in SED LIQ; do
  f="$G/appraise/${c}_appraise.txt"
  [ -s "$f" ] && ! grep -qE '^usage: singlem|unrecognized arguments|error:' "$f" && APOK=0
done
check "appraise SUCCEEDED (not just wrote a file)" "$APOK"

runlog "22_singlem complete"
echo
echo "######## 22 DONE ########"
