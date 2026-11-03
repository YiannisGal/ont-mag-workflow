#!/usr/bin/env bash
# =============================================================================
# 13_readprofile.sh  --  community composition WITHOUT binning
#
# WHY THIS IS THE BIGGEST GAP IN THE CHAPTER
#   Every result so far is genome-resolved, and the genomes cover
#       sediment  ~2%  of the community
#       cultures  ~38% of the community
#   So 98% of the sediment community has not been described at all. A reviewer
#   asking "what actually lives in this lagoon?" cannot be answered from 4 MAGs.
#   Read-based profiling is assembly-free and binning-free, so it sees the rare
#   majority that neither assembles nor bins.
#
# THREE COMPLEMENTARY VIEWS, deliberately not one
#   A. READS      Kraken2/Bracken on the unassembled reads. Sees everything the
#                 sequencer saw, including organisms too rare to assemble.
#   B. CONTIGS    the same classifier on assembled contigs, weighted by
#                 coverage, split into BINNED vs UNBINNED. This says directly
#                 what the MAGs missed — the question step 08d raised and could
#                 not answer.
#   C. rRNA       16S genes recovered from the whole assembly (not just MAGs),
#                 quantified by read coverage. This is the view that is directly
#                 comparable to the amplicon chapter, in the same units.
#
# MEMORY — AND THE FOUR DAYS THIS COST
#   The first run used maxikraken2_1903: a 140 GB database on a 125 GB host, so
#   --memory-mapping was set to read it from disk instead of loading it. That is
#   the documented workaround and it is technically correct. It is also
#   unusable: the run reached ~2% after FOUR DAYS. Memory mapping turns every
#   k-mer lookup into a random disk seek, and a metagenome does hundreds of
#   millions of them. The flag does not make a too-large database slow, it makes
#   it effectively non-terminating.
#
#   The fix is not a flag, it is a smaller database. standard-16 is capped at
#   16 GB, loads entirely into RAM, and finishes in hours. The cost is real and
#   must be stated in Methods: fewer reference genomes, so a lower classified
#   fraction. That trade is the right one here because Kraken2 is NOT this
#   project's rare-taxa method — SILVA SSU (13b) is, and it is built from rRNA
#   rather than genomes, so it sees uncultured lineages Kraken2 structurally
#   cannot. Kraken2's job is the broad, genome-backed overview.
#
#   This script now CHOOSES: it measures the database against available RAM and
#   memory-maps only if the database fits comfortably enough for that to be
#   sane. If the database is larger than RAM it refuses to start rather than
#   silently entering another four-day crawl.
#
# ONT CAVEAT
#   Kraken2 was designed for short reads. Long reads carry more k-mers each, so
#   classification rates are usually higher, but a chimeric or error-rich read
#   can pick up spurious assignments. --confidence 0.1 is applied to suppress
#   that, and the confidence value is reported so it can be stated in Methods.
#
# INPUT    00_input/{sediment,liquid}_reads.fastq.gz, 04_polish/*_asm_final.fasta
#          02_mapping/*_depth_multi.txt, 07_mags/bins/
# OUTPUT   13_profile/
#
# RUNTIME  2-5 h with a RAM-resident 16 GB database. RUN IN tmux.
# USAGE    ws2 && bash scripts/13_readprofile.sh 2>&1 | tee logs/13_profile_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13_readprofile"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
PR="$WS2/13_profile"; mkdir -p "$PR/kraken" "$PR/contigs" "$PR/rrna"

CONF="${KRAKEN_CONFIDENCE:-0.1}"
echo "######## 13 — READ-BASED COMMUNITY PROFILE ########"
echo "  started: $(date -Is)  threads: $THREADS  confidence: $CONF"
runlog "13_readprofile started (confidence=$CONF)"

[ -d "${KRAKEN2_DB:-/nonexistent}" ] || { echo "  !! KRAKEN2_DB not a directory: ${KRAKEN2_DB:-UNSET}"; exit 1; }
echo "  [db] $KRAKEN2_DB  ($(du -sh "$KRAKEN2_DB" 2>/dev/null | cut -f1))"
free -g | awk 'NR<=2{print "  RAM: "$0}'

# ---------------------------------------------------------------- MEMORY GATE
# Decide, do not assume. DB size and available RAM are both measured here.
# Measure ONLY the .k2d files — that is what Kraken2 loads. A directory total
# also counts the source tarball, the taxonomy dumps and the Bracken distribs,
# none of which are loaded, and would overstate the requirement by ~10 GB.
DBM=$(du -cm "$KRAKEN2_DB"/*.k2d 2>/dev/null | tail -1 | cut -f1)
[ -n "${DBM:-}" ] || DBM=$(du -sm "$KRAKEN2_DB" 2>/dev/null | cut -f1)
DBG=$(( ${DBM:-0} / 1024 ))
echo "  k2d index size: ${DBM:-?} MB  (directory total includes non-loaded files)"
AVAILG=$(awk '/MemAvailable/{printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || echo 0)
echo "  database ${DBG} GB   available RAM ${AVAILG} GB"
MMAP=""
if [ "$DBG" -gt "$AVAILG" ]; then
  echo
  echo "  !! THE DATABASE IS LARGER THAN AVAILABLE RAM (${DBG} GB > ${AVAILG} GB)."
  echo "     --memory-mapping would let this start and then never finish: the"
  echo "     previous attempt on a 140 GB database reached ~2% in four days."
  echo "     Refusing to repeat that."
  echo
  echo "     Use a database that fits. standard-16 is capped at 16 GB:"
  echo "       KRAKEN2_DB=${DB_DIR}/kraken2_standard16 \\"
  echo "         bash scripts/13_readprofile.sh"
  echo "     Or set FORCE_MMAP=1 to override, having read the paragraph above."
  [ "${FORCE_MMAP:-0}" = "1" ] || exit 1
  MMAP="--memory-mapping"
  echo "     FORCE_MMAP=1 set — proceeding with memory mapping against advice."
elif [ "$DBG" -gt $(( AVAILG / 2 )) ]; then
  echo "  database is large relative to RAM but fits — loading normally (no mmap)."
else
  echo "  database fits comfortably in RAM — loading normally (no mmap)."
fi
[ -n "$MMAP" ] && echo "  kraken2 memory mode: MEMORY-MAPPED" || echo "  kraken2 memory mode: LOADED"
record_version "kraken2" "${ENV_kraken2:-tax_profile}" kraken2 --version

# Bracken needs kmer-distribution files built for THIS database and read length.
BRACKEN_OK=0
# Bracken's -r is the read length its kmer distribution was built for. ONT reads
# are far longer than any available distribution, so none is a true match; the
# LARGEST available is the closest approximation and the choice must be stated
# in Methods rather than left at a default. Override with BRACKEN_READLEN.
BRK_R="${BRACKEN_READLEN:-}"
if [ -z "$BRK_R" ]; then
  BRK_R=$(ls -1 "$KRAKEN2_DB"/database*mers.kmer_distrib 2>/dev/null \
          | sed 's|.*/database||; s|mers.kmer_distrib||' | sort -n | tail -1)
fi
BRK_R="${BRK_R:-250}"
if ls "$KRAKEN2_DB"/database*mers.kmer_distrib >/dev/null 2>&1; then
  BRACKEN_OK=1
  echo "  Bracken kmer_distrib files present — abundance re-estimation available"
  echo "  Bracken read length: -r $BRK_R (largest distribution available; ONT reads"
  echo "    exceed every built distribution, so this is an approximation — state it)"
else
  echo "  NOTE: no Bracken kmer_distrib in the database."
  echo "        Kraken2 read counts will be reported instead. State this in Methods:"
  echo "        Kraken2 assigns each read to its LCA, so higher ranks accumulate"
  echo "        reads that could not be resolved further — it is NOT a relative"
  echo "        abundance estimate in the way Bracken's is."
fi

# =====================================================================
# A. READS
# =====================================================================
echo
echo "######## A. KRAKEN2 ON READS ($(date +%H:%M)) ########"
kraken_reads() {  # kraken_reads <CODE> <fastq.gz>
  local code="$1" reads="$2"
  local rep="$PR/kraken/${code}_reads.report" out="$PR/kraken/${code}_reads.kraken"
  [ -s "$reads" ] || { echo "  !! missing $reads"; return 1; }
  if [ -s "$rep" ]; then echo "  $code: already done"; return 0; fi
  echo "  $code: classifying $(basename "$reads") ..."
  run_tool kraken2 --db "$KRAKEN2_DB" $MMAP --threads "$THREADS" \
      --confidence "$CONF" --gzip-compressed \
      --report "$rep" --output "$out" "$reads" 2>&1 | sed 's/^/    /'
  # the per-read file is huge and is not needed downstream
  [ -s "$out" ] && gzip -f "$out" &
}
kraken_reads SED "$SED_READS"
kraken_reads LIQ "$LIQ_READS"
wait

if [ "$BRACKEN_OK" -eq 1 ]; then
  echo
  echo "  Bracken re-estimation ..."
  for c in SED LIQ; do
    for lvl in P C O F G S; do
      [ -s "$PR/kraken/${c}_reads.report" ] || continue
      # bracken was never registered in env_resolved.sh, so run_tool had no
      # prefix to invoke and every call died with "no working env resolved".
      # It lives alongside kraken2 in the tax_profile env.
      ENV_bracken="${ENV_bracken:-${ENV_kraken2:-tax_profile}}" \
      run_tool bracken -d "$KRAKEN2_DB" -i "$PR/kraken/${c}_reads.report" \
          -o "$PR/kraken/${c}_bracken_${lvl}.tsv" -r "$BRK_R" -l "$lvl" -t 10 \
          > "$PR/kraken/${c}_bracken_${lvl}.log" 2>&1 \
        || { echo "    bracken $c $lvl failed:"
             tail -3 "$PR/kraken/${c}_bracken_${lvl}.log" | sed 's/^/      /'; }
    done
  done
fi

# =====================================================================
# B. CONTIGS — what did the MAGs miss?
# =====================================================================
echo
echo "######## B. KRAKEN2 ON CONTIGS, SPLIT BINNED vs UNBINNED ($(date +%H:%M)) ########"
cat "$WS2/07_mags/bins"/*.fa 2>/dev/null | grep '^>' | sed 's/^>//;s/ .*//' \
  | sort -u > "$PR/contigs/binned_contigs.txt"
echo "  contigs in MAGs: $(wc -l < "$PR/contigs/binned_contigs.txt")"

for c in SED LIQ; do
  asm="$WS2/04_polish/${c}_asm_final.fasta"
  [ -s "$asm" ] || continue
  rep="$PR/contigs/${c}_contigs.report" out="$PR/contigs/${c}_contigs.kraken"
  if [ ! -s "$out" ]; then
    echo "  $c: classifying contigs ..."
    run_tool kraken2 --db "$KRAKEN2_DB" $MMAP --threads "$THREADS" \
        --confidence "$CONF" --report "$rep" --output "$out" "$asm" 2>&1 | sed 's/^/    /'
  else
    echo "  $c: contigs already classified"
  fi
done

python3 - "$PR" "$WS2" <<'PY'
import sys, os, csv, collections
PR, WS2 = sys.argv[1], sys.argv[2]

binned = set()
# ---- taxid -> scientific name, from the database's own taxonomy dump ----
tax2name, tax2rank = {}, {}
KDB = os.environ.get('KRAKEN2_DB', '')
nd = os.path.join(KDB, 'names.dmp')
if os.path.exists(nd):
    for line in open(nd):
        f = [x.strip() for x in line.split('|')]
        if len(f) >= 4 and f[3] == 'scientific name':
            tax2name[f[0]] = f[1]
no = os.path.join(KDB, 'nodes.dmp')
if os.path.exists(no):
    for line in open(no):
        f = [x.strip() for x in line.split('|')]
        if len(f) >= 3: tax2rank[f[0]] = f[2]
kt = os.path.join(KDB, 'ktaxonomy.tsv')
if not tax2name and os.path.exists(kt):
    for line in open(kt):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 5: tax2name[f[0]] = f[-1].strip()
print(f"  taxid -> name lookup: {len(tax2name)} entries "
      f"({'names.dmp' if os.path.exists(nd) else 'ktaxonomy.tsv' if os.path.exists(kt) else 'NONE — names will print as taxids'})")

p = os.path.join(PR, 'contigs', 'binned_contigs.txt')
if os.path.exists(p):
    binned = {l.strip() for l in open(p) if l.strip()}

for code in ('SED', 'LIQ'):
    kf = os.path.join(PR, 'contigs', f'{code}_contigs.kraken')
    df = os.path.join(WS2, '02_mapping', f'{code}_depth_multi.txt')
    if not os.path.exists(kf): continue

    depth = {}
    if os.path.exists(df):
        with open(df) as fh:
            next(fh, None)
            for line in fh:
                f = line.split('\t')
                if len(f) >= 3:
                    try: depth[f[0]] = float(f[1]) * float(f[2])   # length x depth = bases
                    except ValueError: pass

    tot = collections.Counter(); binned_w = collections.Counter()
    n_un = n_bn = 0; w_un = w_bn = 0.0
    for line in open(kf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 3: continue
        cid, tax = f[1], f[2]
        # kraken2 --output column 3 is a BARE TAXID, not 'name (taxid)' — that
        # format belongs to the --report file. The first version printed raw
        # numbers, which are unreadable and unciteable. Resolve via names.dmp.
        name = tax.split('(')[0].strip()
        if name.isdigit():
            name = tax2name.get(name, f'taxid:{name}')
        w = depth.get(cid, 0.0)
        tot[name] += w
        if cid in binned: binned_w[name] += w; n_bn += 1; w_bn += w
        else: n_un += 1; w_un += w

    total_w = w_bn + w_un
    print()
    print("=" * 92)
    print(f"  {code} — ASSEMBLY COMPOSITION BY COVERAGE-WEIGHTED TAXON")
    print("=" * 92)
    print(f"  contigs: {n_bn} binned / {n_un} unbinned")
    if total_w:
        print(f"  sequencing effort: {100*w_bn/total_w:.1f}% of mapped bases sit in MAGs, "
              f"{100*w_un/total_w:.1f}% unbinned")
    print()
    print(f"  {'taxon':<46}{'% of assembly':>14}{'% of that in MAGs':>19}")
    for name, w in tot.most_common(25):
        if not total_w: break
        inmag = 100*binned_w[name]/w if w else 0
        print(f"  {name[:46]:<46}{100*w/total_w:>13.2f}%{inmag:>18.1f}%")
    print()
    print("  Taxa with a high % of assembly and a LOW % in MAGs are the ones the")
    print("  binning missed — the honest answer to 'what is not in your catalogue'.")

    with open(os.path.join(PR, f'{code}_assembly_composition.tsv'), 'w', newline='') as fh:
        w_ = csv.writer(fh, delimiter='\t')
        w_.writerow(['taxon','weighted_bases','pct_of_assembly','pct_in_mags'])
        for name, w in tot.most_common():
            w_.writerow([name, f"{w:.0f}",
                         f"{100*w/total_w:.4f}" if total_w else '',
                         f"{100*binned_w[name]/w:.2f}" if w else ''])
PY

# =====================================================================
# C. 16S rRNA FROM THE WHOLE ASSEMBLY — comparable to the amplicon chapter
# =====================================================================
echo
echo "######## C. 16S rRNA GENES ACROSS THE FULL ASSEMBLY ($(date +%H:%M)) ########"
for c in SED LIQ; do
  asm="$WS2/04_polish/${c}_asm_final.fasta"
  [ -s "$asm" ] || continue
  gff="$PR/rrna/${c}_rrna.gff"
  if [ ! -s "$gff" ]; then
    echo "  $c: scanning for rRNA (bacteria + archaea) ..."
    { run_tool barrnap --kingdom bac --threads "$THREADS" --quiet "$asm"
      run_tool barrnap --kingdom arc --threads "$THREADS" --quiet "$asm"; } \
      > "$gff" 2>/dev/null
  fi
  n16=$(grep -c '16S_rRNA' "$gff" 2>/dev/null || true)
  printf '  %s: %s 16S rRNA genes in the assembly (MAGs hold only a fraction)\n' "$c" "${n16:-0}"
  # pull the 16S sequences out with their host contig, for later amplicon comparison
  awk -F'\t' '$0!~/^#/ && $9~/16S/ {print $1"\t"$4"\t"$5"\t"$7}' "$gff" \
    > "$PR/rrna/${c}_16S_coords.tsv" 2>/dev/null || true
  if [ -s "$PR/rrna/${c}_16S_coords.tsv" ]; then
    run_tool seqkit subseq --bed <(awk -F'\t' '{print $1"\t"($2-1)"\t"$3"\t16S\t0\t"$4}' \
        "$PR/rrna/${c}_16S_coords.tsv") "$asm" > "$PR/rrna/${c}_16S.fna" 2>/dev/null \
      || echo "    (seqkit subseq failed — coordinates are still in ${c}_16S_coords.tsv)"
    [ -s "$PR/rrna/${c}_16S.fna" ] && \
      echo "    -> $PR/rrna/${c}_16S.fna  ($(grep -c '^>' "$PR/rrna/${c}_16S.fna") sequences)"
  fi
done
echo
echo "  These 16S sequences are the bridge to the amplicon chapter: they come"
echo "  from the same DNA, in the same units, and are NOT limited to binned"
echo "  genomes. Arc_Otu22 can be searched against them directly."

# =====================================================================
# SUMMARY
# =====================================================================
echo
echo "######## CLASSIFIED FRACTIONS ########"
for c in SED LIQ; do
  rep="$PR/kraken/${c}_reads.report"
  [ -s "$rep" ] || continue
  unclass=$(awk -F'\t' '$4=="U"{print $1}' "$rep" | head -1)
  root=$(awk -F'\t' '$4=="R"{print $1}' "$rep" | head -1)
  printf '  %s reads: %s%% classified, %s%% unclassified\n' "$c" "${root:-?}" "${unclass:-?}"
  for rank in D P C O F G S; do
    n=$(awk -F'\t' -v r="$rank" '$4==r' "$rep" | wc -l)
    printf '     distinct taxa at rank %s: %s\n' "$rank" "$n"
  done
done

echo
echo "######## VALIDATION GATE ########"
for c in SED LIQ; do
  check "$c read report"    $([ -s "$PR/kraken/${c}_reads.report" ] && echo 0 || echo 1)
  check "$c contig report"  $([ -s "$PR/contigs/${c}_contigs.report" ] && echo 0 || echo 1)
  check "$c composition table" $([ -s "$PR/${c}_assembly_composition.tsv" ] && echo 0 || echo 1)
done

runlog "13_readprofile complete"
echo
echo "######## 13 DONE ########"
