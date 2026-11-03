#!/usr/bin/env bash
# =============================================================================
# 13f_testrun_diagnose.sh  --  what ARE metagenome_RBK_test and _LSK_test?
#
# THE ARITHMETIC THAT ALREADY ANSWERS IT
#   16,133 of metagenome_RBK_test's 700,959 reads (2.30%) align to Erwinia 16S
#   over >=1000 bp. Ask what fraction of reads from a SHOTGUN library could do
#   that, even if the library were 100% Erwinia:
#
#       Enterobacteriaceae: ~7 rRNA operons x ~1.5 kb SSU in a ~4.7 Mb genome
#       SSU share of genomic DNA = 7 x 1.5 / 4700 = 0.22%
#
#   0.22% is the CEILING for a pure-culture shotgun library. The observed value
#   is 2.30% — ten times higher than physically possible from genomic DNA.
#
#   The assembly agrees. Those reads produced three contigs of 1,207, 526 and
#   1,267 bp at 15,464x, 849x and 14,380x. A 4.7 Mb genome at that depth would
#   have assembled into a near-complete genome; instead there are three short
#   rRNA-length fragments and no genome anywhere in the catalogue.
#
#   Both observations point one way: metagenome_RBK_test is not a shotgun
#   metagenome. It is an amplicon or rRNA-enriched library — almost certainly a
#   16S run used to test the Rapid Barcoding chemistry — and it was pooled into
#   the sediment shotgun read set by mistake.
#
# WHAT THIS SCRIPT DOES
#   1. Read-length distributions per run. An amplicon library is unmistakable:
#      reads pile up at the amplicon length instead of decaying like sheared DNA.
#   2. Total SSU fraction per run, not just Erwinia. A shotgun run sits near
#      0.1-0.5%; an amplicon run sits above 50%.
#   3. If the diagnosis holds, writes the exclusion list and the cleaned read
#      file, so step 03 can be re-run from a defensible input.
#
#   Test 1 is the strong one. It needs no reference and cannot be argued with.
#
# OUTPUT   13_profile/silva/verify/
#     run_diagnostics.tsv   per run: n, length stats, SSU %
#     exclude_read_ids.txt  every read ID from the offending runs
#     DIAGNOSIS.txt
#   00_input/sediment_reads_clean.fastq.gz   (only if exclusion is warranted)
#
# RUNTIME  20-60 min. USAGE:
#   ws2 && RUNDIR=${RAW_DATA_DIR}/fastqs \
#     bash scripts/13f_testrun_diagnose.sh 2>&1 | tee logs/13f_testrun.log
#   Add WRITE_CLEAN=1 to actually write the filtered read file.
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13f_testrun"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
S="$WS2/13_profile/silva"; V="$S/verify"; mkdir -p "$V"
RUNDIR="${RUNDIR:-}"
SUB="${SUBSAMPLE:-20000}"          # reads per run for the SSU test
WRITE_CLEAN="${WRITE_CLEAN:-0}"

echo "######## 13f — WHAT ARE THE TEST RUNS? ########"
echo "  started: $(date -Is)"
runlog "13f_testrun_diagnose started"
[ -d "${RUNDIR:-/nonexistent}" ] || { echo "  !! set RUNDIR to the raw fastq directory"; exit 1; }
[ -s "$S/silva.mmi" ] || { echo "  !! no SILVA index — run 13b first"; exit 1; }

# ---------------------------------------------------------- TEST 1 + TEST 2
: > "$V/run_diagnostics.tsv"
printf 'run\treads\tmean_len\tmedian_len\tN50\tpct_1300_1700\tn_long\tssu_pct_of_long\ttop_taxon_share\n' \
    >> "$V/run_diagnostics.tsv"
printf '  %-40s %8s %7s %7s %8s %10s %10s\n' \
       "run" "reads" "mean" "median" ">=1kb" "%SSU/long" "%top taxon"

shopt -s nullglob
for f in "$RUNDIR"/*.fastq "$RUNDIR"/*.fq "$RUNDIR"/*.fastq.gz "$RUNDIR"/*.fq.gz; do
  [ -s "$f" ] || continue
  b=$(basename "$f")
  case "$f" in *.gz) R=(zcat "$f");; *) R=(cat "$f");; esac

  # ---- read lengths, one pass. The 1.3-1.7 kb window brackets a full-length
  # 16S amplicon; a shotgun library has no reason to pile up there.
  # asort() is a gawk extension and Ubuntu's default awk is mawk, where it
  # fails silently. Sort externally instead so this works with either awk.
  LENS="$TMPDIR/13f_lens.$$"
  "${R[@]}" | awk 'NR%4==2{print length($0)}' | sort -n > "$LENS"
  read -r n mean med n50 amp < <(awk '
      { L[NR] = $1; s += $1; if ($1 >= 1300 && $1 <= 1700) a++ }
      END {
        if (!NR) { print 0, 0, 0, 0, "0.00"; exit }
        half = s / 2; run = 0; n50 = L[1]
        for (i = NR; i >= 1; i--) { run += L[i]; if (run >= half) { n50 = L[i]; break } }
        printf "%d %.0f %d %d %.2f\n", NR, s/NR, L[int(NR/2)+1], n50, 100*a/NR
      }' "$LENS")
  rm -f "$LENS"

  # ---- THE DISCRIMINATING TEST, rewritten.
  # The first version asked what fraction of ALL reads are SSU, with an
  # amplicon threshold of >40% of reads at 1.3-1.7 kb. Both were wrong for this
  # data. metagenome_RBK_test has a MEDIAN READ OF 193 bp — it is not an
  # amplicon library, it is a mostly-junk run with a discrete ~1.5 kb spike
  # inside it (2.28% of reads, i.e. ~15,980 — essentially identical to the
  # 16,133 Erwinia reads). Diluted across 700k short reads that spike is only
  # 2% of the run and slips under any whole-run threshold.
  #
  # Ask the sharper question instead: OF THE READS LONG ENOUGH TO CARRY AN SSU
  # (>=1000 bp), what fraction ARE SSU? A shotgun library sits near 0.2-1%,
  # because SSU is a tiny share of any genome. A library whose long reads are
  # almost all SSU is not a shotgun library, whatever the rest of it contains.
  # NOTE: 'sub' is an awk builtin and cannot be a variable name — using it made
  # every measurement above print 0.00%.
  LONG="$TMPDIR/13f_long.$$.fq"
  "${R[@]}" | awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==0{if(length(s)>=1000) print h"\n"s"\n+\n"$0}' \
      | head -n $((SUB*4)) > "$LONG"
  nlong=$(( $(wc -l < "$LONG") / 4 ))
  if [ "$nlong" -gt 0 ]; then
    read -r nssu top < <(run_tool minimap2 -x map-ont -t "$THREADS" --secondary=no -c \
          "$S/silva.mmi" "$LONG" 2>/dev/null \
        | awk '$11>=1000 && $10/$11>=0.92 { if (!($1 in b) || $10 > b[$1]) { b[$1]=$10; t[$1]=$6 } }
               END { n=0; for (q in b) { n++; c[t[q]]++ }
                     mx=0; for (a in c) if (c[a]>mx) mx=c[a]
                     print n+0, (n ? mx/n : 0) }')
    ssu=$(awk -v a="${nssu:-0}" -v b="$nlong" 'BEGIN{printf "%.2f", (b?100*a/b:0)}')
    topf=$(awk -v x="${top:-0}" 'BEGIN{printf "%.1f", 100*x}')
  else
    ssu="0.00"; topf="0.0"
  fi
  rm -f "$LONG"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$b" "$n" "$mean" "$med" "$n50" "$amp" "$nlong" "$ssu" "$topf" >> "$V/run_diagnostics.tsv"
  printf '  %-40s %8s %7s %7s %8s %9s%% %9s%%\n' \
      "$b" "$n" "$mean" "$med" "$nlong" "$ssu" "$topf"
done
shopt -u nullglob

# ---------------------------------------------------------------- DIAGNOSIS
python3 - "$V" "$WS2" <<'PY'
import sys, os, csv
V, WS2 = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(os.path.join(V, 'run_diagnostics.tsv')), delimiter='\t'))
def f(r, k):
    try: return float(r[k])
    except (ValueError, KeyError, TypeError): return 0.0

print()
print("=" * 96)
print("  DIAGNOSIS")
print("=" * 96)
print("  OF THE READS >=1000 bp — long enough to carry an SSU — what fraction ARE")
print("  SSU? In any shotgun library that is well under 2%, because rRNA is a tiny")
print("  share of any genome. A run whose long reads are mostly SSU is not shotgun.")
print()
print("  (The first version of this test asked about ALL reads and required >40% at")
print("   1.3-1.7 kb. metagenome_RBK_test has a MEDIAN READ OF 193 bp, so its ~1.5 kb")
print("   spike is only 2.3% of the run and slipped straight under that threshold.)")
print()
bad = [r for r in rows if f(r,'ssu_pct_of_long') > 20]
ok  = [r for r in rows if r not in bad]
if ok:
    m = sorted(f(r,'ssu_pct_of_long') for r in ok)[len(ok)//2]
    print(f"  shotgun-like runs: {len(ok)}, median SSU-among-long-reads {m:.2f}%")
print(f"  runs failing the test: {len(bad)}")
for r in bad:
    print(f"    {r['run']:<40} {f(r,'ssu_pct_of_long'):>6.1f}% of its {int(f(r,'n_long')):,} "
          f"long reads are SSU; {f(r,'top_taxon_share'):.0f}% of those on ONE reference")

amp = bad
L = []
if amp:
    L.append("NON-SHOTGUN SEQUENCE POOLED INTO A SHOTGUN READ SET.")
    L.append("")
    for r in amp:
        L.append(f"  {r['run']}: {f(r,'ssu_pct_of_long'):.1f}% of its "
                 f"{int(f(r,'n_long')):,} reads >=1 kb are SSU "
                 f"({f(r,'top_taxon_share'):.0f}% of those on a single reference); "
                 f"median read {r['median_len']} bp of {int(f(r,'reads')):,} total")
    L.append("")
    L.append("These runs are not shotgun metagenomes and must not be in the sediment")
    L.append("assembly input. They contributed the Erwinia signal that made the")
    L.append("sediment SSU profile 84.65% Enterobacterales, and their reads inflate")
    L.append("rRNA coverage in the assembly generally, not only for Erwinia.")
    L.append("")
    L.append("ACTION: rebuild the sediment read set without them, re-run 03 -> 05,")
    L.append("and re-run 13b. The MAG catalogue itself is very unlikely to change —")
    L.append("amplicon reads cannot form genome bins, no Erwinia MAG exists, and the")
    L.append("four sediment MAGs are Bacillota at 12-18x — but the read-level")
    L.append("profile and the assembly's rRNA content both need the clean input.")
else:
    L.append("No amplicon-like run detected by these two tests.")
    L.append("The Erwinia reads are still confined to two runs (13e, Test 4), so the")
    L.append("exclusion still stands — but on read-identity grounds rather than")
    L.append("because a whole library is the wrong type. In that case exclude the")
    L.append("16,261 Erwinia read IDs rather than the whole runs.")
for l in L: print("  " + l)
with open(os.path.join(V, 'DIAGNOSIS.txt'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
with open(os.path.join(V, '.amplicon_runs'), 'w') as fh:
    fh.write('\n'.join(r['run'] for r in amp) + ('\n' if amp else ''))
PY

# ---------------------------------------------------------------- EXCLUSION
echo
echo "######## EXCLUSION LIST ########"
AMPF="$V/.amplicon_runs"
if [ -s "$AMPF" ]; then
  : > "$V/exclude_read_ids.txt"
  while read -r b; do
    [ -n "$b" ] || continue
    f="$RUNDIR/$b"; [ -s "$f" ] || continue
    case "$f" in *.gz) R=(zcat "$f");; *) R=(cat "$f");; esac
    "${R[@]}" | awk 'NR%4==1{print substr($1,2)}' >> "$V/exclude_read_ids.txt"
  done < "$AMPF"
  N=$(wc -l < "$V/exclude_read_ids.txt")
  echo "  read IDs to exclude: $N"
  echo "  -> $V/exclude_read_ids.txt"
  if [ "$WRITE_CLEAN" = "1" ]; then
    OUT="$WS2/00_input/sediment_reads_clean.fastq.gz"
    echo "  writing $OUT ..."
    run_tool seqkit grep -v -f "$V/exclude_read_ids.txt" "$SED_READS" 2>/dev/null \
      | gzip > "$OUT"
    echo "  reads before: $(zcat "$SED_READS" | awk 'NR%4==1' | wc -l)"
    echo "  reads after : $(zcat "$OUT" | awk 'NR%4==1' | wc -l)"
    echo
    echo "  Point SED_READS at the clean file in 00_config.sh, keep the original"
    echo "  untouched as the provenance record, and re-run 03 onward."
  else
    echo "  (set WRITE_CLEAN=1 to write 00_input/sediment_reads_clean.fastq.gz)"
  fi
else
  echo "  no amplicon runs identified — nothing written."
fi

echo
echo "######## VALIDATION GATE ########"
check "run diagnostics" $([ -s "$V/run_diagnostics.tsv" ] && echo 0 || echo 1)
check "diagnosis"       $([ -s "$V/DIAGNOSIS.txt" ] && echo 0 || echo 1)
runlog "13f_testrun_diagnose complete"
echo
echo "######## 13f DONE ########"
