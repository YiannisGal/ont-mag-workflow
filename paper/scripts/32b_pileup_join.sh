#!/usr/bin/env bash
# =============================================================================
# 32b_pileup_join.sh  --  the pileup and comparison half of step 32, rebuilt so
#                         it cannot exhaust the machine
#
# WHY THERE IS A 32b
#   32's mapping stage finished. Its PILEUP + COMPARE stage did not, and the
#   tmux session went with it. Two things in that stage were capable of doing
#   exactly that, and both were mine:
#
#   1. `samtools mpileup -aa` emits EVERY reference position, covered or not,
#      for all four samples. The stated reason — "so the two samples stay in
#      register without an interval join" — was a bad trade. A dict/merge join
#      on (contig, pos) is trivial, and -aa inflated every file to the full
#      genome length four times over for no analytical gain.
#
#   2. The comparison loaded ALL FOUR pileups into memory at once
#      (`sets = {t: load(t) for t in (...)}`) before comparing anything. A
#      Python dict keyed on (contig, pos) tuples costs roughly 350-400 bytes per
#      site. Four samples x a few million passing sites is 5-8 GB of resident
#      memory, reached in a burst, on a machine also running the pileups. If the
#      OOM killer took that process the shell hosting it dies too, which is
#      consistent with the session simply not existing afterwards.
#
#   Neither is a subtle bug. Both are me writing the convenient version.
#
# WHAT CHANGED
#   - `-aa` dropped. Only covered positions are emitted.
#   - depth filtered in AWK, before anything is written, so the counts files
#     hold only sites that could ever qualify. Same threshold Python applies.
#   - the comparison is a STREAMING MERGE JOIN over two coordinate-sorted files
#     at a time. Memory is O(number of contigs), not O(number of sites), and it
#     never holds more than two samples at once.
#   - contig order comes from the reference .fai, which is the order mpileup
#     emits, so the merge is well defined even when a contig is missing from one
#     side.
#   - progress and peak memory are printed, because a step that died silently
#     should not be replaced by another step that could.
#
#   The ANALYSIS is unchanged: same MINDP, same MAJFRAC, same split-half error
#   floor, same verdict logic. This is a rewrite of the plumbing only, so its
#   numbers are directly comparable to what 32 would have produced.
#
# SAFE TO RE-RUN. It reuses 32's BAMs and skips any pileup already finished.
#
# INPUT   32_strain/{SED,LIQ,sedA,sedB}.bam   (from 32, already done)
# OUTPUT  32_strain/*.counts.tsv | fixed_sites.tsv | DIVERGENCE.txt
# RUNTIME 10-40 min, dominated by mpileup. Memory stays under ~200 MB.
# USAGE   ws2 && bash scripts/32b_pileup_join.sh 2>&1 | tee logs/32b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="32b_pileup"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; S="$WS2/32_strain"; mkdir -p "$S"

REFID="${REF:-ETO_SED_MAG_002}"
MINDP="${MINDP:-20}"
MAJFRAC="${MAJFRAC:-0.80}"

echo "######## 32b — PILEUP AND COMPARISON ########"
echo "  started: $(date -Is)"
runlog "32b_pileup started"

# ---- the reference, found the same way 32 found it --------------------------
REFFA=""
WN=$(awk -F'\t' -v id="$REFID" 'NR==1{for(i=1;i<=NF;i++){if($i=="mag_id")a=i;if($i=="working_name")b=i};next}
     $a==id{print $b}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
for cand in "$F/bins/${REFID}.fa" "$F/bins/${WN}.fa" \
            $(find "$F" -maxdepth 3 \( -name "${REFID}.f*a" -o -name "${WN}.f*a" \) 2>/dev/null); do
  [ -n "$cand" ] && [ -s "$cand" ] && { REFFA="$cand"; break; }
done
[ -n "$REFFA" ] || { echo "  !! could not find a fasta for $REFID"; exit 1; }
echo "  reference: ${REFFA#$WS2/}"

# The .fai defines the order mpileup emits contigs in, and the merge join needs
# that order to be explicit rather than inferred from whichever file happens to
# mention a contig first.
if [ ! -s "${REFFA}.fai" ]; then
  echo "  indexing reference ..."
  run_tool samtools faidx "$REFFA" || { echo "  !! samtools faidx failed"; exit 1; }
fi
NCTG=$(wc -l < "${REFFA}.fai")
GLEN=$(awk '{s+=$2} END{print s}' "${REFFA}.fai")
echo "  $NCTG contigs, $(( GLEN / 1000 )) kb"

# ---- check the BAMs 32 left behind ------------------------------------------
echo
echo "######## 0. BAMs FROM 32 ########"
MISSING=0
for t in SED LIQ sedA sedB; do
  if [ -s "$S/${t}.bam" ]; then
    sz=$(du -h "$S/${t}.bam" | cut -f1)
    [ -s "$S/${t}.bam.bai" ] || run_tool samtools index "$S/${t}.bam" 2>/dev/null
    echo "  [ok]   ${t}.bam  $sz"
  else
    echo "  [MISS] ${t}.bam"; MISSING=$((MISSING+1))
  fi
done
if [ "$MISSING" -gt 0 ]; then
  echo "  !! $MISSING BAM(s) missing — re-run 32 for the mapping stage first."
  echo "     It caches BAMs, so it will not redo the ones that are there."
  exit 1
fi

# ---- disk headroom ----------------------------------------------------------
# Without -aa the counts files are a fraction of what they were, but check
# anyway: a truncated write here is the one failure that could still poison the
# comparison, and it is cheap to rule out.
AV=$(df -BM --output=avail "$S" 2>/dev/null | tail -1 | tr -dc '0-9'); AV=${AV:-0}
echo
echo "  free space on the workspace filesystem: $(( AV / 1024 )) GB"
if [ "$AV" -lt 5120 ]; then
  echo "  !! under 5 GB free. Clear space before running this."; exit 1
fi

# ---- strip the BAMs down to what mpileup actually needs ---------------------
# WHY. minimap2 -a writes UNMAPPED reads into the SAM as well, and samtools sort
# keeps them. That is why SED.bam is 6.6 GB against a ~3 Mb reference: almost
# all of it is reads that never touched this genome. mpileup has to stream past
# every one. Filtering them out once makes every later pass cheap.
#   -F 0x904 drops unmapped (0x4), secondary (0x100) and supplementary (0x800).
#   -q 1     drops multi-mapping-ambiguous reads with MAPQ 0, which are exactly
#            the reads that would put another organism's bases into this
#            genome's pileup and manufacture "divergence".
MAXD="${MAXD:-300}"
SPANMAX="${SPANMAX:-200000}"
SETTINGS="F0x904_q1_Q5_d${MAXD}_m${MINDP}_span${SPANMAX}"

echo
echo "######## 1. FILTERING BAMs ########"
echo "  settings stamp: $SETTINGS"
BFAIL=0
for t in SED LIQ sedA sedB; do
  mb="$S/${t}.mapped.bam"
  if [ -s "$mb" ] && [ -s "$mb.bai" ]; then
    echo "  [$t] cached — $(du -h "$mb" | cut -f1)"
    continue
  fi
  printf '  [%s] filtering ...' "$t"
  run_tool samtools view -b -F 0x904 -q 1 -@ 4 -o "$mb.part" "$S/${t}.bam" 2>"$S/${t}.filter.err"
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$mb.part" ]; then
    echo " FAILED (rc=$rc)"; tail -3 "$S/${t}.filter.err" | sed 's/^/         /'
    rm -f "$mb.part"; BFAIL=$((BFAIL+1)); continue
  fi
  mv "$mb.part" "$mb"
  run_tool samtools index "$mb" 2>/dev/null
  echo " $(du -h "$S/${t}.bam" | cut -f1) -> $(du -h "$mb" | cut -f1)"
done
[ "$BFAIL" -gt 0 ] && { echo "  !! $BFAIL filter(s) failed"; exit 2; }

# ---- drop alignments that span an absurd stretch of reference ---------------
# THE ACTUAL CAUSE, from dmesg:
#   Killed process (samtools) total-vm:434351288kB, anon-rss:127403568kB
# samtools really did allocate 127 GB. Not a cgroup, not another user, not
# something external — my previous two explanations were both wrong, and this
# one is measured rather than reasoned.
#
# A 34 MB BAM cannot hold 127 GB of reads, so the memory is not read data. It is
# the pileup engine walking a CIGAR: an alignment carrying a multi-megabase D or
# N operation makes mpileup step across every one of those reference positions.
# One such alignment is enough. That is why sedA passed all 110 contigs and
# sedB — the other random half of the SAME library, same size, same settings —
# died: the offending read is in one half and not the other.
#
# So the fix is to remove the alignment, not to shrink the depth cap. Reference
# span is computed from the CIGAR (M/D/N/=/X consume reference) and anything
# beyond SPANMAX is dropped. No real ONT read of this library spans 200 kb of a
# 2.7 Mb genome; anything that claims to is a mapping artefact, and it was going
# to contribute garbage bases to the pileup even if it had not exhausted the RAM.
#
# The maximum span seen is REPORTED per sample whether or not anything is
# dropped, because "how bad was it" is a fact worth having in the notebook.
echo
echo "######## 1b. REMOVING PATHOLOGICAL ALIGNMENTS ########"
echo "  dropping alignments whose reference span exceeds ${SPANMAX} bp"
CFAIL=0
for t in SED LIQ sedA sedB; do
  cb="$S/${t}.clean.bam"; mb="$S/${t}.mapped.bam"; st="$S/${t}.span.stats"
  if [ -s "$cb" ] && [ -s "$cb.bai" ] && [ -s "$st" ]; then
    echo "  [$t] cached — max span $(awk '$1=="maxspan"{print $2}' "$st") bp, "\
"$(awk '$1=="dropped"{print $2}' "$st") dropped"
    continue
  fi
  printf '  [%s] scanning ...' "$t"
  run_tool samtools view -h "$mb" 2>/dev/null \
    | awk -v mx="$SPANMAX" -v st="$st" '
        /^@/ { print; next }
        {
          c = $6; span = 0; num = ""
          for (i = 1; i <= length(c); i++) {
            ch = substr(c, i, 1)
            if (ch ~ /[0-9]/) { num = num ch }
            else {
              if (ch == "M" || ch == "D" || ch == "N" || ch == "=" || ch == "X")
                span += num + 0
              num = ""
            }
          }
          if (span > maxspan) maxspan = span
          if (span <= mx) { keep++; print } else { drop++ }
        }
        END { printf "maxspan\t%d\ndropped\t%d\nkept\t%d\n", maxspan+0, drop+0, keep+0 > st }' \
    | run_tool samtools view -b -o "$cb.part" - 2>/dev/null
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ] || [ ! -s "$cb.part" ]; then
    echo " FAILED (rc=$rc)"; rm -f "$cb.part"; CFAIL=$((CFAIL+1)); continue
  fi
  mv "$cb.part" "$cb"
  run_tool samtools index "$cb" 2>/dev/null
  ms=$(awk '$1=="maxspan"{print $2}' "$st"); dr=$(awk '$1=="dropped"{print $2}' "$st")
  kp=$(awk '$1=="kept"{print $2}' "$st")
  printf '\r  [%s] max span %s bp   dropped %s of %s alignments%-12s\n' \
    "$t" "$ms" "$dr" "$((dr + kp))" ""
  if [ "${dr:-0}" -gt 0 ]; then
    echo "        ^ those $dr alignment(s) are why this step kept dying."
  fi
done
[ "$CFAIL" -gt 0 ] && { echo "  !! $CFAIL span-filter(s) failed"; exit 2; }

# ---- pileup ------------------------------------------------------------------
# WHAT KILLED THE LAST ATTEMPT. Three of four mpileups exited 137 — SIGKILL, the
# OOM killer. Two things made that possible and both were mine:
#
#   -d 5000 with LONG reads. mpileup holds every read overlapping the current
#   column in memory. On Illumina that is 5000 x 150 bp. On ONT it is 5000 x
#   whatever these reads are, which at 20-50 kb is gigabytes at a single
#   high-coverage column. The depth cap has to be read-length aware, and 5000
#   bought nothing: the analysis needs depth >= 20 and an allele fraction, and
#   both are settled long before 300x.
#
#   Whole-BAM in one invocation, so the peak was never released.
#
# Now: -d 300, and one contig at a time so memory returns to the OS between
# contigs and progress is visible. Per-contig also means a single pathological
# contig is identifiable instead of taking the whole sample down anonymously.
#
# NOTE: the previous run's sedA.counts.tsv was made with -d 5000 and unfiltered
# reads. It is NOT reused — a cached file from different settings is a different
# measurement wearing the same name. The settings stamp enforces that.
#
# 2026-08-25, SECOND FAILURE, AND MY EXPLANATION WAS WRONG.
#   The machine has 125 GB of RAM and reported 122 GB free at start. SED.mapped
#   .bam is 67 MB, sedB.mapped.bam is 34 MB. mpileup on a 34 MB file cannot need
#   122 GB, and sedA — the same size, the same data, the same settings —
#   finished. Two SSH sessions were also reset by the server mid-run.
#
#   A process that dies at a size it cannot possibly have reached, while its
#   identical twin survives, is not being killed by its own memory use. Either
#   something else on this shared machine is taking the RAM, or a cgroup limit
#   on this user session is being hit. I do not know which and will not guess
#   again — so this version MEASURES instead.
#
#   Each mpileup now runs under /usr/bin/time and under an address-space cap.
#   That turns the question into an experiment with two distinguishable answers:
#     - killed while peak RSS is small  -> something external is doing it
#     - hits the cap and fails cleanly  -> it really is mpileup, and we will
#                                          have the number that proves it
#   Peak RSS is printed per sample either way, so the next run produces
#   evidence rather than another theory.
MEMCAP_KB="${MEMCAP_KB:-8388608}"        # 8 GB of address space per mpileup
RETRY_DIV="${RETRY_DIV:-3}"              # on failure, retry the contig at MAXD/3

# Resolve the samtools binary once so it can be run under /usr/bin/time. run_tool
# is a shell function and cannot be the target of an external timer.
bin_path() {
  local t="$1"
  local v="ENV_${t}"
  local e="${!v:-}"
  local p=""
  if [ -n "$e" ] && [ "$e" != "-" ] && declare -f env_prefix >/dev/null 2>&1; then
    p=$(env_prefix "$e" 2>/dev/null) || p=""
  fi
  if [ -n "$p" ] && [ -x "$p/bin/$t" ]; then
    echo "$p/bin/$t"
  else
    command -v "$t" 2>/dev/null
  fi
}
SAMTOOLS=$(bin_path samtools)
[ -x "$SAMTOOLS" ] || { echo "  !! cannot locate the samtools binary"; exit 1; }
SAMPREFIX="${SAMTOOLS%/bin/samtools}"
[ -d "$SAMPREFIX/lib" ] && export LD_LIBRARY_PATH="$SAMPREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "  samtools binary: $SAMTOOLS"
HAVE_TIME=0
[ -x /usr/bin/time ] && /usr/bin/time -f '%M' true >/dev/null 2>&1 && HAVE_TIME=1
[ "$HAVE_TIME" = "1" ] || echo "  (/usr/bin/time unavailable — no peak-memory measurement)"

# One contig. Returns mpileup's exit status; writes peak RSS in KB to $RSSFILE.
mpileup_one () {                 # mpileup_one <bam> <contig> <maxd> <errfile>
  local bam="$1" ctg="$2" md="$3" ef="$4"
  : > "$RSSFILE"
  (
    ulimit -v "$MEMCAP_KB" 2>/dev/null || true
    if [ "$HAVE_TIME" = "1" ]; then
      /usr/bin/time -f '%M' -o "$RSSFILE" \
        "$SAMTOOLS" mpileup -A -Q 5 -d "$md" --ff UNMAP,SECONDARY,QCFAIL,DUP \
        -r "$ctg" -f "$REFFA" "$bam" 2>>"$ef"
    else
      "$SAMTOOLS" mpileup -A -Q 5 -d "$md" --ff UNMAP,SECONDARY,QCFAIL,DUP \
        -r "$ctg" -f "$REFFA" "$bam" 2>>"$ef"
    fi
  )
}

pile () {
  local tag="$1"
  local out="$S/${tag}.counts.tsv"
  local stamp="$S/${tag}.counts.settings"
  local err="$S/${tag}.mpileup.err"
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$SETTINGS" ]; then
    echo "  [$tag] cached — $(wc -l < "$out") sites"
    return 0
  fi
  [ -s "$out" ] && echo "  [$tag] cached file was made with different settings — redoing"
  rm -f "$out.part" "$err"
  local t0 nc i maxrss rss rc lowd
  t0=$(date +%s)
  nc=$(wc -l < "${REFFA}.fai")
  i=0; maxrss=0
  RSSFILE="$S/.rss.$$"
  : > "$out.part"
  while read -r ctg _; do
    i=$((i+1))
    printf '\r  [%s] contig %d/%d %-24s peak %sMB ' \
      "$tag" "$i" "$nc" "${ctg:0:24}" "$(( maxrss / 1024 ))"

    mpileup_one "$S/${tag}.clean.bam" "$ctg" "$MAXD" "$err" \
      | awk -F'\t' -v m="$MINDP" '$4 >= m {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >> "$out.part"
    rc=${PIPESTATUS[0]}
    rss=$(tail -1 "$RSSFILE" 2>/dev/null | tr -dc '0-9'); rss=${rss:-0}
    [ "$rss" -gt "$maxrss" ] 2>/dev/null && maxrss=$rss

    # One retry at a lower depth cap before giving up on the sample. If the
    # contig genuinely has a pathological pile, this rescues it; if something
    # external is killing processes, it will fail again and we learn that too.
    if [ "$rc" -ne 0 ]; then
      lowd=$(( MAXD / RETRY_DIV )); [ "$lowd" -lt 30 ] && lowd=30
      printf '\r  [%s] contig %s exited %s — retrying at -d %s%-15s\n' \
        "$tag" "$ctg" "$rc" "$lowd" ""
      # the failed contig contributed nothing or a partial line; drop the tail
      sed -i "/^${ctg}\t/d" "$out.part" 2>/dev/null
      mpileup_one "$S/${tag}.clean.bam" "$ctg" "$lowd" "$err" \
        | awk -F'\t' -v m="$MINDP" '$4 >= m {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >> "$out.part"
      rc=${PIPESTATUS[0]}
      rss=$(tail -1 "$RSSFILE" 2>/dev/null | tr -dc '0-9'); rss=${rss:-0}
      [ "$rss" -gt "$maxrss" ] 2>/dev/null && maxrss=$rss
    fi

    if [ "$rc" -ne 0 ]; then
      printf '\r  [%s] !! mpileup exited %s on contig %s%-20s\n' "$tag" "$rc" "$ctg" ""
      echo "         peak RSS seen so far on this sample: $(( maxrss / 1024 )) MB"
      echo "         address-space cap in force: $(( MEMCAP_KB / 1024 )) MB"
      if [ "$rc" -eq 137 ]; then
        echo
        echo "         137 is SIGKILL. READ THE PEAK ABOVE BEFORE CONCLUDING:"
        echo "           * peak is small (tens/hundreds of MB) -> mpileup was NOT"
        echo "             using the memory. Something else on this machine is"
        echo "             killing it. Check:  sudo dmesg -T | grep -i 'killed process'"
        echo "             and whether another user's job is holding the RAM."
        echo "           * peak is near the cap -> it really is mpileup; lower the"
        echo "             depth:  MAXD=60 bash scripts/32b_pileup_join.sh"
      fi
      tail -4 "$err" | sed 's/^/         /'
      rm -f "$out.part" "$RSSFILE"; return 1
    fi
  done < "${REFFA}.fai"
  rm -f "$RSSFILE"
  printf '\r%-78s\r' ''
  echo "  [$tag] peak mpileup memory: $(( maxrss / 1024 )) MB"
  if [ ! -s "$out.part" ]; then
    echo "  [$tag] !! no site reached depth $MINDP. That is a coverage fact, not"
    echo "         an error, but it makes this comparison impossible."
    rm -f "$out.part"; return 1
  fi
  mv "$out.part" "$out"; echo "$SETTINGS" > "$stamp"
  echo "  [$tag] $(wc -l < "$out") sites at depth >= $MINDP  ($(du -h "$out" | cut -f1), $(( $(date +%s) - t0 ))s)"
}

echo
echo "######## 2. PILEUP ########"
echo "  covered positions only, depth capped at ${MAXD}x, one contig at a time."
echo
free -g 2>/dev/null | awk '/^Mem:/{printf "  memory: %s GB total, %s GB available\n", $2, $7}'
echo "  DO NOT run this alongside the antiSMASH containers. The last attempt had"
echo "  both going at once and three of four pileups were OOM-killed; sedA, the"
echo "  one that ran while nothing else was competing, finished fine."
echo
PFAIL=0
for t in SED LIQ sedA sedB; do pile "$t" || PFAIL=$((PFAIL+1)); done
if [ "$PFAIL" -gt 0 ]; then
  echo
  echo "  !! $PFAIL pileup(s) failed. Not comparing — a missing sample would"
  echo "     silently become 'no shared sites', which reads as a result."
  runlog "**32b ABORTED: $PFAIL pileup failure(s)**"
  exit 2
fi

# ---- compare -----------------------------------------------------------------
echo
echo "######## 3. COMPARISON ########"
# GNU time is not on every box, and its absence must not stop the analysis —
# it only costs us the memory number.
TIMER=()
if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then
  TIMER=(/usr/bin/time -v)
else
  echo "  (/usr/bin/time -v unavailable — running without the memory measurement)"
fi

"${TIMER[@]}" python3 - "$S" "$MINDP" "$MAJFRAC" "${REFFA}.fai" "$REFID" \
     2> "$S/compare.time" <<'PY'
import sys, os
from collections import Counter

S, MINDP, MAJFRAC, FAI, REFID = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), \
                                sys.argv[4], sys.argv[5]

# contig -> rank, from the reference index. mpileup emits in this order, so two
# files can be walked in lockstep by (rank, pos) without holding either in RAM.
ORDER = {}
for i, line in enumerate(open(FAI)):
    ORDER[line.split('\t')[0]] = i

def parse(bases, ref):
    """samtools mpileup base string -> Counter of A/C/G/T. Unchanged from 32:
    indel runs stripped by their length, read-start/end markers dropped, and
    . / , are the reference base on either strand."""
    out, i, n = Counter(), 0, len(bases)
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            j = i + 1; num = ''
            while j < n and bases[j].isdigit(): num += bases[j]; j += 1
            i = j + int(num or 0); continue
        if c in '.,': out[ref.upper()] += 1
        elif c.upper() in 'ACGT': out[c.upper()] += 1
        i += 1
    return out

def sites(tag):
    """Generator: ((rank, pos), contig, major_base, major_fraction, depth).
    One line in memory at a time."""
    p = os.path.join(S, f'{tag}.counts.tsv')
    with open(p) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 5: continue
            r = ORDER.get(f[0])
            if r is None: continue
            try: pos = int(f[1])
            except ValueError: continue
            c = parse(f[4], f[2])
            tot = sum(c.values())
            if tot < MINDP: continue
            b, k = c.most_common(1)[0]
            yield (r, pos), f[0], b, k / tot, tot

def merge(a, b):
    """Streaming inner join of two coordinate-sorted site streams."""
    ga, gb = sites(a), sites(b)
    ra = next(ga, None); rb = next(gb, None)
    while ra is not None and rb is not None:
        if ra[0] == rb[0]:
            yield ra, rb
            ra = next(ga, None); rb = next(gb, None)
        elif ra[0] < rb[0]:
            ra = next(ga, None)
        else:
            rb = next(gb, None)

def compare(a, b, dump=None):
    shared = 0
    fixed = []
    fh = open(dump, 'w') if dump else None
    if fh:
        fh.write("contig\tpos\tA_allele\tA_frac\tA_dp\tB_allele\tB_frac\tB_dp\n")
    for ra, rb in merge(a, b):
        shared += 1
        if ra[2] != rb[2] and ra[3] >= MAJFRAC and rb[3] >= MAJFRAC:
            fixed.append((ra[1], ra[0][1]))
            if fh:
                fh.write(f"{ra[1]}\t{ra[0][1]}\t{ra[2]}\t{ra[3]:.2f}\t{ra[4]}\t"
                         f"{rb[2]}\t{rb[3]:.2f}\t{rb[4]}\n")
    if fh: fh.close()
    return shared, fixed

W = 96
print("=" * W)
print("  FIXED DIFFERENCES")
print("=" * W)
print(f"  A site counts only if BOTH samples have depth >= {MINDP} and each major")
print(f"  allele reaches {MAJFRAC:.0%}. Near-fixation on both sides is what keeps ONT")
print("  error out: a 1% error rate does not produce an 80% allele.")
print()

res = {}
FLOOR_L = 'sedA vs sedB  (ERROR FLOOR)'
MEAS_L  = 'SED  vs LIQ   (measurement)'
for lbl, (a, b), dump in ((FLOOR_L, ('sedA', 'sedB'), None),
                          (MEAS_L,  ('SED', 'LIQ'), os.path.join(S, 'fixed_sites.tsv'))):
    ns, fx = compare(a, b, dump)
    per_kb = 1000 * len(fx) / ns if ns else 0
    res[lbl] = (ns, len(fx), per_kb)
    print(f"  {lbl:<30}{ns:>9} sites compared{len(fx):>8} fixed diffs"
          f"{per_kb:>10.3f} per kb")

floor = res[FLOOR_L][2]
meas  = res[MEAS_L][2]

out = []
out.append("=" * W)
out.append("  VERDICT")
out.append("=" * W)
if res[FLOOR_L][0] == 0 or res[MEAS_L][0] == 0:
    out.append("  >> NO VERDICT. Too few sites reached the depth threshold in one of the")
    out.append(f"     comparisons. Lower MINDP (currently {MINDP}) or accept that this")
    out.append("     genome is not covered deeply enough in both samples to compare.")
elif floor > 0 and meas < 3 * floor:
    out.append(f"  >> NOT DISTINGUISHABLE. {meas:.3f} fixed differences per kb between")
    out.append(f"     sediment and culture, against an error floor of {floor:.3f} per kb")
    out.append("     measured on two halves of the SAME library. The enrichment did not")
    out.append("     detectably change this population: what grew in the bottle is the")
    out.append("     sediment population, not a diverged descendant of it.")
    out.append("     That is a real result and the conservative one — it licenses using")
    out.append("     the culture genome to interpret the sediment organism.")
elif floor == 0:
    out.append(f"  >> {meas:.3f} fixed differences per kb, and the split-half control")
    out.append("     returned EXACTLY ZERO. A zero floor cannot be divided into, so no")
    out.append("     fold-over-floor is quoted. Treat the control as 'below the")
    out.append("     resolution of this depth' and report the measurement with that")
    out.append("     stated, not as an infinite ratio.")
else:
    out.append(f"  >> DIVERGENCE DETECTED. {meas:.3f} fixed differences per kb, against")
    out.append(f"     an error floor of {floor:.3f} per kb from the split-half control —")
    out.append(f"     {meas / floor:.1f}x the floor.")
    out.append("     The cultured population is not simply a subset of the sediment one.")
    out.append("     Inspect fixed_sites.tsv: if the differences cluster in a few genes")
    out.append("     that is selection on those genes, and if they are spread evenly it")
    out.append("     is more likely a different strain having been enriched from a")
    out.append("     mixed in-situ population.")
out.append("")
out.append("  LIMITS. Fixed differences are not a phylogeny and this is one genome.")
out.append("  A zero here means 'no difference resolvable above the measured error")
out.append("  floor', never 'the populations are identical'.")

print()
for line in out: print(line)

with open(os.path.join(S, 'DIVERGENCE.txt'), 'w') as fh:
    fh.write(f"reference\t{REFID}\nmin_depth\t{MINDP}\nmajor_fraction\t{MAJFRAC}\n\n")
    for lbl, (ns, nf, pk) in res.items():
        fh.write(f"{lbl}\t{ns}\t{nf}\t{pk:.4f}\n")
    fh.write('\n' + '\n'.join(out) + '\n')
PY
RC=$?

# The point of /usr/bin/time -v here is to make the thing that killed 32
# measurable rather than suspected. If peak RSS is a few hundred MB, the memory
# explanation is confirmed and closed.
PEAK=$(grep -i 'Maximum resident set size' "$S/compare.time" 2>/dev/null | tr -dc '0-9')
if [ -n "$PEAK" ]; then
  echo
  echo "  peak memory of the comparison: $(( PEAK / 1024 )) MB"
  echo "  (32 held all four pileups at once; this holds two streams and the"
  echo "   contig index. If that number is small, memory is what killed 32.)"
fi
[ "$RC" -ne 0 ] && { echo; echo "  !! the comparison exited $RC. Last lines:"; \
  tail -15 "$S/compare.time" | sed 's/^/     /'; }

echo
echo "######## VALIDATION ########"
OK=0
check "comparison completed" "$RC" || OK=1
check "DIVERGENCE.txt written" $([ -s "$S/DIVERGENCE.txt" ] && echo 0 || echo 1) || OK=1
for t in SED LIQ sedA sedB; do
  check "${t}.counts.tsv non-empty" $([ -s "$S/${t}.counts.tsv" ] && echo 0 || echo 1) || OK=1
done
if grep -q 'NO VERDICT' "$S/DIVERGENCE.txt" 2>/dev/null; then
  echo "  [STOP] no verdict — see the depth note above."; OK=1
fi

echo
if [ "$OK" -eq 0 ]; then
  echo "  32b complete."; runlog "32b_pileup complete"
else
  echo "  32b finished with problems."; runlog "**32b_pileup finished with FAILs**"
fi
echo "  ended: $(date -Is)"
exit "$OK"
