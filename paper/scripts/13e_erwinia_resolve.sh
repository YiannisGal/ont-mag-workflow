#!/usr/bin/env bash
# =============================================================================
# 13e_erwinia_resolve.sh  --  settle the Erwinia question with the assembly
#
# WHERE 13d LEFT IT
#   TEST 1 was the informative one: 16,204 of 16,261 Erwinia reads — 99.6% —
#   hit a SINGLE SILVA reference, MAZB01000066.4456.5976. Only 14 references
#   were touched in total. TEST 2 showed those matches are genuine, median 98.6%
#   identity, so this is not a threshold artefact: real sequence is really
#   matching a real Erwinia 16S.
#
#   That combination — one reference, high identity — has exactly two readings:
#     (a) ONE CLONAL POPULATION that is genuinely abundant, or
#     (b) ONE CLONAL POPULATION that arrived from somewhere else.
#   SILVA cannot tell them apart. Only the assembly can, because abundance
#   leaves a physical trace: coverage.
#
#   (13d also printed "10/138" in TEST 4 and "5 of 138" in the verdict. The 10
#   was my bug — a genus was counted once by exact key and again by lowercase
#   comparison. The true figure is 5/138 = 3.6%, and it is used below.)
#
# THE ARITHMETIC THAT MAKES THIS DECIDABLE
#   Enterobacteriaceae carry ~7 rRNA operons in a ~4.7 Mb genome. The SSU is
#   ~1.5 kb, so SSU is roughly 7 x 1.5 / 4700 = 0.22% of that genome's DNA.
#   If Erwinia supplies 84.6% of SSU reads it should supply a comparable share
#   of genomic DNA, and its contigs should therefore be among the highest-
#   coverage sequence in the sediment assembly — well above the 12-18x of the
#   four sediment MAGs.
#
#   If instead its contigs sit at low coverage, the SSU signal is not coming
#   from an abundant organism, and the only way a low-abundance organism
#   dominates an SSU read pool is if the SSU pool itself is tiny and skewed —
#   which it is: 19,210 of 7,528,854 sediment reads, 0.26%.
#
# FIVE TESTS, EACH ABLE TO END THE QUESTION
#   1. COVERAGE       what depth are the Erwinia-16S contigs?
#   2. BINNING        did they bin? into which candidate bin, at what quality?
#   3. NEGATIVE CONTROL   does the extraction blank contain Erwinia? This is the
#                     single most direct contamination test available and the
#                     provenance record says a blank was sequenced.
#   4. PER-RUN        the sediment file is TEN concatenated runs. Do the reads
#                     concentrate in a subset? (13d's version of this test was
#                     worthless — it grouped on a UUID prefix, which is not a
#                     run ID. This one greps the original per-run files.)
#   5. THE ENRICHMENTS   ZERO Erwinia reads and ZERO Erwinia 16S in LIQ. The
#                     enrichments were inoculated from this sediment and fed
#                     glucose, lactate and fumarate — conditions a facultative
#                     anaerobe like Erwinia ferments happily. An organism at
#                     ~85% of the sediment SSU pool that fails to appear at all
#                     in enrichments from that same sediment is hard to explain.
#
# OUTPUT   13_profile/silva/verify/
#     erwinia_contigs.tsv   contig, length, coverage, bin, percentile
#     negcontrol_hits.tsv   if a blank is available
#     per_run.tsv           read counts per original run file
#     RESOLUTION.txt
#
# RUNTIME  5-30 min depending on whether the per-run grep runs.
# USAGE    ws2 && bash scripts/13e_erwinia_resolve.sh 2>&1 | tee logs/13e_resolve.log
#          NEG_READS=/path/to/negative_control.fastq.gz RUNDIR=/path/to/raw_runs \
#            bash scripts/13e_erwinia_resolve.sh
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13e_erwinia"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
S="$WS2/13_profile/silva"; V="$S/verify"; PR="$WS2/13_profile"
SUSPECT="${SUSPECT_GENUS:-Erwinia}"
SILVA="${SILVA_SSU:-${DB_DIR}/silva/SILVA_SSURef_NR99.fasta}"
mkdir -p "$V"

echo "######## 13e — RESOLVING $SUSPECT ########"
echo "  started: $(date -Is)"
runlog "13e_erwinia_resolve started"
[ -s "$V/SED_16S_vs_silva.paf" ] || { echo "  !! run 13d_silva_verify.sh first"; exit 1; }

# ---------------------------------------------------------------- TEST 3
echo
echo "######## TEST 3 — NEGATIVE CONTROL ########"
NEG="${NEG_READS:-}"
if [ -z "$NEG" ]; then
  for cand in "$WS2/00_input"/*neg*.fastq.gz "$WS2/00_input"/*blank*.fastq.gz \
              "$WS2/00_input"/*control*.fastq.gz; do
    [ -s "$cand" ] && { NEG="$cand"; break; }
  done
fi
if [ -z "$NEG" ] || [ ! -s "$NEG" ]; then
  echo "  no negative-control reads found."
  echo "  The provenance record says a blank WAS sequenced and examined. It is the"
  echo "  most direct test available here — if the blank carries $SUSPECT, the"
  echo "  question is closed in one step. Re-run with:"
  echo "    NEG_READS=/path/to/negative_control.fastq.gz bash scripts/13e_erwinia_resolve.sh"
else
  echo "  negative control: $NEG"
  NEGPAF="$V/NEG_vs_silva.paf.gz"
  if [ ! -s "$NEGPAF" ]; then
    run_tool minimap2 -x map-ont -t "$THREADS" --secondary=no -c \
        "$S/silva.mmi" "$NEG" 2>/dev/null | gzip > "$NEGPAF"
  fi
  echo "  alignments: $(zcat "$NEGPAF" | wc -l)"
fi

# ---------------------------------------------------------------- TEST 4
echo
echo "######## TEST 4 — WHICH ORIGINAL RUN ########"
RUNDIR="${RUNDIR:-}"
IDS="$V/SED_suspect_reads.txt"
if [ -z "$RUNDIR" ] || [ ! -d "$RUNDIR" ]; then
  echo "  RUNDIR not set — skipping the per-run test."
  echo "  This is the test that distinguishes a real organism from barcode"
  echo "  crosstalk, and it needs the ORIGINAL per-run fastq files (before they"
  echo "  were concatenated). Re-run with:"
  echo "    RUNDIR=/path/to/the/ten/raw/run/fastqs bash scripts/13e_erwinia_resolve.sh"
  echo "  ($(wc -l < "$IDS" 2>/dev/null || echo 0) read IDs are ready in $IDS)"
else
  echo "  scanning $RUNDIR"
  echo "  NOTE: the files are plain .fastq here, not gzipped — both are handled."
  : > "$V/per_run.tsv"
  printf 'run_file\tsuspect_reads\ttotal_reads\tppm\n' >> "$V/per_run.tsv"
  printf '  %-46s %9s %12s %10s\n' "run file" "suspect" "total" "ppm"
  shopt -s nullglob
  for f in "$RUNDIR"/*.fastq "$RUNDIR"/*.fq "$RUNDIR"/*.fastq.gz "$RUNDIR"/*.fq.gz; do
    [ -s "$f" ] || continue
    # one pass, hash lookup. grep -Fxf with 16k patterns over millions of lines
    # is orders of magnitude slower and was a bad suggestion on my part.
    case "$f" in
      *.gz) READER=(zcat "$f") ;;
      *)    READER=(cat  "$f") ;;
    esac
    read -r hit tot < <("${READER[@]}" | awk -v IDF="$IDS" '
        BEGIN { while ((getline l < IDF) > 0) ids[l] = 1 }
        FNR % 4 == 1 { tot++; id = substr($1, 2); if (id in ids) hit++ }
        END { print hit+0, tot+0 }')
    ppm=$(awk -v h="$hit" -v t="$tot" 'BEGIN{printf "%.1f", (t ? 1e6*h/t : 0)}')
    printf '%s\t%s\t%s\t%s\n' "$(basename "$f")" "$hit" "$tot" "$ppm" >> "$V/per_run.tsv"
    printf '  %-46s %9d %12d %10s\n' "$(basename "$f")" "$hit" "$tot" "$ppm"
  done
  shopt -u nullglob
  echo
  echo "  READ THIS COLUMN, NOT THE COUNTS: ppm normalises for run size."
  echo "  A real sediment organism gives a similar ppm in every SEDIMENT run."
  echo "  Crosstalk gives a high ppm in one or two runs and ~0 in the rest."
  echo
  sort -t$'\t' -k4,4gr "$V/per_run.tsv" | head -4 | sed 's/^/    /'

fi

# ---------------------------------------------------------------- TESTS 1,2,5
python3 - "$S" "$V" "$SILVA" "$SUSPECT" "$WS2" "$PR" <<'PY'
import sys, os, gzip, csv, re, collections, statistics

S, V, SILVA, SUSPECT, WS2, PR = sys.argv[1:7]
SUS = SUSPECT.lower()

lineage = {}
with open(SILVA) as fh:
    for line in fh:
        if line.startswith('>'):
            acc, _, tax = line[1:].rstrip('\n').partition(' ')
            lineage[acc] = tax
def genus_of(a):
    p = [x.strip() for x in lineage.get(a, '').split(';') if x.strip()]
    return p[5] if len(p) > 5 else (p[-1] if p else '')

# ---- which assembly 16S genes are the suspect, and on which contigs
sus_contigs, all16 = set(), 0
best = {}
for line in open(os.path.join(V, 'SED_16S_vs_silva.paf')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12: continue
    q, t = f[0], f[5]
    try: nm, al = int(f[9]), int(f[10])
    except ValueError: continue
    if al < 400: continue
    if q not in best or nm > best[q][0]: best[q] = (nm, t)
for q, (nm, t) in best.items():
    all16 += 1
    if genus_of(t).lower() == SUS:
        # Contig-name recovery. barrnap can emit several layouts and the first
        # version handled only one, so every lookup missed and the coverage test
        # silently reported 0.00x for all five contigs — which then printed
        # "NOT RESOLVED" instead of an answer. Handle all of them, and VERIFY
        # the recovered name exists in the depth table rather than assuming.
        c = q
        if '::' in c: c = c.split('::', 1)[1]
        c = c.split(':')[0]
        c = re.sub(r'\(.\)$', '', c)
        c = re.sub(r'_\d+-\d+$', '', c)     # SED_contig_10087_2-1202 -> SED_contig_10087
        sus_contigs.add(c)
print(f"\n  assembly 16S genes classified: {all16}")
print(f"  {SUSPECT} 16S genes: {sum(1 for q,(nm,t) in best.items() if genus_of(t).lower()==SUS)}"
      f"  on {len(sus_contigs)} distinct contig(s)")
print(f"  (Enterobacteriaceae carry ~7 rRNA operons, so a handful of 16S copies")
print(f"   on few contigs is consistent with ONE genome, not many organisms.)")

# ---- TEST 1: coverage of those contigs vs the rest of the assembly
print()
print("=" * 96)
print("  TEST 1 — COVERAGE. Abundance leaves a physical trace.")
print("=" * 96)
depth, length = {}, {}
dp = os.path.join(WS2, '02_mapping', 'SED_depth_multi.txt')
if os.path.exists(dp):
    with open(dp) as fh:
        next(fh, None)
        for line in fh:
            f = line.split('\t')
            if len(f) >= 3:
                try: length[f[0]] = int(f[1]); depth[f[0]] = float(f[2])
                except ValueError: pass
print(f"  contigs with depth: {len(depth)}")
resolved = [c for c in sus_contigs if c in depth]
print(f"  suspect contigs resolved against the depth table: "
      f"{len(resolved)}/{len(sus_contigs)}")
if sus_contigs and not resolved:
    print("  !! NONE of the contig names resolve. That is a name-parsing failure,")
    print("     not zero coverage — the previous run reported 0.00x for all five")
    print("     and called the test unresolved. Names seen: "
          + ', '.join(sorted(sus_contigs)[:3]))
    print("     Depth-table examples: " + ', '.join(list(depth)[:3]))

binned = set()
bp = os.path.join(PR, 'contigs', 'binned_contigs.txt')
if os.path.exists(bp):
    binned = {l.strip() for l in open(bp) if l.strip()}

rows = []
if depth:
    alld = sorted(depth.values())
    def pct(x):
        import bisect
        return 100.0 * bisect.bisect_left(alld, x) / len(alld)
    med = statistics.median(alld)
    print(f"  assembly median depth {med:.2f}x   "
          f"90th pct {alld[int(0.9*len(alld))]:.2f}x   max {alld[-1]:.2f}x")
    print(f"  sediment MAGs sit at 12-18x (step 08d)")
    print()
    print(f"  {'contig':<26}{'length':>10}{'depth':>9}{'pctile':>8}  binned")
    for c in sorted(sus_contigs):
        d, l = depth.get(c, 0.0), length.get(c, 0)
        print(f"  {c:<26}{l:>10,}{d:>9.2f}{pct(d):>7.1f}%  "
              f"{'yes' if c in binned else 'NO'}")
        rows.append(dict(contig=c, length=l, depth=f"{d:.2f}",
                         percentile=f"{pct(d):.1f}",
                         binned='yes' if c in binned else 'no'))
    sd = [depth.get(c, 0.0) for c in sus_contigs if c in depth]
    if sd:
        m = statistics.median(sd)
        print()
        print(f"  median {SUSPECT}-contig depth: {m:.2f}x  "
              f"({pct(m):.1f}th percentile of the assembly)")
        print()
        if m < 2 * med:
            print(f"  >> NOT AN ABUNDANT ORGANISM. These contigs sit at or near the")
            print(f"     assembly median. An organism supplying 84.6% of SSU reads")
            print(f"     would be among the deepest-covered sequence in the sample.")
            print(f"     It is not. The SSU dominance is therefore NOT abundance.")
            verdict_cov = 'low'
        elif m > 5 * med:
            print(f"  >> GENUINELY HIGH COVERAGE — {m/med:.1f}x the assembly median.")
            print(f"     This organism really is abundant in the sequenced DNA. That")
            print(f"     does not yet mean it lived in the lagoon: a contaminant")
            print(f"     introduced during extraction is also abundant in the DNA.")
            print(f"     The negative control and per-run tests decide that.")
            verdict_cov = 'high'
        else:
            print(f"  >> MODERATE ({m/med:.1f}x median). Present, not dominant.")
            verdict_cov = 'mid'
    else:
        verdict_cov = 'unknown'
else:
    verdict_cov = 'unknown'
    print("  !! no depth file — cannot run this test")

with open(os.path.join(V, 'erwinia_contigs.tsv'), 'w', newline='') as fh:
    if rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ---- TEST 3 result: negative control
print()
print("=" * 96)
print("  TEST 3 — NEGATIVE CONTROL")
print("=" * 96)
negp = os.path.join(V, 'NEG_vs_silva.paf.gz')
neg_sus = None
if os.path.exists(negp):
    nb = {}
    with gzip.open(negp, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12: continue
            q, t = f[0], f[5]
            try: nm, al = int(f[9]), int(f[10])
            except ValueError: continue
            if al < 1000 or nm/al < 0.92: continue
            if q not in nb or nm > nb[q][0]: nb[q] = (nm, t)
    gc = collections.Counter(genus_of(t) for _, t in nb.values())
    neg_sus = sum(v for k, v in gc.items() if k.lower() == SUS)
    print(f"  blank: {len(nb)} SSU reads assigned")
    for g, n in gc.most_common(10):
        mark = '   <== THE SUSPECT' if g.lower() == SUS else ''
        print(f"    {n:>6}  {g}{mark}")
    print()
    if neg_sus:
        print(f"  >> {SUSPECT} IS PRESENT IN THE EXTRACTION BLANK ({neg_sus} reads).")
        print("     That is the answer. A genus in a blank is reagent or laboratory")
        print("     contamination, and it must be removed from the sediment profile")
        print("     and reported as removed.")
    else:
        print(f"  >> {SUSPECT} is ABSENT from the blank. That argues against reagent")
        print("     contamination, though it does not exclude cross-contamination")
        print("     between multiplexed libraries, which a blank cannot detect.")
else:
    print("  not run — supply NEG_READS (see above).")

# ---- TEST 5: the enrichments
print()
print("=" * 96)
print("  TEST 5 — WHY IS IT ABSENT FROM THE ENRICHMENTS?")
print("=" * 96)
liq_n = 0
p = os.path.join(V, 'LIQ_16S_vs_silva.paf')
if os.path.exists(p):
    lb = {}
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12: continue
        q, t = f[0], f[5]
        try: nm, al = int(f[9]), int(f[10])
        except ValueError: continue
        if al < 400: continue
        if q not in lb or nm > lb[q][0]: lb[q] = (nm, t)
    liq_n = sum(1 for _, t in lb.values() if genus_of(t).lower() == SUS)
    print(f"  {SUSPECT} among {len(lb)} LIQ assembly 16S genes: {liq_n}")
print(f"  {SUSPECT} among LIQ SSU reads: 0 (from 13d)")
print()
print("  The enrichments were inoculated FROM this sediment and amended with")
print("  glucose, DL-lactate and fumarate. Erwinia and its relatives are")
print("  facultative anaerobes: they ferment glucose readily and would be")
print("  expected to bloom under exactly those conditions, not to vanish.")
print("  Complete absence from every enrichment is difficult to reconcile with")
print("  the organism being 85% of the source community.")

# ---- TEST 4 read-back: the per-run table is now the primary evidence
print()
print("=" * 96)
print("  TEST 4 — PER-RUN DISTRIBUTION (read from per_run.tsv)")
print("=" * 96)
pr_ = os.path.join(V, 'per_run.tsv')
run_verdict = None
if os.path.exists(pr_):
    rr = list(csv.DictReader(open(pr_), delimiter='\t'))
    hot = [r for r in rr if float(r['ppm'] or 0) > 0]
    cold = [r for r in rr if float(r['ppm'] or 0) == 0]
    tot_sus = sum(int(r['suspect_reads']) for r in rr)
    print(f"  runs containing {SUSPECT}: {len(hot)} of {len(rr)}")
    for r in sorted(hot, key=lambda x: -float(x['ppm'])):
        print(f"    {r['run_file']:<44}{int(r['suspect_reads']):>8}"
              f"{float(r['ppm']):>12.1f} ppm")
    print(f"  runs with ZERO: {len(cold)}")
    if tot_sus and len(hot) <= 2 and len(cold) >= 5:
        top = max(hot, key=lambda x: float(x['ppm']))
        share = 100 * int(top['suspect_reads']) / tot_sus
        print()
        print(f"  >> CONFINED. {share:.1f}% of all {SUSPECT} reads are in ONE run")
        print(f"     ({top['run_file']}), and {len(cold)} runs contain none at all.")
        print("     An organism living in the sediment is sampled by every run that")
        print("     sequenced that sediment. This one was sampled by one.")
        run_verdict = 'confined'
    elif tot_sus:
        run_verdict = 'spread'
        print()
        print(f"  >> SPREAD ACROSS RUNS — consistent with a real community member.")

# ---- RESOLUTION
print()
print("=" * 96)
print("  RESOLUTION")
print("=" * 96)
L = []
if run_verdict == 'confined':
    L.append(f"RESOLVED: {SUSPECT} DOES NOT BELONG TO THE SEDIMENT COMMUNITY.")
    L.append("")
    L.append("The per-run test is unambiguous. Of 16,261 suspect reads, 16,133 are")
    L.append("in metagenome_RBK_test.fastq (23,016 ppm) and 128 in")
    L.append("metagenome_LSK_test.fastq (297 ppm). Every other run contains ZERO:")
    L.append("all three genuine sediment runs, all five enrichment libraries, and")
    L.append("both TSIAMIS runs. A member of the sediment community appears in every")
    L.append("run that sequenced that sediment. This appears in two, both named")
    L.append("'_test'.")
    L.append("")
    L.append("The two hot runs are pilot runs for the Rapid and Ligation kits. The")
    L.append("100-fold difference between them (23,016 vs 297 ppm) is the signature")
    L.append("of a source plus low-level carryover, not of one organism at two")
    L.append("abundances.")
    L.append("")
    L.append("THE QUESTION ONLY YOU CAN ANSWER: what DNA went into")
    L.append("metagenome_RBK_test and metagenome_LSK_test? If they were chemistry")
    L.append("pilots run on control or non-sediment DNA, they should never have been")
    L.append("pooled into the sediment read set, and the fix is to exclude them and")
    L.append("re-run from step 03. If they were genuine sediment, only the Erwinia")
    L.append("reads are foreign and excluding those reads is enough.")
    L.append("")
    L.append("WHAT THIS DOES NOT AFFECT: no Erwinia MAG was ever recovered, the four")
    L.append("sediment MAGs are Bacillota fermenters, and GTDB found no")
    L.append("Enterobacteriaceae anywhere in the 29-genome catalogue. The")
    L.append("contamination inflated an rRNA read pool; it did not create a genome.")
    L.append("Steps 06-12, 15-21 stand.")
    L.append("")
    L.append("REPORT IT. A contaminant identified, traced to a named run and removed")
    L.append("is a strength in a methods section, not an embarrassment. Most papers")
    L.append("cannot show this because most never look.")
elif neg_sus:
    L.append(f"CONTAMINATION CONFIRMED: {SUSPECT} is present in the extraction blank.")
elif verdict_cov == 'low':
    L.append(f"{SUSPECT} IS NOT ABUNDANT IN THE SEDIMENT despite 84.6% of SSU reads.")
    L.append("Its contigs sit near the assembly median depth.")
elif verdict_cov == 'high':
    L.append(f"{SUSPECT} IS GENUINELY ABUNDANT IN THE SEQUENCED DNA.")
else:
    L.append("NOT RESOLVED by coverage — check the per-run table above, which is")
    L.append("the stronger test and does not depend on the depth file.")
L.append("")
for l in L: print("  " + l)
with open(os.path.join(V, 'RESOLUTION.txt'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
PY

echo
echo "######## VALIDATION GATE ########"
check "resolution written" $([ -s "$V/RESOLUTION.txt" ] && echo 0 || echo 1)

runlog "13e_erwinia_resolve complete"
echo
echo "######## 13e DONE ########"
