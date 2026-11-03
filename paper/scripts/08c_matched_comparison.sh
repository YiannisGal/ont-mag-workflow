#!/usr/bin/env bash
# =============================================================================
# 08c_matched_comparison.sh  --  did polishing actually change anything?
#
# WHY
#   The workspace_2 rebuild was justified by one inference: workspace_1's MAGs
#   had a 244 aa median gene length, and since `liquid_1` was a CLOSED 4.90 Mb
#   contig at 281 aa, fragmentation could not explain it — so the cause had to
#   be frameshifts from --nano-raw plus missing Medaka polishing.
#
#   That genome has now been reassembled with --nano-hq AND Medaka-polished.
#   Result: 280.78 aa, coding density 0.805 — versus 281 aa / 0.805 before.
#   No change. The inference was wrong: 281 aa is most likely this organism's
#   real mean gene length, and a closed high-coverage genome is precisely where
#   Flye's internal polishing already suffices.
#
#   But that is n=1, and n=1 is what produced the bad inference in the first
#   place. This script does the comparison properly across every genome present
#   in BOTH workspaces.
#
# METHOD
#   FastANI between workspace_2's MAGs and workspace_1's MAGs. Pairs at >=99%
#   ANI are the same organism, so their gene-length difference isolates the
#   effect of the assembly/polishing change — bin composition is held constant
#   because it is literally the same genome either side.
#
# THE TWO OUTCOMES, both of which are publishable
#   Delta ~0 across matched pairs  -> polishing did not fix gene length; the
#     workspace_1 deficit was fragmentation/incompleteness in the bin set, and
#     workspace_2's improvement comes from better binning. The re-run is still
#     justified (correct error model, GUNC screening, provenance, consensus)
#     but the STATED REASON in the audit must be corrected.
#   Delta rises as coverage falls  -> polishing worked where it could, and the
#     closed genome was simply already clean. The original hypothesis survives
#     in a restricted, better-specified form.
#
# workspace_1 IS READ-ONLY. This script only reads from it.
#
# USAGE
#   ws2 && bash scripts/08c_matched_comparison.sh 2>&1 | tee logs/08c_matched.log
#   WS1_BINS="/path/a /path/b" bash scripts/08c_matched_comparison.sh   # override
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="08c_matched"
WS1="${WS1}"
Q="$WS2/06_binqc"; M="$Q/matched"; mkdir -p "$M"

echo "######## 08c — MATCHED-GENOME COMPARISON (workspace_1 vs workspace_2) ########"
echo "  started: $(date -Is)"
runlog "08c_matched_comparison started"

# --------------------------------------------------- locate workspace_1 bins
if [ -n "${WS1_BINS:-}" ]; then
  DIRS="$WS1_BINS"
else
  DIRS=$(find "$WS1" -maxdepth 6 -type d \( -name '*DASTool_bins' -o -name '*dastool*bins*' \) \
         2>/dev/null | tr '\n' ' ')
fi
echo
echo "  workspace_1 bin directories found:"
if [ -z "${DIRS// /}" ]; then
  echo "    NONE. Re-run with an explicit list, e.g."
  echo "    WS1_BINS=\"$WS1/liquid_cultures/binning/dastool_output/x_DASTool_bins\" bash scripts/08c_matched_comparison.sh"
  echo
  echo "  To find them:"
  echo "    find $WS1 -name '*.fa' -path '*bin*' 2>/dev/null | sed 's#/[^/]*\$##' | sort -u | head -30"
  exit 1
fi
for d in $DIRS; do printf '    %s  (%s bins)\n' "$d" "$(ls -1 "$d"/*.fa "$d"/*.fasta 2>/dev/null | wc -l)"; done

# ------------------------------------------------------------- build lists
ls -1 "$Q/union_bins"/*.fa > "$M/ws2_list.txt"
: > "$M/ws1_list.txt"
for d in $DIRS; do ls -1 "$d"/*.fa "$d"/*.fasta 2>/dev/null >> "$M/ws1_list.txt"; done
sort -u "$M/ws1_list.txt" -o "$M/ws1_list.txt"
echo
printf '  workspace_2 genomes: %s\n' "$(wc -l < "$M/ws2_list.txt")"
printf '  workspace_1 genomes: %s\n' "$(wc -l < "$M/ws1_list.txt")"
[ -s "$M/ws1_list.txt" ] || { echo "  !! no workspace_1 bins to compare"; exit 1; }

# ------------------------------------------------------------------ FastANI
echo
echo "######## FastANI ($(date +%H:%M)) ########"
if [ -s "$M/ani.tsv" ]; then
  echo "  already done"
else
  record_version "fastANI" "${ENV_fastANI:-coverm_env}" fastANI --version
  run_tool fastANI --ql "$M/ws2_list.txt" --rl "$M/ws1_list.txt" \
      -o "$M/ani.tsv" -t "$THREADS" --matrix 2>&1 | tail -5 \
    || echo "  !! fastANI failed"
fi
[ -s "$M/ani.tsv" ] || { echo "  !! no ANI output"; exit 1; }
echo "  ANI hits: $(wc -l < "$M/ani.tsv")"

# ------------------------------------------------- collect workspace_1 AGL
echo
echo "######## COLLECTING workspace_1 CheckM2 RESULTS ########"
: > "$M/ws1_checkm2.tsv"
n=0
while IFS= read -r qr; do
  [ -s "$qr" ] || continue
  awk -F'\t' -v src="$qr" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    { print $h["Name"]"\t"$h["Completeness"]"\t"$h["Contamination"]"\t"\
            $h["Average_Gene_Length"]"\t"$h["Coding_Density"]"\t"src }' "$qr" \
    >> "$M/ws1_checkm2.tsv"
  n=$((n+1))
done < <(find "$WS1" -name 'quality_report.tsv' 2>/dev/null)
echo "  merged $n workspace_1 quality_report.tsv files -> $(wc -l < "$M/ws1_checkm2.tsv") rows"

# ------------------------------------------------------------------ COMPARE
python3 - "$M" "$Q/MAG_QC.tsv" <<'PY'
import sys, os, csv, statistics as st

M, qc_path = sys.argv[1], sys.argv[2]
ANI_MIN, AF_MIN = 99.0, 0.50

ws2 = {}
for r in csv.DictReader(open(qc_path), delimiter='\t'):
    ws2[r['bin']] = r

ws1 = {}
for line in open(os.path.join(M, 'ws1_checkm2.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) < 5 or not p[0]: continue
    ws1.setdefault(p[0], p)          # first occurrence wins

def stem(p):
    b = os.path.basename(p)
    for ext in ('.fa', '.fasta', '.fna'):
        if b.endswith(ext): return b[:-len(ext)]
    return b

# best workspace_1 match per workspace_2 genome
best = {}
for line in open(os.path.join(M, 'ani.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) < 5: continue
    q, r, ani, frag, tot = stem(p[0]), stem(p[1]), float(p[2]), int(p[3]), int(p[4])
    af = frag / tot if tot else 0
    if ani < ANI_MIN or af < AF_MIN: continue
    if q not in best or ani > best[q][1]:
        best[q] = (r, ani, af)

print()
print("=" * 78)
print("  MATCHED GENOMES (>=99% ANI, >=50% aligned) — same organism, both workspaces")
print("=" * 78)
if not best:
    print("  No matches. Either the bin sets share no genomes at >=99% ANI, or the")
    print("  workspace_1 directories located were not the MAG directories.")
    raise SystemExit(0)

print(f"  {'workspace_2 MAG':<26} {'workspace_1 MAG':<22} {'ANI':>6} "
      f"{'AGL w1':>7} {'AGL w2':>7} {'delta':>7} {'CD w1':>6} {'CD w2':>6} {'comp w2':>8}")
deltas = []
for q in sorted(best):
    r, ani, af = best[q]
    w1 = ws1.get(r)
    w2 = ws2.get(q)
    if not w1 or not w2: continue
    try:
        a1, a2 = float(w1[3]), float(w2['avg_gene_length'])
    except (ValueError, TypeError, KeyError):
        continue
    d = a2 - a1
    deltas.append((d, float(w2['completeness'] or 0)))
    print(f"  {q:<26} {r:<22} {ani:>6.2f} {a1:>7.1f} {a2:>7.1f} {d:>+7.1f} "
          f"{w1[4][:6]:>6} {(w2['coding_density'] or '')[:6]:>6} "
          f"{(w2['completeness'] or '')[:6]:>8}")

print()
if deltas:
    ds = [d for d, _ in deltas]
    print(f"  matched pairs           : {len(ds)}")
    print(f"  median change in AGL    : {st.median(ds):+.1f} aa")
    print(f"  mean change in AGL      : {st.mean(ds):+.1f} aa")
    print(f"  pairs improved >10 aa   : {sum(1 for d in ds if d > 10)}")
    print(f"  pairs unchanged (+-10)  : {sum(1 for d in ds if abs(d) <= 10)}")
    print(f"  pairs worse   <-10 aa   : {sum(1 for d in ds if d < -10)}")
    print()
    med = st.median(ds)
    print("  " + "-" * 74)
    if abs(med) <= 10:
        print("  >> POLISHING DID NOT CHANGE GENE LENGTH on the same genomes.")
        print("     The workspace_1 244 aa deficit was therefore NOT frameshift error.")
        print("     It was fragmentation and incompleteness in that bin set, and")
        print("     workspace_2's improvement comes from better BINNING.")
        print()
        print("     ACTION: correct the stated rationale in pipeline_audit_and_roadmap.md")
        print("     and in the lab notebook. The re-run remains justified on other")
        print("     grounds (correct error model as best practice, three-binner")
        print("     consensus, GUNC screening, contig provenance, reproducible")
        print("     record) but NOT on the gene-length argument.")
    elif med > 10:
        print("  >> POLISHING IMPROVED GENE LENGTH on matched genomes.")
        print("     The original hypothesis survives. Check the completeness column:")
        print("     if gains concentrate in less-complete (lower-coverage) genomes,")
        print("     state the effect as coverage-dependent, which is the mechanism")
        print("     Medaka would predict.")
    else:
        print("  >> GENE LENGTH DECREASED. Unexpected — investigate before reporting.")
        print("     Check that matched pairs really are the same organism and that")
        print("     the workspace_1 CheckM2 rows were parsed from the right files.")

    lowc  = [d for d, c in deltas if c < 90]
    highc = [d for d, c in deltas if c >= 90]
    if lowc and highc:
        print()
        print(f"  by completeness:  >=90% complete  median {st.median(highc):+.1f} aa (n={len(highc)})")
        print(f"                    <90%  complete  median {st.median(lowc):+.1f} aa (n={len(lowc)})")
        print("  (Medaka's benefit should be larger at lower coverage. If these two")
        print("   are similar, coverage-dependence is not supported either.)")

with open(os.path.join(M, 'matched_pairs.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['ws2_mag','ws1_mag','ani','align_frac','agl_ws1','agl_ws2',
                'delta_agl','cd_ws1','cd_ws2','completeness_ws2'])
    for q in sorted(best):
        r, ani, af = best[q]
        w1, w2 = ws1.get(r), ws2.get(q)
        if not w1 or not w2: continue
        w.writerow([q, r, f"{ani:.2f}", f"{af:.2f}", w1[3],
                    w2['avg_gene_length'],
                    (float(w2['avg_gene_length']) - float(w1[3]))
                    if w1[3] and w2['avg_gene_length'] else '',
                    w1[4], w2['coding_density'], w2['completeness']])
print()
print(f"  written: {os.path.join(M, 'matched_pairs.tsv')}")
PY

runlog "08c_matched_comparison complete"
echo
echo "######## 08c DONE ########"
