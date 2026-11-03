#!/usr/bin/env bash
# =============================================================================
# 08d_finalize_magset.sh  --  freeze the MAG set, fix naming, measure coverage
#
# THREE CLEANUP JOBS, then everything downstream has ONE unambiguous input.
#
#   1. RESTORE BINNER PROVENANCE for DAS_Tool's five `_sub` bins.
#      DAS_Tool appends `_sub` when it MODIFIES a bin by removing contigs to
#      raise its score. The step-08 name parser did not allow the suffix and
#      labelled them binner=unknown. The requirement is that every filename
#      states which binner produced it, so these are renamed
#        <SAMPLE>_dassub_<binner>_<NNN>.fa
#      which keeps the binner readable AND records that DAS_Tool refined it —
#      information the original name lost. Unrefined bins keep
#        <SAMPLE>_das_<binner>_<NNN>.fa
#
#   2. FREEZE THE 29-MAG SET.
#      Selection, fixed at step 08: DAS_Tool score_threshold 0.3, MIMAG high or
#      medium (CheckM2 >=50% complete, <10% contaminated), GUNC clade separation
#      score <= 0.45. Copied to 07_mags/bins/ with a crosswalk carrying
#      publication IDs (ETO_<SAMPLE>_MAG_NNN), working names, source binner,
#      raw binner filename, and every QC statistic.
#      workspace_1 had three conflicting sets (38/22/17) with no recorded
#      criteria. This is the file that prevents a repeat.
#
#   3. PER-MAG COVERAGE with CoverM.
#      Needed for two reasons: completeness proved a weak proxy for coverage and
#      left the coverage-dependence question at step 08c unresolved; and the
#      abundance discussion needs real numbers. Each MAG is quantified against
#      BOTH its own sample's reads AND the other sample's reads, using the self
#      and cross BAMs from step 05 — so a sediment MAG's coverage in culture
#      reads directly answers "did this organism enrich?".
#
# RUNTIME  10-30 min, CoverM dominating.
# USAGE    ws2 && bash scripts/08d_finalize_magset.sh 2>&1 | tee logs/08d_finalize.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="08d_finalize"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
Q="$WS2/06_binqc"; F="$WS2/07_mags"; M="$WS2/02_mapping"
mkdir -p "$F/bins" "$F/coverage"

echo "######## 08d — FINALISE MAG SET ########"
echo "  started: $(date -Is)"
runlog "08d_finalize_magset started"

for f in "$Q/MAG_QC.tsv" "$Q/union_membership.tsv"; do
  [ -s "$f" ] || { echo "  !! missing $f — run 08 first"; exit 1; }
done

# ============================================ 1+2. RENAME AND FREEZE
echo
echo "######## SELECTING AND RENAMING ########"
python3 - "$Q" "$F" <<'PY'
import sys, os, csv, shutil, re

Q, F = sys.argv[1], sys.argv[2]

qc   = {r['bin']: r for r in csv.DictReader(open(f"{Q}/MAG_QC.tsv"), delimiter='\t')}
memb = {r['bin']: r for r in csv.DictReader(open(f"{Q}/union_membership.tsv"), delimiter='\t')}

def fnum(x, d=0.0):
    try: return float(x)
    except: return d

# ---- selection criteria, as fixed at step 08 -------------------------------
sel = []
for b, r in qc.items():
    if r.get('in_s0.3') != 'yes':            continue   # adopted threshold
    if r['mimag'] not in ('high', 'medium'): continue   # MIMAG medium or better
    if r['gunc_chimeric'] == 'yes':          continue   # GUNC CSS > 0.45
    sel.append(r)
print(f"  selected {len(sel)} MAGs "
      f"(score 0.3 + MIMAG high/medium + GUNC clean) from {len(qc)} union bins")

# ---- recover the true binner for DAS_Tool `_sub` bins ----------------------
SRC = re.compile(r'/(SED|LIQ)_(metabat|vamb|semibin)_(\d+)(_sub)?\.fa$')
fixed = 0
for r in sel:
    src = memb[r['bin']].get('source_file', '')
    m = SRC.search(src)
    if not m:
        r['_binner'], r['_sub'], r['_srcnum'] = r['binner'], False, '000'
        continue
    r['_binner'], r['_srcnum'], r['_sub'] = m.group(2), m.group(3), bool(m.group(4))
    if r['binner'] in ('unknown', '') and r['_binner'] != 'unknown':
        fixed += 1
        print(f"    provenance restored: {r['bin']:<24} -> {r['_binner']}"
              f"{' (DAS_Tool-refined subset)' if r['_sub'] else ''}")
print(f"  binner provenance restored for {fixed} bin(s)")

# ---- final names + publication IDs ----------------------------------------
sel.sort(key=lambda r: (r['sample'] != 'SED', -fnum(r['completeness'])))
rows, counters = [], {}
for r in sel:
    s = r['sample']
    stage = 'dassub' if r['_sub'] else 'das'
    counters[s] = counters.get(s, 0) + 1
    pubid = f"ETO_{s}_MAG_{counters[s]:03d}"
    name  = f"{s}_{stage}_{r['_binner']}_{int(r['_srcnum']):03d}"
    n = 1
    while any(x['working_name'] == name for x in rows):     # collision guard
        n += 1; name = f"{s}_{stage}_{r['_binner']}_{int(r['_srcnum']):03d}{chr(96+n)}"
    src = os.path.join(Q, 'union_bins', r['bin'] + '.fa')
    dst = os.path.join(F, 'bins', name + '.fa')
    shutil.copyfile(src, dst)
    rows.append(dict(
        mag_id=pubid, working_name=name, sample=s, binner=r['_binner'],
        das_refined=('yes' if r['_sub'] else 'no'),
        qc_bin_name=r['bin'],
        raw_binner_file=os.path.basename(memb[r['bin']].get('source_file', '')),
        n_contigs=r['n_contigs'], bp=r['bp'],
        completeness=r['completeness'], contamination=r['contamination'],
        mimag=r['mimag'], avg_gene_length=r['avg_gene_length'],
        coding_density=r['coding_density'], gunc_css=r['gunc_css'],
        n_trna_types=r['n_trna_types'], has_5S=r['has_5S'],
        has_16S=r['has_16S'], has_23S=r['has_23S'],
        in_s05=r.get('in_s0.5',''), in_s04=r.get('in_s0.4',''), in_s03=r.get('in_s0.3','')))

with open(os.path.join(F, 'MAG_CROSSWALK.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)

print()
print(f"  {'MAG ID':<18} {'working name':<28} {'MIMAG':<7} {'comp':>6} {'cont':>6} {'Mb':>6}")
for r in rows:
    print(f"  {r['mag_id']:<18} {r['working_name']:<28} {r['mimag']:<7} "
          f"{r['completeness'][:6]:>6} {r['contamination'][:6]:>6} "
          f"{int(r['bp'])/1e6:>6.2f}")
print()
for s in ('SED', 'LIQ'):
    sub = [r for r in rows if r['sample'] == s]
    hi = sum(1 for r in sub if r['mimag'] == 'high')
    print(f"  {s}: {len(sub)} MAGs ({hi} high, {len(sub)-hi} medium), "
          f"{sum(int(r['bp']) for r in sub)/1e6:.1f} Mbp, "
          f"{sum(1 for r in sub if r['das_refined']=='yes')} DAS_Tool-refined")
PY

NMAG=$(ls -1 "$F/bins"/*.fa 2>/dev/null | wc -l)
echo "  frozen MAG set: $NMAG genomes in $F/bins/"

# ============================================ 3. COVERAGE (CoverM)
echo
echo "######## CoverM — per-MAG coverage, self and cross ($(date +%H:%M)) ########"
record_version "coverm" "${ENV_coverm:-coverm_env}" coverm --version

cover() {  # cover <SAMPLE>
  local s="$1"
  local out="$F/coverage/${s}_coverage.tsv"
  [ -s "$out" ] && { echo "  $s: already done"; return 0; }
  # a MAG's contigs exist only in its own sample's BAMs, so quantify per sample
  local d="$TMPDIR/coverm_$s"; rm -rf "$d"; mkdir -p "$d"
  local n=0
  for f in "$F/bins/${s}_"*.fa; do [ -s "$f" ] && { cp "$f" "$d/"; n=$((n+1)); }; done
  [ "$n" -gt 0 ] || { echo "  $s: no MAGs"; return 0; }
  local self="$M/${s}_map_self.sorted.bam" cross="$M/${s}_map_cross.sorted.bam"
  for b in "$self" "$cross"; do
    [ -s "$b" ] || { echo "  !! missing $b — run 05_map.sh"; return 1; }
  done
  echo "  $s: $n MAGs vs self + cross BAMs"
  run_tool coverm genome --genome-fasta-directory "$d" -x fa \
      --bam-files "$self" "$cross" \
      -m mean covered_fraction relative_abundance \
      --min-covered-fraction 0 -t "$THREADS" -o "$out" \
    || echo "    !! coverm failed for $s"
  rm -rf "$d"
}
cover SED
cover LIQ

echo
echo "######## COVERAGE SUMMARY (self vs cross) ########"
python3 - "$F" <<'PY'
import sys, os, csv
F = sys.argv[1]
print(f"  {'MAG':<28} {'mean cov (own)':>15} {'mean cov (other)':>17} "
      f"{'breadth own':>12} {'enriched?':>10}")
for s in ('SED', 'LIQ'):
    p = os.path.join(F, 'coverage', f'{s}_coverage.tsv')
    if not os.path.exists(p): continue
    rd = list(csv.reader(open(p), delimiter='\t'))
    if len(rd) < 2: continue
    hdr, body = rd[0], rd[1:]
    # column order: Genome, then per-BAM blocks (self first, then cross)
    mean_i = [i for i, h in enumerate(hdr) if h.endswith('Mean')]
    br_i   = [i for i, h in enumerate(hdr) if 'Covered Fraction' in h]
    if len(mean_i) < 2: continue
    for row in body:
        if not row or row[0].lower().startswith('unmapped'): continue
        try:
            own, other = float(row[mean_i[0]]), float(row[mean_i[1]])
            br = float(row[br_i[0]]) if br_i else float('nan')
        except (ValueError, IndexError):
            continue
        ratio = (other / own) if own > 0 else 0
        tag = 'both' if ratio > 0.1 else ('sample-specific' if ratio < 0.01 else 'minor')
        print(f"  {row[0]:<28} {own:>15.2f} {other:>17.2f} {br:>12.3f} {tag:>10}")
print()
print("  'own'   = reads from the same sample; 'other' = the opposite sample.")
print("  A SED MAG with high coverage in LIQ reads enriched in culture.")
print("  A LIQ MAG absent from SED reads was below detection in the sediment pool.")
PY

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "MAG set frozen ($NMAG genomes)" $([ "$NMAG" -gt 0 ] && echo 0 || echo 1)
check "crosswalk written" $([ -s "$F/MAG_CROSSWALK.tsv" ] && echo 0 || echo 1)
if [ -s "$F/MAG_CROSSWALK.tsv" ]; then
  nrow=$(( $(wc -l < "$F/MAG_CROSSWALK.tsv") - 1 ))
  check "crosswalk covers every MAG ($nrow/$NMAG)" $([ "$nrow" -eq "$NMAG" ] && echo 0 || echo 1)
  nunk=$(awk -F'\t' 'NR>1 && $4=="unknown"' "$F/MAG_CROSSWALK.tsv" | wc -l)
  check "every MAG has a known binner ($nunk unknown)" $([ "$nunk" -eq 0 ] && echo 0 || echo 1)
fi
bad=$(ls -1 "$F/bins"/*.fa 2>/dev/null \
  | grep -vcE '/(SED|LIQ)_(das|dassub)_(metabat|vamb|semibin)_[0-9]{3}[a-z]?\.fa$' || true)
check "MAG filenames follow the convention (${bad} violations)" \
  $([ "${bad:-0}" -eq 0 ] && echo 0 || echo 1)
for s in SED LIQ; do
  wrong=$(cat "$F/bins/${s}_"*.fa 2>/dev/null | grep '^>' | grep -cv "^>${s}_" || true)
  check "$s MAGs contain only ${s}_ contigs (${wrong} foreign)" \
    $([ "${wrong:-0}" -eq 0 ] && echo 0 || echo 1)
  check "$s coverage table" $([ -s "$F/coverage/${s}_coverage.tsv" ] && echo 0 || echo 1)
done
dup=$(cat "$F/bins"/*.fa 2>/dev/null | grep '^>' | sort | uniq -d | wc -l)
check "no contig appears in two MAGs (${dup} duplicates)" \
  $([ "$dup" -eq 0 ] && echo 0 || echo 1)

runlog "08d_finalize_magset complete ($NMAG MAGs)"
echo
echo "######## 08d DONE — next: 09_taxonomy.sh ########"
echo "  Frozen MAG set : $F/bins/            <- the ONLY input downstream"
echo "  Crosswalk      : $F/MAG_CROSSWALK.tsv"
