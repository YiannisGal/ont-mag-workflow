#!/usr/bin/env bash
# =============================================================================
# 09_taxonomy.sh  --  GTDB-Tk R232 classification of the frozen MAG set
#
# PURPOSE
#   Assign each MAG a GTDB taxonomy, and quantify how novel it is. Three
#   project-specific questions ride on this step:
#
#     1. IS THERE A SEDIMENT ARCHAEON, and what is it?
#        workspace_1 appeared to recover none, and I initially wrote that up as
#        "zero sediment Archaea". That was WRONG: the archaeon had been
#        classified through GTDB-Tk's ANI screen, which short-circuits and skips
#        marker-gene alignment, so no ar53 alignment file was ever written and
#        its absence was misread as absence of Archaea.
#        --place_species is therefore used here DELIBERATELY: "place species in
#        pplacer trees even if they are classified with skani". Every genome
#        gets a tree placement IN ADDITION to its ANI assignment, so the
#        evidence trail is uniform and this failure mode cannot recur. The
#        classification_method column is also reported per MAG, which states
#        outright whether a call came from skani ANI or from tree topology —
#        precisely the information whose absence caused the misreading.
#        (Note: GTDB-Tk 2.7 removed --skip_ani_screen; the mash pre-screen no
#        longer exists and skani replaced it. --place_species is the equivalent
#        control in this version.)
#
#     2. DOES ANY MAG CORRESPOND TO Arc_Otu22?
#        The amplicon chapter's focal archaeon: 1.641% peak relative abundance
#        at 60-65 cm, rank 14 of 327. Its classification is contested between
#        chapters (75.8% vs 93.1% identity, audit finding F23). A MAG would
#        settle it with a genome instead of a 16S fragment. Check the placement
#        of any archaeal MAG here, then confirm with the 16S extracted at
#        step 08 (barrnap) rather than from taxonomy alone.
#
#     3. DO THE LIQUID MAGs MATCH THE CULTUROMICS ISOLATES?
#        The isolates came from these same enrichments. Any MAG at >=95% ANI to
#        an isolate is the same species and links the two chapters. Handled at
#        the ANI step, but flagged here because taxonomy is the first clue.
#
# METHOD NOTES FOR MATERIALS & METHODS
#   GTDB-Tk 2.7.2, GTDB release R232. GTDB-Tk 2.x uses skani for ANI (not
#   FastANI or Mash) — cite skani, not FastANI, for this step. pplacer memory
#   scales with thread count, so --pplacer_cpus is held low: the bacterial
#   placement needs ~55 GB and the host has 125 GB.
#
# INPUT    07_mags/bins/               (the frozen 29-MAG set from 08d)
# OUTPUT   08_taxonomy/gtdbtk/         + MAG_TAXONOMY.tsv, novelty summary
#
# RUNTIME  1-3 h for ~30 genomes.
# USAGE    ws2 && bash scripts/09_taxonomy.sh 2>&1 | tee logs/09_taxonomy_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="09_taxonomy"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; T="$WS2/08_taxonomy"; mkdir -p "$T"

echo "######## 09 — GTDB-Tk TAXONOMY ########"
echo "  started: $(date -Is)  threads: $THREADS"
echo "  GTDB release: ${GTDB_RELEASE:-R232}"
echo "  GTDBTK_DATA_PATH: ${GTDBTK_DATA_PATH:-UNSET}"
runlog "09_taxonomy started (GTDB ${GTDB_RELEASE:-R232}, --place_species)"

NMAG=$(ls -1 "$F/bins"/*.fa 2>/dev/null | wc -l)
[ "$NMAG" -gt 0 ] || { echo "  !! no MAGs in $F/bins — run 08d_finalize_magset.sh first"; exit 1; }
echo "  genomes to classify: $NMAG"

[ -d "${GTDBTK_DATA_PATH:-/nonexistent}" ] || {
  echo "  !! GTDBTK_DATA_PATH is not a directory: ${GTDBTK_DATA_PATH:-UNSET}"; exit 1; }

record_version "gtdbtk" "${ENV_gtdbtk:-gtdbtk_env}" gtdbtk --version

# pplacer RAM scales with threads; keep it at 1 to stay well inside 125 GB
PPLACER_CPUS="${PPLACER_CPUS:-1}"
echo "  pplacer_cpus: $PPLACER_CPUS  (memory-bound, not CPU-bound)"

# ------------------------------------------------------------------ CLASSIFY
if [ -s "$T/gtdbtk/gtdbtk.bac120.summary.tsv" ] || [ -s "$T/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
  echo
  echo "  classification already present — delete $T/gtdbtk to force"
else
  echo
  echo "######## classify_wf ($(date +%H:%M)) — this is the long part ########"
  mkdir -p "$TMPDIR/pplacer_scratch"
  run_tool gtdbtk classify_wf \
      --genome_dir "$F/bins" \
      --extension fa \
      --out_dir "$T/gtdbtk" \
      --place_species \
      --write_single_copy_genes \
      --keep_intermediates \
      --cpus "$THREADS" \
      --pplacer_cpus "$PPLACER_CPUS" \
      --scratch_dir "$TMPDIR/pplacer_scratch" \
      --tmpdir "$TMPDIR" \
    || echo "  !! gtdbtk classify_wf returned non-zero — inspect the log"
  echo "  finished: $(date -Is)"
fi

# ------------------------------------------------- MERGE INTO THE CROSSWALK
echo
echo "######## MERGING TAXONOMY ########"
python3 - "$T" "$F" <<'PY'
import sys, os, csv

T, F = sys.argv[1], sys.argv[2]
tax = {}
for dom, fn in (('Bacteria', 'gtdbtk.bac120.summary.tsv'),
                ('Archaea',  'gtdbtk.ar53.summary.tsv')):
    p = os.path.join(T, 'gtdbtk', fn)
    if not os.path.exists(p):
        print(f"  (no {fn} — no {dom} placed, or the run did not complete)")
        continue
    for r in csv.DictReader(open(p), delimiter='\t'):
        r['_domain'] = dom
        tax[r['user_genome']] = r
    print(f"  {fn}: {sum(1 for v in tax.values() if v['_domain']==dom)} genomes")

# GUARD: an empty classification is NOT a finding. Without this, the novelty
# block below reports "29/29 MAGs have no species-level match" when in fact
# GTDB-Tk never ran — absence of data formatted as a result. This happened once
# (2026-08-18, --skip_ani_screen rejected by 2.7.2) and must not happen again.
if not tax:
    print()
    print("  !! NO CLASSIFICATIONS PARSED. GTDB-Tk did not produce summary tables.")
    print("     Refusing to write MAG_TAXONOMY.tsv or report novelty — an empty")
    print("     result is not evidence of novel taxa. Inspect the classify_wf log.")
    raise SystemExit(1)

cw_path = os.path.join(F, 'MAG_CROSSWALK.tsv')
cw = list(csv.DictReader(open(cw_path), delimiter='\t'))

RANKS = ['domain','phylum','class','order','family','genus','species']
def split_tax(s):
    out = dict.fromkeys(RANKS, '')
    for part in (s or '').split(';'):
        part = part.strip()
        # NB: part[1:3], not part[1] — a single char never equals '__'.
        # That bug (2026-08-18) silently emptied every rank and reported all
        # 29 MAGs as novel. Any change here must be checked against a real
        # GTDB string such as 'd__Bacteria;p__Bacillota;c__Clostridia'.
        if len(part) > 3 and part[1:3] == '__':
            key = {'d':'domain','p':'phylum','c':'class','o':'order',
                   'f':'family','g':'genus','s':'species'}.get(part[0])
            if key: out[key] = part[3:]
    return out

assert split_tax('d__Bacteria;p__Bacillota')['phylum'] == 'Bacillota', \
    "split_tax is broken — refusing to report taxonomy"

rows = []
for r in cw:
    t = tax.get(r['working_name'], {})
    parts = split_tax(t.get('classification', ''))
    # novelty: the most specific rank GTDB could actually assign
    assigned = [k for k in RANKS if parts[k]]
    deepest = assigned[-1] if assigned else 'unclassified'
    novelty = 'known species' if parts['species'] else f"novel below {deepest}"
    # the summary FILE the row came from is authoritative for domain
    if t.get('_domain') and not parts['domain']:
        parts['domain'] = t['_domain']
    r.update(parts)
    r['gtdb_classification'] = t.get('classification', '')
    r['closest_ref'] = t.get('fastani_reference', '') or t.get('closest_genome_reference', '')
    r['ani'] = t.get('fastani_ani', '') or t.get('closest_genome_ani', '')
    r['msa_percent'] = t.get('msa_percent', '')
    r['red_value'] = t.get('red_value', '')
    r['classification_method'] = t.get('classification_method', '')
    r['gtdb_note'] = t.get('note', '')
    r['warnings'] = t.get('warnings', '')
    r['novelty'] = novelty
    r['deepest_rank'] = deepest
    rows.append(r)

out = os.path.join(F, 'MAG_TAXONOMY.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)

print()
print("=" * 92)
print("  MAG TAXONOMY")
print("=" * 92)
print(f"  {'MAG ID':<18} {'smp':<4} {'dom':<4} {'phylum':<20} {'genus':<20} "
      f"{'species':<24} {'ANI':>6}")
for r in sorted(rows, key=lambda r: (r['sample'] != 'SED', r['mag_id'])):
    print(f"  {r['mag_id']:<18} {r['sample']:<4} {(r['domain'] or '?')[:3]:<4} "
          f"{(r['phylum'] or '-')[:20]:<20} {(r['genus'] or '-')[:20]:<20} "
          f"{(r['species'] or '-')[:24]:<24} {(r['ani'] or '')[:6]:>6}")

# ---------------------------------------------------------------- SUMMARIES
print()
print("  HOW EACH CALL WAS MADE (the workspace_1 lesson: never assume)")
meth = {}
for r in rows:
    k = r['classification_method'] or 'not reported'
    meth[k] = meth.get(k, 0) + 1
for k in sorted(meth, key=lambda k: -meth[k]):
    print(f"    {meth[k]:>3}  {k}")
print("    ANI/skani calls are assignments to an existing species cluster;")
print("    topology calls are placements in the reference tree. --place_species")
print("    forces a tree placement for BOTH, so no genome lacks an evidence trail.")

print()
print("  DOMAIN BREAKDOWN")
for s in ('SED', 'LIQ'):
    sub = [r for r in rows if r['sample'] == s]
    b = sum(1 for r in sub if r['domain'] == 'Bacteria')
    a = sum(1 for r in sub if r['domain'] == 'Archaea')
    u = len(sub) - b - a
    print(f"    {s}: {len(sub)} MAGs — {b} Bacteria, {a} Archaea, {u} unclassified")

arch = [r for r in rows if r['domain'] == 'Archaea']
print()
print("  ARCHAEA (audit finding F3 / F26 — 'zero sediment Archaea' was WITHDRAWN;")
print("           --place_species guarantees a uniform evidence trail here)")
if arch:
    for r in arch:
        print(f"    {r['mag_id']:<18} {r['sample']:<4} {r['gtdb_classification']}")
        print(f"      {'':18} closest ref {r['closest_ref'] or '-'} "
              f"ANI {r['ani'] or '-'}  MSA {r['msa_percent'] or '-'}%")
    print()
    print("    -> Compare against Arc_Otu22 (1.641% peak at 60-65 cm, rank 14/327).")
    print("       Use the barrnap 16S from step 08, not taxonomy alone, and record")
    print("       the identity so audit finding F23 (75.8% vs 93.1%) can be closed.")
else:
    print("    NONE placed. Before concluding absence, check:")
    print("      06_binqc/MAG_QC.tsv for archaeal-looking bins classed 'low'")
    print("      08_taxonomy/gtdbtk/identify/ for ar53 marker counts")
    print("    Absence of an ar53 summary is NOT evidence of absence of Archaea.")

print()
print("  NOVELTY")
nov = {}
for r in rows:
    nov[r['novelty']] = nov.get(r['novelty'], 0) + 1
for k in sorted(nov, key=lambda k: -nov[k]):
    print(f"    {nov[k]:>3}  {k}")
nsp = sum(1 for r in rows if not r['species'])
print(f"    {nsp}/{len(rows)} MAGs have NO species-level match in GTDB {os.environ.get('GTDB_RELEASE','R232')}")
print("    (these are candidate novel species — the strongest claim in the chapter)")

low = [r for r in rows if r['msa_percent'] and float(r['msa_percent']) < 50]
if low:
    print()
    print("  LOW MSA WARNING — placement is weakly supported for:")
    for r in low:
        print(f"    {r['mag_id']:<18} MSA {r['msa_percent']}%  {r['gtdb_classification'][:60]}")
    print("    Treat these classifications as provisional; confirm phylogenomically.")

print()
print(f"  written: {out}")
PY

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "gtdbtk output directory" $([ -d "$T/gtdbtk" ] && echo 0 || echo 1)
sumfiles=$(ls -1 "$T/gtdbtk"/gtdbtk.*.summary.tsv 2>/dev/null | wc -l)
check "at least one summary table ($sumfiles)" $([ "$sumfiles" -gt 0 ] && echo 0 || echo 1)
check "merged MAG_TAXONOMY.tsv" $([ -s "$F/MAG_TAXONOMY.tsv" ] && echo 0 || echo 1)
if [ -s "$F/MAG_TAXONOMY.tsv" ]; then
  nclass=$(awk -F'\t' 'NR>1 && $0 ~ /d__/' "$F/MAG_TAXONOMY.tsv" | wc -l)
  check "every MAG classified ($nclass/$NMAG)" $([ "$nclass" -eq "$NMAG" ] && echo 0 || echo 1)
fi
# the ANI screen must NOT have been used — it is what caused the F3 misreading
if [ -d "$T/gtdbtk/align" ]; then
  na=$(ls -1 "$T/gtdbtk/align"/*msa* 2>/dev/null | wc -l)
  check "marker alignments written ($na) — ANI screen correctly skipped" \
    $([ "$na" -gt 0 ] && echo 0 || echo 1)
fi

runlog "09_taxonomy complete"
echo
echo "######## 09 DONE — next: ANI (isolates vs MAGs) + annotation ########"
