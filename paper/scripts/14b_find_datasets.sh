#!/usr/bin/env bash
# =============================================================================
# 14b_find_datasets.sh  --  discover more comparison datasets, reproducibly
#
# WHY A SEARCH TOOL RATHER THAN A LIST
#   Accessions quoted from memory are a reliability hazard — PRJNA404017 turned
#   out to have no SRA records at all, and every "Grevelingen" entry is
#   AMPLICON. This queries SRA directly, filters to SHOTGUN metagenomes
#   (WGS + METAGENOMIC, excluding AMPLICON and RNA), and reports what is
#   actually downloadable. Nothing here is asserted; it is all retrieved.
#
# THE COMPARISON SET SHOULD BE SYSTEMS LIKE ETOLIKO
#   Etoliko is a coastal MEROMICTIC lagoon: permanently stratified, anoxic and
#   sulfidic at depth, sulfate-replete throughout, with no sulfate-methane
#   transition zone. The closest analogues are other meromictic/euxinic basins
#   and anoxic brines, which is why the search terms below target that class
#   rather than "marine sediment" generally.
#
# OUTPUT   14_comparative/candidates.tsv  — paste rows into projects.tsv
#
# RUNTIME  2-5 min (NCBI queries only, no sequence download).
# USAGE    ws2 && bash scripts/14b_find_datasets.sh 2>&1 | tee logs/14b_find.log
#          MIN_SPOTS=5000000 bash scripts/14b_find_datasets.sh
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="14b_find"
C="$WS2/14_comparative"; mkdir -p "$C/meta"
MIN_SPOTS="${MIN_SPOTS:-1000000}"
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

echo "######## 14b — SEARCH SRA FOR COMPARABLE SYSTEMS ########"
echo "  started: $(date -Is)  minimum spots per run: $MIN_SPOTS"
runlog "14b_find_datasets started"

# Search terms, grouped by the property they share with Etoliko.
TERMS=(
  "meromictic AND metagenome"
  "euxinic AND metagenome"
  "anoxic AND sulfidic AND metagenome"
  "sulfidic AND lagoon AND metagenome"
  "coastal lagoon AND sediment AND metagenome"
  "stratified AND lake AND anoxic AND metagenome"
  "deep hypersaline anoxic basin AND metagenome"
  "brine pool AND metagenome"
  "Mediterranean AND anoxic AND sediment AND metagenome"
  "fjord AND anoxic AND metagenome"
  "Cadagno AND metagenome"
  "Ace Lake AND metagenome"
  "Mariager AND metagenome"
  "Framvaren AND metagenome"
  "Lake Faro AND metagenome"
  "blue hole AND metagenome"
  "Guaymas AND sediment AND metagenome"
  "Aarhus Bay AND sediment AND metagenome"
  "White Oak River AND sediment AND metagenome"
  "salt marsh AND sediment AND metagenome"
)

: > "$C/meta/runinfo_all.csv"
for t in "${TERMS[@]}"; do
  enc=$(printf '%s' "$t" | sed 's/ /+/g; s/"/%22/g')
  ids=$(wget -qO- "$EUT/esearch.fcgi?db=sra&term=${enc}&retmax=300" \
        | grep -oE '<Id>[0-9]+</Id>' | grep -oE '[0-9]+' | paste -sd,)
  n=0
  if [ -n "$ids" ]; then
    wget -qO- "$EUT/efetch.fcgi?db=sra&id=${ids}&rettype=runinfo&retmode=text" \
      >> "$C/meta/runinfo_all.csv"
    n=$(echo "$ids" | tr ',' '\n' | wc -l)
  fi
  printf '  %-58s %4s records\n' "$t" "$n"
  sleep 1
done

echo
echo "######## FILTERING TO SHOTGUN METAGENOMES ########"
python3 - "$C" "$MIN_SPOTS" <<'PY'
import sys, os, csv, collections
C, MIN = sys.argv[1], int(sys.argv[2])
path = os.path.join(C, 'meta', 'runinfo_all.csv')

rows, hdr = [], None
with open(path, newline='', errors='replace') as fh:
    for line in fh:
        f = next(csv.reader([line]))
        if f and f[0] == 'Run':
            hdr = f; continue
        if hdr and f and f[0].startswith(('SRR', 'ERR', 'DRR')) and len(f) >= len(hdr):
            rows.append(dict(zip(hdr, f)))

def g(r, k, d=''):
    return (r.get(k) or d).strip()

seen, keep = set(), []
for r in rows:
    run = g(r, 'Run')
    if run in seen: continue
    seen.add(run)
    if g(r, 'LibraryStrategy').upper() != 'WGS': continue
    if g(r, 'LibrarySource').upper() != 'METAGENOMIC': continue
    if 'RNA' in g(r, 'LibrarySelection').upper(): continue
    try:
        if int(g(r, 'spots', '0') or 0) < MIN: continue
    except ValueError:
        continue
    keep.append(r)

print(f"  {len(rows)} run records -> {len(keep)} shotgun metagenome runs "
      f"(>= {MIN:,} spots)")

# already-used projects, so we only report NEW ones
used = set()
p = os.path.join(C, 'meta', 'projects.tsv')
if os.path.exists(p):
    used = {l.split('\t')[0] for l in open(p) if l.strip()}

by_prj = collections.defaultdict(list)
for r in keep:
    by_prj[g(r, 'BioProject')].append(r)

out = os.path.join(C, 'candidates.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['bioproject', 'n_runs', 'total_Gbp', 'platform', 'organism',
                'example_run', 'title'])
    print()
    print(f"  {'BioProject':<14}{'runs':>5}{'Gbp':>8}  {'platform':<12}{'organism / title'}")
    ranked = sorted(by_prj.items(), key=lambda kv: -len(kv[1]))
    for prj, rs in ranked:
        if not prj or prj in used: continue
        gbp = sum(int(g(r, 'bases', '0') or 0) for r in rs) / 1e9
        plat = collections.Counter(g(r, 'Platform') for r in rs).most_common(1)[0][0]
        org = collections.Counter(g(r, 'ScientificName') for r in rs).most_common(1)[0][0]
        title = g(rs[0], 'SampleName') or g(rs[0], 'Experiment')
        w.writerow([prj, len(rs), f"{gbp:.1f}", plat, org, g(rs[0], 'Run'), title])
        if len(ranked) < 200:
            print(f"  {prj:<14}{len(rs):>5}{gbp:>8.1f}  {plat[:12]:<12}{(org+' | '+title)[:58]}")

print()
print(f"  written: {out}")
print()
print("  TO USE: pick the projects whose system genuinely resembles Etoliko —")
print("  permanently stratified, anoxic, sulfidic, sulfate-replete — and append")
print("  them to 14_comparative/meta/projects.tsv as")
print("      <BioProject><TAB><system description><TAB>DNA")
print("  then re-run 14_comparative.sh. Cached datasets are skipped, so only the")
print("  new ones are downloaded.")
print()
print("  Judge by the SYSTEM, not the run count. Ten runs from an oxic estuary")
print("  are worth less than two from a meromictic sulfidic lake.")
PY

echo
echo "######## VALIDATION GATE ########"
check "candidate table written" $([ -s "$C/candidates.tsv" ] && echo 0 || echo 1)
n=$(( $(wc -l < "$C/candidates.tsv") - 1 ))
check "candidate projects found ($n)" $([ "$n" -gt 0 ] && echo 0 || echo 1)

runlog "14b_find_datasets complete ($n candidates)"
echo
echo "######## 14b DONE ########"
