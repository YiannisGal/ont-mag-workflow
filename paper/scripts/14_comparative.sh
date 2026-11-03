#!/usr/bin/env bash
# =============================================================================
# 14_comparative.sh  --  are these organisms Etoliko-endemic or euxinic-basin
#                        generalists? Fragment recruitment against public data.
#
# THE QUESTION
#   29 MAGs, 22 of them candidate novel species. Novel to GTDB is not the same
#   as novel to the world: an organism can be absent from reference databases
#   and still be common in every anoxic sulfidic basin on the planet. Nobody can
#   tell which from a GTDB classification alone.
#
#   So: take public metagenomes from comparable systems and ask whether reads
#   from those sites recruit to YOUR genomes. A MAG that recruits in Cariaco,
#   the Black Sea and Orca Basin is a euxinic-basin generalist that simply had
#   no genome until now. One that recruits nowhere is Etoliko-specific. Both
#   answers are publishable; the second is the stronger novelty claim.
#
# WHY FRAGMENT RECRUITMENT AND NOT FULL DOWNLOADS
#   Full download of eleven projects is terabytes. Recruitment needs only enough
#   reads to detect presence, so ~2M reads per run are STREAMED with
#   fastq-dump -X, which stops after the requested spot count instead of
#   fetching the whole run. ~5 min per run at the measured 1.15 MB/s.
#
# THE METRIC THAT MATTERS: BREADTH, NOT DEPTH
#   Conserved genes (ribosomal proteins, chaperones) recruit reads from distant
#   relatives and inflate mean depth, so depth alone manufactures false
#   presence. Covered FRACTION is the honest test: a genome truly present is
#   covered evenly across its length. The reporting threshold here is
#   >=30% breadth, and both numbers are always shown so the reader can judge.
#
# YOUR OWN SAMPLES ARE INCLUDED as positive controls. If Etoliko reads do not
#   recruit strongly to Etoliko MAGs, the pipeline is wrong and no comparison
#   result should be believed.
#
# OUTPUT   14_comparative/recruitment_matrix.tsv   MAG x dataset breadth+depth
#          14_comparative/summary.tsv              per-MAG biogeography call
#
# RUNTIME  2-4 h. RUN IN tmux.
# USAGE    ws2 && bash scripts/14_comparative.sh 2>&1 | tee logs/14_comparative.log
#          RUNS_PER_PROJECT=1 SPOTS=1000000 bash scripts/14_comparative.sh   # faster
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="14_comparative"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
C="$WS2/14_comparative"; F="$WS2/07_mags"
mkdir -p "$C/reads" "$C/bam" "$C/cov" "$C/meta"

SPOTS="${SPOTS:-2000000}"
NRUN="${RUNS_PER_PROJECT:-2}"
BREADTH_MIN="${BREADTH_MIN:-0.30}"

echo "######## 14 — CROSS-STUDY FRAGMENT RECRUITMENT ########"
echo "  started: $(date -Is)  spots/run: $SPOTS  runs/project: $NRUN"
runlog "14_comparative started (spots=$SPOTS runs=$NRUN)"

# --------------------------------------------------------------- DATASETS
# project<TAB>system<TAB>type   — edit freely; DNA only is used for recruitment
#
# WRITTEN ONLY IF ABSENT. The first version of this script used `>` and
# overwrote the file on every run, silently destroying datasets the user had
# appended (observed 2026-08-19: seven added meromictic-lake projects vanished
# and the run repeated the original nine). Delete the file to restore defaults.
if [ ! -s "$C/meta/projects.tsv" ]; then
cat > "$C/meta/projects.tsv" <<'PRJ'
PRJEB34883	Baltic Sea coastal transect + LMO time series	DNA
PRJNA326482	Cariaco Basin (permanently anoxic, sulfidic)	DNA
PRJNA638805	Black Sea 5/30/150/750 m (euxinic below ~150 m)	DNA
PRJNA247822	Saanich Inlet (seasonally anoxic fjord)	DNA
PRJNA1426066	Orca Basin (hypersaline anoxic brine, Gulf of Mexico)	DNA
PRJNA349072	Santa Barbara Basin (oxygen-deficient, sulfidic sediment)	DNA
PRJNA299097	Red Sea Atlantis II / Discovery / Chain brine pools	DNA
PRJNA404017	Lake Grevelingen saline sediment	DNA
PRJEB26778	Rogoznica Lake chemocline + hypolimnion	DNA
PRJNA428921	Rogoznica / Dragon's Eye salt lake	DNA
PRJNA1413184	Orca Basin metatranscriptome	RNA
PRJ
  echo "  wrote default project list"
else
  echo "  using existing project list (delete $C/meta/projects.tsv to reset)"
fi
# tolerate space-separated rows pasted by hand
if grep -qP '^\S+ +\S' "$C/meta/projects.tsv" 2>/dev/null; then
  awk '{ if ($0 !~ /\t/) { p=$1; $1=""; sub(/^ +/,""); t=$NF; sub(/ +[^ ]+$/,"");
         printf "%s\t%s\t%s\n", p, $0, t } else print }' \
      "$C/meta/projects.tsv" > "$C/meta/projects.tsv.fix" \
    && mv "$C/meta/projects.tsv.fix" "$C/meta/projects.tsv"
  echo "  normalised space-separated rows to tab-separated"
fi
echo "  projects listed: $(wc -l < "$C/meta/projects.tsv")"

# ------------------------------------------------- resolve project -> runs
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
resolve_runs() {  # resolve_runs <bioproject>
  local prj="$1" out="$C/meta/${prj}_runs.txt"
  if [ -s "$out" ]; then cat "$out"; return 0; fi
  local ids
  ids=$(wget -qO- "$EUT/esearch.fcgi?db=sra&term=${prj}&retmax=500" \
        | grep -oE '<Id>[0-9]+</Id>' | grep -oE '[0-9]+' | head -200 | paste -sd,)
  [ -n "$ids" ] || { echo "" > "$out"; return 1; }
  wget -qO- "$EUT/efetch.fcgi?db=sra&id=${ids}&rettype=runinfo&retmode=text" \
    | awk -F',' 'NR>1 && $1 ~ /^[SED]RR/ {print $1"\t"$4"\t"$8}' \
    | sort -k2 -rn > "$out"          # run, spots, bases  (largest first)
  cat "$out"
}

# ---------------------------------------------------------------- MAG REF
REF="$C/all_mags.fna"
if [ ! -s "$REF" ]; then
  for f in "$F/bins"/*.fa; do
    n=$(basename "$f" .fa)
    awk -v m="$n" '/^>/{print ">"m"~"substr($0,2); next}{print}' "$f"
  done > "$REF"
fi
NMAG=$(ls -1 "$F/bins"/*.fa | wc -l)
echo "  reference: $NMAG MAGs, $(grep -c '^>' "$REF") contigs"

MMI="$C/all_mags.mmi"
[ -s "$MMI" ] || run_tool minimap2 -d "$MMI" "$REF" 2>/dev/null

# ------------------------------------------------------------- RECRUIT
recruit() {  # recruit <label> <fastq(.gz)> [preset]
  local lab="$1" fq="$2" preset="${3:-}"
  local bam="$C/bam/${lab}.bam" cov="$C/cov/${lab}.tsv"
  [ -s "$cov" ] && { echo "    $lab: cached"; return 0; }
  [ -s "$fq" ] || { echo "    $lab: no reads"; return 1; }

  if [ -z "$preset" ]; then
    local mlen
    mlen=$(head -400000 <(zcat -f "$fq") | awk 'NR%4==2{s+=length($0); n++} END{print (n? int(s/n):0)}')
    if [ "${mlen:-0}" -gt 500 ]; then preset="map-ont"; else preset="sr"; fi
    echo "    $lab: mean read ${mlen} bp -> minimap2 -x $preset"
  fi

  run_tool minimap2 -ax "$preset" -t "$THREADS" --secondary=no "$MMI" "$fq" 2>/dev/null \
    | run_tool samtools sort -@ 4 -m 1G -o "$bam" - 2>/dev/null
  run_tool samtools index "$bam" 2>/dev/null
  run_tool coverm genome --bam-files "$bam" --separator '~' \
      -m mean covered_fraction relative_abundance \
      --min-covered-fraction 0 -t "$THREADS" -o "$cov" 2>/dev/null \
    || echo "    !! coverm failed for $lab"
  rm -f "$bam" "$bam.bai"          # BAMs are large and not needed downstream
}

# ---- positive controls: our own reads --------------------------------------
echo
echo "######## POSITIVE CONTROLS — Etoliko reads vs Etoliko MAGs ########"
# RANDOM subsampling, not head. The sediment fastq is TEN concatenated
# sequencing runs, so the first N reads come from whichever run is first in the
# file rather than from the community. Measured 2026-08-19: head-based
# subsample recruited 0.04% of reads to the MAGs; a random subsample of the
# same size recruited 4.64% — a 116x difference, and the reason the positive
# control failed. Never use head to subsample a concatenated read set.
for pair in "ETOLIKO_SED:$SED_READS" "ETOLIKO_LIQ:$LIQ_READS"; do
  lab="${pair%%:*}"; fq="${pair#*:}"
  sub="$C/reads/${lab}.fastq.gz"
  if [ ! -s "$sub" ]; then
    echo "  random-subsampling $lab to ~$SPOTS reads ..."
    total=$(zcat -f "$fq" | awk 'END{print int(NR/4)}')
    frac=$(awk -v s="$SPOTS" -v t="$total" 'BEGIN{f=(t>0? s/t : 1); print (f>1?1:f)}')
    echo "    $total reads in file -> sampling fraction $frac"
    run_tool seqkit sample -p "$frac" -s 42 "$fq" -o "$sub" 2>/dev/null \
      || { echo "    !! seqkit sample failed; falling back to head (UNREPRESENTATIVE)"
           zcat -f "$fq" | head -n $((SPOTS*4)) | gzip > "$sub"; }
  fi
  recruit "$lab" "$sub" map-ont
done

# ---- public datasets --------------------------------------------------------
echo
echo "######## PUBLIC DATASETS ($(date +%H:%M)) ########"
while IFS=$'\t' read -r prj system dtype; do
  [ -z "$prj" ] && continue
  echo
  echo "--- $prj  $system  [$dtype]"
  if [ "$dtype" = "RNA" ]; then
    echo "    metatranscriptome — excluded from DNA recruitment (different units;"
    echo "    RNA recruitment measures expression, not presence). Skipping."
    continue
  fi
  runs=$(resolve_runs "$prj" | awk 'NF' | head -"$NRUN" | cut -f1)
  if [ -z "$runs" ]; then echo "    !! could not resolve runs for $prj"; continue; fi
  for run in $runs; do
    lab="${prj}_${run}"
    fq="$C/reads/${run}.fastq.gz"
    if [ ! -s "$fq" ] && [ ! -s "$C/cov/${lab}.tsv" ]; then
      echo "    streaming $run ($SPOTS spots) ..."
      ( cd "$C/reads" && fastq-dump --gzip -X "$SPOTS" --skip-technical \
          --readids --split-spot -O . "$run" > /dev/null 2>&1 ) \
        || echo "      !! fastq-dump failed for $run"
      [ -s "$C/reads/${run}.fastq.gz" ] || \
        mv "$C/reads/${run}_1.fastq.gz" "$fq" 2>/dev/null || true
    fi
    recruit "$lab" "$fq"
    rm -f "$C/reads/${run}"*.fastq.gz     # keep disk bounded
  done
done < "$C/meta/projects.tsv"

# =====================================================================
# MATRIX + BIOGEOGRAPHY CALLS
# =====================================================================
python3 - "$C" "$F" "$BREADTH_MIN" <<'PY'
import sys, os, csv, glob, collections
C, F, BMIN = sys.argv[1], sys.argv[2], float(sys.argv[3])

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

sysmap = {}
p = os.path.join(C, 'meta', 'projects.tsv')
if os.path.exists(p):
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2: sysmap[f[0]] = f[1]

data = collections.defaultdict(dict)      # mag -> dataset -> (breadth, depth)
datasets = []
for cov in sorted(glob.glob(os.path.join(C, 'cov', '*.tsv'))):
    lab = os.path.basename(cov)[:-4]
    rows = list(csv.reader(open(cov), delimiter='\t'))
    if len(rows) < 2: continue
    hdr, body = rows[0], rows[1:]
    mi = [i for i, h in enumerate(hdr) if h.endswith('Mean')]
    bi = [i for i, h in enumerate(hdr) if 'Covered Fraction' in h]
    if not mi or not bi: continue
    datasets.append(lab)
    for r in body:
        if not r or r[0].lower().startswith('unmapped'): continue
        try: data[r[0]][lab] = (float(r[bi[0]]), float(r[mi[0]]))
        except (ValueError, IndexError): pass

if not datasets:
    print("  !! no coverage tables — nothing recruited"); raise SystemExit(1)

pub = [d for d in datasets if not d.startswith('ETOLIKO')]
ctl = [d for d in datasets if d.startswith('ETOLIKO')]

with open(os.path.join(C, 'recruitment_matrix.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['mag_id','working_name','domain','genus']
               + [f"{d}_breadth" for d in datasets] + [f"{d}_depth" for d in datasets])
    for m in sorted(data):
        t = tax.get(m, {})
        w.writerow([t.get('mag_id',''), m, t.get('domain',''), t.get('genus','')]
                   + [f"{data[m].get(d,(0,0))[0]:.3f}" for d in datasets]
                   + [f"{data[m].get(d,(0,0))[1]:.2f}" for d in datasets])

# ---------------------------------------------------- controls first
print()
print("=" * 96)
print("  POSITIVE CONTROLS — Etoliko reads must recruit to Etoliko MAGs")
print("=" * 96)
ok = True
for d in ctl:
    hi = sum(1 for m in data if data[m].get(d, (0, 0))[0] >= BMIN)
    print(f"  {d:<16} {hi}/{len(data)} MAGs at >={BMIN:.0%} breadth")
    if hi == 0: ok = False
if not ok:
    print("  !! CONTROL FAILED — do not interpret the comparison below.")
else:
    print("  Controls pass; the recruitment pipeline is behaving.")

# ---------------------------------------------------- biogeography
print()
print("=" * 96)
print(f"  BIOGEOGRAPHY — public datasets recruiting at >={BMIN:.0%} genome breadth")
print("=" * 96)
print(f"  {'MAG ID':<17}{'genus':<20}{'n sites':>8}   sites")
rows = []
for m in sorted(data, key=lambda m: tax.get(m, {}).get('mag_id', m)):
    t = tax.get(m, {})
    hitsites = [d for d in pub if data[m].get(d, (0, 0))[0] >= BMIN]
    systems = sorted({sysmap.get(d.split('_')[0], d.split('_')[0]) for d in hitsites})
    call = ('cosmopolitan' if len(systems) >= 3 else
            'shared' if systems else 'Etoliko-only')
    rows.append(dict(mag_id=t.get('mag_id',''), working_name=m,
                     domain=t.get('domain',''), genus=t.get('genus',''),
                     n_sites=len(hitsites), n_systems=len(systems),
                     call=call, systems='; '.join(systems)))
    print(f"  {t.get('mag_id',m):<17}{(t.get('genus','-') or '-')[:20]:<20}"
          f"{len(hitsites):>8}   {'; '.join(s[:34] for s in systems[:3]) or '—'}")

with open(os.path.join(C, 'summary.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)

n = collections.Counter(r['call'] for r in rows)
print()
print(f"  cosmopolitan (>=3 systems) : {n['cosmopolitan']}")
print(f"  shared (1-2 systems)       : {n['shared']}")
print(f"  Etoliko-only               : {n['Etoliko-only']}")
print()
print("  'Etoliko-only' means NOT DETECTED in the datasets tested at this depth")
print("  of sampling. It is evidence of restricted distribution, not proof of")
print("  endemism — 2M reads per run cannot detect a genome below ~0.01%")
print("  relative abundance. Say exactly that in the text.")
PY

echo
echo "######## VALIDATION GATE ########"
check "recruitment matrix" $([ -s "$C/recruitment_matrix.tsv" ] && echo 0 || echo 1)
check "summary table"      $([ -s "$C/summary.tsv" ] && echo 0 || echo 1)
nd=$(ls -1 "$C/cov"/*.tsv 2>/dev/null | wc -l)
check "datasets recruited ($nd)" $([ "$nd" -gt 2 ] && echo 0 || echo 1)
nc=$(ls -1 "$C/cov"/ETOLIKO_*.tsv 2>/dev/null | wc -l)
check "positive controls present ($nc/2)" $([ "$nc" -eq 2 ] && echo 0 || echo 1)

runlog "14_comparative complete"
echo
echo "######## 14 DONE ########"
