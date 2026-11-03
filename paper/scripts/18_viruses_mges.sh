#!/usr/bin/env bash
# =============================================================================
# 18_viruses_mges.sh  --  viruses, proviruses and plasmids
#
# THE GAP
#   Step 16 counted phage-like GENES inside MAGs. That is not virus discovery:
#   it cannot find a free viral contig, cannot tell a prophage from a decayed
#   remnant, and cannot link a virus to a host. Viruses are unexamined in this
#   system, and in anoxic sulfidic sediments they typically outnumber cells and
#   are a dominant control on mortality and gene transfer.
#
# WHAT THIS DOES
#   geNomad on BOTH levels, which answers different questions:
#     A. WHOLE ASSEMBLY  -> free viral contigs and plasmids in the community,
#        including the ~98% of sediment that never binned. Quantified with the
#        step-05 coverage tables, so viral abundance is measured, not counted.
#     B. MAG CONTIGS     -> proviruses integrated in the 29 genomes. A provirus
#        found here comes with its host attached, which is the cleanest
#        virus-host link obtainable from a metagenome.
#
#   geNomad rather than VirSorter2/VIBRANT: one tool for viruses AND plasmids,
#   ~3.7 GB database rather than 11 GB, and it reports taxonomy and marker
#   evidence per contig instead of a bare score.
#
# INTERPRETING SCORES — read this before quoting any number
#   geNomad's virus_score is a classifier probability, not proof. The default
#   cutoff is deliberately permissive. This script reports counts at the
#   default AND at a conservative threshold (>=0.8 with >=1 viral hallmark
#   gene), because the difference between those two numbers is the honest
#   uncertainty in "how many viruses are there".
#
# SETUP REQUIRED ONCE (geNomad is not installed on this system):
#   conda create -y -n genomad_env -c conda-forge -c bioconda genomad
#   conda run -n genomad_env genomad download-database ${DB_DIR}
#   echo 'export ENV_genomad="genomad_env"' >> env_resolved.sh
#   echo 'export GENOMAD_DB="${DB_DIR}/genomad_db"' >> 00_config.sh
#   (database ~3.7 GB; at the measured 1.15 MB/s from NCBI expect under an hour)
#
# INPUT    04_polish/{SED,LIQ}_asm_final.fasta, 07_mags/bins/, 02_mapping/*_depth_multi.txt
# OUTPUT   18_viruses/
#
# RUNTIME  1-3 h for the assemblies, minutes for the MAGs. Use tmux.
# USAGE    ws2 && bash scripts/18_viruses_mges.sh 2>&1 | tee logs/18_viruses.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="18_viruses"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
V="$WS2/18_viruses"; F="$WS2/07_mags"; mkdir -p "$V"

echo "######## 18 — VIRUSES, PROVIRUSES AND PLASMIDS ########"
echo "  started: $(date -Is)  threads: $THREADS"
runlog "18_viruses_mges started"

# --------------------------------------------------------------- PREFLIGHT
if ! run_tool genomad --version >/dev/null 2>&1; then
  echo
  echo "  !! geNomad is not callable. Install it once (user-owned env, nothing shared):"
  echo
  echo "     conda create -y -n genomad_env -c conda-forge -c bioconda genomad"
  echo "     conda run -n genomad_env genomad download-database ${DB_DIR}"
  echo "     cd $WS2"
  echo "     echo 'export ENV_genomad=\"genomad_env\"' >> env_resolved.sh"
  echo "     echo 'export GENOMAD_DB=\"${DB_DIR}/genomad_db\"' >> 00_config.sh"
  echo "     source 00_config.sh && run_tool genomad --version"
  echo
  echo "  Then re-run this script."
  exit 1
fi
record_version "genomad" "${ENV_genomad:-genomad_env}" genomad --version
[ -d "${GENOMAD_DB:-/nonexistent}" ] || { echo "  !! GENOMAD_DB not set or missing: ${GENOMAD_DB:-UNSET}"; exit 1; }
echo "  [db] $GENOMAD_DB"

# =====================================================================
# A. WHOLE-ASSEMBLY SURVEY
# =====================================================================
echo
echo "######## A. FREE VIRUSES AND PLASMIDS IN THE ASSEMBLIES ($(date +%H:%M)) ########"
for s in SED LIQ; do
  asm="$WS2/04_polish/${s}_asm_final.fasta"
  [ -s "$asm" ] || { echo "  $s: assembly missing"; continue; }
  out="$V/${s}_assembly"
  if [ -d "$out/${s}_asm_final_summary" ]; then
    echo "  $s: already done"; continue
  fi
  echo "  $s: running geNomad end-to-end (this is the long part) ..."
  run_tool genomad end-to-end --cleanup --threads "$THREADS" \
      "$asm" "$out" "$GENOMAD_DB" > "$V/${s}_genomad.log" 2>&1 \
    || echo "    !! geNomad failed on $s — see $V/${s}_genomad.log"
done

# =====================================================================
# B. PROVIRUSES INSIDE THE MAGs — virus-host links with the host known
# =====================================================================
echo
echo "######## B. PROVIRUSES IN THE 29 MAGs ($(date +%H:%M)) ########"
MAGCAT="$V/all_mags_contigs.fna"
if [ ! -s "$MAGCAT" ]; then
  for f in "$F/bins"/*.fa; do
    n=$(basename "$f" .fa)
    awk -v m="$n" '/^>/{sub(/^>/,""); split($0,a," "); print ">"m"~"a[1]; next}{print}' "$f"
  done > "$MAGCAT"
fi
echo "  MAG contigs: $(grep -c '^>' "$MAGCAT")"
if [ ! -d "$V/mags_genomad/all_mags_contigs_summary" ]; then
  run_tool genomad end-to-end --cleanup --threads "$THREADS" \
      "$MAGCAT" "$V/mags_genomad" "$GENOMAD_DB" > "$V/mags_genomad.log" 2>&1 \
    || echo "  !! geNomad failed on MAG contigs — see $V/mags_genomad.log"
fi

# =====================================================================
# C. REPORT
# =====================================================================
python3 - "$V" "$WS2" "$F" <<'PY'
import sys, os, csv, glob
from collections import defaultdict, Counter

V, WS2, F = sys.argv[1], sys.argv[2], sys.argv[3]
STRICT_SCORE, STRICT_HALLMARK = 0.80, 1

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

def find(pat):
    hits = glob.glob(pat, recursive=True)
    return hits[0] if hits else None

def load(path):
    if not path or not os.path.exists(path): return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def num(r, k, d=0.0):
    try: return float(r.get(k, d) or d)
    except ValueError: return d

# ---------------- assemblies
print()
print("=" * 100)
print("  A. VIRUSES AND PLASMIDS IN THE ASSEMBLIES")
print("=" * 100)
summary_rows = []
for s in ('SED', 'LIQ'):
    vf = find(os.path.join(V, f'{s}_assembly', '**', '*_virus_summary.tsv'))
    pf = find(os.path.join(V, f'{s}_assembly', '**', '*_plasmid_summary.tsv'))
    vir, pla = load(vf), load(pf)
    if not vir and not pla:
        print(f"\n  {s}: no geNomad summary found"); continue

    strict = [r for r in vir
              if num(r, 'virus_score') >= STRICT_SCORE
              and num(r, 'n_hallmarks') >= STRICT_HALLMARK]
    # coverage-weighted abundance from step 05
    depth = {}
    dfile = os.path.join(WS2, '02_mapping', f'{s}_depth_multi.txt')
    if os.path.exists(dfile):
        with open(dfile) as fh:
            next(fh, None)
            for line in fh:
                f_ = line.split('\t')
                if len(f_) >= 3:
                    try: depth[f_[0]] = float(f_[1]) * float(f_[2])
                    except ValueError: pass
    tot_bases = sum(depth.values()) or 1
    vb = sum(depth.get(r.get('seq_name', ''), 0) for r in strict)
    pb = sum(depth.get(r.get('seq_name', ''), 0) for r in pla)

    print(f"\n  {s}")
    print(f"    viral contigs   : {len(vir):>6} at geNomad default | "
          f"{len(strict):>6} at score>={STRICT_SCORE} with >={STRICT_HALLMARK} hallmark")
    print(f"    plasmid contigs : {len(pla):>6}")
    print(f"    viral share of mapped bases   : {100*vb/tot_bases:.2f}%")
    print(f"    plasmid share of mapped bases : {100*pb/tot_bases:.2f}%")
    if strict:
        lens = sorted(num(r, 'length') for r in strict)
        print(f"    viral contig length: median {lens[len(lens)//2]:,.0f} bp, "
              f"max {lens[-1]:,.0f} bp")
        tx = Counter()
        for r in strict:
            t = (r.get('taxonomy') or '').split(';')
            tx[t[-1] if len(t) > 1 else 'unclassified'] += 1
        print(f"    top viral taxa: " +
              ', '.join(f'{k} ({v})' for k, v in tx.most_common(5)))
    summary_rows.append(dict(sample=s, n_virus_default=len(vir),
                             n_virus_strict=len(strict), n_plasmid=len(pla),
                             pct_bases_viral=f"{100*vb/tot_bases:.3f}",
                             pct_bases_plasmid=f"{100*pb/tot_bases:.3f}"))

if summary_rows:
    with open(os.path.join(V, 'assembly_virus_summary.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(summary_rows)

# ---------------- proviruses in MAGs
print()
print("=" * 100)
print("  B. PROVIRUSES INSIDE THE MAGs — virus-host links, host known")
print("=" * 100)
pv = find(os.path.join(V, 'mags_genomad', '**', '*_virus_summary.tsv'))
rows = load(pv)
if not rows:
    print("  no geNomad output for the MAG contigs")
else:
    by_mag = defaultdict(list)
    for r in rows:
        seq = r.get('seq_name', '')
        mag = seq.split('~')[0] if '~' in seq else 'unknown'
        by_mag[mag].append(r)
    print(f"  {'MAG ID':<17}{'genus':<20}{'regions':>8}{'total kb':>10}   top taxonomy")
    out = []
    for m in sorted(by_mag, key=lambda m: -len(by_mag[m])):
        t = tax.get(m, {})
        rs = by_mag[m]
        kb = sum(num(r, 'length') for r in rs) / 1000
        tx = Counter()
        for r in rs:
            parts = (r.get('taxonomy') or '').split(';')
            tx[parts[-1] if len(parts) > 1 else 'unclassified'] += 1
        top = tx.most_common(1)[0][0] if tx else '-'
        print(f"  {t.get('mag_id', m)[:17]:<17}{(t.get('genus','-') or '-')[:20]:<20}"
              f"{len(rs):>8}{kb:>10.1f}   {top}")
        for r in rs:
            out.append(dict(mag_id=t.get('mag_id', m), working_name=m,
                            genus=t.get('genus', ''), domain=t.get('domain', ''),
                            seq_name=r.get('seq_name', ''),
                            length=r.get('length', ''),
                            virus_score=r.get('virus_score', ''),
                            n_hallmarks=r.get('n_hallmarks', ''),
                            taxonomy=r.get('taxonomy', '')))
    if out:
        with open(os.path.join(V, 'provirus_host_links.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter='\t')
            w.writeheader(); w.writerows(out)
    n_inf = len(by_mag)
    print()
    print(f"  {n_inf}/29 MAGs carry at least one predicted provirus.")
    print("  Cross-check against step 16: MAGs with CRISPR arrays AND a provirus")
    print("  have both a record of past infection and a current resident — the")
    print("  strongest single-genome evidence of ongoing virus-host interaction.")

print()
print("  CAVEATS")
print("   - geNomad's default cutoff is permissive. Both counts are reported; the")
print("     gap between them IS the uncertainty. Quote the strict number.")
print("   - A predicted provirus may be a decayed remnant rather than an inducible")
print("     phage. Completeness (CheckV) would separate those; not run here.")
print("   - Viral abundance is coverage-weighted from step 05, so it reflects DNA")
print("     in the extract. It is not a virus:cell ratio — that needs counts.")
PY

echo
echo "######## VALIDATION GATE ########"
check "assembly summary"     $([ -s "$V/assembly_virus_summary.tsv" ] && echo 0 || echo 1)
check "provirus host links"  $([ -s "$V/provirus_host_links.tsv" ] && echo 0 || echo 1)
for s in SED LIQ; do
  check "$s geNomad output" $(ls -d "$V/${s}_assembly"/*_summary >/dev/null 2>&1 && echo 0 || echo 1)
done

runlog "18_viruses_mges complete"
echo
echo "######## 18 DONE ########"
echo "The provirus-host table is the novel part —"
echo "no viral analysis exists for this lagoon at all."
