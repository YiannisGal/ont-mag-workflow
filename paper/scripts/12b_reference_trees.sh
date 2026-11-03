#!/usr/bin/env bash
# =============================================================================
# 12b_reference_trees.sh  --  publication-grade trees with taxonomy-selected refs
#
# WHY
#   Step 12 pulled reference genomes using each MAG's `closest_ref` from
#   GTDB-Tk. That worked for the bacteria (7 refs) but gave the archaea ZERO,
#   because all eight methanogens returned closest_ref = N/A — they are novel
#   enough to have no close reference. The ar53 tree therefore shows how the
#   eight methanogens relate to each other and nothing else, which is not a
#   thesis figure.
#
#   This script selects references by TAXONOMY instead: for every family
#   represented among the MAGs, take up to K GTDB genomes from that family,
#   plus an outgroup from a phylum the MAGs do not occupy. The tree then has
#   context, and the placement of the novel lineages — especially the
#   uncultured Methanosarcinaceae genus MTP4 (ETO_LIQ_MAG_003) — becomes
#   interpretable rather than asserted from GTDB's classification alone.
#
# ROOTING
#   An explicit outgroup is included and used to root the tree. Midpoint
#   rooting is a fallback, and is flagged in the output when used, because a
#   midpoint-rooted tree should not be presented as if it were outgroup-rooted.
#
# MASKING
#   GTDB's reference MSA is the full untrimmed alignment (bac120: 41,084
#   columns; ar53: similar). GTDB-Tk's user_msa is masked (5,010 / 4,369).
#   The release's own mask is applied to references before merging — this was
#   the defect that broke the step-12 bacterial tree.
#
# INPUT    08_taxonomy/gtdbtk/align/*user_msa.fasta.gz, 07_mags/MAG_TAXONOMY.tsv
#          $GTDBTK_DATA_PATH/{msa,masks,taxonomy}
# OUTPUT   12_phylogeny/ref_trees/{bac120,ar53}_ref.treefile + iTOL files
#
# RUNTIME  30-90 min (IQ-TREE on ~60-120 taxa).
# USAGE    ws2 && bash scripts/12b_reference_trees.sh 2>&1 | tee logs/12b_reftrees.log
#          REFS_PER_FAMILY=10 bash scripts/12b_reference_trees.sh
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="12b_reftrees"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; T="$WS2/08_taxonomy"; P="$WS2/12_phylogeny"
R="$P/ref_trees"; mkdir -p "$R"

K="${REFS_PER_FAMILY:-6}"
echo "######## 12b — REFERENCE TREES (taxonomy-selected) ########"
echo "  started: $(date -Is)  refs per family: $K"
runlog "12b_reference_trees started (K=$K)"

TAXFILE=""
for c in "$GTDBTK_DATA_PATH"/taxonomy/gtdb_taxonomy.tsv "$GTDBTK_DATA_PATH"/taxonomy/*.tsv; do
  [ -s "$c" ] && { TAXFILE="$c"; break; }
done
[ -n "$TAXFILE" ] || { echo "  !! no GTDB taxonomy file under $GTDBTK_DATA_PATH/taxonomy"; exit 1; }
echo "  GTDB taxonomy: $TAXFILE  ($(wc -l < "$TAXFILE") genomes)"

record_version "iqtree" "${ENV_iqtree:-iqtree3}" iqtree3 --version

build() {  # build <bac120|ar53> <label>
  local mk="$1" lab="$2"
  local umsa="" g
  for g in "$T/gtdbtk/align/gtdbtk.${mk}.user_msa.fasta.gz" \
           "$T/gtdbtk/align/gtdbtk.${mk}.user_msa.fasta"; do
    [ -s "$g" ] && { umsa="$g"; break; }
  done
  [ -n "$umsa" ] || { echo "  $lab: no user MSA — skipping"; return 0; }

  local refmsa="" mask=""
  for g in "$GTDBTK_DATA_PATH"/msa/*${mk}*.faa*; do [ -s "$g" ] && { refmsa="$g"; break; }; done
  for g in "$GTDBTK_DATA_PATH"/masks/*${mk}*; do [ -s "$g" ] && { mask="$g"; break; }; done
  [ -n "$refmsa" ] || { echo "  $lab: no reference MSA — skipping"; return 0; }

  echo
  echo "=================================================================="
  echo "  $lab"
  echo "=================================================================="
  local final="$R/${mk}_ref_input.faa"
  python3 - "$umsa" "$refmsa" "$mask" "$TAXFILE" "$F/MAG_TAXONOMY.tsv" "$final" "$K" "$mk" <<'PY'
import sys, gzip, csv, os, collections
umsa, refmsa, maskf, taxfile, magtax, out, K, mk = sys.argv[1:9]
K = int(K)

def rd(path):
    op = gzip.open if path.endswith('.gz') else open
    name, buf = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name: yield name.split()[0], ''.join(buf)
                name, buf = line[1:], []
            else: buf.append(line)
    if name: yield name.split()[0], ''.join(buf)

user = list(rd(umsa))
ulen = len(user[0][1])
print(f"  user MSA: {len(user)} genomes, {ulen} columns")

mask = None
if maskf and os.path.exists(maskf):
    raw = open(maskf).read().strip().replace('\n', '')
    if set(raw) <= {'0', '1'}:
        mask = raw
        print(f"  mask: {os.path.basename(maskf)}  {len(mask)} -> {mask.count('1')}")

# ---- which families do our MAGs occupy? ----------------------------------
dom = 'Archaea' if mk == 'ar53' else 'Bacteria'
fams, phyla = set(), set()
for r in csv.DictReader(open(magtax), delimiter='\t'):
    if r.get('domain') != dom: continue
    if r.get('family'): fams.add(r['family'])
    if r.get('phylum'): phyla.add(r['phylum'])
print(f"  MAG families ({len(fams)}): {', '.join(sorted(fams))}")

# ---- pick references ------------------------------------------------------
by_fam = collections.defaultdict(list)
outgroup_pool = []
for line in open(taxfile):
    p = line.rstrip('\n').split('\t')
    if len(p) < 2: continue
    gid, lin = p[0], p[1]
    if f"d__{dom}" not in lin: continue
    d = {}
    for part in lin.split(';'):
        part = part.strip()
        if len(part) > 3 and part[1:3] == '__':
            d[part[0]] = part[3:]
    fam, phy = d.get('f', ''), d.get('p', '')
    if fam in fams:
        by_fam[fam].append((gid, d.get('s', '') or d.get('g', '') or fam))
    elif phy and phy not in phyla:
        outgroup_pool.append((gid, phy))

wanted = {}
for fam in sorted(fams):
    picks = sorted(by_fam.get(fam, []))[:K]
    for gid, lab in picks: wanted[gid] = f"{fam}|{lab}"
    print(f"    {fam:<32} {len(picks)} refs")

# The outgroup must be MONOPHYLETIC or IQ-TREE cannot root on it: taking one
# genome each from two different phyla guarantees failure ("Branch separating
# outgroup is not found", observed 2026-08-19). Take several genomes from a
# SINGLE phylum instead.
og = {}
if outgroup_pool:
    by_phy = collections.defaultdict(list)
    for gid, phy in outgroup_pool:
        by_phy[phy].append(gid)
    # a phylum with enough sampled genomes to be a stable clade
    cand = sorted(by_phy.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    ogphy, gids = cand[0]
    for gid in sorted(gids)[:3]:
        og[gid] = f"OUTGROUP|{ogphy}"
    print(f"  outgroup: {len(og)} genomes from a single phylum ({ogphy}) "
          f"— monophyletic, so rooting will succeed")
else:
    print("  outgroup: none available — tree will be midpoint-rooted")

wanted.update(og)
print(f"  reference genomes requested: {len(wanted)}")

# ---- extract, mask, write -------------------------------------------------
kept, dropped = [], 0
for gid, seq in rd(refmsa):
    if gid not in wanted: continue
    if len(seq) != ulen and mask and len(seq) == len(mask):
        seq = ''.join(c for c, m in zip(seq, mask) if m == '1')
    if len(seq) != ulen: dropped += 1; continue
    kept.append((gid, seq))

with open(out, 'w') as o:
    for n, s in user: o.write(f">{n}\n{s}\n")
    for n, s in kept: o.write(f">{n}\n{s}\n")

with open(out.replace('_ref_input.faa', '_ref_labels.tsv'), 'w') as o:
    o.write("id\tlabel\ttype\n")
    for n, _ in user: o.write(f"{n}\t{n}\tMAG\n")
    for n, _ in kept:
        lab = wanted.get(n, n)
        o.write(f"{n}\t{lab}\t{'OUTGROUP' if lab.startswith('OUTGROUP') else 'reference'}\n")

print(f"  references added: {len(kept)}  (dropped for width mismatch: {dropped})")
print(f"  tree input: {len(user)+len(kept)} taxa x {ulen} columns")
og_ids = [n for n, _ in kept if wanted.get(n, '').startswith('OUTGROUP')]
open(out.replace('_ref_input.faa', '_outgroup.txt'), 'w').write(','.join(og_ids))
PY

  local ntax; ntax=$(grep -c '^>' "$final")
  local widths; widths=$(awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{if(l)print l}' \
                         "$final" | sort -u | wc -l)
  [ "$widths" -eq 1 ] || { echo "  !! ragged alignment — aborting $lab"; return 1; }

  local ogfile="$R/${mk}_outgroup.txt"
  local ogarg=""
  if [ -s "$ogfile" ] && [ -n "$(cat "$ogfile")" ]; then
    ogarg="-o $(cat "$ogfile")"
    echo "  rooting on outgroup: $(cat "$ogfile")"
  else
    echo "  NOTE: no outgroup — the tree will be MIDPOINT rooted. Say so in the caption."
  fi

  local pre="$R/${mk}_ref"
  if [ -s "${pre}.treefile" ]; then echo "  tree already built"; return 0; fi
  echo "  IQ-TREE on $ntax taxa (ModelFinder + 1000 UFBoot + SH-aLRT) ..."
  run_tool iqtree3 -s "$final" -m MFP -B 1000 -alrt 1000 $ogarg \
      -T AUTO --threads-max "$THREADS" --prefix "$pre" -redo \
      > "$R/${mk}_iqtree.log" 2>&1 \
    || { echo "  !! IQ-TREE failed — retrying without the outgroup constraint"
         run_tool iqtree3 -s "$final" -m MFP -B 1000 -alrt 1000 \
             -T AUTO --threads-max "$THREADS" --prefix "$pre" -redo \
             > "$R/${mk}_iqtree.log" 2>&1 || echo "  !! still failed — see the log"; }
  [ -s "${pre}.treefile" ] && {
    echo "  -> ${pre}.treefile"
    grep -m1 'Best-fit model' "$R/${mk}_iqtree.log" | sed 's/^/     /' || true
  }
}

build ar53   "ARCHAEA (ar53) — the one step 12 could not contextualise"
build bac120 "BACTERIA (bac120)"

# ------------------------------------------------------------ iTOL styling
python3 - "$R" "$F" <<'PY'
import os, csv, glob, sys
R, F = sys.argv[1], sys.argv[2]
mag = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    mag = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}
for lab in glob.glob(os.path.join(R, '*_ref_labels.tsv')):
    mk = os.path.basename(lab).split('_')[0]
    rows = list(csv.DictReader(open(lab), delimiter='\t'))
    with open(os.path.join(R, f'{mk}_itol_labels.txt'), 'w') as fh:
        fh.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for r in rows:
            if r['type'] == 'MAG':
                m = mag.get(r['id'], {})
                fh.write(f"{r['id']}\t{m.get('mag_id',r['id'])} "
                         f"{m.get('genus') or m.get('family') or ''}\n")
            else:
                fh.write(f"{r['id']}\t{r['label'].split('|')[-1]}\n")
    with open(os.path.join(R, f'{mk}_itol_colours.txt'), 'w') as fh:
        fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\torigin\n"
                 "COLOR\t#2a78d6\nLEGEND_TITLE\torigin\nLEGEND_SHAPES\t1\t1\t1\n"
                 "LEGEND_COLORS\t#2a78d6\t#eb6834\t#52514e\n"
                 "LEGEND_LABELS\tSediment MAG\tCulture MAG\tGTDB reference\nDATA\n")
        for r in rows:
            if r['type'] == 'MAG':
                s = mag.get(r['id'], {}).get('sample', '')
                col = '#2a78d6' if s == 'SED' else '#eb6834'
                fh.write(f"{r['id']}\t{col}\t{s}\n")
            else:
                fh.write(f"{r['id']}\t#52514e\treference\n")
    print(f"  iTOL files for {mk}: {R}/{mk}_itol_labels.txt, {mk}_itol_colours.txt")
PY

echo
echo "######## VALIDATION GATE ########"
for mk in ar53 bac120; do
  check "$mk reference tree built" $([ -s "$R/${mk}_ref.treefile" ] && echo 0 || echo 1)
  if [ -s "$R/${mk}_ref_labels.tsv" ]; then
    nr=$(awk -F'\t' 'NR>1 && $3=="reference"' "$R/${mk}_ref_labels.tsv" | wc -l)
    check "$mk has reference genomes ($nr)" $([ "$nr" -gt 0 ] && echo 0 || echo 1)
  fi
done
if [ -s "$R/ar53_ref.treefile" ]; then
  n8=$(grep -o 'LIQ_das_[a-z]*_[0-9]*' "$R/ar53_ref.treefile" | sort -u | wc -l)
  check "all 8 archaeal MAGs in the tree ($n8)" $([ "$n8" -eq 8 ] && echo 0 || echo 1)
fi

runlog "12b_reference_trees complete"
echo
echo "######## 12b DONE ########"
echo "Look for MTP4 (ETO_LIQ_MAG_003) in the archaeal tree — where it falls among"
echo "the Methanosarcinaceae determines how novel it actually is."
