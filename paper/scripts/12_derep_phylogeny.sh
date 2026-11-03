#!/usr/bin/env bash
# =============================================================================
# 12_derep_phylogeny.sh  --  dereplication (ANI) + phylogenomic trees
#
# TWO JOBS, both closing open audit findings.
#
#   JOB 1 — DEREPLICATION (ANI within the 29 MAGs)
#     Step 09 recovered Chakrabartyella FOUR times (ETO_SED_MAG_002,
#     ETO_LIQ_MAG_009, _012, _024) and Sedimentibacter TWICE (SED_001,
#     LIQ_025). Are those distinct species, or one genome recovered repeatedly?
#     The catalogue cannot be reported until that is known: "29 MAGs" means
#     something different if six of them are two organisms.
#
#     It also tests the enrichment link directly. A SEDIMENT MAG and a CULTURE
#     MAG at >=95% ANI are the SAME SPECIES — genomic proof that the organism
#     present in the lagoon is the one that grew in the enrichment. That is a
#     stronger statement than the read-mapping coverage result at step 08d,
#     because it is genome-to-genome rather than reads-to-reference.
#
#     Thresholds: 95% ANI = species (Jain et al. 2018), 99% = strain.
#
#   JOB 2 — PHYLOGENOMIC TREES  (audit finding F2)
#     workspace_1's reference tree was built on 22 MAGs, ALL from the liquid
#     cultures, with ZERO sediment genomes — so it could not support any claim
#     about the sediment community. Here both domains and both sample types are
#     included, from the GTDB-Tk marker alignments already computed at step 09
#     (bac120 for 21 bacteria, ar53 for 8 archaea). Nothing is realigned.
#
#     Where the GTDB reference MSA is available on disk, the closest reference
#     genome for each MAG is added, so the tree has context rather than being
#     MAGs floating alone. If it is not available the tree is still built, and
#     the script says plainly that it is a MAG-only topology.
#
# INPUT    07_mags/bins/*.fa, 07_mags/MAG_TAXONOMY.tsv
#          08_taxonomy/gtdbtk/align/*user_msa.fasta.gz
# OUTPUT   12_phylogeny/ani/, 12_phylogeny/trees/, derep_clusters.tsv
#
# RUNTIME  ANI seconds; IQ-TREE 20-60 min. Independent of 11b — safe to run
#          alongside the KOfam download.
# USAGE    ws2 && bash scripts/12_derep_phylogeny.sh 2>&1 | tee logs/12_phylo.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="12_phylogeny"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; T="$WS2/08_taxonomy"; P="$WS2/12_phylogeny"
mkdir -p "$P/ani" "$P/trees" "$P/msa"

echo "######## 12 — DEREPLICATION + PHYLOGENOMICS ########"
echo "  started: $(date -Is)  threads: $THREADS"
runlog "12_derep_phylogeny started"

NMAG=$(ls -1 "$F/bins"/*.fa 2>/dev/null | wc -l)
[ "$NMAG" -gt 0 ] || { echo "  !! no MAGs in $F/bins"; exit 1; }
echo "  MAGs: $NMAG"

# =====================================================================
# JOB 1 — ALL-vs-ALL ANI
# =====================================================================
echo
echo "######## ALL-vs-ALL ANI ($(date +%H:%M)) ########"
ls -1 "$F/bins"/*.fa > "$P/ani/genome_list.txt"
if [ -s "$P/ani/ani_all.tsv" ]; then
  echo "  already computed"
else
  record_version "fastANI" "${ENV_fastANI:-coverm_env}" fastANI --version
  run_tool fastANI --ql "$P/ani/genome_list.txt" --rl "$P/ani/genome_list.txt" \
      -o "$P/ani/ani_all.tsv" -t "$THREADS" 2>&1 | tail -3 \
    || echo "  !! fastANI failed"
fi
[ -s "$P/ani/ani_all.tsv" ] || { echo "  !! no ANI output"; exit 1; }
echo "  pairwise comparisons above the reporting floor: $(wc -l < "$P/ani/ani_all.tsv")"

python3 - "$P" "$F" <<'PY'
import sys, os, csv

P, F = sys.argv[1], sys.argv[2]
SPECIES, STRAIN, MINAF = 95.0, 99.0, 0.30

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}

def stem(x):
    b = os.path.basename(x)
    return b[:-3] if b.endswith('.fa') else b

pairs = []
for line in open(os.path.join(P, 'ani', 'ani_all.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5: continue
    q, r = stem(f[0]), stem(f[1])
    if q == r: continue
    ani = float(f[2]); af = int(f[3]) / int(f[4]) if int(f[4]) else 0
    if af < MINAF: continue
    pairs.append((q, r, ani, af))

# ---------------- single-linkage clustering at the species threshold
names = sorted(tax) or sorted({x for p_ in pairs for x in p_[:2]})
parent = {n: n for n in names}
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

for q, r, ani, af in pairs:
    if ani >= SPECIES and q in parent and r in parent:
        union(q, r)

clusters = {}
for n in names:
    clusters.setdefault(find(n), []).append(n)

print()
print("=" * 96)
print(f"  DEREPLICATION — species clusters at >={SPECIES}% ANI (align. fraction >={MINAF:.0%})")
print("=" * 96)
multi = {k: v for k, v in clusters.items() if len(v) > 1}
print(f"  {len(names)} MAGs  ->  {len(clusters)} species-level clusters "
      f"({len(multi)} contain more than one MAG)")

def q(n):
    t = tax.get(n, {})
    try: return float(t.get('completeness', 0)) - 5*float(t.get('contamination', 0))
    except ValueError: return 0.0

rows = []
if multi:
    print()
    for k, members in sorted(multi.items(), key=lambda kv: -len(kv[1])):
        best = max(members, key=q)
        gen = tax.get(best, {}).get('genus', '-') or '-'
        samples = {tax.get(m, {}).get('sample', '?') for m in members}
        cross = 'SED' in samples and 'LIQ' in samples
        print(f"  cluster: {gen}  ({len(members)} MAGs)"
              + ("   <-- SPANS SEDIMENT AND CULTURES" if cross else ""))
        for m in sorted(members, key=lambda m: -q(m)):
            t = tax.get(m, {})
            anis = [f"{a:.2f}" for x, y, a, _ in pairs
                    if {x, y} <= set(members) and (x == m or y == m)]
            print(f"     {'*' if m == best else ' '} {t.get('mag_id',m):<17} "
                  f"{t.get('sample','?'):<4} comp {t.get('completeness','?'):>6} "
                  f"cont {t.get('contamination','?'):>5}  ANI to cluster: "
                  f"{'/'.join(sorted(set(anis))[:4])}")
        print()
else:
    print("  No two MAGs exceed 95% ANI — all 29 are distinct species.")

for k, members in clusters.items():
    best = max(members, key=q)
    for m in members:
        t = tax.get(m, {})
        rows.append(dict(mag_id=t.get('mag_id',''), working_name=m,
                         sample=t.get('sample',''), genus=t.get('genus',''),
                         cluster=tax.get(best,{}).get('mag_id', best),
                         representative=('yes' if m == best else 'no'),
                         cluster_size=len(members)))
with open(os.path.join(P, 'derep_clusters.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(sorted(rows, key=lambda r: r['mag_id']))

# ---------------- cross-sample species matches: the enrichment link
print("=" * 96)
print("  SEDIMENT <-> CULTURE GENOME PAIRS  (the enrichment link, genome-level)")
print("=" * 96)
cross = [(qn, rn, a, af) for qn, rn, a, af in pairs
         if a >= SPECIES
         and tax.get(qn, {}).get('sample') == 'SED'
         and tax.get(rn, {}).get('sample') == 'LIQ']
if cross:
    for qn, rn, a, af in sorted(cross, key=lambda x: -x[2]):
        print(f"  {tax.get(qn,{}).get('mag_id',qn):<17} (SED)  <->  "
              f"{tax.get(rn,{}).get('mag_id',rn):<17} (LIQ)   ANI {a:.2f}%  AF {af:.2f}")
        print(f"     {tax.get(qn,{}).get('genus','-')} / {tax.get(rn,{}).get('genus','-')}")
    print()
    print("  These are the SAME SPECIES recovered from both the lagoon sediment and")
    print("  the enrichment. Genome-level proof that the cultures grew resident")
    print("  organisms — stronger than the read-mapping evidence at step 08d.")
else:
    print("  No sediment-culture pair reaches 95% ANI.")
    print("  Note this does NOT contradict the coverage result (step 08d): a")
    print("  sediment organism can be abundant in the cultures without a CULTURE")
    print("  MAG having been recovered for it, and vice versa.")

# ---------------- closest relatives generally
print()
print("=" * 96)
print("  NEAREST NEIGHBOUR PER MAG (highest ANI to any other MAG)")
print("=" * 96)
best_hit = {}
for qn, rn, a, af in pairs:
    if qn not in best_hit or a > best_hit[qn][1]: best_hit[qn] = (rn, a, af)
for n in sorted(names, key=lambda n: tax.get(n,{}).get('mag_id','')):
    t = tax.get(n, {})
    if n in best_hit:
        rn, a, af = best_hit[n]
        flag = ' SAME SPECIES' if a >= SPECIES else (' same genus?' if a >= 85 else '')
        print(f"  {t.get('mag_id',n):<17} -> {tax.get(rn,{}).get('mag_id',rn):<17} "
              f"{a:>6.2f}%  AF {af:.2f}{flag}")
    else:
        print(f"  {t.get('mag_id',n):<17} -> (no relative above the ANI reporting floor)")
PY

# =====================================================================
# JOB 2 — PHYLOGENOMIC TREES
# =====================================================================
echo
echo "######## PHYLOGENOMIC TREES ($(date +%H:%M)) ########"
ALIGN="$T/gtdbtk/align"
[ -d "$ALIGN" ] || { echo "  !! $ALIGN missing — step 09 output required"; exit 1; }

record_version "iqtree" "${ENV_iqtree:-iqtree3}" iqtree3 --version

build_tree() {  # build_tree <bac120|ar53> <label>
  local mk="$1" lab="$2"
  local src="" g
  for g in "$ALIGN/gtdbtk.${mk}.user_msa.fasta.gz" "$ALIGN/gtdbtk.${mk}.user_msa.fasta"; do
    [ -s "$g" ] && { src="$g"; break; }
  done
  [ -n "$src" ] || { echo "  $lab: no user MSA found — skipping"; return 0; }

  local msa="$P/msa/${mk}_user.faa"
  if [ ! -s "$msa" ]; then
    case "$src" in *.gz) zcat "$src" > "$msa";; *) cp "$src" "$msa";; esac
  fi
  local n; n=$(grep -c '^>' "$msa")
  echo
  echo "--- $lab: $n genomes in the alignment"
  if [ "$n" -lt 4 ]; then
    echo "    fewer than 4 taxa — a tree is not meaningful; skipping"
    return 0
  fi

  # Opportunistically add GTDB reference genomes for context. The reference MSA
  # ships with the GTDB release; if present, pull the closest reference of each
  # MAG so the tree is not MAGs floating without context.
  local refmsa=""
  for g in "$GTDBTK_DATA_PATH"/msa/*${mk}*.faa "$GTDBTK_DATA_PATH"/msa/*${mk}*.faa.gz; do
    [ -s "$g" ] && { refmsa="$g"; break; }
  done
  local final="$P/msa/${mk}_tree_input.faa"
  if [ -n "$refmsa" ] && [ ! -s "$final" ]; then
    echo "    reference MSA found: $(basename "$refmsa")"
    python3 - "$msa" "$refmsa" "$F/MAG_TAXONOMY.tsv" "$final" "$GTDBTK_DATA_PATH" "$mk" <<'PY'
import sys, gzip, csv, os, glob
msa, refmsa, taxf, out, gtdbdir, mk = sys.argv[1:7]

def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    name, buf = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name: yield name, ''.join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name: yield name, ''.join(buf)

user = list(read_fasta(msa))
if not user:
    print("    !! user MSA empty"); raise SystemExit(1)
ulen = len(user[0][1])
print(f"    user alignment width: {ulen} columns")

wanted = set()
if os.path.exists(taxf):
    for r in csv.DictReader(open(taxf), delimiter='\t'):
        for tok in (r.get('closest_ref') or '').replace(',', ' ').split():
            if tok and tok not in ('N/A', 'none'):
                wanted.add(tok.split('~')[0])

# The GTDB reference MSA is the FULL concatenated alignment; GTDB-Tk's user_msa
# is MASKED to a subset of columns. Adding one to the other unmasked produces
# unequal sequence lengths and IQ-TREE rejects the alignment — this is exactly
# why the bacterial tree failed on 2026-08-19 while the archaeal one (which had
# no references to add) succeeded. Apply GTDB's own mask to the references.
mask = None
for pat in (f"{gtdbdir}/masks/*{mk}*.mask", f"{gtdbdir}/masks/*{mk}*"):
    hits = sorted(glob.glob(pat))
    if hits:
        raw = open(hits[0]).read().strip().replace('\n', '')
        if set(raw) <= {'0', '1'}:
            mask = raw
            print(f"    mask: {os.path.basename(hits[0])} "
                  f"({len(mask)} cols -> {mask.count('1')} kept)")
        break

refs = []
for name, seq in read_fasta(refmsa):
    tag = name.split()[0]
    if not any(w in tag for w in wanted):
        continue
    if len(seq) != ulen and mask and len(seq) == len(mask):
        seq = ''.join(c for c, m in zip(seq, mask) if m == '1')
    if len(seq) != ulen:
        continue                      # still mismatched — drop rather than break
    refs.append((tag.replace(' ', '_'), seq))

with open(out, 'w') as o:
    for n, s in user: o.write(f">{n}\n{s}\n")
    for n, s in refs: o.write(f">{n}\n{s}\n")

print(f"    reference genomes added: {len(refs)}"
      + ("" if refs else "  (none matched the user alignment width — MAG-only tree)"))
PY
  fi
  [ -s "$final" ] || cp "$msa" "$final"
  local ntot; ntot=$(grep -c '^>' "$final")
  # every sequence must be the same width or IQ-TREE aborts
  local widths; widths=$(awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{if(l)print l}' \
                         "$final" | sort -u | tr '\n' ' ')
  echo "    tree input: $ntot sequences, alignment width(s): $widths"
  if [ "$(echo "$widths" | wc -w)" -gt 1 ]; then
    echo "    !! ragged alignment — falling back to the MAG-only alignment"
    cp "$msa" "$final"
    ntot=$(grep -c '^>' "$final")
    echo "    tree input now: $ntot sequences (MAGs only)"
  fi

  local pre="$P/trees/${mk}"
  if [ -s "${pre}.treefile" ]; then echo "    tree already built"; return 0; fi
  echo "    running IQ-TREE (ModelFinder + 1000 UFBoot + SH-aLRT) ..."
  run_tool iqtree3 -s "$final" -m MFP -B 1000 -alrt 1000 \
      -T AUTO --threads-max "$THREADS" --prefix "$pre" -redo \
      > "$P/trees/${mk}_iqtree.log" 2>&1 \
    || echo "    !! IQ-TREE failed — see $P/trees/${mk}_iqtree.log"
  [ -s "${pre}.treefile" ] && echo "    -> ${pre}.treefile"
}

build_tree bac120 "BACTERIA (bac120)"
build_tree ar53   "ARCHAEA (ar53)"

# ---------------------------------------------------- iTOL annotation
python3 - "$F" "$P" <<'PY'
import os, csv, sys
F, P = sys.argv[1], sys.argv[2]
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if not os.path.exists(p): raise SystemExit(0)
rows = list(csv.DictReader(open(p), delimiter='\t'))
out = os.path.join(P, 'trees', 'itol_labels.txt')
with open(out, 'w') as fh:
    fh.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for r in rows:
        lab = f"{r['mag_id']} {r.get('genus') or r.get('family') or r.get('phylum') or ''}".strip()
        fh.write(f"{r['working_name']}\t{lab}\n")
out2 = os.path.join(P, 'trees', 'itol_sample_colours.txt')
with open(out2, 'w') as fh:
    fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tsample\n"
             "COLOR\t#2a78d6\nLEGEND_TITLE\tsample\nLEGEND_SHAPES\t1\t1\n"
             "LEGEND_COLORS\t#2a78d6\t#eb6834\nLEGEND_LABELS\tSediment\tLiquid cultures\nDATA\n")
    for r in rows:
        col = '#2a78d6' if r.get('sample') == 'SED' else '#eb6834'
        fh.write(f"{r['working_name']}\t{col}\t{r.get('sample','')}\n")
print(f"  iTOL annotation files: {out}")
print(f"                         {out2}")
PY

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "ANI table"            $([ -s "$P/ani/ani_all.tsv" ] && echo 0 || echo 1)
check "dereplication table"  $([ -s "$P/derep_clusters.tsv" ] && echo 0 || echo 1)
nc=$(( $(wc -l < "$P/derep_clusters.tsv") - 1 ))
check "every MAG assigned to a cluster ($nc/$NMAG)" $([ "$nc" -eq "$NMAG" ] && echo 0 || echo 1)
tb=$([ -s "$P/trees/bac120.treefile" ] && echo 0 || echo 1)
ta=$([ -s "$P/trees/ar53.treefile" ] && echo 0 || echo 1)
check "bacterial tree built" "$tb"
check "archaeal tree built"  "$ta"
# audit F2: the tree must contain sediment genomes, which workspace_1's did not
if [ -s "$P/trees/bac120.treefile" ]; then
  sed_in=$(grep -o 'SED_[A-Za-z0-9_]*' "$P/trees/bac120.treefile" | sort -u | wc -l)
  check "bacterial tree contains sediment MAGs ($sed_in) — audit F2" \
    $([ "$sed_in" -gt 0 ] && echo 0 || echo 1)
fi

runlog "12_derep_phylogeny complete"
echo
echo "######## 12 DONE ########"
