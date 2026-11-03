#!/usr/bin/env bash
# =============================================================================
# 24_deep_taxonomy.sh  --  how to see MORE taxa in the same data
#
# THE SHORT ANSWER: YOU ALREADY HAVE THEM. SingleM produced 20,919 OTU rows for
# the sediment and 15,591 for the cultures, and step 22 printed the top twenty.
# Everything below the twentieth line has never been looked at. That is by far
# the largest untapped taxonomic resource in this project and it costs nothing
# to mine.
#
# WHY EACH METHOD SEES A DIFFERENT NUMBER OF TAXA
#   Kraken2      0.47% of sediment reads classified. It matches nucleotide
#                k-mers against RefSeq genomes. An organism with no sequenced
#                relative is invisible, and ONT error breaks exact k-mers.
#                Structurally blind to exactly this community.
#   SILVA SSU    1,582 reads. Not a sensitivity problem — a LENGTH problem. The
#                filter required a >=1000 bp alignment, but sediment reads have
#                a MEDIAN OF 364 bp. Most reads are physically too short to
#                produce a 1 kb SSU alignment no matter what lives there.
#   SingleM      20,919 OTUs. Extracts a short conserved window from ~60
#                single-copy ribosomal genes. Needs no reference genome for the
#                organism and no long read. This is why it wins here.
#
# WHAT THIS SCRIPT DOES — three things, all using data already on disk
#   A. MINE SINGLEM PROPERLY. Rank-by-rank taxon counts, how many lineages are
#      unnamed at each rank, richness per marker, and the full ranked table
#      written out rather than the top 20.
#   B. RE-RUN SILVA AT A LENGTH APPROPRIATE TO THE READS. Identity stays strict
#      at 92% — that threshold was hard-won and is not the problem. Alignment
#      length drops 1000 -> 400, which is still long enough to place a read
#      taxonomically but no longer discards the majority of the library.
#      Written to a SEPARATE directory so the strict profile stays intact and
#      both can be reported.
#   C. THE ASSEMBLY'S OWN SSU GENES as a third, independent view (138 SED +
#      323 LIQ full-length genes, not limited to binned genomes).
#
#   It also prints a costed menu of the heavier options that would go further.
#
# OUTPUT   24_deep_taxonomy/
#     singlem_ranks.tsv        taxa per rank, per sample, named vs unnamed
#     singlem_full_SED.tsv     every lineage ranked, not the top 20
#     singlem_full_LIQ.tsv
#     silva_len400/            the relaxed SSU profile
#     COMPARISON.txt           how many taxa each method sees, and why
#
# RUNTIME  1-2 h (the SILVA re-map dominates).
# USAGE    ws2 && bash scripts/24_deep_taxonomy.sh 2>&1 | tee logs/24_deeptax.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="24_deeptax"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
G="$WS2/22_singlem"; S="$WS2/13_profile/silva"; PR="$WS2/13_profile"
D="$WS2/24_deep_taxonomy"; mkdir -p "$D"

echo "######## 24 — DEEPER TAXONOMY FROM THE SAME DATA ########"
echo "  started: $(date -Is)"
runlog "24_deep_taxonomy started"

# =====================================================================
# B. SILVA at a length the reads can actually satisfy
# =====================================================================
echo
echo "######## B. SILVA SSU at MINLEN=400 (identity stays 92%) ########"
echo "  Sediment median read is 364 bp. Requiring a 1000 bp SSU alignment"
echo "  discards most of the library before taxonomy is even attempted."
S4="$D/silva_len400"; mkdir -p "$S4"
cp -n "$S/silva.mmi" "$S4/" 2>/dev/null || true
for pair in "SED:$SED_READS" "LIQ:$LIQ_READS"; do
  code="${pair%%:*}"; reads="${pair#*:}"
  src="$S/${code}_vs_silva.paf.gz"
  dst="$S4/${code}_vs_silva.paf.gz"
  if [ -s "$src" ]; then
    # the PAF already holds every alignment; only the FILTER changes, so there
    # is no need to re-map. Re-using it also guarantees the two profiles differ
    # only in the threshold, which is the point of the comparison.
    echo "  $code: reusing existing alignments (only the filter changes)"
    cp -n "$src" "$dst" 2>/dev/null || true
  else
    echo "  $code: mapping ($(date +%H:%M)) ..."
    run_tool minimap2 -x map-ont -t "$THREADS" --secondary=no -c \
        "$S4/silva.mmi" "$reads" 2>/dev/null | gzip > "$dst"
  fi
done
# The strict profile was overwritten by the first version of this script, so
# rebuild it from the cached PAF (seconds) before writing the relaxed one.
echo "  rebuilding the STRICT (>=1000 bp) profile that was clobbered ..."
MINID=0.92 MINLEN=1000 SILVA_OUT="$S" \
  bash "$(dirname "$0")/13b_silva_ssu.sh" >/dev/null 2>&1 || true
echo "  writing the RELAXED (>=400 bp) profile to $S4 ..."
MINID=0.92 MINLEN=400 SILVA_OUT="$S4" \
  bash "$(dirname "$0")/13b_silva_ssu.sh" >/dev/null 2>&1 || true

# =====================================================================
python3 - "$G" "$S" "$S4" "$PR" "$D" "$WS2" <<'PY'
import sys, os, csv, gzip, glob
from collections import Counter, defaultdict

G, S, S4, PR, D, WS2 = sys.argv[1:7]
RANKS = ['domain','phylum','class','order','family','genus','species']
PFX   = ['d__','p__','c__','o__','f__','g__','s__']

# =====================================================================
# A. MINE SINGLEM
# =====================================================================
print()
print("=" * 100)
print("  A. SINGLEM — THE FULL PICTURE, NOT THE TOP 20")
print("=" * 100)

rank_rows = []
for code in ('SED', 'LIQ'):
    p = os.path.join(G, f'{code}_otu_table.tsv')
    if not os.path.exists(p):
        print(f"  {code}: no OTU table"); continue
    rows = list(csv.DictReader(open(p), delimiter='\t'))
    def cov(r):
        for k in ('coverage', 'num_hits'):
            try: return float(r.get(k) or 0)
            except (TypeError, ValueError): pass
        return 0.0

    # distinct taxa at each rank, and how many are UNNAMED (alphanumeric codes
    # like g__JAOZPL01 are real lineages with no cultured representative)
    seen = {r: set() for r in RANKS}
    unnamed = {r: set() for r in RANKS}
    covr = {r: Counter() for r in RANKS}
    markers = Counter(); seqs = set()
    for r in rows:
        t = [x.strip() for x in (r.get('taxonomy') or '').split(';') if x.strip()]
        markers[r.get('gene', '?')] += 1
        if r.get('sequence'): seqs.add(r['sequence'])
        for i, pf in enumerate(PFX):
            v = next((x[len(pf):] for x in t if x.startswith(pf)), '')
            if not v: continue
            seen[RANKS[i]].add(v)
            covr[RANKS[i]][v] += cov(r)
            # a named lineage has a lowercase letter and no long digit run
            if v[:1].isupper() and any(ch.isdigit() for ch in v) and len(v) > 5:
                unnamed[RANKS[i]].add(v)
            elif v.startswith(('UBA','JA','CG','GCA','GW','RBG','B','AB','FEB','DAO','CAY','WJM','BS')) and any(ch.isdigit() for ch in v):
                unnamed[RANKS[i]].add(v)

    print(f"\n  {code}: {len(rows)} OTU rows, {len(seqs)} distinct marker windows, "
          f"{len(markers)} marker genes")
    print(f"  {'rank':<10}{'distinct taxa':>15}{'no cultured name':>19}{'% unnamed':>12}")
    for r in RANKS:
        n, u = len(seen[r]), len(unnamed[r])
        print(f"  {r:<10}{n:>15}{u:>19}{100*u/max(1,n):>11.1f}%")
        rank_rows.append(dict(sample=code, rank=r, distinct=n, unnamed=u,
                              pct_unnamed=f"{100*u/max(1,n):.1f}"))

    # write the FULL ranked table
    out = os.path.join(D, f'singlem_full_{code}.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['rank','taxon','coverage','pct_of_rank'])
        for r in RANKS:
            tot = sum(covr[r].values()) or 1
            for tx, c in covr[r].most_common():
                w.writerow([r, tx, f"{c:.4f}", f"{100*c/tot:.4f}"])
    print(f"    full ranked table -> {out}")

    g = covr['genus']
    gt = sum(g.values()) or 1
    print(f"\n    genera 21-45 by coverage (the ones step 22 never printed):")
    for tx, c in g.most_common()[20:45]:
        print(f"      {100*c/gt:>6.2f}%  {tx}")

with open(os.path.join(D, 'singlem_ranks.tsv'), 'w', newline='') as fh:
    if rank_rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rank_rows[0]))
        w.writeheader(); w.writerows(rank_rows)

# =====================================================================
# B/C. HOW MANY TAXA DOES EACH METHOD SEE?
# =====================================================================
def silva_genera(d, code):
    p = os.path.join(d, f'{code}_ssu_profile.tsv')
    if not os.path.exists(p): return None, 0
    gs, reads = set(), 0
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r.get('rank') == 'genus':
            gs.add(r['taxon'])
            try: reads += int(r['reads'])
            except (TypeError, ValueError): pass
    return gs, reads

print()
print("=" * 100)
print("  HOW MANY TAXA EACH METHOD SEES — same DNA, same sample")
print("=" * 100)
lines = []
for code in ('SED', 'LIQ'):
    strict, sr = silva_genera(S, code)
    loose, lr = silva_genera(S4, code)
    sm = os.path.join(G, f'{code}_otu_table.tsv')
    smg = set()
    if os.path.exists(sm):
        for r in csv.DictReader(open(sm), delimiter='\t'):
            t = [x.strip() for x in (r.get('taxonomy') or '').split(';') if x.strip()]
            v = next((x[3:] for x in t if x.startswith('g__')), '')
            if v: smg.add(v)
    print(f"\n  {code}")
    print(f"    SILVA  >=1000 bp : {len(strict or []):>6} genera from {sr:>8,} reads")
    if loose is not None:
        print(f"    SILVA  >= 400 bp : {len(loose):>6} genera from {lr:>8,} reads"
              f"   ({lr/max(1,sr):.1f}x more reads)")
        new = loose - (strict or set())
        print(f"      genera visible ONLY at 400 bp: {len(new)}")
        if new:
            print("      e.g. " + ', '.join(sorted(new)[:12]))
    print(f"    SingleM          : {len(smg):>6} genera")
    lines.append((code, len(strict or []), len(loose or []), len(smg)))

print()
print("=" * 100)
print("  WHAT ELSE COULD BE DONE, AND WHAT EACH WOULD COST")
print("=" * 100)
print("""
  1. KAIJU — protein-level read classification.            HIGH VALUE
     Kraken2 matches nucleotide k-mers; a 15% divergent organism shares almost
     no exact 31-mers with its nearest relative and ONT indels break them
     further. Kaiju translates reads and does protein-level maximum-exact-match
     search, which tolerates both. On environmental samples it typically
     classifies several-fold more reads than Kraken2.
       conda create -y -p ENVS/kaiju_env -c conda-forge -c bioconda kaiju
       kaiju-makedb -s progenomes     (~15 GB, or nr_euk ~65 GB)
     Expect the sediment's 0.47% to rise substantially. Worth doing.

  2. A GTDB-BUILT KRAKEN2 DATABASE instead of RefSeq standard-16.  HIGH VALUE
     standard-16 contains RefSeq complete genomes — cultured organisms. GTDB
     r220/r226 includes hundreds of thousands of environmental MAGs, which is
     precisely the reference material this community needs. Struo2 builds one,
     or pre-built GTDB Kraken2 databases are distributed. Check it fits RAM
     before starting (step 13 will refuse if it does not).

  3. MMseqs2 TAXONOMY ON CONTIGS — protein-level LCA.       MEDIUM VALUE
     MMseqs2 18.8cc5c is already installed. Assigns taxonomy to every contig by
     protein-level LCA, far more sensitive than the nucleotide Kraken2 pass that
     classified only 1.07% of sediment contigs. Needs a target database
     (GTDB or UniRef90, tens of GB).

  4. 18S / EUKARYOTES.                                      ALREADY PARTLY DONE
     SILVA SSURef contains 18S as well as 16S, so eukaryotes are already in the
     profile, and SingleM independently reports Viridiplantae at 2.4-2.8%.
     Report it rather than adding a tool.

  5. WHAT WILL NOT HELP.
     Re-assembling deeper, or a fourth binner. Recovery is limited by read
     length (sediment median 364 bp), which is a property of the DNA. Step 06b
     measured this: above 10 kb, sediment bins as well as the cultures do.
     Nothing downstream fixes fragmented input.
""")
with open(os.path.join(D, 'COMPARISON.txt'), 'w') as fh:
    for code, a, b, c in lines:
        fh.write(f"{code}\tSILVA1000\t{a}\nSILVA400\t{b}\nSingleM\t{c}\n")
PY

echo
echo "######## VALIDATION GATE ########"
check "singlem rank table"  $([ -s "$D/singlem_ranks.tsv" ] && echo 0 || echo 1)
check "SED full table"      $([ -s "$D/singlem_full_SED.tsv" ] && echo 0 || echo 1)
check "LIQ full table"      $([ -s "$D/singlem_full_LIQ.tsv" ] && echo 0 || echo 1)

runlog "24_deep_taxonomy complete"
echo
echo "######## 24 DONE ########"
