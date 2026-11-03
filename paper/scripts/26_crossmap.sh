#!/usr/bin/env bash
# =============================================================================
# 26_crossmap.sh  --  are the organisms we cultured actually IN the sediment?
#
# THE GAP THIS FILLS
#   Step 05 mapped each read set to its OWN assembly. So the catalogue carries
#   cov_SED for the 4 sediment MAGs, cov_LIQ for the 25 culture MAGs, and
#   NOTHING off the diagonal. The most important question about an enrichment
#   experiment has therefore never been asked with data:
#
#       Did the enrichments recover abundant residents of the lagoon, or did
#       they bloom something rare that was barely there to begin with?
#
#   Both look identical in a self-mapped catalogue. This script maps SEDIMENT
#   reads against all 29 genomes and CULTURE reads against all 29, filling the
#   off-diagonal.
#
# THE METHODOLOGICAL POINT, WHICH IS THE WHOLE ANALYSIS
#   PRESENCE MUST BE CALLED ON BREADTH, NOT DEPTH.
#
#   A genome at 0.5x with 5% of its length covered is reads piling onto
#   ribosomal proteins and other conserved regions from RELATIVES — it is
#   cross-mapping, not presence. The same 0.5x spread over 60% of the genome is
#   a real low-abundance population. Mean depth cannot distinguish them, and
#   reporting depth alone would let every culture MAG appear "detectable" in
#   sediment purely through conserved genes.
#
#   Thresholds used, both stated rather than buried:
#     BREADTH_PRESENT = 0.30   (30% of the genome covered)
#     BREADTH_STRONG  = 0.60
#   These are conventions, not measurements. The script prints the full
#   distribution so you can see where the natural break falls in YOUR data and
#   move them if it disagrees.
#
#   CoverM is run with --min-covered-fraction 0. Its default of 10% silently
#   reports zero for anything below, which would erase precisely the
#   low-abundance signal this analysis exists to detect.
#
# CONTROLS
#   Self-mapping is the positive control and it is free: the 4 sediment MAGs
#   MUST be present in sediment and the 25 culture MAGs MUST be present in
#   culture. If any fails, the mapping is broken and no off-diagonal number
#   means anything. The script refuses to interpret if the controls fail.
#
# READS
#   Sediment uses 00_input/sediment_reads_clean.fastq.gz if it exists — the
#   amplicon-free set from 13f. If it does not, the script says so loudly and
#   uses $SED_READS, which still contains the RBK_test contamination.
#
# OUTPUT   26_crossmap/cross_abundance.tsv   one row per MAG, both samples
#          26_crossmap/{SED,LIQ}_vs_all.tsv  raw CoverM output
# RUNTIME  30-90 min depending on read count. Mapping dominates.
# USAGE    ws2 && bash scripts/26_crossmap.sh 2>&1 | tee logs/26.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="26_crossmap"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
X="$WS2/26_crossmap"; F="$WS2/07_mags"; mkdir -p "$X"

BREADTH_PRESENT="${BREADTH_PRESENT:-0.30}"
BREADTH_STRONG="${BREADTH_STRONG:-0.60}"

echo "######## 26 — CROSS-MAPPING SEDIMENT <-> CULTURES ########"
echo "  started: $(date -Is)"
runlog "26_crossmap started"

# ---- locate the MAG fasta directory rather than assuming it -----------------
MAGDIR=""; MAGEXT=""
for ext in fa fna fasta; do
  d=$(find "$F" -maxdepth 3 -name "*.${ext}" -printf '%h\n' 2>/dev/null \
      | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
  [ -n "$d" ] || continue
  n=$(find "$d" -maxdepth 1 -name "*.${ext}" | wc -l)
  if [ "$n" -ge 20 ]; then MAGDIR="$d"; MAGEXT="$ext"; break; fi
done
[ -n "$MAGDIR" ] || { echo "  !! could not find a directory under 07_mags holding >=20 genome fastas"; exit 1; }
NMAG=$(find "$MAGDIR" -maxdepth 1 -name "*.${MAGEXT}" | wc -l)
echo "  genomes: $NMAG x .${MAGEXT} in ${MAGDIR#$WS2/}"
[ "$NMAG" -eq 29 ] || echo "  !! expected 29 genomes, found $NMAG — check before trusting the table"

# ---- reads -----------------------------------------------------------------
CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
if [ -s "$CLEAN" ]; then
  SED_IN="$CLEAN"; echo "  sediment reads: 00_input/sediment_reads_clean.fastq.gz (amplicon removed)"
else
  SED_IN="$SED_READS"
  echo "  !! CLEAN SEDIMENT READS NOT FOUND."
  echo "     Falling back to \$SED_READS, which STILL CONTAINS the RBK_test"
  echo "     amplicon library identified at step 13f. Every sediment abundance"
  echo "     below is then computed against a denominator inflated by ~700k"
  echo "     contaminant reads. Run 13f with WRITE_CLEAN=1 first."
fi
LIQ_IN="$LIQ_READS"
for f in "$SED_IN" "$LIQ_IN"; do
  [ -s "$f" ] || { echo "  !! missing reads: $f"; exit 1; }
done

# ---- map -------------------------------------------------------------------
# --min-covered-fraction 0 is REQUIRED. CoverM's default of 10 reports a hard
# zero below that, which would delete the low-abundance signal this whole
# analysis is looking for and make absence indistinguishable from "under 10%".
run_cov () {
  # NEVER reference a variable in the same `local` statement that declares it.
  # `local tag="$1" out="$X/${tag}_..."` looks fine and is not: bash declares
  # every name in the statement as an unset local BEFORE evaluating the
  # assignments, so ${tag} on the same line resolves to the shadowing local
  # that has not been assigned yet — and under `set -u` that is a fatal
  # "unbound variable" on a line where the variable is visibly being set.
  local tag="$1"
  local reads="$2"
  local out="$X/${tag}_vs_all.tsv"
  if [ -s "$out" ]; then echo "  [$tag] cached"; return 0; fi
  echo "  [$tag] mapping $(basename "$reads") against all $NMAG genomes ..."
  run_tool coverm genome \
      --single "$reads" \
      --genome-fasta-directory "$MAGDIR" -x "$MAGEXT" \
      -p minimap2-ont \
      -m mean covered_fraction relative_abundance \
      --min-covered-fraction 0 \
      -t "$THREADS" -o "$out" 2>&1 | tail -3
  [ -s "$out" ] || { echo "  !! $tag mapping produced nothing"; return 1; }
}
run_cov SED "$SED_IN" || exit 1
run_cov LIQ "$LIQ_IN" || exit 1

python3 - "$X" "$F" "$BREADTH_PRESENT" "$BREADTH_STRONG" <<'PY'
import sys, os, csv, re

X, F, BP, BS = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
by = {}
for wn, r in tax.items():
    for k in (wn, r.get('mag_id',''), r.get('qc_bin_name',''), r.get('raw_binner_file','')):
        if k:
            b = os.path.splitext(os.path.basename(k))[0]
            by[b.lower()] = wn
def resolve(name):
    b = os.path.splitext(os.path.basename((name or '').strip()))[0]
    return by.get(b.lower())
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus','') or '-')
def phy(m): return ((tax.get(m, {}) or {}).get('phylum','') or '')
def origin(m): return 'SED' if '_SED_' in mid(m) else 'LIQ'

def read_coverm(p):
    """CoverM names value columns after the input file, so a fixed field map
    cannot work — that is what lost the coverage join at step 25. Find columns
    by what their header CONTAINS."""
    rows = list(csv.DictReader(open(p), delimiter='\t'))
    if not rows: return {}
    hdr = list(rows[0])
    gcol = next((h for h in hdr if h.strip().lower() in ('genome','bin')), hdr[0])
    def pick(*pats):
        for h in hdr:
            hl = h.lower()
            if all(x in hl for x in pats): return h
        return None
    dcol = pick('mean'); bcol = pick('covered', 'fraction'); rcol = pick('relative')
    out = {}
    for r in rows:
        g = (r.get(gcol) or '').strip()
        if not g or g.lower() == 'unmapped': continue
        k = resolve(g)
        if not k: continue
        def fl(c):
            try: return float(r.get(c) or 0)
            except ValueError: return 0.0
        out[k] = dict(depth=fl(dcol), breadth=fl(bcol), relab=fl(rcol))
    return out

sed = read_coverm(os.path.join(X, 'SED_vs_all.tsv'))
liq = read_coverm(os.path.join(X, 'LIQ_vs_all.tsv'))
mags = sorted(set(sed) | set(liq), key=mid)
if not mags:
    print("  !! no CoverM row resolved to a MAG — check genome file naming"); sys.exit(1)

# CoverM reports covered_fraction as a fraction (0-1) in most versions but as a
# percentage in some. Detect rather than assume.
vals = [v['breadth'] for v in list(sed.values()) + list(liq.values())]
scale = 100.0 if max(vals, default=0) > 1.5 else 1.0
if scale == 100.0:
    print("  note: covered_fraction is on a 0-100 scale in this CoverM build; "
          "thresholds rescaled")
def br(d, m): return (d.get(m, {}).get('breadth', 0.0)) / scale

print()
print("=" * 112)
print("  1. EVERY GENOME IN BOTH SAMPLES")
print("=" * 112)
print(f"  present = breadth >= {BP:.2f}   strong = breadth >= {BS:.2f}")
print("  Breadth, not depth, decides presence. Depth is shown so you can see the")
print("  pairs where they disagree — those are the conserved-region artefacts.")
print()
print(f"  {'MAG ID':<17}{'genus':<19}{'from':<5}"
      f"{'SED depth':>10}{'SED brdth':>10}{'SED%':>7}"
      f"{'LIQ depth':>10}{'LIQ brdth':>10}{'LIQ%':>7}  call")
rows = []
for m in mags:
    s, l = sed.get(m, {}), liq.get(m, {})
    sb, lb = br(sed, m), br(liq, m)
    inS = sb >= BP; inL = lb >= BP
    call = ('BOTH' if inS and inL else
            'sediment only' if inS else
            'culture only' if inL else 'neither')
    print(f"  {mid(m):<17}{gen(m)[:18]:<19}{origin(m):<5}"
          f"{s.get('depth',0):>10.2f}{sb:>10.2f}{s.get('relab',0):>7.2f}"
          f"{l.get('depth',0):>10.2f}{lb:>10.2f}{l.get('relab',0):>7.2f}  {call}")
    rows.append(dict(mag_id=mid(m), genus=gen(m), phylum=phy(m), recovered_from=origin(m),
                     sed_depth=f"{s.get('depth',0):.3f}", sed_breadth=f"{sb:.4f}",
                     sed_relabund=f"{s.get('relab',0):.4f}",
                     liq_depth=f"{l.get('depth',0):.3f}", liq_breadth=f"{lb:.4f}",
                     liq_relabund=f"{l.get('relab',0):.4f}", call=call))
with open(os.path.join(X, 'cross_abundance.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

# ---------------------------------------------------------------- CONTROLS
print()
print("=" * 112)
print("  2. CONTROLS — SELF-MAPPING")
print("=" * 112)
fs = [m for m in mags if origin(m) == 'SED' and br(sed, m) < BP]
fl = [m for m in mags if origin(m) == 'LIQ' and br(liq, m) < BP]
OK = True
if fs:
    OK = False
    print(f"  [FAIL] {len(fs)} sediment MAG(s) not detected in the sediment they came")
    print(f"         from: {', '.join(mid(x) for x in fs)}")
elif [m for m in mags if origin(m) == 'SED']:
    print(f"  [PASS] all {sum(1 for m in mags if origin(m)=='SED')} sediment MAGs "
          f"recovered from sediment reads.")
if fl:
    OK = False
    print(f"  [FAIL] {len(fl)} culture MAG(s) not detected in the culture reads: "
          f"{', '.join(mid(x) for x in fl)}")
else:
    print(f"  [PASS] all {sum(1 for m in mags if origin(m)=='LIQ')} culture MAGs "
          f"recovered from culture reads.")
if not OK:
    print()
    print("  >> A genome must be found in the sample it was assembled from. Until")
    print("     that holds, NOTHING off the diagonal can be read. Stop here.")

# ------------------------------------------------- THE QUESTION THIS ANSWERS
print()
print("=" * 112)
print("  3. DID THE ENRICHMENTS PULL REAL RESIDENTS?")
print("=" * 112)
cult = [m for m in mags if origin(m) == 'LIQ']
seen = [m for m in cult if br(sed, m) >= BP]
strong = [m for m in cult if br(sed, m) >= BS]
if not OK:
    print("  >> NO VERDICT — controls failed.")
elif not cult:
    print("  no culture MAGs to assess.")
else:
    print(f"  {len(seen)} of {len(cult)} culture genomes are detectable in the "
          f"sediment ({len(strong)} strongly).")
    # THE THRESHOLD IS MINE; THE GAP IS THE DATA'S. Rather than defend 0.30,
    # find the largest jump in the sorted breadths and say whether the cut sits
    # in it. A threshold that lands inside a real gap is robust to being moved;
    # one that slices through a continuum is not, and the reader should be told
    # which of those they are looking at.
    bs = sorted((br(sed, m) for m in cult), reverse=True)
    gaps = [(bs[i] - bs[i+1], bs[i+1], bs[i]) for i in range(len(bs) - 1)]
    g, lo, hi = max(gaps) if gaps else (0, 0, 0)
    print(f"  Sorted sediment breadths: {', '.join(f'{b:.2f}' for b in bs)}")
    print(f"  Largest gap: {g:.2f}, between {lo:.2f} and {hi:.2f}.")
    if lo <= BP <= hi:
        print(f"  The {BP:.2f} cut falls INSIDE that gap, so the split is a feature of")
        print("  the data and moving the threshold within the gap changes nothing.")
    else:
        print(f"  The {BP:.2f} cut does NOT fall in that gap — it slices a continuum.")
        print(f"  Quote the breadth for each genome rather than the count above; the")
        print(f"  natural break in these data is at {lo:.2f}-{hi:.2f}, not at {BP:.2f}.")
    print()
    for m in sorted(cult, key=lambda z: -br(sed, z)):
        b = br(sed, m)
        mark = 'STRONG' if b >= BS else ('present' if b >= BP else 'below threshold')
        print(f"    {mid(m):<17}{gen(m)[:20]:<21}sediment breadth {b:>6.2f}   {mark}")
    print()
    frac = len(seen) / len(cult)
    if frac >= 0.6:
        print("  >> THE ENRICHMENTS SELECTED FROM THE RESIDENT COMMUNITY. Most of")
        print("     what grew was already there in detectable amounts, so the")
        print("     culture genomes can be used to interpret the sediment. That is")
        print("     the permission this chapter needs to link the two halves.")
    elif frac >= 0.25:
        print("  >> MIXED. A substantial minority of the cultured organisms are not")
        print("     detectable in situ. They are not necessarily absent — the")
        print("     sediment catalogue holds ~2% of the community and detection")
        print("     needs coverage — but each one must be argued individually, not")
        print("     assumed to represent the lagoon.")
    else:
        print("  >> THE ENRICHMENTS DID NOT RECOVER THE ABUNDANT COMMUNITY. Most")
        print("     culture genomes are undetectable in the sediment. They are real")
        print("     organisms that grew from that inoculum, but they were rare, and")
        print("     no statement about in-situ function may rest on them.")
    print()
    print("  NOTE ON ABSENCE. Sediment sequencing depth sets a detection floor. A")
    print("  genome below threshold is 'not detected at this depth', never 'absent")
    print("  from the lagoon'. Quote the floor with the claim.")

print()
print("=" * 112)
print("  4. THE CONVERSE — sediment MAGs in the enrichments")
print("=" * 112)
# THE ENRICHMENT FACTOR MUST COME FROM RELATIVE ABUNDANCE, NOT RAW DEPTH.
# The first version divided culture depth by sediment depth. The two libraries
# are not the same size — 6.83M sediment reads against 2.82M culture reads —
# so a depth ratio silently carries the library-size ratio inside it and
# UNDERSTATES enrichment by roughly 2.4x. Relative abundance is already
# normalised to its own library, which is the whole point of it.
print("  Enrichment = culture relative abundance / sediment relative abundance.")
print("  NOT a depth ratio: the two libraries differ in size, and a depth ratio")
print("  would carry that difference inside it.")
print()
for m in [z for z in mags if origin(z) == 'SED']:
    sb, lb = br(sed, m), br(liq, m)
    sr = sed.get(m, {}).get('relab', 0.0); lr = liq.get(m, {}).get('relab', 0.0)
    ef = (lr / sr) if sr else float('inf')
    efs = '   inf' if ef == float('inf') else f"{ef:>6.1f}"
    tag = ('enriched' if ef > 2 else
           'unchanged' if ef > 0.5 else 'DEPLETED by the enrichment')
    print(f"  {mid(m):<17}{gen(m)[:20]:<21}sed {sb:>5.2f} ({sr:>5.2f}%)  "
          f"cult {lb:>5.2f} ({lr:>5.2f}%)  x{efs}  {tag}")
print()
print("  ETO_SED_MAG_002 and ETO_LIQ_MAG_012 are the same Chakrabartyella strain")
print("  at 99.8% ANI, so this row is the cleanest measurement of what enrichment")
print("  did to one organism, with the organism held constant. inStrain (step 32)")
print("  would take it further and ask whether the cultured population is a")
print("  SUBSET of the sediment population or a diverged one.")
PY

echo
echo "######## VALIDATION GATE ########"
check "cross abundance table" $([ -s "$X/cross_abundance.tsv" ] && echo 0 || echo 1)
check "SED mapping"           $([ -s "$X/SED_vs_all.tsv" ] && echo 0 || echo 1)
check "LIQ mapping"           $([ -s "$X/LIQ_vs_all.tsv" ] && echo 0 || echo 1)
runlog "26_crossmap complete"
echo
echo "######## 26 DONE ########"
