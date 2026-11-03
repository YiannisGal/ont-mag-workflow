#!/usr/bin/env bash
# =============================================================================
# 32_strain_divergence.sh  --  what did the enrichment do to one organism?
#
# THE QUESTION
#   ETO_SED_MAG_002 (sediment) and ETO_LIQ_MAG_012 (culture) are the same
#   Chakrabartyella strain at 99.80% ANI — the only conspecific pair in a
#   catalogue of 29 MAGs / 28 species. Step 26 confirmed both are present in
#   BOTH samples (sed breadth 0.99, cult 0.98). That makes this the one place in
#   the dataset where "what does enrichment select for" can be asked with the
#   ORGANISM HELD CONSTANT.
#
#   Is the cultured population a SUBSET of the sediment population — the same
#   strains, just fewer — or has it diverged?
#
# WHY NOT inStrain
#   inStrain is the standard tool and it is built for ILLUMINA. Its SNV model
#   assumes per-base error around 0.1%; ONT R10.4.1 with Dorado sup is nearer
#   1%. Run unadjusted it would call thousands of false variants from sequencing
#   error and popANI would be meaningless. Rather than force a short-read tool
#   onto long reads and quietly hope, this does the comparison natively with an
#   EMPIRICAL error floor.
#
# THE CONTROL, AND IT IS THE WHOLE DESIGN
#   Split the sediment reads in half at random and run the two halves through
#   the identical pipeline. Any "difference" between two halves of the SAME
#   library is pure sequencing and mapping error. That is the floor.
#
#       SED-half-A vs SED-half-B  ->  the error floor, measured not assumed
#       SED (all)   vs LIQ (all)  ->  the measurement
#
#   The comparison means something only if it clears the floor. If it does not,
#   the answer is "not distinguishable at this depth", which is a real answer.
#
# WHAT COUNTS AS A DIFFERENCE
#   A site is a FIXED DIFFERENCE only when both samples have depth >= MINDP, the
#   major allele differs, AND each major allele is at >= MAJFRAC of its sample's
#   reads. Requiring near-fixation on both sides is what keeps ONT error out: a
#   1% error does not reach 80% frequency, and a site that is 50/50 in one
#   sample is polymorphism, not divergence.
#
# OUTPUT  32_strain/<ref>.divergence.tsv   per-comparison summary
#         32_strain/<ref>.fixed_sites.tsv  the sites themselves
# RUNTIME 20-40 min: four mappings against one small genome.
# USAGE   ws2 && bash scripts/32_strain_divergence.sh 2>&1 | tee logs/32.log
#         REF=ETO_SED_MAG_002 to choose a different reference genome
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="32_strain"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; S="$WS2/32_strain"; mkdir -p "$S"
rm -f "$S"/*.divergence.tsv

REFID="${REF:-ETO_SED_MAG_002}"
MINDP="${MINDP:-20}"          # minimum depth in BOTH samples at a site
MAJFRAC="${MAJFRAC:-0.80}"    # major allele must reach this fraction in both
SEED="${SEED:-42}"

echo "######## 32 — STRAIN DIVERGENCE, SEDIMENT vs ENRICHMENT ########"
echo "  started: $(date -Is)"
runlog "32_strain started"
echo "  reference: $REFID   min depth $MINDP   major-allele fraction $MAJFRAC"

# ---- find the reference genome fasta ---------------------------------------
REFFA=""
WN=$(awk -F'\t' -v id="$REFID" 'NR==1{for(i=1;i<=NF;i++){if($i=="mag_id")a=i;if($i=="working_name")b=i};next}
     $a==id{print $b}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
for cand in "$F/bins/${REFID}.fa" "$F/bins/${WN}.fa" \
            $(find "$F" -maxdepth 3 \( -name "${REFID}.f*a" -o -name "${WN}.f*a" \) 2>/dev/null); do
  [ -n "$cand" ] && [ -s "$cand" ] && { REFFA="$cand"; break; }
done
[ -n "$REFFA" ] || { echo "  !! could not find a fasta for $REFID (working name '${WN:-?}')"; exit 1; }
echo "  fasta: ${REFFA#$WS2/}  ($(grep -c '^>' "$REFFA") contigs)"

# ---- reads ------------------------------------------------------------------
CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
if [ -s "$CLEAN" ]; then SED_IN="$CLEAN"
else SED_IN="$SED_READS"; echo "  !! using contaminated \$SED_READS — run 13f WRITE_CLEAN=1"; fi

# SPLIT-HALF CONTROL. Deterministic, seeded, and built from the SAME library, so
# every difference it produces is error rather than biology.
HA="$S/sedA.fq.gz"; HB="$S/sedB.fq.gz"
if [ ! -s "$HA" ] || [ ! -s "$HB" ]; then
  echo "  building split-half control (seed $SEED) ..."
  zcat -f "$SED_IN" | awk -v s="$SEED" 'BEGIN{srand(s)}
    NR%4==1{keep=(rand()<0.5)}
    {if(keep) print > "/dev/stdout"; else print > "/dev/stderr"}' \
    2> >(gzip > "$HB") | gzip > "$HA"
  echo "    A: $(zcat "$HA" | awk 'NR%4==1' | wc -l) reads   B: $(zcat "$HB" | awk 'NR%4==1' | wc -l) reads"
fi

# ---- map --------------------------------------------------------------------
domap () {
  local tag="$1"
  local reads="$2"
  local bam="$S/${tag}.bam"
  if [ -s "$bam" ]; then echo "  [$tag] cached"; return 0; fi
  echo "  [$tag] mapping ..."
  run_tool minimap2 -ax map-ont -t "$THREADS" --secondary=no "$REFFA" "$reads" 2>/dev/null \
    | run_tool samtools sort -@ 4 -o "$bam.part" - 2>/dev/null
  local rc=$?
  [ "$rc" -eq 0 ] && [ -s "$bam.part" ] || { rm -f "$bam.part"; echo "  !! $tag mapping failed"; return 1; }
  mv "$bam.part" "$bam"; run_tool samtools index "$bam" 2>/dev/null
}
echo
echo "######## MAPPING ########"
domap SED  "$SED_IN"     || exit 1
domap LIQ  "$LIQ_READS"  || exit 1
domap sedA "$HA"         || exit 1
domap sedB "$HB"         || exit 1

# ---- pileup -----------------------------------------------------------------
pile () {
  local tag="$1"
  local out="$S/${tag}.counts.tsv"
  if [ -s "$out" ]; then return 0; fi
  # -aa gives every position including zero-depth, so the two samples stay in
  # register without an interval join.
  run_tool samtools mpileup -aa -A -Q 5 -d 5000 -f "$REFFA" "$S/${tag}.bam" 2>/dev/null \
    | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > "$out.part" && mv "$out.part" "$out"
}
echo
echo "######## PILEUP ########"
for t in SED LIQ sedA sedB; do echo "  [$t]"; pile "$t"; done

python3 - "$S" "$MINDP" "$MAJFRAC" <<'PY'
import sys, os, re
from collections import Counter
S, MINDP, MAJFRAC = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])

CLEAN = re.compile(r'\^.|\$|[+-](\d+)')
def parse(bases, ref):
    """samtools mpileup base string -> Counter of A/C/G/T. Indel runs are
    stripped by length, read-start/end markers dropped; . and , are the
    reference base on either strand."""
    out, i, n = Counter(), 0, len(bases)
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            j = i + 1; num = ''
            while j < n and bases[j].isdigit(): num += bases[j]; j += 1
            i = j + int(num or 0); continue
        if c in '.,': out[ref.upper()] += 1
        elif c.upper() in 'ACGT': out[c.upper()] += 1
        i += 1
    return out

def load(tag):
    p = os.path.join(S, f'{tag}.counts.tsv')
    d = {}
    if not os.path.exists(p): return d
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5: continue
        try: dp = int(f[3])
        except ValueError: continue
        if dp < MINDP: continue
        c = parse(f[4], f[2])
        tot = sum(c.values())
        if tot < MINDP: continue
        b, n = c.most_common(1)[0]
        d[(f[0], f[1])] = (b, n / tot, tot)
    return d

sets = {t: load(t) for t in ('SED', 'LIQ', 'sedA', 'sedB')}
for t, d in sets.items():
    print(f"  {t}: {len(d)} sites at depth >= {MINDP}")

def compare(a, b):
    A, B = sets[a], sets[b]
    shared = A.keys() & B.keys()
    fixed = [k for k in shared
             if A[k][0] != B[k][0] and A[k][1] >= MAJFRAC and B[k][1] >= MAJFRAC]
    return len(shared), fixed

print()
print("=" * 96)
print("  FIXED DIFFERENCES")
print("=" * 96)
print(f"  A site counts only if BOTH samples have depth >= {MINDP} and each major")
print(f"  allele reaches {MAJFRAC:.0%}. Near-fixation on both sides is what keeps ONT")
print("  error out: a 1% error rate does not produce an 80% allele.")
print()
res = {}
for lbl, (a, b) in (('sedA vs sedB  (ERROR FLOOR)', ('sedA', 'sedB')),
                    ('SED  vs LIQ   (measurement)', ('SED', 'LIQ'))):
    ns, fx = compare(a, b)
    per_kb = 1000 * len(fx) / ns if ns else 0
    res[lbl] = (ns, len(fx), per_kb)
    print(f"  {lbl:<30}{ns:>9} sites compared{len(fx):>8} fixed diffs"
          f"{per_kb:>10.3f} per kb")
    if a == 'SED':
        with open(os.path.join(S, 'fixed_sites.tsv'), 'w') as fh:
            fh.write("contig\tpos\tSED_allele\tSED_frac\tSED_dp\t"
                     "LIQ_allele\tLIQ_frac\tLIQ_dp\n")
            for k in sorted(fx):
                sa, sf, sd = sets['SED'][k]; la, lf, ld = sets['LIQ'][k]
                fh.write(f"{k[0]}\t{k[1]}\t{sa}\t{sf:.2f}\t{sd}\t{la}\t{lf:.2f}\t{ld}\n")

floor = res['sedA vs sedB  (ERROR FLOOR)'][2]
meas  = res['SED  vs LIQ   (measurement)'][2]
print()
print("=" * 96)
print("  VERDICT")
print("=" * 96)
if res['sedA vs sedB  (ERROR FLOOR)'][0] == 0 or res['SED  vs LIQ   (measurement)'][0] == 0:
    print("  >> NO VERDICT. Too few sites reached the depth threshold in one of the")
    print(f"     comparisons. Lower MINDP (currently {MINDP}) or accept that this")
    print("     genome is not covered deeply enough in both samples to compare.")
elif floor > 0 and meas < 3 * floor:
    print(f"  >> NOT DISTINGUISHABLE. {meas:.3f} fixed differences per kb between")
    print(f"     sediment and culture, against an error floor of {floor:.3f} per kb")
    print("     measured on two halves of the SAME library. The enrichment did not")
    print("     detectably change this population: what grew in the bottle is the")
    print("     sediment population, not a diverged descendant of it.")
    print("     That is a real result and the conservative one — it licenses using")
    print("     the culture genome to interpret the sediment organism.")
else:
    print(f"  >> DIVERGENCE DETECTED. {meas:.3f} fixed differences per kb, against")
    print(f"     an error floor of {floor:.3f} per kb from the split-half control —")
    print(f"     {meas/max(floor,1e-9):.1f}x the floor.")
    print("     The cultured population is not simply a subset of the sediment one.")
    print("     Inspect fixed_sites.tsv: if the differences cluster in a few genes")
    print("     that is selection on those genes, and if they are spread evenly it")
    print("     is more likely a different strain having been enriched from a")
    print("     mixed in-situ population.")
print()
print("  LIMITS. Fixed differences are not a phylogeny and this is one genome.")
print("  A zero here means 'no difference resolvable above the measured error",
      "floor', never 'the populations are identical'.")
PY

echo
echo "######## VALIDATION GATE ########"
check "SED bam"       $([ -s "$S/SED.bam" ] && echo 0 || echo 1)
check "LIQ bam"       $([ -s "$S/LIQ.bam" ] && echo 0 || echo 1)
check "split-half A"  $([ -s "$S/sedA.bam" ] && echo 0 || echo 1)
check "split-half B"  $([ -s "$S/sedB.bam" ] && echo 0 || echo 1)
runlog "32_strain complete"
echo
echo "######## 32 DONE ########"
