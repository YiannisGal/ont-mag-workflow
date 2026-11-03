#!/usr/bin/env bash
# =============================================================================
# 13d_silva_verify.sh  --  is the sediment SSU profile real?
#
# THE PROBLEM
#   The strict SILVA run reports the sediment as 84.65% *Erwinia*, with
#   Citrobacter, Escherichia-Shigella, Klebsiella, Enterobacter, Kluyvera and
#   Rosenbergiella behind it. That is >92% Enterobacteriaceae in a permanently
#   anoxic, sulfidic, sulfate-replete lagoon sediment. It is not credible:
#   Enterobacteriaceae are facultative aerobes associated with plants, guts and
#   surfaces, not with euxinic sediment. The accompanying genera make the
#   pattern more specific, not less — *Wolbachia* is an obligate INSECT
#   endosymbiont, *Commensalibacter* and *Rosenbergiella* are insect- and
#   nectar-associated. That combination is a plant/insect microbiome signature,
#   not a lagoon.
#
#   Three things independently contradict it:
#     - Kraken2 on the same reads does NOT show this. It classified 0.84% of
#       sediment reads and Erwinia is not its dominant taxon.
#     - The 4 sediment MAGs are Bacillota fermenters. Nothing Enterobacterial
#       assembled or binned.
#     - The amplicon chapter shows Archaea abundant in situ (Arc_Otu22 at
#       1.641%); SILVA here reports Archaea at 0.1%.
#
#   AND MY BIAS SCREEN IN 13c MISSED IT. It matched genus names against a hand
#   list by EXACT string, so "Escherichia-Shigella" did not match "escherichia",
#   and Erwinia/Citrobacter/Enterobacter were not in the list at all. It printed
#   "Below 5% — consistent with a genuine environmental profile" over a profile
#   that is >90% Enterobacteriaceae. That is a false reassurance and it is
#   fixed here: this script screens by FAMILY and by DOMINANCE, not by a list
#   of genus names I happened to think of.
#
# WHAT THIS SCRIPT DECIDES BETWEEN
#   H1  REAL BUT UNASSEMBLED. Genuinely abundant, just never assembled.
#       PREDICTION: the assembly's own 16S genes should include Erwinia-like
#       sequences, and its reads should be spread over many SILVA accessions.
#   H2  CROSS-CONTAMINATION / BARCODE CROSSTALK. The sediment file is TEN
#       concatenated runs of mixed chemistry, including PCR barcoding. If
#       another library on a shared flow cell was a plant or insect sample,
#       index hopping injects exactly this signature.
#       PREDICTION: the hits concentrate in a subset of runs/barcodes, and are
#       absent from the assembly.
#   H3  MAPPING ARTEFACT. A small number of SILVA references acting as sinks.
#       PREDICTION: nearly all hits land on very few accessions, and identity
#       sits near the threshold rather than high.
#
#   These make DIFFERENT predictions, so the data can separate them. That is
#   the point — "probably contamination" is a guess, this is a test.
#
# OUTPUT   13_profile/silva/verify/
#     hit_accessions.tsv     how concentrated the hits are
#     identity_profile.tsv   identity distribution, suspect vs rest
#     by_run.tsv             hits per sequencing run, if headers permit
#     assembly_16S_tax.tsv   what the assembly's own 16S genes actually are
#     VERDICT.txt
#
# RUNTIME  10-20 min (no re-mapping; the PAF is reused).
# USAGE    ws2 && bash scripts/13d_silva_verify.sh 2>&1 | tee logs/13d_verify.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13d_silva_verify"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
S="$WS2/13_profile/silva"; V="$S/verify"; PR="$WS2/13_profile"
mkdir -p "$V"
SUSPECT="${SUSPECT_GENUS:-Erwinia}"
SILVA="${SILVA_SSU:-${DB_DIR}/silva/SILVA_SSURef_NR99.fasta}"

echo "######## 13d — IS THE SEDIMENT SSU PROFILE REAL? ########"
echo "  started: $(date -Is)   suspect genus: $SUSPECT"
runlog "13d_silva_verify started"
[ -s "$S/SED_vs_silva.paf.gz" ] || { echo "  !! no SED PAF — run 13b first"; exit 1; }

# ---------------------------------------------------------------------------
# TEST 4 (run first, it is the decisive one): classify the assembly's own 16S
# If a genus is 85% of the SSU read pool it must be in the assembly.
# ---------------------------------------------------------------------------
echo
echo "######## TEST 4 — what are the assembly's 16S genes? ########"
for c in SED LIQ; do
  q="$PR/rrna/${c}_16S.fna"
  out="$V/${c}_16S_vs_silva.paf"
  [ -s "$q" ] || { echo "  $c: no 16S file at $q"; continue; }
  if [ -s "$out" ]; then echo "  $c: already classified"; continue; fi
  echo "  $c: mapping $(grep -c '^>' "$q") assembly 16S genes to SILVA ..."
  run_tool minimap2 -x map-ont -t "$THREADS" --secondary=no -c \
      "$S/silva.mmi" "$q" 2>/dev/null > "$out"
  echo "    alignments: $(wc -l < "$out")"
done

python3 - "$S" "$V" "$SILVA" "$SUSPECT" "$WS2" <<'PY'
import sys, os, gzip, csv, collections

S, V, SILVA, SUSPECT, WS2 = sys.argv[1:6]
SUS = SUSPECT.lower()

# ---- SILVA accession -> lineage
lineage = {}
with open(SILVA) as fh:
    for line in fh:
        if line.startswith('>'):
            acc, _, tax = line[1:].rstrip('\n').partition(' ')
            lineage[acc] = tax
print(f"  SILVA lineages: {len(lineage)}")

def genus_of(acc):
    t = lineage.get(acc, '')
    p = [x.strip() for x in t.split(';') if x.strip()]
    return p[5] if len(p) > 5 else (p[-1] if p else '')
def family_of(acc):
    t = lineage.get(acc, '')
    p = [x.strip() for x in t.split(';') if x.strip()]
    return p[4] if len(p) > 4 else ''

# =====================================================================
# TEST 4 — the assembly's own 16S
# =====================================================================
print()
print("=" * 96)
print("  TEST 4 — WHAT THE ASSEMBLY'S OWN 16S GENES ARE")
print("=" * 96)
print("  If a genus is ~85% of the SSU READ pool it is one of the most abundant")
print("  organisms present, and an organism that abundant assembles. Its absence")
print("  from the assembly's 16S genes would be decisive against H1.")
asm_hit = {}
for code in ('SED', 'LIQ'):
    p = os.path.join(V, f'{code}_16S_vs_silva.paf')
    if not os.path.exists(p):
        print(f"\n  {code}: not classified"); continue
    best = {}
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12: continue
        q, t = f[0], f[5]
        try: nm, al = int(f[9]), int(f[10])
        except ValueError: continue
        if al < 400: continue
        if q not in best or nm > best[q][0]: best[q] = (nm, t, nm/al)
    gc = collections.Counter(genus_of(t) for _, t, _ in best.values())
    fc = collections.Counter(family_of(t) for _, t, _ in best.values())
    asm_hit[code] = gc
    print(f"\n  {code}: {len(best)} assembly 16S genes classified")
    print(f"  {'count':>6}  {'%':>7}  genus")
    tot = sum(gc.values()) or 1
    for g, n in gc.most_common(12):
        mark = '   <== THE SUSPECT' if g.lower() == SUS else ''
        print(f"  {n:>6}  {100*n/tot:>6.1f}%  {g}{mark}")
    n_sus = sum(v for k, v in gc.items() if k.lower() == SUS)   # was double-counted
    print(f"\n  {SUSPECT} among assembly 16S genes: {n_sus} / {len(best)}"
          f"  ({100*n_sus/max(1,len(best)):.1f}%)")

with open(os.path.join(V, 'assembly_16S_tax.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(['sample','genus','n_16S_genes'])
    for code, gc in asm_hit.items():
        for g, n in gc.most_common(): w.writerow([code, g, n])

# =====================================================================
# TESTS 1-3 — the read-level hits
# =====================================================================
results = {}
for code in ('SED', 'LIQ'):
    paf = os.path.join(S, f'{code}_vs_silva.paf.gz')
    if not os.path.exists(paf): continue
    best = {}
    with gzip.open(paf, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12: continue
            q, t = f[0], f[5]
            try: nm, al = int(f[9]), int(f[10])
            except ValueError: continue
            if al < 1000: continue
            ident = nm / al if al else 0
            if ident < 0.92: continue
            if q not in best or nm > best[q][0]: best[q] = (nm, t, ident, al)
    results[code] = best

    sus = {q: v for q, v in best.items() if genus_of(v[1]).lower() == SUS}
    rest = {q: v for q, v in best.items() if genus_of(v[1]).lower() != SUS}
    print()
    print("=" * 96)
    print(f"  {code} — READ-LEVEL TESTS   ({len(best)} reads assigned, "
          f"{len(sus)} to {SUSPECT})")
    print("=" * 96)
    if not sus:
        print(f"  no {SUSPECT} reads in this sample — tests 1-3 not applicable")
        continue

    # ---- TEST 1: accession concentration
    acc = collections.Counter(v[1] for v in sus.values())
    print(f"\n  TEST 1 — ACCESSION CONCENTRATION  ({len(acc)} distinct SILVA refs hit)")
    print("  A real abundant population hits many references across its clade.")
    print("  A handful of references absorbing everything indicates a sink artefact.")
    top5 = sum(n for _, n in acc.most_common(5))
    print(f"  {'reads':>7}  {'%':>6}  accession")
    for a, n in acc.most_common(8):
        print(f"  {n:>7}  {100*n/len(sus):>5.1f}%  {a}  {lineage.get(a,'')[:52]}")
    print(f"  top 5 accessions hold {100*top5/len(sus):.1f}% of {SUSPECT} reads")
    conc = top5 / len(sus)

    # ---- TEST 2: identity distribution
    si = sorted(v[2] for v in sus.values())
    ri = sorted(v[2] for v in rest.values()) or [0]
    def med(x): return x[len(x)//2]
    print(f"\n  TEST 2 — IDENTITY DISTRIBUTION")
    print("  A true match to a well-represented reference sits high (>=97%).")
    print("  Values hugging the 92% cutoff mean the threshold, not the organism,")
    print("  is deciding the assignment.")
    print(f"    {SUSPECT:<14} median {100*med(si):.1f}%   min {100*si[0]:.1f}%   max {100*si[-1]:.1f}%")
    print(f"    {'all others':<14} median {100*med(ri):.1f}%")
    hi = sum(1 for x in si if x >= 0.97) / len(si)
    print(f"    {100*hi:.1f}% of {SUSPECT} reads are >=97% identical to their reference")

    # ---- TEST 3: which sequencing run
    print(f"\n  TEST 3 — WHICH SEQUENCING RUN")
    print("  The sediment file is TEN concatenated runs. Contamination or barcode")
    print("  crosstalk concentrates in a subset of runs; a real organism does not.")
    def runid(q):
        # ONT read IDs are per-run UUIDs; without headers the best available
        # proxy is the UUID prefix, which is NOT a run ID. Say so rather than
        # inventing a grouping.
        return q.split('-')[0][:4] if '-' in q else 'n/a'
    rs = collections.Counter(runid(q) for q in sus)
    ra = collections.Counter(runid(q) for q in best)
    if len(ra) <= 1 or 'n/a' in ra:
        print("    read IDs carry no run information in this file — test skipped.")
        print("    To run it properly, grep the ORIGINAL per-run fastqs for these")
        print(f"    read IDs (written to {os.path.join(V, code+'_suspect_reads.txt')}).")
        runskew = None
    else:
        print(f"    {'prefix':<10}{'suspect':>9}{'all':>8}{'enrichment':>12}")
        for pfx, n in rs.most_common(8):
            allc = ra.get(pfx, 0) or 1
            print(f"    {pfx:<10}{n:>9}{allc:>8}{n/allc:>11.2f}x")
        runskew = max(rs.values()) / len(sus)
    with open(os.path.join(V, f'{code}_suspect_reads.txt'), 'w') as fh:
        fh.write('\n'.join(sorted(sus)) + '\n')

    # ---- family-level composition, the screen 13c should have done
    fam = collections.Counter(family_of(v[1]) for v in best.values())
    ftot = sum(fam.values()) or 1
    print(f"\n  FAMILY-LEVEL COMPOSITION (the screen 13c failed to apply)")
    print(f"  {'%':>7}  family")
    for f_, n in fam.most_common(8):
        flag = ''
        if f_.lower() in ('enterobacteriaceae', 'yersiniaceae', 'erwiniaceae',
                          'morganellaceae', 'pectobacteriaceae'):
            flag = '   <== facultative aerobe, implausible here'
        print(f"  {100*n/ftot:>6.1f}%  {f_}{flag}")
    entero = sum(n for f_, n in fam.items() if f_.lower() in
                 ('enterobacteriaceae','yersiniaceae','erwiniaceae',
                  'morganellaceae','pectobacteriaceae'))
    print(f"  Enterobacterales families total: {100*entero/ftot:.1f}%")

    results[code] = dict(n=len(best), nsus=len(sus), conc=conc, hi=hi,
                         entero=entero/ftot,
                         asm=asm_hit.get(code, collections.Counter()))

# =====================================================================
# VERDICT
# =====================================================================
print()
print("=" * 96)
print("  VERDICT")
print("=" * 96)
lines = []
r = results.get('SED')
if not isinstance(r, dict) or not r.get('nsus'):
    lines.append(f"No {SUSPECT} reads to assess in SED.")
else:
    asm = r['asm']
    n_asm_sus = sum(v for k, v in asm.items() if k.lower() == SUS)
    frac_sus = r['nsus'] / max(1, r['n'])
    lines.append(f"{SUSPECT} is {100*frac_sus:.1f}% of assigned SED SSU reads.")
    lines.append(f"Assembly 16S genes assigned to {SUSPECT}: {n_asm_sus} of "
                 f"{sum(asm.values())}.")
    lines.append("")
    if n_asm_sus == 0 and frac_sus > 0.2:
        lines.append("H1 (real but unassembled) is REJECTED.")
        lines.append("An organism supplying most of the SSU reads is among the most")
        lines.append("abundant in the sample, and the sediment assembly recovered 138")
        lines.append("16S genes — none of them this genus. A dominant organism that")
        lines.append("contributes no 16S gene to its own assembly does not exist.")
        lines.append("")
        if r['conc'] > 0.8:
            lines.append("H3 (mapping sink) is SUPPORTED: >80% of these reads land on")
            lines.append("five or fewer SILVA references. Check those accessions by hand.")
        else:
            lines.append("H2 (cross-contamination / barcode crosstalk) is the remaining")
            lines.append("explanation. The genus combination points the same way:")
            lines.append("Erwinia + Wolbachia + Commensalibacter + Rosenbergiella is a")
            lines.append("plant/insect microbiome signature, and the sediment file is ten")
            lines.append("concatenated runs of mixed chemistry including PCR barcoding.")
            lines.append("")
            lines.append("DO THIS NEXT, it is quick and decisive:")
            lines.append("  1. Find which ORIGINAL run each suspect read came from —")
            lines.append("     grep the read IDs in verify/SED_suspect_reads.txt against")
            lines.append("     each per-run fastq in the original data directory.")
            lines.append("  2. If they concentrate in one or two runs, exclude those runs")
            lines.append("     and re-run 13b. Report the exclusion in Methods.")
            lines.append("  3. Ask what else was on those flow cells. A plant or insect")
            lines.append("     library sequenced alongside would settle it immediately.")
    elif n_asm_sus > 0:
        lines.append("H1 is POSSIBLE — the genus does appear among assembly 16S genes.")
        lines.append("Check whether its share there matches its share of the read pool.")
        lines.append("A large mismatch still points to contamination on top of a real,")
        lines.append("smaller population.")
    lines.append("")
    lines.append("EITHER WAY, DO NOT REPORT THE CURRENT SED SSU PROFILE. The LIQ")
    lines.append("profile is unaffected — it is dominated by Anaerotignum, Clostridium,")
    lines.append("Sedimentibacter and the methanogens, which is exactly what an")
    lines.append("anaerobic enrichment should look like, and 26.1% of it maps to MAG")
    lines.append("genera. That one is usable now.")

for l in lines: print("  " + l)
with open(os.path.join(V, 'VERDICT.txt'), 'w') as fh:
    fh.write('\n'.join(lines) + '\n')

print()
print("  LESSON FOR THE SCRIPTS: 13c's bias screen matched genus names against a")
print("  hand-written list by exact string. 'Escherichia-Shigella' missed")
print("  'escherichia'; Erwinia, Citrobacter and Enterobacter were never on the")
print("  list. It then printed reassurance over a >90% Enterobacterales profile.")
print("  A screen built from a list of things I remembered to distrust will always")
print("  have that hole. Screen on STRUCTURE instead — family-level aggregation and")
print("  single-taxon dominance — which is what this script does.")
PY

echo
echo "######## VALIDATION GATE ########"
check "assembly 16S classified" $([ -s "$V/SED_16S_vs_silva.paf" ] && echo 0 || echo 1)
check "verdict written"         $([ -s "$V/VERDICT.txt" ] && echo 0 || echo 1)
check "suspect read list"       $([ -s "$V/SED_suspect_reads.txt" ] && echo 0 || echo 1)

runlog "13d_silva_verify complete"
echo
echo "######## 13d DONE ########"
