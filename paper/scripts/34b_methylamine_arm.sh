#!/usr/bin/env bash
# =============================================================================
# 34b_methylamine_arm.sh  --  the methylamine arm, done with the right KOs and
#                             with the pyrrolysine problem handled
#
# WHY THIS EXISTS — I WAS WRONG, AND THE ko_list GREP SHOWS IT
#   Step 34 tried to resolve mtmB with an AUTO keyword search for
#   'monomethylamine;methyltransferase'. It found nothing and aborted, and I
#   read that as "KEGG has no KO for mtmB". That was wrong. KEGG has one; it
#   simply does not use the word *monomethylamine* in its name:
#
#       K16176  methylamine---corrinoid protein Co-methyltransferase [2.1.1.248]
#
#   EC 2.1.1.248 IS MtmB. The only entry carrying the word *monomethylamine* is
#   K16177, the corrinoid protein MtmC — the promiscuous B12-binding partner
#   that absorbed 1,978 reads in step 34 and had to be thrown out. The guard
#   was right to abort; my reading of why it aborted was not.
#
#   Rule 24 said: when a resolution step replaces your assertion, read what it
#   resolved to. The failure here is the mirror image — when a resolution step
#   finds NOTHING, that is a statement about the search string, not about KEGG.
#   Rule 25, below.
#
# THE THREE METHYLAMINE ROUTES, CORRECTLY NUMBERED
#   monomethylamine   MtmB  K16176   +  MtmC  K16177 (corrinoid)
#   dimethylamine     MtbB  K16178   +  MtbC  K16179 (corrinoid)
#   trimethylamine    MttB  K14083   +  MttC  K14084 (corrinoid)
#   shared second step MtbA K14082   -- methylated corrinoid -> coenzyme M
#
#   methanol          MtaB  K04480   +  MtaC  K14081 (corrinoid)  + MtaA K14080
#
#   The methanol KOs were confirmed against ko_list on 2026-08-27 and step 34's
#   label for K04480 was CORRECT: KEGG calls it 'methanol---5-hydroxybenzimid-
#   azolylcobamide Co-methyltransferase [EC:2.1.1.90]', which is MtaB. I had
#   doubted it; the doubt was wrong. K14080 is MtaA, the second step, not MtaB.
#
#   MtaB IS NOT A PYRROLYSINE PROTEIN. It is a zinc-dependent methyltransferase
#   with the same job in the same organisms as MtmB/MtbB/MttB and no amber
#   codon — which makes it the single best negative control the amber test
#   could have, and it also gives the methanol arm a genome-level test at last.
#
# THE PROBLEM NOBODY'S PIPELINE HANDLES: PYRROLYSINE
#   MtmB, MtbB and MttB each carry an IN-FRAME UAG AMBER CODON that encodes
#   pyrrolysine, the 22nd amino acid. Prodigal, Bakta and every standard caller
#   run translation table 11, where UAG is a stop — so they TRUNCATE these three
#   genes, usually mid-protein. Step 34's level-2 search ran on prodigal
#   proteins. Its assembly-level view of exactly these genes was compromised by
#   construction, and no threshold could have fixed it.
#
#   This step re-translates the candidate contigs READING THROUGH UAG and
#   searches those. The amber codon then stops being an obstacle and becomes
#   the discriminator:
#
# THE AMBER TEST, WITH ITS CONTROL INSIDE THE SAME PATHWAY
#   MtmB / MtbB / MttB       pyrrolysine proteins    -> expect >=1 in-frame UAG
#   MtbA / MtaB / MtaA       same pathway, no Pyl     -> expect 0
#   MtmC / MtbC / MttC/MtaC  corrinoid partners, no Pyl -> expect 0
#
#   MtaB is the sharpest of these: a corrinoid-dependent methyltransferase, in
#   the same operonic neighbourhoods, in the same methanogens, doing the same
#   chemistry on a different substrate — and carrying no pyrrolysine. If the
#   amber separates MtmB/MtbB/MttB from MtaB, it is not separating
#   methyltransferases from non-methyltransferases; it is separating
#   pyrrolysine proteins from everything else, which is the claim.
#
#   If MtbA and the corrinoids show as many in-frame ambers as MtmB does, the
#   test is measuring random UAGs in spurious reading frames and is void. That
#   check runs automatically and blocks the amber verdict.
#
#   This also settles the MttB superfamily trap. Step 34's mttB decoy absorbed
#   312,347 sediment reads because non-pyrrolysine MttB homologues (MtgB and
#   relatives, glycine betaine methyltransferases) are everywhere. A true
#   trimethylamine MttB has the amber codon; MtgB does not.
#
# ORTHOGONAL CORROBORATION FROM A DIFFERENT GENE SET
#   Pyrrolysine has to be made and charged: pylB/pylC/pylD plus pylS/pylT. A
#   genome that encodes methylamine methyltransferases should encode the
#   pyrrolysine machinery, and that machinery shares no homology with the
#   methyltransferases at all. Agreement between two unrelated gene sets is
#   worth more than any score threshold.
#
# WHAT THIS DOES NOT DO
#   It does not redo step 34's read level. Step 34's read numbers stand, with
#   corrected labels. mtmB (K16176) was never searched at read level by any
#   step — that is stated as a gap, not filled with a guess.
#
# OUTPUT  34b_methylamine/markers.tsv · amber_test.tsv · pyl_machinery.tsv
#         VERDICT.txt
# RUNTIME 20-60 min
# USAGE   ws2 && mkdir -p logs && bash scripts/34b_methylamine_arm.sh 2>&1 | tee logs/34b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="34b_methylamine"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"

A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"; U="$WS2/34_unbinned"
D="$WS2/34b_methylamine"; mkdir -p "$D/profiles" "$D/nom" "$D/rt" "$D/hits"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 34b — THE METHYLAMINE ARM, PYRROLYSINE-AWARE ########"
echo "  started: $(date -Is)"
runlog "34b_methylamine started"

# ------------------------------------------------------------------ binaries
bin_env() {
  local tool="$1"; local var="ENV_${tool//-/_}"
  local cands="${!var:-} ${ENV_hmmsearch:-} ${ENV_hmmscan:-} ${ENV_hmmer:-} ${ENV_bakta:-} -"
  local e prefix root d
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }; continue; fi
    if declare -f env_prefix >/dev/null 2>&1; then
      prefix=$(env_prefix "$e" 2>/dev/null) || continue
      [ -x "$prefix/bin/$tool" ] && { echo "$e"; return 0; }
    fi
  done
  for root in ${CONDA_ROOTS:-}; do
    [ -d "$root" ] || continue
    for d in "$root"/*; do [ -x "$d/bin/$tool" ] && { echo "$(basename "$d")"; return 0; }; done
  done
  return 1
}
HS=$(bin_env hmmsearch || true)
[ -n "$HS" ] || { echo "  !! hmmsearch not found"; runlog "**34b ABORTED: no hmmsearch**"; exit 2; }
export ENV_hmmsearch="$HS"; echo "  hmmsearch: env '$HS'"

# ------------------------------------------------------------------- ko_list
for src in "$K/ko_list" "$U/ko_list" "$WS2/36_gapclose/ko_list" "$WS2/19_auxotrophy/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$D/ko_list" 2>/dev/null || true; break; }
done
[ -s "$D/ko_list" ] || { wget -q -O "$D/ko_list.gz" "$BASE/ko_list.gz" 2>/dev/null && gunzip -f "$D/ko_list.gz" 2>/dev/null || true; }
[ -s "$D/ko_list" ] || { echo "  !! no ko_list"; runlog "**34b ABORTED: no ko_list**"; exit 3; }
echo "  ko_list: $(wc -l < "$D/ko_list") families"

# ============================================================ 0. KO GUARD
echo
echo "######## 0. KO IDENTITY GUARD ########"
# Every KO is asserted BY NUMBER and checked against KEGG's own definition.
# No keyword resolution for the markers — that is the step that went wrong.
cat > "$D/targets.tsv" <<'KOS'
K16176	MtmB	PYL	monomethylamine methyltransferase	2.1.1.248
K16178	MtbB	PYL	dimethylamine methyltransferase	2.1.1.249
K14083	MttB	PYL	trimethylamine methyltransferase	2.1.1.250
K14082	MtbA	NONPYL_CTRL	methylated corrinoid:coenzyme M methyltransferase	2.1.1.247
K04480	MtaB	NONPYL_CTRL	methanol methyltransferase	2.1.1.90
K14080	MtaA	NONPYL_CTRL	methanol/betaine corrinoid:coenzyme M methyltransferase	2.1.1.246
K16177	MtmC	CORRINOID_DECOY	monomethylamine corrinoid protein	-
K16179	MtbC	CORRINOID_DECOY	dimethylamine corrinoid protein	-
K14084	MttC	CORRINOID_DECOY	trimethylamine corrinoid protein	-
K14081	MtaC	CORRINOID_DECOY	methanol corrinoid protein	-
KOS
python3 - "$D" <<'PYG'
import sys, os
D = sys.argv[1]
defs = {}
with open(os.path.join(D, 'ko_list'), errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    di = hdr.index('definition') if 'definition' in hdr else len(hdr) - 1
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > di:
            defs[f[0]] = f[di]
bad = 0
print("  asserted by KO NUMBER, verified against KEGG's own definition:")
for line in open(os.path.join(D, 'targets.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5:
        continue
    ko, gene, role, claim, ec = f
    d = defs.get(ko)
    if d is None:
        print(f"  [FAIL] {ko} ({gene}) is not in ko_list at all"); bad += 1; continue
    # the EC number is the identity check; a name can drift, an EC cannot
    if ec != '-' and ec not in d:
        print(f"  [FAIL] {ko} claimed as {gene} (EC {ec}) but KEGG says: {d}")
        bad += 1; continue
    if ec == '-' and 'corrinoid' not in d.lower():
        print(f"  [FAIL] {ko} claimed as corrinoid protein {gene} but KEGG says: {d}")
        bad += 1; continue
    print(f"  [ok]   {ko}  {gene:<5} {role:<16} {d[:52]}")
# Pyrrolysine machinery. The first version searched the single word
# 'pyrrolysin' and found ZERO KOs — then said so and moved on, which is the
# very failure Rule 25 was written about three sections above this one. KEGG
# names pylB/pylC/pylD by their CHEMISTRY, not by the amino acid they build:
#   pylB  methylornithine synthase                    [EC:5.4.99.58]
#   pylC  3-methylornithine---L-lysine ligase         [EC:6.3.2.-]
#   pylD  3-methylornithine-N6-L-lysine dehydrogenase [EC:1.5.1.-]
#   pylS  pyrrolysine---tRNA(Pyl) ligase              [EC:6.1.1.26]
# So the search is over the reaction names as well, and every candidate it
# finds is PRINTED WITH ITS DEFINITION so the resolution can be read rather
# than trusted.
TERMS = ('pyrrolysin', 'methylornithine', '5.4.99.58', '6.1.1.26')
pyl = {k: v for k, v in defs.items()
       if any(x in v.lower() for x in TERMS)}
print()
print(f"  pyrrolysine machinery, searched on {len(TERMS)} terms: {len(pyl)} KOs")
for k in sorted(pyl):
    print(f"    {k}  {pyl[k][:70]}")
with open(os.path.join(D, 'pyl_kos.tsv'), 'w') as fh:
    for k in sorted(pyl):
        fh.write(f"{k}\tPYL_MACHINERY\t{pyl[k]}\n")
if not pyl:
    print("  !! NONE of these terms matches anything in ko_list:")
    print(f"     {', '.join(TERMS)}")
    print("     That is a statement about ko_list, not about these genomes — KEGG")
    print("     may simply not carry orthologue families for the pyl operon. The")
    print("     orthogonal control cannot run, and section 4 will say so rather")
    print("     than report an absence. To look for other wording:")
    print("       grep -iE 'ornithine|lysine ligase|tRNA ligase' <ko_list> | cut -f1,12")
if bad:
    print(f"  !! {bad} KO(s) disagree with KEGG. NOTHING IS RUN.")
    sys.exit(1)
sys.exit(0)
PYG
[ $? -eq 0 ] || { runlog "**34b ABORTED: KO guard failed**"; exit 4; }

# ============================================================ 1. PROFILES
echo
echo "######## 1. PROFILES ########"
mapfile -t KOS < <( { cut -f1 "$D/targets.tsv"; cut -f1 "$D/pyl_kos.tsv" 2>/dev/null; } | sort -u )
got=0; missing=()
for ko in "${KOS[@]}"; do
  [ -s "$D/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$U/profiles/${ko}.hmm" "$WS2/36_gapclose/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$D/profiles/"; break; }
  done
  if [ ! -s "$D/profiles/${ko}.hmm" ]; then
    wget -q -O "$D/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" 2>/dev/null || true
    grep -q '^HMMER3' "$D/profiles/${ko}.hmm" 2>/dev/null || rm -f "$D/profiles/${ko}.hmm"
  fi
  [ -s "$D/profiles/${ko}.hmm" ] && got=$((got+1)) || missing+=("$ko")
done
echo "  $got/${#KOS[@]} profiles"
printf '%s\n' "${missing[@]:-}" | grep -v '^$' > "$D/not_searched.txt" || true
if [ "${#missing[@]}" -gt 0 ]; then
  echo "  !! NO PROFILE for: ${missing[*]}"
  echo "     These are carried through the report as NOT SEARCHED. Nothing about"
  echo "     them is reported as an absence — that distinction is the one step 36"
  echo "     lost, and it turned an unrun search into 190 downgraded network edges."
  echo "     To get them: see runbook step 3 (the KOfam profile tarball)."
fi
cat "$D"/profiles/*.hmm > "$D/all.hmm" 2>/dev/null

# EVERY CACHE BELOW IS KEYED ON THE PROFILE SET, NOT ON FILE EXISTENCE.
# The run of 2026-09-04 added the three pyrrolysine profiles and then reused a
# nomination and a segment search built from the ten-profile set — so the new
# profiles were never searched against anything, and section 4 reported "no pyl
# machinery found" when the truth was "not searched". That is the step-36
# defect for the third time (Rule 34), and existence-based caching is what
# keeps producing it. A stamp makes the cache honest.
PSTAMP=$(md5sum "$D/all.hmm" 2>/dev/null | cut -c1-12)
echo "  profile-set stamp: $PSTAMP"
# An ABSENT stamp is a mismatch, not a match. The first stamped run had no
# .pstamp, so this block did not fire and the segment search stayed cached from
# the ten-profile set — the same defect one layer down.
if [ "$(cat "$D/.pstamp" 2>/dev/null || echo NONE)" != "$PSTAMP" ]; then
  echo "  !! the profile set CHANGED since the last run."
  echo "     Nomination, read-through segments and the segment search are all"
  echo "     downstream of it, so they are being rebuilt. Nothing is reused."
  rm -f "$D"/nom/*.tbl "$D"/nom/*.stamp \
        "$D/rt/segments.faa" "$D/rt/amber_pos.tsv" \
        "$D/hits/rt.tbl" "$D/hits/rt.domtbl"
fi
printf '%s\n' "$PSTAMP" > "$D/.pstamp"
[ -s "$D/all.hmm" ] || { echo "  !! no profiles"; runlog "**34b ABORTED: no profiles**"; exit 5; }

# ============================================================ 2. NOMINATE
echo
echo "######## 2. NOMINATING CONTIGS (permissive, E<1e-5) ########"
# Nomination is deliberately permissive: a Pyl gene truncated by the gene caller
# leaves a FRAGMENT, and a fragment can score below a curated full-length
# threshold. Nominating at E<1e-5 and adjudicating at the curated threshold on
# the read-through translation is intentional, and the two cutoffs are never
# mixed in the reporting.
declare -A FAA
for d in "$A"/*/; do n=$(basename "$d"); [ -s "$d/$n.faa" ] && FAA["$n"]="$d/$n.faa"; done
for code in SED LIQ; do
  [ -s "$U/prot/${code}_all.faa" ] && FAA["ASM_${code}"]="$U/prot/${code}_all.faa"
done
echo "  proteomes to scan: ${#FAA[@]} (${#FAA[@]} = 29 MAGs + 2 assemblies when both exist)"
if [ "${#FAA[@]}" -eq 0 ]; then
  echo "  !! no proteomes. MAG proteomes come from step 10; assembly proteins from"
  echo "     step 34 stage A ($U/prot/{SED,LIQ}_all.faa)."
  runlog "**34b ABORTED: no proteomes**"; exit 6
fi
i=0
for name in "${!FAA[@]}"; do
  i=$((i+1)); tb="$D/nom/${name}.tbl"; st="$D/nom/${name}.stamp"
  if [ -s "$tb" ] && [ "$(cat "$st" 2>/dev/null)" = "$PSTAMP" ]; then
    printf '\r  [%2d/%2d] %-32s cached' "$i" "${#FAA[@]}" "$name"; continue
  fi
  printf '\r  [%2d/%2d] %-32s      ' "$i" "${#FAA[@]}" "$name"
  run_tool hmmsearch --cpu "$THREADS" -E 1e-5 --domE 1e-5 --tblout "$tb.part" \
      "$D/all.hmm" "${FAA[$name]}" >/dev/null 2>"$D/nom/${name}.err"
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$tb.part" ]; then
    printf '\r  [%2d/%2d] %-32s FAILED rc=%s\n' "$i" "${#FAA[@]}" "$name" "$rc"
    head -3 "$D/nom/${name}.err" | sed 's/^/       /'; rm -f "$tb.part"; continue
  fi
  mv "$tb.part" "$tb"; echo "$PSTAMP" > "$st"
done
echo
echo "  nomination tables: $(ls -1 "$D/nom"/*.tbl 2>/dev/null | wc -l)/${#FAA[@]}"

# ============================================================ 3. READ-THROUGH
echo
echo "######## 3. RE-TRANSLATING CANDIDATE CONTIGS THROUGH THE AMBER CODON ########"
python3 - "$D" "$F" "$WS2" <<'PYR'
import sys, os, re, gzip
D, F, WS2 = sys.argv[1], sys.argv[2], sys.argv[3]

# which contigs were nominated, per source
nom = {}
for fn in sorted(os.listdir(os.path.join(D, 'nom'))):
    if not fn.endswith('.tbl'):
        continue
    src = fn[:-4]
    s = set()
    for line in open(os.path.join(D, 'nom', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if not f:
            continue
        p = f[0]
        # Bakta/prodigal protein ids are <contig>_<n>; strip the final _n
        s.add(p.rsplit('_', 1)[0] if re.search(r'_\d+$', p) else p)
    if s:
        nom[src] = s
print(f"  sources with a nomination: {len(nom)}")
tot = sum(len(v) for v in nom.values())
print(f"  candidate contigs: {tot}")
if tot == 0:
    print("  nothing nominated — no methylamine methyltransferase candidate anywhere.")
    print("  That is a RESULT, but only if section 4's control passes.")

CODON = {}
B = 'TCAG'
AAS = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
i = 0
for a in B:
    for b in B:
        for c in B:
            CODON[a + b + c] = AAS[i]; i += 1
# READ THROUGH THE AMBER CODON. TAA and TGA stay real stops.
#
# The amber positions are tracked SEPARATELY rather than written into the
# sequence as a marker character. A first version translated TAG to 'X' and
# counted X's later — but any codon containing an N also translates to X, so
# an assembly gap would have been counted as a pyrrolysine site and inflated
# the very signal this step rests on. Same class of error as step 32d comparing
# tuples lexicographically: the shortcut looked equivalent and was not.
CODON['TAG'] = 'X'
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def translate(seq):
    """-> (protein, set of indices that came from a real in-frame TAG)"""
    aa = []
    amb = set()
    for j in range(0, len(seq) - 2, 3):
        c = seq[j:j + 3]
        if c == 'TAG':
            amb.add(j // 3)
            aa.append('X')
        else:
            aa.append(CODON.get(c, 'X'))
    return ''.join(aa), amb

def frames(seq):
    rc = seq.translate(COMP)[::-1]
    for f in range(3):
        p, a = translate(seq[f:]); yield f + 1, p, a
    for f in range(3):
        p, a = translate(rc[f:]); yield -(f + 1), p, a

def readfa(path):
    op = gzip.open if path.endswith('.gz') else open
    name, buf = None, []
    with op(path, 'rt', errors='replace') as fh:
        for line in fh:
            if line[:1] == '>':
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, ''.join(buf)

# where each source's contigs live
srcfa = {}
bd = os.path.join(F, 'bins')
if os.path.isdir(bd):
    for fn in os.listdir(bd):
        if fn.endswith(('.fa', '.fna', '.fasta')):
            srcfa[os.path.splitext(fn)[0]] = os.path.join(bd, fn)
for code in ('SED', 'LIQ'):
    for cand in (os.path.join(WS2, '04_polish', f'{code}_asm_final.fasta'),):
        if os.path.exists(cand):
            srcfa[f'ASM_{code}'] = cand

MINAA = 60
nseg = 0
namber = 0
out = open(os.path.join(D, 'rt', 'segments.faa'), 'w')
posf = open(os.path.join(D, 'rt', 'amber_pos.tsv'), 'w')
posf.write('segment\tpositions_1based_in_segment\n')
for src, contigs in sorted(nom.items()):
    fa = srcfa.get(src)
    if not fa:
        print(f"  !! no contig fasta for {src} — its candidates cannot be re-translated")
        continue
    want = contigs
    for cname, seq in readfa(fa):
        if cname not in want:
            continue
        seq = seq.upper()
        for fr, prot, amb in frames(seq):
            for m in re.finditer(r'[^*]+', prot):
                s = m.group(0)
                if len(s) < MINAA:
                    continue
                # ambers inside THIS segment only, counted from tracked positions
                # positions are recorded RELATIVE TO THE SEGMENT, 1-based, so
                # they can be tested against hmmsearch's alignment envelope.
                # Presence anywhere in the segment is not the test: in a 100-aa
                # stretch of arbitrary frame a TAG turns up by chance about four
                # times in five. Inside a scored alignment of real coding
                # sequence it essentially never does, because TAG is a stop and
                # selection removes it.
                pos = sorted(p - m.start() + 1 for p in amb
                             if m.start() <= p < m.end())
                nseg += 1
                namber += len(pos)
                sid = f"{src}|{cname}|f{fr}|{m.start()}|{len(s)}|amb{len(pos)}"
                out.write(f">{sid}\n")
                for j in range(0, len(s), 60):
                    out.write(s[j:j + 60] + '\n')
                if pos:
                    posf.write(sid + '\t' + ','.join(str(x) for x in pos) + '\n')
out.close()
posf.close()
print(f"  read-through segments written: {nseg}")
print(f"  in-frame TAG codons inside those segments: {namber}")
print("  (TAG read through; TAA and TGA remain stops; N-containing codons")
print("   translate to X but are NOT counted as ambers)")
PYR

if [ -s "$D/rt/segments.faa" ]; then
  echo
  echo "  searching the read-through segments at CURATED thresholds ..."
  # keyed on BOTH inputs: the profile set and the segments actually searched.
  SSTAMP=$(md5sum "$D/rt/segments.faa" | cut -c1-12)
  RTSTAMP="${PSTAMP}_${SSTAMP}"
  if [ -s "$D/hits/rt.tbl" ] \
     && [ "$(cat "$D/hits/rt.stamp" 2>/dev/null || echo NONE)" = "$RTSTAMP" ]; then
    echo "  cached (profiles and segments both unchanged: $RTSTAMP)"
  else
    rm -f "$D/hits/rt.tbl" "$D/hits/rt.domtbl"
    run_tool hmmsearch --cpu "$THREADS" -E 1e-5 --tblout "$D/hits/rt.tbl.part" \
        --domtblout "$D/hits/rt.domtbl.part" "$D/all.hmm" "$D/rt/segments.faa" \
        >/dev/null 2>"$D/hits/rt.err"
    rc=$?
    if [ "$rc" -ne 0 ] || [ ! -s "$D/hits/rt.tbl.part" ]; then
      echo "  !! hmmsearch on segments exited $rc"; head -3 "$D/hits/rt.err" | sed 's/^/     /'
      rm -f "$D/hits/rt.tbl.part" "$D/hits/rt.domtbl.part"
    else
      mv "$D/hits/rt.tbl.part" "$D/hits/rt.tbl"
      mv "$D/hits/rt.domtbl.part" "$D/hits/rt.domtbl"
      echo "$RTSTAMP" > "$D/hits/rt.stamp"
      echo "  searched $(grep -c -v '^#' "$D/hits/rt.tbl") hit rows"
    fi
  fi
else
  echo "  no segments to search."
fi

# ============================================================ 4. INTERPRET
echo
python3 - "$D" "$F" "$U" <<'PY'
import sys, os, re, csv
from collections import defaultdict, Counter

D, F, U = sys.argv[1], sys.argv[2], sys.argv[3]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

# ---- KO metadata -------------------------------------------------------
role, gene = {}, {}
for line in open(os.path.join(D, 'targets.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 4:
        gene[f[0]] = f[1]; role[f[0]] = f[2]
pylko = {}
pp = os.path.join(D, 'pyl_kos.tsv')
if os.path.exists(pp):
    for line in open(pp):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 3:
            pylko[f[0]] = f[2]
            gene.setdefault(f[0], f[0]); role.setdefault(f[0], 'PYL_MACHINERY')

# A KO with no profile was NOT SEARCHED. Same distinction step 36 lost.
NOTSEARCHED = set()
_ns = os.path.join(D, 'not_searched.txt')
if os.path.exists(_ns):
    NOTSEARCHED = {l.strip() for l in open(_ns) if l.strip()}
else:
    NOTSEARCHED = {k for k in list(role) + list(pylko)
                   if not os.path.exists(os.path.join(D, 'profiles', f'{k}.hmm'))}

thr = {}
with open(os.path.join(D, 'ko_list'), errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    ti = hdr.index('threshold') if 'threshold' in hdr else 1
    si = hdr.index('score_type') if 'score_type' in hdr else 2
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > max(ti, si):
            try:
                thr[f[0]] = (float(f[ti]), f[si].strip().lower())
            except ValueError:
                thr[f[0]] = (None, f[si].strip().lower())

# ---- read-through hits, adjudicated ------------------------------------
# tblout: 0 target 2 query 4 E(full) 5 score(full) 7 E(dom) 8 score(dom)
# amber positions per segment, 1-based within the segment
ampos = defaultdict(list)
app = os.path.join(D, 'rt', 'amber_pos.tsv')
if os.path.exists(app):
    for line in open(app, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if len(f) == 2 and f[0] != 'segment' and f[1]:
            ampos[f[0]] = [int(x) for x in f[1].split(',') if x.isdigit()]

# alignment envelope per (segment, KO), from --domtblout
# columns (0-based): 0 target 3 query 13 domain score 17 ali_from 18 ali_to
env = {}
dp = os.path.join(D, 'hits', 'rt.domtbl')
if os.path.exists(dp):
    for line in open(dp, errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 19:
            continue
        try:
            sc = float(f[13]); a, b = int(f[17]), int(f[18])
        except ValueError:
            continue
        k = (f[0], f[3])
        if k not in env or sc > env[k][0]:
            env[k] = (sc, a, b)

EDGE = 0.05      # an amber within 5% of either end of the envelope is TERMINAL
passed = []      # (ko, src, contig, seg, score, n_amber, ali_len, n_internal, rels)
nothr = set(); below = Counter()
tb = os.path.join(D, 'hits', 'rt.tbl')
if os.path.exists(tb):
    best = {}
    for line in open(tb, errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 9:
            continue
        seg, ko = f[0], f[2]
        if ko not in gene:
            ko = f[3] if len(f) > 3 and f[3] in gene else None
            if ko is None:
                continue
        t, st = thr.get(ko, (None, ''))
        try:
            sc = float(f[8]) if st and st.startswith('dom') else float(f[5])
        except ValueError:
            continue
        if t is None:
            nothr.add(ko)
        elif sc < t:
            below[ko] += 1
            continue
        # keep the best-scoring KO per segment: a segment is one protein
        if seg not in best or sc > best[seg][1]:
            best[seg] = (ko, sc)
    for seg, (ko, sc) in best.items():
        parts = seg.split('|')
        src = parts[0] if parts else '?'
        contig = parts[1] if len(parts) > 1 else '?'
        # WHERE IN THE ALIGNMENT THE AMBER SITS DECIDES WHAT IT IS.
        #
        # Reading through TAG has a systematic cost I did not account for: TAG is
        # also an ordinary stop codon, used by roughly a fifth of bacterial genes.
        # Every such gene yields a segment whose TRUE C-terminal stop appears as an
        # internal position, and if hmmsearch's envelope overshoots the real end by
        # even a few residues, that stop is counted as a pyrrolysine site.
        #
        # A pyrrolysine residue sits MID-PROTEIN (MtmB's is near residue 202 of
        # ~460). A read-through stop sits at the very END. So the position within
        # the envelope separates them, and the verdict uses INTERNAL ambers only.
        e = env.get((seg, ko))
        if e:
            _s, a, b = e
            inside = [p for p in ampos.get(seg, ()) if a <= p <= b]
            alen = max(1, b - a + 1)
            rels = [(p - a) / max(1, b - a) for p in inside]
            namb = len(inside)
            nint = sum(1 for r in rels if EDGE <= r <= 1.0 - EDGE)
        else:
            namb = len(ampos.get(seg, ()))
            alen = 0; rels = []; nint = 0
        passed.append((ko, src, contig, seg, sc, namb, alen, nint, rels))

rule('1. WHAT THE CORRECTED KOs FIND')
say("  Step 34 aborted trying to resolve mtmB by keyword and I concluded KEGG had")
say("  no KO for it. KEGG has one — K16176, which it calls 'methylamine---corrinoid")
say("  protein Co-methyltransferase', EC 2.1.1.248. The word 'monomethylamine'")
say("  appears in KEGG only on K16177, the corrinoid protein. The search string was")
say("  wrong, not the database.")
say()
bysrc = defaultdict(lambda: defaultdict(int))
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    bysrc[ko][src] += 1
say(f"  {'KO':<9}{'gene':<7}{'role':<18}{'segments':>10}{'sources':>9}   KEGG threshold")
say('-' * W)
# ORDER COMES FROM targets.tsv, NOT FROM A LIST RETYPED HERE.
# The first version hardcoded seven KOs and then three more were added to the
# targets — MtaB, MtaA and MtaC — so section 1 silently omitted them, including
# the methanol marker this step exists to test. A display list that has to be
# kept in sync with a data file will fall out of sync.
order = []
for line in open(os.path.join(D, 'targets.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 2 and f[0] not in order:
        order.append(f[0])
order += [k for k in sorted(pylko) if k not in order]
nsflag = False
for ko in order:
    tv = thr.get(ko, (None, ''))[0]
    if ko in NOTSEARCHED:
        nsflag = True
        say(f"  {ko:<9}{gene.get(ko, '?'):<7}{role.get(ko, ''):<18}"
            f"{'NOT SEARCHED':>19}   {('%.1f' % tv) if tv else 'none'}")
        continue
    n = sum(bysrc[ko].values())
    say(f"  {ko:<9}{gene.get(ko, '?'):<7}{role.get(ko, ''):<18}{n:>10}"
        f"{len(bysrc[ko]):>9}   {('%.1f' % tv) if tv else 'none'}")
if nsflag:
    say()
    say("  NOT SEARCHED means no HMM profile could be fetched for that family. It is")
    say("  not a zero. Nothing below reports those KOs as absent.")
if below:
    say()
    say("  segments that hit but failed their family's curated threshold: "
        + ', '.join(f'{k}:{v}' for k, v in below.most_common()))
if nothr:
    say(f"  no curated threshold for {', '.join(sorted(nothr))} — kept on E-value alone")

# ---- the amber test ----------------------------------------------------
rule('2. THE AMBER TEST — AND ITS CONTROL INSIDE THE SAME PATHWAY')
say("  MtmB, MtbB and MttB carry an in-frame UAG encoding pyrrolysine. MtbA does")
say("  not, and neither do the three corrinoid partners. All seven were searched")
say("  the same way on the same re-translated segments, so if the amber signal is")
say("  real it must separate the first group from the other two.")
say()
say("  The amber is counted ONLY INSIDE hmmsearch's alignment envelope, never")
say("  merely somewhere in the segment. Across 100 residues of arbitrary reading")
say("  frame a TAG turns up by chance about four times in five, so 'contains one'")
say("  would be a coin-flip dressed as evidence. Inside a scored alignment of real")
say("  coding sequence it is close to impossible — TAG is a stop, and selection")
say("  has removed it everywhere except where pyrrolysine is actually used.")
say()
# [segments, with_internal, total_internal, ali_aa, with_any, total_terminal]
grp = defaultdict(lambda: [0, 0, 0, 0, 0, 0])
noenv = 0
allrel = defaultdict(list)
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    r = role.get(ko, '?')
    grp[r][0] += 1
    if nint:
        grp[r][1] += 1
    grp[r][2] += nint
    grp[r][3] += alen
    if namb:
        grp[r][4] += 1
    grp[r][5] += (namb - nint)
    allrel[r].extend(rels)
    if alen == 0:
        noenv += 1
say("  An amber within 5% of either end of the alignment is counted as TERMINAL —")
say("  most likely the gene's own TAG stop codon read through — and excluded. The")
say("  verdict uses INTERNAL ambers only.")
say()
say(f"  {'class':<22}{'segs':>7}{'any UAG':>9}{'TERMINAL':>10}{'INTERNAL':>10}"
    f"{'% internal':>12}{'int/1000aa':>12}")
say('-' * W)
for r in ('PYL', 'NONPYL_CTRL', 'CORRINOID_DECOY', 'PYL_MACHINERY'):
    n, wi, ti, aa, wa, term = grp.get(r, [0, 0, 0, 0, 0, 0])
    if not n:
        say(f"  {r:<22}{0:>7}{'-':>9}{'-':>10}{'-':>10}{'-':>12}{'-':>12}")
        continue
    rate = (1000.0 * ti / aa) if aa else 0.0
    say(f"  {r:<22}{n:>7}{wa:>9}{term:>10}{ti:>10}{100.0 * wi / n:>11.0f}%{rate:>12.2f}")
say()
if noenv:
    say(f"  ({noenv} segment(s) had no domain row; their ambers were counted over the")
    say("   whole segment and are the weaker evidence.)")
say()

# ---- A SEGMENT IS NOT A GENE, AND THE FIRST GATE FORGOT THAT -----------
# v1 required at least 60% of PYL segments to carry an in-alignment amber
# before it would pass. On the real data that rejected a 24-fold separation
# (31% vs 1%, Fisher p = 7e-9) as a failure.
#
# The 60% was ill-founded for two reasons:
#   1. hmmsearch's envelope need not span the pyrrolysine position. An
#      alignment covering the N-terminal half of MtmB misses a residue that
#      sits mid-protein, and the segment is then a true hit with no amber.
#   2. 653 candidate contigs searched in six frames produce many partial and
#      spurious segments per real gene. The denominator is segments, not genes.
#
# What the test is actually asking is whether the amber SEPARATES the classes.
# That is a comparison, so it gets a comparison's statistics — and a per-gene
# view alongside the per-segment one.
gene_amber = defaultdict(bool)      # (source, gene) -> any in-alignment amber
gene_role = {}
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    k = (src, gene.get(ko, ko))
    gene_role[k] = role.get(ko, '?')
    if nint:                       # INTERNAL ambers only, same rule as the verdict
        gene_amber[k] = True
    else:
        gene_amber.setdefault(k, False)
gclass = defaultdict(lambda: [0, 0])
for k, has in gene_amber.items():
    r = gene_role.get(k, '?')
    gclass[r][0] += 1
    if has:
        gclass[r][1] += 1
say(f"  {'class':<22}{'genes hit':>11}{'with UAG':>10}{'% of genes':>12}")
say('-' * W)
for r in ('PYL', 'NONPYL_CTRL', 'CORRINOID_DECOY', 'PYL_MACHINERY'):
    n, wa = gclass.get(r, [0, 0])
    say(f"  {r:<22}{n:>11}{wa if n else '-':>10}"
        + (f"{100.0 * wa / n:>11.0f}%" if n else f"{'-':>12}"))
say()
say("  A gene counts as amber-bearing if ANY of its segments carries an INTERNAL")
say("  amber. That is the right unit: the question is whether the gene has a")
say("  pyrrolysine site, not whether every partial alignment happens to cover it.")
say()

def fisher_greater(a, b, c, d):
    """one-sided Fisher exact, P(as extreme or more) for a/(a+b) > c/(c+d)"""
    from math import comb
    r1, r2, c1 = a + b, c + d, a + c
    n = r1 + r2
    if min(r1, r2, c1) < 0 or n == 0:
        return 1.0
    tot = comb(n, c1)
    if tot == 0:
        return 1.0
    s = 0
    for x in range(a, min(r1, c1) + 1):
        s += comb(r1, x) * comb(r2, c1 - x)
    return s / tot

pn, pw, pt, pa = grp.get('PYL', [0, 0, 0, 0, 0, 0])[:4]
cn, cw, ct, ca = grp.get('NONPYL_CTRL', [0, 0, 0, 0, 0, 0])[:4]
dn, dw, dt, da = grp.get('CORRINOID_DECOY', [0, 0, 0, 0, 0, 0])[:4]
neg_n = cn + dn; neg_w = cw + dw; neg_t = ct + dt; neg_a = ca + da
AMBER = None
if pn == 0:
    say("  [NO TEST] no pyrrolysine-protein segments passed threshold, so there is")
    say("            nothing to test. Section 1 already says so.")
elif neg_n < 10 or pn < 10:
    say(f"  [UNDERPOWERED] {pn} pyrrolysine and {neg_n} control segments is not enough")
    say("                 to separate two rates. The amber evidence is not usable")
    say("                 either way — that is a power statement, not a negative.")
else:
    pr = pw / pn; nr = neg_w / neg_n
    prate = (1000.0 * pt / pa) if pa else 0.0
    nrate = (1000.0 * neg_t / neg_a) if neg_a else 0.0
    P = fisher_greater(pw, pn - pw, neg_w, neg_n - neg_w)
    say(f"  pyrrolysine proteins with an INTERNAL UAG     : {pw}/{pn}  ({100 * pr:.0f}%)"
        f"   {prate:.2f} per 1000 aligned aa")
    say(f"  non-pyrrolysine controls with one             : {neg_w}/{neg_n}  ({100 * nr:.0f}%)"
        f"   {nrate:.2f} per 1000 aligned aa")
    say(f"  Fisher exact, one-sided                      : p = {P:.3g}")
    say()
    say("  Three things must agree: the proportion of segments carrying one, the")
    say("  rate per aligned residue, and that the difference is not chance. The")
    say("  proportion alone rewards long alignments; the rate alone can be carried")
    say("  by one protein with several; a ratio alone can come from small numbers.")
    sep = pr >= 3 * max(nr, 1e-9)
    ratesep = prate >= 3 * max(nrate, 1e-9)
    if sep and ratesep and P < 1e-3:
        AMBER = True
        say(f"  [PASS] the amber separates the classes {pr / max(nr, 1e-9):.0f}-fold by proportion and"
            f" {prate / max(nrate, 1e-9):.0f}-fold by")
        say("         rate, at p < 0.001. A segment carrying both the profile hit and")
        say("         an in-alignment UAG is a real pyrrolysine methyltransferase, and")
        say("         that is a sequence-level identification rather than a score.")
        say()
        say("         NOTE the absolute level: most PYL segments do NOT carry one, and")
        say("         that is expected — an alignment covering only part of the protein")
        say("         misses a residue that sits mid-sequence. The evidence is in the")
        say("         CONTRAST with the controls, not in the absolute fraction.")
    else:
        AMBER = False
        say("  [FAIL] the in-frame UAG does not separate proteins that can contain")
        say("         pyrrolysine from proteins that cannot"
            + ("" if P < 1e-3 else f" (p = {P:.3g})") + ".")
        say("         It is being read out of spurious frames. The amber evidence is")
        say("         void; the profile hits still stand on their thresholds, but")
        say("         nothing here confirms them as pyrrolysine proteins.")
    if neg_w:
        say()
        say(f"  {neg_w} control segment(s) DID carry an in-alignment UAG. Those are listed")
        say("  in markers.tsv with in_frame_UAG > 0 and a NONPYL or DECOY role; a real")
        say("  read-through, an assembly error and a spurious frame all look like this,")
        say("  and at this count they do not change the verdict either way.")
say()
say("  This is also the MttB superfamily test. Step 34's mttB decoy took 312,347")
say("  sediment reads because non-pyrrolysine MttB homologues — MtgB and kin,")
say("  acting on glycine betaine — are common. Those carry no amber codon.")

# ---- confirmed callers -------------------------------------------------
rule('3. WHICH GENOMES CARRY A CONFIRMED METHYLAMINE METHYLTRANSFERASE')
conf = defaultdict(set)     # src -> {gene}
allhit = defaultdict(set)
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    if role.get(ko) == 'PYL':
        allhit[src].add(gene[ko])
        if nint:
            conf[src].add(gene[ko])
say(f"  {'source':<30}{'profile hit':<22}{'with in-frame UAG':<24}")
say('-' * W)
for src in sorted(set(list(allhit) + list(conf))):
    say(f"  {src:<30}{','.join(sorted(allhit.get(src, ()))) or '-':<22}"
        f"{','.join(sorted(conf.get(src, ()))) or '-':<24}")
if not allhit:
    say("  none.")

# ---- the methanol arm --------------------------------------------------
rule('3b. THE METHANOL ARM — its first genome-level test')
say("  Step 34 called this THE UNTESTED ARM of pectin -> methanol -> CH4 and could")
say("  only reach it through reads, where MtaB came back at 2.2x mcrA and broke the")
say("  ceiling every methanogen imposes. Here it is searched in the genomes and the")
say("  assemblies, where a hit is a gene and not a read share.")
say()
say("  Confirmed against ko_list 2026-08-27: K04480 is 'methanol---5-hydroxybenz-")
say("  imidazolylcobamide Co-methyltransferase [EC:2.1.1.90]', which IS MtaB, so")
say("  step 34's label was right. K14080 is MtaA, the second step. K14081 is MtaC.")
say()
METHANOL = {'K04480': 'MtaB (methanol -> corrinoid)',
            'K14081': 'MtaC (corrinoid partner)',
            'K14080': 'MtaA (corrinoid -> coenzyme M)'}
msrc = defaultdict(set)
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    if ko in METHANOL:
        msrc[ko].add(src)
say(f"  {'KO':<9}{'protein':<34}{'sources':>9}   where")
say('-' * W)
for ko, lab in METHANOL.items():
    s = sorted(msrc.get(ko, ()))
    say(f"  {ko:<9}{lab:<34}{len(s):>9}   {', '.join(s)[:38] or '-'}")
say()
mtab = msrc.get('K04480', set())
if not mtab:
    say("  MtaB is absent from every MAG and both assemblies. Step 34's read-level")
    say("  signal, if real, therefore comes from organisms below the assembly")
    say("  threshold — the same conclusion dsrAB reached, reached the same way.")
elif msrc.get('K14081') & mtab:
    say("  MtaB and its corrinoid partner MtaC occur in the same source(s). A")
    say("  methyltransferase without its corrinoid is a lone hit; the pair is a")
    say("  pathway, and that is the difference worth reporting.")
else:
    say("  MtaB is present but MtaC is not found alongside it. Treat that as a lone")
    say("  hit rather than a functioning methanol route until the partner is found —")
    say("  the same operon logic that made dsrA and dsrB check each other.")
say()
say("  NOTE ON THE DECOY ROLE: MtaC is also a promiscuous B12-binding protein, the")
say("  same family that let MtmC absorb 1,978 reads in step 34. It counts here as")
say("  corroboration only when MtaB is already present, never on its own.")

# ---- orthogonal control ------------------------------------------------
rule('4. ORTHOGONAL CORROBORATION — PYRROLYSINE MACHINERY')
say("  pylB/pylC/pylD make pyrrolysine; pylS/pylT charge it. They share no")
say("  homology with the methyltransferases, so agreement between the two gene")
say("  sets is independent evidence — worth more than either score alone.")
say()
pylhit = defaultdict(set)
for ko, src, contig, seg, sc, namb, alen, nint, rels in passed:
    if role.get(ko) == 'PYL_MACHINERY':
        pylhit[src].add(ko)
PYL_SEARCHED = [k for k in pylko if k not in NOTSEARCHED]
if pylko and not PYL_SEARCHED:
    say(f"  ko_list DOES carry the machinery — {len(pylko)} KOs resolved:")
    for k in sorted(pylko):
        say(f"    {k}  {pylko[k][:66]}")
    say()
    say("  But NO PROFILE could be fetched for any of them, so NONE WAS SEARCHED.")
    say("  This control did not run. It did not fail, and it found nothing because")
    say("  it never looked — those are different things, and collapsing them is the")
    say("  defect that gave step 36 a table of 190 fictional edge downgrades.")
    say()
    say("  Nothing about the pyrrolysine machinery is claimed here in either")
    say("  direction. The amber test in section 2 stands on its within-pathway")
    say("  controls without it.")
    say()
    say("  To close this: fetch the profiles for "
        + ", ".join(sorted(pylko)) + " (runbook step 3)")
    say("  and re-run — everything else is cached, so it costs about a minute.")
elif not pylko:
    say("  This control COULD NOT RUN: no KO in ko_list matched any of the terms")
    say("  searched (pyrrolysin, methylornithine, EC 5.4.99.58, EC 6.1.1.26). KEGG")
    say("  names pylB/pylC/pylD by their chemistry — methylornithine synthase and")
    say("  kin — and may not carry orthologue families for them at all.")
    say()
    say("  NOT SEARCHED is not ABSENT. Nothing about the pyrrolysine machinery is")
    say("  claimed here in either direction, and the amber test in section 2 stands")
    say("  on its own controls without it. To go looking for other wording:")
    say("    grep -iE 'ornithine|lysine ligase|tRNA ligase' 34b_methylamine/ko_list \\")
    say("      | cut -f1,12")
else:
    if len(PYL_SEARCHED) < len(pylko):
        say(f"  NOT SEARCHED (no profile): "
            + ", ".join(sorted(set(pylko) - set(PYL_SEARCHED))))
        say("  Only the searched families below can contribute either way.")
        say()
    say(f"  {'source':<30}{'MT genes':>10}{'pyl machinery KOs':>20}   agreement")
    say('-' * W)
    agree = disagree = 0
    for src in sorted(set(list(allhit) + list(pylhit))):
        nm = len(allhit.get(src, ()))
        npl = len(pylhit.get(src, ()))
        if nm and npl:
            v = 'both'; agree += 1
        elif nm:
            v = 'MT only — no pyl machinery found'; disagree += 1
        else:
            v = 'pyl machinery only'
        say(f"  {src:<30}{nm:>10}{npl:>20}   {v}")
    say()
    say(f"  {agree} source(s) carry both; {disagree} carry methyltransferases with no")
    say("  pyrrolysine machinery detected. The second class is not disproof — a")
    say("  partial MAG loses genes — but it is not corroborated either, and the")
    say("  completeness of those genomes has to travel with the claim.")

# ---- what step 34 said, relabelled -------------------------------------
rule('5. STEP 34 RESTATED WITH THE CORRECT LABELS')
say("  Step 34's read-level numbers stand. What changes is what they are called.")
say()
say(f"  {'step 34 called it':<22}{'it is':<8}{'KO':<9}   status")
say('-' * W)
say(f"  {'mtbB':<22}{'MtbB':<8}{'K16178':<9}   correct; 146 sediment reads, not reportable")
say(f"  {'mttB (decoy)':<22}{'MttB':<8}{'K14083':<9}   superfamily magnet, 312,347 reads")
say(f"  {'mtaB':<22}{'MtaB':<8}{'K04480':<9}   CORRECT as labelled; EC 2.1.1.90 confirms it")
say(f"  {'mtmB (withdrawn)':<22}{'MtmB':<8}{'K16176':<9}   EXISTS. Never searched at read level.")
say()
say("  >> GAP, stated as a gap: K16176 has not been run against the reads by any")
say("     step. Step 34's read pipeline can do it — the marker row is simply")
say("     K16176 with MtmC (K16177) as its decoy rather than as its target.")
say()
say("  The methanol arm is no longer unnamed: MtaB K04480, MtaC K14081, MtaA")
say("  K14080, all confirmed against ko_list. Section 3b tests it in the genomes.")
say("  What remains open at READ level is MtmB (K16176) only.")

rule('WHAT TO REPORT')
if not passed:
    say("  >> NOTHING FOUND. No methylamine methyltransferase passed a curated")
    say("     threshold in any MAG or either assembly. Given step 34 found 146")
    say("     mtbB reads, the likely reading is that these organisms are below the")
    say("     assembly threshold — the same conclusion dsrAB reached.")
elif AMBER is False:
    say("  >> The profile hits stand on their thresholds. The pyrrolysine evidence")
    say("     is void (section 2), so call them 'MttB-superfamily methyltransferases'")
    say("     and not trimethylamine methyltransferases.")
else:
    say("  1. Report methylamine methyltransferases only where the profile hit and")
    say("     the in-frame UAG coincide. That pair is a sequence-level")
    say("     identification and it is what separates MtmB/MtbB/MttB from the")
    say("     MttB superfamily that swallowed 312,347 reads in step 34.")
    say("  2. Say in Methods that standard gene calling truncates these three genes")
    say("     at the amber codon, and that this step re-translated candidate contigs")
    say("     reading through UAG. Any pipeline that did not do this under-detects")
    say("     them by construction — including step 34's own level 2.")
    if pylko and [k for k in pylko if k not in NOTSEARCHED]:
        say("  3. Carry the pyrrolysine-machinery agreement from section 4 with the")
        say("     claim; it is independent of every score in section 1.")
    else:
        say("  3. The orthogonal pyrrolysine-machinery control DID NOT RUN (section 4),")
        say("     so the amber evidence rests on its within-pathway controls alone —")
        say("     MtbA, MtaB, MtaA and the four corrinoid proteins. That is still a")
        say("     controlled comparison, but say in the text that the independent")
        say("     corroboration was not available rather than implying it was.")
say()
say("  LIMITS:")
say("   1. Nomination ran at E<1e-5 and adjudication at KEGG's curated thresholds.")
say("      Those are different cutoffs used for different jobs, deliberately, and")
say("      the nomination number must never be quoted as a count of genes.")
say("   2. A read-through segment is an open reading frame, not a called gene. It")
say("      has no start codon test and no RBS. It is evidence a gene is there.")
say("   3. Contigs not nominated were never re-translated. A methyltransferase")
say("      whose truncated fragment scored below E<1e-5 is invisible to this step.")
say("   4. Read level is untouched here; step 34's numbers stand with the labels")
say("      corrected in section 5.")

with open(os.path.join(D, 'markers.tsv'), 'w') as fh:
    fh.write('ko\tgene\trole\tsource\tcontig\tsegment\tscore\tin_frame_UAG\t'
             'internal_UAG\trel_positions\n')
    for ko, src, contig, seg, sc, namb, alen, nint, rels in sorted(passed):
        fh.write(f"{ko}\t{gene.get(ko,'')}\t{role.get(ko,'')}\t{src}\t{contig}\t"
                 f"{seg}\t{sc:.1f}\t{namb}\t{nint}\t"
                 + ','.join(f"{r:.3f}" for r in rels) + "\n")
open(os.path.join(D, 'VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'markers.tsv')}")
print(f"  written: {os.path.join(D, 'VERDICT.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "nomination ran" "$([ "$(ls -1 "$D/nom"/*.tbl 2>/dev/null | wc -l)" -ge 1 ] && echo 0 || echo 1)" || OK=1
check "verdict written" "$([ -s "$D/VERDICT.txt" ] && echo 0 || echo 1)" || OK=1
[ "$OK" -eq 0 ] && runlog "34b_methylamine complete" || runlog "**34b_methylamine FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
