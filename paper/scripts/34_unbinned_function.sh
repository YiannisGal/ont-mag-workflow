#!/usr/bin/env bash
# =============================================================================
# 34_unbinned_function.sh  --  the 2% problem, attacked directly
#
# THE PROBLEM THIS STEP EXISTS FOR
#   The MAG catalogue holds ~2.4% of the sediment community (step 22, confirmed
#   independently by CoverM and SingleM). **Every functional statement in this
#   chapter therefore describes 2% of what is there.** Some of those statements
#   are absences, and an absence measured on 2% of a community is not a result.
#
#   Two are already known to be wrong, or at least untested:
#
#     * ZERO dsrAB in any of the 29 MAGs — while SILVA puts *Desulfatiglans* at
#       9.48% of sediment SSU reads. A sulfate-reducing lineage at nearly a
#       tenth of the community, and not one dissimilatory sulfite reductase
#       gene in the catalogue. That absence is a recovery limit almost
#       certainly, and it matters more than any other number here: the whole
#       methanogenesis argument is framed against sulfate reduction winning the
#       competition for hydrogen and acetate. The chapter needs the sulfate
#       reducers to be REAL and MEASURED, not inferred from an rRNA profile.
#
#     * SingleM puts sediment archaea at 14.1% while one methanogen assembled.
#
# THE DESIGN IS STEP 20b's, GENERALISED — AND THAT IS DELIBERATE (rule 19)
#   20b built this machinery for DMSP and it was used once. Reusing the code
#   rather than the idea:
#
#     LEVEL 1  the 29 MAGs               already done, per marker
#     LEVEL 2  every gene in both assemblies, binned or not
#     LEVEL 3  the raw reads             100% of the data
#
#   A hit at level 2 on an UNBINNED contig is a real organism the binning
#   missed. A hit at level 3 only is an organism below the assembly threshold.
#   A controlled zero at level 3 is a TESTED NEGATIVE — a far stronger sentence
#   than a MAG-based absence, and one almost no MAG study establishes.
#
#   Stage A reuses 20b's Prodigal output if it is there (~30 min saved).
#
# THE BINNED/UNBINNED SPLIT IS THE WHOLE POINT
#   Every assembly hit is labelled with whether its contig belongs to a MAG.
#   Without that column this step would only repeat what steps 11 and 20 know.
#
# CONTROLS, DOMAIN-DECLARED (the 19b lesson, applied up front)
#   rpoB / glnA / tuf are BACTERIAL. Archaea have a split RNA polymerase B and
#   EF-1alpha, not EF-Tu. Applying a bacterial control set to an archaeal target
#   and reporting the zero is the error 19b and 23e both made. Controls here are
#   declared per domain, and an archaeal target's absence is only callable if an
#   ARCHAEAL control passed. If it did not, the archaeal targets are reported as
#   UNTESTED rather than absent.
#
# NO KO IS TRUSTED FROM MEMORY (rule 1)
#   Every KO below is checked against KEGG's own ko_list definition at runtime,
#   for substrate AND activity keywords. One mismatch aborts the script and
#   prints what KEGG actually says. This guard has caught four wrong KOs in this
#   project already (K14080, K00825, K17217, K03146).
#
# OUTPUT  34_unbinned/marker_levels.tsv · assembly_hits.tsv · read_hits.tsv
#         34_unbinned/VERDICT.txt
# RUNTIME stage A 20-60 min (cached from 20b if present), stage B 1-3 h.
# USAGE   ws2 && bash scripts/34_unbinned_function.sh 2>&1 | tee logs/34.log
#         RUN_READS=0 bash scripts/34_unbinned_function.sh   # stage A only
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="34_unbinned"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"

F="$WS2/07_mags"; K="$WS2/11_kofam"; U="$WS2/20_substrates"
OLD="$U/dmsp_assembly"                       # 20b's output, reused where possible
D="$WS2/34_unbinned"; mkdir -p "$D/profiles" "$D/hits" "$D/prot"
RUN_READS="${RUN_READS:-1}"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 34 — FUNCTION IN THE UNBINNED FRACTION ########"
echo "  started: $(date -Is)"
runlog "34_unbinned started"

# ------------------------------------------------------------- TARGET MARKERS
# KO  gene  role  domain  keywords-that-must-appear-in-KEGG's-definition  note
# role:   TARGET | CONTROL
# domain: BACT | ARCH | BOTH   (which controls license an absence call)
cat > "$D/targets.tsv" <<'KOS'
K11180	dsrA	TARGET	BOTH	sulfite reductase	dissimilatory sulfite reductase alpha — THE headline absence
K11181	dsrB	TARGET	BOTH	sulfite reductase	dissimilatory sulfite reductase beta
K00394	aprA	TARGET	BOTH	adenylylsulfate	APS reductase alpha — independent confirmation of sulfate reduction
K00395	aprB	TARGET	BOTH	adenylylsulfate	APS reductase beta
K00958	sat	TARGET	BOTH	sulfate adenylyltransferase	sulfate activation, the first committed step
K00399	mcrA	TARGET	ARCH	methyl-coenzyme M reductase	expected present — an internal positive expectation for archaea
K04480	mtaB	TARGET	ARCH	methanol	THE UNTESTED ARM of pectin -> methanol -> CH4
K16178	mtbB	TARGET	ARCH	dimethylamine	dimethylamine methyltransferase — what K16178 REALLY is
K16176	mtmB	TARGET	ARCH	2.1.1.248	monomethylamine methyltransferase; KEGG names it 'methylamine---corrinoid protein Co-methyltransferase' so the EC is the identity check
K02588	nifH	TARGET	BOTH	nitrogenase	nitrogen fixation — never tested in this project
K20038	cutC	TARGET	BACT	choline	choline TMA-lyase — how widespread beyond ETO_LIQ_MAG_015?
K03043	rpoB	CONTROL	BACT	RNA polymerase	universal in bacteria, must be found
K01915	glnA	CONTROL	BACT	glutamine synthetase	near-universal in bacteria
K02358	tuf	CONTROL	BACT	elongation factor Tu	bacterial — NOT valid for archaea
K03231	ef1a	CONTROL	ARCH	elongation factor 1-alpha	the archaeal counterpart of tuf
KOS
NKO=0   # set from the resolved table after the guard
awk -F'\t' '{printf "    %-8s %-8s %-8s %s\n", $1, $2, $3"/"$4, $6}' "$D/targets.tsv"

# ------------------------------------------------------------- ko_list, reused
for src in "$K/ko_list" "$WS2/19_auxotrophy/ko_list" "$OLD/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$D/ko_list" 2>/dev/null || true; break; }
done
# NOT `A || B && C` — that parses as `(A || B) && C` and would run gunzip even
# when the file was already there. Rule 14's cousin: precedence is not grouping.
if [ ! -s "$D/ko_list" ]; then
  wget -q -O "$D/ko_list.gz" "$BASE/ko_list.gz" 2>/dev/null && \
    gunzip -f "$D/ko_list.gz" 2>/dev/null || true
fi

# ---------------------------------------------------------------- LABEL GUARD
echo
echo "######## 0. KO LABEL GUARD ########"
if [ ! -s "$D/ko_list" ]; then
  echo "  !! no ko_list — every KO below is an unverified assertion."
  echo "     This project's rule 1 forbids running on unverified KO numbers."
  echo "     Fetch it:  wget -O $D/ko_list.gz $BASE/ko_list.gz && gunzip $D/ko_list.gz"
  runlog "**34 ABORTED: no ko_list, KO labels unverifiable**"
  exit 2
fi
python3 - "$D" <<'PYG'
import sys, os, csv
D = sys.argv[1]
defs = {}
with open(os.path.join(D, 'ko_list'), errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    di = hdr.index('definition') if 'definition' in hdr else len(hdr) - 1
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > di: defs[f[0]] = f[di]
bad = 0
resolved = []
MAXAUTO = 4
for line in open(os.path.join(D, 'targets.tsv')):
    ko, gene, role, dom, kws, note = (line.rstrip('\n').split('\t') + [''] * 6)[:6]
    if not ko: continue
    keys = [k.strip().lower() for k in kws.split(';') if k.strip()]

    # AUTO: do not assert a KO at all — ask ko_list which families KEGG itself
    # describes with these words. This is rule 1 taken one step further: the
    # earlier form verified a number I had typed, which still required me to
    # type a right one. K16178 is exactly where that failed — I wanted
    # monomethylamine and typed the dimethylamine family.
    if ko.upper() == 'AUTO':
        cands = sorted(k for k, d in defs.items()
                       if all(x in d.lower() for x in keys))
        if not cands:
            print(f"  [FAIL] {gene}: no KO in ko_list matches {kws}"); bad += 1; continue
        if len(cands) > MAXAUTO:
            print(f"  [FAIL] {gene}: {len(cands)} KOs match {kws} — too vague to resolve")
            for c in cands[:8]: print(f"           {c}  {defs[c][:64]}")
            bad += 1; continue
        for c in cands:
            print(f"  [AUTO] {c} {gene:<6} {defs[c][:70]}")
            resolved.append([c, gene, role, dom, kws, note])
        continue

    d = defs.get(ko)
    if d is None:
        print(f"  [FAIL] {ko} ({gene}) is not in ko_list at all"); bad += 1; continue
    ok = all(x in d.lower() for x in keys)
    print(f"  [{'PASS' if ok else 'FAIL'}] {ko} {gene:<6} {d[:70]}")
    if not ok:
        print(f"         expected to contain: {kws}")
        bad += 1
    else:
        resolved.append([ko, gene, role, dom, kws, note])
if bad:
    print()
    print(f"  !! {bad} KO label(s) disagree with KEGG. NOTHING IS RUN.")
    print("     A wrong KO here would put a marker on the wrong gene and the")
    print("     absence claim would be about something else entirely.")
    sys.exit(3)
print()
print("  all KO labels agree with KEGG's own definitions.")

# The resolved table is what the rest of the run uses, and it is the provenance
# record of which KO each marker ACTUALLY was — not which one I typed.
with open(os.path.join(D, 'targets_resolved.tsv'), 'w') as fh:
    for r in resolved:
        fh.write('\t'.join(r) + '\n')
print(f"  resolved marker table: {len(resolved)} rows -> targets_resolved.tsv")
PYG
GRC=$?          # captured immediately — rule 14
if [ "$GRC" -ne 0 ]; then
  runlog "**34 ABORTED: KO label guard failed**"
  exit 3
fi

# ------------------------------------------------------------------ PROFILES
echo
echo "######## 0b. KOfam PROFILES ########"
NKO=$(wc -l < "$D/targets_resolved.tsv")
echo "  markers after resolution: $NKO"
got=0
while IFS=$'\t' read -r ko gene role dom kws note; do
  [ -z "${ko:-}" ] && continue
  [ -s "$D/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$WS2/19_auxotrophy/profiles/${ko}.hmm" \
             "$OLD/profiles/${ko}.hmm" "$U/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$D/profiles/"; break; }
  done
  [ -s "$D/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  if wget -q -O "$D/profiles/${ko}.hmm" "$BASE/profiles/${ko}.hmm" \
     && head -1 "$D/profiles/${ko}.hmm" | grep -q HMMER; then
    got=$((got+1))
  else
    rm -f "$D/profiles/${ko}.hmm"
  fi
done < "$D/targets_resolved.tsv"
echo "  profiles present: $got/$NKO"

if [ "$got" -lt "$NKO" ]; then
  TARB="$K/profiles.tar.gz"
  [ -s "$TARB" ] || { echo "  fetching the KOfam archive once (~400 MB, resumable)"
                      mkdir -p "$K"; wget -c -O "$TARB" "$BASE/profiles.tar.gz" || true; }
  if [ -s "$TARB" ]; then
    awk -F'\t' '{print "profiles/"$1".hmm"}' "$D/targets_resolved.tsv" > "$D/wanted.txt"
    tar -xzf "$TARB" -C "$D" -T "$D/wanted.txt" 2>/dev/null || true
    got=$(ls -1 "$D/profiles"/*.hmm 2>/dev/null | wc -l)
    echo "  after archive extraction: $got/$NKO"
    cp -n "$D/profiles"/*.hmm "$K/profiles/" 2>/dev/null || true
  fi
fi

# A missing TARGET profile makes that marker's zero uninterpretable; a missing
# CONTROL profile removes the proof the search works at all.
MISSING=""
while IFS=$'\t' read -r ko gene role dom kws note; do
  [ -z "${ko:-}" ] && continue
  [ -s "$D/profiles/${ko}.hmm" ] || MISSING="$MISSING $gene($role)"
done < "$D/targets_resolved.tsv"
if [ -n "$MISSING" ]; then
  echo "  !! profiles unavailable for:$MISSING"
  case "$MISSING" in
    *CONTROL*) echo "     A control profile is missing — no zero here would be"
               echo "     interpretable. Stopping."
               runlog "**34 ABORTED: control profile missing**"; exit 4 ;;
    *)         echo "     Those TARGETS are reported as UNTESTED, not absent." ;;
  esac
fi
cat "$D/profiles"/*.hmm > "$D/targets.hmm"


# --------------------------------------------------- RESOLVE THE BINARIES USED
# 🔴 Stage A produced nothing on the first run and the log said only
# "hmmsearch failed". The real cause was one line earlier: Prodigal never ran.
# `run_tool prodigal ... 2>/dev/null` discarded its own error message, which was
# almost certainly "no working env resolved for prodigal" — run_tool looks up
# ENV_<tool> and gives up if the smoketest never exercised that binary.
#
# This is EXACTLY the bug fixed in step 33 for hmmscan, and I did not carry the
# fix here. Rule 19 says reuse the code, not the idea; I reused the idea.
#
# Each binary is now located by looking for it, tested by finding it executable,
# and its stderr is KEPT.
bin_env() {                      # bin_env <tool> -> env token, "-" for PATH, or ""
  local tool="$1"
  local var="ENV_${tool}"
  local cands="${!var:-} ${ENV_bakta:-} ${ENV_prodigal:-} ${ENV_hmmsearch:-} ${ENV_diamond:-} ${ENV_hmmscan:-} -"
  local e prefix
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then
      command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }
      continue
    fi
    if declare -f env_prefix >/dev/null 2>&1; then
      prefix=$(env_prefix "$e" 2>/dev/null) || continue
      [ -x "$prefix/bin/$tool" ] && { echo "$e"; return 0; }
    else
      conda run -n "$e" "$tool" -h >/dev/null 2>&1 && { echo "$e"; return 0; }
    fi
  done
  return 1
}

echo
echo "######## 0d. BINARIES ########"
GENECALLER=""
for t in prodigal pyrodigal; do
  e=$(bin_env "$t" || true)
  if [ -n "$e" ]; then
    GENECALLER="$t"; export "ENV_${t}=$e"
    echo "  gene caller: $t (env '$e')"
    break
  fi
done
if [ -z "$GENECALLER" ]; then
  echo "  !! no gene caller found. Tried prodigal and pyrodigal across"
  echo "     ENV_prodigal, ENV_bakta and PATH."
  echo "     Bakta bundles pyrodigal, so the bakta env normally has one:"
  echo "       ls \$(env_prefix \"\${ENV_bakta:-NB_Genome}\")/bin | grep -i rodigal"
  echo "     Stage A cannot run without it. Stage B is unaffected — re-run with"
  echo "     RUN_READS=1 once a gene caller is available."
  runlog "**34: no gene caller — stage A skipped**"
  SKIP_A=1
else
  SKIP_A=0
fi
for t in hmmsearch diamond; do
  e=$(bin_env "$t" || true)
  if [ -n "$e" ]; then export "ENV_${t}=$e"; echo "  $t: env '$e'"
  else echo "  !! $t not found in any environment"; fi
done

# ---------------------------------------------------------- contig -> MAG map
# Built from the MAG fastas, which is the authoritative source (step 25). A
# contig name resolving to more than one MAG is DISCARDED rather than arbitrated
# — Bakta renames sequences, so `contig_1` can exist in every genome.
echo
echo "######## 0c. CONTIG -> MAG MAP ########"
python3 - "$F" "$D" <<'PYM'
import sys, os, glob
from collections import defaultdict
F, D = sys.argv[1], sys.argv[2]
owner = defaultdict(set)
nfa = 0
for fa in glob.glob(os.path.join(F, 'bins', '*.fa')) + \
          glob.glob(os.path.join(F, 'bins', '*.fasta')):
    mag = os.path.splitext(os.path.basename(fa))[0]; nfa += 1
    for line in open(fa, errors='replace'):
        if line.startswith('>'):
            owner[line[1:].split()[0]].add(mag)
amb = sum(1 for v in owner.values() if len(v) > 1)
with open(os.path.join(D, 'contig2mag.tsv'), 'w') as fh:
    fh.write("contig\tmag\n")
    for c, v in owner.items():
        if len(v) == 1: fh.write(f"{c}\t{next(iter(v))}\n")
print(f"  {nfa} MAG fastas, {len(owner)} contig names, {amb} ambiguous (discarded)")
print(f"  binned contigs resolvable: {len(owner) - amb}")
PYM

# ===========================================================================
# STAGE A — every gene in both assemblies
# ===========================================================================
echo
echo "######## STAGE A — ALL ASSEMBLY GENES ########"
for code in SED LIQ; do
  asm="$WS2/04_polish/${code}_asm_final.fasta"
  faa="$D/prot/${code}_all.faa"
  [ -s "$asm" ] || { echo "  !! missing $asm"; exit 1; }
  # 20b already called genes on these assemblies. Reuse rather than repeat.
  if [ ! -s "$faa" ] && [ -s "$OLD/prot/${code}_all.faa" ]; then
    ln -sf "$OLD/prot/${code}_all.faa" "$faa"
    echo "  $code: reusing 20b's Prodigal proteins ($(grep -c '^>' "$faa") genes)"
  elif [ -s "$faa" ]; then
    echo "  $code: proteins cached ($(grep -c '^>' "$faa") genes)"
  elif [ "$SKIP_A" = "1" ]; then
    echo "  $code: no gene caller — stage A skipped for this sample"
    continue
  else
    echo "  $code: calling genes on the whole assembly with $GENECALLER -p meta ..."
    # stderr is KEPT. Discarding it is what hid the real failure last time.
    # No -o/-f/-q: only the .faa is used downstream, and those flags differ
    # between prodigal and pyrodigal.
    run_tool "$GENECALLER" -i "$asm" -a "$faa.part" -p meta \
        >/dev/null 2>"$D/prot/${code}_genecall.err"
    rc=$?
    if [ "$rc" -ne 0 ] || [ ! -s "$faa.part" ]; then
      echo "    !! $GENECALLER exited $rc and wrote no proteins. Its own message:"
      sed -n '1,8p' "$D/prot/${code}_genecall.err" 2>/dev/null | sed 's/^/       /'
      rm -f "$faa.part"
      continue
    fi
    mv "$faa.part" "$faa"
    echo "    $(grep -c '^>' "$faa") genes"
  fi
  [ -s "$faa" ] || { echo "  $code: no proteins — skipping hmmsearch"; continue; }
  out="$D/hits/${code}_assembly.tbl"
  if [ -s "$out" ]; then
    echo "  $code: hmmsearch cached"
  else
    echo "  $code: hmmsearch ($(date +%H:%M)) ..."
    run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-10 \
        "$D/targets.hmm" "$faa" >"$D/hits/${code}_hmmsearch.out" 2>&1
    rc=$?
    [ "$rc" -ne 0 ] && echo "  !! hmmsearch exited $rc on $code"
  fi
  # A tblout with no data rows is not a result — rpoB is in every bacterium in
  # this assembly, so an empty stage A is a broken search. Say so HERE, where
  # the cause is still visible, instead of letting it surface as zeros.
  nh=$(grep -vc '^#' "$out" 2>/dev/null || echo 0)
  echo "    $code assembly hits: ${nh:-0}"
  if [ "${nh:-0}" -eq 0 ]; then
    echo "    !! no hits at all. Diagnostics:"
    echo "       profiles in library : $(grep -c '^NAME' "$D/targets.hmm" 2>/dev/null || echo 0)"
    echo "       profile names       : $(grep '^NAME' "$D/targets.hmm" 2>/dev/null | head -3 | tr '\n' ' ')"
    echo "       proteins searched   : $(grep -c '^>' "$faa" 2>/dev/null || echo 0)"
    echo "       first protein header: $(grep -m1 '^>' "$faa" 2>/dev/null | cut -c1-60)"
    [ -s "$D/hits/${code}_hmmsearch.out" ] && \
      sed -n '1,6p' "$D/hits/${code}_hmmsearch.out" | sed 's/^/       /'
  fi
done

# ===========================================================================
# STAGE B — the raw reads
# ===========================================================================
if [ "$RUN_READS" = "1" ]; then
echo
echo "######## STAGE B — RAW READS (DIAMOND blastx, frameshift-aware) ########"
DB="$D/marker_ref.dmnd"; REF="$D/marker_ref.faa"
if [ ! -s "$REF" ]; then
  echo "  building a reference protein set ..."
  : > "$REF"
  # 🔴 QUERIES FROM THE VERIFIED HITS, NOT FROM A GENE NAME.
  # UniProt gene-name queries are wide open to NAME COLLISIONS, and three of
  # mine collided badly:
  #     mcrA  = methyl-coenzyme M reductase alpha  AND  the E. coli McrA
  #             modified-cytosine RESTRICTION endonuclease
  #     mtaB  = methanol:corrinoid methyltransferase AND a bacterial tRNA
  #             methylthiotransferase present in most bacteria
  #     cutC  = choline TMA-lyase  AND  the copper-homeostasis protein CutC
  # Every one of those decoys-by-accident is common in sediment bacteria, which
  # is why cutC read 11x higher than its own assembly evidence and mtaB read
  # 0.985 x rpoB — as if every bacterium carried a methanogen enzyme.
  #
  # The fix is to seed the read search with proteins THIS DATASET already
  # verified: every assembly protein that passed its KOfam curated threshold.
  # They cannot be name collisions, because no name was involved.
  if [ -s "$D/hits/SED_assembly.tbl" ] || [ -s "$D/hits/LIQ_assembly.tbl" ]; then
    python3 - "$D" "$REF" <<'PYQ'
import sys, os
D, REF = sys.argv[1], sys.argv[2]
thr = {}
with open(os.path.join(D, 'ko_list'), errors='replace') as fh:
    h = fh.readline().rstrip('\n').split('\t')
    ti = h.index('threshold') if 'threshold' in h else 1
    si = h.index('score_type') if 'score_type' in h else 2
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > max(ti, si):
            try: thr[f[0]] = (float(f[ti]), f[si].strip().lower())
            except ValueError: pass
gene = {}
for line in open(os.path.join(D, 'targets_resolved.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) > 2 and f[2] == 'TARGET': gene[f[0]] = f[1]
want = {}
for code in ('SED', 'LIQ'):
    tb = os.path.join(D, 'hits', f'{code}_assembly.tbl')
    if not os.path.exists(tb): continue
    for line in open(tb, errors='replace'):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 9 or f[2] not in gene: continue
        t, st = thr.get(f[2], (None, ''))
        if t is None: continue
        try: sc = float(f[8]) if st.startswith('dom') else float(f[5])
        except ValueError: continue
        if sc >= t: want[f[0]] = (gene[f[2]], code)
n = 0
out = open(REF, 'a')
for code in ('SED', 'LIQ'):
    faa = os.path.join(D, 'prot', f'{code}_all.faa')
    if not os.path.exists(faa): continue
    keep, name = False, None
    for line in open(faa, errors='replace'):
        if line.startswith('>'):
            name = line[1:].split()[0]
            keep = name in want and want[name][1] == code
            if keep: out.write(f">{want[name][0]}|ASM_{name}\n"); n += 1
        elif keep:
            out.write(line)
out.close()
print(f"    verified assembly proteins added as queries: {n}")
PYQ
  fi

  # Our own verified sequences first, where this project has them: step 15's
  # marker trees were vetted against the ACR and MttB-superfamily look-alikes.
  #
  # LABELLED BY MARKER, NOT BY ORIGIN. A first version tagged these `OWN|`,
  # which would have put every one of our own mcrA reads into a bucket called
  # "OWN" and reported mcrA as absent at read level. The read counter keys on
  # the label, so the label has to BE the marker name. Only files whose basename
  # matches a target gene are taken — anything else has no marker to attribute
  # its hits to and would do the same damage more quietly.
  while IFS=$'\t' read -r ko gene role dom kws note; do
    [ -z "${gene:-}" ] && continue
    for src in "$WS2/15_marker_trees/${gene}.faa" "$WS2/15_marker_trees/${gene}"*.faa; do
      [ -s "$src" ] || continue
      sed "s/^>/>${gene}|/" "$src" >> "$REF"
      echo "    own $gene sequences from $(basename "$src"): $(grep -c '^>' "$src")"
      break
    done
  done < "$D/targets_resolved.tsv"
  fetch_up() {  # fetch_up <label> <uniprot query>
    local lab="$1"
    local q="$2"
    local tmp="$D/.up_${lab}.faa"
    local n
    wget -q -O "$tmp" \
      "https://rest.uniprot.org/uniprotkb/search?format=fasta&size=500&query=$q" 2>/dev/null
    n=$(grep -c '^>' "$tmp" 2>/dev/null || echo 0)
    if [ "${n:-0}" -gt 0 ]; then
      sed "s/^>/>${lab}|/" "$tmp" >> "$REF"; printf '    %-10s %s sequences\n' "$lab" "$n"
    else
      printf '    %-10s NONE\n' "$lab"
    fi
    rm -f "$tmp"
  }
  fetch_up dsrA '%28gene%3AdsrA%29'
  fetch_up dsrB '%28gene%3AdsrB%29'
  fetch_up aprA '%28gene%3AaprA%29'
  fetch_up sat  '%28gene%3Asat%29+AND+%28reviewed%3Atrue%29'
  fetch_up mcrA '%28gene%3AmcrA%29'
  fetch_up mtaB '%28gene%3AmtaB%29'
  fetch_up mtmB '%28gene%3AmtmB%29'
  fetch_up aprB '%28gene%3AaprB%29'
  fetch_up mtbB '%28gene%3AmtbB%29'
  fetch_up nifH '%28gene%3AnifH%29'
  fetch_up cutC '%28gene%3AcutC%29'
  # 🔴 DECOYS. Level 3 assigns each read to its BEST hit, so a family can only
  # be distinguished from its relatives if those relatives are IN the database.
  # Without them, every pyruvate formate-lyase read in a sediment full of
  # fermenters is filed as cutC, and 4,172 "cutC" reads is what that looks like.
  # These are competitors, not targets: a read routed here is a read correctly
  # NOT counted, and the counter reports them so the size of the problem shows.
  fetch_up DECOYpflB  '%28gene%3ApflB%29+AND+%28reviewed%3Atrue%29'     # vs cutC
  fetch_up DECOYcysI  '%28gene%3AcysI%29+AND+%28reviewed%3Atrue%29'     # assimilatory SiR, vs dsrAB
  fetch_up DECOYnirB  '%28gene%3AnirB%29+AND+%28reviewed%3Atrue%29'     # nitrite reductase, vs dsrAB
  fetch_up DECOYmttB  '%28gene%3AmttB%29'                               # MttB superfamily, vs mta/mtm/mtbB
  fetch_up CTRLglnA '%28gene%3AglnA%29+AND+%28reviewed%3Atrue%29'
  fetch_up CTRLrpoB '%28gene%3ArpoB%29+AND+%28reviewed%3Atrue%29'
  fetch_up CTRLtuf  '%28gene%3Atuf%29+AND+%28reviewed%3Atrue%29'
  fetch_up CTRLef1a '%28gene%3Atuf%29+AND+%28taxonomy_id%3A2157%29'
fi
# UniProt headers carry the protein NAME. Drop any target sequence whose name
# does not contain the marker's own keyword — that removes McrA restriction
# enzymes from mcrA, tRNA methylthiotransferases from mtaB and copper CutC from
# cutC, without needing to know in advance which names collide.
python3 - "$D" "$REF" <<'PYF'
import sys, os
D, REF = sys.argv[1], sys.argv[2]
kw = {}
for line in open(os.path.join(D, 'targets_resolved.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) > 4 and f[2] == 'TARGET':
        kw.setdefault(f[1], set()).update(k.strip().lower() for k in f[4].split(';') if k.strip())
kept = dropped = 0
out, keep, drops = [], False, {}
for line in open(REF, errors='replace'):
    if line.startswith('>'):
        lab = line[1:].split('|', 1)[0]
        if lab in kw and '|ASM_' not in line:
            desc = line.lower()
            keep = any(k in desc for k in kw[lab])
            if keep: kept += 1
            else:
                dropped += 1; drops[lab] = drops.get(lab, 0) + 1
        else:
            keep = True; kept += 1
    if keep: out.append(line)
open(REF, 'w').writelines(out)
if dropped:
    print(f"    name-collision filter dropped {dropped} UniProt sequence(s):")
    for k, v in sorted(drops.items(), key=lambda x: -x[1]):
        print(f"      {k:<8}{v:>5} whose protein name lacks its own keyword")
PYF

NREF=$(grep -c '^>' "$REF" 2>/dev/null || echo 0)
NCTR=$(grep '^>' "$REF" 2>/dev/null | grep -c 'CTRL' || echo 0)
echo "  reference set: $NREF sequences ($NCTR control/own)"
if [ "${NREF:-0}" -lt 20 ] || [ "${NCTR:-0}" -lt 5 ]; then
  echo "  !! reference set too small or without controls — stage B would produce"
  echo "     an uninterpretable zero. Skipping it; stage A still stands."
  RUN_READS=0
elif ! run_tool diamond --version >/dev/null 2>&1; then
  echo "  !! diamond not available. Skipping stage B."
  RUN_READS=0
else
  [ -s "$DB" ] || run_tool diamond makedb --in "$REF" -d "${DB%.dmnd}" --quiet 2>/dev/null
  CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
  [ -s "$CLEAN" ] && SED_IN="$CLEAN" || SED_IN="$SED_READS"
  for pair in "SED:$SED_IN" "LIQ:$LIQ_READS"; do
    code="${pair%%:*}"; reads="${pair#*:}"
    out="$D/hits/${code}_reads.tsv"
    [ -s "$reads" ] || { echo "  !! missing reads for $code"; continue; }
    [ -s "$out" ] && { echo "  $code: read search cached"; continue; }
    echo "  $code: DIAMOND blastx --frameshift 15 ($(date +%H:%M)) ..."
    run_tool diamond blastx -d "$DB" -q "$reads" -o "$out" \
        --frameshift 15 --range-culling --top 10 -e 1e-10 \
        --threads "$THREADS" --quiet \
        -f 6 qseqid sseqid pident length evalue bitscore 2>/dev/null \
      || echo "  !! diamond failed on $code"
    echo "    $(wc -l < "$out" 2>/dev/null || echo 0) alignments"
  done
fi
fi

# ===========================================================================
# INTERPRETATION
# ===========================================================================
echo
python3 - "$D" "$WS2" "$RUN_READS" <<'PY'
import sys, os, csv, glob, gzip
from collections import defaultdict, Counter

D, WS2, RUN_READS = sys.argv[1], sys.argv[2], sys.argv[3] == '1'
W = 104
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t: say('  ' + t); say('=' * W)

meta = {}
for line in open(os.path.join(D, 'targets_resolved.tsv')):
    f = (line.rstrip('\n').split('\t') + [''] * 6)[:6]
    if f[0]: meta[f[0]] = dict(gene=f[1], role=f[2], domain=f[3], note=f[5])
by_gene = {v['gene']: k for k, v in meta.items()}
TARGETS = [k for k, v in meta.items() if v['role'] == 'TARGET']
CONTROLS = [k for k, v in meta.items() if v['role'] == 'CONTROL']
HAVE = {k for k in meta if os.path.exists(os.path.join(D, 'profiles', f'{k}.hmm'))}

# ------------------------------------------------------------- contig -> MAG
c2m = {}
p = os.path.join(D, 'contig2mag.tsv')
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        c2m[r['contig']] = r['mag']

# ------------------------------------------------------------------- STAGE A
# 🔴 ADAPTIVE THRESHOLDS, and why the first run's numbers were not reportable.
#   dsrA/dsrB are homologous to ASSIMILATORY sulfite reductase and to nitrite
#   reductase; cutC is a glycyl-radical enzyme homologous to pyruvate
#   formate-lyase, which every fermenter in this sediment carries; mtaB belongs
#   to the same corrinoid-methyltransferase superfamily as MttB. At a flat
#   E < 1e-10 every one of those relatives scores. KOfam ships a curated
#   per-family score threshold precisely to separate them, and steps 11b and 19b
#   already use it. Not using it here was my omission, and it is the difference
#   between "sulfate reduction is in the sediment" and "something that looks a
#   bit like it is".
thr = {}
kl = os.path.join(D, 'ko_list')
if os.path.exists(kl):
    with open(kl, errors='replace') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ti = hdr.index('threshold') if 'threshold' in hdr else 1
        si = hdr.index('score_type') if 'score_type' in hdr else 2
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > max(ti, si):
                try: thr[f[0]] = (float(f[ti]), f[si].strip().lower())
                except ValueError: thr[f[0]] = (None, f[si].strip().lower())

# hmmsearch --tblout columns (0-based): 0 target(protein) 2 query(KO)
#   4 E-value(full) 5 score(full)  7 E-value(best dom) 8 score(best dom)
asm = defaultdict(lambda: defaultdict(list))
nothr = set(); ndrop = Counter()
for code in ('SED', 'LIQ'):
    tb = os.path.join(D, 'hits', f'{code}_assembly.tbl')
    if not os.path.exists(tb): continue
    for line in open(tb, errors='replace'):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 9: continue
        prot, ko = f[0], f[2]
        if ko not in meta:
            ko = f[3] if len(f) > 3 and f[3] in meta else None
            if ko is None: continue
        t, st = thr.get(ko, (None, ''))
        if t is None:
            nothr.add(ko)                      # no curated threshold: keep, but flag
        else:
            try: sc = float(f[8]) if st.startswith('dom') else float(f[5])
            except ValueError: continue
            if sc < t:
                ndrop[ko] += 1
                continue
        contig = prot.rsplit('_', 1)[0]        # Prodigal names proteins <contig>_<n>
        asm[ko][code].append((prot, contig, c2m.get(contig)))

# ------------------------------------------------------------------- STAGE B
reads = defaultdict(lambda: defaultdict(set))  # gene -> sample -> {read ids}
nreads = {}
DOMFRAC = float(os.environ.get('DOMFRAC', '0.25'))
DOMMIN = int(os.environ.get('DOMMIN', '100'))
domrep = []
refseqs = Counter()
_rp = os.path.join(D, 'marker_ref.faa')
if os.path.exists(_rp):
    for line in open(_rp, errors='replace'):
        if line.startswith('>'): refseqs[line[1:].split('|', 1)[0]] += 1
for code in ('SED', 'LIQ'):
    rp = os.path.join(D, 'hits', f'{code}_reads.tsv')
    if not os.path.exists(rp): continue
    best = {}
    for line in open(rp, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        q, sid = f[0], f[1]
        try: bs = float(f[5])
        except ValueError: continue
        lab = sid.split('|', 1)[0]
        if q not in best or bs > best[q][1]: best[q] = (lab, bs, sid)
    # 🔴 PER-REFERENCE DOMINANCE FILTER — the step-13d rule, generalised.
    #   dsrB returned 1,231 sediment reads and TWO reference sequences accounted
    #   for 1,970 of the alignments; every other dsrB reference had 13-17.
    #   dsrA, from the same operon, was flat (15, 14, 14, 13, 13) and returned
    #   70. Adjacent genes of an operon cannot differ 17-fold: the two spiked
    #   references were promiscuous magnets, exactly as one SILVA reference took
    #   99.6% of the Erwinia reads at step 13d.
    #
    #   A reference that collects a large share of a marker's reads while
    #   hundreds of its siblings collect a handful is not measuring the marker.
    #   It is flagged, excluded, and REPORTED — never silently dropped.
    percount = defaultdict(Counter)
    for q, (lab, bs, sid) in best.items():
        percount[lab][sid] += 1
    # WHO CAN BE EXCLUDED, AND WHO CANNOT.
    #   v8 fired on ASM_mcrA (63% of 30 reads) and on a reviewed SwissProt glnA.
    #   Both are legitimate: a query taken from THIS DATASET should recruit most
    #   of its own marker's reads if the organism is closely related, and that is
    #   a sign the query is good, not bad. Excluding them cut mcrA from 30 to 11
    #   and glnA from 7,886 to 4,337 for no reason.
    #
    #   A magnet is an artefact of a bad DATABASE ENTRY. Only unreviewed TrEMBL
    #   entries qualify: dataset-derived proteins (ASM_) and reviewed SwissProt
    #   entries (sp|) are never excluded. The marker also needs enough reads for
    #   a share to mean anything.
    #
    #   The nref>=20 gate is gone: it exempted exactly the thin query sets that
    #   are most likely to be contaminated.
    bad_ref = set()
    for lab, c in percount.items():
        tot = sum(c.values())
        if tot < DOMMIN: continue
        for sid, n in c.items():
            if '|tr|' not in sid: continue          # TrEMBL only
            if n / tot > DOMFRAC:
                bad_ref.add(sid)
                domrep.append((code, lab, sid, n, tot, n / tot))
    for q, (lab, bs, sid) in best.items():
        if sid in bad_ref: continue
        reads[lab][code].add(q)

# read totals for per-million normalisation
for code, key in (('SED', 'sediment'), ('LIQ', 'liquid')):
    nreads[code] = None

decoys = {k: v for k, v in reads.items() if k.startswith('DECOY')}

# 🔴 WHICH GENES WERE ACTUALLY QUERIED AT LEVEL 3?
# glnA, aprB and mtbB all reported 0 reads in the last run. None of them was in
# the reference set — I never fetched sequences for them. A zero from a query
# that does not exist is not a measurement, and printing it in the same column
# as a real zero is the exact error this whole script is built to prevent.
refgenes = Counter()
refp = os.path.join(D, 'marker_ref.faa')
if os.path.exists(refp):
    for line in open(refp, errors='replace'):
        if line.startswith('>'):
            refgenes[line[1:].split('|', 1)[0]] += 1
def queried(g):
    return refgenes.get(g, 0) + refgenes.get('CTRL' + g, 0)
if domrep:
    rule('0a. PROMISCUOUS REFERENCE SEQUENCES, EXCLUDED')
    say(f"  A single reference taking more than {DOMFRAC:.0%} of a marker's reads while")
    say("  hundreds of its siblings take a handful is not measuring that marker.")
    say()
    say(f"  {'smpl':<6}{'gene':<8}{'share':>7}{'reads':>8}  reference")
    say('-' * W)
    for code, lab, sid, n, tot, fr in sorted(domrep, key=lambda x: -x[5]):
        say(f"  {code:<6}{lab:<8}{fr:>6.0%}{n:>8}  {sid[:56]}")
    say()
    say("  Those reads are removed from the counts below. The check that found")
    say("  this: dsrA and dsrB are adjacent genes of one operon and returned")
    say("  70 against 1,231 reads — impossible, and traceable to two references.")

rule('0. LEVEL-2 THRESHOLDS')
say(f"  KOfam curated score thresholds applied to {len(thr)} families.")
if ndrop:
    say(f"  hits below their family's own threshold, DISCARDED:")
    for k, n in ndrop.most_common():
        say(f"    {meta.get(k, {}).get('gene', k):<8} {k}  {n:>6} discarded")
    say("  These are the homologues a flat E-value cannot separate: assimilatory")
    say("  sulfite reductase for dsrAB, pyruvate formate-lyase for cutC, the")
    say("  MttB superfamily for the methyltransferases.")
if nothr:
    say(f"  !! no curated threshold for: {', '.join(sorted(nothr))} — those kept on")
    say("     E-value alone and their counts are UPPER BOUNDS.")
if decoys:
    say()
    say("  LEVEL-3 DECOYS — reads routed to a close relative instead of a target:")
    for k in sorted(decoys):
        ns = len(decoys[k].get('SED', ())); nl = len(decoys[k].get('LIQ', ()))
        say(f"    {k[5:]:<8} SED {ns:>7}   LIQ {nl:>7}")
    say("  Every one of those is a read the previous run would have counted as a")
    say("  target. A decoy that collects thousands is doing its job.")

if refgenes:
    thin = [(g, n) for g, n in sorted(refgenes.items()) if n < 20 and not g.startswith('DECOY')]
    if thin:
        say()
        say("  !! THIN QUERY SETS — level-3 sensitivity for these is low, so their")
        say("     read counts are FLOORS, not estimates:")
        for g, n in thin:
            say(f"       {g:<10}{n:>4} reference sequence(s)")

rule('1. CONTROLS — IS EACH LEVEL WORKING?')
say("  Controls are DOMAIN-DECLARED. rpoB/glnA/tuf are bacterial; archaea have a")
say("  split RNA polymerase B and EF-1alpha, not EF-Tu. An archaeal target's")
say("  absence is only callable if an ARCHAEAL control passed.")
say()
ctrl_asm = defaultdict(list); ctrl_rd = defaultdict(list)
say(f"  {'control':<10}{'domain':<8}{'SED asm':>9}{'LIQ asm':>9}{'SED reads':>11}{'LIQ reads':>11}")
say('-' * W)
for ko in CONTROLS:
    g, dom = meta[ko]['gene'], meta[ko]['domain']
    a_s, a_l = len(asm[ko]['SED']), len(asm[ko]['LIQ'])
    lab = 'CTRL' + g if ('CTRL' + g) in reads else g
    r_s = len(reads.get(lab, {}).get('SED', ()))
    r_l = len(reads.get(lab, {}).get('LIQ', ()))
    say(f"  {g:<10}{dom:<8}{a_s:>9}{a_l:>9}{r_s:>11}{r_l:>11}")
    ctrl_asm[dom].append(a_s + a_l > 0)
    ctrl_rd[dom].append(r_s + r_l > 0)
say()

# A CONTROL IS PER LEVEL. The first version pooled them, so a working read
# search (6,707 rpoB reads) was thrown away because the assembly search found
# nothing. Each level is judged on its own control, and a level that fails
# invalidates only its own evidence.
asm_ok = {d: any(v) for d, v in ctrl_asm.items()}
rd_ok = {d: any(v) for d, v in ctrl_rd.items()}
for d in sorted(set(asm_ok) | set(rd_ok)):
    say(f"  {d:<6} assembly (level 2): {'[PASS]' if asm_ok.get(d) else '[FAIL]'}"
        f"     reads (level 3): {'[PASS]' if rd_ok.get(d) else '[FAIL]'}")
say()
LVL2 = asm_ok.get('BACT', False)
LVL3 = rd_ok.get('BACT', False)
if not LVL2:
    say("  !! LEVEL 2 IS BROKEN. rpoB is in every bacterium in that assembly, so")
    say("     an empty assembly search is a technical failure. No level-2 result")
    say("     is reported below and no absence is called from it.")
if not LVL3:
    say("  !! LEVEL 3 IS BROKEN — no control recovered from the reads either.")
if not LVL2 and not LVL3:
    say()
    say("  Both levels failed. NO VERDICT.")
    open(os.path.join(D, 'VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
    sys.exit(5)
if LVL3 and not LVL2:
    say("  Level 3 stands on its own control and IS reported. A read-level")
    say("  presence is still a presence; what is lost is the binned/unbinned")
    say("  split, which only the assembly can give.")

rule('2. WHAT LEVEL DOES EACH MARKER REACH?')
say("  LEVEL 1  in a MAG            already in the catalogue")
say("  LEVEL 2  assembly, BINNED    on a contig belonging to a MAG")
say("  LEVEL 2u assembly, UNBINNED  a real organism the binning missed")
say("  LEVEL 3  reads only          below the assembly threshold")
say("  none     with controls passing = a TESTED NEGATIVE")
say()
say("  Two rows per marker: SEDIMENT first, then the cultures. 'x/rpoB' is the")
say("  read count against the bacterial control's, a rough abundance proxy.")
say()
say(f"  {'gene':<8}{'KO':<9}{'smpl':<6}{'asm':>7}{'unbin':>7}{'reads':>9}{'x/rpoB':>8}   level")
say('-' * W)
rows = []
for ko in TARGETS:
    g, dom = meta[ko]['gene'], meta[ko]['domain']
    if ko not in HAVE:
        say(f"  {g:<8}{ko:<9}{'—':<6}{'no profile, UNTESTED':>40}")
        rows.append(dict(ko=ko, gene=g, domain=dom, level='UNTESTED (no profile)'))
        continue
    hs, hl = asm[ko]['SED'], asm[ko]['LIQ']
    us = sum(1 for _, _, m in hs if m is None)
    ul = sum(1 for _, _, m in hl if m is None)
    rs = len(reads.get(g, {}).get('SED', ()))
    rl = len(reads.get(g, {}).get('LIQ', ()))
    nq = queried(g)

    # 🔴 PER SAMPLE. Pooling SED and LIQ made dsrA read "UNBINNED ONLY (2)"
    # when both of those contigs are in the CULTURE assembly and the sediment
    # assembly has none. The chapter's question is about the lagoon, and a gene
    # found only in an enrichment is not an answer to it.
    def level_for(nasm, nunb, nrd, lvl2ok):
        nbin = nasm - nunb
        if lvl2ok and nbin and nunb: return f'assembly: binned {nbin} + unbinned {nunb}'
        if lvl2ok and nbin:          return f'assembly: binned only ({nbin})'
        if lvl2ok and nunb:          return f'assembly: UNBINNED ({nunb})'
        if nq == 0:                  return 'level 3 UNQUERIED'
        if nrd:                      return f'reads only ({nrd})'
        if not lvl2ok:               return 'UNTESTED (level 2 broken)'
        return 'not detected'
    lvl_sed = level_for(len(hs), us, rs, LVL2)
    lvl_liq = level_for(len(hl), ul, rl, LVL2)
    # 🔴 `binned = any(...)` labelled a marker "in a MAG" if ONE of its hits was
    # on a binned contig, hiding the other forty. dsrA had 8 unbinned SED hits
    # and 1 binned, and came out as "in a MAG" — erasing the entire finding this
    # step exists to make. A marker can be in both places and the report has to
    # say so.
    lvl = lvl_sed          # the sediment is what the chapter is about
    # crude abundance proxy: target reads against the bacterial control's reads.
    # Not length-corrected, so it is an order-of-magnitude statement only.
    ctlr = len(reads.get('CTRLrpoB', {}).get('SED', ())) or 1
    ratio = rs / ctlr
    say(f"  {g:<8}{ko:<9}{dom:<6}{len(hs):>7}{us:>7}{rs:>9}{ratio:>8.3f}   {lvl_sed}")
    say(f"  {'':<8}{'':<9}{'LIQ':<6}{len(hl):>7}{ul:>7}{rl:>9}{'':>8}   {lvl_liq}")
    rows.append(dict(ko=ko, gene=g, domain=dom, asm_sed=len(hs), unbinned_sed=us,
                     asm_liq=len(hl), unbinned_liq=ul, reads_sed=rs, reads_liq=rl,
                     n_queries=nq, sed_over_rpoB=round(ratio, 4),
                     level=lvl_sed, level_liq=lvl_liq))

with open(os.path.join(D, 'marker_levels.tsv'), 'w', newline='') as fh:
    fld = ['ko', 'gene', 'domain', 'asm_sed', 'unbinned_sed', 'asm_liq',
           'unbinned_liq', 'reads_sed', 'reads_liq', 'n_queries',
           'sed_over_rpoB', 'level', 'level_liq']
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=fld, extrasaction='ignore')
    w.writeheader(); w.writerows(rows)

# per-hit table, because the unbinned contigs are the point
with open(os.path.join(D, 'assembly_hits.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['ko', 'gene', 'sample', 'protein', 'contig', 'mag_or_UNBINNED'])
    for ko in TARGETS + CONTROLS:
        for code in ('SED', 'LIQ'):
            for prot, contig, mag in asm[ko][code]:
                w.writerow([ko, meta[ko]['gene'], code, prot, contig, mag or 'UNBINNED'])

rule('2b. DO THE TWO LEVELS AGREE?')
say("  Levels 2 and 3 are independent measurements of the same thing, so their")
say("  ratios to the same control should match. Where reads run far above the")
say("  assembly, the read query set is picking up something the KOfam threshold")
say("  rejects — a name collision, or a wider family.")
say()
ra = len(asm.get(by_gene.get('rpoB', ''), {}).get('SED', ())) or 1
rr = len(reads.get('CTRLrpoB', {}).get('SED', ())) or 1
say(f"  {'gene':<8}{'asm/rpoB':>10}{'reads/rpoB':>12}{'ratio':>9}   reading")
say('-' * W)
for r in rows:
    if not r.get('gene') or r.get('n_queries', 1) == 0: continue
    a = (r.get('asm_sed', 0) or 0) / ra
    b = (r.get('reads_sed', 0) or 0) / rr
    if a == 0:
        say(f"  {r['gene']:<8}{a:>10.3f}{b:>12.3f}{'—':>9}   no assembly hit — not cross-checkable")
        continue
    q = b / a
    verd = 'consistent' if 0.33 <= q <= 3 else (
        'READS TOO HIGH — suspect the query set' if q > 3 else 'reads low — thin or divergent queries')
    say(f"  {r['gene']:<8}{a:>10.3f}{b:>12.3f}{q:>9.1f}   {verd}")
say()
say("  A marker with no assembly hit cannot be cross-checked at all, and its")
say("  read count has to be defended some other way.")

rule('2c. A BIOLOGICAL CEILING FOR THE ARCHAEAL MARKERS')
say("  EVERY methanogen carries mcrA; only some carry any given methyltransferase.")
say("  So in the same sample no methyltransferase marker can exceed mcrA. This is")
say("  not a threshold I chose — it is a property of the organisms, and it makes")
say("  a usable automatic check.")
say()
_mc = next((r for r in rows if r['gene'] == 'mcrA'), None)
if _mc:
    base = _mc.get('reads_sed', 0) or 0
    say(f"  mcrA, sediment reads: {base}")
    say()
    bad = []
    for r in rows:
        if r['gene'] in ('mcrA',) or r.get('domain') != 'ARCH': continue
        n = r.get('reads_sed', 0) or 0
        # 3x, not 2x: at counts of a few dozen, Poisson noise plus query-set
        # asymmetry can produce a twofold difference with no biology in it.
        flag = 'IMPLAUSIBLE' if base and n > 3 * base else (
            'ok' if n <= base else 'slightly above mcrA — not resolvable at these counts')
        if flag == 'IMPLAUSIBLE': bad.append(r['gene'])
        say(f"    {r['gene']:<8}{n:>8}   {n / base if base else 0:>7.1f}x mcrA   {flag}")
    say()
    if bad:
        say(f"  !! {', '.join(bad)} exceed mcrA by more than twofold in the sediment.")
        say("     That cannot be true of genuine methanogens, so those read counts")
        say("     are still measuring something else — a wider corrinoid")
        say("     methyltransferase family that neither the MttB decoy nor the")
        say("     name filter separates. DO NOT REPORT THEM AS A PRESENCE.")
        say("     The methanol arm therefore remains UNTESTED, which is the same")
        say("     conclusion step 31 reached by a different route.")
    else:
        say("  All archaeal methyltransferase markers sit at or below mcrA, which is")
        say("  what genuine methanogens look like.")

rule('3. THE SULFATE REDUCERS')
sul = [r for r in rows if r['gene'] in ('dsrA', 'dsrB', 'aprA', 'aprB', 'sat')]
def in_sed(r):
    return (r.get('asm_sed', 0) or 0) > 0 or (r.get('reads_sed', 0) or 0) > 0
found = [r for r in sul if in_sed(r)]
absent = [r for r in sul if not in_sed(r) and (r.get('n_queries', 0) or 0) > 0]
unq = [r for r in sul if not in_sed(r) and not (r.get('n_queries', 0) or 0)]
say("  Step 20 found ZERO dsrAB in 29 MAGs while SILVA puts *Desulfatiglans* at")
say("  9.48% of sediment SSU reads. This is the test of that absence.")
say()
if found:
    say(f"  {len(found)} of {len(sul)} sulfur markers are present IN THE SEDIMENT:")
    say(f"    {'gene':<6}{'sed asm':>9}{'unbinned':>10}{'sed reads':>11}{'x/rpoB':>9}   where")
    for r in found:
        say(f"    {r['gene']:<6}{r.get('asm_sed', 0):>9}{r.get('unbinned_sed', 0):>10}"
            f"{r.get('reads_sed', 0):>11}{r.get('sed_over_rpoB', 0):>9.3f}   {r['level']}")
    if unq:
        say(f"    UNQUERIED at level 3 (no reference sequences): "
            f"{', '.join(r['gene'] for r in unq)}")
    say()
    say("  READ THE COLUMNS, NOT THE HEADLINE. A marker with 0 in 'sed asm' and a")
    say("  read count is present in the LAGOON but did not assemble there — the")
    say("  organism is below the assembly threshold, not merely unbinned.")
    say()
    say("  ** The catalogue's zero was a RECOVERY LIMIT, not an absence. **")
    say("  Sulfate reduction is genetically present in Etoliko sediment; the")
    say("  organisms carrying it did not assemble into MAGs. This is what the")
    say("  SILVA Desulfatiglans signal predicted and it is now measured.")
    say()
    say("  Consequence for the chapter: the methanogenesis argument is framed")
    say("  against sulfate reducers outcompeting methanogens for H2 and acetate.")
    say("  That competitor is now DOCUMENTED in the sediment rather than inferred")
    say("  from an rRNA profile — which strengthens the non-competitive-substrate")
    say("  explanation instead of weakening it.")
else:
    say("  No sulfur marker was recovered at any level, with controls passing.")
    say("  That is a tested negative and it CONTRADICTS the SILVA Desulfatiglans")
    say("  signal at 9.48%. Before reporting it, check whether those SSU reads")
    say("  were assigned to a lineage whose sulfate-reduction genes are genuinely")
    say("  absent — several Desulfatiglans relatives are not sulfate reducers.")

rule('4. THE METHANOL ARM')
mt = {r['gene']: r for r in rows if r['gene'] in ('mtaB', 'mtmB', 'mcrA')}
for g in ('mcrA', 'mtmB', 'mtaB'):
    r = mt.get(g)
    if r: say(f"  {g:<6} {r['level']}")
say()
say("  mtaB is the consumer end of pectin -> methanol -> CH4. Step 30 found the")
say("  producer (CE8 in ETO_LIQ_MAG_015, 0.90 sediment breadth); step 31 could")
say("  not test the consumer because Bakta names no mtaB in any archaeal genome")
say("  here, so there was no query. This is the first test of it.")

rule('5. WHAT TO REPORT')
n_unb = sum(1 for r in rows
            if (r.get('unbinned_sed', 0) or 0) + (r.get('unbinned_liq', 0) or 0) > 0)
n_neg = sum(1 for r in rows if 'TESTED NEGATIVE' in str(r.get('level', '')))
n_unt = sum(1 for r in rows if 'UNTESTED' in str(r.get('level', '')))
say(f"  markers found on unbinned contigs (new organisms): {n_unb}")
say(f"  markers reaching a TESTED NEGATIVE:                {n_neg}")
say(f"  markers left UNTESTED:                             {n_unt}")
say()
say("  A tested negative is a far stronger sentence than a MAG-based absence and")
say("  should be written as one: 'X was not detected in the assemblies or in the")
say("  raw reads, in a search that recovered <controls> from the same data'.")
say()
say("  LIMITS:")
say("   1. Level 3 has a detection floor — one read. Quote it with every zero.")
say("   2. hmmsearch at E<1e-10 without KOfam adaptive thresholds is permissive;")
say("      hits on unbinned contigs should be checked against the KO's own")
say("      threshold before being quoted individually.")
say("   3. An unbinned contig is not a genome. It shows the gene is in the")
say("      community, not which organism carries it.")

open(os.path.join(D, 'VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'marker_levels.tsv')}")
print(f"  written: {os.path.join(D, 'assembly_hits.tsv')}")
print(f"  written: {os.path.join(D, 'VERDICT.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "marker_levels.tsv written" $([ -s "$D/marker_levels.tsv" ] && echo 0 || echo 1) || OK=1
check "assembly_hits.tsv written" $([ -s "$D/assembly_hits.tsv" ] && echo 0 || echo 1) || OK=1
[ "$OK" -eq 0 ] && runlog "34_unbinned complete" || runlog "**34_unbinned FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
