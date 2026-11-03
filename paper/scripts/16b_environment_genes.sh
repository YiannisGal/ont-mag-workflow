#!/usr/bin/env bash
# =============================================================================
# 16b_environment_genes.sh  --  SURVIVING ETOLIKO: sulfide, salt, metals, oxygen
#
# WHY 16b AND NOT 16
#   `16_environment_genes.sh` exists and was never run to completion — the
#   notebook has carried a placeholder for it since Volume II. It was written
#   before this project learned three things the hard way, and re-running it
#   would reproduce all three:
#
#     - it matches KEYWORDS against Bakta product strings. Step 30 replaced that
#       method for CAZymes and step 34's failure 8 showed why: a gene NAME is
#       not a gene, and product strings collide across unrelated proteins.
#     - it predates `00_config_patch.sh`, so its `run_tool` calls resolve
#       `ENV_<tool>` variables that are only set for tools the smoketest probed.
#       That is what broke steps 33, 34 and 35.
#     - it has no domain-declared control, and steps 19b and 23e both once read
#       a guaranteed absence as a result.
#
#   16b asks the same questions with KOfam profiles at KEGG's curated
#   thresholds, the method steps 34/36/36b converged on. 16 is left un-renamed
#   as the provenance of the question, exactly as 32 was left for its BAMs.
#
# THE QUESTION, AND WHY IT IS NOT GENERIC
#   Etoliko is permanently stratified, anoxic, SULFIDIC and sulfate-replete.
#   Four stressor axes, each with a specific expectation:
#
#   1. SULFIDE. H2S at millimolar concentrations poisons haem-copper oxidases
#      and attacks Fe-S clusters. Organisms living in it carry sulfide:quinone
#      oxidoreductase (SQR) and persulfide dioxygenase (ETHE1) to oxidise it,
#      rhodanese-family sulfurtransferases to shuttle it, and cytochrome bd
#      (CydAB) — the terminal oxidase that keeps working under sulfide when
#      cytochrome c oxidase does not.
#   2. SALT. A coastal lagoon with a halocline. Compatible-solute synthesis
#      (ectoine, glycine betaine, trehalose) and the ABC/BCCT uptake systems.
#      Step 29 already tested the betaine TRANSPORTERS; this tests SYNTHESIS,
#      which is the different question of whether they can make their own.
#   3. METALS. Sulfidic sediment concentrates As and Hg and mobilises them on
#      redox swings. ars operon, merA, and the P-type ATPase efflux pumps.
#   4. OXYGEN. And here is the trap this step exists to avoid:
#
# 🔴 THE ABSENCE-OF-THE-WRONG-FAMILY TRAP, WRITTEN IN BEFORE IT HAPPENS
#   Step 36 reported zero B12 and zero thiamine ABC transporters across 29
#   genomes and the answer was that Firmicutes use ECF transporters instead —
#   Rule 35, an absence is an absence of what you searched for. The identical
#   trap is waiting here: **strict anaerobes largely do not use superoxide
#   dismutase and catalase.** They use superoxide REDUCTASE (SOR /
#   desulfoferrodoxin), rubrerythrin and rubredoxin — a non-haem, non-radical-
#   producing route. Searching SOD/KatE alone in a catalogue of anaerobic
#   Clostridia and methanogens would produce a table of zeros and a completely
#   wrong conclusion.
#
#   So both families are searched, and the report is explicitly conditional on
#   which one is found. That is the whole design of SECTION 6.
#
#   RESULT, recorded here so the header is not left predicting: the anaerobe
#   route leads (SOR 23/29, rubrerythrin 18/29) but classical defence is
#   widespread too (SOD 20/29) -- 18/29 carry BOTH and 0/29 carry neither. The
#   prediction that the anaerobe route would dominate starkly was half right,
#   and searching only SOD/catalase would still have understated the catalogue.
#
# CONTROLS
#   1. POSITIVE, DOMAIN-DECLARED: tuf (K02358) in >=80% of Bacteria, EF-1alpha
#      (K03231) in >=80% of Archaea. Archaea have no EF-Tu. A failure here
#      blocks every verdict.
#   2. LABEL GUARD: every KO's expected identity is checked against KEGG's own
#      ko_list definition BEFORE anything is searched. Step 34 lost six versions
#      to a KO that was silently the wrong gene.
#   3. NOT SEARCHED != ABSENT: a KO with no profile is carried through every
#      table as NOT SEARCHED and can never be reported as a zero. This defect
#      has now appeared in steps 36, 34b and 34b again (Rule 34).
#   4. THRESHOLDS: hits clearing E<1e-10 but failing the curated threshold are
#      counted, not discarded silently. Step 36 measured that at 2,895%; this
#      step measures 511% on its own target set.
#   5. THE JOIN IS DETERMINISTIC. Genus/domain come from a row-local one-hop
#      propagation through sorted() ids, after two identical runs of an earlier
#      version returned DIFFERENT genera for the same genome (Rule 46).
#
# THE COMPARISON WORTH THE STEP
#   Sediment-derived vs enrichment-derived genomes, and in-situ (step 26 breadth
#   >= 0.5) vs not. If sulfide-handling genes concentrate in the genomes actually
#   confirmed in the lagoon, that is a selection signal measured within one
#   catalogue by one method. Fisher exact, two-sided. It can come back null.
#
# OUTPUT  16_environment/env_matrix.tsv · ENVIRONMENT.txt · not_searched.txt
#         guard.tsv (the KO list AS RESOLVED, including AUTO rows)
# RUNTIME 10-25 min first run, seconds when cached (29 proteomes x 31 profiles)
# USAGE   ws2 && bash scripts/16b_environment_genes.sh 2>&1 | tee logs/16b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="16b_environment"
A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"
D="$WS2/16_environment"; mkdir -p "$D/profiles" "$D/hits"

echo "######## 16b — ENVIRONMENTAL STRESS GENES, BY PROFILE ########"
echo "  started: $(date -Is)"

bin_env() {
  local tool="$1" p
  if command -v "$tool" >/dev/null 2>&1; then echo ""; return 0; fi
  for root in ${CONDA_ROOTS:-${CONDA_ENVS_SHARED} /opt/conda/envs}; do
    for p in "$root"/*/bin/"$tool"; do
      [ -x "$p" ] && { echo "${p%/bin/$tool}"; return 0; }
    done
  done
  return 1
}
HS_ENV="$(bin_env hmmsearch)" || { echo "  !! no hmmsearch anywhere. Stop."; exit 2; }
HMMSEARCH="${HS_ENV:+$HS_ENV/bin/}hmmsearch"
echo "  hmmsearch: ${HS_ENV:-<on PATH>}"

KOLIST=""
for c in "$K/ko_list" "$K/db/ko_list" "$WS2/11_kofam/ko_list"; do
  [ -s "$c" ] && { KOLIST="$c"; break; }
done
[ -n "$KOLIST" ] || { echo "  !! ko_list not found under $K — curated thresholds are the method here."; exit 2; }
echo "  ko_list  : $KOLIST  ($(( $(wc -l < "$KOLIST") - 1 )) families)"

PROFDIR=""
for c in "$K/profiles" "$K/db/profiles" "$WS2/11_kofam/profiles"; do
  [ -d "$c" ] && { PROFDIR="$c"; break; }
done
echo "  profiles : ${PROFDIR:-<none>}  ($(ls -1 "$PROFDIR" 2>/dev/null | wc -l) files)"

# ------------------------------------------------------------ the targets
# Asserted BY KO NUMBER with the expected identity beside it. The guard checks
# each against KEGG's own definition and refuses to run on a mismatch, because
# a name can drift and step 34 spent six versions on a KO that was the wrong
# gene. Columns: KO | short | axis | family | expect (regex, case-insensitive)
cat > "$D/targets.tsv" <<'KOS'
K17218	sqr	SULFIDE	OXIDATION	sulfide:?quinone
K17725	ethe1	SULFIDE	OXIDATION	persulfide dioxygenase|sulfur dioxygenase
K01011	tst/sseA	SULFIDE	TRANSFER	thiosulfate sulfurtransferase|rhodanese|3-mercaptopyruvate
K00425	cydA	SULFIDE	TOLERANT_OXIDASE	cytochrome (bd|d) ubiquinol oxidase
K00426	cydB	SULFIDE	TOLERANT_OXIDASE	cytochrome (bd|d) ubiquinol oxidase
K08738	cycM	SULFIDE	CTRL_OXIDASE	cytochrome c
K06718	ectA	SALT	ECTOINE	diaminobutyr(ate|ic).*acetyltransferase
K00836	ectB	SALT	ECTOINE	diaminobutyr(ate|ic).*(transaminase|aminotransferase)
K06720	ectC	SALT	ECTOINE	ectoine synthase
K00108	betA	SALT	BETAINE_SYN	choline dehydrogenase
K00130	betB	SALT	BETAINE_SYN	betaine.*aldehyde dehydrogenase
K00697	otsA	SALT	TREHALOSE	trehalose.*phosphate synthase
K01087	otsB	SALT	TREHALOSE	trehalose.*phosphatase
K00147	proA	SALT	PROLINE	glutamate.*5.*semialdehyde dehydrogenase
K00537	arsC	METAL	ARSENIC	arsenate reductase
K03325	arsB	METAL	ARSENIC	arsenite|arsenical
K03741	arsC2	METAL	ARSENIC	arsenate reductase
K00520	merA	METAL	MERCURY	mercuric reductase
K01534	zntA	METAL	ZN_CD_PB	Zn2?\+?.*(Cd|transport|export)|cadmium|zinc.*(export|ATPase)
K17686	copA	METAL	COPPER	Cu\+?.*(export|transport|ATPase)|copper.*(export|transport|ATPase)
K15726	czcA/cusA	METAL	RND_EFFLUX	heavy metal efflux|cobalt.*zinc.*cadmium|CzcA|CusA
K04564	sodA	OXYGEN	AEROBE_ROUTE	superoxide dismutase
K03781	katE	OXYGEN	AEROBE_ROUTE	catalase
K05919	sor	OXYGEN	ANAEROBE_ROUTE	superoxide reductase|desulfoferrodoxin|neelaredoxin
AUTO	rbr	OXYGEN	ANAEROBE_ROUTE	rubrerythrin
K03386	ahpC/prx	OXYGEN	PEROXIREDOXIN	peroxiredoxin|alkyl hydroperoxide reductase
K03671	trxA	OXYGEN	THIOL_REDOX	thioredoxin
K00384	trxB	OXYGEN	THIOL_REDOX	thioredoxin reductase
K03564	bcp	OXYGEN	PEROXIREDOXIN	peroxiredoxin|thiol peroxidase
K02358	tuf	CONTROL	BACTERIA	elongation factor Tu
K03231	ef1a	CONTROL	ARCHAEA	elongation factor 1
KOS
NT=$(( $(wc -l < "$D/targets.tsv") ))
echo "  targets  : $NT KOs across 4 stressor axes + 2 domain controls"

echo
echo "######## LABEL GUARD ########"
python3 - "$D" "$KOLIST" <<'PY'
import sys, os, re
D, KOLIST = sys.argv[1], sys.argv[2]
defn = {}
with open(KOLIST, errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    try:
        di = hdr.index('definition')
    except ValueError:
        di = -1
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > max(1, di):
            defn[f[0]] = f[di] if di >= 0 else f[-1]
def candidates(exp, limit=8):
    """every KO whose KEGG definition matches the expectation"""
    try:
        rx = re.compile(exp, re.I)
    except re.error:
        return []
    return [(k, v) for k, v in defn.items() if rx.search(v)][:limit + 1]

bad = 0
rows = []
for line in open(os.path.join(D, 'targets.tsv'), errors='replace'):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5:
        continue
    ko, short, axis, fam, exp = f[:5]

    # ---- AUTO: resolve the KO from ko_list by its definition ---------------
    # Used where asserting a number failed. It prints what it resolved to,
    # because a resolution nobody reads is how step 34 spent six versions on
    # the corrinoid protein (Rule 24); and finding nothing is recorded as a
    # statement about the SEARCH STRING, not about the genomes (Rule 25).
    if ko == 'AUTO':
        cs = candidates(exp)
        if not cs:
            print(f"  [AUTO ] {short:<10} /{exp}/ matches no KEGG definition —")
            print(f"          recorded as NOT SEARCHED. That is a statement about")
            print(f"          the search string, not about these genomes.")
            rows.append(('AUTO_' + short, short, axis, fam, 'NOKO'))
            continue
        if len(cs) > 8:
            print(f"  [FAIL] {short:<10} /{exp}/ matches >8 KOs — too broad to resolve")
            bad += 1
            continue
        # all matches are searched. Ambiguity is handled by searching the whole
        # family, never by silently picking one of them.
        for k, v in cs:
            print(f"  [AUTO ] {k}  {short:<10} {v[:60]}")
            rows.append((k, short, axis, fam, 'OK'))
        continue

    d = defn.get(ko)
    if d is None:
        print(f"  [ABSENT] {ko}  {short:<10} not in ko_list — cannot be searched")
        rows.append((ko, short, axis, fam, 'NOKO'))
        continue
    if re.search(exp, d, re.I):
        print(f"  [ok]   {ko}  {short:<10} {d[:62]}")
        rows.append((ko, short, axis, fam, 'OK'))
    else:
        print(f"  [FAIL] {ko}  {short:<10} expected /{exp}/")
        print(f"         KEGG says: {d[:80]}")
        # v1 stopped here and made you re-run once per wrong KO. Print the KOs
        # that DO match the expectation, so the fix is one edit, not a guess.
        cs = candidates(exp, limit=5)
        if cs:
            print(f"         KOs whose definition DOES match /{exp}/:")
            for k, v in cs[:5]:
                print(f"           {k}  {v[:66]}")
        else:
            print(f"         ...and no KO matches /{exp}/ either, so the")
            print(f"            expectation string is the thing to fix.")
        rows.append((ko, short, axis, fam, 'BAD'))
        bad += 1
with open(os.path.join(D, 'guard.tsv'), 'w') as oh:
    for r in rows:
        oh.write('\t'.join(r) + '\n')
if bad:
    print(f"\n  !! {bad} KO label(s) disagree with KEGG. NOTHING IS RUN.")
    print("     TWO things can cause this and they need different fixes:")
    print("       (a) the KO number is wrong  -> find the right KO")
    print("       (b) the expectation string is too narrow for KEGG's wording")
    print("           -> widen the regex in targets.tsv")
    print("     Read the printed definition to tell which. The guard cannot, and")
    print("     it aborts either way rather than running on an unverified label:")
    print("     a KO that was silently the wrong gene cost step 34 six versions.")
    sys.exit(9)
print("\n  [PASS] every KO's identity agrees with KEGG's own definition.")
PY
[ $? -eq 0 ] || { runlog "**16b_environment guard FAILED**"; exit 9; }

# ------------------------------------------------------------ profiles
echo
echo "######## PROFILES ########"
: > "$D/not_searched.txt"
NEED=(); HAVE=0
while IFS=$'\t' read -r ko short axis fam st; do
  [ "$st" = "NOKO" ] && { echo -e "$ko\t$short\tnot in ko_list" >> "$D/not_searched.txt"; continue; }
  p="$D/profiles/$ko.hmm"
  if [ -s "$p" ]; then HAVE=$((HAVE+1)); continue; fi
  for src in "$PROFDIR/$ko.hmm" "$PROFDIR/$ko.HMM" "$WS2/34b_methylamine/profiles/$ko.hmm" \
             "$WS2/36_gapclose/profiles/$ko.hmm" "$WS2/34_unbinned/profiles/$ko.hmm"; do
    [ -s "$src" ] && { cp "$src" "$p"; break; }
  done
  if [ -s "$p" ]; then HAVE=$((HAVE+1)); else NEED+=("$ko"); fi
done < "$D/guard.tsv"

if [ "${#NEED[@]}" -gt 0 ]; then
  # v1 tried to fetch each KO from KEGG and got 404 twenty-six times. That URL
  # does not serve individual profiles -- KOfam is distributed as one tarball,
  # which is already on this machine (1.46 GB, downloaded for step 36). So:
  # extract from the local archive, and only then consider the network. An
  # unrun search that 404s is still an unrun search, but there is no reason for
  # it to be unrun when the file is on disk.
  echo "  ${#NEED[@]} profile(s) not in $PROFDIR — looking for the KOfam archive"
  TARB=""
  for c in "$K/profiles.tar.gz" "$K/db/profiles.tar.gz" "$WS2/profiles.tar.gz" \
           "$K/ko_profiles.tar.gz"; do
    [ -s "$c" ] && { TARB="$c"; break; }
  done
  if [ -z "$TARB" ]; then
    for r in "$WS2" "${DBDIR:-${DB_DIR}}" "$HOME" ${SCRATCH_ROOT}; do
      [ -d "$r" ] || continue
      TARB="$(find "$r" -maxdepth 6 -type f \( -name 'profiles.tar.gz' \
             -o -name 'kofam*profiles*.tar.gz' -o -name 'ko_profiles*.tar.gz' \) \
             -size +100M 2>/dev/null | head -1)"
      [ -n "$TARB" ] && break
    done
  fi
  if [ -n "$TARB" ] && [ -s "$TARB" ]; then
    echo "  archive: $TARB"
    # extract only the ones needed; the archive nests them under profiles/
    PATS=(); for ko in "${NEED[@]}"; do PATS+=("*/$ko.hmm" "$ko.hmm"); done
    tar -xzf "$TARB" -C "$D/profiles" --strip-components=1 --wildcards \
        "${PATS[@]}" 2>/dev/null
    tar -xzf "$TARB" -C "$D/profiles" --wildcards "${PATS[@]}" 2>/dev/null
    find "$D/profiles" -mindepth 2 -name '*.hmm' -exec mv -n {} "$D/profiles"/ \; 2>/dev/null
    STILL=()
    for ko in "${NEED[@]}"; do
      if [ -s "$D/profiles/$ko.hmm" ]; then
        HAVE=$((HAVE+1)); echo "    $ko  extracted"
      else
        STILL+=("$ko"); echo "    $ko  NOT IN ARCHIVE"
      fi
    done
    NEED=("${STILL[@]:-}")
    [ -z "${NEED[0]:-}" ] && NEED=()
  else
    echo "  !! no local KOfam archive found. Searched \$WS2, \$DBDIR, \$HOME,"
    echo "     ${SCRATCH_ROOT} for profiles.tar.gz (>100 MB)."
    echo "     If it is elsewhere, point KOFAM_TAR at it and re-run:"
    echo "       KOFAM_TAR=/path/to/profiles.tar.gz bash scripts/16b_environment_genes.sh"
  fi
  # explicit override, checked first on a re-run
  if [ -n "${KOFAM_TAR:-}" ] && [ -s "${KOFAM_TAR}" ] && [ "${#NEED[@]}" -gt 0 ]; then
    echo "  using KOFAM_TAR=$KOFAM_TAR"
    for ko in "${NEED[@]}"; do
      tar -xzf "$KOFAM_TAR" -C "$D/profiles" --strip-components=1 --wildcards \
          "*/$ko.hmm" 2>/dev/null || tar -xzf "$KOFAM_TAR" -C "$D/profiles" \
          --wildcards "$ko.hmm" 2>/dev/null
      [ -s "$D/profiles/$ko.hmm" ] && HAVE=$((HAVE+1))
    done
  fi
  # whatever is still missing is NOT SEARCHED, recorded by KO
  for ko in "${NEED[@]:-}"; do
    [ -n "$ko" ] || continue
    [ -s "$D/profiles/$ko.hmm" ] && continue
    echo -e "$ko\t-\tno profile on disk or in archive" >> "$D/not_searched.txt"
  done
fi
echo "  profiles available: $HAVE/$NT"
NS=$(wc -l < "$D/not_searched.txt")
[ "$NS" -gt 0 ] && echo "  NOT SEARCHED: $NS  (never reported as absence — see $D/not_searched.txt)"

# the domain controls must be searchable or nothing can be validated
for ko in K02358 K03231; do
  [ -s "$D/profiles/$ko.hmm" ] || { echo "  !! control profile $ko missing. Stop."; exit 3; }
done
if [ "$HAVE" -lt $((NT / 2)) ]; then
  echo "  !! fewer than half the profiles are present. Step 36 v1 manufactured 190"
  echo "     fictional downgrades from exactly this state. Stopping instead."
  exit 3
fi

cat "$D/profiles"/*.hmm > "$D/all.hmm" 2>/dev/null
PSTAMP=$(md5sum "$D/all.hmm" | cut -c1-12)
echo "  profile-set stamp: $PSTAMP"

# ------------------------------------------------------------ search
echo
echo "######## SEARCHING ########"
mapfile -t MAGS < <(ls -1 "$A" 2>/dev/null)
n=0
for m in "${MAGS[@]}"; do
  faa="$A/$m/$m.faa"
  [ -s "$faa" ] || faa="$(ls -1 "$A/$m"/*.faa 2>/dev/null | head -1)"
  [ -n "$faa" ] && [ -s "$faa" ] || continue
  n=$((n+1))
  out="$D/hits/$m.tbl"
  if [ -s "$out" ] && [ "$(cat "$D/hits/$m.stamp" 2>/dev/null)" = "$PSTAMP" ]; then
    printf "\r  [%d/%d] %-34s cached" "$n" "${#MAGS[@]}" "$m"; continue
  fi
  printf "\r  [%d/%d] %-34s        " "$n" "${#MAGS[@]}" "$m"
  "$HMMSEARCH" --cpu "${THREADS:-8}" -E 1e-10 --tblout "$out.part" \
      "$D/all.hmm" "$faa" >/dev/null 2>&1
  if [ -s "$out.part" ]; then
    mv "$out.part" "$out"; echo "$PSTAMP" > "$D/hits/$m.stamp"
  else
    rm -f "$out.part"
  fi
done
echo
echo "  tables: $(ls -1 "$D/hits"/*.tbl 2>/dev/null | wc -l)/$n"

echo
rm -f "$D/ENVIRONMENT.txt" "$D/env_matrix.tsv"
python3 - "$D" "$A" "$F" "$KOLIST" "$WS2" <<'PY'
import sys, os, re, csv, glob, math
from collections import defaultdict, Counter

D, A, F, KOLIST, WS2 = sys.argv[1:6]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)
STATUS = ['INCOMPLETE']

def dump(code=0):
    with open(os.path.join(D, 'ENVIRONMENT.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'ENVIRONMENT.txt')}")
    sys.exit(code)

def _lf(n):
    return math.lgamma(n + 1)

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, r2, c1 = a + b, c + d, a + c
    lo, hi = max(0, c1 - r2), min(r1, c1)
    def lp(x):
        return (_lf(r1) + _lf(r2) + _lf(c1) + _lf(n - c1)
                - _lf(n) - _lf(x) - _lf(r1 - x) - _lf(c1 - x) - _lf(r2 - c1 + x))
    ps = {x: math.exp(lp(x)) for x in range(lo, hi + 1)}
    tot = sum(ps.values()); obs = ps.get(a, 0.0)
    return min(1.0, sum(v for v in ps.values() if v <= obs * (1 + 1e-9)) / tot) if tot else 1.0

# ------------------------------------------------------- targets & thresholds
TGT = {}
ORDER = []
# guard.tsv, not targets.tsv: AUTO rows are resolved to real KO numbers by the
# guard, and reading the asserted list here would lose that resolution.
for line in open(os.path.join(D, 'guard.tsv'), errors='replace'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 4 and f[0] not in TGT:
        TGT[f[0]] = (f[1], f[2], f[3])
        ORDER.append(f[0])

NOTSEARCHED = {}
p = os.path.join(D, 'not_searched.txt')
if os.path.exists(p):
    for line in open(p, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if f and f[0]:
            NOTSEARCHED[f[0]] = f[-1]
SEARCHED = [k for k in ORDER if k not in NOTSEARCHED]

thr = {}
with open(KOLIST, errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    ti = hdr.index('threshold') if 'threshold' in hdr else 1
    si = hdr.index('score_type') if 'score_type' in hdr else 2
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[0] in TGT:
            try:
                thr[f[0]] = (float(f[ti]), f[si] if len(f) > si else '')
            except (ValueError, IndexError):
                thr[f[0]] = (None, '')

# ------------------------------------------------------------ read hits
hits = defaultdict(set)     # ko -> {mag}
nprot = defaultdict(int)
sub = Counter()
nothr = set()
mags = []
for fn in sorted(os.listdir(os.path.join(D, 'hits'))):
    if not fn.endswith('.tbl'):
        continue
    mag = re.sub(r'\.tbl$', '', fn)
    mags.append(mag)
    for line in open(os.path.join(D, 'hits', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 9:
            continue
        ko = f[2] if f[2] in TGT else (f[3] if len(f) > 3 and f[3] in TGT else None)
        if ko is None:
            continue
        t, st = thr.get(ko, (None, ''))
        if t is None:
            nothr.add(ko)
            hits[ko].add(mag); nprot[ko] += 1
            continue
        try:
            sc = float(f[8]) if st.startswith('dom') else float(f[5])
        except ValueError:
            continue
        if sc >= t:
            hits[ko].add(mag); nprot[ko] += 1
        else:
            sub[ko] += 1

# ------------------------------------------------------ taxonomy / breadth
# 🔴 v1's join was NON-DETERMINISTIC. Its token filter accepted any string of
# >=4 characters, so 'Bacteria', 'Methanosarcina' and '99.85' all became alias
# keys; a GLOBAL alias merge then linked unrelated genomes through them, and the
# genus was picked by iterating an unordered set. Two identical runs of 16b gave
# LIQ_das_metabat_001 -> Methanosarcina and -> JBAIYS01. A join whose answer
# changes between runs of the same script on the same data is not a join.
#
# Fixed three ways: a genome id in this project always contains an underscore
# AND a digit (ETO_LIQ_MAG_002, LIQ_das_metabat_075) while genus and domain
# names do not (Lentimicrobium, GCA-2746605, Methanobacterium_C, JBAIYS01);
# propagation is ONE HOP and ROW-LOCAL so nothing merges transitively; and every
# choice goes through sorted() so it cannot depend on hash order.
IDTOK = re.compile(r'^(?=.*_)(?=.*\d)[A-Za-z0-9._-]{6,}$')
RANKP = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

def id_tokens(r):
    out = set()
    for v in r.values():
        v = str(v or '').strip()
        if v and ' ' not in v and ';' not in v and not v.startswith(RANKP) \
                and IDTOK.match(v) and not re.fullmatch(r'[0-9._-]+', v):
            out.add(v)
    return out

def _rows(path):
    if not os.path.exists(path):
        return []
    with open(path, errors='replace') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

dom, gen, alias = {}, {}, {}
INTAX = set()
for r in _rows(os.path.join(F, 'MAG_TAXONOMY.tsv')):
    blob = ' '.join(str(v or '') for v in r.values())
    d = next((str(r[c]).strip() for c in r if c and c.strip().lower() == 'domain'), '')
    if not d:
        m = re.search(r'd__([A-Za-z]+)', blob); d = m.group(1) if m else ''
    g = next((str(r[c]).strip() for c in r if c and c.strip().lower() == 'genus'), '')
    if not g:
        m = re.search(r'g__([^;\t]+)', blob); g = m.group(1).strip() if m else ''
    toks = id_tokens(r)
    INTAX |= toks
    for i in toks:
        if d:
            dom[i] = d
        if g:
            gen[i] = g
AMBIG = 0
for r in _rows(os.path.join(F, 'MAG_CROSSWALK.tsv')):
    ids = sorted(id_tokens(r))
    if len(ids) < 2:
        continue
    for i in ids:
        alias[i] = set(ids)          # row-local, never merged across rows
    gs = sorted({gen[i] for i in ids if i in gen})
    ds = sorted({dom[i] for i in ids if i in dom})
    if len(gs) > 1:
        AMBIG += 1
    if any(i in INTAX for i in ids):
        INTAX |= set(ids)
    for i in ids:
        if ds:
            dom.setdefault(i, ds[0])
        if gs:
            gen.setdefault(i, gs[0])

# ---- step 26's breadth, located rather than guessed -----------------------
# v1 globbed for '*breadth*'/'*crossmap*' and found nothing, because step 26's
# files are cross_abundance.tsv, SED_vs_all.tsv and LIQ_vs_all.tsv. v2 scanned
# columns instead and picked 12_results/MAG_CATALOGUE.tsv -- whose breadth_SED
# column only carries the FOUR sediment-derived MAGs. The contrast then ran as
# "4 SED MAGs vs 25 LIQ MAGs", which is a comparison of ASSEMBLY ORIGIN wearing
# the label of in-situ presence. A table that covers a sixth of the catalogue
# cannot answer this question, so coverage is now checked and a thin table is
# refused.
#
# Also: CoverM names its column "Covered Fraction" WITH A SPACE, and its genome
# ids carry a .fa extension the crosswalk records but the proteome names do not.
# Both are handled rather than assumed away.
def _norm(x):
    return re.sub(r'\.(fa|fna|fasta|fa\.gz|fna\.gz)$', '', str(x or '').strip())

def find_breadth(ws2, alias, need):
    files = sorted(glob.glob(os.path.join(ws2, '26_crossmap', '*.tsv'))) + \
            sorted(glob.glob(os.path.join(ws2, '12_results', '*.tsv')))
    files.sort(key=lambda f: (0 if re.search(r'sed', os.path.basename(f), re.I) else 1,
                              os.path.basename(f)))
    # normalised alias index, so LIQ_metabat_001.fa finds LIQ_das_metabat_001
    ax = defaultdict(set)
    for k, grp in alias.items():
        ax[_norm(k)].update(grp)
    seen, best = [], None
    for c2 in files:
        try:
            rows = list(csv.DictReader(open(c2, errors='replace'), delimiter='\t'))
        except Exception:
            continue
        if not rows:
            continue
        cols = [x for x in rows[0].keys() if x]
        bcol = next((x for x in cols if re.search(r'breadth', x, re.I)
                     and re.search(r'sed', x, re.I)), None) \
            or next((x for x in cols if re.search(r'breadth|covered[ _-]?frac|'
                                                  r'coverage[ _-]?frac|frac[ _-]?cov', x, re.I)), None)
        seen.append((os.path.relpath(c2, ws2), cols, bcol is not None))
        if not bcol:
            continue
        idc = next((x for x in cols if re.search(r'mag|genome|bin|^name$', x, re.I)), cols[0])
        out = {}
        for r in rows:
            mid = _norm(r.get(idc))
            try:
                v = float(str(r.get(bcol) or '').strip() or 'nan')
            except ValueError:
                continue
            if v != v:
                continue
            for a2 in sorted(ax.get(mid) or {mid}):
                out[a2] = v
        cov = len([m for m in need if m in out])
        if best is None or cov > best[2]:
            best = (out, (os.path.relpath(c2, ws2), bcol, idc), cov)
        if cov >= 0.9 * len(need):
            break
    if best is None:
        return {}, None, seen, 0
    return best[0], best[1], seen, best[2]

breadth, bsrc, bseen, BCOV = find_breadth(WS2, alias, mags)

# ============================================================ 1. CONTROL
rule('1. CONTROL — CAN THIS SEARCH FIND WHAT MUST BE THERE?')
bact = [m for m in mags if dom.get(m, '').lower().startswith('bacteri')]
arch = [m for m in mags if dom.get(m, '').lower().startswith('archae')]
unk = [m for m in mags if m not in bact and m not in arch]
say(f"  {len(bact)} bacterial, {len(arch)} archaeal, {len(unk)} undeclared proteomes")
nogen = [m for m in mags if m not in gen]
say(f"  genus joined for {len(mags) - len(nogen)}/{len(mags)} proteomes; "
    f"{AMBIG} crosswalk row(s) carried conflicting genera")
if nogen:
    nojoin = sorted(m for m in nogen if m not in INTAX)
    nogtdb = sorted(m for m in nogen if m in INTAX)
    if nogtdb:
        say(f"  no GTDB genus: {', '.join(nogtdb)}")
        say("     These join fine; GTDB assigned them no genus (g__ is empty —")
        say("     'novel below family'). That is a RESULT, not a join failure, and")
        say("     it is why the pre-fix join had to invent one from a neighbour.")
    if nojoin:
        say(f"  NOT JOINED: {', '.join(nojoin)}")
        say("     Absent from MAG_CROSSWALK.tsv or carrying only one id token.")
if AMBIG:
    say("  !! a crosswalk row maps one genome to more than one genus. The genus")
    say("     column below is unreliable until that is resolved.")
say()
say("  Domain-declared, because ARCHAEA HAVE NO EF-Tu. Steps 19b and 23e both")
say("  applied a bacterial control to archaeal genomes and read the guaranteed")
say("  absence as a result.")
say()
CTRL_OK = True
for ko, label, pool in (('K02358', 'tuf (Bacteria)', bact),
                        ('K03231', 'EF-1a (Archaea)', arch)):
    if not pool:
        say(f"  {ko} {label:<18}: no genomes in this domain — not applicable")
        continue
    n = len([m for m in pool if m in hits.get(ko, ())])
    frac = n / len(pool)
    ok = frac >= 0.8
    CTRL_OK = CTRL_OK and ok
    say(f"  {ko} {label:<18}: {n}/{len(pool)} genomes ({100 * frac:.0f}%)"
        f"  [{'PASS' if ok else 'FAIL'}]")
say()
if not CTRL_OK:
    say("  [FAIL] the search cannot recover a universal single-copy gene in every")
    say("         domain present. NO VERDICT — every absence below would be")
    say("         uninterpretable.")
    STATUS[0] = 'CONTROL_FAILED'
    dump(0)
say("  [PASS] absences below are informative.")

# ============================================================ 2. thresholds
rule('2. THRESHOLDS — WHAT A FLAT E-VALUE WOULD HAVE ADDED')
kept = sum(nprot[k] for k in nprot)
rej = sum(sub.values())
say(f"  passed the curated KOfam threshold : {kept:,}")
say(f"  cleared E<1e-10 but FAILED it      : {rej:,}")
if kept:
    say(f"  a flat E<1e-10 would have inflated this table by {100 * rej / kept:.0f}%")
if nothr:
    say(f"  !! no curated threshold for: {', '.join(sorted(nothr))} — those were")
    say("     kept on E-value alone and are flagged in the matrix.")
say()
if NOTSEARCHED:
    say(f"  NOT SEARCHED ({len(NOTSEARCHED)}): "
        f"{', '.join(sorted(NOTSEARCHED))}")
    say("  These appear as NOT SEARCHED everywhere below and NEVER as zeros.")
    say()

# ============================================================ 3. by axis
NG = len(mags)
def line_for(ko):
    short, axis, fam = TGT[ko]
    if ko in NOTSEARCHED:
        return f"  {ko:<8}{short:<12}{fam:<18}{'NOT SEARCHED':>14}{'':>9}"
    n = len(hits.get(ko, ()))
    return (f"  {ko:<8}{short:<12}{fam:<18}{n:>14}"
            f"{100 * n / max(1, NG):>8.0f}%")

for axis, title, note in (
    ('SULFIDE', '3. SULFIDE — the defining stressor of this lagoon',
     "H2S poisons haem-copper oxidases and attacks Fe-S clusters. SQR and ETHE1\n"
     "  oxidise it; rhodanese-family enzymes shuttle it; cytochrome bd keeps\n"
     "  working under sulfide where cytochrome c oxidase does not."),
    ('SALT', '4. SALT — a coastal lagoon with a halocline',
     "Step 29 tested the betaine TRANSPORTERS. This tests SYNTHESIS, which is the\n"
     "  different question of whether an organism can make its own osmolyte or\n"
     "  must scavenge one."),
    ('METAL', '5. METALS — sulfidic sediment concentrates As and Hg',
     "Metal resistance is also a proxy for exposure history, since these operons\n"
     "  are readily lost when the pressure is absent."),
):
    rule(title)
    for l in note.split('\n'):
        say('  ' + l.strip() if not l.startswith('  ') else l)
    say()
    say(f"  {'KO':<8}{'gene':<12}{'role':<18}{'genomes':>14}{'%':>9}")
    say('-' * W)
    for ko in ORDER:
        if TGT[ko][1] == axis:
            say(line_for(ko))
    say()

# ============================================================ 6. OXYGEN
rule('6. OXIDATIVE STRESS — and the absence-of-the-wrong-family trap')
say("  Strict anaerobes largely do NOT use superoxide dismutase and catalase.")
say("  They use superoxide REDUCTASE (SOR / desulfoferrodoxin), rubrerythrin and")
say("  rubredoxin — a non-haem route that produces no O2 and no hydroxyl radical.")
say("  Searching SOD/KatE alone in a catalogue of anaerobic Clostridia and")
say("  methanogens would return zeros and support a completely wrong conclusion.")
say("  That is Rule 35, and it is why both families are here.")
say()
say(f"  {'KO':<8}{'gene':<12}{'route':<18}{'genomes':>14}{'%':>9}")
say('-' * W)
for ko in ORDER:
    if TGT[ko][1] == 'OXYGEN':
        say(line_for(ko))
say()
AER = [k for k in ORDER if TGT[k][2] == 'AEROBE_ROUTE' and k not in NOTSEARCHED]
ANA = [k for k in ORDER if TGT[k][2] == 'ANAEROBE_ROUTE' and k not in NOTSEARCHED]
aer_g = set().union(*[hits.get(k, set()) for k in AER]) if AER else set()
ana_g = set().union(*[hits.get(k, set()) for k in ANA]) if ANA else set()
say(f"  genomes with the AEROBE route (SOD or catalase)      : {len(aer_g)}/{NG}")
say(f"  genomes with the ANAEROBE route (SOR or rubrerythrin): {len(ana_g)}/{NG}")
say(f"  genomes with BOTH                                    : {len(aer_g & ana_g)}/{NG}")
say(f"  genomes with NEITHER                                 : "
    f"{NG - len(aer_g | ana_g)}/{NG}")
say()
if not AER or not ANA:
    say("  >> One of the two families was NOT SEARCHED, so this comparison did not")
    say("     run. Do not read either number as a route being absent.")
elif len(ana_g) > len(aer_g):
    say("  >> The anaerobe route dominates. Had this step searched SOD and catalase")
    say("     alone, it would have reported oxidative-stress defence as largely")
    say("     ABSENT from an anoxic-sediment catalogue — the step-36 vitamin-")
    say("     transporter error, in a different pathway.")
elif len(aer_g) > len(ana_g):
    say("  >> The aerobe route dominates, which is worth a sentence in a catalogue")
    say("     from permanently anoxic sediment: these organisms retain classical")
    say("     ROS defence, and enrichment cultures see more O2 than the lagoon does.")
else:
    say("  >> Neither route dominates. Report both and claim neither.")
say()

# ============================================================ 7. the contrast
rule('7. DOES SULFIDE HANDLING TRACK BEING IN THE LAGOON?')
say("  The comparison worth the step: genomes CONFIRMED in the sediment (step 26")
say("  breadth >= 0.50) against the rest of the catalogue. Same tool, same")
say("  thresholds, same run — an internally controlled contrast.")
say()
SULF = [k for k in ORDER if TGT[k][1] == 'SULFIDE' and TGT[k][2] != 'CTRL_OXIDASE'
        and k not in NOTSEARCHED]
if not breadth:
    say("  !! no breadth column found, so this contrast did not run. Recorded as")
    say("     NOT LOOKED UP, not as a null result. Headers actually present:")
    for f2, cols in bseen[:8]:
        say(f"       {f2}: {', '.join(cols[:12])}")
elif BCOV < 0.6 * len(mags):
    say(f"  !! the best breadth table covers only {BCOV}/{len(mags)} genomes"
        f" ({bsrc[0] if bsrc else '-'}).")
    say("     Running the contrast on it would split the catalogue by ASSEMBLY")
    say("     ORIGIN — the 4 sediment-derived MAGs against the 25 culture-derived")
    say("     ones — and label that split 'in situ vs not'. Those are different")
    say("     questions and the second is the one step 26 exists to answer.")
    say("     NOT RUN. Tables scanned, and whether each had a breadth column:")
    for f2, cols, hasb in bseen[:8]:
        say(f"       {'[b] ' if hasb else '[ ] '}{f2}: {', '.join(cols[:10])}")
elif not SULF:
    say("  !! no sulfide KO was searched. Nothing to contrast.")
else:
    ins = [m for m in mags if breadth.get(m, 0.0) >= 0.5]
    out = [m for m in mags if m not in ins]
    say(f"  source: {bsrc[0]}   breadth column: {bsrc[1]}   id column: {bsrc[2]}")
    say(f"  breadth values for {BCOV}/{len(mags)} genomes")
    say(f"  in situ (breadth >= 0.50): {len(ins)}   other: {len(out)}")
    say(f"  (breadth found for {len([m for m in mags if m in breadth])}/{NG} genomes)")
    say()
    if len(ins) < 3 or len(out) < 3:
        say("  Too few genomes on one side to test. Reporting counts only.")
    say(f"  {'KO':<8}{'gene':<12}{'in situ':>10}{'other':>10}{'fold':>9}{'p':>12}")
    say('-' * W)
    for ko in SULF:
        h = hits.get(ko, set())
        a = len([m for m in ins if m in h]); b = len(ins) - a
        c = len([m for m in out if m in h]); d = len(out) - c
        fa = a / max(1, len(ins)); fb = c / max(1, len(out))
        fold = (fa / fb) if fb else (float('inf') if fa else float('nan'))
        pv = fisher_two_sided(a, b, c, d)
        fs = ('    inf' if fold == float('inf')
              else ('      -' if fold != fold else f"{fold:>9.2f}"))
        say(f"  {ko:<8}{TGT[ko][0]:<12}{a:>10}{c:>10}{fs}{pv:>12.2e}")
    say()
    anyh = set().union(*[hits.get(k, set()) for k in SULF])
    a = len([m for m in ins if m in anyh]); b = len(ins) - a
    c = len([m for m in out if m in anyh]); d = len(out) - c
    pv = fisher_two_sided(a, b, c, d)
    say(f"  ANY sulfide-handling gene: {a}/{len(ins)} in situ vs {c}/{len(out)} other,"
        f" p = {pv:.2e}")
    if pv < 0.05:
        say("  >> The difference is significant. Note what it is NOT: 29 MAGs is a")
        say("     small sample and these genomes are not independent draws — they")
        say("     share phylogeny. Report it as an association within this")
        say("     catalogue, never as selection demonstrated.")
    else:
        say("  >> No significant difference. With 29 genomes this test has little")
        say("     power, so this is 'not detected here', not 'no difference'.")
say()

# ============================================================ per genome
rule('8. PER GENOME')
say(f"  {'genome':<26}{'genus':<20}{'sulf':>6}{'salt':>6}{'metal':>7}"
    f"{'O2-an':>7}{'O2-ae':>7}{'breadth':>9}")
say('-' * W)
def cnt(m, axis):
    return len([k for k in ORDER if TGT[k][1] == axis and k not in NOTSEARCHED
                and m in hits.get(k, ())])
for m in sorted(mags, key=lambda x: (-cnt(x, 'SULFIDE'), x)):
    b = breadth.get(m)
    say(f"  {m:<26}{(gen.get(m, '-') or '-')[:19]:<20}"
        f"{cnt(m, 'SULFIDE'):>6}{cnt(m, 'SALT'):>6}{cnt(m, 'METAL'):>7}"
        f"{len([k for k in ANA if m in hits.get(k, ())]):>7}"
        f"{len([k for k in AER if m in hits.get(k, ())]):>7}"
        f"{(f'{b:.2f}' if b is not None else '-'):>9}")

with open(os.path.join(D, 'env_matrix.tsv'), 'w') as oh:
    oh.write('mag\tgenus\tdomain\tbreadth\t' + '\t'.join(ORDER) + '\n')
    for m in sorted(mags):
        b = breadth.get(m)
        cells = []
        for ko in ORDER:
            cells.append('NOT_SEARCHED' if ko in NOTSEARCHED
                         else ('1' if m in hits.get(ko, ()) else '0'))
        oh.write(f"{m}\t{gen.get(m, '')}\t{dom.get(m, '')}\t"
                 f"{'' if b is None else f'{b:.3f}'}\t" + '\t'.join(cells) + '\n')

STATUS[0] = 'OK' if CTRL_OK else 'CONTROL_FAILED'
rule('WHAT TO REPORT')
say("  1. This closes the step-16 placeholder that has been in the notebook since")
say("     Volume II — and it closes it by PROFILE SEARCH, not by the keyword")
say("     matching the original script used. Say which method produced the")
say("     numbers, because the chapter has both kinds in it.")
say("  2. Section 6 is the methodological point worth a Methods sentence: the")
say("     anaerobe ROS route had to be searched alongside the aerobe one, or an")
say("     anoxic-sediment catalogue would have been reported as defenceless.")
_n = [3]
def item(*lines):
    say(f"  {_n[0]}. {lines[0]}")
    for l in lines[1:]:
        say(f"     {l}")
    _n[0] += 1
if NOTSEARCHED:
    item(f"{len(NOTSEARCHED)} KO(s) had no profile and were NOT SEARCHED. They are",
         "absent from every count above as a category, not as a zero.")
item("Presence is gene presence. None of this is expression, and metal",
     "resistance operons in particular are readily lost and readily acquired,",
     "so a genome without one is not thereby sensitive.")
item("Section 7's contrast is an association inside one catalogue of 29",
     "genomes that share phylogeny. It is not a demonstration of selection.")

dump(0)
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "report written" "$([ -s "$D/ENVIRONMENT.txt" ] && echo 0 || echo 1)" || OK=1
check "matrix written" "$([ -s "$D/env_matrix.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/ENVIRONMENT.txt" 2>/dev/null | awk '{print $2}')
echo "  status line: ${ST:-<none>}"
[ "${ST:-MISSING}" = "OK" ] || { echo "  [STOP] ${ST:-no status line} — not quotable."; OK=1; }
[ "$OK" -eq 0 ] && runlog "16b_environment complete" || runlog "**16b_environment FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
