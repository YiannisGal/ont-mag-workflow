#!/usr/bin/env bash
# =============================================================================
# 36_hmm_gap_closure.sh  --  run the two profile searches earlier steps deferred
#
# WHAT THIS REPLACES
#   The roadmap's step 36 was "dN/dS + pangenome openness". Step 32c already
#   produced dN/dS = 0.08 for the one strain pair we have, and step 32d showed
#   there is NO RESOLVABLE DIVERGENCE between those genomes at all. Running a
#   selection analysis on a pair with 22 differences in 385,644 sites would
#   produce a number with no information in it. Pangenome openness needs
#   multiple genomes per species and we have 28 species in 29 MAGs.
#   That step is SUPERSEDED and is recorded as such, not silently skipped.
#
#   What is genuinely outstanding is two HMM searches that earlier steps
#   identified, costed, and explicitly left undone:
#
#   GAP A — step 29 could not test 193 network edges. Its transporter calls came
#           from the step-11 annotation, and where a transporter family was
#           detected in <20% of genomes it invoked the step-19 rule ("that is
#           marker-set failure, not absence") and left the edge intact but
#           UNTESTED. The rule was right. The remedy was never run: search the
#           proteomes DIRECTLY with the KOfam profiles for those families, at
#           KEGG's own curated thresholds. Then the thin detection is either
#           confirmed as real absence or exposed as an annotation gap — and
#           either answer closes the edges.
#
#   GAP B — step 23d could not put a floor under T9SS cargo. It said so in its
#           own text: "WHAT WOULD SETTLE IT PROPERLY: hmmsearch with TIGR04183
#           and TIGR04131". The C-terminal sorting domain is what marks a
#           protein for Por secretion; without it, "secreted" was inference.
#
# THE TRANSPORTER LIST IS READ OUT OF STEP 29, NOT RETYPED HERE
#   Step 33 asserted an antiSMASH image tag and a container path and was wrong
#   about both; the fix was to read them out of the wrapper. Same principle: the
#   metabolite -> KO map is parsed from 29_transporter_check.sh itself, so this
#   step cannot drift from the step it is closing. If the parse fails it says so
#   and falls back to a literal list, flagged as a fallback.
#
# CONTROLS
#   POSITIVE, domain-declared (the 19b/23e error was a bacterial control applied
#     to archaea): K02358 tuf for Bacteria, K03231 EF-1alpha for Archaea. Both
#     are universal within their domain. If the search cannot find them, the
#     search is broken and NOTHING is reported.
#   THRESHOLD, internal: every hit that clears E<1e-10 but fails its family's
#     curated score threshold is counted. That number is the false-positive rate
#     a flat cutoff would have had — step 34 measured it at 1,878 hits.
#   GAP B, structural: a T9SS C-terminal domain should occur ONLY in genomes
#     that carry the T9SS core machinery (step 23d's t9ss_machinery.tsv). CTDs
#     scattered across genomes with no machinery are cross-hits, not cargo.
#
# OUTPUT  36_gapclose/transporter_matrix.tsv · edges_closed.tsv · ctd_cargo.tsv
#         VERDICT.txt
# RUNTIME 10-25 min (29 proteomes x ~17 profiles; gap B adds 2 more)
# USAGE   ws2 && bash scripts/36_hmm_gap_closure.sh 2>&1 | tee logs/36.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="36_gapclose"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"

A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"
T29="$WS2/29_transport"; T9="$WS2/23_extras/t9ss"
D="$WS2/36_gapclose"; mkdir -p "$D/profiles" "$D/hits" "$D/ctd"
SRC29="$(dirname "$0")/29_transporter_check.sh"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 36 — CLOSING THE TWO DEFERRED HMM SEARCHES ########"
echo "  started: $(date -Is)"
runlog "36_gapclose started"

# ------------------------------------------------------------------ binaries
# Resolve by finding the binary. ENV_<tool> is only set for tools the smoketest
# happened to probe; assuming it exists has now broken steps 33, 34 and 35.
bin_env() {
  local tool="$1"; local var="ENV_${tool//-/_}"
  local cands="${!var:-} ${ENV_hmmsearch:-} ${ENV_hmmscan:-} ${ENV_hmmer:-} ${ENV_bakta:-} -"
  local e prefix root d
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then
      command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }; continue
    fi
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
if [ -z "$HS" ]; then
  echo "  !! hmmsearch not found in PATH or under: ${CONDA_ROOTS:-<unset>}"
  runlog "**36 ABORTED: no hmmsearch**"; exit 2
fi
export ENV_hmmsearch="$HS"; echo "  hmmsearch: env '$HS'"
HF=$(bin_env hmmfetch || true); [ -n "$HF" ] && export ENV_hmmfetch="$HF"
HP=$(bin_env hmmpress  || true); [ -n "$HP" ] && export ENV_hmmpress="$HP"
record_version "hmmsearch" "$HS" hmmsearch -h 2>/dev/null || true

# ------------------------------------------------------------------ proteomes
mapfile -t MAGS < <(for d in "$A"/*/; do
  n=$(basename "$d"); [ -s "$d/$n.faa" ] && echo "$n"; done | sort)
echo "  proteomes: ${#MAGS[@]}"
if [ "${#MAGS[@]}" -lt 5 ]; then
  echo "  !! fewer than 5 proteomes under $A — nothing to search."
  runlog "**36 ABORTED: no proteomes**"; exit 3
fi

# ------------------------------------------------------------------- ko_list
for src in "$K/ko_list" "$WS2/34_unbinned/ko_list" "$WS2/19_auxotrophy/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$D/ko_list" 2>/dev/null || true; break; }
done
if [ ! -s "$D/ko_list" ]; then
  wget -q -O "$D/ko_list.gz" "$BASE/ko_list.gz" 2>/dev/null && gunzip -f "$D/ko_list.gz" 2>/dev/null || true
fi
if [ ! -s "$D/ko_list" ]; then
  echo "  !! no ko_list. Without curated thresholds every hit below is an"
  echo "     unverified assertion — step 34 showed a flat E<1e-10 produced 1,878"
  echo "     homologue hits that the curated cutoffs rejected."
  echo "     wget -O $D/ko_list.gz $BASE/ko_list.gz && gunzip $D/ko_list.gz"
  runlog "**36 ABORTED: no ko_list**"; exit 4
fi
echo "  ko_list: $(wc -l < "$D/ko_list") families"

# =========================================================================
echo
echo "######## A0. THE TRANSPORTER LIST, READ OUT OF STEP 29 ########"
# Parse the TRANS = { ... } literal from 29_transporter_check.sh. Do not retype
# it: a retyped list that drifts from step 29 would close edges step 29 never
# opened, and nobody would notice.
python3 - "$SRC29" "$D" <<'PYT'
import sys, os, re, ast
src, D = sys.argv[1], sys.argv[2]
rows, how = [], 'FALLBACK'
txt = ''
if os.path.exists(src):
    txt = open(src, errors='replace').read()
m = re.search(r'^TRANS\s*=\s*(\{.*?^\})', txt, re.S | re.M)
if m:
    try:
        trans = ast.literal_eval(m.group(1))
        for key, kos in trans.items():
            for ko, (label, kws) in kos.items():
                rows.append((ko, key, label))
        how = f'parsed from {os.path.basename(src)}'
    except Exception as e:
        print(f"  !! TRANS block found but would not parse: {e}")
if not rows:
    # literal fallback — flagged, never silent
    FB = {'acetate': {'K14393': 'acetate symporter',
                      'K02106': 'short-chain fatty acid transporter'},
          'formate': {'K06212': 'formate transporter'},
          'choline': {'K02168': 'choline/betaine transporter',
                      'K05020': 'glycine betaine transporter'},
          'b12':     {'K06858': 'B12 transport, binding',
                      'K06073': 'B12 transport, permease',
                      'K06074': 'B12 transport, ATPase'},
          'amino':   {'K01999': 'branched-chain AA transport',
                      'K01997': 'branched-chain AA transport',
                      'K01998': 'branched-chain AA transport',
                      'K10005': 'glutamate transport'},
          'biotin':  {'K03523': 'biotin transport'},
          'thiamin': {'K02064': 'thiamine transport',
                      'K02062': 'thiamine transport'}}
    for key, kos in FB.items():
        for ko, label in kos.items():
            rows.append((ko, key, label))
with open(os.path.join(D, 'transporter_kos.tsv'), 'w') as fh:
    fh.write('ko\tmetabolite_key\tlabel\n')
    for ko, key, label in sorted(set(rows)):
        fh.write(f'{ko}\t{key}\t{label}\n')
print(f"  source: {how}")
print(f"  {len(set(r[0] for r in rows))} transporter KOs across "
      f"{len(set(r[1] for r in rows))} metabolite families")
if how == 'FALLBACK':
    print("  !! step 29's own table could not be read. The list below is a literal")
    print("     copy and MAY HAVE DRIFTED from what step 29 actually tested.")
PYT

# ---- add the domain-declared positive controls --------------------------
{ cat "$D/transporter_kos.tsv"
  printf 'K02358\t__CONTROL_BACT\telongation factor Tu\n'
  printf 'K03231\t__CONTROL_ARCH\telongation factor 1-alpha\n'
} > "$D/kos_all.tsv"
mapfile -t KOS < <(awk -F'\t' 'NR>1{print $1}' "$D/kos_all.tsv" | sort -u)
echo "  searching ${#KOS[@]} profiles (transporters + 2 domain controls)"

# ---- LABEL GUARD: every KO must mean what this script says it means -----
# Step 34 labelled K16178 mtmB when KEGG calls it dimethylamine methyltransferase.
# The guard caught it before anything ran. Same guard here.
echo
echo "######## A1. KO LABEL GUARD ########"
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
KEY = {'acetate': ('acetate', 'fatty acid'), 'formate': ('formate',),
       'choline': ('choline', 'betaine'), 'b12': ('cobalamin', 'b12', 'btu'),
       'amino': ('amino acid', 'glutamate', 'branched', 'liv'),
       'biotin': ('biotin', 'bio'), 'thiamin': ('thiamin', 'thi'),
       '__CONTROL_BACT': ('elongation factor tu',),
       '__CONTROL_ARCH': ('elongation factor 1',)}
bad = miss = 0
for line in open(os.path.join(D, 'kos_all.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 3 or f[0] == 'ko':
        continue
    ko, key, label = f[0], f[1], f[2]
    d = defs.get(ko)
    if d is None:
        print(f"  [MISS] {ko} ({label}) is not in ko_list — no profile, no threshold")
        miss += 1
        continue
    kws = KEY.get(key, ())
    if kws and not any(k in d.lower() for k in kws):
        print(f"  [FAIL] {ko} is filed under metabolite family '{key}' as '{label}',")
        print(f"         but KEGG's own definition is: {d}")
        bad += 1
    else:
        print(f"  [ok]   {ko}  {d[:66]}")
if bad:
    print(f"  !! {bad} KO label(s) disagree with KEGG. NOTHING IS RUN.")
    sys.exit(1)
if miss:
    print(f"  !! {miss} KO(s) absent from ko_list; they are dropped, not guessed.")
sys.exit(0)
PYG
if [ $? -ne 0 ]; then runlog "**36 ABORTED: KO label guard failed**"; exit 5; fi

# ---- fetch profiles -----------------------------------------------------
echo
echo "######## A2. PROFILES ########"
# THE FIRST VERSION HID WHY FETCHES FAILED.
#   It ran `wget -q ... 2>/dev/null` and reported only that a profile was
#   absent. The run that followed downloaded 2 of 17 — both from a local cache —
#   and every network fetch failed silently. That is the same defect that hid
#   prodigal's failure in step 34: an error sent to /dev/null becomes an absence,
#   and an absence becomes a result. Every fetch now reports its HTTP status.
fetch_one() {  # fetch_one <KO> -> 0 on a real HMMER3 file
  local ko="$1" dst="$D/profiles/$1.hmm" code
  code=$(wget -S -O "$dst" "$BASE/profiles/${ko}.hmm" 2>&1 \
         | sed -n 's/^ *HTTP\/[0-9.]* \([0-9]*\).*/\1/p' | tail -1)
  if grep -q '^HMMER3' "$dst" 2>/dev/null; then return 0; fi
  rm -f "$dst"; echo "${code:-no-response}"; return 1
}
got=0; missing=(); codes=()
for ko in "${KOS[@]}"; do
  [ -s "$D/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$WS2/34_unbinned/profiles/${ko}.hmm" \
             "$WS2/34b_methylamine/profiles/${ko}.hmm" \
             "$WS2/19_auxotrophy/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$D/profiles/"; break; }
  done
  if [ ! -s "$D/profiles/${ko}.hmm" ]; then
    c=$(fetch_one "$ko") || { missing+=("$ko"); codes+=("$ko:$c"); continue; }
  fi
  got=$((got+1))
done
echo "  $got/${#KOS[@]} profiles available"
if [ "${#missing[@]}" -gt 0 ]; then
  echo "  !! not fetched: ${missing[*]}"
  echo "     HTTP status per KO: ${codes[*]}"
fi

# ---- fallback: the whole KOfam profile set, once ------------------------
# If per-KO fetching is not working (KEGG rate-limits, the path moved, the
# machine is behind a proxy), one tarball settles it. ~360 MB, kept in $D so a
# re-run never downloads it twice.
if [ "${#missing[@]}" -gt 0 ] && [ "${SKIP_TARBALL:-0}" != "1" ]; then
  TB="$D/profiles.tar.gz"
  if [ ! -s "$TB" ]; then
    AV=$(df -BG --output=avail "$D" 2>/dev/null | tail -1 | tr -dc '0-9'); AV=${AV:-0}
    if [ "$AV" -lt 5 ]; then
      echo "  !! only ${AV} GB free at $D — not downloading the 360 MB profile set."
    else
      echo "  fetching the full KOfam profile set (~360 MB, once) ..."
      wget -O "$TB" "$BASE/profiles.tar.gz" 2>&1 | tail -2 | sed 's/^/     /'
    fi
  fi
  if [ -s "$TB" ]; then
    echo "  extracting only the ${#missing[@]} missing profile(s) ..."
    for ko in "${missing[@]}"; do
      tar xzf "$TB" -C "$D/profiles" --strip-components=1 "profiles/${ko}.hmm" 2>/dev/null \
        || tar xzf "$TB" -C "$D/profiles" "${ko}.hmm" 2>/dev/null || true
    done
    still=(); got=0
    for ko in "${KOS[@]}"; do
      if [ -s "$D/profiles/${ko}.hmm" ]; then got=$((got+1)); else still+=("$ko"); fi
    done
    missing=("${still[@]}")
    echo "  after the tarball: $got/${#KOS[@]} profiles"
  fi
fi
if [ "${#missing[@]}" -gt 0 ]; then
  echo "  !! NO PROFILE for: ${missing[*]}"
  echo "     These are carried through the whole report as NOT SEARCHED. They are"
  echo "     never counted as absent, and no edge is decided from them."
fi
printf '%s\n' "${missing[@]}" | grep -v '^$' > "$D/not_searched.txt" || true

cat "$D"/profiles/*.hmm > "$D/all_profiles.hmm" 2>/dev/null
if [ ! -s "$D/all_profiles.hmm" ]; then
  echo "  !! no profiles at all."; runlog "**36 ABORTED: no profiles**"; exit 6
fi

# ---- refuse to run a closure that cannot close anything -----------------
# The first run searched 2 of 17 profiles and then reported 190 edges
# "DOWNGRADED — direct search finds nothing". Nothing had been searched. Only
# the positive control stopped that from being written down as a result.
NTR=0; NTRGOT=0
while IFS=$'\t' read -r ko key label; do
  [ "$ko" = "ko" ] && continue
  NTR=$((NTR+1)); [ -s "$D/profiles/${ko}.hmm" ] && NTRGOT=$((NTRGOT+1))
done < "$D/transporter_kos.tsv"
echo "  transporter profiles: $NTRGOT/$NTR"
if [ "$NTR" -gt 0 ] && [ "$NTRGOT" -lt $(( (NTR + 1) / 2 )) ]; then
  echo
  echo "  !! fewer than half the transporter profiles are available."
  echo "     Gap A cannot be closed and a run would manufacture absences."
  echo "     Get the profiles first:"
  echo "       wget -O $D/profiles.tar.gz $BASE/profiles.tar.gz"
  echo "       (then re-run; extraction and everything else is cached)"
  runlog "**36 ABORTED: transporter profile coverage $NTRGOT/$NTR**"
  exit 7
fi

# a control KO with no profile means the control cannot run — say so loudly
for c in K02358 K03231; do
  [ -s "$D/profiles/${c}.hmm" ] || echo "  !! CONTROL profile $c is missing — its domain will have NO control."
done

# ---- cache stamp --------------------------------------------------------
# The hits cache MUST be keyed on the profile set. The second run of v1 reused
# tables built from 2 profiles and reported them as a 17-profile search.
PSTAMP=$(md5sum "$D/all_profiles.hmm" | cut -c1-12)
echo "  profile-set stamp: $PSTAMP"

# ---- search -------------------------------------------------------------
echo
echo "######## A3. SEARCHING ${#MAGS[@]} PROTEOMES ########"
n=0
for name in "${MAGS[@]}"; do
  n=$((n+1)); tb="$D/hits/${name}.tbl"; st="$D/hits/${name}.stamp"
  if [ -s "$tb" ] && [ "$(cat "$st" 2>/dev/null)" = "$PSTAMP" ]; then
    printf '\r  [%2d/%2d] %-30s cached ' "$n" "${#MAGS[@]}" "$name"; continue
  fi
  printf '\r  [%2d/%2d] %-30s        ' "$n" "${#MAGS[@]}" "$name"
  run_tool hmmsearch --cpu "$THREADS" -E 1e-10 --domE 1e-10 \
      --tblout "$tb.part" "$D/all_profiles.hmm" "$A/$name/$name.faa" \
      >/dev/null 2>"$D/hits/${name}.err"
  rc=$?
  # exit 0 with no table is not success — step 33 accepted 29 of those
  if [ "$rc" -ne 0 ] || [ ! -s "$tb.part" ]; then
    printf '\r  [%2d/%2d] %-30s FAILED rc=%s\n' "$n" "${#MAGS[@]}" "$name" "$rc"
    head -3 "$D/hits/${name}.err" 2>/dev/null | sed 's/^/       /'
    rm -f "$tb.part"; continue
  fi
  mv "$tb.part" "$tb"; echo "$PSTAMP" > "$st"
done
echo
echo "  tables written: $(ls -1 "$D/hits"/*.tbl 2>/dev/null | wc -l)/${#MAGS[@]}"

# =========================================================================
echo
echo "######## B. T9SS C-TERMINAL SORTING DOMAIN ########"
# TIGR04183 (type A CTD) and TIGR04131 (type B). These are TIGRFAM/NCBIfam
# models and are NOT in Pfam-A, which is why step 23d could not run this.
CTD_OK=0
CTDH="$D/ctd/ctd.hmm"
if [ -s "$CTDH" ]; then
  CTD_OK=1; echo "  models: cached"
else
  # look for any local HMM library that contains them, then fetch by accession
  LIBS=$(ls -1 ${DB_DIR}/*/*.{hmm,HMM} \
                ${DB_DIR}/*/*/*.hmm 2>/dev/null | head -40)
  for lib in $LIBS; do
    grep -qE '^ACC +TIGR04183' "$lib" 2>/dev/null || continue
    if [ -n "${ENV_hmmfetch:-}" ]; then
      run_tool hmmfetch -o "$D/ctd/a.hmm" "$lib" TIGR04183 >/dev/null 2>&1
      run_tool hmmfetch -o "$D/ctd/b.hmm" "$lib" TIGR04131 >/dev/null 2>&1
      cat "$D/ctd"/[ab].hmm > "$CTDH" 2>/dev/null
      [ -s "$CTDH" ] && { CTD_OK=1; echo "  models: fetched from $lib"; break; }
    fi
  done
fi
if [ "$CTD_OK" -eq 0 ]; then
  echo "  !! TIGR04183 / TIGR04131 not found locally. GAP B IS SKIPPED — it is"
  echo "     not reported as 'no cargo found', which would be a fabricated"
  echo "     negative. Gap A below is unaffected and stands on its own."
  echo
  echo "     To get them (NCBIfam carries both):"
  echo "       mkdir -p $D/ctd && cd $D/ctd"
  echo "       wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz"
  echo "       tar xzf hmm_PGAP.HMM.tgz"
  echo "       cat hmm_PGAP/TIGR04183.HMM hmm_PGAP/TIGR04131.HMM > ctd.hmm"
  echo "     then re-run this script; everything else is cached."
else
  m=0
  for name in "${MAGS[@]}"; do
    m=$((m+1)); ct="$D/ctd/${name}.tbl"
    [ -s "$ct" ] && continue
    printf '\r  [%2d/%2d] %-30s' "$m" "${#MAGS[@]}" "$name"
    run_tool hmmsearch --cpu "$THREADS" --cut_nc --tblout "$ct.part" \
        "$CTDH" "$A/$name/$name.faa" >/dev/null 2>&1
    # --cut_nc needs the model to carry NC lines; fall back to an E-value
    [ -s "$ct.part" ] || run_tool hmmsearch --cpu "$THREADS" -E 1e-5 \
        --tblout "$ct.part" "$CTDH" "$A/$name/$name.faa" >/dev/null 2>&1
    [ -s "$ct.part" ] && mv "$ct.part" "$ct" || rm -f "$ct.part"
  done
  echo
  echo "  CTD tables: $(ls -1 "$D/ctd"/*.tbl 2>/dev/null | wc -l)/${#MAGS[@]}"
fi

# =========================================================================
echo
python3 - "$D" "$F" "$T29" "$T9" "$WS2" "$CTD_OK" <<'PY'
import sys, os, re, csv
from collections import defaultdict, Counter

# 🔴 THE TOKEN FILTER, TIGHTENED -- see Rule 46.
# A filter accepting any string of >=4 characters made 'Bacteria',
# 'Methanosarcina' and '99.85' into alias keys; a global alias merge then linked
# unrelated genomes through them, and reading the resulting set gave a DIFFERENT
# answer per process because set order follows hash randomisation. Two identical
# runs of 16b returned different genera for the same genome.
# A genome id in this project always contains an underscore AND a digit
# (ETO_LIQ_MAG_002, LIQ_das_metabat_075); a genus never does in a way that
# collides (Lentimicrobium, GCA-2746605, Methanobacterium_C, JBAIYS01).
IDTOK = re.compile(r'^(?=.*_)(?=.*\d)[A-Za-z0-9._-]{6,}$')

D, F, T29, T9, WS2, CTD_OK = (sys.argv[1], sys.argv[2], sys.argv[3],
                              sys.argv[4], sys.argv[5], sys.argv[6] == '1')
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)

# ------------------------------------------------------------ thresholds
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

komap = {}   # ko -> (metabolite_key, label)
for line in open(os.path.join(D, 'kos_all.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 3 and f[0] != 'ko':
        komap[f[0]] = (f[1], f[2])

# A KO with no profile was NOT SEARCHED. It is not absent. The first version
# treated the two the same and reported 190 edges "DOWNGRADED — direct search
# finds nothing" after searching 2 profiles of 17. Only the positive control
# stopped that from being written down.
NOTSEARCHED = set()
nsp = os.path.join(D, 'not_searched.txt')
if os.path.exists(nsp):
    NOTSEARCHED = {l.strip() for l in open(nsp) if l.strip()}
else:
    NOTSEARCHED = {k for k in komap
                   if not os.path.exists(os.path.join(D, 'profiles', f'{k}.hmm'))}
SEARCHED = {k for k in komap if k not in NOTSEARCHED}

# ------------------------------------------------------------ domains
def parse_taxonomy(paths):
    """Key the domain/genus map on EVERY id-like field in each row.

    The proteome directories are named by BIN id (SED_das_vamb_004) while the
    catalogue may key on the ETO working name. v1 picked one column, matched
    nothing, and declared 29 proteomes 'undeclared' — which turned a metadata
    join failure into a failed positive control. Index every plausible id and
    the join cannot miss."""
    dom, gen = {}, {}
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path, errors='replace') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not rows:
            continue
        cols = [c for c in rows[0].keys() if c]
        dc = next((c for c in cols if c.strip().lower() == 'domain'), None)
        gc = next((c for c in cols if c.strip().lower() == 'genus'), None)
        for r in rows:
            blob = ' '.join(str(v or '') for v in r.values())
            d = (r.get(dc) or '').strip() if dc else ''
            if not d:
                m = re.search(r'd__([A-Za-z]+)', blob)
                d = m.group(1) if m else ''
            g = (r.get(gc) or '').strip() if gc else ''
            if not g:
                m = re.search(r'g__([^;\t]+)', blob)
                g = m.group(1).strip() if m else ''
            ids = set()
            for c in cols:
                v = (r.get(c) or '').strip()
                # an id: no spaces, not a lineage string, plausibly a name
                if v and ' ' not in v and ';' not in v and IDTOK.match(v) \
                        and not re.fullmatch(r'[0-9._-]+', v) \
                        and not v.startswith(('d__', 'p__', 'c__', 'o__',
                                              'f__', 'g__', 's__')):
                    ids.add(v)
            for i in ids:
                if d and i not in dom:
                    dom[i] = d
                if g and i not in gen:
                    gen[i] = g
    return dom, gen

def id_fields(row, cols):
    out = set()
    for c in cols:
        v = (row.get(c) or '').strip()
        if v and ' ' not in v and ';' not in v and IDTOK.match(v) \
                and not re.fullmatch(r'[0-9._-]+', v) \
                and not v.startswith(('d__', 'p__', 'c__', 'o__',
                                      'f__', 'g__', 's__')):
            out.add(v)
    return out

def propagate_aliases(dom, gen, path):
    """MAG_CROSSWALK.tsv links the ETO working name to the bin id but carries no
    taxonomy of its own. Reading it as a taxonomy table stores nothing, so the
    join still fails — which is what a first fix did. Treat each of its rows as
    an ALIAS GROUP and push the domain across the group instead."""
    if not os.path.exists(path):
        return 0
    with open(path, errors='replace') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if not rows:
        return 0
    cols = [c for c in rows[0].keys() if c]
    n = 0
    for r in rows:
        ids = id_fields(r, cols)
        d = next((dom[i] for i in sorted(ids) if i in dom), None)
        g = next((gen[i] for i in sorted(ids) if i in gen), None)
        for i in ids:
            if d and i not in dom:
                dom[i] = d; n += 1
            if g and i not in gen:
                gen[i] = g
    return n

dom, gen = parse_taxonomy([os.path.join(F, 'MAG_TAXONOMY.tsv')])
_al = propagate_aliases(dom, gen, os.path.join(F, 'MAG_CROSSWALK.tsv'))

# ------------------------------------------------------------ read hits
# hmmsearch --tblout: 0 target(protein) 2 query(profile NAME = KO)
#                     4 E(full) 5 score(full) 7 E(dom) 8 score(dom)
hits = defaultdict(set)       # ko -> {mag}
nprot = defaultdict(int)      # ko -> protein count (passing)
sub = Counter()               # ko -> hits that cleared E but failed threshold
nothr = set()
mags = []
for fn in sorted(os.listdir(os.path.join(D, 'hits'))):
    if not fn.endswith('.tbl'):
        continue
    mag = fn[:-4]
    mags.append(mag)
    for line in open(os.path.join(D, 'hits', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 9:
            continue
        ko = f[2] if f[2] in komap else (f[3] if len(f) > 3 and f[3] in komap else None)
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

rule('A1. CONTROL — CAN THIS SEARCH FIND WHAT MUST BE THERE?')
bact = [m for m in mags if dom.get(m, '').lower().startswith('bacteri')]
arch = [m for m in mags if dom.get(m, '').lower().startswith('archae')]
unk = [m for m in mags if m not in bact and m not in arch]
say(f"  {len(bact)} bacterial, {len(arch)} archaeal, {len(unk)} undeclared proteomes")
if _al:
    say(f"  ({_al} name(s) resolved through MAG_CROSSWALK.tsv — the proteome"
        f" directories and the catalogue do not use the same ids)")
say()
say("  A control is DOMAIN-DECLARED. Archaea have no EF-Tu — they use EF-1alpha —")
say("  and steps 19b and 23e both once applied a bacterial control to archaeal")
say("  genomes and read the guaranteed absence as a result.")
say()
CTRL = {}
for ko, grp, gname in (('K02358', bact, 'Bacteria'), ('K03231', arch, 'Archaea')):
    have = len([m for m in grp if m in hits.get(ko, ())])
    prof = os.path.exists(os.path.join(D, 'profiles', f'{ko}.hmm'))
    if not grp:
        say(f"  {ko} ({gname:<8}) : no genomes of this domain — control not applicable")
        CTRL[gname] = None
    elif not prof:
        say(f"  {ko} ({gname:<8}) : PROFILE MISSING — this domain has NO control")
        CTRL[gname] = False
    else:
        r = 100.0 * have / len(grp)
        ok = r >= 80.0
        CTRL[gname] = ok
        say(f"  {ko} ({gname:<8}) : {have}/{len(grp)} genomes  ({r:.0f}%)  "
            f"[{'PASS' if ok else 'FAIL'}]")
say()
live = [v for v in CTRL.values() if v is not None]
CTRL_OK = bool(live) and all(live)
DEGRADED = False
if not live:
    # No domain could be joined to any proteome. That is a metadata failure, not
    # a search failure, and blocking on it would discard a working search. Fall
    # back to a weaker but still valid control: EVERY genome must carry one
    # elongation factor or the other. It does not test domain assignment, and
    # the report says so wherever it is used.
    DEGRADED = True
    have = len([m for m in mags
                if m in hits.get('K02358', ()) or m in hits.get('K03231', ())])
    r = 100.0 * have / max(1, len(mags))
    CTRL_OK = r >= 80.0
    say("  !! no domain could be joined to any proteome — check that")
    say(f"     {os.path.join(F, 'MAG_TAXONOMY.tsv')} keys on the same names as")
    say(f"     {os.path.join(WS2, '09_annotation')}/*/ (bin id vs ETO working name).")
    say()
    say("  DEGRADED CONTROL, union form: every genome must carry tuf OR EF-1alpha.")
    say(f"  {have}/{len(mags)} genomes do ({r:.0f}%)  [{'PASS' if CTRL_OK else 'FAIL'}]")
    say("  This does not check the domain call — only that the search works.")
elif CTRL_OK:
    say("  [PASS] the profile search recovers universal single-copy genes in every")
    say("         domain present. Absences below are therefore informative.")
if not CTRL_OK:
    say("  [FAIL] the search cannot find genes that must be in these genomes. Every")
    say("         'absent' below is uninterpretable. NO VERDICT.")

rule('A2. THRESHOLDS — WHAT A FLAT E-VALUE WOULD HAVE ADDED')
say(f"  computed over the {len(SEARCHED)} families that were actually searched"
    + (f"; {len(NOTSEARCHED)} had no profile and contribute nothing"
       if NOTSEARCHED else ''))
say()
tot_sub = sum(v for k, v in sub.items() if k in SEARCHED)
tot_keep = sum(v for k, v in nprot.items() if k in SEARCHED)
say(f"  passed the family's curated KOfam threshold : {tot_keep}")
say(f"  cleared E<1e-10 but FAILED that threshold   : {tot_sub}")
if tot_keep + tot_sub:
    say(f"  a flat E<1e-10 would have inflated this table by "
        f"{100.0 * tot_sub / max(1, tot_keep):.0f}%")
say()
if sub:
    say(f"  {'KO':<9}{'kept':>7}{'rejected':>10}   family")
    say('-' * W)
    for ko, c in sorted(sub.items(), key=lambda kv: (-kv[1], kv[0]))[:12]:
        say(f"  {ko:<9}{nprot.get(ko, 0):>7}{c:>10}   {komap.get(ko, ('', ''))[1][:52]}")
if nothr:
    say()
    say(f"  no curated threshold for {', '.join(sorted(nothr))} — those are kept on")
    say("  the E-value alone and are marked as such in the matrix.")

rule('A3. TRANSPORTER PRESENCE — DIRECT SEARCH vs STEP 29')
say(f"  {'KO':<9}{'metabolite':<12}{'genomes':>9}{'%':>7}   family")
say('-' * W)
matrix = os.path.join(D, 'transporter_matrix.tsv')
with open(matrix, 'w') as mh:
    mh.write('ko\tmetabolite_key\tlabel\tn_genomes\tpct_genomes\tproteins\t'
             'sub_threshold\tgenomes\n')
    for ko in sorted(komap):
        key, label = komap[ko]
        if key.startswith('__CONTROL'):
            continue
        if ko in NOTSEARCHED:
            say(f"  {ko:<9}{key:<12}{'NOT SEARCHED':>16}   {label[:46]}")
            mh.write(f"{ko}\t{key}\t{label}\tNOT_SEARCHED\t-\t-\t-\t\n")
            continue
        g = sorted(hits.get(ko, ()))
        pct = 100.0 * len(g) / max(1, len(mags))
        say(f"  {ko:<9}{key:<12}{len(g):>9}{pct:>6.0f}%   {label[:50]}")
        mh.write(f"{ko}\t{key}\t{label}\t{len(g)}\t{pct:.1f}\t{nprot.get(ko, 0)}\t"
                 f"{sub.get(ko, 0)}\t{','.join(g)}\n")
say()
say("  Step 29 called a family DETECTED-THIN when the step-11 annotation found it")
say("  in <20% of genomes, and by the step-19 rule refused to read that as")
say("  absence. Any family at or above 20% here that step 29 had below it was an")
say("  ANNOTATION gap; any family still below 20% is a real absence, now measured")
say("  directly rather than inferred from a pipeline that might not have looked.")

# ------------------------------------------------- close the open edges
rule('A4. EDGES CLOSED')
ep = os.path.join(T29, 'edge_transport.tsv')
by_key = defaultdict(list)
for ko, (key, _l) in komap.items():
    if not key.startswith('__CONTROL'):
        by_key[key].append(ko)


# ---- NAME JOIN: edge tables key on ETO working names, proteome dirs on bin ids
# Step 36's edge closure looked up `consumer` (ETO_LIQ_MAG_003) in a table keyed
# on proteome directory names (LIQ_das_metabat_001). The intersection is EMPTY,
# so `found` was always empty and every openable edge became DOWNGRADED without
# the transporter ever being looked up. The domain control needed the crosswalk
# and got it; the edge closure needed it and did not.
def build_alias(paths):
    """canonical name -> set of all names for the same genome, both directions"""
    groups = []
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path, errors='replace') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not rows:
            continue
        cols = [c for c in rows[0].keys() if c]
        for r in rows:
            ids = set()
            for c in cols:
                v = (r.get(c) or '').strip()
                if v and ' ' not in v and ';' not in v and IDTOK.match(v) \
                        and not re.fullmatch(r'[0-9._-]+', v) \
                        and not v.startswith(('d__', 'p__', 'c__', 'o__',
                                              'f__', 'g__', 's__')):
                    ids.add(v)
            if len(ids) > 1:
                groups.append(ids)
    alias = {}
    for g in groups:
        for i in g:
            alias[i] = set(g)          # row-local; never merged
    return alias

ALIAS = build_alias([os.path.join(F, 'MAG_CROSSWALK.tsv'),
                     os.path.join(F, 'MAG_TAXONOMY.tsv')])

def resolve(nm, known):
    """map a name onto whichever alias of it appears in `known`"""
    if nm in known:
        return nm
    for a in sorted(ALIAS.get(nm, ())):
        if a in known:
            return a
    return None

KNOWN = set(mags)
_res = sum(1 for r in {(r.get('consumer') or '').strip()
                       for r in csv.DictReader(open(ep, errors='replace'), delimiter='\t')}
           if r and resolve(r, KNOWN)) if os.path.exists(ep) else 0

closed = reopened = still = nofam = unsearched = unjoined = 0
# a metabolite family is decidable only if at least one of its KOs was searched
decidable = {k for k in by_key if any(ko in SEARCHED for ko in by_key[k])}
rows_out = []
if not os.path.exists(ep):
    say(f"  !! {ep} not found — step 29's edge table is what this closes.")
    say("     Re-run step 29, or this section has nothing to act on.")
else:
    with open(ep, errors='replace') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    say(f"  step 29 edges: {len(rows)}")
    vc = Counter((r.get('verdict') or '').strip() for r in rows)
    for v, c in vc.most_common():
        say(f"    {v:<26}{c:>6}")
    say()
    for r in rows:
        v = (r.get('verdict') or '').strip()
        cons = (r.get('consumer') or '').strip()
        met = (r.get('metabolite') or '').strip()
        if v not in ('untested', 'DOWNGRADED'):
            continue
        cons = resolve(cons, KNOWN) or cons
        if cons not in KNOWN:
            # the consumer has no proteome under any of its names, so the
            # transporter cannot be looked up for it AT ALL. That is not a
            # downgrade, it is an unanswerable question.
            unjoined += 1
            rows_out.append((cons, met, v, 'NOT JOINED (no proteome for this consumer)', ''))
            continue
        key = next((k for k in by_key if k in met.lower()), None)
        if key is None:
            nofam += 1
            continue
        if key not in decidable:
            # NO PROFILE WAS SEARCHED FOR THIS FAMILY. Calling the edge
            # DOWNGRADED here would report a search that never happened.
            unsearched += 1
            rows_out.append((cons, met, v, 'NOT SEARCHED (no profile for this family)', ''))
            continue
        found = [ko for ko in by_key[key] if ko in SEARCHED and cons in hits.get(ko, ())]
        if found:
            new = 'supported'
            closed += 1
        else:
            new = 'DOWNGRADED (direct search, transporter absent)'
            if v == 'untested':
                reopened += 1
            else:
                still += 1
        rows_out.append((cons, met, v, new, ','.join(found)))
    with open(os.path.join(D, 'edges_closed.tsv'), 'w') as oh:
        oh.write('consumer\tmetabolite\told_verdict\tnew_verdict\ttransporter_kos\n')
        for t in rows_out:
            oh.write('\t'.join(t) + '\n')
    say(f"  {closed:>5}  now SUPPORTED — the transporter is there; step 11 missed it")
    say(f"  {reopened:>5}  now DOWNGRADED — was untested, direct search finds nothing")
    say(f"  {still:>5}  DOWNGRADED confirmed by direct search")
    say(f"  {unsearched:>5}  NOT SEARCHED — no profile for that family; still open")
    say(f"  {unjoined:>5}  NOT JOINED — no proteome matches that consumer's name")
    say(f"  {nofam:>5}  still untestable — no transporter family for that metabolite")
    say()
    if unjoined:
        say()
        say(f"  !! {unjoined} edges name a consumer with no proteome under any alias.")
        say("     Those were NOT decided — the lookup could not run. If this number")
        say("     is large the crosswalk is incomplete and the closure is not usable.")
    tot = closed + reopened + still
    den = tot + nofam + unsearched + unjoined
    if tot:
        say(f"  {tot} of step 29's open edges are now decided "
            f"({100.0 * tot / max(1, den):.0f}% of the open ones).")
    if unsearched:
        say()
        say(f"  >> {unsearched} edges could not be touched because their transporter")
        say("     family had no profile. They are NOT downgraded. Fetch the profiles")
        say("     and re-run; the proteome search re-keys on the profile set, so it")
        say("     will re-search rather than reuse tables built from a smaller set.")
    say()
    say("  A DOWNGRADED edge is not a deleted edge. It means step 21 predicted an")
    say("  exchange whose uptake route we cannot find, in a genome of stated")
    say("  completeness. Report the number, and report the completeness with it.")

# ------------------------------------------------------------ gap B
rule('B. T9SS CARGO — DOES THE C-TERMINAL SORTING DOMAIN CONFINE ITSELF?')
CTD_VERDICT = 'NOTSEARCHED'
if not CTD_OK:
    say("  SKIPPED — TIGR04183/TIGR04131 not available locally. This is recorded as")
    say("  NOT SEARCHED. Step 23d's cargo floor remains open, exactly as 23d left")
    say("  it. It is not reported as zero cargo.")
else:
    ctd = defaultdict(set)
    for fn in sorted(os.listdir(os.path.join(D, 'ctd'))):
        if not fn.endswith('.tbl'):
            continue
        mag = fn[:-4]
        for line in open(os.path.join(D, 'ctd', fn), errors='replace'):
            if line.startswith('#'):
                continue
            f = line.split()
            if len(f) >= 3:
                ctd[mag].add(f[0])
    # ---- machinery call -------------------------------------------------
    # v3 read the id column straight out of t9ss_machinery.tsv and tested
    # `m in mach` against proteome directory names. t9ss_machinery.tsv keys on
    # ETO working names (ETO_LIQ_MAG_003); `mags` are bin ids
    # (LIQ_das_metabat_001). The intersection is EMPTY, so the header line said
    # "3 genomes with T9SS" and the table beneath it said 0 — the same crosswalk
    # bug already fixed twice in this script, in the domain control and in the
    # edge closure, and never applied here. A contradiction between a count and
    # the table built from it is a join failure until proven otherwise.
    mach_raw, mach = set(), set()
    mp = os.path.join(T9, 't9ss_machinery.tsv')
    if os.path.exists(mp):
        with open(mp, errors='replace') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if rows:
            cols = list(rows[0].keys())
            say(f"  t9ss_machinery.tsv columns: {', '.join(c for c in cols if c)}")
            idc = next((c for c in cols if c and re.search(r'mag|genome|bin', c, re.I)),
                       cols[0])
            cnc = next((c for c in cols if c and re.fullmatch(r'n_core|core', c.strip(), re.I)),
                       None)
            if cnc is None:
                cnc = next((c for c in cols if c and re.search(r'core|count|n_', c, re.I)
                            and not re.search(r'accessor|gene', c, re.I)), None)
            for r in rows:
                mid = (r.get(idc) or '').strip()
                if not mid:
                    continue
                if cnc:
                    try:
                        if float(re.sub(r'[^0-9.]', '', r.get(cnc) or '0') or 0) > 0:
                            mach_raw.add(mid)
                    except ValueError:
                        mach_raw.add(mid)
                else:
                    mach_raw.add(mid)
            say(f"  machinery call taken from column '{cnc or idc}' "
                f"-> {len(mach_raw)} rows called T9SS-positive")
            # ---- CONTROL B0: does the id namespace join onto the proteomes?
            allids = {(r.get(idc) or '').strip() for r in rows}
            allids.discard('')
            direct = len(allids & set(mags))
            joined = {i: resolve(i, KNOWN) for i in allids}
            viacw = sum(1 for i, v in joined.items() if v)
            say(f"  id join: {direct}/{len(allids)} table ids match a proteome name "
                f"directly, {viacw}/{len(allids)} after crosswalk resolution")
            for i in sorted(mach_raw):
                v = resolve(i, KNOWN)
                if v:
                    mach.add(v)
            if mach_raw and not mach:
                say(f"  !! none of the {len(mach_raw)} T9SS-positive ids resolve onto a")
                say(f"     proteome. example table id: {sorted(mach_raw)[0]}")
                say(f"     example proteome name    : {sorted(mags)[0] if mags else '-'}")
                say("     That is a JOIN FAILURE, not an absence of machinery.")
            elif len(mach) != len(mach_raw):
                say(f"  note: {len(mach_raw) - len(mach)} T9SS-positive rows have no"
                    f" proteome in this set (kept out of both classes)")
    else:
        say(f"  !! {mp} not found — the confinement control cannot run.")
    say()
    JOINFAIL = bool(mach_raw) and not mach
    withm = [m for m in mags if m in mach]
    without = [m for m in mags if m not in mach]
    cw = sum(len(ctd.get(m, ())) for m in withm)
    co = sum(len(ctd.get(m, ())) for m in without)
    say(f"  {'genome class':<34}{'genomes':>9}{'CTD proteins':>14}{'per genome':>12}")
    say('-' * W)
    say(f"  {'T9SS machinery present':<34}{len(withm):>9}{cw:>14}"
        f"{cw / max(1, len(withm)):>12.1f}")
    say(f"  {'no T9SS machinery':<34}{len(without):>9}{co:>14}"
        f"{co / max(1, len(without)):>12.1f}")
    say()
    if JOINFAIL:
        CTD_VERDICT = 'JOINFAIL'
        say("  [NO CONTROL — JOIN FAILED] the machinery table and the proteome set")
        say("         do not share a name space, so every genome landed in the")
        say("         'no machinery' class by construction. The confinement test")
        say("         did not run. These counts are neither validated nor refuted;")
        say("         they are unjoined. Do not read the FAIL that this table would")
        say("         otherwise print as biology.")
    elif not mach:
        CTD_VERDICT = 'NOCTRL'
        say("  [NO CONTROL] no machinery call available; counts stand unvalidated.")
    elif cw and (co / max(1, len(without))) * 5 <= (cw / max(1, len(withm))):
        CTD_VERDICT = 'PASS'
        say("  [PASS] the domain confines itself to genomes that carry the secretion")
        say("         system it feeds. These counts are a real cargo floor.")
    else:
        CTD_VERDICT = 'FAIL'
        say("  [FAIL] CTDs appear about as often in genomes with no T9SS as in")
        say("         genomes with one. That is a cross-hitting model, not cargo.")
        say("         Do not quote these counts.")
    say()
    # ---- phylum sanity check, independent of the machinery table ---------
    # T9SS is a Bacteroidota system. If the CTD carriers are Bacteroidota the
    # model is behaving even when the machinery join is broken; if they are
    # scattered across unrelated phyla the model is a cross-hitter regardless.
    phy = {}
    tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
    if os.path.exists(tp):
        with open(tp, errors='replace') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                blob = ' '.join(str(v or '') for v in r.values())
                m = re.search(r'p__([^;\t]+)', blob)
                p = m.group(1).strip() if m else ''
                if not p:
                    pc = next((c for c in r if c and c.strip().lower() == 'phylum'), None)
                    p = (r.get(pc) or '').strip() if pc else ''
                if not p:
                    continue
                for i in id_fields(r, [c for c in r.keys() if c]):
                    phy.setdefault(i, p)
        for grp in ({v for v in ALIAS.get(k, ())} for k in list(ALIAS)):
            p = next((phy[i] for i in grp if i in phy), None)
            if p:
                for i in grp:
                    phy.setdefault(i, p)
    ctdtot = sum(len(v) for v in ctd.values())
    if phy and ctdtot:
        byphy = Counter()
        for m, ps in ctd.items():
            byphy[phy.get(m, 'unassigned')] += len(ps)
        say(f"  {'phylum of the CTD carriers':<40}{'CTD proteins':>14}{'% of all':>12}")
        say('-' * W)
        for p, n in sorted(byphy.items(), key=lambda kv: (-kv[1], kv[0]))[:8]:
            say(f"  {p[:39]:<40}{n:>14}{100.0 * n / ctdtot:>11.1f}%")
        bact_ct = byphy.get('Bacteroidota', 0)
        say()
        say(f"  Bacteroidota share of all CTD hits: "
            f"{100.0 * bact_ct / ctdtot:.1f}%")
        if bact_ct > 0.6 * ctdtot:
            say("  >> T9SS is a Bacteroidota system. The domain is landing where the")
            say("     secretion system lives, which is independent evidence that the")
            say("     model is not cross-hitting — it does not depend on the")
            say("     machinery table joining.")
        else:
            say("  >> The domain is NOT concentrated in Bacteroidota. Treat these")
            say("     counts as unvalidated whatever the machinery join says.")
        say()
    say(f"  {'genome':<30}{'genus':<24}{'CTD proteins':>14}{'T9SS':>8}")
    say('-' * W)
    for m in sorted(mags, key=lambda x: (-len(ctd.get(x, ())), x))[:15]:
        if not ctd.get(m):
            continue
        say(f"  {m:<30}{(gen.get(m, '') or '-')[:23]:<24}"
            f"{len(ctd[m]):>14}{('yes' if m in mach else 'no'):>8}")
    # cross with the CAZyme call from 23d, if it is there
    cz = os.path.join(T9, 'cazyme_endo_exo.tsv')
    if os.path.exists(cz):
        czp = set()
        with open(cz, errors='replace') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                for v in r.values():
                    v = str(v or '')
                    if re.match(r'^[A-Za-z0-9_.|-]{6,}$', v) and '_' in v:
                        czp.add(v)
        allctd = set().union(*ctd.values()) if ctd else set()
        both = czp & allctd
        say()
        # CONTROL B1: the two protein-id name spaces must be comparable at all.
        # v3 printed "0" for this intersection and I read it as a finding. Zero
        # overlap between two id sets that never shared a prefix is a join
        # failure with the same shape as the machinery bug above.
        def _pref(s):
            return Counter(re.split(r'[_.|]', x)[0] for x in s if x)
        pz, pc = _pref(czp), _pref(allctd)
        shared_pref = set(pz) & set(pc)
        say(f"  23d CAZyme protein ids: {len(czp):,}   CTD protein ids: {len(allctd):,}")
        if not czp or not allctd:
            say("  [NO CROSS] one of the two id sets is empty; nothing to intersect.")
        elif not shared_pref:
            say("  [JOIN FAILED] the two id sets share no id prefix at all "
                f"(e.g. {sorted(pz)[0]!r} vs {sorted(pc)[0]!r}).")
            say("  They are different name spaces, so the intersection is 0 by")
            say("  construction. This is NOT 'no CAZyme carries a CTD' — it is an")
            say("  unrun comparison. Report it as not tested.")
        else:
            say(f"  CAZyme proteins listed by 23d that carry a CTD: {len(both)}")
            say("  Those are the ones with an actual export tag. A CAZyme without a")
            say("  CTD may still be secreted by another route — this is a FLOOR,")
            say("  which is the quantity 23d said it could not put a number on.")
    with open(os.path.join(D, 'ctd_cargo.tsv'), 'w') as oh:
        oh.write('mag\tgenus\tt9ss_machinery\tctd_proteins\tprotein_ids\n')
        for m in sorted(mags):
            ps = sorted(ctd.get(m, ()))
            oh.write(f"{m}\t{gen.get(m, '')}\t{'yes' if m in mach else 'no'}\t"
                     f"{len(ps)}\t{','.join(ps)}\n")

rule('WHAT TO REPORT')
if not CTRL_OK:
    say("  >> NO VERDICT. The positive control failed; see A1.")
elif NOTSEARCHED and len(SEARCHED) < 2:
    say("  >> NO VERDICT. Almost nothing was searched; see A2/A3.")
else:
    if DEGRADED:
        say("  0. The control that passed was the DEGRADED union form — the domain")
        say("     metadata would not join. Fix that join before this goes in a")
        say("     methods section, because the domain-declared control is the one")
        say("     that catches the 19b/23e error.")
    if NOTSEARCHED:
        say(f"  0b. {len(NOTSEARCHED)} families had no profile and were NOT SEARCHED.")
        say("     Nothing about them appears as an absence anywhere above.")
    if closed + reopened + still == 0:
        say("  1. GAP A CLOSED NOTHING. No transporter family was both searched and")
        say("     matched to an open edge, so step 29's 193 edges stand exactly as")
        say("     step 29 left them. Fetch the missing profiles and re-run before")
        say("     any of this is quoted — a zero here is an unrun search, not a")
        say("     finding.")
    else:
        say("  1. Step 29's transporter calls were annotation-limited. Direct profile")
        say("     search at KEGG's curated thresholds decides edges that step 29 could")
        say("     only leave open, and the numbers in A4 are the closure.")
    say("  2. Report DOWNGRADED edges as downgraded, not deleted. Step 21's")
    say("     prediction stands; what is missing is the uptake route, in a genome")
    say("     of stated completeness.")
    say("  3. The rejected-hit count in A2 is worth a Methods sentence: it is the")
    say("     size of the error a flat E-value cutoff would have introduced, in")
    say("     this dataset, measured rather than asserted.")
    # item 4 must follow section B's verdict. v3 printed "controlled by
    # confinement" whenever the models were merely PRESENT, so it said the
    # counts were controlled in the same run where section B printed FAIL.
    if CTD_VERDICT == 'PASS':
        say("  4. The CTD counts are a FLOOR on T9SS cargo, controlled by confinement")
        say("     to machinery-bearing genomes.")
    elif CTD_VERDICT == 'JOINFAIL':
        say("  4. T9SS CTDs were searched and found, but the confinement control did")
        say("     NOT run: the machinery table and the proteome set do not share a")
        say("     name space. Quote the CTD counts only alongside the phylum")
        say("     breakdown in section B, and say the machinery control is pending.")
    elif CTD_VERDICT == 'FAIL':
        say("  4. T9SS CTDs do not confine themselves to machinery-bearing genomes.")
        say("     The model is cross-hitting; report the search as done and the cargo")
        say("     floor as NOT established. Do not quote the counts.")
    elif CTD_VERDICT == 'NOCTRL':
        say("  4. T9SS CTDs were searched, but no machinery call was available, so")
        say("     the counts are unvalidated. Report them as uncontrolled or not at")
        say("     all.")
    else:
        say("  4. T9SS cargo remains unmeasured. Say 'not searched', never 'none'.")
say()
say("  SUPERSEDED, and recorded so it is not re-run: the roadmap's dN/dS and")
say("  pangenome-openness step. 32c measured dN/dS = 0.08 on the only strain pair")
say("  available and 32d then showed that pair has no resolvable divergence at")
say("  all; with 28 species in 29 MAGs there is no within-species set to open a")
say("  pangenome on. Both are limits of the sample, not omissions.")

open(os.path.join(D, 'VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(D, 'VERDICT.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "hmm tables written" "$([ "$(ls -1 "$D/hits"/*.tbl 2>/dev/null | wc -l)" -ge 5 ] && echo 0 || echo 1)" || OK=1
check "transporter matrix" "$([ -s "$D/transporter_matrix.tsv" ] && echo 0 || echo 1)" || OK=1
check "verdict written" "$([ -s "$D/VERDICT.txt" ] && echo 0 || echo 1)" || OK=1
if grep -q 'NO VERDICT' "$D/VERDICT.txt" 2>/dev/null; then
  echo "  [STOP] control failed — nothing in $D may be reported."; OK=1
fi
[ "$OK" -eq 0 ] && runlog "36_gapclose complete" || runlog "**36_gapclose FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
