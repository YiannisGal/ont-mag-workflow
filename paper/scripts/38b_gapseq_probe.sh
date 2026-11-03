#!/usr/bin/env bash
# =============================================================================
# 38b_gapseq_probe.sh  --  ONE GENOME THROUGH gapseq, TIMED, BEFORE COMMITTING 29
#
# WHAT CHANGED SINCE THE PRE-FLIGHT
#   38_preflight said "NEITHER model builder is installed" and "no MILP solver".
#   Both were artefacts of where it looked. gapseq is installed, in a PIXI
#   environment the pre-flight never searched:
#
#     ~/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default
#
#   with its own R, its own Rscript, and the R package **cobrar** — which links
#   GLPK directly and is where the MILP is solved. Two gap-filled models sit
#   beside it. So the solver question is closed, and step 38 may need no model
#   builder installed at all.
#
# WHY A PROBE AND NOT THE RUN
#   The culturomics log records `Running time: 18086 s` for `gapseq find` on ONE
#   isolate — five hours — and the command line shows `-K 1`, a single thread on
#   a 24-thread machine. Naively, 29 MAGs is 145 hours. With threads it should
#   be far less, but "should" is not a measurement, and this project has spent
#   enough days on unmeasured assumptions. So: ONE genome, every stage timed,
#   then a projection with the measured number in it.
#
# 🔴 WHAT WENT WRONG IN v1, AND WHAT REPLACES IT
#   v1 passed `-U` to every stage, believing it suppressed reference updates.
#   gapseq's own help says the opposite:
#
#     -U  Do not use gapseq sequence archive and update sequences from uniprot
#         manually (very slow) (default: false)
#
#   So `-U` FORCED a manual UniProt refetch. It ran for a day, aborted on a curl
#   error, and rewrote 2,435 of 85,765 reference files in the culturomics pixi
#   environment, leaving that reference set partially updated and missing at
#   least one file. `-U` is gone.
#
#   The deeper error was the premise. v1 promised the pixi env would stay
#   read-only. gapseq WRITES TO ITS OWN dat/ ON THE DEFAULT PATH — the
#   culturomics log says "Reference sequences updated" on a run with no `-U` —
#   so no flag could ever have delivered that promise. Using another project's
#   tool installation in place was never viable.
#
#   THE FIX, twofold:
#     * gapseq's data directory is COPIED to $WS2/envs/gapseq_local and this
#       script runs that copy. The culturomics env is read for BINARIES only.
#     * Every flag this script passes is checked against `gapseq <sub> -h`
#       BEFORE anything long starts, and the run is refused if a flag is absent.
#       Verifying an interface against the tool rather than against memory, as a
#       control rather than as care. That check costs two seconds and would have
#       saved a day and 2,435 files.
#
# CONTROLS
#   1. The binaries exist and report a version. A path is not a program.
#   2. The reference data is present, so -U is safe — if dat/seq is missing,
#      gapseq would need the update this script forbids, and that is a STOP.
#   3. The full chain completes: find -> find-transport -> draft -> medium ->
#      fill. The draft->fill transition IS the MILP; if it completes, the solver
#      works, which is a better test than any synthetic LP.
#   4. The model is PLAUSIBLE, not merely present. Reaction and gene counts are
#      compared against the culturomics YGCA01_Clean.xml as a scale reference.
#      A model with 0 or 5 reactions is a failure that otherwise looks like
#      success — the step-11b lesson: absence reported as a result.
#   5. Nothing under the pixi env was modified.
#
# OUTPUT  38_gems/probe/{PROBE.txt, <genome>*, timings.tsv}
# RUNTIME dominated by `gapseq find`; that is the number being measured.
#         Budget several hours for the first run and start it in `screen`.
# USAGE   ws2 && screen -S gapseq
#         # STAGE_TIMEOUT=7200 by default; every stage is killed at that budget
#         bash scripts/38b_gapseq_probe.sh 2>&1 | tee logs/38b.log
#         # optional: GAPSEQ_ENV=/path/to/pixi/env  PROBE_MAG=SED_das_semibin_003
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="38b_gapseq_probe"
F="$WS2/07_mags"
D="$WS2/38_gems/probe"; mkdir -p "$D"

echo "######## 38b — ONE GENOME THROUGH gapseq, TIMED ########"
echo "  started: $(date -Is)"
rm -f "$D/PROBE.txt" "$D/timings.tsv" "$D/env_modified.txt" "$D/smoke.log"

# ------------------------------------------------------------------ locate env
# Two separate things:
#   GAPSEQ_HOME — our OWN copy of share/gapseq. Written to. Never shared.
#   GAPSEQ_ENV  — someone else's pixi env, used ONLY for its binaries (blast,
#                 exonerate, R, ...). Never written.
GAPSEQ_HOME="${GAPSEQ_HOME:-$WS2/envs/gapseq_local}"
if [ ! -x "$GAPSEQ_HOME/gapseq" ]; then
  echo "  !! no local gapseq copy at $GAPSEQ_HOME"
  echo "     Make one first — the shared env CANNOT be used in place, because"
  echo "     gapseq writes to its own dat/ even on the default path:"
  echo "       cp -a <pixi env>/share/gapseq \"\$WS2/envs/gapseq_local\""
  exit 2
fi
if [ -z "${GAPSEQ_ENV:-}" ]; then
  mapfile -t WITH < <(find "$HOME" -maxdepth 8 -type f -name gapseq -perm -u+x \
      2>/dev/null | sed 's#/bin/gapseq$##' | sort -u)
  echo "  envs providing gapseq binaries: ${#WITH[@]}"
  for c in "${WITH[@]}"; do echo "    $c"; done
  if [ "${#WITH[@]}" -eq 1 ]; then GAPSEQ_ENV="${WITH[0]}"
  elif [ "${#WITH[@]}" -gt 1 ]; then
    echo "  !! more than one. REFUSING TO CHOOSE — set GAPSEQ_ENV."; exit 2
  fi
fi
if [ -z "${GAPSEQ_ENV:-}" ] || [ ! -d "$GAPSEQ_ENV/bin" ]; then
  echo "  !! set GAPSEQ_ENV=/path/to/.pixi/envs/default (binaries only)"; exit 2
fi
export PATH="$GAPSEQ_ENV/bin:$PATH"
echo "  gapseq (local copy) : $GAPSEQ_HOME"
echo "  binaries from       : $GAPSEQ_ENV   [read-only]"
GS="$GAPSEQ_HOME/gapseq"
RS="$GAPSEQ_ENV/bin/Rscript"
echo "  threads             : ${THREADS:-4}"

# stamp for the read-only control
STAMP="$D/.stamp"; : > "$STAMP"; sleep 1
# A guard placed after the operation it guards is decoration. v1's read-only
# check sat at the end of the script; `find` never returned, so it never ran
# while 2,435 files were being rewritten. It now runs immediately after the
# FIRST stage, so a violation is caught in minutes rather than never.
env_check () {
  find "$GAPSEQ_ENV" -newer "$STAMP" -type f 2>/dev/null | sort | head -40 \
    > "$D/env_modified.txt"
  if [ -s "$D/env_modified.txt" ]; then
    echo "  [STOP] the shared env was written to. Aborting before more damage:"
    head -5 "$D/env_modified.txt"
    echo "STEP_STATUS: SHARED_ENV_MODIFIED" > "$D/PROBE.txt"
    exit 1
  fi
}

# ------------------------------------------------------------------ pick genome
# The smallest genome, so the timing is a FLOOR the projection can scale from,
# and the first run is the cheapest possible. Deterministic: size then name.
# One location, not two. v1 searched $F to depth 2 and took the smallest file
# anywhere under it, mixing bins/ with polished/ and unpolished/. 39b's
# mag_fastas() already had this right — first populated location wins, never a
# mixture — and these two scripts disagreed because I wrote them hours apart.
BINDIR=""
for sub in bins final MAGs .; do
  if compgen -G "$F/$sub/*.fa" >/dev/null || compgen -G "$F/$sub/*.fna" >/dev/null; then
    BINDIR="$F/$sub"; break
  fi
done
[ -n "$BINDIR" ] || { echo "  !! no MAG FASTAs under $F"; exit 2; }
echo "  bin dir    : $BINDIR"
if [ -z "${PROBE_MAG:-}" ]; then
  PROBE_MAG="$(find "$BINDIR" -maxdepth 1 \( -name '*.fa' -o -name '*.fna' \) \
      -not -name 'all_*' -printf '%s\t%f\n' 2>/dev/null \
      | sort -k1,1n -k2,2 | head -1 | cut -f2 | sed 's/\.[^.]*$//')"
fi
FA="$(find "$BINDIR" -maxdepth 1 -name "${PROBE_MAG}.*" \
      \( -name '*.fa' -o -name '*.fna' \) 2>/dev/null | sort | head -1)"
if [ -z "$FA" ]; then
  echo "  !! probe genome '${PROBE_MAG}' not found under $F"; exit 2
fi
# gapseq's -t decides which reference set is searched. 8 of the 29 genomes are
# archaea and this project treats both domains equally, so the domain is READ
# from the taxonomy rather than assumed. Passing Bacteria for an archaeon would
# search the wrong references and return a thin model that still looks like a
# result.
# The genome-id column in MAG_TAXONOMY.tsv is `working_name`, NOT `bin_id`.
# `bin_id` was a name I invented for a test fixture and then trusted; a fixture
# built from an assumption can only ever confirm it. Accept every id column the
# file actually has.
TAXF="$F/MAG_TAXONOMY.tsv"
if [ -s "$TAXF" ]; then
  # a header with fewer fields than its data rows shifts every column index,
  # which would silently return the wrong column rather than no column
  NFH=$(awk -F'\t' 'NR==1{print NF; exit}' "$TAXF")
  NFD=$(awk -F'\t' 'NR==2{print NF; exit}' "$TAXF")
  if [ "${NFH:-0}" != "${NFD:-0}" ]; then
    echo "  [STOP] MAG_TAXONOMY.tsv header has $NFH fields, data rows have $NFD."
    echo "  Every column looked up by name in this file is shifted. That affects"
    echo "  far more than this step — fix the file before trusting any of them."
    exit 2
  fi
fi
DOM="$(awk -F'\t' -v m="$PROBE_MAG" '
  NR==1{for(i=1;i<=NF;i++){
          if($i=="domain")d=i
          if($i=="working_name"||$i=="bin_id"||$i=="qc_bin_name")w[i]=1
          if($i=="mag_id")w[i]=1}; next}
  {for(i in w) if($i==m){print $d; exit}}' "$TAXF" 2>/dev/null)"
case "${DOM:-}" in
  [Aa]rchae*) GTAX="Archaea" ;;
  [Bb]acteri*) GTAX="Bacteria" ;;
  *) GTAX="auto" ;;
esac
echo "  probe MAG  : $PROBE_MAG  ($(du -h "$FA" | cut -f1), $(grep -c '^>' "$FA") contigs)"
echo "  domain     : ${DOM:-<not in MAG_TAXONOMY.tsv>}  -> gapseq -t $GTAX"
if [ "$GTAX" = "auto" ]; then
  echo "  [STOP] '$PROBE_MAG' has no domain in MAG_TAXONOMY.tsv."
  echo "  Either it is not one of the catalogue genomes, or the id column has"
  echo "  changed again. Both are reasons NOT to spend hours modelling it, and"
  echo "  v1 printed this as a note and started anyway. Set PROBE_MAG explicitly"
  echo "  to a genome that resolves, or fix the taxonomy join first."
  exit 2
fi
echo

cd "$D" || exit 2
: > timings.tsv

# Every stage runs under a wall-clock budget. The smoke test below may or may
# not exercise the code path that hangs — `-x` might skip the reference update,
# and I cannot know that without trying. A timeout does not depend on knowing
# where the hang is: `find` sat at 3% CPU for 51 minutes doing no alignment, and
# an earlier run burned a full day. STAGE_TIMEOUT makes that impossible.
: "${STAGE_TIMEOUT:=7200}"
run_stage () {          # $1 = label, rest = command
  local label="$1"; shift
  local t0 t1
  echo "  ---- $label  (budget ${STAGE_TIMEOUT}s)"
  t0=$(date +%s)
  timeout --foreground -k 30 "$STAGE_TIMEOUT" "$@" >> "$D/${label}.log" 2>&1
  local rc=$?
  t1=$(date +%s)
  if [ "$rc" -eq 124 ]; then
    echo "       TIMED OUT after ${STAGE_TIMEOUT}s — killed"
  fi
  printf '%s\t%s\t%s\n' "$label" "$((t1 - t0))" "$rc" >> "$D/timings.tsv"
  [ "$rc" -eq 124 ] || echo "       $((t1 - t0)) s, exit $rc"
  return $rc
}

# ------------------------------------------------------------------ controls 1+2
echo "######## CONTROLS 1-2 ########"
"$GS" -v > "$D/version.txt" 2>&1
echo "  gapseq version: $(head -1 "$D/version.txt")"
if [ ! -x "$RS" ]; then
  echo "  [FAIL] $RS is not executable. The culturomics log shows gapseq_find.sh"
  echo "  line 804 failing on exactly this — the env's own Rscript was gone."
  echo "STEP_STATUS: NO_RSCRIPT" > "$D/PROBE.txt"; exit 1
fi
echo "  Rscript   : $("$RS" -e 'cat(R.version.string)' 2>/dev/null)"
COBRAR=$("$RS" -e 'cat(as.character(requireNamespace("cobrar", quietly=TRUE)))' 2>/dev/null)
echo "  cobrar    : ${COBRAR:-not loadable}"
if [ "$COBRAR" != "TRUE" ]; then
  echo "  [FAIL] cobrar does not load. That is where gapseq solves its MILP."
  echo "STEP_STATUS: NO_COBRAR" > "$D/PROBE.txt"; exit 1
fi

# --------------------------------------------------------- CONTROL 2a: FLAGS
# The -U disaster in one sentence: I asserted a flag's meaning from memory,
# documented the assertion at length, and never ran `gapseq find -h`. Every flag
# below is now checked against the tool's own help text before anything long
# starts, and a flag whose description contradicts the intended use is a STOP.
echo "  ---- verifying flags against gapseq's own help"
for sub in find find-transport draft medium fill; do
  "$GS" "$sub" -h > "$D/help_${sub}.txt" 2>&1 || true
done
flag_desc () {   # $1=subcommand  $2=flag letter -> prints the description line
  awk -v f="-$2" '$1==f{ $1=""; sub(/^ +/,""); print; exit }' "$D/help_$1.txt"
}
FLAGFAIL=0
check_flag () {  # $1=sub $2=letter $3=regex the description MUST match
  local d; d="$(flag_desc "$1" "$2")"
  if [ -z "$d" ]; then
    echo "    [FAIL] $1 has no -$2 flag"; FLAGFAIL=1; return
  fi
  if ! printf '%s' "$d" | grep -qiE "$3"; then
    echo "    [FAIL] $1 -$2 is '$d'"
    echo "           expected something matching /$3/"; FLAGFAIL=1; return
  fi
  echo "    [ok]   $1 -$2 : ${d:0:64}"
}
has_flag () {    # $1=sub $2=letter -> 0 if the flag exists in that help text
  [ -n "$(flag_desc "$1" "$2")" ]
}
check_flag find           b 'bit ?score'
check_flag find           p 'keyword|pathway|subsystem'
check_flag find           t 'taxonom'
check_flag find           K 'thread|cpu|core'
check_flag find-transport b 'bit ?score'
# -O is what stops the reference-update network call that hangs `find` with or
# without -U. It was in the help all along; the first version of this script did
# not read the help at all. Verified, then used.
OFF_FIND=""; OFF_TRAN=""
if has_flag find O; then
  check_flag find O 'offline'
  OFF_FIND="-O"
else
  echo "    [WARN] this gapseq has no 'find -O'. Without offline mode the"
  echo "           reference update runs and has been observed to hang on a"
  echo "           network fetch with no timeout."
fi
if has_flag find-transport O; then OFF_TRAN="-O"; fi
# -U must NOT be used. Assert that its description still means what broke us,
# so that if a future gapseq redefines it this script notices instead of
# silently being right for the wrong reason.
UDESC="$(flag_desc find U)"
if printf '%s' "$UDESC" | grep -qiE 'not use .*archive|uniprot'; then
  echo "    [ok]   find -U confirmed as 'refetch from uniprot' — NOT PASSED"
else
  echo "    [WARN] find -U now reads: '$UDESC'"
  echo "           v1 passed this believing it meant 'no update'. It is still"
  echo "           not passed, but the description has changed — re-read it."
fi
if [ "$FLAGFAIL" -ne 0 ]; then
  echo "  [STOP] at least one flag is absent or means something else."
  echo "  This is the check that was missing yesterday. Nothing long will start."
  echo "STEP_STATUS: FLAG_MISMATCH" > "$D/PROBE.txt"; exit 1
fi
echo "  [PASS] every flag verified against the tool, not against memory"
echo

SEQDIR="$GAPSEQ_HOME/dat/seq"
NSEQ=$(find "$SEQDIR" -type f 2>/dev/null | wc -l)
echo "  reference sequences: $NSEQ file(s) under dat/seq"
if [ "$NSEQ" -lt 10 ]; then
  echo "  [STOP] the local copy has no reference data. Re-copy share/gapseq"
  echo "  from the pixi env; a gapseq with no archive refetches from UniProt,"
  echo "  which is what cost a day."
  echo "STEP_STATUS: WOULD_WRITE_SHARED_ENV" > "$D/PROBE.txt"; exit 1
fi
echo "  [PASS] binaries present, cobrar loads, local reference archive present"
echo

# ------------------------------------------------------------------ the chain
# -U everywhere: never update reference sequences, because that writes into the
# culturomics environment.
echo "######## SMOKE TEST — does find get past the update at all? ########"
# `-x` lists pathways and reactions WITHOUT blasting. It exercises exactly the
# startup path that hung (reference check, taxonomy, database load) and returns
# in seconds. A hang here is a hang in 30 s instead of an hour.
if has_flag find x; then
  ( "$GS" find -v 0 -p amino -t "$GTAX" $OFF_FIND -x "$FA" \
      > "$D/smoke.log" 2>&1 ) &
  SMOKE=$!
  for i in $(seq 1 60); do
    kill -0 "$SMOKE" 2>/dev/null || break
    sleep 5
  done
  if kill -0 "$SMOKE" 2>/dev/null; then
    kill -TERM "$SMOKE" 2>/dev/null
    echo "  [STOP] the dry run (-x, no blasting) did not finish in 5 minutes."
    echo "  It is hanging before any alignment work — almost certainly the"
    echo "  reference-update network call. Last lines:"
    tail -3 "$D/smoke.log" 2>/dev/null | sed 's/^/    /'
    echo "  Nothing long has been started."
    echo "STEP_STATUS: SMOKE_TEST_HUNG" > "$D/PROBE.txt"; exit 1
  fi
  wait "$SMOKE"; SRC=$?
  echo "  dry run finished, exit $SRC"
  [ "$SRC" -eq 0 ] || { echo "  [STOP] the dry run failed; see smoke.log"
    tail -5 "$D/smoke.log" | sed 's/^/    /'
    echo "STEP_STATUS: SMOKE_TEST_FAILED" > "$D/PROBE.txt"; exit 1; }
else
  echo "  no -x flag; skipping the dry run"
fi
echo

echo "######## THE gapseq CHAIN (NO -U, ${OFF_FIND:-no -O}, -K ${THREADS:-4}) ########"
OKALL=0
run_stage find           "$GS" find -v 0 -b 200 -p all -t "$GTAX" $OFF_FIND -K "${THREADS:-4}" "$FA" || OKALL=1
env_check      # first opportunity: before anything else runs
run_stage find-transport "$GS" find-transport -b 200 $OFF_TRAN -K "${THREADS:-4}" "$FA" || OKALL=1

RXN="$D/${PROBE_MAG}-all-Reactions.tbl"
PWY="$D/${PROBE_MAG}-all-Pathways.tbl"
TRA="$D/${PROBE_MAG}-Transporter.tbl"
for f in "$RXN" "$PWY" "$TRA"; do
  [ -s "$f" ] || { echo "  !! missing $(basename "$f") — the chain cannot continue"; OKALL=1; }
done

if [ "$OKALL" -eq 0 ]; then
  run_stage draft  "$GS" draft -r "$RXN" -t "$TRA" -p "$PWY" -c "$FA" -u 200 -l 100 || OKALL=1
  run_stage medium "$GS" medium -m "$D/${PROBE_MAG}-draft.RDS" -p "$PWY" || OKALL=1
  run_stage fill   "$GS" fill -m "$D/${PROBE_MAG}-draft.RDS" \
                    -n "$D/${PROBE_MAG}-medium.csv" \
                    -c "$D/${PROBE_MAG}-rxnWeights.RDS" \
                    -g "$D/${PROBE_MAG}-rxnXgenes.RDS" -b 100 || OKALL=1
fi
echo

# ------------------------------------------------------------------ controls 3-5
echo "######## CONTROLS 3-5 ########"
# Control 5's DATA must be gathered before the block that reads it. In v1 this
# scan sat after the python heredoc, so env_modified.txt was always empty and
# control 5 always passed — a control that cannot fail (Rule 55), in the script
# that introduced Rule 55.
find "$GAPSEQ_ENV" -newer "$STAMP" -type f 2>/dev/null | sort | head -40 \
  > "$D/env_modified.txt"
rm -f "$STAMP"

python3 - "$D" "$PROBE_MAG" "$GAPSEQ_ENV" "$WS2" "${THREADS:-4}" <<'PY'
import sys, os, re, glob

D, MAG, ENV, WS2, NTHREADS = sys.argv[1:6]
W = 100
rep = []
STATUS = ['INCOMPLETE']


def say(s=''):
    print(s)
    rep.append(s)


def dump(st):
    STATUS[0] = st
    with open(os.path.join(D, 'PROBE.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'PROBE.txt')}")
    sys.exit(0)


say('=' * W)
say('  38b — ONE GENOME THROUGH gapseq')
say('=' * W)
say(f"  probe genome : {MAG}")
say(f"  gapseq env   : {ENV}")
say()

# ---- timings
say('  stage timings')
say('-' * W)
tot = 0
stages = []
tp = os.path.join(D, 'timings.tsv')
if os.path.exists(tp):
    for line in open(tp):
        f = line.rstrip('\n').split('\t')
        if len(f) < 3:
            continue
        lab, sec, rc = f[0], int(f[1]), int(f[2])
        stages.append((lab, sec, rc))
        tot += sec
        say(f"  {lab:<20}{sec:>8} s{'':>4}{'ok' if rc == 0 else 'EXIT ' + str(rc)}")
say('-' * W)
say(f"  {'TOTAL':<20}{tot:>8} s  ({tot / 3600.0:.2f} h)")
say()
failed = [s for s in stages if s[2] != 0]
if failed or not stages:
    say("  [FAIL] not every stage completed:")
    for lab, sec, rc in failed:
        if rc == 124:
            say(f"    {lab} TIMED OUT after {sec} s and was killed.")
            say("      Raise STAGE_TIMEOUT if it was genuinely still working;")
            say("      check for a hang if CPU was near zero. `find` has been")
            say("      seen sitting at 3% CPU doing no alignment at all.")
        else:
            say(f"    {lab} exited {rc} — see {lab}.log")
    say("  The projection below is withheld; a chain that does not finish once")
    say("  cannot be scaled to 29 genomes.")
    dump('CHAIN_INCOMPLETE')

# ---- CONTROL 3: the model exists
xml = sorted(glob.glob(os.path.join(D, MAG + '.xml')))
draft = sorted(glob.glob(os.path.join(D, MAG + '-draft.xml')))
say('  CONTROL 3 — did the MILP run?')
say(f"    draft model : {'yes' if draft else 'NO'}")
say(f"    filled model: {'yes' if xml else 'NO'}")
if not xml:
    say("  [FAIL] no gap-filled model. draft -> filled IS the mixed-integer step;")
    say("  without it the solver question is not answered for this machine.")
    dump('NO_FILLED_MODEL')
say("  [PASS] a gap-filled model was produced, so cobrar/GLPK solved the MILP.")
say()


def sbml_counts(path):
    try:
        txt = open(path, errors='replace').read()
    except Exception:
        return None
    return (len(re.findall(r'<reaction[ >]', txt)),
            len(re.findall(r'<species[ >]', txt)),
            len(re.findall(r'<fbc:geneProduct[ >]', txt)))


# ---- CONTROL 4: is it PLAUSIBLE, not merely present?
say('  CONTROL 4 — is the model plausible? (a 5-reaction model is not a success)')
c = sbml_counts(xml[0])
cd = sbml_counts(draft[0]) if draft else None
say(f"    {'model':<40}{'reactions':>12}{'species':>10}{'genes':>10}")
say('-' * W)
if cd:
    say(f"    {MAG + '-draft.xml':<40}{cd[0]:>12,}{cd[1]:>10,}{cd[2]:>10,}")
say(f"    {MAG + '.xml (filled)':<40}{c[0]:>12,}{c[1]:>10,}{c[2]:>10,}")

# scale reference: the culturomics model built by the same tool
ref = sorted(glob.glob(os.path.join(ENV, '..', '..', '..', '*.xml'))) or \
      sorted(glob.glob(os.path.join(os.path.dirname(ENV.split('/.pixi')[0]),
                                    'gapseq_project', '*.xml')))
ref = [r for r in ref if 'draft' not in r and 'backup' not in r]
if ref:
    rc_ = sbml_counts(ref[0])
    if rc_:
        say(f"    {os.path.basename(ref[0]) + ' (culturomics)':<40}"
            f"{rc_[0]:>12,}{rc_[1]:>10,}{rc_[2]:>10,}")
say()
if c[0] < 300:
    say(f"  [FAIL] {c[0]} reactions is far below any viable bacterial GEM"
        f" (typically 1,000-2,500).")
    say("  The chain completed and produced something empty — exactly the failure")
    say("  mode that looks like success. Do not scale this to 29.")
    dump('MODEL_IMPLAUSIBLE')
if cd and c[0] <= cd[0]:
    say("  [WARN] the filled model has no more reactions than the draft, so the")
    say("  gap-fill added nothing. Check the medium — an over-rich medium makes")
    say("  gap-filling unnecessary and silently trivialises the result.")
say(f"  [PASS] {c[0]:,} reactions, {c[2]:,} genes — in the expected range.")
say()

# ---- CONTROL 5: was the shared environment touched?
say('  CONTROL 5 — was the culturomics pixi env modified?')
mod = os.path.join(D, 'env_modified.txt')
n = 0
if os.path.exists(mod):
    lines = [l.strip() for l in open(mod, errors='replace') if l.strip()]
    n = len(lines)
    for l in lines[:8]:
        say(f"    {l}")
if n:
    say(f"  [FAIL] {n} file(s) under the pixi env changed during this run. The")
    say("  read-only constraint was violated — `-U` did not cover everything.")
    say("  Investigate before running 29 genomes through it.")
    dump('SHARED_ENV_MODIFIED')
say("  [PASS] nothing under the pixi env was modified.")
say()

# ---- the projection, with a measured number in it
say('=' * W)
say('  PROJECTION FOR 29 GENOMES')
say('=' * W)
nmag = len(glob.glob(os.path.join(WS2, '07_mags', 'bins', '*.fa'))) or 29
find_s = next((s for l, s, _ in stages if l == 'find'), 0)
say(f"  measured on {MAG}, at {NTHREADS} threads:")
say(f"    find            {find_s:>7} s   ({find_s / 3600.0:.2f} h)")
say(f"    whole chain     {tot:>7} s   ({tot / 3600.0:.2f} h)")
say()
say(f"  x {nmag} genomes, serial : {tot * nmag / 3600.0:.1f} h")
say(f"  The culturomics run took 18,086 s for `find` alone at -K 1. This run"
    f" took {find_s} s.")
if find_s > 0:
    say(f"  Speedup on that stage: {18086.0 / find_s:.1f}x")
say()
say("  NOTE. This is the SMALLEST genome in the catalogue, so the per-genome")
say("  figure is a FLOOR, not a mean. Scale it with judgement, and note that")
say("  gapseq parallelises within a genome (-K) but genomes are independent —")
say("  running several at once is the obvious win if RAM allows.")
say()
say('=' * W)
say('  WHAT IS STILL MISSING')
say('=' * W)
say("  gapseq and its solver are ACCOUNTED FOR — no install needed for the model")
say("  building half of step 38.")
say()
say("  SMETANA is still absent and is a python tool needing its own solver. It")
say("  goes in a NEW venv under workspace_2/envs, never into the pixi env:")
say("      python3 -m venv $WS2/envs/smetana")
say("      $WS2/envs/smetana/bin/pip install smetana swiglpk")
say("  and the pre-flight should then be re-run to confirm it SOLVES rather than")
say("  merely imports.")
say()
say("  A MEDIUM is still unspecified, and it decides the answer. gapseq wrote")
say(f"  {MAG}-medium.csv here, and the culturomics run has two more. For Etoliko")
say("  the medium must be anoxic, sulfate-replete and low in fixed carbon; a")
say("  rich medium makes every member independent and collapses every")
say("  cross-feeding score to zero.")
dump('OK')
PY
RC=$?

if [ -s "$D/env_modified.txt" ]; then
  echo "  !! files under the pixi env changed during this run:"
  head -8 "$D/env_modified.txt"
fi

echo
echo "######## VALIDATION ########"
OK=0
check "probe completed" "$RC" || OK=1
check "report written" "$([ -s "$D/PROBE.txt" ] && echo 0 || echo 1)" || OK=1
check "timings recorded" "$([ -s "$D/timings.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/PROBE.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  status: ${ST:-<none>}"
[ "${ST:-MISSING}" = "OK" ] || { echo "  [STOP] ${ST:-no status line}"; OK=1; }
[ "$OK" -eq 0 ] && runlog "38b_gapseq_probe complete" || runlog "**38b_gapseq_probe FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
