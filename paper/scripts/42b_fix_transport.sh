#!/usr/bin/env bash
# 42b -- find-transport, draft and fill failed to start: exit=1 in 0s on all ten.
#
# I passed `-O -K 3` to find-transport because those are `find`'s flags. 38b verified
# `find`'s flags against its own -h and I then assumed the other subcommands take the
# same options. They plainly do not. So: read -h for EVERY remaining subcommand, build
# each command line only from flags that subcommand actually documents, and run.
#
# `find` is already complete for all ten (3.5-3.7 h each) and is never re-run here.
set -uo pipefail
WS2="${WS2:-${WS2}}"
GS="${GS:-$WS2/envs/gapseq_local}"; OUT="$WS2/42_gems"
FILLMED="$WS2/41_env_checks/ETO_FILL-medium.csv"
export TMPDIR="${TMPDIR42:-${SCRATCH_ROOT}/tmp42}"; mkdir -p "$TMPDIR"
STAGE_TIMEOUT="${STAGE_TIMEOUT:-21600}"; JOBS="${JOBS:-8}"; KTHREADS="${KTHREADS:-3}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "42b  $(date -Is)"
PIXI=$(find "$HOME/Documents/UniqueELT_Solid" -maxdepth 6 -type d -name default -path '*envs*' \
       2>/dev/null | grep -v damaged | while read -r e; do [ -x "$e/bin/gapseq" ] && echo "$e"; done)
[ -n "$PIXI" ] && export PATH="$PIXI/bin:$PATH"

say "A. what the failure actually was"
G=$(head -1 "$OUT/targets.tsv" | cut -f1); FA=$(head -1 "$OUT/targets.tsv" | cut -f3)
note "reproducing on $G:"
( cd "$OUT/$G" && timeout 120 "$GS/gapseq" find-transport -O -K "$KTHREADS" "$FA" 2>&1 | head -6 | sed 's/^/           /' )
note "exit was $?"

say "B. read -h for every remaining subcommand"
for sub in find-transport draft fill; do
  printf '\n  --- gapseq %s -h ---\n' "$sub"
  "$GS/gapseq" "$sub" -h 2>&1 | head -30 | sed 's/^/           /'
done

say "C. build each command line from flags the subcommand documents"
HELP_T=$("$GS/gapseq" find-transport -h 2>&1)
HELP_D=$("$GS/gapseq" draft -h 2>&1)
HELP_F=$("$GS/gapseq" fill -h 2>&1)
has(){ printf '%s' "$1" | grep -qE "(^|[[:space:],])-$2([[:space:]]|,|$)"; }
TARGS=""
has "$HELP_T" O && { TARGS="$TARGS -O"; ok "find-transport accepts -O"; } || note "find-transport has NO -O"
has "$HELP_T" K && { TARGS="$TARGS -K $KTHREADS"; ok "find-transport accepts -K"; } || note "find-transport has NO -K"
has "$HELP_T" b && note "find-transport accepts -b (bitscore) -- not passing, default kept"
note "find-transport args:${TARGS:- (none)}"
for f in r t c u l; do has "$HELP_D" $f && printf '           draft accepts -%s\n' "$f" \
  || printf '           draft does NOT accept -%s  <-- would break\n' "$f"; done
for f in m n c g; do has "$HELP_F" $f && printf '           fill accepts -%s\n' "$f" \
  || printf '           fill does NOT accept -%s  <-- would break\n' "$f"; done

say "D. run transport -> draft -> fill for all ten"
cat > "$OUT/rest.sh" <<'INNER'
#!/usr/bin/env bash
set -uo pipefail
W="$1"; FA="$2"; GS="$3"; OUT="$4"; MED="$5"; T="$6"; TARGS="$7"
D="$OUT/$W"; cd "$D" || exit 1; L="$D/gapseq.log"
B=$(basename "$FA"); B="${B%.*}"
run(){ local nm="$1" chk="$2"; shift 2
  if [ -s "$chk" ]; then echo "[$W] skip $nm" | tee -a "$L"; return 0; fi
  echo "[$W] $nm start $(date -Is)" | tee -a "$L"; local t0=$SECONDS
  timeout "$T" "$@" >>"$L" 2>&1; local rc=$?
  echo "[$W] $nm exit=$rc $((SECONDS-t0))s" | tee -a "$L"
  [ $rc -ne 0 ] && return $rc
  [ -s "$chk" ] || { echo "[$W] $nm wrote no $chk" | tee -a "$L"; return 3; }; return 0; }
# shellcheck disable=SC2086
run transport "$D/$B-Transporter.tbl" "$GS/gapseq" find-transport $TARGS "$FA" || exit 11
run draft "$D/$B-draft.RDS" "$GS/gapseq" draft -r "$D/$B-all-Reactions.tbl" \
    -t "$D/$B-Transporter.tbl" -c "$FA" -u 200 -l 100 || exit 12
model(){ find "$D" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1; }
if [ -n "$(model)" ]; then echo "[$W] skip fill" | tee -a "$L"; else
  echo "[$W] fill start $(date -Is)" | tee -a "$L"; t0=$SECONDS
  timeout "$T" "$GS/gapseq" fill -m "$D/$B-draft.RDS" -n "$MED" \
      -c "$D/$B-rxnWeights.RDS" -g "$D/$B-rxnXgenes.RDS" >>"$L" 2>&1; rc=$?
  echo "[$W] fill exit=$rc $((SECONDS-t0))s" | tee -a "$L"
  [ $rc -ne 0 ] && exit 13
  [ -z "$(model)" ] && { echo "[$W] fill wrote no .xml; dir:" | tee -a "$L"; ls -la "$D" | tee -a "$L"; exit 13; }
fi
echo "[$W] COMPLETE model=$(basename "$(model)")" | tee -a "$L"
INNER
chmod +x "$OUT/rest.sh"
awk -F'\t' -v g="$GS" -v o="$OUT" -v m="$FILLMED" -v t="$STAGE_TIMEOUT" -v a="$TARGS" \
  '{print $1"\t"$3"\t"g"\t"o"\t"m"\t"t"\t"a}' "$OUT/targets.tsv" \
  | xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -ra x <<< "{}"; '"$OUT"'/rest.sh "${x[@]}"'

say "E. result"
D=0; T=0
while IFS=$'\t' read -r w dom f; do
  T=$((T+1)); MX=$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1)
  if [ -n "$MX" ]; then D=$((D+1))
    R=$(grep -c '<reaction ' "$MX" 2>/dev/null || echo 0)
    printf '  [OK]   %-24s %-8s %s reactions\n' "$w" "$dom" "$R"
  else
    printf '  [FAIL] %-24s %-8s no model:\n' "$w" "$dom"
    tail -3 "$OUT/$w/gapseq.log" | sed 's/^/             /'
  fi
done < "$OUT/targets.tsv"
echo "           $D/$T"
note "compare the Archaea row against the Bacteria: a much smaller reaction count would"
note "mean gapseq's archaeal pathway coverage is thin, as Bakta's was (41f, 95%)."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
