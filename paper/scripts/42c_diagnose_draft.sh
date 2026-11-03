#!/usr/bin/env bash
# 42c -- three genomes failed at `draft`: exit=0 in ~10s, no -draft.RDS. One of them is
# the only archaeon. Diagnose before guessing.
#
# Two errors of mine in 42b, for the record:
#  * the flag detector reported draft/fill as rejecting flags they accept. `draft -h`
#    is not a help handler (it crashes), so I parsed an R traceback as help text; and
#    fill spells options `-m|--model`, which my regex could not match. The detector
#    only ever gated find-transport. draft and fill ran on hardcoded flags that were
#    right by luck.
#  * an exit code of 0 with no output file is the failure mode this whole project keeps
#    hitting. The step check caught it, which is why nothing downstream used a
#    half-built model.
set -uo pipefail
WS2="${WS2:-${WS2}}"
GS="${GS:-$WS2/envs/gapseq_local}"; OUT="$WS2/42_gems"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "42c  $(date -Is)"
PIXI=$(find "$HOME/Documents/UniqueELT_Solid" -maxdepth 6 -type d -name default -path '*envs*' \
       2>/dev/null | grep -v damaged | while read -r e; do [ -x "$e/bin/gapseq" ] && echo "$e"; done)
[ -n "$PIXI" ] && export PATH="$PIXI/bin:$PATH"
FAILED=$(while IFS=$'\t' read -r w d f; do
  [ -z "$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null)" ] \
    && echo "$w"$'\t'"$d"$'\t'"$f"; done < "$OUT/targets.tsv")
note "genomes without a model:"; printf '%s\n' "$FAILED" | sed 's/^/           /'

say "A. draft's REAL option list, from its source"
SRC=$(find "$GS/src" -name 'generate_GSdraft.R' 2>/dev/null | head -1)
if [ -n "$SRC" ]; then
  ok "$SRC"
  grep -nE "make_option|c\(\"-" "$SRC" | head -40 | sed 's/^/           /'
else
  no "generate_GSdraft.R not found; falling back to running draft with no args"
  "$GS/gapseq" draft 2>&1 | head -25 | sed 's/^/           /'
fi

say "B. what the failing runs actually printed"
printf '%s\n' "$FAILED" | while IFS=$'\t' read -r w d f; do
  [ -z "$w" ] && continue
  printf '\n  --- %s (%s) ---\n' "$w" "$d"
  awk '/draft start/{c=1} c' "$OUT/$w/gapseq.log" 2>/dev/null | grep -vE '^\[' | head -25 \
    | sed 's/^/           /'
  note "inputs present:"
  ls -la "$OUT/$w"/*-all-Reactions.tbl "$OUT/$w"/*-Transporter.tbl 2>/dev/null \
    | awk '{printf "           %10s  %s\n",$5,$9}'
  R=$(find "$OUT/$w" -name '*-all-Reactions.tbl' | head -1)
  [ -n "$R" ] && printf '           Reactions.tbl rows: %s\n' "$(( $(wc -l < "$R") - 1 ))"
  T=$(find "$OUT/$w" -name '*-Transporter.tbl' | head -1)
  [ -n "$T" ] && printf '           Transporter.tbl rows: %s\n' "$(( $(wc -l < "$T") - 1 ))"
done

say "C. compare against a genome that WORKED"
GOOD=$(while IFS=$'\t' read -r w d f; do
  [ -n "$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null)" ] \
    && { echo "$w"; break; }; done < "$OUT/targets.tsv")
note "working genome: $GOOD"
for g in $GOOD $(printf '%s\n' "$FAILED" | cut -f1); do
  [ -z "$g" ] && continue
  R=$(find "$OUT/$g" -name '*-all-Reactions.tbl' | head -1)
  T=$(find "$OUT/$g" -name '*-Transporter.tbl' | head -1)
  printf '           %-24s reactions=%-7s transporters=%-6s draft=%s\n' "$g" \
    "$( [ -n "$R" ] && wc -l < "$R" || echo -)" \
    "$( [ -n "$T" ] && wc -l < "$T" || echo -)" \
    "$( [ -n "$(find "$OUT/$g" -name '*-draft.RDS' -size +0 2>/dev/null)" ] && echo yes || echo NO)"
done

say "D. rerun draft verbosely on one failure"
W=$(printf '%s\n' "$FAILED" | head -1 | cut -f1)
F=$(printf '%s\n' "$FAILED" | head -1 | cut -f3)
if [ -n "$W" ]; then
  B=$(basename "$F"); B="${B%.*}"
  note "cd $OUT/$W && gapseq draft -r ... -t ... -c ... -u 200 -l 100"
  ( cd "$OUT/$W" && timeout 600 "$GS/gapseq" draft \
      -r "$OUT/$W/$B-all-Reactions.tbl" -t "$OUT/$W/$B-Transporter.tbl" \
      -c "$F" -u 200 -l 100 2>&1 | tail -30 | sed 's/^/           /' )
  note "exit=$?"
fi

say "E. are the seven existing models sound?"
note "two overlapping runs wrote into these directories after the SSH drop; verify."
PY="$WS2/envs/smetana/bin/python"
[ -x "$PY" ] && while IFS=$'\t' read -r w d f; do
  MX=$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null | head -1)
  [ -z "$MX" ] && continue
  "$PY" - "$MX" "$w" "$d" <<'PY' 2>/dev/null | grep -v Warning
import sys,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
p,w,d=sys.argv[1],sys.argv[2],sys.argv[3]
try:
    m=load_cbmodel(p,flavor='fbc2')
    ex=[r for r in m.reactions if 'EX_' in r]
    s=FBA(m); g=getattr(s,'fobj',None)
    print(f"  [{'OK' if (g or 0)>1e-9 else 'FAIL'}]   {w:<24}{d:<9}{len(m.reactions):>5} rxn "
          f"{len(ex):>4} exch   growth={g!r}")
except Exception as e:
    print(f"  [FAIL] {w:<24}{d:<9}unreadable: {type(e).__name__}")
PY
done < "$OUT/targets.tsv"
note "growth here is on an UNCONSTRAINED model, so it only shows the model is not broken."
note "Growth on ETO_SED is a separate question and comes after all ten exist."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
