#!/usr/bin/env bash
# make_figures.sh -- build every thesis figure from MASTER.tsv.
#
#   bash scripts/make_figures.sh              # all figures found
#   bash scripts/make_figures.sh fig02        # just one
#
# Puts PDF (vector, for submission) and PNG (600 dpi, for looking at) in
# 50_figures/. Re-runnable: a figure is rebuilt every time, so when a number
# changes you re-run and every panel that uses it is current.
#
# HOUSE STYLE lives in etoliko_figstyle.py. Change a colour or a font size there
# and every panel changes together -- which is the point. A thesis whose figures
# are each styled separately reads as separate documents.
#
# The palette was validated, not chosen by eye:
#   Bacteria #2a78d6 · Archaea #eb6834 on a #fcfcfb surface
#   CVD separation dE 24.7 (protan) / 32.7 (tritan)   -- target >=8
#   normal-vision dE 33.6                             -- floor  >=15
#   contrast vs surface >=3:1, both slots
# Colour encodes DOMAIN only. Nothing in any panel is identified by colour alone.
set -uo pipefail
WS2="${WS2:-${WS2}}"
SRC="${MASTER:-$WS2/47_master/MASTER.tsv}"
OUT="${OUTDIR:-$WS2/50_figures}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
echo "figures  $(date -Is)"

say "A. a python with matplotlib"
PY=""
for c in "$WS2/envs/smetana/bin/python" \
         "$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default/bin/python" \
         python3 python; do
  command -v "$c" >/dev/null 2>&1 || [ -x "$c" ] || continue
  if "$c" -c "import matplotlib" 2>/dev/null; then
    PY="$c"; V=$("$c" -c "import matplotlib;print(matplotlib.__version__)")
    ok "$c  (matplotlib $V)"; break
  else
    no "$c  -- no matplotlib"
  fi
done
if [ -z "$PY" ]; then
  bad "no interpreter with matplotlib."
  note "install into the venv we own, never a shared conda env:"
  note "  $WS2/envs/smetana/bin/pip install matplotlib"
  echo; echo "STEP_STATUS: FAILED"; exit 1
fi

say "B. inputs"
if [ ! -s "$SRC" ]; then
  bad "no MASTER.tsv at $SRC -- run 47 first"; echo "STEP_STATUS: FAILED"; exit 1
fi
ok "$SRC  ($(( $(wc -l < "$SRC") - 1 )) genomes, $(head -1 "$SRC" | tr '\t' '\n' | wc -l) columns)"
mkdir -p "$OUT"; ok "output: $OUT"

say "C. build"
WANT="${1:-}"
N=0; FAIL=0
for f in "$HERE"/fig*.py; do
  [ -s "$f" ] || continue
  b=$(basename "$f" .py)
  if [ -n "$WANT" ] && [[ "$b" != "$WANT"* ]]; then continue; fi
  printf '  --- %s ---\n' "$b"
  if ( cd "$HERE" && "$PY" "$f" "$SRC" "$OUT" ); then N=$((N+1)); else FAIL=$((FAIL+1)); bad "$b"; fi
done
echo
[ "$N" -gt 0 ] && ok "$N figure(s) built" || no "nothing built"
[ "$FAIL" -gt 0 ] && bad "$FAIL failed"

say "D. what is in $OUT"
ls -lh "$OUT"/*.pdf "$OUT"/*.png 2>/dev/null | awk '{printf "     %-10s %s\n",$5,$9}' || no "empty"
note "PDF for the thesis (vector, fonts embedded as TrueType subsets, not Type 3)."
note "PNG at 600 dpi for checking. Look at the PNG before using the PDF: the palette"
note "is validated by script, but label collisions are only visible to the eye."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
