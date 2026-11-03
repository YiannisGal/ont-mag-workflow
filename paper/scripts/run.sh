#!/usr/bin/env bash
# run.sh -- run a script only if its contents match the md5 in its own filename.
#
#   bash scripts/run.sh scripts/53_mcra_nomination_a860e20.sh
#   bash scripts/run.sh scripts/53_mcra_nomination_a860e20.sh REBUILD=1
#
# WHY. Every script in this project is named <name>_<first 7 of its md5>.sh. That
# convention only helps if something checks it. Five times in this project a
# script has been transferred or pasted incompletely and the WRONG FILE RAN --
# most recently 52's contents sitting in a file named 53, which produced a
# perfectly plausible-looking output for the wrong analysis. An output you cannot
# tell is from the wrong script is worse than no output.
#
# This refuses to run when the file's md5 does not match its filename, prints
# both, and exits non-zero. It adds two seconds and removes an entire class of
# wasted run.

set -uo pipefail

F="${1:-}"
if [ -z "$F" ]; then
  echo "usage: bash scripts/run.sh <script_name_xxxxxxx.sh> [VAR=value ...]"
  exit 2
fi
shift

if [ ! -f "$F" ]; then
  echo "  [FAIL] no such file: $F"
  exit 1
fi

WANT=$(basename "$F" | sed -n 's/.*_\([0-9a-f]\{7\}\)\.sh$/\1/p')
GOT=$(md5sum "$F" | cut -c1-7)
FULL=$(md5sum "$F" | cut -d' ' -f1)
SIZE=$(stat -c %s "$F")
LINES=$(wc -l < "$F")

echo "================================================================"
echo "  file   : $F"
echo "  size   : $SIZE bytes, $LINES lines"
echo "  md5    : $FULL"
if [ -z "$WANT" ]; then
  echo "  stamp  : (filename carries no md5 stamp)"
  echo "================================================================"
  echo "  [--]   cannot verify this file. Running anyway."
else
  echo "  stamp  : $WANT   (from the filename)"
  echo "================================================================"
  if [ "$WANT" != "$GOT" ]; then
    echo "  [FAIL] CONTENTS DO NOT MATCH THE FILENAME."
    echo "         filename claims  $WANT"
    echo "         file actually is $GOT"
    echo
    echo "         The transfer was incomplete, or an older script is sitting"
    echo "         under this name. Re-transfer it and run this again."
    echo "         Nothing has been executed."
    exit 1
  fi
  echo "  [OK]   verified. Running."
fi

# first banner line of the script, so the log says what it thinks it is
echo "  first line of its own banner:"
grep -m1 -E '^echo "[0-9]+ ' "$F" | sed 's/^/    /' || true
echo "================================================================"
echo

exec env "$@" bash "$F"
