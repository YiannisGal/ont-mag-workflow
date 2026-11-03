#!/usr/bin/env bash
# 78_housekeeping_fix.sh -- correct three wrong assumptions from step 77.
#
# Step 77 refused both deletions. The refusals were correct and my reasons for
# proposing them were not:
#
#   A. chunks/ is NOT a copy of all_hits.tbl. The files are named c00.hmm and
#      hold 51.4 M lines against all_hits.tbl's 1.5 M; 49.9 M of them are not
#      in it. That is the SPLIT KOfam PROFILE DATABASE, i.e. the search's
#      INPUT, not its output. Comparing its lines to a hit table was comparing
#      two unrelated things -- the gate held only because it demanded identity
#      rather than plausibility.
#
#   B. all_proteins.faa is NOT a concatenation of the per-genome proteins.
#      919,806 sequences against 159,523 in 06_binqc, and the surplus carries
#      ids like ACR|OFV66072.1 -- an external namespace. It holds reference
#      proteins as well as ours, so the per-genome files cannot rebuild it.
#
#   C. the curl conversion it wrote is BROKEN. `wget -qO-` became
#      `curl -fsSLO -qO-`, which is not a command: curl's -O means "write to a
#      remote-named file", the opposite of -O- "write to stdout", and -q has a
#      different meaning again. A converter that emits something unrunnable is
#      worse than no converter, because the breakage surfaces at fetch time.
#      It also tested reachability against UNIPROT while the script fetches
#      from NCBI E-utilities -- a connectivity test against a host the script
#      never contacts.
#
# AND THE PREMISE BEHIND ALL OF IT WAS WRONG. ${SCRATCH_ROOT} has 1.8 TB free.
# There is no disk pressure. Nothing here needs deleting; this step only
# IDENTIFIES what the two large things are, so the next person is not told a
# 7.6 GB directory is a duplicate when it is a database.
set -u
WS2="${WS2}"
K="$WS2/45_kofam_full"
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. what chunks/ actually is, and whether it is regenerable"
echo "=============================================================="
CH="$K/chunks"
if [ -d "$CH" ]; then
  echo "  $(ls -1 "$CH" | wc -l) files, $(du -sh "$CH" | cut -f1)"
  ls -1 "$CH" | head -4 | sed 's/^/     /'
  F=$(ls -1 "$CH" | head -1)
  echo "  first three lines of $F:"
  head -3 "$CH/$F" | cut -c1-90 | sed 's/^/     /'
  # An HMMER3 profile file starts with HMMER3/f. A hit table starts with '#'
  # and target-name columns. One line tells them apart.
  if head -1 "$CH/$F" | grep -q "^HMMER3"; then
    echo "  >> HMMER3 profile database. This is the SEARCH INPUT."
    echo "     Profiles in this chunk: $(grep -c '^HMMER3' "$CH/$F")"
    echo "     Profiles across chunks: $(cat "$CH"/*.hmm 2>/dev/null | grep -c '^HMMER3')"
  else
    echo "  >> not an HMMER3 profile file; identify before assuming anything"
  fi
  echo
  echo "  is it regenerable? looking for the KOfam source it was split from:"
  for c in "$K/profiles" "$K/db/profiles" "$WS2/db/kofam/profiles" \
           "$HOME/kofam/profiles" "$K/ko_list" "$WS2/db/kofam/ko_list"; do
    [ -e "$c" ] && printf '     %-12s %s\n' "$(du -sh "$c" 2>/dev/null | cut -f1)" "${c/#$WS2/\$WS2}"
  done
  echo "     (a chunks/ directory with the source profiles still on disk is"
  echo "      re-splittable and therefore disposable; without them it is the"
  echo "      only copy of the database this project searched against, and"
  echo "      deleting it would make step 45 unreproducible)"
else
  echo "  chunks/ not present"
fi
echo

echo "=============================================================="
echo "B. what is actually inside all_proteins.faa"
echo "=============================================================="
FAA="$K/all_proteins.faa"
if [ -s "$FAA" ]; then
  "$PY" - "$FAA" <<'PYEOF'
import sys, collections
ns = collections.Counter()
ex = {}
n = 0
with open(sys.argv[1], "r", errors="replace") as fh:
    for ln in fh:
        if ln[:1] != ">":
            continue
        n += 1
        i = ln[1:].split()[0]
        # the namespace is whatever precedes the first separator; our own MAG
        # proteins carry bin~contig_gene, imported references carry SET|acc
        if "|" in i:
            k = i.split("|", 1)[0] + "|"
        elif "~" in i:
            k = i.split("~", 1)[0][:3] + "…~  (per-genome)"
        else:
            k = "(no separator)"
        ns[k] += 1
        ex.setdefault(k, i)
print("   %d sequences, %d distinct namespaces" % (n, len(ns)))
for k, c in ns.most_common(14):
    print("     %-28s %8d   e.g. %s" % (k, c, ex[k]))
PYEOF
  echo
  echo "  >> anything in an external namespace is an IMPORTED reference set."
  echo "     It is only rebuildable by re-fetching, which needs the network"
  echo "     and would not return byte-identical sequences later. Treat this"
  echo "     file as primary data, not as a cache."
else
  echo "  not present"
fi
echo

echo "=============================================================="
echo "C. a curl conversion that actually runs"
echo "=============================================================="
for S in "$WS2/scripts/15_marker_trees.sh" "$WS2/15_marker_trees.sh"; do
  [ -f "$S" ] || continue
  NEW="${S%.sh}_curl.sh"
  rm -f "$NEW"
  # Order matters: the stdout forms must be rewritten BEFORE the bare form,
  # or the bare rule eats "wget " and leaves "-qO-" dangling after curl.
  sed -e 's#wget[[:space:]]\{1,\}-qO-#curl -fsSL#g' \
      -e 's#wget[[:space:]]\{1,\}-O-[[:space:]]\{1,\}-q#curl -fsSL#g' \
      -e 's#wget[[:space:]]\{1,\}-O-#curl -fsSL#g' \
      -e 's#wget[[:space:]]\{1,\}-q[[:space:]]\{1,\}-O[[:space:]]\{1,\}\([^[:space:]]\{1,\}\)#curl -fsSL -o \1#g' \
      -e 's#wget[[:space:]]\{1,\}-O[[:space:]]\{1,\}\([^[:space:]]\{1,\}\)#curl -fsSL -o \1#g' \
      -e 's#wget[[:space:]]\{1,\}-q[[:space:]]#curl -fsSLO #g' \
      -e 's#wget[[:space:]]#curl -fsSLO #g' "$S" > "$NEW"
  echo "  ${S#$WS2/}  ->  ${NEW#$WS2/}"
  diff "$S" "$NEW" | sed 's/^/     /'
  echo
  # VERIFY THE OUTPUT. The previous version shipped an unrunnable command
  # because it checked nothing after substituting.
  BAD=$(grep -n "curl" "$NEW" | grep -E '\-qO|O-|--quiet|[[:space:]]-nv|[[:space:]]-nc|--tries' || true)
  if [ -n "$BAD" ]; then
    echo "  [--] leftover wget flags survived into curl commands:"
    printf '%s\n' "$BAD" | sed 's/^/     /'
    echo "       NOT usable. Fix by hand."
  else
    echo "  [OK] no wget flags remain on any curl line"
    echo "  syntax check:"
    if bash -n "$NEW" 2>&1 | sed 's/^/     /'; then echo "     [OK] parses"; fi
  fi
  echo
  echo "  reachability -- against the hosts THIS SCRIPT actually uses:"
  # extracted from the script, not assumed. Step 77 tested UniProt; this
  # script fetches from NCBI E-utilities.
  grep -oE 'https?://[A-Za-z0-9.-]+' "$S" | sort -u | while read -r h; do
    if curl -fsS --max-time 20 -o /dev/null "$h" 2>/dev/null; then
      echo "     [OK] $h"
    else
      echo "     [--] $h  unreachable from this machine"
    fi
  done
  # $EUT is a variable; resolve it too
  E=$(grep -m1 -E '^[[:space:]]*(export[[:space:]]+)?EUT=' "$S" | sed 's/.*EUT=//; s/["'"'"']//g')
  if [ -n "${E:-}" ]; then
    echo "     EUT resolves to: $E"
    curl -fsS --max-time 20 -o /dev/null "$E/einfo.fcgi" 2>/dev/null \
      && echo "     [OK] E-utilities answers" \
      || echo "     [--] E-utilities did not answer -- do not re-run the fetch yet"
  fi
  echo
  echo "  NOTE. NCBI throttles unauthenticated E-utilities to 3 requests/second"
  echo "  and will return 429s in bursts. -f turns those into a non-zero exit"
  echo "  instead of an empty FASTA that the tree step reads as 'no references'."
done
echo

echo "=============================================================="
echo "DISK -- for the record"
echo "=============================================================="
df -h ${SCRATCH_ROOT} | sed 's/^/   /'
echo "   45_kofam_full: $(du -sh "$K" | cut -f1)"
echo
echo "   1.8 TB free. NOTHING here needs deleting. The two large items are:"
echo "     chunks/          the profile database step 45 searched against"
echo "     all_proteins.faa our proteins plus imported references"
echo "   Neither is a duplicate of anything. Keep both."
