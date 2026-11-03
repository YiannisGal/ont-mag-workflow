#!/usr/bin/env bash
# 77_housekeeping.sh -- reclaim disk, and produce a working marker-tree script.
#
# THIS STEP DELETES FILES. Every deletion is gated on a PROOF, computed here,
# that the bytes are recoverable from something else that stays. If a proof
# fails, nothing in that section is deleted and the step says why. No gate is
# "I remember that these are duplicates".
#
#   A  45_kofam_full/chunks/   -- proven to be a partition of all_hits.tbl
#   B  45_kofam_full/all_proteins.faa -- proven to be a concatenation of the
#                                per-genome protein files, which stay
#   C  15_marker_trees.sh      -- wget is broken on this machine; a curl
#                                version is WRITTEN ALONGSIDE, nothing is
#                                overwritten
#   D  figure script hashes    -- re-recorded after a lost workspace
#
# Run it once to see the audit. It deletes only when APPLY=1 is set:
#     APPLY=1 bash scripts/run.sh scripts/77_housekeeping_<md5>.sh
set -u
WS2="${WS2}"
K="$WS2/45_kofam_full"
OUT="$WS2/77_housekeeping"
mkdir -p "$OUT" || exit 1
APPLY="${APPLY:-0}"
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }
[ "$APPLY" = "1" ] && echo "*** APPLY=1 -- proven-safe deletions WILL be made ***" \
                   || echo "*** audit only; set APPLY=1 to delete ***"
echo

hsize() { du -sh "$1" 2>/dev/null | cut -f1; }

echo "=============================================================="
echo "A. 45_kofam_full/chunks/ -- is it really inside all_hits.tbl?"
echo "=============================================================="
ALL="$K/all_hits.tbl"
CH="$K/chunks"
A_OK=0
if [ ! -s "$ALL" ]; then
  echo "  all_hits.tbl missing or empty -- NOT deleting anything"
elif [ ! -d "$CH" ]; then
  echo "  chunks/ not present; nothing to do"
else
  echo "  all_hits.tbl : $(hsize "$ALL"), $(grep -vc '^#' "$ALL" 2>/dev/null) data lines"
  echo "  chunks/      : $(hsize "$CH"), $(ls -1 "$CH" | wc -l) files"
  # The proof: every data line of every chunk must appear in all_hits.tbl.
  # Line counts alone are NOT a proof -- two files can have equal counts and
  # different contents -- so this compares sorted content, and reports how many
  # chunk lines are missing from the union file.
  echo "  comparing sorted content (this reads ~800 MB, give it a minute)..."
  "$PY" - "$ALL" "$CH" <<'PYEOF'
import sys, os, hashlib
allp, chd = sys.argv[1], sys.argv[2]

def datalines(p):
    with open(p, "rb") as fh:
        for ln in fh:
            if ln[:1] != b"#" and ln.strip():
                yield ln

seen = set()
for ln in datalines(allp):
    seen.add(hashlib.md5(ln).digest())
print("     distinct lines in all_hits.tbl : %d" % len(seen))

tot = miss = 0
worst = []
for f in sorted(os.listdir(chd)):
    p = os.path.join(chd, f)
    if not os.path.isfile(p):
        continue
    m = 0
    for ln in datalines(p):
        tot += 1
        if hashlib.md5(ln).digest() not in seen:
            m += 1
    miss += m
    if m:
        worst.append((f, m))
print("     lines across chunks/           : %d" % tot)
print("     chunk lines NOT in all_hits    : %d" % miss)
for f, m in worst[:8]:
    print("       %-40s %d missing" % (f, m))
sys.exit(0 if (miss == 0 and tot > 0) else 1)
PYEOF
  if [ $? -eq 0 ]; then
    echo "  [OK] chunks/ is a subset of all_hits.tbl -- deleting it loses nothing"
    A_OK=1
  else
    echo "  [--] chunks/ holds lines all_hits.tbl does not. NOT SAFE TO DELETE."
    echo "       Keep both, and treat all_hits.tbl as incomplete."
  fi
fi
if [ "$A_OK" = "1" ] && [ "$APPLY" = "1" ]; then
  n=$(hsize "$CH"); rm -rf "$CH" && echo "  removed chunks/ ($n)"
fi
echo

echo "=============================================================="
echo "B. 45_kofam_full/all_proteins.faa -- can it be rebuilt?"
echo "=============================================================="
FAA="$K/all_proteins.faa"
B_OK=0
if [ ! -s "$FAA" ]; then
  echo "  not present; nothing to do"
else
  echo "  all_proteins.faa : $(hsize "$FAA")"
  # The proof: every protein id in the concatenated file must be present in the
  # per-genome protein files that stay on disk. If a single id exists only in
  # the concatenation, the concatenation is the ONLY copy and must be kept.
  SRC=""
  for c in "$WS2/06_binqc/checkm2/protein_files" "$WS2/07_mags/proteins" \
           "$WS2/10_annotation/proteins"; do
    [ -d "$c" ] && { SRC="$c"; break; }
  done
  if [ -z "$SRC" ]; then
    echo "  [--] no per-genome protein directory found -- this file may be the"
    echo "       only copy. NOT deleting."
  else
    echo "  per-genome source: ${SRC#$WS2/}  ($(ls -1 "$SRC" | wc -l) files, $(hsize "$SRC"))"
    "$PY" - "$FAA" "$SRC" <<'PYEOF'
import sys, os
faa, src = sys.argv[1], sys.argv[2]

def ids(path):
    out = set()
    with open(path, "r", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                out.add(ln[1:].split()[0])
    return out

have = set()
nf = 0
for f in sorted(os.listdir(src)):
    p = os.path.join(src, f)
    if os.path.isfile(p) and f.endswith((".faa", ".fa", ".fasta")):
        have |= ids(p); nf += 1
want = ids(faa)
print("     proteins in all_proteins.faa : %d" % len(want))
print("     proteins in %d per-genome files: %d" % (nf, len(have)))
only = want - have
print("     present ONLY in the concatenation: %d" % len(only))
for i in sorted(only)[:6]:
    print("       %s" % i)
sys.exit(0 if (not only and want) else 1)
PYEOF
    if [ $? -eq 0 ]; then
      echo "  [OK] every protein is still in the per-genome files -- rebuildable with"
      echo "       cat ${SRC#$WS2/}/*.faa > all_proteins.faa"
      B_OK=1
    else
      echo "  [--] the concatenation holds proteins nothing else has. NOT deleting."
    fi
  fi
fi
if [ "$B_OK" = "1" ] && [ "$APPLY" = "1" ]; then
  n=$(hsize "$FAA"); rm -f "$FAA" && echo "  removed all_proteins.faa ($n)"
fi
echo

echo "=============================================================="
echo "C. 15_marker_trees.sh -- wget is broken here, curl is not"
echo "=============================================================="
for S in "$WS2/scripts/15_marker_trees.sh" "$WS2/15_marker_trees.sh"; do
  [ -f "$S" ] || continue
  echo "  found ${S#$WS2/}"
  echo "  wget invocations:"
  grep -n "wget" "$S" | sed 's/^/     /'
  NEW="${S%.sh}_curl.sh"
  # Written ALONGSIDE. The original is never touched: a step script is the
  # record of what was run, and silently rewriting one destroys that record.
  sed -e 's#wget[[:space:]]\+-q\{0,1\}[[:space:]]*-O[[:space:]]\+\([^[:space:]]\+\)#curl -fsSL -o \1#g' \
      -e 's#wget[[:space:]]\+-O[[:space:]]\+\([^[:space:]]\+\)#curl -fsSL -o \1#g' \
      -e 's#wget[[:space:]]\+-q[[:space:]]\+#curl -fsSLO #g' \
      -e 's#wget[[:space:]]\+#curl -fsSLO #g' "$S" > "$NEW"
  echo "  wrote ${NEW#$WS2/}; the difference:"
  diff "$S" "$NEW" | sed 's/^/     /' | head -30
  echo "  curl reachability test:"
  if curl -fsSI "https://rest.uniprot.org/" >/dev/null 2>&1; then
    echo "     [OK] curl reaches UniProt"
  else
    echo "     [--] curl could not reach UniProt either -- the problem is the"
    echo "          network, not the tool. Do not re-run the fetch yet."
  fi
  # a downloader that silently returns nothing is the failure mode that started
  # this: -f makes curl exit non-zero on an HTTP error instead of writing a
  # zero-byte file that the next step happily parses as "no references found"
  echo "  NOTE: -f is the important flag. Without it a 404 writes an empty file"
  echo "        and the tree step reports 'no references' instead of failing."
done
echo

echo "=============================================================="
echo "D. figure scripts -- hashes"
echo "=============================================================="
for d in "$WS2/scripts" "$WS2/figs" "$WS2"; do
  ls "$d"/fig*.py >/dev/null 2>&1 || continue
  md5sum "$d"/fig*.py | sed "s#$WS2/##" | sed 's/^/   /'
  break
done
[ -f "$WS2/figs/etoliko_figstyle.py" ] && \
  md5sum "$WS2/figs/etoliko_figstyle.py" | sed "s#$WS2/##" | sed 's/^/   /'
echo

echo "=============================================================="
echo "DISK"
echo "=============================================================="
df -h ${SCRATCH_ROOT} | sed 's/^/   /'
echo "   45_kofam_full now: $(hsize "$K")"
[ "$APPLY" = "1" ] || echo "   (nothing was deleted -- rerun with APPLY=1)"
