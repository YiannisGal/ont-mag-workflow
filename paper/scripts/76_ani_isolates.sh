#!/usr/bin/env bash
# 76_ani_isolates.sh -- did the ANI run ever contain a cultured isolate?
#
# READ-ONLY except its own output directory.
#
# WHY THIS STEP EXISTS. Step 75 found the table: 40_integration/skani_raw.tsv,
# and 40_integration/INTEGRATION.txt quotes the control at 99.75 %. But every
# path visible in that table so far points into 07_mags/bins -- the 29 MAGs.
# The chapter's claim is "no MAG matches a CULTURED ISOLATE at >=95 % ANI", and
# those are not the same statement:
#
#     a MAG scoring below 95 % against isolates that WERE in the comparison
#     is a result;
#     a MAG scoring below 95 % against isolates that were NEVER in the
#     comparison is not a result, it is an empty set.
#
# That is the same trap as before -- "absent" and "not tested" printed in one
# column. So this step does not compute anything new. It asks one question:
# WHICH GENOMES WERE ACTUALLY COMPARED. If the answer is "the 29 MAGs and
# nothing else", the figure cannot be drawn and the sentence in the chapter has
# to change, not be illustrated.
set -u
WS2="${WS2}"
OUT="$WS2/76_ani_isolates"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }
BINS="$WS2/07_mags/bins"

echo "=============================================================="
echo "A. what step 40 says it did, in its own words"
echo "=============================================================="
if [ -s "$WS2/40_integration/INTEGRATION.txt" ]; then
  # print the ANI section only, up to the next top-level heading
  awk '/^[[:space:]]*A\./ {p=1} p && /^[[:space:]]*B\./ {exit} p' \
      "$WS2/40_integration/INTEGRATION.txt" | head -60
else
  echo "  INTEGRATION.txt missing"
fi
echo

echo "=============================================================="
echo "B. the genome sets actually present in the ANI table"
echo "=============================================================="
"$PY" - "$WS2/40_integration/skani_raw.tsv" "$BINS" <<'PYEOF'
import sys, os, collections
src, bins = sys.argv[1], sys.argv[2]
if not os.path.exists(src):
    sys.exit("  skani_raw.tsv missing")

mags = set()
if os.path.isdir(bins):
    mags = {f for f in os.listdir(bins) if f.endswith((".fa", ".fna", ".fasta"))}
print("  genomes in 07_mags/bins            : %d" % len(mags))

dirs = collections.Counter()
seen = set()
rows = 0
with open(src, encoding="utf-8", errors="replace") as fh:
    head = fh.readline()
    # skani's own header starts with Ref_file; a headerless file must not lose
    # its first row, so put it back if it is not a header
    lines = [] if head.startswith("Ref_file") else [head]
    lines += fh.readlines()
for ln in lines:
    f = ln.rstrip("\n").split("\t")
    if len(f) < 3:
        continue
    rows += 1
    for p in (f[0], f[1]):
        seen.add(os.path.basename(p))
        dirs[os.path.dirname(p)] += 1

print("  pairs in the table                 : %d" % rows)
print("  distinct genome files in the table : %d" % len(seen))
print()
print("  every directory the table draws from:")
for d, n in dirs.most_common():
    print("    %6d  %s" % (n, d))
print()

outside = sorted(g for g in seen if g not in mags)
print("  GENOMES IN THE TABLE THAT ARE NOT ONE OF THE 29 MAGs: %d" % len(outside))
if outside:
    for g in outside[:40]:
        print("    %s" % g)
    if len(outside) > 40:
        print("    ... and %d more" % (len(outside) - 40))
else:
    print("    none. The run compared the MAGs against each other and against")
    print("    nothing else, so it cannot answer a question about isolates.")
PYEOF
echo

echo "=============================================================="
echo "C. how step 40 built its genome list"
echo "=============================================================="
# The list file is the evidence. Whatever the prose says, the ANI tool only saw
# what was written into --ql / --rl.
grep -n -A4 -B4 -E -- "--ql|--rl|-q |--query|genome_list|iso.*list|list.*iso" \
     "$WS2/scripts/40_integration.sh" 2>/dev/null | head -50
echo
for f in "$WS2"/40_integration/*list*.txt "$WS2"/40_integration/*.txt; do
  case "$f" in *INTEGRATION.txt) continue ;; esac
  [ -s "$f" ] || continue
  echo "--- ${f#$WS2/}  ($(wc -l < "$f") lines)"
  head -4 "$f" | sed 's/^/     /'
done
echo

echo "=============================================================="
echo "D. do isolate genome assemblies exist at all, anywhere reachable?"
echo "=============================================================="
# If they exist the comparison is a cheap job and the figure becomes drawable.
# If they do not, no amount of re-running produces the number.
for R in "$HOME/Documents/UniqueELT_Solid" \
         "$HOME/Documents/Metagenome_analysis_PhD" ; do
  [ -d "$R" ] || { echo "  (no such tree: $R)"; continue; }
  echo "--- under ${R/#$HOME/~}"
  find "$R" -maxdepth 6 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) \
       -size +500k ! -path "*/07_mags/*" ! -path "*/06_binqc/*" \
       ! -path "*/.pixi/*" ! -path "*/envs/*" 2>/dev/null \
    | head -40 | while read -r g; do
        printf '   %8s  %s\n' "$(du -h "$g" | cut -f1)" "${g/#$HOME/~}"
      done
  n=$(find "$R" -maxdepth 6 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) \
       -size +500k ! -path "*/07_mags/*" ! -path "*/06_binqc/*" \
       ! -path "*/.pixi/*" ! -path "*/envs/*" 2>/dev/null | wc -l)
  echo "   total candidate genome-sized FASTAs: $n"
done
echo

echo "=============================================================="
echo "E. the two controls step 40 set for itself"
echo "=============================================================="
# POSITIVE: the Chakrabartyella pair recovered from both sediment and culture
# must come back near 99.8 %. NEGATIVE: an archaeal MAG against a bacterial
# genome must produce NO ANI value at all -- skani reports nothing below its
# screening threshold, and a number there would mean the run is not measuring
# what it claims.
"$PY" - "$WS2/40_integration/skani_raw.tsv" <<'PYEOF'
import sys, os
src = sys.argv[1]
if not os.path.exists(src):
    sys.exit("  table missing")
best = []
with open(src, encoding="utf-8", errors="replace") as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3 or f[0] == "Ref_file":
            continue
        try: ani = float(f[2])
        except ValueError: continue
        a, b = os.path.basename(f[0]), os.path.basename(f[1])
        if a == b:
            continue                      # self-hit, not evidence of anything
        best.append((ani, a, b))
best.sort(reverse=True)
print("  highest non-self ANI pairs in the table:")
for ani, a, b in best[:8]:
    print("    %6.2f  %s  x  %s" % (ani, a, b))
print()
n95 = sum(1 for ani, a, b in best if ani >= 95.0)
print("  non-self pairs at >=95 %% ANI: %d" % n95)
print("  (each such pair is two genome files that are the same species)")
PYEOF
echo

echo "=============================================================="
echo "VERDICT TO READ OFF THIS"
echo "=============================================================="
echo "B says 'none'  + D finds no isolate FASTAs"
echo "     -> the comparison was never possible. The chapter cannot claim a"
echo "        measured 0/29 against isolates; it can only say the isolate"
echo "        genomes were not sequenced. No figure."
echo "B says 'none'  + D finds isolate FASTAs"
echo "     -> the comparison is runnable and cheap. I write step 77, it runs,"
echo "        and fig11 is drawn from its output."
echo "B lists isolate genomes"
echo "     -> the number already exists; fig11 is drawn from B and E"
echo "        directly."
