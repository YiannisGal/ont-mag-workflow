#!/usr/bin/env bash
# 80_mtbB_mtaB_finish.sh -- can the decoy labels be quoted as gene counts?
#
# WHAT STEP 79 ESTABLISHED. Its control reproduced the published mtmB result
# exactly (18 / 24 / 0.75) and every decoy count matches Volume XV verbatim, so
# the best-hit table is being read correctly. The two numbers wanted are on it:
#
#     DECOYmtbB   283 sediment / 296 cultures
#     DECOYmtaB    17 sediment / 194 cultures
#
# WHAT STEP 79 GOT WRONG, IN THE SCRIPT AFTER THE ONE THAT MADE THE RULE.
# Section C ran `ls *.faa | head -1`, which is ctrl.faa -- the 13-sequence
# archaeal control -- and so reported the reference-set sizes of the control
# instead of the panel. 61_mtmB holds FOUR fasta files (ctrl, own, panel,
# refs). That is rule 125 again: a directory is not one kind of file, and
# taking the first name in it is not sampling. Section A below reads all four.
#
# THE QUESTION THIS STEP ACTUALLY DECIDES. Best-hit assignment is symmetric --
# a read goes to the highest-scoring reference whether that reference was
# labelled target or decoy -- so DECOYmtbB counts ARE mtbB counts *provided the
# mtbB reference set was built like the mtmB one*. If the decoy sets were
# assembled from broad superfamily queries, they are deliberately greedy and
# their catch is not a gene count. Section A is what tells the two apart, and
# 283 sediment reads against 24 mcrA -- 11.8x, worse than the 4.9x that made
# the old number unreportable -- is the reason to ask before quoting.
#
# THE CEILING IS DERIVED, THEN VALIDATED ON A KNOWN ANSWER. Expected reads for
# a gene scale as (copies in the community) x (gene length + read length),
# because a read hits a gene by overlapping it. Section D computes that from
# the panel's own sequence lengths and the library's own read length, and
# checks it reproduces mtmB's published 0.78 BEFORE applying it to anything
# else. A formula that cannot recover the case it was validated on does not
# get to adjudicate the cases it was not.
set -u
WS2="${WS2}"
P61="$WS2/61_mtmB"
P63="$WS2/63_mtmB_relabel"
OUT="$WS2/80_mtbB_finish"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. the panel, all four fasta files, per label"
echo "=============================================================="
"$PY" - "$P61" "$P63" <<'PYEOF'
import sys, os, collections, statistics
P61, P63 = sys.argv[1], sys.argv[2]
bp = os.path.join(P63, "bacterial.txt")
bact = {l.strip() for l in open(bp)} if os.path.exists(bp) else set()

def lab(i):
    return "DECOYbact" if i.split("|")[-1] in bact else i.split("|")[0]

allp = collections.defaultdict(list)
for f in sorted(os.listdir(P61)):
    if not f.endswith((".faa", ".fa", ".fasta")):
        continue
    per = collections.defaultdict(list)
    cur, n = None, 0
    for ln in open(os.path.join(P61, f), errors="replace"):
        if ln.startswith(">"):
            if cur:
                per[cur].append(n); allp[cur].append(n)
            cur = lab(ln[1:].split()[0]); n = 0
        else:
            n += len(ln.strip())
    if cur:
        per[cur].append(n); allp[cur].append(n)
    print("   %s" % f)
    for k in sorted(per):
        v = per[k]
        print("      %-14s %4d seq   median %4d aa   range %d-%d"
              % (k, len(v), int(statistics.median(v)), min(v), max(v)))
print()
print("   COMBINED, every label in the panel:")
print("   %-14s %6s %10s %12s" % ("label", "refs", "median aa", "range"))
for k in sorted(allp):
    v = allp[k]
    print("   %-14s %6d %10d %12s"
          % (k, len(v), int(statistics.median(v)), "%d-%d" % (min(v), max(v))))
import json
json.dump({k: v for k, v in allp.items()}, open(os.path.join(os.path.dirname(P61.rstrip("/")), "80_mtbB_finish", "panel_lengths.json"), "w"))
print()
print("   Read this before believing any decoy count: a label whose references")
print("   are few and long behaves differently from one that is many and short,")
print("   and a superfamily label is greedy by construction.")
PYEOF
echo
echo "   what the panel says about its own provenance:"
[ -s "$P61/panel_kos.tsv" ] && head -20 "$P61/panel_kos.tsv" | sed 's/^/     /'
echo

echo "=============================================================="
echo "B. how step 61 itself computed the 0.78 prediction"
echo "=============================================================="
# The arithmetic exists in the run. Reconstructing it from memory is how the
# wrong ceiling got asserted in the first place (rule 110).
S=$(ls "$WS2"/scripts/61_mtmB*.sh "$WS2"/61_mtmB*.sh 2>/dev/null | head -1)
if [ -n "${S:-}" ]; then
  echo "   ${S#$WS2/}"
  grep -n -i -E "0\.78|predict|span|ceiling|1864|2914|gene_len|readlen" "$S" \
    | head -30 | sed 's/^/     /'
else
  echo "   !! step 61's script not found; section D derives the formula instead"
fi
echo

echo "=============================================================="
echo "C. catalogue gene counts ABOVE the curated threshold"
echo "=============================================================="
# Step 79 printed any-hit counts. Any-hit and above-threshold differ: mcrA came
# back as 13 any-hit where the published above-threshold count is 9. The
# threshold table is the arbiter (rule 97).
KO="$WS2/45_kofam_full/all_hits.tbl"
THR=""
for c in "$WS2/36_gapclose/ko_list" "$WS2/45_kofam_full/ko_list" \
         "$WS2/db/kofam/ko_list"; do
  [ -s "$c" ] && { THR="$c"; break; }
done
echo "   hits    : ${KO#$WS2/}"
echo "   threshold table: ${THR:-NOT FOUND}"
if [ -n "$THR" ] && [ -s "$KO" ]; then
  "$PY" - "$KO" "$THR" "$OUT/ko_counts.tsv" <<'PYEOF'
import sys
hits, thr = sys.argv[1], sys.argv[2]
WANT = {"K00399": "mcrA", "K16176": "mtmB", "K16178": "mtbB", "K04480": "mtaB"}
cut = {}
with open(thr, errors="replace") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    try:
        ik, it = head.index("knum"), head.index("threshold")
    except ValueError:
        ik, it = 0, 1
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) > max(ik, it) and f[ik] in WANT:
            try:
                cut[f[ik]] = float(f[it])
            except ValueError:
                cut[f[ik]] = None
print("   curated thresholds: %s"
      % ", ".join("%s %s" % (k, cut.get(k)) for k in sorted(WANT)))

# hmmsearch --tblout: field 0 target, field 2 query(KO), field 5 full-seq score.
# The description field carries prodigal coordinates and is NOT a score
# (rule 114) -- that mistake produced an identical 96142.0 across eight
# profiles in step 71.
above, any_ = {k: set() for k in WANT}, {k: set() for k in WANT}
for ln in open(hits, errors="replace"):
    if ln.startswith("#"):
        continue
    f = ln.split()
    if len(f) < 6 or f[2] not in WANT:
        continue
    any_[f[2]].add(f[0])
    t = cut.get(f[2])
    try:
        sc = float(f[5])
    except ValueError:
        continue
    if t is None or sc >= t:
        above[f[2]].add(f[0])
print()
print("   %-8s %-6s %10s %16s" % ("KO", "gene", "any hit", "above threshold"))
for k in sorted(WANT):
    print("   %-8s %-6s %10d %16d" % (k, WANT[k], len(any_[k]), len(above[k])))
print()
print("   CONTROL: mtmB must be 11 and mcrA 9 above threshold. If they are not,")
print("   this is not the count the published prediction was built from and the")
print("   ceiling below is not comparable to it.")
with open(sys.argv[3], "w") as fh:
    for k in sorted(WANT):
        fh.write("%s\t%s\t%d\n" % (k, WANT[k], len(above[k])))
PYEOF
else
  echo "   !! cannot count above threshold without both files."
  echo "      Using any-hit counts here would repeat rule 97: the curated"
  echo "      threshold is the point of KOfam."
fi
echo

echo "=============================================================="
echo "D. the ceiling, derived and validated on mtmB"
echo "=============================================================="
"$PY" - "$P61" "$OUT" "$WS2" <<'PYEOF'
import sys, os, json, statistics
P61, OUT, WS2 = sys.argv[1], sys.argv[2], sys.argv[3]

lens = {}
p = os.path.join(OUT, "panel_lengths.json")
if os.path.exists(p):
    lens = {k: statistics.median(v) * 3.0 for k, v in json.load(open(p)).items()}

ko = {}
kp = os.path.join(OUT, "ko_counts.tsv")
if os.path.exists(kp):
    for ln in open(kp):
        k, g, n = ln.rstrip("\n").split("\t")
        ko[g] = int(n)
if not ko or not lens:
    sys.exit("   sections A or C did not produce their inputs; no ceiling computed")

counts = {"mtmB": (18, 112), "mcrA": (24, 355),
          "mtbB": (283, 296), "mtaB": (17, 194)}
LBL = {"mtmB": "mtmB", "mcrA": "mcrA", "mtbB": "DECOYmtbB", "mtaB": "DECOYmtaB"}

# mean read length of the library the sediment counts came from
R = None
for c in ("00_input/sediment_reads_NanoStats.txt", "00_input/SED_NanoStats.txt"):
    f = os.path.join(WS2, c)
    if os.path.exists(f):
        for ln in open(f, errors="replace"):
            if "mean read length" in ln.lower():
                try:
                    R = float(ln.split(":")[1].strip().replace(",", ""))
                except (IndexError, ValueError):
                    pass
        break
if R is None:
    R = 8.44e9 / 7528854            # bases / reads, from the settled library size
    print("   read length not found in NanoStats; using bases/reads = %.0f bp" % R)
else:
    print("   mean read length from NanoStats: %.0f bp" % R)

def span(g):
    L = lens.get(LBL[g])
    return None if L is None else L + R

print()
print("   %-6s %6s %10s %10s" % ("gene", "genes", "median bp", "span bp"))
for g in ("mcrA", "mtmB", "mtbB", "mtaB"):
    s = span(g)
    print("   %-6s %6s %10s %10s"
          % (g, ko.get(g, "?"),
             "%.0f" % lens[LBL[g]] if LBL[g] in lens else "?",
             "%.0f" % s if s else "?"))

sm, nm = span("mcrA"), ko.get("mcrA")
if not sm or not nm:
    sys.exit("   mcrA span or gene count unavailable -- no prediction possible")

print()
print("   predicted read ratio = (genes_x / genes_mcrA) x (span_x / span_mcrA)")
print()
print("   %-6s %10s %10s %8s %s" % ("gene", "predicted", "measured", "excess", "verdict"))
val = None
for g in ("mtmB", "mtbB", "mtaB"):
    s, n = span(g), ko.get(g)
    if not s or not n:
        print("   %-6s  inputs missing" % g); continue
    pred = (n / float(nm)) * (s / sm)
    meas = counts[g][0] / float(counts["mcrA"][0])
    exc = meas / pred if pred else float("inf")
    if g == "mtmB":
        val = abs(pred - 0.78) <= 0.08
        v = "CONTROL: %s published 0.78" % ("reproduced" if val else "NOT reproduced")
    else:
        v = ("within the ceiling -- reportable" if exc <= 1.5 else
             "%.1fx over the ceiling -- NOT reportable" % exc)
    print("   %-6s %10.2f %10.2f %8.1fx %s" % (g, pred, meas, exc, v))

print()
if val:
    print("   The formula recovers the published mtmB prediction, so the two")
    print("   verdicts above are made on the same basis as the number already in")
    print("   the report. A label that still exceeds its ceiling is measuring a")
    print("   family, not a gene, and stays out of the chapter.")
else:
    print("   THE FORMULA DOES NOT RECOVER mtmB's PUBLISHED 0.78. Nothing above")
    print("   is adjudicated. Compare sections A, C and D to find the")
    print("   difference -- the constants are in step 61.")
PYEOF
