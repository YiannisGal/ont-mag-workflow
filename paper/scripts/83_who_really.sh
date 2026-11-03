#!/usr/bin/env bash
# 83_who_really.sh -- resolve the organism behind every reference, properly.
#
# STEP 82 FAILED, AND ITS OWN CONTROL CAUGHT IT. It read the organism out of
# the fasta header, and step 61's panel has no organism in any header -- the
# ids are bare accessions (TQS84672.1) and bare contig names. So every
# reference came back "NOT a methanogen", INCLUDING mcrA's, whose references
# are our own mcrA-bearing contigs and are methanogens by construction.
#
#     mcrA -- 24 reads ... methanogen references hold 0 reads; others hold 24
#
# A test that returns the same verdict for every input, control included, has
# measured nothing. The "methanogen-only 0.00" column in step 82 is an artefact
# and must not be quoted.
#
# THE TWO KINDS OF REFERENCE NEED TWO DIFFERENT LOOKUPS, which is why one
# header rule could never have worked:
#
#   OURS      LIQ_contig_522_668 -- a gene on one of our own contigs. Its
#             organism is whichever MAG that contig belongs to, and that is on
#             disk in 07_mags/bins/. No network.
#
#   IMPORTED  TQS84672.1 -- an NCBI accession. Its organism is at NCBI, and
#             E-utilities is reachable from this machine (step 78 confirmed it
#             against the host this code actually uses).
#
# CONTROLS, BOTH DIRECTIONS. mcrA's references must resolve to archaea -- if
# they do not, the contig-to-MAG join is broken and nothing else is
# trustworthy. And the AAM* accessions should return Methanosarcina
# acetivorans, which Volume XV already identified by hand; if NCBI says
# otherwise, the accession lookup is wrong.
#
# Anything that resolves to NOTHING is reported as unknown and counted into
# neither column. "Absent" and "not tested" do not share a column (rule 113).
set -u
WS2="${WS2}"
P82="$WS2/82_mtbB_who/reference_taxonomy.tsv"
OUT="$WS2/83_who_really"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
[ -s "$P82" ] || { echo "   need $P82 -- run step 82 first"; exit 1; }

echo "=============================================================="
echo "A. split the references by where their taxonomy lives"
echo "=============================================================="
"$PY" - "$P82" "$OUT" <<'PYEOF'
import sys, csv, re, os
src, OUT = sys.argv[1], sys.argv[2]
own, ext = set(), set()
for r in csv.DictReader(open(src), delimiter="\t"):
    ref = r["reference"].split("|")[-1]
    (own if re.match(r"^(LIQ|SED)_contig_", ref) else ext).add(ref)
open(os.path.join(OUT, "own.txt"), "w").write("\n".join(sorted(own)) + "\n")
open(os.path.join(OUT, "ext.txt"), "w").write("\n".join(sorted(ext)) + "\n")
print("   %d reference(s) are our own contigs -> 07_mags/bins" % len(own))
print("   %d reference(s) are NCBI accessions -> E-utilities" % len(ext))
PYEOF
echo

echo "=============================================================="
echo "B. our own contigs -> which MAG -> which domain"
echo "=============================================================="
BINS="$WS2/07_mags/bins"
: > "$OUT/own_tax.tsv"
if [ -d "$BINS" ]; then
  grep -H "^>" "$BINS"/*.fa 2>/dev/null \
    | sed 's#^.*/##; s/\.fa:>/\t/' | awk -F'\t' '{print $2"\t"$1}' \
    | sort -u > "$OUT/contig2bin.tsv"
  echo "   $(wc -l < "$OUT/contig2bin.tsv") contig->bin entries"
  "$PY" - "$OUT" "$WS2" <<'PYEOF'
import sys, os, csv, re
OUT, WS2 = sys.argv[1], sys.argv[2]
c2b = {}
for ln in open(os.path.join(OUT, "contig2bin.tsv"), errors="replace"):
    f = ln.rstrip("\n").split("\t")
    if len(f) == 2:
        c2b[f[0].split()[0]] = f[1]

dom = {}
m = os.path.join(WS2, "47_master", "MASTER.tsv")
if os.path.exists(m):
    for r in csv.DictReader(open(m, errors="replace"), delimiter="\t"):
        for k in ("working_name", "mag_id", "bin", "qc_bin_name"):
            v = (r.get(k) or "").strip()
            if v:
                dom[v] = ((r.get("domain") or "?").strip(),
                          (r.get("genus") or "?").strip())
else:
    print("   !! MASTER.tsv not found -- domain cannot be assigned")

rows, miss = [], 0
for ref in open(os.path.join(OUT, "own.txt")):
    ref = ref.strip()
    if not ref:
        continue
    mm = re.match(r"^((?:LIQ|SED)_contig_\d+)_\d+$", ref)
    contig = mm.group(1) if mm else ref
    b = c2b.get(contig)
    d, g = dom.get(b, ("?", "?")) if b else ("?", "?")
    if not b:
        miss += 1
    rows.append((ref, contig, b or "(unbinned)", d, g))
with open(os.path.join(OUT, "own_tax.tsv"), "w") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\n")
print("   %d resolved to a bin, %d unbinned" % (len(rows) - miss, miss))
for r in rows[:10]:
    print("     %-22s %-22s %-10s %s" % (r[0][:22], r[2][:22], r[3], r[4]))
if len(rows) > 10:
    print("     ... %d more in own_tax.tsv" % (len(rows) - 10))
PYEOF
else
  echo "   !! $BINS not found"
fi
echo

echo "=============================================================="
echo "C. NCBI accessions -> organism"
echo "=============================================================="
: > "$OUT/ext_tax.tsv"
rm -f "$OUT/esummary.json.part"
N=$(grep -c . "$OUT/ext.txt" 2>/dev/null || echo 0)
echo "   $N accession(s) to resolve"
if [ "$N" -gt 0 ]; then
  split -l 150 "$OUT/ext.txt" "$OUT/batch_"
  for b in "$OUT"/batch_*; do
    ids=$(paste -sd, "$b")
    curl -fsSL --max-time 60 \
      "$EUT/esummary.fcgi?db=protein&retmode=json&id=$ids" \
      >> "$OUT/esummary.json.part" 2>/dev/null \
      || echo "   !! request failed for $(basename "$b")"
    echo >> "$OUT/esummary.json.part"
    sleep 0.4
  done
  rm -f "$OUT"/batch_*
  "$PY" - "$OUT" <<'PYEOF'
import sys, os, json, re
OUT = sys.argv[1]
got = {}
p = os.path.join(OUT, "esummary.json.part")
if os.path.exists(p):
    txt = open(p, errors="replace").read()
    for t in re.findall(r"\{.*?\n\}", txt, re.S) or [txt]:
        try:
            d = json.loads(t)
        except ValueError:
            continue
        res = d.get("result", {})
        for k in res.get("uids", []):
            r = res.get(k, {})
            org = r.get("organism") or ""
            if not org:
                m = re.search(r"\[([^\]]+)\]\s*$", r.get("title", ""))
                org = m.group(1) if m else ""
            acc = r.get("accessionversion") or r.get("caption") or k
            if org:
                got[acc] = org
                got[acc.split(".")[0]] = org
want = [l.strip() for l in open(os.path.join(OUT, "ext.txt")) if l.strip()]
with open(os.path.join(OUT, "ext_tax.tsv"), "w") as fh:
    for a in want:
        fh.write("%s\t%s\n" % (a, got.get(a) or got.get(a.split(".")[0], "")))
res = sum(1 for a in want if got.get(a) or got.get(a.split(".")[0]))
print("   resolved %d of %d (%.0f%%)" % (res, len(want), 100.0 * res / len(want)))
if res < 0.9 * len(want):
    print("   !! under 90 percent resolved. Treat the split below as provisional.")
print()
print("   CONTROL -- Volume XV identified the AAM* entries by hand as")
print("   Methanosarcina acetivorans. NCBI says:")
for a in want:
    if a.startswith("AAM"):
        print("     %-16s %s" % (a, got.get(a) or got.get(a.split(".")[0], "(unresolved)")))
print()
print("   the references holding the most reads:")
for a in ("TQS84672.1", "WWV68927.1", "ATW28134.1", "RZN62963.1",
          "TFG55810.1", "AAM03976.1", "AAM05811.1", "AAM04366.1"):
    print("     %-16s %s" % (a, got.get(a) or got.get(a.split(".")[0], "(unresolved)")))
PYEOF
fi
echo

echo "=============================================================="
echo "D. recount each label by carrier type"
echo "=============================================================="
"$PY" - "$P82" "$OUT" <<'PYEOF'
import sys, os, csv, re, collections
src, OUT = sys.argv[1], sys.argv[2]

own = {}
p = os.path.join(OUT, "own_tax.tsv")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 5:
            own[f[0]] = (f[3], f[4])

ext = {}
p = os.path.join(OUT, "ext_tax.tsv")
if os.path.exists(p):
    for ln in open(p, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) == 2:
            ext[f[0]] = f[1]

MG = re.compile(r"^(Methano|Methermicoccus|Candidatus Methano)", re.I)


def carrier(ref):
    r = ref.split("|")[-1]
    if r in own:
        d, g = own[r]
        if d.lower().startswith("arch"):
            return "archaeon"
        if d.lower().startswith("bact"):
            return "bacterium"
        return "unknown"
    o = ext.get(r, "")
    if not o:
        return "unknown"
    return "archaeon" if MG.match(o) else "bacterium"


agg = collections.defaultdict(collections.Counter)
for r in csv.DictReader(open(src), delimiter="\t"):
    agg[(r["sample"], r["label"])][carrier(r["reference"])] += int(r["reads"])

for tag in ("SED", "LIQ"):
    mcra = sum(agg[(tag, "mcrA")].values())
    print("   %s   mcrA = %d reads" % (tag, mcra))
    print("      %-13s %7s %9s %10s %8s %10s %14s"
          % ("label", "total", "archaeal", "bacterial", "unknown",
             "all/mcrA", "archaeal/mcrA"))
    for lab in ("mtmB", "DECOYmtaB", "DECOYmtbB", "mcrA"):
        c = agg.get((tag, lab))
        if not c:
            continue
        t = sum(c.values())
        a, b, u = c["archaeon"], c["bacterium"], c["unknown"]
        print("      %-13s %7d %9d %10d %8d %10s %14s"
              % (lab, t, a, b, u,
                 "%.2f" % (t / float(mcra)) if mcra else "--",
                 "%.2f" % (a / float(mcra)) if mcra else "--"))
    print()

print("   CONTROL: mcrA's own row must be 100 percent archaeal. If it is not,")
print("   the contig-to-MAG join in section B is broken and every other row")
print("   here is void -- that is exactly how step 82 produced its zeros.")
print()
print("   The ceiling 'mtbB cannot exceed mcrA' constrains ARCHAEAL carriers")
print("   only. Read the archaeal/mcrA column against it, and treat the")
print("   bacterial column as a separate observation about MttB-superfamily")
print("   genes in this sediment -- which this project has already met once,")
print("   in the two Desulfosporosinus meridiei proteins step 63 removed.")
PYEOF
