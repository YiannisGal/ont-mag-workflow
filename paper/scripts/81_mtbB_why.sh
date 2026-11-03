#!/usr/bin/env bash
# 81_mtbB_why.sh -- why mtbB fails where mtaB passes.
#
# STEP 80 SETTLED THE VERDICTS. The ceiling formula reproduced the published
# mtmB prediction (0.76 against 0.78, measured 0.75), so it is the same basis
# as the number already in the report:
#
#     gene   genes  predicted  measured  excess
#     mtmB      11       0.76      0.75    1.0x   control, reproduced
#     mtaB       9       0.89      0.71    0.8x   within the ceiling
#     mtbB      11       0.77     11.79   15.4x   over, by fifteenfold
#
# THE ANOMALY WORTH EXPLAINING. mtbB and mtmB have the SAME catalogue
# abundance -- 11 proteins each above their curated thresholds -- and almost
# the same predicted ratio, 0.77 against 0.76. Yet mtbB collects 283 sediment
# reads where mtmB collects 18. Fifteen times the reads for the same number of
# genes is not biology; it is the reference set. A thesis sentence saying "mtbB
# could not be measured" is much stronger if it says what went wrong.
#
# TWO TESTS, BOTH ALREADY PAID FOR. Every alignment is in step 61's hit table,
# so neither needs a new search:
#
#   A. PROMISCUITY, per reference, at 4x the expected share 1/n (rule 108 --
#      a flat 25% threshold is a different test at every set size, and it is
#      what missed mtmB's promiscuous entries originally). Step 79 never ran
#      this on the decoys because it read the wrong fasta file.
#
#   B. RUNNER-UP, the test step 62 used to separate genuine from spurious.
#      For every read won by a label, what is the best hit OUTSIDE that label
#      and how far behind is it? Sequence only that family can see has no
#      runner-up at all. Reads whose runner-up is a hair behind were won by
#      an accident of scoring and belong to no gene in particular.
set -u
WS2="${WS2}"
P61="$WS2/61_mtmB"
P63="$WS2/63_mtmB_relabel"
OUT="$WS2/81_mtbB_why"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$P61" "$P63" "$OUT" <<'PYEOF'
import sys, os, collections, statistics
P61, P63, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MINLEN = 30
FOCUS = ("mtmB", "DECOYmtaB", "DECOYmtbB", "DECOYmttB", "DECOYmttBsf", "mcrA")

bp = os.path.join(P63, "bacterial.txt")
bact = {l.strip() for l in open(bp)} if os.path.exists(bp) else set()
if not bact:
    sys.exit("   bacterial.txt missing -- refusing (step 63 removed 7 imported "
             "bacterial proteins holding 67 sediment reads)")


def lab(i):
    return "DECOYbact" if i.split("|")[-1] in bact else i.split("|")[0]


# --- unique references per label. The four fasta files overlap: panel.faa is
# largely own.faa + refs.faa, so adding their counts double-counts and makes
# 1/n too small. Dedupe by id.
ids = collections.defaultdict(set)
for f in sorted(os.listdir(P61)):
    if not f.endswith((".faa", ".fa", ".fasta")):
        continue
    for ln in open(os.path.join(P61, f), errors="replace"):
        if ln.startswith(">"):
            i = ln[1:].split()[0]
            ids[lab(i)].add(i)
NREF = {k: len(v) for k, v in ids.items()}
print("=" * 62)
print("UNIQUE references per label (deduplicated across all four fasta files)")
print("=" * 62)
print("   %-14s %6s" % ("label", "unique"))
for k in sorted(NREF):
    print("   %-14s %6d" % (k, NREF[k]))
print()
print("   Step 80's combined table added the files together and so counted the")
print("   same sequences more than once. Medians were unaffected -- duplicates")
print("   are identical -- but 1/n was not, and 1/n is the promiscuity test.")
print()


def alignments(tag):
    """every alignment, not just the best -- the runner-up test needs them"""
    per = collections.defaultdict(list)
    for ln in open(os.path.join(P61, "%s_hits.tsv" % tag), errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            L, s = int(f[3]), float(f[5])
        except ValueError:
            continue
        if L < MINLEN:
            continue
        per[f[0]].append((s, f[1]))
    return per


for tag in ("SED", "LIQ"):
    per = alignments(tag)
    best = {r: max(v) for r, v in per.items()}

    print("=" * 62)
    print("%s" % tag)
    print("=" * 62)

    # ---- A. promiscuity, per reference, against 1/n
    print("   A. references taking more than 4x their expected share")
    byref = collections.Counter(ref for _, ref in best.values())
    bylab = collections.defaultdict(collections.Counter)
    for ref, n in byref.items():
        bylab[lab(ref)][ref] = n
    flagged = collections.Counter()
    for L in FOCUS:
        cnt = bylab.get(L)
        if not cnt:
            continue
        tot = sum(cnt.values())
        n = NREF.get(L, 0)
        if not n or not tot:
            continue
        exp = 1.0 / n
        lines = []
        for ref, c in cnt.most_common():
            sh = c / float(tot)
            if sh > 4 * exp and c >= 3:
                flagged[L] += c
                lines.append("        %-40s %4d/%4d = %4.1f%%  (%.1fx expected)"
                             % (ref.split("|")[-1][:40], c, tot, 100 * sh, sh / exp))
        print("      %-13s n=%-4d total=%-5d  %s"
              % (L, n, tot, "%d reference(s) flagged, %d reads"
                 % (len(lines), flagged[L]) if lines else "none flagged"))
        for x in lines[:5]:
            print(x)

    # ---- B. runner-up outside the winning label
    print()
    print("   B. runner-up outside the winning label")
    print("      %-13s %6s %8s %8s %8s"
          % ("label", "reads", "no other", "gap<10%", "median gap"))
    for L in FOCUS:
        rs = [r for r, (s, ref) in best.items() if lab(ref) == L]
        if not rs:
            continue
        alone = close = 0
        gaps = []
        for r in rs:
            top_s, top_ref = best[r]
            other = [s for s, ref in per[r] if lab(ref) != L]
            if not other:
                alone += 1
                continue
            g = (top_s - max(other)) / top_s if top_s else 0.0
            gaps.append(g)
            if g < 0.10:
                close += 1
        med = statistics.median(gaps) if gaps else float("nan")
        print("      %-13s %6d %7d%% %7d%% %9s"
              % (L, len(rs), round(100 * alone / float(len(rs))),
                 round(100 * close / float(len(rs))),
                 "--" if gaps != gaps or not gaps else "%.0f%%" % (100 * med)))
    print()
    print("      'no other' = the read hits nothing outside this label at all,")
    print("      i.e. sequence the rest of the family cannot see. High is good.")
    print("      'gap<10%' = a sibling scored within a tenth of the winner; the")
    print("      assignment turned on very little. High is bad.")
    print()

print("=" * 62)
print("HOW TO READ THIS")
print("=" * 62)
print("   mtmB is the calibrated case: its measured 0.75 sits on a predicted")
print("   0.76, so whatever its 'no other' and 'gap<10%' figures are, they are")
print("   the values of a marker that DOES measure its gene. Compare mtaB and")
print("   mtbB against mtmB's row, not against an absolute standard.")
print()
print("   If DECOYmtbB shows a low 'no other' and a high 'gap<10%' relative to")
print("   mtmB, its 283 reads are superfamily sequence that happened to land on")
print("   an mtbB reference -- the count measures the family, and the chapter")
print("   says exactly that, with this as the evidence.")
PYEOF
