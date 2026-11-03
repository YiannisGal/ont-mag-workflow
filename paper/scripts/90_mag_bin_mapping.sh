#!/usr/bin/env bash
# 90_mag_bin_mapping.sh -- two files disagree about which genome is which, and
# one of them is the source of a table in the report.
#
# STEP 89 PRINTED BOTH AND THEY DO NOT MATCH.
#
#   11_kofam/verify/mcrA_K00399_candidates.faa   and the 105 mcrA tree tips
#       ETO_LIQ_MAG_005 | LIQ_das_vamb_011
#       ETO_LIQ_MAG_007 | LIQ_das_semibin_028
#       ETO_LIQ_MAG_010 | LIQ_das_vamb_005
#
#   54_mcr_fragmentation/hits_with_coords.tsv
#       ETO_LIQ_MAG_005 | LIQ_das_vamb_005
#       ETO_LIQ_MAG_007 | LIQ_das_vamb_011
#       ETO_LIQ_MAG_010 | LIQ_das_semibin_028
#
# The other five archaea agree exactly (MAG_003 vamb_003, MAG_011 vamb_017,
# MAG_016 metabat_047, MAG_018 metabat_030, and MAG_017 appears only in the
# second file). So this is not a naming convention difference -- three genomes
# have their bin assignments PERMUTED between the two files, and everything
# else lines up.
#
# WHY IT MATTERS RATHER THAN BEING A COSMETIC MISMATCH. If the mag column of
# hits_with_coords.tsv is wrong, then the mcr scores this project has been
# quoting per genome are attached to the wrong genomes -- including the
# 993,9-1046,7 range in FINDINGS_REPORT section 5 and, more importantly, the
# identity of the FRAGMENTED genome. Section 5 says that genome is
# ETO_LIQ_MAG_017, a Methanogenium. That claim is only as good as this join.
#
# If instead the tree tips are wrong, then nine tip labels in a figure name the
# wrong organisms, and fig14 must not be drawn until they are fixed.
#
# THREE INDEPENDENT ROUTES ARE TAKEN HERE AND NONE OF THEM IS EITHER OF THE TWO
# FILES IN DISPUTE:
#
#   A. 47_master/MASTER.tsv -- the curated table, working_name <-> mag_id.
#   D. 07_mags/bins/*.fa    -- which bin FASTA physically contains the contig
#                              each mcr protein sits on. A contig is in exactly
#                              one bin file, so this is a fact about bytes on
#                              disk and cannot be a bookkeeping error.
#   E. the locus tags       -- ETOLIQMAG018_08440 carries its own mag id, and
#                              it was written by the annotation of one genome.
#
# Nothing here is written outside this step's own directory, and nothing is
# deleted. Rule 121: a deletion needs a proof, and no deletion is proposed.
set -u
WS2="${WS2}"
OUT="$WS2/90_mag_bin_mapping"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, csv, collections
WS2, OUT = sys.argv[1], sys.argv[2]

# ----------------------------------------------------------------- A. MASTER
print("=" * 62)
print("A. 47_master/MASTER.tsv -- the curated mapping")
print("=" * 62)
m = os.path.join(WS2, "47_master", "MASTER.tsv")
MASTER = {}          # working_name -> mag_id
GENUS, BREADTH, DOM = {}, {}, {}
if not os.path.exists(m):
    print("   MISSING %s -- nothing below can be settled" % m)
    sys.exit(1)
rows = list(csv.DictReader(open(m, errors="replace"), delimiter="\t"))
for r in rows:
    wn = (r.get("working_name") or "").strip()
    mg = (r.get("mag_id") or "").strip()
    if wn and mg:
        MASTER[wn] = mg
        GENUS[mg] = (r.get("genus") or "?").strip()
        BREADTH[mg] = (r.get("sed_breadth") or "").strip()
        DOM[mg] = (r.get("domain") or "?").strip()
print("   %d genome(s); the eight archaea:" % len(MASTER))
print("   %-22s %-24s %-20s %s" % ("mag_id", "working_name", "genus", "sed_breadth"))
ARC = sorted([mg for mg in GENUS if DOM.get(mg) == "Archaea"])
REV = {v: k for k, v in MASTER.items()}
for mg in ARC:
    print("   %-22s %-24s %-20s %s" % (mg, REV.get(mg, "?"), GENUS[mg], BREADTH.get(mg, "")))
if len(MASTER) != len(REV):
    print("   !! working_name -> mag_id is not one-to-one; the reverse lookup")
    print("      above is unsafe and every join below inherits that.")
print()

# ------------------------------------------------------ B. the tree / the faa
print("=" * 62)
print("B. the mcrA tree tips and the nominated FASTA")
print("=" * 62)
PAIRS_TREE = collections.OrderedDict()


def add(d, mag, binn, src):
    d.setdefault((mag, binn), []).append(src)


tp = None
for c in ("mcrA.contree", "mcrA.treefile"):
    p = os.path.join(WS2, "15_marker_trees", "mcrA", c)
    if os.path.exists(p):
        tp = p
        break
if tp:
    txt = open(tp, errors="replace").read()
    for nm in re.findall(r"QUERY\|[^:,()\s]+", txt):
        f = nm.split("|")
        if len(f) >= 4:
            add(PAIRS_TREE, f[1], f[2], "tree")
    print("   %s: %d QUERY tip(s)" % (os.path.basename(tp),
                                      sum(len(v) for v in PAIRS_TREE.values())))
cand = os.path.join(WS2, "11_kofam", "verify", "mcrA_K00399_candidates.faa")
PAIRS_FAA = collections.OrderedDict()
if os.path.exists(cand):
    for ln in open(cand, errors="replace"):
        if ln.startswith(">"):
            f = ln[1:].strip().split("|")
            if len(f) >= 3:
                add(PAIRS_FAA, f[0], f[1], "faa")
    print("   candidates.faa: %d sequence(s)"
          % sum(len(v) for v in PAIRS_FAA.values()))
print()
print("   %-22s %-24s %-8s %-8s %s" % ("mag_id", "bin", "in tree", "in faa", "MASTER says"))
keys = sorted(set(PAIRS_TREE) | set(PAIRS_FAA))
bad_tree = []
for mag, binn in keys:
    says = MASTER.get(binn, "(bin not in MASTER)")
    ok = (says == mag)
    if not ok:
        bad_tree.append((mag, binn, says))
    print("   %-22s %-24s %-8s %-8s %s %s"
          % (mag, binn,
             len(PAIRS_TREE.get((mag, binn), [])),
             len(PAIRS_FAA.get((mag, binn), [])),
             says, "" if ok else "  <-- DISAGREES"))
print()

# ------------------------------------------------- C. hits_with_coords.tsv
print("=" * 62)
print("C. 54_mcr_fragmentation/hits_with_coords.tsv")
print("=" * 62)
hc = os.path.join(WS2, "54_mcr_fragmentation", "hits_with_coords.tsv")
PAIRS_HC = collections.OrderedDict()
widths = collections.Counter()
hc_rows = []
if os.path.exists(hc):
    with open(hc, errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            widths[len(f)] += 1
            if len(f) >= 5:
                hc_rows.append(f)
                add(PAIRS_HC, f[0], f[1], "hc")
    print("   header has %d column(s): %s" % (len(hdr), ", ".join(hdr)))
    print("   row widths seen: %s"
          % ", ".join("%d fields x%d" % (k, v) for k, v in sorted(widths.items())))
    if len(widths) > 1:
        print("   !! this file is NOT rectangular. Step 89 showed rows such as")
        print("        ... LIQ_contig_1360934  41489  43201  +  570")
        print("      where the contig name and the gene index have run together,")
        print("      i.e. a separator was lost when the file was written. The")
        print("      COORDINATE columns of this file therefore cannot be parsed")
        print("      positionally. Fields 0-5 (mag, bin, ko, protein, score,")
        print("      evalue) and the LAST field (aa_len) are unaffected, and")
        print("      nothing below uses any of the others.")
    print()
    print("   %-22s %-24s %s" % ("mag_id", "bin", "MASTER says"))
    for mag, binn in sorted(PAIRS_HC):
        says = MASTER.get(binn, "(bin not in MASTER)")
        print("   %-22s %-24s %s %s" % (mag, binn, says,
                                        "" if says == mag else "  <-- DISAGREES"))
else:
    print("   MISSING %s" % hc)
print()

# ------------------------------------------------------ D. the bin FASTAs
print("=" * 62)
print("D. which bin FASTA actually contains each mcr contig")
print("=" * 62)
print("   A contig is in exactly one bin file. This route is a fact about")
print("   bytes on disk and cannot be a bookkeeping mistake.")
print()
bd = os.path.join(WS2, "07_mags", "bins")
contigs = set()
for f in hc_rows:
    p = f[3]
    mm = re.match(r"^((?:LIQ|SED)_contig_\d+)_\d+$", p)
    if mm:
        contigs.add(mm.group(1))
print("   %d distinct mcr-bearing contig(s) to locate" % len(contigs))
WHERE = {}
if os.path.isdir(bd):
    fas = sorted(f for f in os.listdir(bd) if f.endswith((".fa", ".fna", ".fasta")))
    print("   %d bin file(s) in 07_mags/bins" % len(fas))
    for fn in fas:
        stem = re.sub(r"\.(fa|fna|fasta)$", "", fn)
        try:
            for ln in open(os.path.join(bd, fn), errors="replace"):
                if ln.startswith(">"):
                    c = ln[1:].split()[0]
                    c = c.rpartition("~")[2] if "~" in c else c
                    if c in contigs:
                        WHERE.setdefault(c, []).append(stem)
        except OSError as e:
            print("   %s unreadable (%s)" % (fn, e))
else:
    print("   MISSING %s" % bd)
print()
print("   %-22s %-10s %-26s %-24s %s"
      % ("contig", "ko", "bin file it is in", "mag by MASTER", "mag in hits_with_coords"))
seen = set()
for f in hc_rows:
    mm = re.match(r"^((?:LIQ|SED)_contig_\d+)_\d+$", f[3])
    if not mm:
        continue
    c = mm.group(1)
    if (c, f[2]) in seen:
        continue
    seen.add((c, f[2]))
    bins = WHERE.get(c, [])
    b = bins[0] if len(bins) == 1 else ("|".join(bins) if bins else "(in no bin file)")
    says = MASTER.get(b, "?")
    flag = "" if says == f[0] else "   <-- hits_with_coords disagrees"
    print("   %-22s %-10s %-26s %-24s %s%s" % (c, f[2], b, says, f[0], flag))
    if len(bins) > 1:
        print("      !! this contig is in %d bin files: %s" % (len(bins), ", ".join(bins)))
print()

# ------------------------------------------------------ E. the locus tags
print("=" * 62)
print("E. the locus tags carry their own genome id")
print("=" * 62)
print("   ETOLIQMAG018_08440 was written by the annotation of ONE genome, so")
print("   the tag itself is a fourth opinion, independent of every table.")
print()
if os.path.exists(cand):
    for ln in open(cand, errors="replace"):
        if not ln.startswith(">"):
            continue
        f = ln[1:].strip().split("|")
        if len(f) < 3:
            continue
        mag, binn, tag = f[0], f[1], f[2]
        canon = mag.replace("_", "")
        tagmag = tag.rsplit("_", 1)[0]
        print("   %-22s %-24s %-24s tag-vs-header %s"
              % (mag, binn, tag, "AGREE" if tagmag == canon else "DISAGREE"))
print()

# ------------------------------------------------------ F. the real scores
print("=" * 62)
print("F. per-genome mcr scores, joined the way the ids are actually written")
print("=" * 62)
print("   Step 89 section F failed because it looked for '~' or an ETO_ prefix")
print("   in 45_kofam_full/all_hits.tbl, and that file writes")
print("       LIQ_das_metabat_030|LIQ_contig_18995_2")
print("   i.e. BIN|contig_index. The bin is the join key and MASTER turns it")
print("   into a mag id. Rule 111: the join is printed before it is believed.")
print()
THR = {}
kl = os.path.join(WS2, "36_gapclose", "ko_list")
if os.path.exists(kl):
    for ln in open(kl, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 1 and re.match(r"^K\d{5}$", f[0]):
            try:
                THR[f[0]] = float(f[1])
            except ValueError:
                pass
tb = os.path.join(WS2, "45_kofam_full", "all_hits.tbl")
KOS = ("K00399", "K00401", "K00402")
hits = []
if os.path.exists(tb):
    for ln in open(tb, errors="replace"):
        if ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 6 or f[2] not in KOS:
            continue
        try:
            hits.append((f[0], f[2], float(f[5])))
        except ValueError:
            pass
else:
    print("   MISSING %s" % tb)
print("   %d hit row(s) for %s" % (len(hits), "/".join(KOS)))
unmapped = set()
per = collections.defaultdict(lambda: collections.defaultdict(list))
for pid, ko, sc in hits:
    binn = pid.split("|")[0]
    mag = MASTER.get(binn)
    if not mag:
        unmapped.add(binn)
        continue
    per[mag][ko].append((sc, pid))
print("   %d bin(s) mapped, %d unmapped%s"
      % (len({p.split("|")[0] for p, _, _ in hits}) - len(unmapped), len(unmapped),
         (": " + ", ".join(sorted(unmapped))) if unmapped else ""))
print()
print("   %-22s %-20s %-10s %8s %10s %6s  %s"
      % ("mag_id", "genus", "KO", "best", "threshold", "n ORF", "all scores"))
for mag in sorted(per, key=lambda g: (DOM.get(g, ""), g)):
    for ko in KOS:
        v = sorted(per[mag].get(ko, []), reverse=True)
        if not v:
            continue
        t = THR.get(ko)
        mark = "" if (t and v[0][0] >= t) else "   BELOW THRESHOLD"
        print("   %-22s %-20s %-10s %8.1f %10s %6d  %s%s"
              % (mag, GENUS.get(mag, "?")[:20], ko, v[0][0],
                 ("%.2f" % t) if t else "?", len(v),
                 ", ".join("%.1f" % s for s, _ in v[:6]), mark))
print()
print("   and the same thing counted the way the report states it:")
above = [g for g in per if per[g].get("K00399")
         and max(s for s, _ in per[g]["K00399"]) >= THR.get("K00399", 1e9)]
arch_above = [g for g in above if DOM.get(g) == "Archaea"]
bact_any = [g for g in per if DOM.get(g) == "Bacteria"]
print("      archaea with a K00399 ORF above %.2f : %d  (%s)"
      % (THR.get("K00399", 0), len(arch_above), ", ".join(sorted(arch_above))))
below = sorted(g for g in per if DOM.get(g) == "Archaea" and g not in arch_above)
print("      archaea with K00399 hits but none above: %s"
      % (", ".join(below) if below else "(none)"))
print("      BACTERIAL genomes with any mcr hit at all: %d  (%s)"
      % (len(bact_any), ", ".join(sorted(bact_any)) if bact_any else "none"))
print()
with open(os.path.join(OUT, "mcr_scores_by_mag.tsv"), "w") as fh:
    fh.write("mag_id\tworking_name\tgenus\tko\tprotein\tscore\tthreshold\tabove\n")
    for mag in sorted(per):
        for ko in KOS:
            for sc, pid in sorted(per[mag].get(ko, []), reverse=True):
                t = THR.get(ko, 0.0)
                fh.write("%s\t%s\t%s\t%s\t%s\t%.1f\t%.2f\t%s\n"
                         % (mag, REV.get(mag, "?"), GENUS.get(mag, "?"), ko,
                            pid, sc, t, "yes" if sc >= t else "no"))
print("   -> 90_mag_bin_mapping/mcr_scores_by_mag.tsv")
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   If MASTER, the bin FASTAs and the locus tags all agree with the tree"
echo "   tips, then hits_with_coords.tsv's mag column is wrong for three"
echo "   genomes, FINDINGS_REPORT section 5's per-genome mcr scores inherit"
echo "   that error, and the identity of the fragmented genome has to be"
echo "   restated from section F's table instead."
echo
echo "   If they agree with hits_with_coords.tsv instead, then nine tip labels"
echo "   of the mcrA tree name the wrong organisms and fig14 cannot be drawn"
echo "   until the tree is relabelled."
echo
echo "   Section F's table is the one panel (a) of fig14 will be built from,"
echo "   whichever way D and E come out."
