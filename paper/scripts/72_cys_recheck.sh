#!/usr/bin/env bash
# 72_cys_recheck.sh -- step 71 section C read the wrong column. Redo it.
#
#   bash scripts/run.sh scripts/72_cys_recheck_XXXXXXX.sh
#
# 🔴 THE NUMBERS IN STEP 71 SECTION C ARE NOT SCORES.
#
#     LIQ_contig_1777_86   score 96142.0   thr 524.0
#     LIQ_contig_1779_49   score 64975.0   thr 442.23
#
# KOfam bit scores do not reach 96,000. And the same gene carries the SAME
# number against eight different profiles:
#
#     LIQ_contig_1777_86 = 96142,0 under K01766, K04487, K06868, K11325,
#                                        K11717, K18913, K20247 and K22207
#
# A bit score is a property of the (gene, profile) PAIR. A number that is
# identical across eight profiles is a property of the GENE alone -- a
# coordinate out of the prodigal header, which `45_kofam_full` carries in its
# description field. My parser took "the largest plausible float on the line",
# and on that file the largest float is a genome coordinate.
#
# So every "above threshold" in that section is unverified, and the conclusion
# drawn from it -- "MTP4 CARRIES A CYSTEINE KO ABOVE THRESHOLD" -- is withdrawn.
# This is rule 97 (use the curated threshold) defeated by not knowing which column
# the score was in. I wrote a heuristic instead of reading the format, in the step
# whose whole purpose was to stop exactly this kind of error.
#
#   > Rule 114. Identify the column before you compare its contents to anything.
#   > "The largest number on the line" is not a measurement, it is a guess with a
#   > decimal point.
#
#   > Rule 115. A score that does not change when the profile changes is not a
#   > score. Test for it -- one line of code, and it catches the whole class.
#
# WHAT SURVIVES FROM STEP 71, because it never depended on that parse:
#
#   B. 19b's cysteine route is K01738 (cysK), K00640 (cysE), K10150 (cysM). All
#      three are the BACTERIAL O-acetylserine route. The archaeal tRNA-dependent
#      route -- K07587 O-phosphoseryl-tRNA synthetase and K06868 Sep-tRNA:Cys-tRNA
#      synthetase -- is not in routes.tsv at all, though both profiles exist in
#      ko_list with thresholds of 658,73 and 490,73. For an archaeon, 19b's Cys is
#      0,00 BY CONSTRUCTION. MTP4's zero means "no bacterial route", not "no
#      cysteine".
#   D. MTP4's model is 1.125 reactions in all TEN directories it appears in. It was
#      never re-filled. MAG_011 moved 880 -> 923, MAG_025 1.187 -> 1.217.
#   E. ETO_SED holds no amino acid: water, phosphate, CO2, ammonia, glucose,
#      xylose, cellobiose, galacturonate, sulfate, sulfide, and ions. Exactly as
#      CHAPTER_STRUCTURE §5 specifies.
#   A. And the archaeal deficit is NOT one-sided: archaea are lower on 12 of 20
#      pathways and HIGHER on 8, including B12 by +0,387. That is a pattern, not a
#      tool failure, and my Rule 66 hypothesis is not supported for the amino acid
#      set. Cys at -0,269 is mid-pack, not an outlier.
#
# WHAT THIS STEP DOES. Reads the column layout instead of guessing it, tests the
# parse against a control every organism must have, and then answers one narrow
# question for all eight archaea: is the ARCHAEAL cysteine route present?
#
# READ-ONLY except 72_cys/. Seconds. No network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/72_cys"; mkdir -p "$OUT"
exec > >(tee "$OUT/72.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "72_cys_recheck -- the cysteine question, with the score column identified"
echo "WS2 = $WS2"
echo "started $(date -Is)"

HITS="$WS2/11_kofam/hits"
KOLIST="$WS2/36_gapclose/ko_list"
MASTER="$WS2/47_master/MASTER.tsv"
PC="$WS2/19b_auxotrophy/pathway_completeness.tsv"

say "A. THE FORMAT, READ NOT GUESSED"
sub "the comment header of one table, verbatim"
F=$(find "$HITS" -name '*.tbl' | sort | head -1)
[ -n "$F" ] || { no "no tables under ${HITS#$WS2/}"; echo "STEP_STATUS: FAILED"; exit 1; }
echo "    $(basename "$F")"
grep -m6 '^#' "$F" | sed 's/^/    /'
sub "first two data lines, every field numbered"
grep -v '^#' "$F" | head -2 | awk '{for(i=1;i<=NF && i<=12;i++) printf "    %d:%s\n", i-1, $i; print "    ---"}'

say "B. THE PARSE, AND A CONTROL THAT PROVES IT"
python3 - "$HITS" "$KOLIST" "$MASTER" "$PC" "$OUT" <<'PYEOF'
import sys, os, re, csv, glob, collections
hits, kolist, master, pcf, out = sys.argv[1:6]

# ---- ko_list: threshold and definition, keyed by KO
thr, defn = {}, {}
with open(kolist) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    it = head.index("threshold") if "threshold" in head else 1
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) <= it:
            continue
        try: thr[p[0]] = float(p[it])
        except ValueError: thr[p[0]] = None
        defn[p[0]] = p[-1]
print("    ko_list: %d KO(s)" % len(thr))

# ---- hmmsearch --tblout: target 0, query(KO) 2, E-value 4, SCORE 5
KO_I, SC_I = 2, 5


def rows(p):
    with open(p, errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            if len(f) <= SC_I:
                continue
            ko = f[KO_I]
            if not re.fullmatch(r"K\d{5}", ko):
                continue
            try:
                sc = float(f[SC_I])
            except ValueError:
                continue
            yield f[0], ko, sc


files = sorted(glob.glob(os.path.join(hits, "*.tbl")))
print("    %d table(s)" % len(files))

# ---- RULE 115 SELF-TEST: a score must change when the profile changes
same = 0; pairs = 0
byg = collections.defaultdict(set)
for p in files[:4]:
    for gene, ko, sc in rows(p):
        byg[gene].add(round(sc, 3))
multi = [g for g in byg if len(byg[g]) >= 1]
const = sum(1 for g in byg if len(byg[g]) == 1)
seen = collections.defaultdict(list)
for p in files[:4]:
    for gene, ko, sc in rows(p):
        seen[gene].append((ko, sc))
suspect = [g for g, v in seen.items()
           if len(v) > 2 and len({round(s, 3) for _, s in v}) == 1]
print("    rule-115 self-test: %d gene(s) hit 3+ profiles with an IDENTICAL value"
      % len(suspect))
if suspect:
    print("    !! that is the signature of a non-score column. STOP and re-read A.")
    for g in suspect[:3]:
        print("       %s -> %s" % (g, seen[g][:4]))
else:
    print("    scores vary across profiles, as a bit score must. Column %d is the"
          % SC_I)
    print("    full-sequence score.")

# ---- CONTROL: cysteinyl-tRNA synthetase. Every organism that makes protein has
# one. If the parse cannot find it, the parse is wrong, not the biology.
CTRL = "K01883"
print()
print("    control %s (%s), threshold %s" % (CTRL, defn.get(CTRL, "?"), thr.get(CTRL)))

M = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            M[r["working_name"]] = r

dom = {}
if os.path.exists(pcf):
    for r in csv.DictReader(open(pcf), delimiter="\t"):
        wn = r.get("working_name")
        if wn:
            dom[wn] = r.get("domain", "")

best = {}
for p in files:
    bin_ = os.path.basename(p).replace(".tbl", "")
    b = {}
    for gene, ko, sc in rows(p):
        if sc > b.get(ko, -1):
            b[ko] = sc
    best[bin_] = b

npass = sum(1 for b in best.values()
            if b.get(CTRL, -1) >= (thr.get(CTRL) or 0))
print("    %d of %d genome(s) pass the control." % (npass, len(best)))
if npass < len(best) * 0.7:
    print("    !! the control fails in most genomes. The parse is still wrong and")
    print("       nothing below may be believed.")
else:
    print("    the parse finds a housekeeping enzyme where one must be. It works.")

# ---- the question
BACT = [("K01738", "cysK"), ("K00640", "cysE"), ("K10150", "cysM")]
ARCH = [("K07587", "sepS / O-phosphoseryl-tRNA synthetase"),
        ("K06868", "sepCysS / Sep-tRNA:Cys-tRNA synthetase")]
print()
print("================================================================")
print("C. THE BACTERIAL ROUTE AND THE ARCHAEAL ROUTE, PER GENOME")
print("================================================================")
print("    19b's Cys pathway uses only the first three. The last two are the")
print("    archaeal route and are in ko_list but not in routes.tsv.")
print()
hdr = "      %-22s %-9s" % ("genome", "domain")
for k, _ in BACT + ARCH:
    hdr += "%10s" % k
print(hdr)
rowsout = []
for bin_ in sorted(best, key=lambda b: (dom.get(b, "") != "Archaea", b)):
    mag = (M.get(bin_, {}).get("mag_id", bin_) or bin_).replace("ETO_", "")
    d = dom.get(bin_, "?")
    line = "      %-22s %-9s" % (mag, d[:9])
    rec = {"mag": mag, "domain": d}
    for k, _ in BACT + ARCH:
        sc = best[bin_].get(k, -1.0)
        t = thr.get(k)
        hit = (t is not None and sc >= t)
        line += "%10s" % (("%.0f" % sc) if hit else ("(%.0f)" % sc if sc > 0 else "-"))
        rec[k] = "%.1f" % sc if sc > 0 else ""
        rec[k + "_pass"] = "yes" if hit else "no"
    rowsout.append(rec)
    print(line)
print()
print("    a bare number passed its curated threshold; (a number in brackets) hit")
print("    the profile but below threshold; a dash means no hit at all.")

arc = [r for r in rowsout if str(r["domain"]).startswith("Arch")]
bac = [r for r in rowsout if not str(r["domain"]).startswith("Arch")]


def frac(rs, k):
    return sum(1 for r in rs if r[k + "_pass"] == "yes"), len(rs)


print()
for k, nm in BACT + ARCH:
    a, na = frac(arc, k)
    b, nb = frac(bac, k)
    print("      %-9s %-44s archaea %d/%d   bacteria %d/%d"
          % (k, nm[:44], a, na, b, nb))

mtp = next((r for r in rowsout if r["mag"] == "LIQ_MAG_003"), None)
print()
if mtp:
    print("    MTP4 (LIQ_MAG_003):")
    for k, nm in BACT + ARCH:
        print("      %-9s %-44s %s" % (k, nm[:44],
                                       "PASS" if mtp[k + "_pass"] == "yes" else "no"))
    archok = any(mtp[k + "_pass"] == "yes" for k, _ in ARCH)
    bactok = any(mtp[k + "_pass"] == "yes" for k, _ in BACT)
    print()
    if archok and not bactok:
        print("    >>> MTP4 CARRIES THE ARCHAEAL ROUTE AND NOT THE BACTERIAL ONE.")
        print("        19b's Cys = 0,00 is a marker-set gap, exactly as section B of")
        print("        step 71 predicted from routes.tsv. The model needing no")
        print("        cysteine and the genome making its own AGREE.")
    elif bactok:
        print("    >>> MTP4 carries the bacterial route after all. Then 19b's 0,00 is")
        print("        a scoring question, not a route question -- check what score")
        print("        cysK got against its threshold above.")
    else:
        print("    >>> NEITHER route passes. The honest line is that cysteine")
        print("        biosynthesis in MTP4 is unresolved: the bacterial route is")
        print("        absent and the archaeal one does not pass its threshold")
        print("        either. fig06's empty MTP4 column then rests on the model")
        print("        alone and the caption must say so.")

with open(os.path.join(out, "cys_routes.tsv"), "w", newline="") as fh:
    cols = ["mag", "domain"] + [k for k, _ in BACT + ARCH] + \
           [k + "_pass" for k, _ in BACT + ARCH]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rowsout:
        w.writerow(r)
print()
print("    wrote %s" % os.path.join(out, "cys_routes.tsv"))
PYEOF

say "D. WHAT THIS SETTLES"
cat <<'EOT'
  If the archaeal route passes in MTP4 and the bacterial one does not, then the
  chain closes:

    19b says Cys 0,00        because routes.tsv only knows cysK/cysE/cysM
    the genome says otherwise because K07587/K06868 are there and pass
    the model says no uptake  because it can make cysteine
    fig06's empty column      is a result, for MTP4

  and MAG_025's empty column is still NOT a result, because its genome is 52,5%
  complete. Two zeros, two different meanings, which is why the panel now carries
  completeness beside every genome.

  If NEITHER route passes, the honest sentence is that MTP4's cysteine
  biosynthesis is unresolved, and fig06's MTP4 column rests on the model alone.

  Either way one thing has to be written into the record: 19b's amino acid routes
  are bacterial-canonical, and every archaeal auxotrophy call in this project
  inherits that. Not a catastrophe -- section A of step 71 shows the archaeal
  profile is mixed rather than uniformly depressed, so this is not Rule 66 at the
  scale of the Bakta failure -- but Cys is demonstrably affected and the others
  have not been checked one by one.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
