#!/usr/bin/env bash
# 73_cys_fullset.sh -- the control failed because the KOs were never searched.
#
#   bash scripts/run.sh scripts/73_cys_fullset_XXXXXXX.sh
#
# WHAT STEP 72 ACTUALLY FOUND. "0 of 29 genomes pass the control." K01883 is
# cysteinyl-tRNA synthetase: every organism that makes protein has one. Twenty-nine
# genomes cannot all lack it. And the parse was right this time -- the header is
# there in black and white, hmmsearch --tblout, query name in column 2, full
# sequence score in column 5, and the rule-115 self-test came back clean.
#
# So the tables do not contain K01883. `11_kofam/hits/*.tbl` is the output of
# step 11, which searched a TARGETED profile set -- `kofam_targets.hmm` -- not all
# 28.339 KOfam profiles. Every KO I asked about was outside it.
#
#   > Rule 56. A KO absent from the searched profile set is not evidence of
#   > absence.
#
# The control did its job. It refused to let five "-" be read as five absences,
# which is exactly the reading step 72 would otherwise have invited and exactly
# the reading that has gone wrong twice already in this chain. A control that
# fails loudly is worth more than one that passes quietly.
#
# WHERE THE ANSWER ACTUALLY LIVES. `45_kofam_full/` -- the name is the clue, it is
# the full-database run, and step 71 found all five KOs in it. The only reason
# step 71's numbers were junk is that I read a coordinate out of the description
# field instead of the score.
#
# 🔴 AND THAT DIRECTORY HAS A SECOND TRAP, ALREADY RECORDED. Rule 101: it holds
# BOTH `all_hits.tbl` AND `chunks/c*.tbl`, and globbing the lot counts every hit
# twice -- that is how "72 mcr hits" became 36 in step 54. This step reads
# all_hits.tbl and nothing else, and says so.
#
# THIS STEP THEREFORE:
#   A. measures what step 11 actually searched, so rule 56 is quantified and not
#      just asserted -- how many KOs, and whether the six of interest are in it
#   B. reads 45_kofam_full's own format the same way, column by column
#   C. runs the SAME control there. If K01883 passes in ~29 genomes the parse and
#      the source are both sound and only then does anything below get believed
#   D. answers the cysteine question for all 29 genomes
#
# READ-ONLY except 73_cys/. Seconds. No network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/73_cys"; mkdir -p "$OUT"
exec > >(tee "$OUT/73.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "73_cys_fullset -- the same question, against the full-database run"
echo "WS2 = $WS2"
echo "started $(date -Is)"

HITS="$WS2/11_kofam/hits"
FULL="$WS2/45_kofam_full"
KOLIST="$WS2/36_gapclose/ko_list"
MASTER="$WS2/47_master/MASTER.tsv"

# ---------------------------------------------------------------- A. rule 56
say "A. WHAT STEP 11 ACTUALLY SEARCHED -- rule 56, quantified"
python3 - "$HITS" <<'PYEOF'
import sys, os, glob, re, collections
hits = sys.argv[1]
files = sorted(glob.glob(os.path.join(hits, "*.tbl")))
kos = collections.Counter()
per = {}
for p in files:
    s = set()
    with open(p, errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            if len(f) > 2 and re.fullmatch(r"K\d{5}", f[2]):
                s.add(f[2])
    per[os.path.basename(p)] = len(s)
    kos.update(s)
print("    %d table(s); %d DISTINCT KO(s) appear anywhere in them"
      % (len(files), len(kos)))
print("    per genome: %d .. %d KO(s)"
      % (min(per.values()) if per else 0, max(per.values()) if per else 0))
print()
WANT = {"K01883": "cysteinyl-tRNA synthetase (the control)",
        "K01738": "cysK", "K00640": "cysE", "K10150": "cysM",
        "K07587": "sepS, archaeal", "K06868": "sepCysS, archaeal"}
for k in sorted(WANT):
    print("      %-9s %-34s %s" % (k, WANT[k],
          "searched, in %d genome(s)" % kos[k] if kos[k]
          else "NOT IN THIS PROFILE SET AT ALL"))
print()
print("    That is the whole of step 72's result: five dashes and a failed control,")
print("    because none of these profiles was in the targeted set. Nothing about")
print("    cysteine was measured there, in either direction.")
PYEOF

# ---------------------------------------------------------------- B. format
say "B. 45_kofam_full -- what is in it, and its format"
sub "files present (rule 101: all_hits.tbl AND chunks would double-count)"
ls -la "$FULL" 2>/dev/null | sed 's/^/    /' | head -20
[ -d "$FULL/chunks" ] && echo "    chunks/: $(find "$FULL/chunks" -name '*.tbl' | wc -l) table(s) -- DELIBERATELY NOT READ"
AH="$FULL/all_hits.tbl"
if [ ! -s "$AH" ]; then
  no "no all_hits.tbl; falling back to the chunks, read ONCE each"
  AH=""
fi
[ -n "$AH" ] && ok "reading ${AH#$WS2/} only  ($(wc -l < "$AH") line(s))"
sub "its comment header, verbatim"
[ -n "$AH" ] && grep -m4 '^#' "$AH" | sed 's/^/    /'
sub "first two data lines, every field numbered"
[ -n "$AH" ] && grep -v '^#' "$AH" | head -2 | \
  awk '{for(i=1;i<=NF && i<=12;i++) printf "    %d:%s\n", i-1, $i; print "    ---"}'

# ---------------------------------------------------------------- C+D
say "C. THE CONTROL, THEN THE ANSWER"
python3 - "$FULL" "$KOLIST" "$MASTER" "$OUT" <<'PYEOF'
import sys, os, re, csv, glob, collections
full, kolist, master, out = sys.argv[1:5]

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

src = os.path.join(full, "all_hits.tbl")
if not os.path.exists(src) or os.path.getsize(src) == 0:
    cand = sorted(glob.glob(os.path.join(full, "chunks", "*.tbl")))
    print("    all_hits.tbl unusable; reading %d chunk(s), each exactly once"
          % len(cand))
else:
    cand = [src]
    print("    reading %s only (rule 101)" % os.path.relpath(src, full))

M, dom = {}, {}
for r in csv.DictReader(open(master), delimiter="\t"):
    wn = r.get("working_name")
    if wn:
        M[wn] = (r.get("mag_id", "") or "").replace("ETO_", "")
        dom[wn] = r.get("domain", "")

KO_I, SC_I = 2, 5
CTRL = "K01883"
BACT = [("K01738", "cysK"), ("K00640", "cysE"), ("K10150", "cysM")]
ARCH = [("K07587", "sepS / O-phosphoseryl-tRNA synthetase"),
        ("K06868", "sepCysS / Sep-tRNA:Cys-tRNA synthetase")]
WANT = {CTRL} | {k for k, _ in BACT + ARCH}

best = collections.defaultdict(dict)
seen = collections.defaultdict(list)
nline = 0
for p in cand:
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
            nline += 1
            try:
                sc = float(f[SC_I])
            except ValueError:
                continue
            gene = f[0]
            bin_ = gene.split("|")[0] if "|" in gene else None
            if bin_ is None:
                continue
            seen[gene].append((ko, sc))
            if sc > best[bin_].get(ko, -1.0):
                best[bin_][ko] = sc
print("    %d hit line(s) parsed, %d genome prefix(es)" % (nline, len(best)))

# rule 115 again, on this file
susp = [g for g, v in seen.items()
        if len(v) > 2 and len({round(s, 3) for _, s in v}) == 1]
print("    rule-115 self-test: %d gene(s) with an identical value across 3+ profiles"
      % len(susp))
if susp:
    print("    !! still reading a non-score column. Stop.")
    for g in susp[:3]:
        print("       %s -> %s" % (g, seen[g][:4]))
    sys.exit(0)

npass = sum(1 for b in best.values() if b.get(CTRL, -1) >= (thr.get(CTRL) or 0))
print()
print("    control %s (%s) threshold %s" % (CTRL, defn.get(CTRL, "?"), thr.get(CTRL)))
print("    %d of %d genome(s) pass." % (npass, len(best)))
if npass < len(best) * 0.7:
    print("    !! the control still fails. Nothing below may be believed, and the")
    print("       next thing to check is whether this file is per-genome at all.")
    sys.exit(0)
print("    a housekeeping enzyme is found where one must be. Source and parse hold.")

print()
print("================================================================")
print("D. CYSTEINE, ALL GENOMES, AGAINST CURATED THRESHOLDS")
print("================================================================")
hdr = "      %-22s %-9s" % ("genome", "domain")
for k, _ in BACT + ARCH:
    hdr += "%10s" % k
print(hdr)
rowsout = []
for bin_ in sorted(best, key=lambda b: (dom.get(b, "") != "Archaea", M.get(b, b))):
    mag = M.get(bin_, bin_)
    d = dom.get(bin_, "?")
    line = "      %-22s %-9s" % (mag, d[:9])
    rec = {"mag": mag, "domain": d, "bin": bin_}
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
print("    bare = passed its curated threshold; (brackets) = hit but below it;")
print("    dash = no hit. Rule 97: the threshold is the ko_list value, not 1e-5.")

arc = [r for r in rowsout if str(r["domain"]).startswith("Arch")]
bac = [r for r in rowsout if not str(r["domain"]).startswith("Arch")]
print()
for k, nm in BACT + ARCH:
    a = sum(1 for r in arc if r[k + "_pass"] == "yes")
    b = sum(1 for r in bac if r[k + "_pass"] == "yes")
    print("      %-9s %-42s archaea %2d/%-3d bacteria %2d/%-3d"
          % (k, nm[:42], a, len(arc), b, len(bac)))

mtp = next((r for r in rowsout if r["mag"] == "LIQ_MAG_003"), None)
if mtp:
    archok = any(mtp[k + "_pass"] == "yes" for k, _ in ARCH)
    bactok = any(mtp[k + "_pass"] == "yes" for k, _ in BACT)
    print()
    print("    MTP4: bacterial route %s, archaeal route %s"
          % ("PRESENT" if bactok else "absent", "PRESENT" if archok else "absent"))
    if archok and not bactok:
        print("    >>> 19b's Cys = 0,00 is a MARKER-SET GAP. routes.tsv knows only")
        print("        cysK/cysE/cysM; MTP4 makes cysteine the archaeal way. The")
        print("        model needing none and the genome making its own AGREE, and")
        print("        fig06's empty MTP4 column is a result.")
    elif bactok:
        print("    >>> the bacterial route is there after all, so 19b's zero is a")
        print("        scoring question. Compare the score above with routes.tsv's")
        print("        own cutoff before writing anything.")
    else:
        print("    >>> neither route passes even in the full search. Then cysteine")
        print("        biosynthesis in MTP4 is genuinely unresolved and fig06's")
        print("        empty column rests on the model alone -- say so in the text.")

with open(os.path.join(out, "cys_routes_full.tsv"), "w", newline="") as fh:
    cols = (["mag", "domain", "bin"] + [k for k, _ in BACT + ARCH]
            + [k + "_pass" for k, _ in BACT + ARCH])
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rowsout:
        w.writerow(r)
print()
print("    wrote %s" % os.path.join(out, "cys_routes_full.tsv"))
PYEOF

say "E. AND ONE THING THAT IS TRUE WHATEVER THIS RETURNS"
cat <<'EOT'
  19b's cysteine pathway is cysK / cysE / cysM and nothing else. Those are the
  bacterial route. For any archaeon in this catalogue, 19b's Cys column reports
  the absence of a route archaea are not expected to have. That belongs in the
  record regardless of what MTP4 turns out to carry, because the thesis promises
  equal treatment of bacteria and archaea and this is a place where the pipeline
  does not deliver it.

  It is NOT the Bakta-scale failure of ARCHAEAL_ANNOTATION_GAP.md. Step 71 section
  A showed archaea higher than bacteria on 8 of 20 pathways, including B12 by
  +0,387 -- a pattern, not a uniform one-sided deficit. What is established is
  narrow: the Cys route is bacterial-only, and the other 19 have not been checked
  one by one.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
