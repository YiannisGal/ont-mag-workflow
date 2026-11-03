#!/usr/bin/env bash
# 68_network_control.sh -- is the direction of the network biology, or gap-filling?
#
#   bash scripts/run.sh scripts/68_network_control_XXXXXXX.sh
#
# WHY THIS RUNS BEFORE THE FIGURE IS DRAWN.
#
# Step 67's matrix has a very strong shape, and it is exactly the shape an
# ARTEFACT would have:
#
#     ETO_LIQ_MAG_011  (Methanofollis)  receives from all ten others, 28.78 total
#     ETO_LIQ_MAG_003  (MTP4)           receives from nobody at all, donates 39.45
#
# Two methanogens, opposite poles of the same axis. Read as biology that is a
# headline. But SMETANA scores a transfer only when the RECEIVER's model cannot
# make the compound itself. So "receives from everyone" is also what the least
# self-sufficient model looks like, and "receives nothing" is what the most
# self-sufficient model looks like -- and self-sufficiency in a draft GEM is set
# by how much gap-filling it was given, not by the organism.
#
# CHAPTER_STRUCTURE.md §7 already commits to this in the Methods: "With gap-filling
# in the models as well, the honest claim is about CAPACITY, not activity." A
# figure whose whole structure might be a gap-filling gradient has to be tested
# against that, not merely caveated under it.
#
# THE TEST IS AN INDEPENDENT ONE, AND IT ALREADY EXISTS. Step 19b called
# auxotrophies from annotation -- KOfam markers and pathway completeness, no flux
# model, no gap-filling, a completely separate line of evidence. If MAG_011 is
# independently short of the amino-acid pathways whose products it receives, the
# two lines agree and the network's direction is a result. If 19b says it can make
# them all, then what the figure would be showing is the gap-filler's arbitrary
# choices and the panel must change.
#
#   > Rule 112. A structure that is also the signature of your method's known
#   > artefact is not a finding until it has been tested against a line of evidence
#   > that does not share the artefact.
#
# It also settles the community's identity. §5 of CHAPTER_STRUCTURE planned run 1
# on THE IN-SITU TEN. This run has eleven genomes. Which eleven, and whether the
# eleventh is MTP4, decides the caption -- "the in-situ ten plus the methylotrophic
# methanogen" is a designed community; "eleven of twenty-nine" is an unexplained
# subset.
#
# READ-ONLY except 68_control/. Seconds. No environment touched, no network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/68_control"; mkdir -p "$OUT"
exec > >(tee "$OUT/68.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "68_network_control -- gap-filling control on the exchange network"
echo "WS2 = $WS2"
echo "started $(date -Is)"

MAT="$WS2/67_network/fig06_matrix.tsv"
MASTER="$WS2/47_master/MASTER.tsv"
MB="$WS2/43_smetana/run_43o/models_bounded"
AUX="$WS2/19b_auxotrophy"
[ -s "$MAT" ] || { no "run step 67 first -- $MAT missing"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. community
say "A. WHICH ELEVEN GENOMES, AND WHY THOSE"
python3 - "$MAT" "$MASTER" <<'PYEOF'
import sys, csv, os
mat, master = sys.argv[1], sys.argv[2]
E = list(csv.DictReader(open(mat), delimiter="\t"))
inrun = sorted({r["donor_mag"] for r in E} | {r["receiver_mag"] for r in E})
if not os.path.exists(master):
    print("    MASTER.tsv missing; cannot check membership"); sys.exit(0)
M = list(csv.DictReader(open(master), delimiter="\t"))


def short(r):
    return (r.get("mag_id") or "").replace("ETO_", "")


ten = [r for r in M if str(r.get("in_situ_ten", "")).lower().startswith(("y", "t", "1"))]
tenset = {short(r) for r in ten}
runset = set(inrun)
print("    the run holds %d genome(s); MASTER marks %d as in-situ"
      % (len(runset), len(tenset)))
print()
print("      %-16s %-22s %-9s %-10s %s" % ("mag", "genus", "domain", "in-situ?", "breadth"))
by = {short(r): r for r in M}
for g in sorted(runset, key=lambda x: (by.get(x, {}).get("domain", ""), x)):
    r = by.get(g, {})
    print("      %-16s %-22s %-9s %-10s %s"
          % (g, (r.get("genus") or "-")[:22], (r.get("domain") or "?")[:9],
             "yes" if g in tenset else "NO", r.get("sed_breadth", "?")))
print()
extra = sorted(runset - tenset)
absent = sorted(tenset - runset)
print("    in the run but NOT in the in-situ ten: %s" % (", ".join(extra) or "none"))
print("    in the in-situ ten but NOT in the run: %s" % (", ".join(absent) or "none"))
print()
if len(extra) == 1 and not absent:
    print("    >>> the community is THE IN-SITU TEN PLUS ONE. That is a designed")
    print("        community and the caption can say so. CHAPTER_STRUCTURE §5")
    print("        planned run 1 on the ten; the eleventh was added deliberately.")
elif absent:
    print("    >>> the run is NOT the in-situ ten. The caption cannot call it that,")
    print("        and the missing genome(s) above need an explanation in the text.")
PYEOF

# ---------------------------------------------------------------- B. model size
say "B. MODEL SIZE AGAINST POSITION IN THE NETWORK"
sub "counting reactions and metabolites in the models SMETANA consumed"
if [ ! -d "$MB" ]; then
  no "no models at ${MB#$WS2/}"
else
python3 - "$MB" "$MAT" "$MASTER" "$WS2/43_smetana" <<'PYEOF'
import sys, os, re, csv, glob, collections
mb, mat, master, smroot = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

M = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")


def count(p):
    txt = open(p, encoding="utf-8", errors="replace").read()
    sp = len(re.findall(r"<species\b", txt))
    rx = len(re.findall(r"<reaction\b", txt))
    ex = len(re.findall(r'<reaction[^>]*id="R?_?EX_', txt))
    return sp, rx, ex


stat = {}
for p in sorted(glob.glob(os.path.join(mb, "*.xml"))):
    b = os.path.basename(p).replace(".xml", "")
    stat[M.get(b, b)] = count(p)

E = list(csv.DictReader(open(mat), delimiter="\t"))
don = collections.defaultdict(float)
rec = collections.defaultdict(float)
nrec = collections.Counter()
for r in E:
    don[r["donor_mag"]] += float(r["sum_smetana"])
    rec[r["receiver_mag"]] += float(r["sum_smetana"])
    nrec[r["receiver_mag"]] += 1

names = sorted(stat, key=lambda g: -rec.get(g, 0.0))
print("      %-16s %8s %8s %8s %10s %10s %8s"
      % ("mag", "species", "rxns", "exch", "received", "donated", "donors"))
for g in names:
    sp, rx, ex = stat[g]
    print("      %-16s %8d %8d %8d %10.2f %10.2f %8d"
          % (g, sp, rx, ex, rec.get(g, 0.0), don.get(g, 0.0), nrec.get(g, 0)))

# the correlation that matters: does a bigger model receive less?
xs = [stat[g][1] for g in names]
ys = [rec.get(g, 0.0) for g in names]
n = len(xs)
if n > 2:
    mx = sum(xs) / n; my = sum(ys) / n
    sx = (sum((x - mx) ** 2 for x in xs) / n) ** 0.5
    sy = (sum((y - my) ** 2 for y in ys) / n) ** 0.5
    r = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n / (sx * sy)) if sx and sy else 0.0
    print()
    print("    reactions vs score RECEIVED:  r = %+.3f  (n = %d)" % (r, n))
    print("      a strong negative r means the biggest models receive least, i.e. the")
    print("      network's direction is largely a model-size gradient.")
    print("      |r| < 0.4 means model size does not explain who receives.")

# how much gap-filling: bounded models against any pre-fill draft of the same name
print()
print("    GAP-FILLING, measured as the reaction count difference against any")
print("    earlier draft of the same genome elsewhere under 43_smetana:")
found = False
for other in ("models", "models_sed", "models_fillrich"):
    d = os.path.join(smroot, other)
    if not os.path.isdir(d):
        continue
    hits = []
    for p in sorted(glob.glob(os.path.join(d, "*.xml"))):
        b = os.path.basename(p).replace(".xml", "")
        g = M.get(b, b)
        if g in stat:
            hits.append((g, count(p)[1], stat[g][1]))
    if hits:
        found = True
        print("      vs %s:" % other)
        print("        %-16s %10s %10s %10s" % ("mag", other, "bounded", "added"))
        for g, a, b2 in sorted(hits, key=lambda h: -(h[2] - h[1])):
            print("        %-16s %10d %10d %+10d" % (g, a, b2, b2 - a))
if not found:
    print("      no comparable draft set found -- gap-filling cannot be quantified")
    print("      this way, and section C carries the control alone.")
PYEOF
fi

# ---------------------------------------------------------------- C. auxotrophy
say "C. THE INDEPENDENT LINE -- 19b's ANNOTATION-BASED AUXOTROPHIES"
sub "what 19b left behind"
for f in "$AUX/auxotrophies.tsv" "$AUX/pathway_completeness.tsv" "$AUX/marker_audit.tsv"; do
  if [ -s "$f" ]; then
    printf '    %-54s %6s line(s)\n' "${f#$WS2/}" "$(wc -l < "$f")"
  else
    printf '    %-54s MISSING\n' "${f#$WS2/}"
  fi
done
if [ -s "$AUX/auxotrophies.tsv" ]; then
  sub "auxotrophies.tsv header and first rows"
  head -1 "$AUX/auxotrophies.tsv" | tr '\t' '\n' | nl -ba | sed 's/^/    /'
  sed -n '2,6p' "$AUX/auxotrophies.tsv" | sed 's/^/    /'
fi
if [ -s "$AUX/pathway_completeness.tsv" ]; then
  sub "pathway_completeness.tsv header and first rows"
  head -1 "$AUX/pathway_completeness.tsv" | tr '\t' '\n' | nl -ba | sed 's/^/    /'
  sed -n '2,6p' "$AUX/pathway_completeness.tsv" | sed 's/^/    /'
fi

sub "do the two lines agree"
python3 - "$MAT" "$MASTER" "$AUX" "$WS2/67_network/compound_names.tsv" \
          "$WS2/43_smetana/run_43q/pairs/pairs_detailed.tsv" <<'PYEOF'
import sys, os, csv, re, collections
mat, master, aux, mapf, det = sys.argv[1:6]

if not os.path.isdir(aux):
    print("    19b output missing -- the independent line is unavailable, and the")
    print("    figure must be drawn as a model result only, with §7's caveat in the")
    print("    caption rather than a test behind it.")
    sys.exit(0)

M = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")

nm = {}
if os.path.exists(mapf):
    nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(mapf), delimiter="\t")}


def bare(c):
    m = re.search(r"(cpd\d+)", c or "")
    return m.group(1) if m else c


# what each genome RECEIVES, as compound names
rows = list(csv.DictReader(open(det), delimiter="\t")) if os.path.exists(det) else []
low = {c.lower(): c for c in rows[0]} if rows else {}
gets = collections.defaultdict(set)
for r in rows:
    try:
        if float(r[low["smetana"]]) <= 0:
            continue
    except (KeyError, ValueError):
        continue
    g = M.get(r[low["receiver"]], r[low["receiver"]])
    gets[g].add(nm.get(bare(r[low["compound"]]), bare(r[low["compound"]])))

# 19b's calls, keyed by whatever its first column is
def load(p):
    if not os.path.exists(p):
        return [], []
    rs = list(csv.DictReader(open(p), delimiter="\t"))
    return rs, (list(rs[0]) if rs else [])


A, ac = load(os.path.join(aux, "auxotrophies.tsv"))
P, pc = load(os.path.join(aux, "pathway_completeness.tsv"))
print("    auxotrophies.tsv columns     : %s" % ", ".join(ac))
print("    pathway_completeness columns : %s" % ", ".join(pc))
if not A:
    print("    no auxotrophy rows to compare")
    sys.exit(0)

gcol = ac[0]
pcol = next((c for c in ac if "path" in c.lower() or "compound" in c.lower()
             or "nutrient" in c.lower()), ac[-1])
print("    joining on '%s' (genome) and reading '%s' as the pathway" % (gcol, pcol))

byg = collections.defaultdict(set)
for r in A:
    k = r[gcol]
    byg[M.get(k, k.replace("ETO_", ""))].add(str(r[pcol]))

print()
print("      %-16s %-10s %-12s %s" % ("mag", "receives", "19b auxo", "overlap (name match)"))
for g in sorted(gets, key=lambda x: -len(gets[x])):
    au = byg.get(g, set())
    ov = sorted({a for a in au
                 for c in gets[g]
                 if a and a.split("-")[0].lower()[:5] in c.lower()})
    print("      %-16s %-10d %-12d %s"
          % (g, len(gets[g]), len(au), ", ".join(ov[:6]) or "-"))
print()
print("    READ IT THIS WAY. A genome that receives many compounds AND carries many")
print("    independent auxotrophies is corroborated: two methods, one conclusion.")
print("    A genome that receives many compounds and 19b calls prototrophic is a")
print("    gap-filling artefact wearing a result's clothes.")
PYEOF

# ---------------------------------------------------------------- D. verdict
say "D. WHAT THIS DECIDES"
cat <<'EOT'
  Section B, the model-size correlation:

    |r| >= 0.7 negative  -> the network's direction is mostly a model-size
                            gradient. The matrix still goes in the thesis, but as
                            a METHOD result, and the two methanogens' opposite
                            poles cannot be reported as ecology
    |r| <  0.4           -> model size does not explain who receives, and the
                            direction survives its most obvious confound

  Section C is the one that settles it. MAG_011 receives from all ten partners. If
  19b independently calls it auxotrophic for the amino acids it is being handed,
  two unrelated methods agree and that is as close to confirmation as a
  reconstruction gets. If 19b calls it prototrophic, the figure is measuring the
  gap-filler.

  Either way the panel gets drawn -- but the caption, and the sentence in the
  chapter it supports, are different sentences. That is why this runs first.

  One more thing this output will show, for free: whether the eleventh genome is
  MTP4. If it is, fig06's subject is "the in-situ ten plus the methylotrophic
  methanogen the chapter argues about", which is a much stronger caption than
  "eleven genomes".
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
