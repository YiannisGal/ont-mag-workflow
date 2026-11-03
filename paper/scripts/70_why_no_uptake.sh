#!/usr/bin/env bash
# 70_why_no_uptake.sh -- why do MTP4 and MAG_025 receive nothing?
#
#   bash scripts/run.sh scripts/70_why_no_uptake_XXXXXXX.sh
#
# THE QUESTION. fig06 has two empty columns. SMETANA scores a transfer as
# SCS x MUS x MPS, and MUS -- the metabolite uptake score -- is zero unless the
# RECEIVER's model actually needs to import that compound. So an empty column
# means: across every simulation, this model never needed anything a partner
# makes. There are four ways that happens, and they mean completely different
# things:
#
#   1. THE MEDIUM ALREADY SUPPLIES IT. If ETO_SED contains free amino acids, no
#      model needs a partner for them, and the empty column says nothing about the
#      organism. CHAPTER_STRUCTURE §5 says the sediment medium deliberately has no
#      free amino acids and no B12 -- but that is the PLAN. This reads the file.
#   2. THE MODEL CANNOT IMPORT IT AT ALL. If a model has no exchange reaction for
#      cysteine, SMETANA can never score it as a cysteine receiver, however
#      auxotrophic the organism is. That is a STRUCTURAL zero and it is invisible
#      in the output -- it looks exactly like biological independence.
#   3. THE BIOMASS OBJECTIVE DEMANDS LESS. A model whose biomass equation has
#      fewer precursors needs fewer things. 42d and 43c both fought biomass
#      problems in this project; a methanogen biomass rebuilt in 43k may simply
#      ask for less than the bacterial ones do.
#   4. THE ORGANISM REALLY IS PROTOTROPHIC ON THIS MEDIUM. The interesting case,
#      and the only one that is a finding.
#
# Only the fourth belongs in the chapter as biology. Three and four are
# distinguishable, two is decisive, and one would end the discussion. None of them
# can be told apart from fig06 itself.
#
# 🔴 WHAT STEP 69 ALREADY SETTLED, AND IT IS NOT WHAT I EXPECTED.
#
#   MTP4        genome completeness 99,84   mean amino-acid pathway 0,642
#   MAG_011     genome completeness 90,48   mean amino-acid pathway 0,532
#
# MTP4 is the most complete genome in the set and is better equipped for amino
# acid biosynthesis than the methanogen that receives from everyone. That points
# at explanation 4. But its Cys pathway is 0,00 -- no cysteine route at all -- and
# cysteine is the single most-transferred compound in the whole network. A genome
# with no cysteine pathway whose model needs no cysteine is a contradiction, and
# it is exactly what explanation 2 produces.
#
# AND THE PAIRED TEST CAME BACK NULL. 4 of 9 genomes negative, mean +0,014. At
# COMPOUND level the two methods do not agree. Step 68's apparent corroboration
# for MAG_011 does not survive: it receives Arg and Leu, whose pathways are 1,00.
# Whatever fig06 shows about which compounds move is the model's claim alone.
#
# READ-ONLY except 70_uptake/. Seconds. No solver needed, no network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/70_uptake"; mkdir -p "$OUT"
exec > >(tee "$OUT/70.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "70_why_no_uptake -- the two empty columns in fig06"
echo "WS2 = $WS2"
echo "started $(date -Is)"

MB="$WS2/43_smetana/run_43o/models_bounded"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
NAMES="$WS2/67_network/compound_names.tsv"
SUM="$WS2/67_network/compound_summary.tsv"
PC="$WS2/19b_auxotrophy/pathway_completeness.tsv"
MASTER="$WS2/47_master/MASTER.tsv"

# ---------------------------------------------------------------- A. medium
say "A. WHAT IS ACTUALLY IN THE MEDIUM"
if [ ! -s "$MDB" ]; then
  no "media file not at ${MDB#$WS2/}"
  find "$WS2/41_env_checks" -name '*media*' 2>/dev/null | sed 's/^/    /'
else
  ok "${MDB#$WS2/}  ($(($(wc -l < "$MDB")-1)) row(s))"
  sub "header"
  head -1 "$MDB" | tr '\t' '\n' | nl -ba | sed 's/^/    /'
  sub "media defined in the file"
  awk -F'\t' 'NR>1{print $1}' "$MDB" | sort | uniq -c | sort -rn | sed 's/^/    /'
  python3 - "$MDB" "$NAMES" <<'PYEOF'
import sys, csv, re, collections
mdb, namesf = sys.argv[1], sys.argv[2]
try:
    nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(namesf), delimiter="\t")}
except OSError:
    nm = {}
rows = list(csv.DictReader(open(mdb), delimiter="\t"))
if not rows:
    print("    empty"); sys.exit(0)
cols = list(rows[0])
mcol = cols[0]
ccol = next((c for c in cols if "compound" in c.lower() or c.lower() in
             ("cpd", "metabolite", "id")), cols[-1])
print()
print("    reading '%s' as the medium and '%s' as the compound" % (mcol, ccol))


def bare(c):
    m = re.search(r"(cpd\d+)", str(c) or "")
    return m.group(1) if m else str(c)


AA = ("cystein methionin arginin aspart glutam valin leucin isoleucin threonin "
      "tryptophan tyrosin alanin glycin serin prolin histidin lysin phenylalanin "
      "asparagin ornithin").split()
for med in sorted({r[mcol] for r in rows}):
    cp = sorted({bare(r[ccol]) for r in rows if r[mcol] == med})
    hits = [(c, nm.get(c, "")) for c in cp
            if any(k in (nm.get(c, "") or "").lower() for k in AA)]
    print()
    print("    MEDIUM %-16s %d compound(s)" % (med, len(cp)))
    if hits:
        print("      🔴 CONTAINS FREE AMINO ACIDS -- %d of them:" % len(hits))
        for c, n in hits:
            print("         %-12s %s" % (c, n))
        print("      If the medium supplies them, no model needs a partner for them,")
        print("      and an empty column says nothing about the organism.")
    else:
        print("      no free amino acid among the named compounds"
              " (%d of %d compounds have a name here)"
              % (sum(1 for c in cp if nm.get(c)), len(cp)))
    named = [(c, nm.get(c, "")) for c in cp if nm.get(c)]
    if named:
        print("      named compounds present: %s"
              % ", ".join("%s" % (n or c) for c, n in named[:25]))
PYEOF
fi

# ---------------------------------------------------------------- B. structure
say "B. CAN THE MODEL IMPORT IT AT ALL -- the structural zero"
if [ ! -d "$MB" ]; then
  no "no models at ${MB#$WS2/}"
else
python3 - "$MB" "$MASTER" "$SUM" "$NAMES" "$OUT" <<'PYEOF'
import sys, os, re, csv, glob, collections
mb, master, sumf, namesf, out = sys.argv[1:6]

M = {}
for r in csv.DictReader(open(master), delimiter="\t"):
    if r.get("working_name"):
        M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")

top = []
if os.path.exists(sumf):
    rs = list(csv.DictReader(open(sumf), delimiter="\t"))
    rs.sort(key=lambda r: -float(r.get("sum_smetana", 0) or 0))
    top = [(r["cpd"], r.get("name", "")) for r in rs[:16]]


def parse(p):
    txt = open(p, encoding="utf-8", errors="replace").read()
    ex = set()
    for rid in re.findall(r'<reaction[^>]*id="([^"]+)"', txt):
        m = re.search(r"EX_(cpd\d+)", rid)
        if m:
            ex.add(m.group(1))
    # the objective: the reaction the model maximises, and what it consumes
    bio = None
    for blk in re.findall(r"<reaction\b.*?</reaction>", txt, re.S):
        rid = re.search(r'id="([^"]+)"', blk)
        if rid and re.search(r"bio", rid.group(1), re.I):
            reac = re.search(r"<listOfReactants>(.*?)</listOfReactants>", blk, re.S)
            spp = re.findall(r'species="([^"]+)"', reac.group(1)) if reac else []
            if bio is None or len(spp) > bio[1]:
                bio = (rid.group(1), len(spp))
    return ex, bio


stat = {}
for p in sorted(glob.glob(os.path.join(mb, "*.xml"))):
    b = os.path.basename(p).replace(".xml", "")
    stat[M.get(b, b)] = parse(p)

print("    %d model(s) parsed" % len(stat))
print()
print("      %-16s %10s %-22s %12s" % ("mag", "exchanges", "objective", "precursors"))
for g in sorted(stat):
    ex, bio = stat[g]
    print("      %-16s %10d %-22s %12s"
          % (g, len(ex), bio[0] if bio else "(none found)",
             bio[1] if bio else "-"))
print()
print("    A SMALLER PRECURSOR COUNT MEANS A MODEL THAT NEEDS LESS. If MTP4's")
print("    biomass asks for fewer precursors than the bacterial ones, its empty")
print("    column is a property of its objective function, not of the organism.")

if top:
    print()
    print("    CAN EACH MODEL EVEN EXCHANGE THE MOST-TRADED COMPOUNDS?")
    print("    A missing exchange reaction is a STRUCTURAL zero: SMETANA cannot")
    print("    score an uptake that the model has no reaction for.")
    print()
    hdr = "      %-16s" % "mag"
    for c, n in top[:10]:
        hdr += "%9s" % (n or c)[:8]
    print(hdr)
    for g in sorted(stat):
        ex = stat[g][0]
        line = "      %-16s" % g
        for c, n in top[:10]:
            line += "%9s" % ("yes" if c in ex else "NO")
        print(line)
    print()
    miss = {g: [c for c, n in top if c not in stat[g][0]] for g in stat}
    for g in sorted(miss):
        if miss[g]:
            print("      %-16s cannot exchange: %s"
                  % (g, ", ".join(dict(top).get(c, c) or c for c in miss[g])))
    if not any(miss.values()):
        print("      every model has an exchange for every one of the top compounds,")
        print("      so no empty column is explained by a missing reaction.")

with open(os.path.join(out, "model_structure.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["mag", "n_exchange", "objective", "n_precursors"])
    for g in sorted(stat):
        ex, bio = stat[g]
        w.writerow([g, len(ex), bio[0] if bio else "", bio[1] if bio else ""])
print()
print("    wrote %s" % os.path.join(out, "model_structure.tsv"))
PYEOF
fi

# ---------------------------------------------------------------- C. the two
say "C. THE TWO NON-RECEIVERS, AGAINST 19b"
python3 - "$PC" "$MASTER" <<'PYEOF'
import sys, csv, os
pcf, master = sys.argv[1], sys.argv[2]
if not os.path.exists(pcf):
    print("    19b pathway_completeness.tsv missing"); sys.exit(0)
PC = list(csv.DictReader(open(pcf), delimiter="\t"))
pcols = [c for c in PC[0] if c not in
         ("mag_id", "working_name", "sample", "domain", "genus", "completeness")]
AAC = [c for c in pcols if not c.startswith("B")]
by = {(r.get("mag_id", "") or "").replace("ETO_", ""): r for r in PC}

WANT = ["LIQ_MAG_003", "LIQ_MAG_025", "LIQ_MAG_011"]
have = [g for g in WANT if g in by]
if not have:
    print("    none of %s found in 19b" % ", ".join(WANT)); sys.exit(0)

print("      %-18s %s" % ("pathway", "  ".join("%-13s" % g for g in have)))
for p in pcols:
    line = "      %-18s" % p
    for g in have:
        try: line += "  %-13.2f" % float(by[g][p])
        except (KeyError, TypeError, ValueError): line += "  %-13s" % "-"
    print(line)
print()
for g in have:
    vals = []
    for c in AAC:
        try: vals.append(float(by[g][c]))
        except (KeyError, TypeError, ValueError): pass
    zeros = [c for c in AAC if str(by[g].get(c, "")).strip() in ("0", "0.0", "0.00")]
    print("      %-16s genome %6s   mean amino-acid %.3f   pathways at ZERO: %s"
          % (g, by[g].get("completeness", "?"),
             sum(vals) / len(vals) if vals else float("nan"),
             ", ".join(zeros) or "none"))
print()
print("    A pathway at 0,00 means 19b found no route at all. If a model whose")
print("    genome has NO route for a compound still needs none of it, the model and")
print("    the annotation disagree about that genome -- and the model is the one")
print("    with gap-filling in it.")
PYEOF

# ---------------------------------------------------------------- D. provenance
say "D. WHERE MTP4'S MODEL CAME FROM"
sub "which directories hold a model for each of the eleven"
python3 - "$WS2/43_smetana" "$MASTER" <<'PYEOF'
import sys, os, glob, csv, collections
root, master = sys.argv[1], sys.argv[2]
M = {}
for r in csv.DictReader(open(master), delimiter="\t"):
    if r.get("working_name"):
        M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")
# TWO LEVELS DEEP, because models_bounded lives under run_43o and a one-level
# glob silently reports that the bounded set does not exist
dirs = sorted(d for d in (glob.glob(os.path.join(root, "*"))
                          + glob.glob(os.path.join(root, "*", "*")))
              if os.path.isdir(d) and glob.glob(os.path.join(d, "*.xml")))
seen = collections.defaultdict(set)
for d in dirs:
    tag = os.path.relpath(d, root)
    for p in glob.glob(os.path.join(d, "*.xml")):
        b = os.path.basename(p).replace(".xml", "")
        seen[M.get(b, b)].add(tag)
print("      %-16s %s" % ("mag", "directories holding a model for it"))
for g in sorted(seen):
    print("      %-16s %s" % (g, ", ".join(sorted(seen[g]))))
print()
print("    A genome present in ONLY the bounded set has no earlier draft to measure")
print("    gap-filling against -- step 68's +0 rows were real comparisons; a genome")
print("    missing from that comparison was never compared at all.")
PYEOF

# ---------------------------------------------------------------- E. verdict
say "E. HOW TO READ THE ANSWER"
cat <<'EOT'
  Section A first, because it can end the discussion. If ETO_SED contains free
  amino acids, nobody needs a partner for them and the empty columns are a
  property of the medium. CHAPTER_STRUCTURE §5 says it should not -- that
  exclusion is described there as "the whole point" -- so agreement here also
  confirms the medium was built as designed.

  Section B is the decisive structural check. An exchange reaction that does not
  exist cannot carry an uptake, so a missing EX_ reaction produces an empty column
  that looks exactly like independence. If MTP4 has an exchange for cysteine and
  still takes none, that is a real prediction. If it has none, fig06's empty
  column for MTP4 is a modelling gap and the caption must say so.

  The precursor count in the same table answers the other half: a biomass equation
  asking for fewer precursors needs fewer partners, and the archaeal biomass was
  rebuilt separately in 43k.

  Section C is the contradiction to resolve. MTP4's Cys pathway is 0,00 and
  cysteine is the most-transferred compound in the network. Either its model can
  make cysteine because something was filled in, or it cannot take cysteine up.
  Section B says which.

  WHAT IS ALREADY DECIDED, whatever this returns: the paired test in step 69 came
  back null, so fig06's compound identities are the model's claim and not
  corroborated by annotation. The panel is a capacity statement. The caption
  already says so; the chapter text has to as well.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
