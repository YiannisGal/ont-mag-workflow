#!/usr/bin/env bash
# 67_network_compounds.sh -- what is actually being traded, and in which direction.
#
#   bash scripts/run.sh scripts/67_network_compounds_XXXXXXX.sh
#
# 🔴 FIRST, A CORRECTION TO STEP 66.
#
# Section C of step 66 printed:
#
#     "of those, the ones this chapter argues about:
#        none of H2, formate, acetate, methanol or methylamine appears at all"
#
# THAT IS FALSE, and the proof is eleven lines further down its own output:
# fig06_input.tsv row 11 is M_cpd00116_e, which is methanol. My lookup table held
# bare ModelSEED ids (cpd00116) and the file holds SBML species ids (M_cpd00116_e).
# Two id spaces for the same thing, and I compared across them without converting.
#
# This is the SAME mistake as the two KOfam id spaces in Volume XII -- Bakta locus
# tags in 11_kofam against BIN|PRODIGAL_ID in 45_kofam_full -- and it produced the
# same shape of wrong answer: a confident NONE that was really a failed join.
# 43q's own awk was immune because it used a regex match, not equality.
#
#   > Rule 111. When a lookup returns nothing, test the join before believing the
#   > answer. An empty intersection is evidence about the KEYS at least as often as
#   > it is evidence about the data.
#
# SECOND, A GAP. 64 distinct compounds moved, and step 66 reported them as
# accession numbers. cpd00084 sums to 15.55 and is the single most-transferred
# thing in the reconstruction -- and nobody can read that. The chapter's claim is
# about WHAT is exchanged, so the accessions have to become names, and they have
# to come from this machine rather than from recollection of ModelSEED.
#
# The SBML models SMETANA actually consumed carry the names in their species
# elements. That is the authoritative local source: same files, same run, no
# network call, no recollection.
#
# WHAT THIS STEP DECIDES. Whether fig06 is a syntrophy figure or an amino-acid
# figure. Those are different panels making different claims, and step 66's top
# fifteen -- if the accessions mean what I suspect -- point hard at the second.
# The chapter currently argues about H2 and acetate. If the reconstruction is
# trading amino acids instead, that is the finding, and it changes the text and
# not just the picture.
#
# READ-ONLY except 67_network/. Seconds. No environment touched, no network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/67_network"; mkdir -p "$OUT"
exec > >(tee "$OUT/67.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "67_network_compounds -- naming the traded metabolites and orienting the edges"
echo "WS2 = $WS2"
echo "started $(date -Is)"

DET="$WS2/43_smetana/run_43q/pairs/pairs_detailed.tsv"
ALT="$WS2/43_smetana/run_43p/pairs/pairs_detailed.tsv"
MASTER="$WS2/47_master/MASTER.tsv"
[ -s "$DET" ] || { no "missing $DET"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. the names
say "A. WHERE THE COMPOUND NAMES COME FROM"
sub "SBML models available"
MDIR=""
for d in "$WS2/43_smetana/run_43o/models_bounded" \
         "$WS2/43_smetana/run_43q/pairs" \
         "$WS2/43_smetana/models_sed" "$WS2/43_smetana/models"; do
  N=$(find "$d" -maxdepth 1 -name '*.xml' 2>/dev/null | wc -l)
  printf '    %-46s %2s model(s)\n' "${d#$WS2/}" "$N"
  [ -z "$MDIR" ] && [ "$N" -gt 0 ] && MDIR="$d"
done
[ -n "$MDIR" ] && ok "reading names from ${MDIR#$WS2/}" || no "no SBML found -- names will be blank"

python3 - "$DET" "${MDIR:-/nonexistent}" "$OUT/compound_names.tsv" <<'PYEOF'
import sys, os, re, csv, glob, collections
det, mdir, outmap = sys.argv[1], sys.argv[2], sys.argv[3]

rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
cC, cS = low["compound"], low["smetana"]


def bare(cid):
    """M_cpd00116_e -> cpd00116. The join step 66 failed to do."""
    m = re.search(r"(cpd\d+)", cid or "")
    return m.group(1) if m else (cid or "")


ids = sorted({r[cC] for r in rows})
print("    %d distinct compound id(s) in the transfer table" % len(ids))
print("    example, as stored: %s   -> bare id %s" % (ids[0], bare(ids[0])))

name = {}
nfile = 0
for p in sorted(glob.glob(os.path.join(mdir, "*.xml"))):
    nfile += 1
    try:
        txt = open(p, encoding="utf-8", errors="replace").read()
    except OSError:
        continue
    for sid, nm in re.findall(r'<species[^>]*?id="([^"]+)"[^>]*?name="([^"]+)"', txt):
        b = bare(sid)
        if b.startswith("cpd"):
            nm = re.sub(r"_[ce]\d*$", "", nm).replace("_", " ").strip()
            name.setdefault(b, nm)
    for nm, sid in re.findall(r'<species[^>]*?name="([^"]+)"[^>]*?id="([^"]+)"', txt):
        b = bare(sid)
        if b.startswith("cpd"):
            nm = re.sub(r"_[ce]\d*$", "", nm).replace("_", " ").strip()
            name.setdefault(b, nm)
print("    parsed %d SBML file(s); %d compound name(s) learned" % (nfile, len(name)))

hit = [i for i in ids if bare(i) in name]
print("    %d of %d transfer compounds resolve to a name" % (len(hit), len(ids)))
missing = [i for i in ids if bare(i) not in name]
if missing:
    print("    unresolved: %s" % ", ".join(missing[:12]))

with open(outmap, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["compound_id", "cpd", "name"])
    for i in ids:
        w.writerow([i, bare(i), name.get(bare(i), "")])
print("    wrote %s" % outmap)
PYEOF

# ---------------------------------------------------------------- B. corrected
say "B. THE CORRECTED TEST -- ARE THE ENERGETIC SUBSTRATES THERE"
python3 - "$DET" "$OUT/compound_names.tsv" <<'PYEOF'
import sys, csv, re, collections
det, mapf = sys.argv[1], sys.argv[2]
nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(mapf), delimiter="\t")}
rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
cC, cS, cD, cR = low["compound"], low["smetana"], low["donor"], low["receiver"]


def bare(c):
    m = re.search(r"(cpd\d+)", c or "")
    return m.group(1) if m else c


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


# the five the chapter argues about, by ModelSEED accession
WANT = [("cpd11640", "H2"), ("cpd00047", "formate"), ("cpd00029", "acetate"),
        ("cpd00116", "methanol"), ("cpd00441", "methylamine"),
        ("cpd00011", "CO2"), ("cpd00239", "H2S")]
print("    tested on the BARE id, which is the join step 66 got wrong:")
tot = collections.Counter(); cnt = collections.Counter()
for r in rows:
    b = bare(r[cC]); tot[b] += f(r[cS]); cnt[b] += 1
for cpd, label in WANT:
    if cnt[cpd]:
        print("      %-10s %-12s %-26s sum %7.3f   in %3d transfer(s)"
              % (cpd, label, nm.get(cpd, "")[:26], tot[cpd], cnt[cpd]))
    else:
        print("      %-10s %-12s absent from the table" % (cpd, label))
print()
grand = sum(tot.values())
ener = sum(tot[c] for c, _ in WANT)
print("    energetic + one-carbon substrates carry %.3f of %.3f total score (%.1f%%)"
      % (ener, grand, 100.0 * ener / grand if grand else 0.0))
PYEOF

# ---------------------------------------------------------------- C. classes
say "C. WHAT IS ACTUALLY BEING TRADED"
python3 - "$DET" "$OUT/compound_names.tsv" "$OUT/compound_summary.tsv" <<'PYEOF'
import sys, csv, re, collections
det, mapf, outf = sys.argv[1], sys.argv[2], sys.argv[3]
nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(mapf), delimiter="\t")}
rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
cC, cS = low["compound"], low["smetana"]


def bare(c):
    m = re.search(r"(cpd\d+)", c or "")
    return m.group(1) if m else c


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


AA = ("alanine glycine serine threonine cysteine valine leucine isoleucine "
      "proline phenylalanine tyrosine tryptophan aspart asparagine glutam "
      "lysine arginine histidine methionine ornithine citrulline").split()
VIT = ("thiamin riboflavin niacin nicotin pantothen pyridox biotin folate folic "
       "cobalamin cobinamide b12 menaquinone ubiquinone heme siroheme lipoate "
       "queuosine molybdopterin coa tetrahydrofolate").split()
NUC = ("adenine guanine cytosine uracil thymine adenosine guanosine cytidine "
       "uridine thymidine inosine xanthine hypoxanthine ribose deoxy").split()
ENE = ("h2 hydrogen formate acetate ethanol lactate succinate propionate "
       "butyrate pyruvate methanol methylamine co2 sulfide").split()


def klass(cpd, n):
    s = (n or "").lower()
    if cpd in ("cpd11640", "cpd00047", "cpd00029", "cpd00116", "cpd00441",
               "cpd00011", "cpd00239"):
        return "energy/C1"
    for k in ENE:
        if k in s: return "energy/C1"
    for k in AA:
        if k in s: return "amino acid"
    for k in VIT:
        if k in s: return "vitamin/cofactor"
    for k in NUC:
        if k in s: return "nucleotide/base"
    return "other" if s else "unnamed"


tot = collections.Counter(); cnt = collections.Counter()
for r in rows:
    b = bare(r[cC]); tot[b] += f(r[cS]); cnt[b] += 1

cls_t = collections.Counter(); cls_c = collections.Counter()
with open(outf, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["cpd", "name", "class", "sum_smetana", "n_transfers"])
    for b, s in tot.most_common():
        k = klass(b, nm.get(b, ""))
        cls_t[k] += s; cls_c[k] += cnt[b]
        w.writerow([b, nm.get(b, ""), k, "%.4f" % s, cnt[b]])

grand = sum(tot.values())
print("    BY CLASS -- the classifier is a name-substring rule, so every compound")
print("    is listed below and the call can be audited line by line.")
print()
print("      %-18s %9s %8s  %8s" % ("class", "sum", "share", "transfers"))
for k, s in cls_t.most_common():
    print("      %-18s %9.3f %7.1f%%  %8d" % (k, s, 100.0 * s / grand if grand else 0, cls_c[k]))
print()
print("    EVERY COMPOUND, most-transferred first:")
print("      %-10s %-34s %-17s %8s %6s" % ("cpd", "name", "class", "sum", "n"))
for b, s in tot.most_common():
    print("      %-10s %-34s %-17s %8.3f %6d"
          % (b, (nm.get(b, "") or "(no name in SBML)")[:34], klass(b, nm.get(b, "")), s, cnt[b]))
print()
print("    wrote %s" % outf)
PYEOF

# ---------------------------------------------------------------- D. direction
say "D. THE DIRECTED MATRIX -- WHO FEEDS WHOM"
python3 - "$DET" "$MASTER" "$OUT/fig06_matrix.tsv" <<'PYEOF'
import sys, csv, re, os, collections
det, master, outf = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
cD, cR, cC, cS = low["donor"], low["receiver"], low["compound"], low["smetana"]

meta = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            meta[r["working_name"]] = r


def tag(b):
    r = meta.get(b)
    if not r: return b
    return (r.get("mag_id", b) or b).replace("ETO_", "")


def dom(b):
    r = meta.get(b)
    return "A" if r and str(r.get("domain", "")).lower().startswith("arch") else "B"


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


names = sorted({r[cD] for r in rows} | {r[cR] for r in rows}, key=lambda b: (dom(b), tag(b)))
S = collections.defaultdict(float)
N = collections.Counter()
for r in rows:
    S[(r[cD], r[cR])] += f(r[cS])
    N[(r[cD], r[cR])] += 1

print("    rows DONATE, columns RECEIVE. Cell = summed SMETANA score.")
print("    (A) marks an archaeon. A dash is no predicted transfer at all.")
print()
hdr = "      %-20s" % "donor \\ receiver"
for b in names:
    hdr += "%13s" % tag(b)
print(hdr)
for a in names:
    line = "      %-20s" % (tag(a) + (" (A)" if dom(a) == "A" else ""))
    for b in names:
        v = S.get((a, b), 0.0)
        line += "%13s" % ("-" if (a == b or v <= 0) else "%.2f" % v)
    print(line)

print()
print("    the same matrix as NUMBER OF COMPOUNDS moved:")
print(hdr)
for a in names:
    line = "      %-20s" % (tag(a) + (" (A)" if dom(a) == "A" else ""))
    for b in names:
        v = N.get((a, b), 0)
        line += "%13s" % ("-" if (a == b or v <= 0) else str(v))
    print(line)

print()
print("    TOTALS")
print("      %-20s %9s %9s %8s %8s" % ("genome", "donated", "received", "partners", "dom"))
for b in names:
    d = sum(v for (x, y), v in S.items() if x == b)
    q = sum(v for (x, y), v in S.items() if y == b)
    p = len({y for (x, y) in S if x == b} | {x for (x, y) in S if y == b})
    print("      %-20s %9.2f %9.2f %8d %8s" % (tag(b), d, q, p, dom(b)))

src = [tag(b) for b in names if sum(v for (x, y), v in S.items() if y == b) == 0]
snk = [tag(b) for b in names if sum(v for (x, y), v in S.items() if x == b) == 0]
print()
print("    receives nothing at all: %s" % (", ".join(src) or "none"))
print("    donates nothing at all : %s" % (", ".join(snk) or "none"))

with open(outf, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["donor_mag", "donor_domain", "receiver_mag", "receiver_domain",
                "sum_smetana", "n_compounds"])
    for (a, b), v in sorted(S.items(), key=lambda kv: -kv[1]):
        w.writerow([tag(a), dom(a), tag(b), dom(b), "%.4f" % v, N[(a, b)]])
print("    wrote %s  (%d directed edge(s))" % (outf, len(S)))
PYEOF

# ---------------------------------------------------------------- E. archaea
say "E. THE TWO ARCHAEA, IN FULL"
python3 - "$DET" "$MASTER" "$OUT/compound_names.tsv" <<'PYEOF'
import sys, csv, re, os, collections
det, master, mapf = sys.argv[1], sys.argv[2], sys.argv[3]
nm = {r["cpd"]: r["name"] for r in csv.DictReader(open(mapf), delimiter="\t")}
rows = list(csv.DictReader(open(det), delimiter="\t"))
low = {c.lower(): c for c in rows[0]}
cD, cR, cC, cS = low["donor"], low["receiver"], low["compound"], low["smetana"]
meta = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            meta[r["working_name"]] = r
arc = [b for b in meta if str(meta[b].get("domain", "")).lower().startswith("arch")]
inrun = {r[cD] for r in rows} | {r[cR] for r in rows}
arc = [b for b in arc if b in inrun]


def bare(c):
    m = re.search(r"(cpd\d+)", c or "")
    return m.group(1) if m else c


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


def tag(b):
    return (meta.get(b, {}).get("mag_id", b) or b).replace("ETO_", "")


for b in arc:
    g = meta.get(b, {}).get("genus", "-")
    print("\n    %s (%s)  [%s]" % (tag(b), g, b))
    rec = [r for r in rows if r[cR] == b and f(r[cS]) > 0]
    don = [r for r in rows if r[cD] == b and f(r[cS]) > 0]
    print("      RECEIVES %d transfer(s):" % len(rec))
    agg = collections.defaultdict(float)
    for r in rec:
        agg[(bare(r[cC]), r[cD])] += f(r[cS])
    for (c, d), v in sorted(agg.items(), key=lambda kv: -kv[1])[:20]:
        print("        %-10s %-30s from %-22s %.3f"
              % (c, (nm.get(c, "") or "?")[:30], tag(d), v))
    if not rec:
        print("        nothing. On this medium the model predicts it is fed by no one.")
    print("      DONATES %d transfer(s); top compounds:" % len(don))
    agg = collections.defaultdict(float)
    for r in don:
        agg[bare(r[cC])] += f(r[cS])
    for c, v in sorted(agg.items(), key=lambda kv: -kv[1])[:12]:
        print("        %-10s %-30s %.3f" % (c, (nm.get(c, "") or "?")[:30], v))
PYEOF

# ---------------------------------------------------------------- F. 43p vs 43q
say "F. 43p AND 43q BOTH RAN PAIRWISE -- ARE THEY THE SAME RUN"
if [ -s "$ALT" ]; then
  python3 - "$DET" "$ALT" <<'PYEOF'
import sys, csv, re
a, b = sys.argv[1], sys.argv[2]


def load(p):
    rs = list(csv.DictReader(open(p), delimiter="\t"))
    low = {c.lower(): c for c in rs[0]}
    g = {r[low["donor"]] for r in rs} | {r[low["receiver"]] for r in rs}
    med = {r[low["medium"]] for r in rs} if "medium" in low else set()
    return rs, g, med, low


ra, ga, ma, la = load(a)
rb, gb, mb, lb = load(b)
print("    43q: %5d row(s), %2d genome(s), medium %s" % (len(ra), len(ga), sorted(ma)))
print("    43p: %5d row(s), %2d genome(s), medium %s" % (len(rb), len(gb), sorted(mb)))
print("    genomes only in 43q: %s" % (", ".join(sorted(ga - gb)) or "none"))
print("    genomes only in 43p: %s" % (", ".join(sorted(gb - ga)) or "none"))
ca = {r[la["compound"]] for r in ra}
cb = {r[lb["compound"]] for r in rb}
print("    compounds: 43q %d, 43p %d, shared %d" % (len(ca), len(cb), len(ca & cb)))
print()
print("    43q is the run this step and the figure use. It is the larger of the two")
print("    and the one 43q_real_test built on the bounded models. If the genome sets")
print("    differ above, the two are not alternative views of one experiment and the")
print("    figure must name which it shows.")
PYEOF
else
  no "no 43p pairwise output to compare"
fi

# ---------------------------------------------------------------- G. next
say "G. WHAT FIG06 BECOMES"
cat <<'EOT'
  Section B decides the panel's subject.

    energetic substrates carry a large share  -> the syntrophy figure the chapter
                                                 currently promises: H2 and
                                                 formate edges, methanogens as
                                                 receivers
    they carry almost none                    -> the figure is about amino acid
                                                 and cofactor exchange, and the
                                                 chapter text has to change with
                                                 it. That is not a worse result.
                                                 A sulfate-replete lagoon where
                                                 the model predicts no H2
                                                 syntrophy is consistent with
                                                 every other line in this thesis:
                                                 the sulfate reducers take the H2

  Section D decides the form. 11 nodes and 32 directed edges is at the boundary
  where a graph stops being readable, and the matrix printed there is already a
  legible figure. The choice between a directed graph and a matrix follows what that
  matrix looks like, not before.

  Section E is the one that matters most for the argument: if the two methanogens
  receive nothing, then the reconstruction says they are not being fed by this
  community, and the methylotrophic route in fig04 is not a preference but the
  only door left open.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
