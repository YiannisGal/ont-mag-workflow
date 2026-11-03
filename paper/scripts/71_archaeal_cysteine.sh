#!/usr/bin/env bash
# 71_archaeal_cysteine.sh -- is MTP4's Cys 0,00 an auxotrophy or a marker-set failure?
#
#   bash scripts/run.sh scripts/71_archaeal_cysteine_XXXXXXX.sh
#
# WHERE STEP 70 LEFT IT. Three of the four explanations for fig06's empty columns
# are dead:
#
#   the medium supplies it   -> NO. None of the 64 traded compounds appears among
#                               ETO_SED's 24, so the medium explains nothing
#   a missing exchange       -> NO. All eleven models carry EX_ for cysteine and
#                               for every other top compound bar pantothenate
#   a smaller biomass        -> NO, not for this. Both archaea ask for 58
#                               precursors and all nine bacteria ask for 85, so
#                               the objective cannot explain why one archaeon
#                               receives from everyone and the other from no one
#
# What remains is that MTP4's model really does make everything it needs on
# ETO_SED. Its genome is 99,84% complete -- the most complete in the set -- and
# its mean amino-acid pathway completeness, 0,642, is the highest of the archaea.
# That is consistent.
#
# EXCEPT CYSTEINE. 19b gives MTP4 Cys = 0,00: no route at all. Cysteine is the
# single most-transferred compound in the whole network, 15,55. A genome with no
# cysteine route whose model needs no cysteine is a contradiction, and there are
# only two ways out:
#
#   (a) gap-filling handed the model a cysteine route the genome does not have
#   (b) 19b's Cys marker set cannot see the ARCHAEAL route
#
# (b) is not a guess. ARCHAEAL_ANNOTATION_GAP.md records three independent
# instances already:
#
#   > Rule 66. A tool trained on one domain, applied to another, fails silently
#   > and in one direction.
#
# Bakta called 240 of 253 KOfam-confirmed archaeal markers "hypothetical" -- 95%,
# including mcrA at score 1000,3. DeepSig's Gram-positive model under-called the
# archaeal secretome. Both archaea here sit low on Cys (0,00 and 0,33). That is
# the signature of Rule 66, not of two unrelated auxotrophies.
#
# AND THE THESIS CARRIES A STANDING REQUIREMENT OF EQUAL TREATMENT OF BACTERIA AND
# ARCHAEA. If 19b's amino-acid routes are bacterial-canonical, every archaeal
# auxotrophy call in this project inherits the same one-sided error, and fig06 is
# the least of it.
#
# THE TESTS, none of which assumes a KO accession I have not read off this machine:
#
#   A. Is the deficit SYSTEMATIC? All 29 genomes, every 19b pathway, archaea
#      against bacteria. One low pathway in one genome is biology; a one-sided
#      deficit across the board is the tool.
#   B. WHICH KOs does 19b actually use for Cys? Read from routes.tsv, not recalled.
#   C. Does MTP4 carry ANY cysteine-related KO? Found by reading ko_list's own
#      definitions for "cystein" and intersecting with the KOfam hits already on
#      disk -- above each KO's curated threshold, per rule 97.
#   D. HOW MUCH GAP-FILLING did MTP4's model get? Step 68 compared against
#      models / models_sed / models_fillrich and MTP4 was in none of them, so it
#      was never compared. Step 70 shows it IS in run_43l, 43m and 43n. Those are
#      the comparisons that were available and that I failed to make.
#
# READ-ONLY except 71_cys/. Seconds. No network.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/71_cys"; mkdir -p "$OUT"
exec > >(tee "$OUT/71.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "71_archaeal_cysteine -- auxotrophy, or Rule 66 again?"
echo "WS2 = $WS2"
echo "started $(date -Is)"

PC="$WS2/19b_auxotrophy/pathway_completeness.tsv"
ROUTES="$WS2/19b_auxotrophy/routes.tsv"
KOLIST="$WS2/36_gapclose/ko_list"
MASTER="$WS2/47_master/MASTER.tsv"

# ---------------------------------------------------------------- A. systematic
say "A. IS THE ARCHAEAL DEFICIT SYSTEMATIC -- all 29 genomes, every pathway"
python3 - "$PC" <<'PYEOF'
import sys, csv, os
pcf = sys.argv[1]
if not os.path.exists(pcf):
    print("    missing %s" % pcf); sys.exit(0)
R = list(csv.DictReader(open(pcf), delimiter="\t"))
cols = [c for c in R[0] if c not in
        ("mag_id", "working_name", "sample", "domain", "genus", "completeness")]
arc = [r for r in R if str(r.get("domain", "")).lower().startswith("arch")]
bac = [r for r in R if not str(r.get("domain", "")).lower().startswith("arch")]
print("    %d archaea, %d bacteria in 19b's table" % (len(arc), len(bac)))
if not arc or not bac:
    sys.exit(0)


def mean(rs, c):
    v = []
    for r in rs:
        try: v.append(float(r[c]))
        except (TypeError, ValueError, KeyError): pass
    return (sum(v) / len(v)) if v else float("nan")


print()
print("      %-18s %9s %9s %9s" % ("pathway", "archaea", "bacteria", "diff"))
lower = 0
rows = []
for c in cols:
    a, b = mean(arc, c), mean(bac, c)
    rows.append((c, a, b, a - b))
    if a < b:
        lower += 1
for c, a, b, d in sorted(rows, key=lambda x: x[3]):
    flag = "  <-- Cys" if c == "Cys" else ""
    print("      %-18s %9.3f %9.3f %+9.3f%s" % (c, a, b, d, flag))
print()
print("    archaea score lower on %d of %d pathways." % (lower, len(cols)))
if lower >= len(cols) * 0.75:
    print("    >>> THAT IS ONE-SIDED. A genuine biological difference would not put")
    print("        archaea below bacteria on nearly every pathway at once. This is")
    print("        the shape of Rule 66: a bacterial-canonical marker set applied to")
    print("        archaea, failing silently and in one direction.")
elif lower <= len(cols) * 0.4:
    print("    >>> not one-sided. The Cys result is specific and is not explained by")
    print("        a domain-wide marker failure.")
else:
    print("    >>> mixed. Section C is then the one that decides Cys specifically.")

# genome completeness is the obvious confound -- check it
def gc(rs):
    v = []
    for r in rs:
        try: v.append(float(r["completeness"]))
        except (TypeError, ValueError, KeyError): pass
    return sum(v) / len(v) if v else float("nan")


print()
print("    the obvious confound: mean genome completeness, archaea %.2f vs"
      " bacteria %.2f" % (gc(arc), gc(bac)))
print("    if the archaea are the MORE complete genomes and still score lower on")
print("    the pathways, the marker set is the only explanation left.")
PYEOF

# ---------------------------------------------------------------- B. the routes
say "B. WHICH KOs 19b USES FOR CYSTEINE -- read, not recalled"
if [ -s "$ROUTES" ]; then
  ok "${ROUTES#$WS2/}  ($(wc -l < "$ROUTES") line(s))"
  sub "every route row mentioning Cys"
  awk -F'\t' 'tolower($0) ~ /cys/' "$ROUTES" | sed 's/^/    /'
  sub "the full set of pathways routes.tsv defines"
  awk -F'\t' '$3!="CONTROL"{print $3}' "$ROUTES" | sort -u | tr '\n' ' ' | fold -w 72 -s | sed 's/^/    /'
  echo
else
  no "routes.tsv not at ${ROUTES#$WS2/} -- 19b writes it into its own output dir"
  find "$WS2/19b_auxotrophy" -maxdepth 1 -type f | sed 's/^/    /'
fi

# ---------------------------------------------------------------- C. the hits
say "C. DOES MTP4 CARRY ANY CYSTEINE KO AT ALL"
sub "KOfam tables available"
for d in "$WS2/11_kofam/hits" "$WS2/45_kofam_full" "$WS2/36_gapclose"; do
  printf '    %-40s %s file(s)\n' "${d#$WS2/}" "$(find "$d" -maxdepth 1 -type f 2>/dev/null | wc -l)"
done
if [ ! -s "$KOLIST" ]; then
  no "ko_list not at ${KOLIST#$WS2/} -- cannot resolve KO definitions locally"
else
  ok "ko_list: $(($(wc -l < "$KOLIST")-1)) KO(s) with curated thresholds"
  sub "every KO whose own definition matches the cysteine rule"
  echo "    rule: cystein | phosphoseryl | seryl-tRNA | the bare token cys"
  echo "    the last three are there because THE ARCHAEAL ROUTE IS tRNA-DEPENDENT"
  echo "    and its enzymes are not named after cysteine. Matching only the string"
  echo "    'cystein' would search for the bacterial route and call its absence an"
  echo "    auxotrophy -- rule 66 committed inside the test for rule 66."
  awk -F'\t' 'NR>1 && tolower($0) ~ /cystein|phosphoseryl|seryl-trna|[^a-z]cys[^a-z]/{printf "    %-10s thr %-10s %s\n", $1, $2, $NF}' \
    "$KOLIST" | head -40
  NCYS=$(awk -F'\t' 'NR>1 && tolower($0) ~ /cystein|phosphoseryl|seryl-trna|[^a-z]cys[^a-z]/' "$KOLIST" | wc -l)
  echo "    $NCYS cysteine-related KO(s) in ko_list"

  sub "which of them MTP4 hits, above each KO's own curated threshold"
  python3 - "$WS2" "$KOLIST" "$MASTER" <<'PYEOF'
import sys, os, csv, glob, re
ws2, kolist, master = sys.argv[1], sys.argv[2], sys.argv[3]

thr, defn = {}, {}
with open(kolist) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    ik = 0
    it = head.index("threshold") if "threshold" in head else 1
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) <= it:
            continue
        try: thr[p[ik]] = float(p[it])
        except ValueError: thr[p[ik]] = None
        defn[p[ik]] = p[-1]
CYSRE = re.compile(r"cystein|phosphoseryl|seryl-trna|(?<![a-z])cys(?![a-z])")
cys = {k for k, d in defn.items() if CYSRE.search((d or "").lower())}
print("    testing %d cysteine KO(s)" % len(cys))

# find MTP4's bin / mag names
mag, bin_ = "ETO_LIQ_MAG_003", None
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("mag_id") == mag:
            bin_ = r.get("working_name")
print("    MTP4 = %s / %s" % (mag, bin_))

cand = []
for pat in ("11_kofam/hits/*.tbl", "45_kofam_full/*.tbl", "11b_kofam_archaea/**/*.tbl"):
    cand += glob.glob(os.path.join(ws2, pat), recursive=True)
print("    %d KOfam table(s) visible" % len(cand))
for p in cand[:1]:
    print("    format of %s, first 3 data lines with field numbers:"
          % os.path.relpath(p, ws2))
    shown = 0
    for ln in open(p, errors="replace"):
        if ln.startswith("#") or not ln.strip():
            continue
        print("      " + " | ".join("%d:%s" % (i, x)
                                    for i, x in enumerate(ln.split()[:10])))
        shown += 1
        if shown >= 3:
            break

hits = {}
for p in cand:
    base = os.path.basename(p)
    whole = "all_hits" in base or "chunks" in p
    interesting = whole or (bin_ and bin_ in base) or ("MAG_003" in base)
    if not interesting:
        continue
    try:
        fh = open(p, errors="replace")
    except OSError:
        continue
    with fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split()
            if len(f) < 3:
                continue
            gene, ko = f[0], f[1] if f[1].startswith("K") else None
            if ko is None:
                for x in f[:4]:
                    if re.fullmatch(r"K\d{5}", x):
                        ko = x; break
            if ko is None or ko not in cys:
                continue
            # THE SCORE, NOT THE E-VALUE. A KOfam table carries both, and the
            # E-value comes first on the line -- "take the first float" would
            # compare 4,4e-302 against a threshold of 180 and reject every real
            # hit. Skip anything in scientific notation and take the largest
            # remaining plausible number.
            vals = []
            for x in f[1:]:
                if "e-" in x.lower() or "e+" in x.lower():
                    continue
                try:
                    v = float(x)
                except ValueError:
                    continue
                if 0.0 < v < 1e5:
                    vals.append(v)
            sc = max(vals) if vals else None
            if bin_ and not (bin_ in gene or bin_ in base or "MAG_003" in gene
                             or "MAG_003" in base):
                continue
            t = thr.get(ko)
            if sc is not None and t is not None and sc < t:
                continue
            hits.setdefault(ko, []).append((gene, sc, t, os.path.relpath(p, ws2)))

if not hits:
    print()
    print("    NO cysteine KO passes threshold for MTP4 in the tables searched.")
    print("    That does not settle it either: rule 56 says a KO absent from the")
    print("    searched profile set is not evidence of absence. Section B's list")
    print("    says which KOs were ever searched at all.")
else:
    print()
    for ko in sorted(hits):
        print("    %-10s %-46s %d hit(s)" % (ko, (defn.get(ko) or "")[:46], len(hits[ko])))
        for g, sc, t, src in hits[ko][:4]:
            print("        %-34s score %-10s thr %-8s %s" % (g, sc, t, src))
    print()
    print("    >>> MTP4 CARRIES A CYSTEINE KO ABOVE THRESHOLD. 19b's Cys = 0,00 is")
    print("        then a marker-set failure, not an auxotrophy, and the model and")
    print("        the genome agree after all.")
PYEOF
fi

# ---------------------------------------------------------------- D. gapfilling
say "D. MTP4'S GAP-FILLING, MEASURED AGAINST THE DRAFTS THAT DO EXIST"
python3 - "$WS2/43_smetana" "$MASTER" <<'PYEOF'
import sys, os, re, csv, glob, collections
root, master = sys.argv[1], sys.argv[2]
M = {}
for r in csv.DictReader(open(master), delimiter="\t"):
    if r.get("working_name"):
        M[r["working_name"]] = (r.get("mag_id", "") or "").replace("ETO_", "")


def nrxn(p):
    return len(re.findall(r"<reaction\b",
                          open(p, encoding="utf-8", errors="replace").read()))


dirs = sorted(d for d in (glob.glob(os.path.join(root, "*"))
                          + glob.glob(os.path.join(root, "*", "*")))
              if os.path.isdir(d) and glob.glob(os.path.join(d, "*.xml")))
tab = collections.defaultdict(dict)
for d in dirs:
    tag = os.path.relpath(d, root)
    for p in glob.glob(os.path.join(d, "*.xml")):
        b = os.path.basename(p).replace(".xml", "")
        tab[M.get(b, b)][tag] = nrxn(p)

for g in ("LIQ_MAG_003", "LIQ_MAG_011", "LIQ_MAG_025"):
    if g not in tab:
        print("    %s: no model anywhere" % g); continue
    print()
    print("    %s" % g)
    for tag in sorted(tab[g], key=lambda t: tab[g][t]):
        print("      %-34s %6d reaction(s)" % (tag, tab[g][tag]))
    v = sorted(tab[g].values())
    print("      spread %d .. %d  (%+d across every version on disk)"
          % (v[0], v[-1], v[-1] - v[0]))
print()
print("    A genome whose model is the SAME SIZE in every directory was not being")
print("    re-filled between runs. A wide spread means reactions were added, and")
print("    the bounded set used by 43q is the one at the end of that chain.")
PYEOF

# ---------------------------------------------------------------- E. medium
say "E. THE MEDIUM, READ WITH ITS OWN NAMES THIS TIME"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
if [ -s "$MDB" ]; then
  echo "    Section A of step 70 looked ETO_SED's compounds up in the 64-compound"
  echo "    map from the transfer table and reported '0 of 24 named'. That is a real"
  echo "    result -- none of the traded compounds is in the medium -- but it is NOT"
  echo "    a reading of the medium, because the file carries its own name column"
  echo "    and I never looked at it."
  sub "ETO_SED, every compound, with the file's own names"
  awk -F'\t' '$1=="ETO_SED"{printf "    %-14s %s\n", $3, $4}' "$MDB"
  sub "any amino acid among them"
  awk -F'\t' '$1=="ETO_SED" && tolower($4) ~ /cystein|methionin|arginin|aspart|glutam|valin|leucin|threonin|tryptophan|tyrosin|alanin|glycin|serin|prolin|histidin|lysin|phenylalanin|asparagin/{printf "    🔴 %-14s %s\n", $3, $4}' "$MDB"
  echo "    (nothing printed above = no free amino acid in ETO_SED, as designed)"
else
  no "media file missing"
fi

say "F. WHAT EACH ANSWER MEANS"
cat <<'EOT'
  A one-sided across nearly every pathway -> 19b's amino-acid calls are not safe
                                             for archaea, and that is a finding
                                             bigger than fig06: it touches every
                                             archaeal auxotrophy in the project,
                                             and the thesis's standing requirement
                                             of equal treatment says it has to be
                                             written down, not worked around

  C finds a cysteine KO above threshold   -> MTP4's Cys = 0,00 is a marker-set
                                             failure. Model and genome agree. The
                                             empty column is a result and fig06's
                                             MTP4 column can be discussed as
                                             biology -- while MAG_025's, at 52,5%
                                             genome completeness, cannot

  C finds none, and B shows the archaeal  -> absence of evidence (rule 56). The
  route was never in the profile set         honest line is that cysteine
                                             biosynthesis in MTP4 is untested

  D shows MTP4's model grew between runs  -> quantify it and put the number in the
                                             caption next to the empty column
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
