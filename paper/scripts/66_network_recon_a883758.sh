#!/usr/bin/env bash
# 66_network_recon.sh -- what is actually in the SMETANA output, and a figure
#                        input table built from it.
#
#   bash scripts/run.sh scripts/66_network_recon_XXXXXXX.sh
#
# WHY A RECON STEP AND NOT THE FIGURE. The interaction network is the one
# remaining panel, and it is the only one whose input I have never seen. I know
# from 43q's own code that the pairwise run wrote pairs_detailed.tsv and that
# the awk in section D reads $3 as receiver, $4 as donor, $5 as compound and $9
# as the SMETANA score. That is the script's belief about its own output. It is
# not a measurement of the file, and the difference between those two things is
# what produced Volume XV's 280-fold mtmB discrepancy.
#
#   > Rule 106. Before running an analysis for the first time, check that it has
#   > not already been run twice.
#
# So this step measures, and it settles four things the figure cannot be written
# without:
#
#   1. WHICH FILE. 43p and 43q both ran pairwise. If both left output, the figure
#      must be told which, and the two must be compared rather than assumed
#      identical.
#   2. HOW MANY EDGES SURVIVE. A network panel is only worth drawing if there is
#      a network. If the non-zero rows are three transfers between two genomes,
#      the honest figure is a table.
#   3. THE NAMES. SMETANA ran on bin names (LIQ_das_semibin_003). The thesis
#      speaks MAG ids (ETO_LIQ_MAG_003). MASTER.tsv carries both, in working_name
#      and mag_id. Every bin in the community set has to resolve, and any that
#      does not has to be named here rather than silently dropped from a panel.
#   4. WHICH GENOMES WERE IN IT AT ALL. The pairwise run used the bounded models
#      in run_43o, which is ELEVEN genomes, not 29. A network panel that does not
#      say so implies the whole catalogue was tested.
#
# It also writes 66_network/fig06_input.tsv -- one row per transfer, already
# carrying mag_id, genus and domain for both partners -- so the figure script
# reads a table with known columns instead of re-deriving the join. If the join
# is wrong, it is wrong HERE, in a file you can look at, not inside a renderer.
#
# READ-ONLY except for 66_network/. Runs in seconds. Touches no environment.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/66_network"; mkdir -p "$OUT"
exec > >(tee "$OUT/66.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "66_network_recon -- SMETANA output structure and the figure input table"
echo "WS2 = $WS2"
echo "started $(date -Is)"

MASTER="$WS2/47_master/MASTER.tsv"

# ---------------------------------------------------------------- A. find them
say "A. EVERY SMETANA OUTPUT IN THE WORKSPACE"
sub "anything named *detailed* or *smetana* under 43_smetana"
FOUND=0
while read -r f; do
  [ -z "$f" ] && continue
  FOUND=$((FOUND+1))
  printf '    %-58s %9s bytes  %7s lines  md5 %s\n' \
    "${f#$WS2/}" "$(stat -c %s "$f")" "$(wc -l < "$f")" "$(md5sum "$f" | cut -c1-8)"
done < <(find "$WS2/43_smetana" -type f \( -name '*detailed*' -o -name '*smetana*' \) 2>/dev/null | sort)
[ "$FOUND" -eq 0 ] && no "none found -- the pairwise run left nothing, and that is the finding"
echo "    $FOUND file(s)"

sub "the run directories that exist"
for d in "$WS2"/43_smetana/*/; do
  [ -d "$d" ] || continue
  printf '    %-34s %s file(s)\n' "${d#$WS2/}" "$(find "$d" -type f | wc -l)"
done

# choose the file the figure will use, but PRINT the choice and the alternatives
DET=""
for c in "$WS2/43_smetana/run_43q/pairs/pairs_detailed.tsv" \
         "$WS2/43_smetana/run_43q/pairs_detailed.tsv" \
         "$WS2/43_smetana/run_43p/pairs/pairs_detailed.tsv"; do
  [ -s "$c" ] && { DET="$c"; break; }
done
if [ -z "$DET" ]; then
  DET=$(find "$WS2/43_smetana" -type f -name '*detailed*' -size +1c 2>/dev/null | head -1)
fi
sub "the file this step will read"
if [ -n "$DET" ] && [ -s "$DET" ]; then
  ok "${DET#$WS2/}"
else
  no "NO non-empty detailed output anywhere. Sections B-E cannot run."
  echo
  echo "  This is not a failure of this script. It means the pairwise SMETANA run"
  echo "  left no interaction table, and the interaction-network panel has no"
  echo "  input. Say so and the figure is replaced, not debugged."
  echo
  echo "STEP_STATUS: NO_INPUT $(date -Is)"
  exit 0
fi

# ---------------------------------------------------------------- B. structure
say "B. WHAT THE COLUMNS ACTUALLY ARE"
sub "header, verbatim, one column per line with its index"
head -1 "$DET" | tr '\t' '\n' | nl -ba | sed 's/^/    /'
NC=$(head -1 "$DET" | awk -F'\t' '{print NF}')
NR_=$(($(wc -l < "$DET")-1))
echo
echo "    $NC columns, $NR_ data row(s)"

sub "first 8 data rows, verbatim"
sed -n '2,9p' "$DET" | sed 's/^/    /'

sub "43q's own awk read \$3 receiver, \$4 donor, \$5 compound, \$9 smetana"
echo "    column 3 header: $(head -1 "$DET" | cut -f3)"
echo "    column 4 header: $(head -1 "$DET" | cut -f4)"
echo "    column 5 header: $(head -1 "$DET" | cut -f5)"
echo "    column 9 header: $(head -1 "$DET" | cut -f9)"
echo
echo "    If those four names are not receiver / donor / compound / smetana, then"
echo "    every number section D of 43q printed was read out of the wrong column,"
echo "    and Volume XIII's interaction findings need re-checking before anything"
echo "    is drawn."

# ---------------------------------------------------------------- C. the edges
say "C. IS THERE A NETWORK TO DRAW"
python3 - "$DET" <<'PYEOF'
import sys, csv, collections
p = sys.argv[1]
rows = list(csv.DictReader(open(p), delimiter="\t"))
if not rows:
    print("    no data rows"); sys.exit(0)
cols = list(rows[0])


def pick(*names):
    low = {c.lower(): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


cR, cD, cC = pick("receiver"), pick("donor"), pick("compound")
cS = pick("smetana")
print("    resolved by NAME (not by position): receiver=%s donor=%s compound=%s score=%s"
      % (cR, cD, cC, cS))
if not all((cR, cD, cC, cS)):
    print("    !! a required column is missing by name. Stopping before guessing.")
    sys.exit(0)


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


pos = [r for r in rows if f(r[cS]) > 0]
print("    %d row(s) total, %d with smetana > 0" % (len(rows), len(pos)))
if not pos:
    print("    >>> EVERY score is zero. There is no network. The panel becomes a")
    print("        statement, not a graph: the reconstruction predicts no transfer.")
    sys.exit(0)

sc = sorted(f(r[cS]) for r in pos)
q = lambda t: sc[min(len(sc) - 1, int(t * (len(sc) - 1)))]
print("    score range %.4g .. %.4g   median %.4g   p90 %.4g" % (sc[0], sc[-1], q(.5), q(.9)))

genomes = sorted({r[cD] for r in pos} | {r[cR] for r in pos})
print("    %d genome(s) appear in a non-zero transfer:" % len(genomes))
for g in genomes:
    d = sum(1 for r in pos if r[cD] == g)
    x = sum(1 for r in pos if r[cR] == g)
    print("      %-28s donates %3d   receives %3d" % (g, d, x))

pair = collections.Counter()
for r in pos:
    pair[(r[cD], r[cR])] += 1
print("    %d ordered donor->receiver pair(s) carry at least one transfer" % len(pair))

cmp_ = collections.Counter()
tot = collections.Counter()
for r in pos:
    cmp_[r[cC]] += 1
    tot[r[cC]] += f(r[cS])
print("    %d distinct compound(s); the 15 most transferred by summed score:" % len(cmp_))
for c, s in tot.most_common(15):
    print("      %-18s sum %8.3f   in %3d row(s)" % (c, s, cmp_[c]))

KNOWN = {"cpd11640": "H2", "cpd00047": "formate", "cpd00029": "acetate",
         "cpd00011": "CO2", "cpd00067": "H+", "cpd00001": "H2O",
         "cpd00013": "NH3", "cpd00009": "phosphate", "cpd00048": "sulfate",
         "cpd00239": "H2S", "cpd00116": "methanol", "cpd00441": "methylamine"}
print("    of those, the ones this chapter argues about:")
hit = False
for c, s in tot.most_common():
    if c in KNOWN:
        hit = True
        print("      %-18s %-10s sum %8.3f   in %3d row(s)" % (c, KNOWN[c], s, cmp_[c]))
if not hit:
    print("      none of H2, formate, acetate, methanol or methylamine appears at all")
PYEOF

# ---------------------------------------------------------------- D. the names
say "D. DO THE SMETANA NAMES RESOLVE TO MAG IDS"
if [ ! -s "$MASTER" ]; then
  no "MASTER.tsv not at ${MASTER#$WS2/} -- the join cannot be checked"
else
  ok "MASTER.tsv: $(($(wc -l < "$MASTER")-1)) row(s)"
  python3 - "$DET" "$MASTER" <<'PYEOF'
import sys, csv
det, master = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(det), delimiter="\t"))
cols = list(rows[0]) if rows else []
low = {c.lower(): c for c in cols}
cD, cR = low.get("donor"), low.get("receiver")
if not (cD and cR):
    print("    donor/receiver columns not found by name; skipping the join")
    sys.exit(0)
names = sorted({r[cD] for r in rows} | {r[cR] for r in rows})

m = list(csv.DictReader(open(master), delimiter="\t"))
mcols = list(m[0]) if m else []
print("    MASTER.tsv columns: %s" % ", ".join(mcols))
wn = "working_name" if "working_name" in mcols else None
if not wn:
    print("    !! MASTER.tsv has no working_name column -- nothing to join on")
    sys.exit(0)
bin2 = {r[wn]: r for r in m if r.get(wn)}

print("    %d distinct genome name(s) in the SMETANA output:" % len(names))
miss = []
for n in names:
    r = bin2.get(n)
    if r:
        print("      %-28s -> %-18s %-22s %s"
              % (n, r.get("mag_id", "?"), (r.get("genus") or "-")[:22],
                 r.get("domain", "?")[:8]))
    else:
        miss.append(n)
        print("      %-28s -> NO MATCH in MASTER.tsv" % n)
if miss:
    print("    !! %d name(s) do not resolve. A panel drawn now would drop them" % len(miss))
    print("       silently. Fix the join before the figure, not inside it.")
else:
    print("    all %d resolve." % len(names))

arc = [bin2[n] for n in names if n in bin2
       and str(bin2[n].get("domain", "")).lower().startswith("arch")]
print("    %d of the %d are archaea -- %s"
      % (len(arc), len(names),
         ", ".join(a.get("mag_id", "?") for a in arc) or "none"))
print("    THE CATALOGUE HAS 29 GENOMES. This run tested %d. Whatever the panel"
      % len(names))
print("    shows, its caption has to say which %d and why the other %d are absent."
      % (len(names), 29 - len(names)))
PYEOF
fi

# ---------------------------------------------------------------- E. the table
say "E. THE FIGURE INPUT TABLE"
python3 - "$DET" "$MASTER" "$OUT/fig06_input.tsv" <<'PYEOF'
import sys, csv, os
det, master, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(det), delimiter="\t"))
if not rows:
    print("    nothing to write"); sys.exit(0)
low = {c.lower(): c for c in rows[0]}
cD, cR, cC, cS = (low.get("donor"), low.get("receiver"),
                  low.get("compound"), low.get("smetana"))
if not all((cD, cR, cC, cS)):
    print("    required columns missing by name; no table written"); sys.exit(0)

meta = {}
if os.path.exists(master):
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("working_name"):
            meta[r["working_name"]] = r


def g(n, k, d="?"):
    r = meta.get(n)
    return (r.get(k) or d) if r else d


def f(v):
    try: return float(v)
    except (TypeError, ValueError): return 0.0


hdr = ["donor_bin", "donor_mag", "donor_genus", "donor_domain",
       "receiver_bin", "receiver_mag", "receiver_genus", "receiver_domain",
       "compound", "smetana"]
n = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(hdr)
    for r in rows:
        s = f(r[cS])
        if s <= 0:
            continue
        d_, r_ = r[cD], r[cR]
        w.writerow([d_, g(d_, "mag_id"), g(d_, "genus", "-"), g(d_, "domain"),
                    r_, g(r_, "mag_id"), g(r_, "genus", "-"), g(r_, "domain"),
                    r[cC], "%.6g" % s])
        n += 1
print("    wrote %s" % out)
print("    %d non-zero transfer row(s), %d column(s)" % (n, len(hdr)))
PYEOF

sub "fig06_input.tsv, first 20 rows"
if [ -s "$OUT/fig06_input.tsv" ]; then
  if command -v column >/dev/null 2>&1; then
    head -21 "$OUT/fig06_input.tsv" | column -t -s $'\t' | sed 's/^/    /'
  else
    head -21 "$OUT/fig06_input.tsv" | awk -F'\t' '{for(i=1;i<=NF;i++)printf "%-20s",$i; print ""}' | sed 's/^/    /'
  fi
  echo
  echo "    $(($(wc -l < "$OUT/fig06_input.tsv")-1)) row(s) total, md5 $(md5sum "$OUT/fig06_input.tsv" | cut -c1-8)"
else
  no "no table written"
fi

# ---------------------------------------------------------------- F. next
say "F. WHAT THE ANSWER DECIDES"
cat <<'EOT'
  Section C decides whether fig06 is a network at all.

    many edges, several compounds  -> the panel is a directed graph: genomes as
                                      nodes, one edge per donor->receiver pair,
                                      edge weight the summed SMETANA score, and
                                      the compounds that matter coloured
    a handful of edges             -> a graph with six lines on it is a table
                                      pretending to be a figure. It becomes a
                                      small matrix or a labelled list instead
    every score zero               -> there is no panel. The result is that the
                                      reconstruction predicts no cross-feeding
                                      on the sediment medium, which is a finding
                                      and belongs in the text, not in a diagram

  Section D decides whether the nodes can be named at all, and forces the caption
  to state how many genomes were tested rather than implying all 29 were.

  Nothing is drawn until this output is read. The last four panels were written
  against real column names and real ranges, and that is the only reason they
  rendered first time on MASTER.tsv.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
