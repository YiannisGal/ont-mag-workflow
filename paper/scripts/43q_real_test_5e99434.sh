#!/usr/bin/env bash
# 43q_the_real_test.sh -- the medium that can actually show syntrophy.
#
# RUN IN tmux.
#
# 43p's RESULT, and why it could not have answered the question
#     MTP4 as RECEIVER across all 55 pairs :   0 rows
#     MTP4 as DONOR                        : 110 rows
#     H2 / formate / acetate, anywhere     :   0 rows
#     most transferred: L-Cysteine 13.82 over 25 transfers
#   Not one energy-currency transfer. Everything moving is amino acids and cofactors,
#   and the methanogen is a NET DONOR.
#
#   That is not evidence against syntrophy. It is a consequence of MY medium choice.
#   ETO_SED_FERM SUPPLIES H2, formate and acetate. SMETANA scores DEPENDENCY: if the
#   environment already provides a compound, nobody depends on a neighbour for it and
#   the score is correctly zero. In 43n I added those compounds so the methanogen would
#   have something to eat, and thereby removed the only reason it would ever need the
#   fermenters. Rule 63, and I walked into it.
#
# THE MEDIUM THAT ASKS THE QUESTION: ETO_SED.
#   Sugars present, so the fermenters grow and produce H2 and formate (43i measured
#   169 and 111 maximum secretion). No H2, no formate, no acetate supplied. MTP4's own
#   sugar uptake is blocked by the declared constraint, so it CANNOT grow there alone --
#   43n measured exactly 0.0000. Its only possible route to growth is a neighbour.
#
#   If MTP4 grows in a pair on ETO_SED, something was transferred, and section D names
#   it. If no pair works, section C says whether that is because the pair cannot grow at
#   all, which is a different statement from "no interaction".
#
# ALSO CORRECTED: 43p reported "largest working community: 4". n=5 TIMED OUT at the
# 15-minute cap I set; it did not fail. 468s at n=3, 685s at n=4. Only the eleven
# crashed outright. Section F re-probes n=5-7 with a real time budget.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43q"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp43q"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
BD="$WS2/43_smetana/run_43o/models_bounded"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017; MED=ETO_SED
PROBE="${PROBE:-0}"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43q  $(date -Is)"; echo "out: $OUT"; echo "medium: $MED"
[ -d "$BD" ] || { bad "no bounded models at $BD"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. premise
say "A. the premise this test rests on -- checked, not assumed"
note "MTP4 must be UNABLE to grow alone on ETO_SED, or a rescue proves nothing."
"$PY" - "$BD" "$MDB" "$MTP4" "$HYD" <<'PYEOF'
import sys, os, csv, glob, warnings
from reframed import load_cbmodel, Environment, FBA
bd, mdb, mtp4, hyd = sys.argv[1:5]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
print("   ETO_SED: %d compounds" % len(med))
for k in ('cpd11640','cpd00047','cpd00029','cpd00116','cpd00441'):
    print("     %s in medium: %s" % (k, k in med))
fmt=lambda x: "R_EX_{}_e".format(x)
alive=0
for p in sorted(glob.glob(os.path.join(bd,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    t=m.copy(); Environment.from_compounds(med, fmt_func=fmt).apply(t,inplace=True,warning=False)
    s=FBA(t); g=s.fobj if s and s.fobj else 0.0
    if g>1e-6: alive+=1
    tag = '  <-- grows alone' if g>1e-6 else ''
    print("   %-24s alone on ETO_SED = %.4f%s" % (n, g, tag))
print()
print("   %d of 11 can grow alone. The rest can only grow if something is transferred." % alive)
PYEOF

# ---------------------------------------------------------------- B. viability
say "B. which PAIRS can grow at all on ETO_SED?"
note "a pair that cannot grow reports no interaction for a reason that has nothing to"
note "do with interaction. Separating those two is the whole point of this section."
"$PY" - "$BD" "$MDB" "$MTP4" "$OUT" <<'PYEOF'
import sys, os, csv, glob, itertools, json, warnings
from reframed import load_cbmodel, Environment, FBA
from smetana.interface import Community
bd, mdb, mtp4, out = sys.argv[1:5]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
models={}
for p in sorted(glob.glob(os.path.join(bd,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        models[n]=load_cbmodel(p, flavor='fbc2')
env=lambda: Environment.from_compounds(med, fmt_func=lambda x:"R_EX_M_{}_e_pool".format(x))
viable=[]; mtp4_viable=[]
for a,b in itertools.combinations(sorted(models),2):
    try:
        c=Community('p',[models[a].copy(), models[b].copy()])
        mm=c.merged.copy(); env().apply(mm, inplace=True, warning=False)
        s=FBA(mm); g=s.fobj if s and s.fobj else 0.0
    except Exception as e:
        g=0.0
    if g>1e-6:
        viable.append((a,b,g))
        if mtp4 in (a,b): mtp4_viable.append((a,b,g))
print("   viable pairs (community growth > 0 on ETO_SED): %d of 55" % len(viable))
for a,b,g in sorted(viable, key=lambda x:-x[2])[:15]:
    print("     %-24s + %-24s  %.4f" % (a,b,g))
print()
print("   pairs containing MTP4 that grow: %d" % len(mtp4_viable))
for a,b,g in mtp4_viable:
    print("     %-24s + %-24s  %.4f" % (a,b,g))
json.dump({'viable':[[a,b,g] for a,b,g in viable]}, open(os.path.join(out,'viable_pairs.json'),'w'), indent=2)
if not mtp4_viable:
    print()
    print("   >>> NO pair containing MTP4 grows on ETO_SED. A pair is not enough; the")
    print("       fermenter cannot both feed itself and support the methanogen. Section E")
    print("       tries larger communities, which is where syntrophy would appear.")
PYEOF

# ---------------------------------------------------------------- C. pairwise
say "C. all 55 pairs on ETO_SED"
PD="$OUT/pairs"; mkdir -p "$PD"; cp "$BD"/*.xml "$PD/"
"$PY" - "$PD" > "$PD/communities.tsv" <<'PYEOF'
import sys, os, glob, itertools
d=sys.argv[1]
names=sorted(os.path.basename(p).replace('.xml','') for p in glob.glob(os.path.join(d,'*.xml')))
for i,(a,b) in enumerate(itertools.combinations(names,2),1):
    cid="P%02d_%s__%s" % (i,a,b)
    print("%s\t%s" % (cid,a)); print("%s\t%s" % (cid,b))
PYEOF
echo "   $(cut -f1 "$PD/communities.tsv" | sort -u | wc -l) pairs"
note "started $(date -Is)"
if [ -s "$PD/pairs_detailed.tsv" ]; then ok "already present"; else
  ( cd "$PD" && "$CLI" -o pairs -c communities.tsv -m "$MED" --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | tail -6 | sed 's/^/   /'
fi
DET="$PD/pairs_detailed.tsv"
[ -s "$DET" ] && ok "$(($(wc -l < "$DET")-1)) interaction row(s)" || no "no output"

# ---------------------------------------------------------------- D. answer
say "D. THE QUESTION: is anything energetic transferred on ETO_SED?"
if [ -s "$DET" ]; then
  echo "   --- H2 (cpd11640), formate (cpd00047), acetate (cpd00029) ---"
  tail -n +2 "$DET" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/{printf "     %-16s %-24s -> %-24s smetana=%s\n",$5,$4,$3,$9}' | head -30
  H=$(tail -n +2 "$DET" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/' | wc -l)
  echo "     total: $H row(s)"
  echo
  echo "   --- MTP4 as RECEIVER ---"
  tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%-8s (%s)\n",$5,$4,$9,$1}' | head -25
  echo "     total: $(tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0' | wc -l) row(s)"
  echo
  echo "   --- $HYD as RECEIVER ---"
  tail -n +2 "$DET" | awk -F'\t' -v a="$HYD" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%s\n",$5,$4,$9}' | head -15
  echo "     total: $(tail -n +2 "$DET" | awk -F'\t' -v a="$HYD" '$3==a && $9>0' | wc -l) row(s)"
  echo
  echo "   --- most-transferred overall ---"
  tail -n +2 "$DET" | awk -F'\t' '$9>0{s[$5]+=$9;n[$5]++} END{for(x in s) printf "     %-18s %7.2f  (%d)\n",x,s[x],n[x]}' | sort -k2,2gr | head -15
fi

# ---------------------------------------------------------------- E. quads
say "E. four-member communities built around MTP4 and the best H2 exporters"
note "43i measured maximum secretion on ETO_SED: LIQ_das_semibin_003 H2 169, formate 111;"
note "LIQ_das_semibin_016 H2 157; LIQ_das_semibin_006 H2 145; SED_das_semibin_002 H2 87."
note "If H2 syntrophy exists anywhere in this reconstruction, it is in these groups."
QD="$OUT/quads"; mkdir -p "$QD"; cp "$BD"/*.xml "$QD/"
cat > "$QD/communities.tsv" <<EOF
Q1_H2donors	$MTP4
Q1_H2donors	LIQ_das_semibin_003
Q1_H2donors	LIQ_das_semibin_016
Q1_H2donors	LIQ_das_semibin_006
Q2_mixed	$MTP4
Q2_mixed	LIQ_das_semibin_003
Q2_mixed	SED_das_semibin_002
Q2_mixed	SED_das_vamb_004
Q3_botharchaea	$MTP4
Q3_botharchaea	$HYD
Q3_botharchaea	LIQ_das_semibin_003
Q3_botharchaea	LIQ_das_semibin_016
EOF
sed 's/^/     /' "$QD/communities.tsv" | head -4
note "started $(date -Is) -- roughly 12 min per quad"
if [ -s "$QD/quads_detailed.tsv" ]; then ok "already present"; else
  ( cd "$QD" && "$CLI" -o quads -c communities.tsv -m "$MED" --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | tail -6 | sed 's/^/   /'
fi
QD_T="$QD/quads_detailed.tsv"
if [ -s "$QD_T" ]; then
  echo "   rows: $(($(wc -l < "$QD_T")-1))"
  echo "   --- energetic transfers in the quads ---"
  tail -n +2 "$QD_T" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/{printf "     %-12s %-16s %-24s -> %-24s %s\n",$1,$5,$4,$3,$9}' | head -20
  echo "     total: $(tail -n +2 "$QD_T" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/' | wc -l) row(s)"
  echo "   --- MTP4 as receiver in the quads ---"
  tail -n +2 "$QD_T" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     %-12s %-16s from %-24s %s\n",$1,$5,$4,$9}' | head -20
  echo "     total: $(tail -n +2 "$QD_T" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0' | wc -l) row(s)"
fi

# ---------------------------------------------------------------- F. ceiling
if [ "$PROBE" = "1" ]; then
  say "F. re-probe the size ceiling with a real time budget (PROBE=1)"
  note "43p's n=5 hit MY 15-minute cap, it did not fail. 45 min each here."
  SZ="$OUT/sizes"; mkdir -p "$SZ"
  ORDER=("$MTP4" LIQ_das_semibin_003 LIQ_das_semibin_006 SED_das_semibin_002 LIQ_das_semibin_016 SED_das_vamb_004 "$HYD")
  for N in 5 6 7; do
    D="$SZ/n$N"; rm -rf "$D"; mkdir -p "$D"
    for i in $(seq 0 $((N-1))); do cp "$BD/${ORDER[$i]}.xml" "$D/" 2>/dev/null; done
    START=$(date +%s)
    ( cd "$D" && timeout 2700 "$CLI" -o s -m "$MED" --mediadb "$MDB" --flavor fbc2 -d "*.xml" >/dev/null 2>"$D/err.txt" )
    RC=$?; EL=$(( $(date +%s) - START ))
    R=0; [ -s "$D/s_detailed.tsv" ] && R=$(($(wc -l < "$D/s_detailed.tsv")-1))
    if [ "$R" -gt 0 ]; then printf '   n=%-3d %5ds  %5d rows  OK\n' "$N" "$EL" "$R"
    elif [ "$RC" -eq 124 ]; then printf '   n=%-3d %5ds  timed out at 45 min\n' "$N" "$EL"; break
    else printf '   n=%-3d %5ds  FAILED: %s\n' "$N" "$EL" "$(grep -oE 'Exception:.*' "$D/err.txt" | tail -1)"; break; fi
  done
else
  note "PROBE=1 to re-probe the community-size ceiling with 45-minute budgets"
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
