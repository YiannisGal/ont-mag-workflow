#!/usr/bin/env bash
# 43p_pairwise.sh -- SCIP cannot take the eleven. Ask the question pairwise instead.
#
# RUN IN tmux.
#
# WHAT 43o ESTABLISHED
#   - reframed 1.6.0's solver registry is ['scip']. highs and glpk import, but reframed
#     does not register them, so --solver has nothing to switch to. One solver, take it
#     or leave it.
#   - 21,959 infinite bounds replaced by +/-1000 across the eleven models, and growth on
#     ETO_SED_FERM was UNCHANGED for every one. The transformation is safe.
#   - and the eleven STILL crash. So the infinities that break SCIP are introduced
#     downstream by SMETANA's own merge, not by the models. Patching a library's solver
#     backend to get a thesis result is not something I am willing to do quietly.
#   - 2 models: 59 rows. 3 real models: 36 rows, including MTP4 donating CoA (0.28),
#     formaldehyde (0.26) and serine (0.06) to LIQ_das_semibin_006.
#
# THE FAILURE SCALES WITH COMMUNITY SIZE, NOT WITH THE MODELS. So this script finds the
# largest size that works and asks the question at that size. Pairwise SMETANA is
# standard practice, gives a full interaction matrix rather than one aggregate number,
# and is provably tractable here.
#
# Section A first checks something 43n may only have appeared to do.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43p"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp43p"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
BD="$WS2/43_smetana/run_43o/models_bounded"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017; MED=ETO_SED_FERM
export PYTHONWARNINGS="ignore::UserWarning"
echo "43p  $(date -Is)"; echo "out: $OUT"
[ -d "$BD" ] || { bad "no bounded models at $BD -- run 43o first"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. constraint
say "A. did the sugar constraint actually PERSIST in the files?"
note "43o section C showed MTP4 at 0.5421 on ETO_SED_FERM -- its PRE-constraint value."
note "That test applies the environment, which resets exchange bounds, so it could not"
note "have shown the stored bound either way. This reads the bound off the file."
"$PY" - "$BD" "$MTP4" "$HYD" <<'PYEOF'
import sys, os, warnings
from reframed import load_cbmodel
bd, mtp4, hyd = sys.argv[1:4]
SUGAR={'cpd00027':'D-Glucose','cpd00154':'Xylose','cpd00158':'Cellobiose','cpd00280':'D-Galacturonate'}
bad=0
for g in (mtp4, hyd):
    p=os.path.join(bd, g+'.xml')
    if not os.path.exists(p): print("   %s absent" % g); continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    print("\n   [%s]" % g)
    for c,n in SUGAR.items():
        r="R_EX_%s_e" % c
        if r not in m.reactions: continue
        lb=m.reactions[r].lb
        state = "uptake BLOCKED" if (lb is not None and lb >= 0) else "uptake OPEN"
        print("     %-14s %-22s lb=%-9s %s" % (n, r, lb, state))
        if lb is None or lb < 0: bad+=1
print()
if bad:
    print("CONSTRAINT_STATUS: NOT PERSISTED (%d open) -- re-applying to the bounded copies" % bad)
else:
    print("CONSTRAINT_STATUS: PERSISTED")
sys.exit(1 if bad else 0)
PYEOF
if [ $? -ne 0 ]; then
  note "re-applying the declared constraint to the bounded models"
  "$PY" - "$BD" "$MTP4" "$HYD" <<'PYEOF'
import sys, os, warnings
from reframed import load_cbmodel, save_cbmodel
bd, mtp4, hyd = sys.argv[1:4]
SUGAR=['cpd00027','cpd00154','cpd00158','cpd00280']
for g in (mtp4, hyd):
    p=os.path.join(bd,g+'.xml')
    if not os.path.exists(p): continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    n=0
    for c in SUGAR:
        r="R_EX_%s_e" % c
        if r in m.reactions: m.reactions[r].lb=0.0; n+=1
    save_cbmodel(m, p, flavor='fbc2')
    print("   %s: %d sugar exchange(s) closed and written" % (g, n))
PYEOF
fi

# ---------------------------------------------------------------- B. size limit
say "B. how large a community can SCIP actually take?"
note "detailed mode on subsets of increasing size, 15 min each. This tells us the"
note "ceiling instead of guessing it, and every subset includes MTP4."
SZ="$OUT/sizes"; mkdir -p "$SZ"
ORDER=("$MTP4" LIQ_das_semibin_003 LIQ_das_semibin_006 SED_das_semibin_002 \
       LIQ_das_semibin_016 SED_das_vamb_004 "$HYD" LIQ_das_semibin_004 \
       SED_das_metabat_018 SED_das_semibin_003 LIQ_dassub_semibin_017)
MAXOK=0
for N in 3 4 5 6 7 8 9 10 11; do
  D="$SZ/n$N"; rm -rf "$D"; mkdir -p "$D"
  for i in $(seq 0 $((N-1))); do cp "$BD/${ORDER[$i]}.xml" "$D/" 2>/dev/null; done
  C=$(ls "$D"/*.xml 2>/dev/null | wc -l)
  [ "$C" -eq "$N" ] || { no "n=$N: only $C models available"; continue; }
  START=$(date +%s)
  ( cd "$D" && timeout 900 "$CLI" -o s -m "$MED" --mediadb "$MDB" --flavor fbc2 -d "*.xml" >/dev/null 2>"$D/err.txt" )
  RC=$?; EL=$(( $(date +%s) - START ))
  R=0; [ -s "$D/s_detailed.tsv" ] && R=$(($(wc -l < "$D/s_detailed.tsv")-1))
  if [ "$R" -gt 0 ]; then
    printf '   n=%-3d %4ds  %5d rows  OK\n' "$N" "$EL" "$R"; MAXOK=$N
  elif [ "$RC" -eq 124 ]; then
    printf '   n=%-3d %4ds  TIMED OUT at 15 min\n' "$N" "$EL"; break
  else
    printf '   n=%-3d %4ds  FAILED: %s\n' "$N" "$EL" "$(grep -oE 'Exception:.*' "$D/err.txt" | tail -1)"; break
  fi
done
echo
if [ "$MAXOK" -ge 3 ]; then ok "largest working community: $MAXOK members"
else bad "nothing above 2 members works"; fi

# ---------------------------------------------------------------- C. pairwise
say "C. all 55 pairs, via SMETANA's own -c subcommunity mechanism"
PD="$OUT/pairs"; mkdir -p "$PD"; cp "$BD"/*.xml "$PD/"
"$PY" - "$PD" > "$PD/communities.tsv" <<'PYEOF'
import sys, os, glob, itertools
d=sys.argv[1]
names=sorted(os.path.basename(p).replace('.xml','') for p in glob.glob(os.path.join(d,'*.xml')))
for i,(a,b) in enumerate(itertools.combinations(names,2), 1):
    cid = "P%02d_%s__%s" % (i, a, b)
    print("%s\t%s" % (cid, a)); print("%s\t%s" % (cid, b))
PYEOF
NP=$(cut -f1 "$PD/communities.tsv" | sort -u | wc -l)
echo "   $NP pairs from $(ls "$PD"/*.xml | wc -l) models"
head -4 "$PD/communities.tsv" | sed 's/^/     /'
note "started $(date -Is) -- roughly a minute a pair"
if [ -s "$PD/pairs_detailed.tsv" ]; then ok "already present"; else
  ( cd "$PD" && "$CLI" -o pairs -c communities.tsv -m "$MED" --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | tail -8 | sed 's/^/   /'
fi

DET="$PD/pairs_detailed.tsv"
if [ ! -s "$DET" ]; then
  bad "no pairwise output"; echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
N=$(($(wc -l < "$DET")-1)); ok "$N interaction row(s) across $NP pairs"

# ---------------------------------------------------------------- D. read it
say "D. what the pairwise matrix says"
echo "   --- pairs with any interaction ---"
tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$1]++} END{for(x in c) printf "     %-46s %d\n",x,c[x]}' | sort -k2,2nr | head -20
echo
echo "   --- MTP4 ($MTP4) as RECEIVER ---"
tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%-8s (%s)\n",$5,$4,$9,$1}' | head -30
echo "     total: $(tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0' | wc -l) row(s)"
echo
echo "   --- MTP4 as DONOR ---"
tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$4==a && $9>0{printf "     %-16s to %-26s smetana=%-8s\n",$5,$3,$9}' | sort -u | head -25
echo "     total: $(tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$4==a && $9>0' | wc -l) row(s)"
echo
echo "   --- H2 (cpd11640), formate (cpd00047), acetate (cpd00029) transfers, anywhere ---"
tail -n +2 "$DET" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/{printf "     %-16s %-24s -> %-24s smetana=%s\n",$5,$4,$3,$9}' | head -25
H=$(tail -n +2 "$DET" | awk -F'\t' '$9>0 && $5 ~ /cpd11640|cpd00047|cpd00029/' | wc -l)
echo "     total: $H row(s)"
[ "$H" -eq 0 ] && note "no H2/formate/acetate transfer anywhere, with a passing detector."
echo
echo "   --- most-transferred compounds overall ---"
tail -n +2 "$DET" | awk -F'\t' '$9>0{s[$5]+=$9; n[$5]++} END{for(x in s) printf "     %-18s %7.2f  (%d transfers)\n",x,s[x],n[x]}' | sort -k2,2gr | head -20

say "E. the largest workable community, if larger than a pair"
if [ "$MAXOK" -ge 4 ]; then
  BIG="$SZ/n$MAXOK"
  echo "   $MAXOK members:"; ls "$BIG"/*.xml | xargs -n1 basename | sed 's/^/     /'
  if [ -s "$BIG/s_detailed.tsv" ]; then
    echo "   --- MTP4 in the $MAXOK-member community ---"
    tail -n +2 "$BIG/s_detailed.tsv" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     receives %-16s from %-24s %s\n",$5,$4,$9}' | head -15
    tail -n +2 "$BIG/s_detailed.tsv" | awk -F'\t' -v a="$MTP4" '$4==a && $9>0{printf "     donates  %-16s to   %-24s %s\n",$5,$3,$9}' | head -15
  fi
else note "no community above 3 completed; the pairwise matrix is the result"; fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
