#!/usr/bin/env bash
# 43o_solver_fix.sh -- make the detailed run survive, and prove it on a small case first.
#
# RUN IN tmux. Section F is the long run and only starts if D and E pass.
#
# 43n's RESULTS, and two corrections to what I said about them.
#
#   mcr CARRIES FLUX IN BOTH ARCHAEA:
#       LIQ_das_vamb_003  2.500 on ETO_METH   51.250 on ETO_SED_FERM
#       LIQ_das_vamb_017  2.391 on ETO_METH    5.446 on ETO_SED_FERM
#   So "pathway does not carry flux" was WRONG. Withdrawn. The true statement is
#   narrower and unremarkable: methanogenesis is NOT COUPLED TO BIOMASS at the FBA
#   optimum in vamb_017 -- another route to growth exists, so maximising biomass leaves
#   mcr at zero. Alternate optima, not a reconstruction gap. (Third time in this stretch
#   I have called a gap that was my reading; see CORRECTIONS_INDEX.)
#
#   THE CONSTRAINT WORKS, and it changed the biology the right way:
#       LIQ_das_vamb_003 ETO_SED_FERM  0.5421 -> 0.0944 growth, methane 16.88, no sugar
#   MTP4 now grows by making methane from fermentation products instead of eating
#   glucose. That is the first time in this exercise the organism has behaved as a
#   methanogen inside the community's environment.
#
#   Control 59 rows. MRO on ETO_SED_FERM = 0.6034 for the eleven. MIP n/a everywhere,
#   which 43h already explained: it needs every member at min_growth simultaneously.
#
# THE FAILURE THIS SCRIPT FIXES
#     SCIP: error in input data!   -- one hour into the detailed run.
#   The GLOBAL runs succeeded on the same models, so the problem is the MILP that
#   detailed mode builds. Infinite bounds are the usual cause. But I am not fixing that
#   blind and spending three more hours finding out: sections D and E prove the fix on a
#   2-member and a 3-member case, which take minutes, BEFORE section F commits.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43o"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp43o"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
SRC="$WS2/43_smetana/run_43n/models"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
export PYTHONWARNINGS="ignore::UserWarning"
echo "43o  $(date -Is)"; echo "out: $OUT"
[ -d "$SRC" ] || { bad "no constrained models at $SRC -- run 43n first"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. solvers
say "A. which solvers does this environment actually have?"
"$PY" - <<'PYEOF'
import importlib
avail=[]
for name, mod in (('gurobi','gurobipy'), ('cplex','cplex'), ('scip','pyscipopt'),
                  ('highs','highspy'), ('glpk','swiglpk')):
    try:
        importlib.import_module(mod); avail.append(name)
        print("   [OK]   %-8s (%s)" % (name, mod))
    except Exception as e:
        print("   [--]   %-8s (%s) %s" % (name, mod, type(e).__name__))
try:
    import reframed.solvers as S
    print("   reframed knows: %s" % [k for k in dir(S) if 'solver' in k.lower()][:8])
    from reframed.solvers import solvers as _s
    print("   reframed solver registry: %s" % list(_s.keys()))
except Exception as e:
    print("   could not read reframed's registry: %s: %s" % (type(e).__name__, e))
print("\n   usable: %s" % (avail or 'NONE -- that is the real problem'))
PYEOF

# ---------------------------------------------------------------- B. bounds
say "B. audit the models for what upsets a MILP: inf, nan, lb>ub"
"$PY" - "$SRC" <<'PYEOF'
import sys, os, glob, math, warnings
from reframed import load_cbmodel
src=sys.argv[1]
tot_inf=tot_nan=tot_bad=0
for p in sorted(glob.glob(os.path.join(src,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    inf=nan=bad=0
    for r,rx in m.reactions.items():
        for b in (rx.lb, rx.ub):
            if b is None: continue
            if math.isinf(b): inf+=1
            elif math.isnan(b): nan+=1
        if rx.lb is not None and rx.ub is not None and rx.lb > rx.ub: bad+=1
    tot_inf+=inf; tot_nan+=nan; tot_bad+=bad
    print("   %-24s %5d rxn  inf=%-5d nan=%-3d lb>ub=%d" % (n, len(m.reactions), inf, nan, bad))
print("\n   totals: inf=%d  nan=%d  lb>ub=%d" % (tot_inf, tot_nan, tot_bad))
if tot_inf: print("   >>> infinite bounds are the prime suspect for 'error in input data'")
if tot_nan or tot_bad: print("   >>> nan or inverted bounds would break any solver")
PYEOF

# ---------------------------------------------------------------- C. bounded copies
say "C. bounded copies -- inf replaced by +/-1000, nothing else changed"
BD="$OUT/models_bounded"; mkdir -p "$BD"
"$PY" - "$SRC" "$BD" "$MDB" <<'PYEOF'
import sys, os, glob, math, csv, warnings
from reframed import load_cbmodel, save_cbmodel, Environment, FBA
src, dst, mdb = sys.argv[1:4]
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
fmt=lambda x: "R_EX_{}_e".format(x)
BIG=1000.0
print("   %-24s %10s %10s   growth before/after on ETO_SED_FERM" % ('model','inf_lb','inf_ub'))
for p in sorted(glob.glob(os.path.join(src,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    a=m.copy(); Environment.from_compounds(media['ETO_SED_FERM'], fmt_func=fmt).apply(a,inplace=True,warning=False)
    sa=FBA(a); ga=sa.fobj if sa and sa.fobj else 0.0
    il=iu=0
    for r,rx in m.reactions.items():
        if rx.lb is not None and math.isinf(rx.lb): rx.lb=-BIG; il+=1
        if rx.ub is not None and math.isinf(rx.ub): rx.ub= BIG; iu+=1
    out=os.path.join(dst, n+'.xml')
    save_cbmodel(m, out, flavor='fbc2')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m2=load_cbmodel(out, flavor='fbc2')
    b=m2.copy(); Environment.from_compounds(media['ETO_SED_FERM'], fmt_func=fmt).apply(b,inplace=True,warning=False)
    sb=FBA(b); gb=sb.fobj if sb and sb.fobj else 0.0
    flag = '' if abs(ga-gb) < 1e-6 else '   <-- CHANGED'
    print("   %-24s %10d %10d   %.4f -> %.4f%s" % (n, il, iu, ga, gb, flag))
print("\n   bounding is only safe if growth is unchanged. Any CHANGED line above means")
print("   1000 was a binding limit for that organism and the number must be raised.")
PYEOF

# ---------------------------------------------------------------- D. small test
say "D. does detailed mode survive on the 2-member control? (minutes)"
T2="$OUT/t2"; mkdir -p "$T2"
cp "$WS2/43_smetana/run_43g"/AUXO_*.xml "$T2/" 2>/dev/null || { bad "no control pair"; exit 1; }
RUN_OK=0
for SOLV in "" scip glpk gurobi cplex highs; do
  [ -n "$SOLV" ] && ARG="--solver $SOLV" || ARG=""
  rm -f "$T2"/t2_detailed.tsv
  note "solver: ${SOLV:-default}"
  ( cd "$T2" && "$CLI" -o t2 -m ETO_SED --mediadb "$MDB" --flavor fbc2 $ARG -d AUXO_*.xml 2>&1 ) | tail -4 | sed 's/^/     /'
  R=0; [ -s "$T2/t2_detailed.tsv" ] && R=$(($(wc -l < "$T2/t2_detailed.tsv")-1))
  echo "     rows: $R"
  if [ "$R" -gt 0 ]; then ok "works with ${SOLV:-default}"; GOODSOLV="$SOLV"; RUN_OK=1; break; fi
done
[ "$RUN_OK" -eq 1 ] || { bad "no solver completes even the 2-member control"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- E. three real models
say "E. and on three REAL models, bounded, including MTP4? (minutes)"
T3="$OUT/t3"; mkdir -p "$T3"; rm -f "$T3"/*.xml
cp "$BD/$MTP4.xml" "$BD/LIQ_das_semibin_003.xml" "$BD/LIQ_das_semibin_006.xml" "$T3/" 2>/dev/null
ls "$T3" | sed 's/^/     /'
ARG=""; [ -n "${GOODSOLV:-}" ] && ARG="--solver $GOODSOLV"
( cd "$T3" && "$CLI" -o t3 -m ETO_SED_FERM --mediadb "$MDB" --flavor fbc2 $ARG -d "*.xml" 2>&1 ) | tail -6 | sed 's/^/     /'
R3=0; [ -s "$T3/t3_detailed.tsv" ] && R3=$(($(wc -l < "$T3/t3_detailed.tsv")-1))
echo "   rows: $R3"
if [ "$R3" -gt 0 ]; then
  head -6 "$T3/t3_detailed.tsv" | column -t -s$'\t' | sed 's/^/     /'
  ok "three real models complete with ${GOODSOLV:-default} on bounded inputs"
else
  bad "three real models still fail -- the fix is not the bounds or the solver."
  note "section F is NOT started. Three hours would be wasted."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi

# ---------------------------------------------------------------- F. full run
say "F. the eleven, detailed, on ETO_SED_FERM -- ~3 h, resumable"
TD="$OUT/community"; mkdir -p "$TD"; cp "$BD"/*.xml "$TD/"
echo "   $(ls "$TD"/*.xml | wc -l) models, solver=${GOODSOLV:-default}"
note "started $(date -Is)"
if [ -s "$TD/d_FERM_detailed.tsv" ]; then ok "already present"; else
  ( cd "$TD" && "$CLI" -o d_FERM -m ETO_SED_FERM --mediadb "$MDB" --flavor fbc2 $ARG -d "*.xml" 2>&1 ) | tail -8 | sed 's/^/   /'
fi
DET="$TD/d_FERM_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET")-1)); echo "   rows: $N"
  if [ "$N" -gt 0 ]; then
    echo "   --- strongest 20 ---"
    { head -1 "$DET"; tail -n +2 "$DET" | sort -t$'\t' -k9,9gr | head -20; } | column -t -s$'\t' | sed 's/^/     /'
    echo "   --- donors ---";    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    echo "   --- receivers ---"; tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    for a in "$MTP4" "$HYD"; do
      echo "   --- what $a receives ---"
      tail -n +2 "$DET" | awk -F'\t' -v a="$a" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%s\n",$5,$4,$9}' | head -20
      echo "     total: $(tail -n +2 "$DET" | awk -F'\t' -v a="$a" '$3==a && $9>0' | wc -l) row(s)"
    done
    echo "   --- H2 / formate / acetate transfers specifically ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0 && ($5 ~ /cpd11640|cpd00047|cpd00029/){printf "     %-16s %-24s -> %-24s smetana=%s\n",$5,$4,$3,$9}' | head -20
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
