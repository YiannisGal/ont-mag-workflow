#!/usr/bin/env bash
# 43n_fix_constraint_and_trace.sh -- verify the constraint properly, find the real gap.
#
# RUN IN tmux. Section E is the ~3 h detailed run and is resumable.
#
# TWO CORRECTIONS TO 43m, both mine.
#
# 1. THE CONSTRAINT TEST WAS SELF-DEFEATING.
#    Environment.apply() rewrites the bounds of every exchange in the medium. I set
#    lb=0 on the sugar exchanges and then applied ETO_SED, which overwrote it. Growth
#    was identical before and after because nothing was actually constrained during
#    the test.
#    In the COMMUNITY this matters less than it appears: SMETANA applies the medium to
#    the POOL reactions (R_EX_M_<cpd>_e_pool), not to each member's own exchange, so a
#    member-level lb=0 survives there. Section B verifies it the way the community will
#    see it -- apply the environment FIRST, then close the member's uptake -- instead of
#    in a way that erases it.
#
# 2. I MISREAD SECTION A AND CALLED A GAP THAT IS NOT THERE.
#    R_rxn03127_c, "2-(Methylthio)ethanesulfonate:N-(7-thioheptanoyl)-3-O-phospho-
#    threonine S-(2-sulfoethyl)thiotransferase", IS methyl-CoM reductase under its
#    systematic name. BOTH models carry it, both carry CoM, Methyl-CoM and CoM-S-S-CoB,
#    and both genomes carry mcrA/B/G:
#        LIQ_das_vamb_003  K00399 1e-299   K00401 1.6e-226  K00402 7.9e-124
#        LIQ_das_vamb_017  K00399 0        K00401 3.1e-227  K00402 4.8e-135
#    So "LIQ_das_vamb_017 has no mcr" was wrong. Withdrawn.
#
#    The real question is narrower: vamb_017 HAS mcr and still makes no methane on
#    H2/CO2 while growing at 0.0845. Something upstream in the CO2 -> methyl-CoM route
#    cannot carry flux, and biomass is coming from somewhere else. Section C locates it.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43n"; mkdir -p "$OUT/models"
export TMPDIR="${SCRATCH_ROOT}/tmp43n"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
SRC="$WS2/43_smetana/run_43l/models"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017; MCR=R_rxn03127_c
SKIP_DETAILED="${SKIP_DETAILED:-0}"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43n  $(date -Is)"; echo "out: $OUT"
[ -d "$SRC" ] || { bad "no models at $SRC"; echo "STEP_STATUS: FAILED"; exit 1; }
cp "$SRC"/*.xml "$OUT/models/"

# ---------------------------------------------------------------- A. mcr flux
say "A. can mcr carry flux at all, in each archaeon?"
note "objective set to $MCR directly. If it maximises to zero, methanogenesis is"
note "blocked upstream and the growth we see is coming from another route entirely."
"$PY" - "$OUT/models" "$MDB" "$MTP4" "$HYD" "$MCR" <<'PYEOF'
import sys, os, csv, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, mdb, mtp4, hyd, mcr = sys.argv[1:6]
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
fmt=lambda x: "R_EX_{}_e".format(x)
for g in (mtp4, hyd):
    p=os.path.join(mdir,g+'.xml')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    print("\n   [%s]" % g)
    if mcr not in m.reactions:
        print("     %s not in this model" % mcr); continue
    rx=m.reactions[mcr]
    print("     %s  lb=%s ub=%s" % (mcr, rx.lb, rx.ub))
    print("     stoichiometry:")
    for k,v in rx.stoichiometry.items():
        print("       %+6.1f  %-18s %s" % (v, k, m.metabolites[k].name))
    for med in ('ETO_METH','ETO_SED_FERM'):
        if med not in media: continue
        t=m.copy()
        Environment.from_compounds(media[med], fmt_func=fmt).apply(t,inplace=True,warning=False)
        t.set_objective({mcr: 1})
        s=FBA(t); v=s.fobj if s and s.fobj else 0.0
        print("     max flux through mcr on %-14s = %.6f" % (med, v))
    # unlimited uptake: is it the medium or the network?
    t=m.copy(); Environment.complete(t, max_uptake=1000.0).apply(t,inplace=True,warning=False)
    t.set_objective({mcr: 1}); s=FBA(t); v=s.fobj if s and s.fobj else 0.0
    print("     max flux with UNLIMITED uptake            = %.6f" % v)
    if v <= 1e-9:
        print("     >>> blocked by the NETWORK, not the medium. Tracing substrates:")
        for k,c in rx.stoichiometry.items():
            if c >= 0: continue
            t2=m.copy(); Environment.complete(t2, max_uptake=1000.0).apply(t2,inplace=True,warning=False)
            try: t2.add_reaction_from_str('SINK_X: %s --> ' % k)
            except Exception: continue
            t2.set_objective({'SINK_X': 1}); s2=FBA(t2)
            vv=s2.fobj if s2 and s2.fobj else 0.0
            print("        %-18s %-40s producible=%.4f %s"
                  % (k, m.metabolites[k].name[:40], vv, '<-- BLOCKED' if vv<=1e-9 else ''))
PYEOF

# ---------------------------------------------------------------- B. constraint
say "B. the sugar constraint, verified the way the COMMUNITY will see it"
note "environment applied FIRST, then the member's own uptake closed. That is the"
note "order SMETANA produces: medium on the pool, member bounds from the model file."
"$PY" - "$OUT/models" "$MDB" "$MTP4" "$HYD" <<'PYEOF'
import sys, os, csv, re, warnings
from reframed import load_cbmodel, save_cbmodel, Environment, FBA
mdir, mdb, mtp4, hyd = sys.argv[1:5]
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
SUGAR={'cpd00027':'D-Glucose','cpd00154':'Xylose','cpd00158':'Cellobiose','cpd00280':'D-Galacturonate'}
fmt=lambda x: "R_EX_{}_e".format(x)
for g in (mtp4, hyd):
    p=os.path.join(mdir,g+'.xml')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    present=[c for c in SUGAR if fmt(c) in m.reactions]
    print("\n   [%s] sugar exchanges present: %s" % (g, [SUGAR[c] for c in present]))
    for med in ('ETO_SED','ETO_SED_FERM','ETO_METH'):
        if med not in media: continue
        # unconstrained
        a=m.copy(); Environment.from_compounds(media[med], fmt_func=fmt).apply(a,inplace=True,warning=False)
        sa=FBA(a); ga=sa.fobj if sa and sa.fobj else 0.0
        # constrained: environment FIRST, then close uptake -- the community's order
        b=m.copy(); Environment.from_compounds(media[med], fmt_func=fmt).apply(b,inplace=True,warning=False)
        for c in present: b.reactions[fmt(c)].lb = 0.0
        sb=FBA(b); gb=sb.fobj if sb and sb.fobj else 0.0
        exch=set(b.get_exchange_reactions())
        sug=[round(v,2) for r,v in (sb.values.items() if sb and sb.values else [])
             if v<-1e-6 and r in exch and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in b.reactions[r].stoichiometry)]
        ch4=[v for r,v in (sb.values.items() if sb and sb.values else [])
             if r in exch and any(re.match(r'^methane',(b.metabolites[k].name or ''),re.I) for k in b.reactions[r].stoichiometry)]
        print("     %-14s growth %.4f -> %.4f   sugar_uptake=%s  methane=%.4f"
              % (med, ga, gb, sug or 'none', max(ch4) if ch4 else 0.0))
    # persist the constraint in the file for the community run
    for c in present: m.reactions[fmt(c)].lb = 0.0
    save_cbmodel(m, p, flavor='fbc2')
    print("     constraint written to %s" % os.path.basename(p))
PYEOF

# ---------------------------------------------------------------- C. gate
say "C. GATE -- what may be claimed, per organism"
"$PY" - "$OUT/models" "$MDB" "$MTP4" "$HYD" "$MCR" <<'PYEOF'
import sys, os, csv, re, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, mdb, mtp4, hyd, mcr = sys.argv[1:6]
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
fmt=lambda x: "R_EX_{}_e".format(x)
ok_any=False
for g in (mtp4, hyd):
    p=os.path.join(mdir,g+'.xml')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    t=m.copy()
    Environment.from_compounds(media['ETO_METH'], fmt_func=fmt).apply(t,inplace=True,warning=False)
    s=FBA(t); gr=s.fobj if s and s.fobj else 0.0
    exch=set(t.get_exchange_reactions())
    ch4=[v for r,v in (s.values.items() if s and s.values else [])
         if r in exch and any(re.match(r'^methane',(t.metabolites[k].name or ''),re.I) for k in t.reactions[r].stoichiometry)]
    f=max(ch4) if ch4 else 0.0
    verdict = "METHANOGENIC" if f>1e-6 else "mcr present but pathway does not carry flux"
    print("   %-24s ETO_METH growth=%-9.4f methane=%-9.4f  %s" % (g, gr, f, verdict))
    if f>1e-6: ok_any=True
print()
if ok_any:
    print("GATE_STATUS: PASSED  at least one archaeon is a working methanogen")
    print("  Any organism marked 'pathway does not carry flux' stays in the community as")
    print("  an archaeon, but NO methanogenesis claim may rest on it. Section A says why.")
else:
    print("GATE_STATUS: FAILED  no working methanogen in the community")
import sys as _s; _s.exit(0 if ok_any else 3)
PYEOF
[ $? -eq 0 ] || { bad "no working methanogen -- stopping"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- D. control + global
say "D. control, then global scores"
CD="$OUT/control"; mkdir -p "$CD"
if ls "$WS2/43_smetana/run_43g"/AUXO_*.xml >/dev/null 2>&1; then
  cp "$WS2/43_smetana/run_43g"/AUXO_*.xml "$CD/"
  ( cd "$CD" && "$CLI" -o ctrl -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d AUXO_*.xml 2>&1 ) | tail -4 | sed 's/^/   /'
  R=0; [ -s "$CD/ctrl_detailed.tsv" ] && R=$(($(wc -l < "$CD/ctrl_detailed.tsv")-1))
  echo "   control rows: $R"
  [ "$R" -gt 0 ] || { bad "control zero -- nothing below is interpretable"; echo "STEP_STATUS: FAILED"; exit 1; }
  ok "detector confirmed on this run"
else no "no control pair -- zeros below are NOT interpretable"; fi
TD="$OUT/community"; mkdir -p "$TD"; cp "$OUT/models"/*.xml "$TD/"
echo "   $(ls "$TD"/*.xml | wc -l) models"
for MED in ETO_SED_FERM ETO_SED ETO_METH; do
  awk -F'\t' -v m="$MED" 'NR==1{next} $1==m{n++} END{exit (n+0==0)}' "$MDB" || { no "$MED absent"; continue; }
  note "$MED"
  ( cd "$TD" && "$CLI" -o "g_$MED" -m "$MED" --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | tail -4 | sed 's/^/     /'
  [ -s "$TD/g_${MED}_global.tsv" ] && column -t -s$'\t' "$TD/g_${MED}_global.tsv" | sed 's/^/     /'
done
[ "$SKIP_DETAILED" = "1" ] && { note "SKIP_DETAILED=1"; say DONE; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0; }

# ---------------------------------------------------------------- E. detailed
say "E. detailed on ETO_SED_FERM -- ~3 h, resumable"
note "started $(date -Is)"
if [ -s "$TD/d_FERM_detailed.tsv" ]; then ok "already present"; else
  ( cd "$TD" && "$CLI" -o d_FERM -m ETO_SED_FERM --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | tail -6 | sed 's/^/   /'
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
      C=$(tail -n +2 "$DET" | awk -F'\t' -v a="$a" '$3==a && $9>0' | wc -l)
      echo "     total: $C row(s)"
    done
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
