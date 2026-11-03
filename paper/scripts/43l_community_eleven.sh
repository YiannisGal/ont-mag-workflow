#!/usr/bin/env bash
# 43l_community_eleven.sh -- the eleven-member community, at last.
#
# RUN IN tmux. Section G is the ~3 h detailed run and is resumable.
#     tmux new -s 43l
#     bash scripts/43l_community_eleven.sh 2>&1 | tee logs/43l.log
#
# WHERE WE ARE
#   MTP4 (LIQ_das_vamb_003) built and gated:
#       ETO_METH  growth 0.0110, H2 -10.00 and CO2 -2.72 in, Methane 2.27 out, no sugar.
#   That is hydrogenotrophic methanogenesis. dat/seq/Archaea is installed and the
#   pathway table mentions archaea, so the annotation was not bounded by a bacterial
#   database -- that worry is closed.
#
#   CAVEAT TO CARRY, not to fix: on ETO_SED MTP4 takes up glucose at -10.00 and emits
#   methane at 21.3. Gapfilling on ETO_METH added only 19 reactions, so that glycolysis
#   came from the GENOME ANNOTATION, not from my medium. It is still not real biology,
#   and it makes ETO_SED the wrong environment for this organism. ETO_SED_FERM is the
#   meaningful one. Both are reported.
#
# THE TRAP THIS SCRIPT EXISTS TO AVOID
#   The new archaeal models use `_e0`. The nine bacterial models in 43_smetana/models
#   were relabelled to `_e`. Merge them as they are and cpd00009_e0 and cpd00009_e are
#   TWO DIFFERENT pool metabolites: nothing is shared, every score is zero, and nothing
#   in the output says why. That is the 43c failure with a new cause. Section C
#   relabels and then ASSERTS the identifier spaces actually overlap.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43l"; mkdir -p "$OUT/models"
export TMPDIR="${SCRATCH_ROOT}/tmp43l"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
GAPSEQ="$WS2/envs/gapseq_local/gapseq"
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONWARNINGS="ignore::UserWarning"
ARCH="$WS2/43k_archaea"; BACT="$WS2/43_smetana/models"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
METHCSV="$WS2/41_env_checks/ETO_METH-medium.csv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
SKIP_DETAILED="${SKIP_DETAILED:-0}"
echo "43l  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- A. archaea
say "A. the two archaeal models"
for g in "$MTP4" "$HYD"; do
  m=$(ls "$ARCH/$g"/*.xml 2>/dev/null | grep -v draft | head -1)
  if [ -n "$m" ]; then ok "$g: $(basename "$m")"
  else no "$g: no filled model yet"; fi
done
HM=$(ls "$ARCH/$HYD"/*.xml 2>/dev/null | grep -v draft | head -1)
if [ -z "$HM" ]; then
  note "refilling $HYD on ETO_METH (43k2 section D may not have run)"
  SRC="$WS2/42_gems/$HYD"; DRAFT="$SRC/$HYD-draft.RDS"
  if [ -s "$DRAFT" ]; then
    mkdir -p "$ARCH/$HYD"
    WT="$SRC/$HYD-rxnWeights.RDS"; XG="$SRC/$HYD-rxnXgenes.RDS"
    ( cd "$ARCH/$HYD" && "$GAPSEQ" fill -m "$DRAFT" -n "$METHCSV" \
        ${WT:+-c "$WT"} ${XG:+-g "$XG"} ) 2>&1 | tail -12 | sed 's/^/     /'
    RC=${PIPESTATUS[0]:-1}; [ "$RC" -eq 0 ] && ok "fill exited 0" || no "fill exited $RC"
    HM=$(ls "$ARCH/$HYD"/*.xml 2>/dev/null | grep -v draft | head -1)
  else no "no draft at $DRAFT"; fi
fi
[ -n "$HM" ] && ok "$HYD model present" || \
  note "$HYD unavailable -- the community will run with MTP4 as the sole archaeon"

# ---------------------------------------------------------------- B. relabel
say "B. relabel the archaeal models to the bacterial convention"
note "_e0 -> _e and _c0 -> _c, leaving _p0 alone, exactly as the bacteria were done."
for g in "$MTP4" "$HYD"; do
  src=$(ls "$ARCH/$g"/*.xml 2>/dev/null | grep -v draft | head -1)
  [ -n "$src" ] || continue
  dst="$OUT/models/$g.xml"
  sed -e 's/_e0/_e/g' -e 's/_c0/_c/g' "$src" > "$dst"
  printf '     %-24s %s -> %s bytes\n' "$g" "$(wc -c < "$src")" "$(wc -c < "$dst")"
done
note "and the nine bacteria, copied unchanged (excluding the old $HYD):"
n=0
for f in "$BACT"/*.xml; do
  b=$(basename "$f" .xml)
  [ "$b" = "$HYD" ] && continue
  cp "$f" "$OUT/models/$b.xml"; n=$((n+1))
done
ok "$n bacterial models copied"
ls "$OUT/models" | sed 's/^/     /'

# ---------------------------------------------------------------- C. assert
say "C. ASSERT the identifier spaces overlap -- the check 43c did not have"
"$PY" - "$OUT/models" "$MTP4" <<'PYEOF'
import sys, os, glob, re, warnings
from reframed import load_cbmodel
mdir, mtp4 = sys.argv[1], sys.argv[2]
sets={}
for p in sorted(glob.glob(os.path.join(mdir,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    ex=set(m.get_exchange_reactions())
    sfx=set()
    for r in ex:
        t=re.match(r'R_EX_(cpd\d+)_([A-Za-z]\d*)$', r)
        if t: sfx.add(t.group(2))
    sets[n]=(ex, sfx, len(m.reactions))
    print("   %-24s %5d rxn  %4d exch  suffixes=%s" % (n, len(m.reactions), len(ex), sorted(sfx)))
# EVERY gapseq model carries R_EX_cpd11416_c, the biomass DRAIN, which is a cytosolic
# exchange. So {'c','e'} is the normal state, not a mixture. v1 of this check demanded
# a single suffix and failed all eleven identically-named models. The question is not
# "how many suffixes" but "does every model use the SAME convention".
shapes = {}
for n,(_,s,_2) in sets.items(): shapes.setdefault(frozenset(s), []).append(n)
print("\n   naming conventions present: %d" % len(shapes))
for shape, members in shapes.items():
    print("     %-14s %d model(s): %s" % (sorted(shape), len(members),
          ', '.join(members[:3]) + (' ...' if len(members) > 3 else '')))
if len(shapes) > 1:
    print("   [FAIL] MIXED EXCHANGE NAMING. Merged like this, the same compound in two")
    print("          conventions becomes two pool metabolites and nothing is shared.")
    sys.exit(3)
# and the metabolite-side convention, which is what the pool is actually keyed on
extsfx = set()
for n,(ex,_,_2) in sets.items():
    for r in ex:
        t = re.match(r'R_EX_(cpd\d+)_(e\d*)$', r)
        if t: extsfx.add(t.group(2))
print("   external-compartment suffix in use: %s" % sorted(extsfx))
if len(extsfx) != 1:
    print("   [FAIL] more than one external naming -- the pool cannot be shared")
    sys.exit(3)
arc = sets.get(mtp4)
if arc:
    others=[k for k in sets if k!=mtp4]
    ov=[(k, len(arc[0] & sets[k][0])) for k in others]
    ov.sort(key=lambda x:-x[1])
    print("\n   exchange reactions %s shares with each bacterium:" % mtp4)
    for k,v in ov: print("     %-24s %d" % (k,v))
    worst=min(v for _,v in ov)
    print("\n   smallest overlap: %d" % worst)
    if worst < 20:
        print("   [FAIL] too little shared vocabulary for a community to mean anything")
        sys.exit(3)
print("\n   [OK] one naming convention, and the archaeon shares exchanges with every member")
PYEOF
[ $? -eq 0 ] || { bad "identifier check failed -- not running SMETANA on this"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- D. gate
say "D. the archaea still behave after relabelling"
"$PY" - "$OUT/models" "$MDB" "$MTP4" "$HYD" <<'PYEOF'
import sys, os, csv, re, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, mdb, mtp4, hyd = sys.argv[1:5]
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
SUGAR={'cpd00027','cpd00154','cpd00158','cpd00280'}
fails=0
for g in (mtp4, hyd):
    p=os.path.join(mdir,g+'.xml')
    if not os.path.exists(p): print("   %s: absent" % g); continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    fmt=lambda x: "R_EX_{}_e".format(x)
    for med in ('ETO_METH','ETO_SED','ETO_SED_FERM'):
        if med not in media: continue
        mm=m.copy()
        Environment.from_compounds(media[med], fmt_func=fmt).apply(mm,inplace=True,warning=False)
        s=FBA(mm); gr=s.fobj if s and s.fobj else 0.0
        exch=set(mm.get_exchange_reactions())
        rows=[(r,v) for r,v in (s.values.items() if s and s.values else []) if r in exch and abs(v)>1e-6]
        ch4=[v for r,v in rows if any(re.match(r'^methane',(mm.metabolites[k].name or ''),re.I) for k in mm.reactions[r].stoichiometry)]
        sug=[round(v,2) for r,v in rows if v<0 and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in mm.reactions[r].stoichiometry)]
        print("   %-22s %-14s growth=%-9.4f methane=%-9.4f sugar=%s"
              % (g, med, gr, max(ch4) if ch4 else 0.0, sug or 'none'))
        if med=='ETO_METH' and (gr<=1e-6 or not ch4):
            fails+=1
print()
print("ARCH_GATE: %s" % ("PASSED" if fails==0 else "FAILED"))
sys.exit(0 if fails==0 else 3)
PYEOF
[ $? -eq 0 ] || { bad "an archaeon lost its methanogenesis in relabelling"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "relabelling did not change the biology"

# ---------------------------------------------------------------- E. control
say "E. the positive control, re-run against THIS model set"
note "43g's control passed on a different set of files. A control does not transfer."
CD="$OUT/control"; mkdir -p "$CD"
if [ -s "$WS2/43_smetana/run_43g/AUXO_cpd00220.xml" ] && [ -s "$WS2/43_smetana/run_43g/AUXO_cpd00039.xml" ]; then
  cp "$WS2/43_smetana/run_43g"/AUXO_*.xml "$CD/"
  ( cd "$CD" && "$CLI" -o ctrl -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d AUXO_*.xml 2>&1 ) | tail -6 | sed 's/^/   /'
  R=0; [ -s "$CD/ctrl_detailed.tsv" ] && R=$(($(wc -l < "$CD/ctrl_detailed.tsv")-1))
  echo "   control rows: $R"
  [ "$R" -gt 0 ] && ok "detector still registers cross-feeding" || { bad "control now gives zero -- stop"; echo "STEP_STATUS: FAILED"; exit 1; }
else
  no "43g's auxotroph pair not found -- cannot re-verify the detector"
  note "proceeding, but a zero below is then NOT interpretable"
fi

# ---------------------------------------------------------------- F. global
say "F. global scores (MIP/MRO) on every medium -- fast"
TD="$OUT/community"; mkdir -p "$TD"; cp "$OUT/models"/*.xml "$TD/"
echo "   $(ls "$TD"/*.xml | wc -l) models"
for MED in ETO_SED ETO_SED_FERM ETO_METH; do
  awk -F'\t' -v m="$MED" 'NR==1{next} $1==m{n++} END{if(n+0==0) exit 1}' "$MDB" || { no "$MED not in media db, skipping"; continue; }
  note "$MED"
  ( cd "$TD" && "$CLI" -o "g_$MED" -m "$MED" --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | tail -6 | sed 's/^/     /'
  [ -s "$TD/g_${MED}_global.tsv" ] && column -t -s$'\t' "$TD/g_${MED}_global.tsv" | sed 's/^/     /'
done

if [ "$SKIP_DETAILED" = "1" ]; then note "SKIP_DETAILED=1"; say DONE; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0; fi

# ---------------------------------------------------------------- G. detailed
say "G. detailed on ETO_SED_FERM -- the ~3 h step, resumable"
note "ETO_SED_FERM, not ETO_SED: MTP4 eats glucose on ETO_SED, which is an annotation"
note "artefact, so ETO_SED would ask the methanogen the wrong question."
note "started $(date -Is)"
if [ -s "$TD/d_FERM_detailed.tsv" ]; then
  ok "already present -- delete $TD/d_FERM_detailed.tsv to force a re-run"
else
  ( cd "$TD" && "$CLI" -o d_FERM -m ETO_SED_FERM --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | tail -8 | sed 's/^/   /'
fi
DET="$TD/d_FERM_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET")-1)); echo "   rows: $N"
  if [ "$N" -gt 0 ]; then
    echo "   --- strongest 20 ---"
    { head -1 "$DET"; tail -n +2 "$DET" | sort -t$'\t' -k9,9gr | head -20; } | column -t -s$'\t' | sed 's/^/     /'
    echo "   --- donors ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    echo "   --- receivers ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    echo "   --- what MTP4 receives, and from whom ---"
    tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%s\n",$5,$4,$9}' | sort | head -25
    C=$(tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0' | wc -l)
    echo "   MTP4 as receiver: $C row(s)"
    [ "$C" -eq 0 ] && note "MTP4 receives nothing -- with a passing control, that is a result"
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
