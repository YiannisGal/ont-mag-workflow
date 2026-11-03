#!/usr/bin/env bash
# 43m_mcr_audit_and_run.sh -- audit methanogenesis, declare the constraints, then run.
#
# RUN IN tmux. Section F is the ~3 h detailed run and is resumable.
#
# WHAT 43l's GATE CAUGHT -- correctly
#     LIQ_das_vamb_003  ETO_METH  growth 0.0110  methane 2.2701   proper methanogen
#     LIQ_das_vamb_017  ETO_METH  growth 0.0845  methane 0.0000   grows, no methane
#
#   LIQ_das_vamb_017 grows on H2/CO2 without methanogenesis. 43j already reported
#   "methyl-CoM reductase reactions: NONE FOUND BY NAME" for it and I did not follow it
#   up. Its reconstruction has no mcr on any medium. The 0.973 methane flux 43i saw came
#   from some other route.
#
#   KOfam found mcr in 8/8 archaea (41f). So the gene is in the genome and gapseq left
#   the reaction out. That is a measurable instance of the archaeal annotation gap, and
#   section B checks it against the KOfam evidence rather than asserting it.
#
# THE SECOND PROBLEM, and the declared fix
#   MTP4 takes up glucose at -10.0 on ETO_SED AND on ETO_SED_FERM (which is ETO_SED plus
#   fermentation products, so the sugars are still present). Only ETO_METH keeps it
#   honest, and ETO_METH starves the fermenters. The community needs sugars; MTP4's model
#   eats them.
#
#   Section C applies a DECLARED CONSTRAINT: sugar uptake blocked in the archaea, on the
#   grounds that Methanosarcinaceae and Methanomicrobiales are not saccharolytic. Uptake
#   only -- secretion stays open. This is a modelling decision, it goes in the methods in
#   one sentence, and section C prints its effect so it is auditable. Silently allowing a
#   methanogen to ferment glucose is the alternative, and it is worse.
#
# THE GATE IS NOW CONDITIONAL, not relaxed
#   mcr present in the model  -> methane flux on ETO_METH is REQUIRED
#   mcr absent                -> reported as a reconstruction gap; the organism may stay
#                                in the community, but NO methanogenesis claim rests on it
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43m"; mkdir -p "$OUT/models"
export TMPDIR="${SCRATCH_ROOT}/tmp43m"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
SRC="$WS2/43_smetana/run_43l/models"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
SKIP_DETAILED="${SKIP_DETAILED:-0}"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43m  $(date -Is)"; echo "out: $OUT"
[ -d "$SRC" ] || { bad "no relabelled models at $SRC -- run 43l first"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. mcr in models
say "A. is methyl-CoM reductase in either archaeal model?"
note "searched by reaction NAME, by EC 2.8.4.1, and by the CoM/CoB metabolites the"
note "reaction must consume. Three independent handles, because one name can miss."
"$PY" - "$SRC" "$MTP4" "$HYD" "$OUT" <<'PYEOF'
import sys, os, re, json, warnings
from reframed import load_cbmodel
src, mtp4, hyd, out = sys.argv[1:5]
res={}
for g in (mtp4, hyd):
    p=os.path.join(src, g+'.xml')
    if not os.path.exists(p): print("   %s absent" % g); continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    txt=open(p, encoding='utf-8', errors='replace').read()
    byname=[(r,rx.name) for r,rx in m.reactions.items()
            if re.search(r'methyl-?coenzyme\s*m\s*reductase|coenzyme-?B.*coenzyme-?M|\bmcr[ABG]?\b',
                         (rx.name or ''), re.I)]
    byec = re.findall(r'2\.8\.4\.1', txt)
    coms = {mid: met.name for mid, met in m.metabolites.items()
            if re.search(r'coenzyme\s*m|coenzyme\s*b|CoM-|CoB-|methyl-?com', (met.name or ''), re.I)}
    ch4  = [mid for mid, met in m.metabolites.items() if re.match(r'^methane', (met.name or ''), re.I)]
    # any reaction producing methane
    prod=[]
    for r, rx in m.reactions.items():
        for k in ch4:
            v = rx.stoichiometry.get(k, 0)
            if (v > 0 and rx.ub > 0) or (v < 0 and rx.lb < 0):
                prod.append((r, rx.name)); break
    print("\n   [%s]  %d reactions" % (g, len(m.reactions)))
    print("     by name (mcr / methyl-CoM reductase) : %s" % (byname or 'NONE'))
    print("     EC 2.8.4.1 occurrences in the SBML   : %d" % len(byec))
    print("     CoM / CoB metabolites present        : %d %s" % (len(coms), list(coms.values())[:3]))
    print("     methane metabolites                  : %s" % ch4)
    print("     reactions that can PRODUCE methane   : %d" % len(prod))
    for r, nme in prod[:6]: print("        %-22s %s" % (r, nme))
    has = bool(byname) or bool(byec) or bool(prod)
    res[g]=dict(byname=[r for r,_ in byname], ec=len(byec), n_com=len(coms),
                producers=[r for r,_ in prod], has_methanogenesis=has)
    print("     -> methanogenesis representable: %s" % has)
json.dump(res, open(os.path.join(out,'mcr_audit.json'),'w'), indent=2)
PYEOF

# ---------------------------------------------------------------- B. kofam
say "B. does the GENOME carry mcr, whatever the model says?"
note "the KO numbers are looked up in ko_list by description, never from memory."
KOF="$WS2/11_kofam"
if [ -d "$KOF" ]; then
  KOLIST=$(find "$KOF" -maxdepth 2 -name ko_list | head -1)
  if [ -n "$KOLIST" ]; then
    note "KOs whose description mentions methyl-coenzyme M reductase:"
    awk -F'\t' 'tolower($NF) ~ /methyl-coenzyme m reductase/ {printf "     %-10s %s\n", $1, $NF}' "$KOLIST" | head -12
    KOS=$(awk -F'\t' 'tolower($NF) ~ /methyl-coenzyme m reductase/ {print $1}' "$KOLIST" | tr '\n' '|' | sed 's/|$//')
    if [ -n "$KOS" ]; then
      note "hits for those KOs in each archaeal MAG:"
      for g in "$MTP4" "$HYD"; do
        f=$(find "$KOF" -type f \( -name "*${g}*" \) 2>/dev/null | grep -viE 'ko_list|profiles' | head -1)
        if [ -n "$f" ]; then
          c=$(grep -cE "\b(${KOS})\b" "$f" 2>/dev/null || echo 0)
          printf '     %-24s %-52s %s hit(s)\n' "$g" "$(basename "$f")" "$c"
          grep -E "\b(${KOS})\b" "$f" 2>/dev/null | head -4 | cut -c1-140 | sed 's/^/        /'
        else
          printf '     %-24s no KOfam result file found\n' "$g"
        fi
      done
    fi
  else no "no ko_list under $KOF"; fi
else no "no $KOF -- cannot cross-check the genome evidence"; fi
note "if the genome has mcr and the model does not, that is gapseq's archaeal gap,"
note "and it is a result for ARCHAEAL_ANNOTATION_GAP.md, not something to patch here."

# ---------------------------------------------------------------- C. constraint
say "C. DECLARED CONSTRAINT -- archaea may not take up sugars"
note "Methanosarcinaceae and Methanomicrobiales are not saccharolytic. Uptake of"
note "glucose, xylose, cellobiose and galacturonate is blocked in the two archaeal"
note "models only. Secretion is left open. Bacteria are untouched."
note "This is a modelling decision and belongs in the methods in one sentence."
cp "$SRC"/*.xml "$OUT/models/"
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
    if not os.path.exists(p): continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    before={}
    for med in ('ETO_SED','ETO_SED_FERM','ETO_METH'):
        if med not in media: continue
        t=m.copy(); Environment.from_compounds(media[med], fmt_func=fmt).apply(t,inplace=True,warning=False)
        s=FBA(t); before[med]=s.fobj if s and s.fobj else 0.0
    blocked=[]
    for cpd,nm in SUGAR.items():
        r=fmt(cpd)
        if r in m.reactions:
            m.reactions[r].lb=0.0          # no uptake; secretion (ub) untouched
            blocked.append(nm)
    save_cbmodel(m, p, flavor='fbc2')
    after={}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m2=load_cbmodel(p, flavor='fbc2')
    for med in ('ETO_SED','ETO_SED_FERM','ETO_METH'):
        if med not in media: continue
        t=m2.copy(); Environment.from_compounds(media[med], fmt_func=fmt).apply(t,inplace=True,warning=False)
        s=FBA(t); after[med]=s.fobj if s and s.fobj else 0.0
    print("\n   [%s] blocked uptake of: %s" % (g, ', '.join(blocked) or 'none present'))
    for med in before:
        print("     %-14s growth %.4f -> %.4f" % (med, before[med], after[med]))
PYEOF

# ---------------------------------------------------------------- D. gate
say "D. CONDITIONAL GATE -- methane required only where mcr is representable"
"$PY" - "$OUT/models" "$MDB" "$MTP4" "$HYD" "$OUT" <<'PYEOF'
import sys, os, csv, re, json, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, mdb, mtp4, hyd, out = sys.argv[1:6]
audit=json.load(open(os.path.join(out,'mcr_audit.json')))
media={}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'],[]).append(row['compound'].strip())
SUGAR={'cpd00027','cpd00154','cpd00158','cpd00280'}
fmt=lambda x: "R_EX_{}_e".format(x)
fails=0
for g in (mtp4, hyd):
    p=os.path.join(mdir,g+'.xml')
    if not os.path.exists(p): continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    can = audit.get(g,{}).get('has_methanogenesis', False)
    print("\n   [%s]  mcr representable in the model: %s" % (g, can))
    for med in ('ETO_METH','ETO_SED','ETO_SED_FERM'):
        if med not in media: continue
        t=m.copy(); Environment.from_compounds(media[med], fmt_func=fmt).apply(t,inplace=True,warning=False)
        s=FBA(t); gr=s.fobj if s and s.fobj else 0.0
        exch=set(t.get_exchange_reactions())
        rows=[(r,v) for r,v in (s.values.items() if s and s.values else []) if r in exch and abs(v)>1e-6]
        ch4=[v for r,v in rows if any(re.match(r'^methane',(t.metabolites[k].name or ''),re.I) for k in t.reactions[r].stoichiometry)]
        sug=[round(v,2) for r,v in rows if v<0 and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in t.reactions[r].stoichiometry)]
        print("     %-14s growth=%-9.4f methane=%-9.4f sugar_uptake=%s" % (med, gr, max(ch4) if ch4 else 0.0, sug or 'none'))
        if med=='ETO_METH' and can and (gr<=1e-6 or not ch4):
            print("       [FAIL] mcr is present but no methane on a methanogenic medium"); fails+=1
        if sug:
            print("       [FAIL] sugar uptake survived the declared constraint"); fails+=1
print()
print("GATE_STATUS: %s (%d failure(s))" % ('PASSED' if fails==0 else 'FAILED', fails))
print("NOTE: an archaeon without mcr is NOT failed here, but no methanogenesis claim")
print("      may rest on it. That is recorded, not silently passed.")
sys.exit(0 if fails==0 else 3)
PYEOF
[ $? -eq 0 ] || { bad "gate failed -- see above"; echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1; }
ok "gate passed"

# ---------------------------------------------------------------- E. control + global
say "E. control, then global scores"
CD="$OUT/control"; mkdir -p "$CD"
if ls "$WS2/43_smetana/run_43g"/AUXO_*.xml >/dev/null 2>&1; then
  cp "$WS2/43_smetana/run_43g"/AUXO_*.xml "$CD/"
  ( cd "$CD" && "$CLI" -o ctrl -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d AUXO_*.xml 2>&1 ) | tail -5 | sed 's/^/   /'
  R=0; [ -s "$CD/ctrl_detailed.tsv" ] && R=$(($(wc -l < "$CD/ctrl_detailed.tsv")-1))
  echo "   control rows: $R"
  [ "$R" -gt 0 ] || { bad "control gives zero -- nothing below would be interpretable"; echo "STEP_STATUS: FAILED"; exit 1; }
  ok "detector confirmed"
else no "no control pair available -- zeros below are NOT interpretable"; fi

TD="$OUT/community"; mkdir -p "$TD"; cp "$OUT/models"/*.xml "$TD/"
echo "   $(ls "$TD"/*.xml | wc -l) models"
for MED in ETO_SED ETO_SED_FERM ETO_METH; do
  awk -F'\t' -v m="$MED" 'NR==1{next} $1==m{n++} END{exit (n+0==0)}' "$MDB" || { no "$MED absent"; continue; }
  note "$MED"
  ( cd "$TD" && "$CLI" -o "g_$MED" -m "$MED" --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | tail -5 | sed 's/^/     /'
  [ -s "$TD/g_${MED}_global.tsv" ] && column -t -s$'\t' "$TD/g_${MED}_global.tsv" | sed 's/^/     /'
done
[ "$SKIP_DETAILED" = "1" ] && { note "SKIP_DETAILED=1"; say DONE; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0; }

# ---------------------------------------------------------------- F. detailed
say "F. detailed on ETO_SED_FERM -- the ~3 h step, resumable"
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
    echo "   --- what MTP4 receives ---"
    tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0{printf "     %-16s from %-26s smetana=%s\n",$5,$4,$9}' | head -25
    C=$(tail -n +2 "$DET" | awk -F'\t' -v a="$MTP4" '$3==a && $9>0' | wc -l)
    echo "   MTP4 as receiver: $C row(s)"
    [ "$C" -eq 0 ] && note "with the control passing, MTP4 receiving nothing is a RESULT"
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
