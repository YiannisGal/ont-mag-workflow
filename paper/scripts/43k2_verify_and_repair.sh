#!/usr/bin/env bash
# 43k2_verify_and_repair.sh -- the model was fine; my gate was wrong.
#
# 43k's own gapseq output for MTP4 (LIQ_das_vamb_003):
#     Added reactions:   19
#     Final growth rate: 0.01096382
#     Uptake at limit:   H2:10
#     Top produced:      H2O:5.358, Methane:2.27, H+:0.061, ...
# It takes up H2 at the limit and makes methane. It is a methanogen.
#
# My section F then reported growth 0.0000 and failed it. THE GATE WAS WRONG:
# the new model is raw gapseq output with `_e0` exchanges, and I built the medium with
# `R_EX_{}_e`, the RELABELLED convention used in 43_smetana/models*. The model was fed
# nothing. I hardcoded a naming convention instead of reading it off the model -- the
# same error as assuming SMETANA's pool names in 43c.
#
# Also wrong in 43k:
#   - section E passed `-rxnXgenes.RDS` to `gapseq fill`. My find used `-o` without
#     grouping so head -1 took the wrong file. A data.table is not a model, hence
#     "no applicable method for `@` applied to an object of class data.table".
#   - `cmd | tail || fallback` cannot work: the pipe's exit status is tail's, always 0.
#     So the fallback never ran and gapseq's own failures were invisible.
#
# TO CHECK, NOT ASSUME: the banner said `Sequence DB: b470543 (2026-05-29, Bacteria)`.
# If find ran against the bacterial database for an archaeon, 1125 reactions is an
# undercount. Section A looks for evidence either way and does NOT silently re-run
# a multi-hour step on a guess.
#
# Nothing here takes hours. Section E is the only gapseq call and it is a fill.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43k_archaea"; mkdir -p "$OUT/models"
export TMPDIR="${SCRATCH_ROOT}/tmp43k"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"
GAPSEQ="$WS2/envs/gapseq_local/gapseq"
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
export PATH="$ENVDIR/bin:$PATH"
METHCSV="$WS2/41_env_checks/ETO_METH-medium.csv"
SEDCSV="$WS2/41_env_checks/ETO_SED-medium.csv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
echo "43k2  $(date -Is)"

# ---------------------------------------------------------------- A. what was built
say "A. what 43k actually produced, and against which database"
MD="$OUT/$MTP4"
[ -d "$MD" ] || { bad "no $MD"; echo "STEP_STATUS: FAILED"; exit 1; }
ls -la "$MD" | sed 's/^/     /'
echo
note "pathway/reaction tables and their sizes:"
for f in "$MD"/*.tbl; do [ -s "$f" ] && printf '     %-52s %8s lines\n' "$(basename "$f")" "$(wc -l < "$f")"; done
echo
note "evidence of which sequence database was used:"
for f in "$MD"/*Pathways.tbl "$MD"/*Reactions.tbl; do
  [ -s "$f" ] || continue
  echo "     $(basename "$f") header:"; head -1 "$f" | cut -c1-200 | sed 's/^/       /'
done
grep -rilE 'archae' "$MD" 2>/dev/null | head -5 | sed 's/^/     mentions archaea: /'
note "gapseq's installed databases:"
ls -d "$WS2"/envs/gapseq_local/dat/seq/* 2>/dev/null | sed 's/^/     /' || \
  note "(no dat/seq under the local checkout)"
note "if only Bacteria is present, MTP4's annotation is bounded by that and the"
note "reaction count is a floor, not an estimate. That is a stated limitation, not"
note "a reason to re-run four hours on a guess."

# ---------------------------------------------------------------- B. the real gate
say "B. the gate, re-run with the exchange suffix READ FROM THE MODEL"
"$PY" - "$MD" "$METHCSV" "$SEDCSV" "$OUT" <<'PYEOF'
import sys, os, csv, glob, re, json, warnings
from reframed import load_cbmodel, Environment, FBA
md, methcsv, sedcsv, out = sys.argv[1:5]

def cpds(path):
    with open(path, newline='') as fh:
        rd = csv.DictReader(fh)
        idc = next((c for c in rd.fieldnames if c.lower() in ('compounds','compound','id','cpd')), rd.fieldnames[0])
        return [r[idc].strip() for r in rd]
meth, sed = cpds(methcsv), cpds(sedcsv)
SUGAR = {'cpd00027','cpd00154','cpd00158','cpd00280'}

def suffix(m):
    """read the exchange naming off the model instead of assuming it"""
    for r in m.get_exchange_reactions():
        mt = re.match(r'R_EX_(cpd\d+)_([A-Za-z]\d*)$', r)
        if mt: return mt.group(2)
    return None

cand = [p for p in glob.glob(os.path.join(md,'*.xml')) if 'draft' not in os.path.basename(p)]
if not cand:
    print("   no non-draft .xml in %s" % md); sys.exit(3)
p = cand[0]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_cbmodel(p, flavor='fbc2')
sfx = suffix(m)
print("   model            : %s" % os.path.basename(p))
print("   reactions        : %d" % len(m.reactions))
print("   exchanges        : %d" % len(m.get_exchange_reactions()))
print("   exchange suffix  : %r   <-- 43k assumed 'e' and fed the model nothing" % sfx)
if sfx is None:
    print("   cannot determine exchange naming"); sys.exit(3)
fmt = lambda x: "R_EX_{}_{}".format(x, sfx)
print("   example          : %s" % fmt('cpd00009'))
present = [c for c in meth if fmt(c) in set(m.get_exchange_reactions())]
print("   ETO_METH compounds the model can take up: %d of %d" % (len(present), len(meth)))
missing = [c for c in meth if fmt(c) not in set(m.get_exchange_reactions())]
if missing: print("     no exchange for: %s" % ', '.join(missing))

res = {}
for label, cs in (('ETO_METH', meth), ('ETO_SED', sed)):
    mm = m.copy()
    Environment.from_compounds(cs, fmt_func=fmt).apply(mm, inplace=True, warning=False)
    s = FBA(mm); g = s.fobj if s and s.fobj else 0.0
    print("\n   --- %s ---" % label)
    print("   growth: %.6f" % g)
    res[label] = g
    if not (s and s.values): continue
    exch = set(mm.get_exchange_reactions())
    rows = [(r,v) for r,v in s.values.items() if r in exch and abs(v) > 1e-6]
    def nm(r):
        st = mm.reactions[r].stoichiometry
        k = list(st)[0] if st else None
        return mm.metabolites[k].name if k in mm.metabolites else r
    up = sorted([x for x in rows if x[1] < 0], key=lambda x: x[1])
    pr = sorted([x for x in rows if x[1] > 0], key=lambda x: -x[1])
    print("   consumed: %s" % ', '.join('%s %.2f' % (nm(r), v) for r,v in up[:8]))
    print("   produced: %s" % ', '.join('%s %.2f' % (nm(r), v) for r,v in pr[:8]))
    ch4 = [v for r,v in rows if any(re.match(r'^methane', (mm.metabolites[k].name or ''), re.I)
                                     for k in mm.reactions[r].stoichiometry)]
    sug = [(nm(r), v) for r,v in rows if v < 0 and any(
             k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in mm.reactions[r].stoichiometry)]
    print("   methane flux : %s" % (max(ch4) if ch4 else 0.0))
    print("   sugar uptake : %s" % (sug or 'none'))
    res[label+'_ch4'] = max(ch4) if ch4 else 0.0
    res[label+'_sugar'] = sug

json.dump(res, open(os.path.join(out,'mtp4_gate.json'),'w'), indent=2, default=str)
print()
okg = res.get('ETO_METH',0) > 1e-6
okc = res.get('ETO_METH_ch4',0) > 1e-6
oks = not res.get('ETO_METH_sugar')
print("   1. grows on ETO_METH     : %.6f  %s" % (res.get('ETO_METH',0), 'OK' if okg else 'FAIL'))
print("   2. methane flux          : %.6f  %s" % (res.get('ETO_METH_ch4',0), 'OK' if okc else 'FAIL'))
print("   3. no sugar on ETO_METH  : %-9s %s" % (res.get('ETO_METH_sugar') or 'none', 'OK' if oks else 'FAIL'))
print("MTP4_GATE: %s" % ('PASSED' if (okg and okc and oks) else 'FAILED'))
sys.exit(0 if (okg and okc and oks) else 3)
PYEOF
MRC=$?
[ "$MRC" -eq 0 ] && ok "MTP4 passes on its own terms" || no "MTP4 gate failed -- see above"

# ---------------------------------------------------------------- C. the right RDS
say "C. the draft RDS for $HYD -- the file 43k should have picked"
SRC="$WS2/42_gems/$HYD"
ls -la "$SRC"/*.RDS 2>/dev/null | sed 's/^/     /'
DRAFT=$(ls "$SRC"/*-draft.RDS 2>/dev/null | head -1)
if [ -z "$DRAFT" ]; then
  DRAFT=$(ls "$SRC"/*.RDS 2>/dev/null | grep -viE 'rxnWeights|rxnXgenes' | head -1)
fi
[ -n "$DRAFT" ] && ok "draft: $(basename "$DRAFT")" || no "no draft RDS distinguishable from the data tables"

# ---------------------------------------------------------------- D. refill hyd
say "D. refill $HYD on ETO_METH, with the correct inputs"
HD="$OUT/$HYD"; mkdir -p "$HD"
if ls "$HD"/*.xml 2>/dev/null | grep -qv draft; then
  ok "already refilled"
elif [ -n "$DRAFT" ]; then
  WT=$(ls "$SRC"/*-rxnWeights.RDS 2>/dev/null | head -1)
  XG=$(ls "$SRC"/*-rxnXgenes.RDS 2>/dev/null | head -1)
  note "gapseq fill -m $(basename "$DRAFT") -n ETO_METH ${WT:+-c $(basename "$WT")} ${XG:+-g $(basename "$XG")}"
  ( cd "$HD" && "$GAPSEQ" fill -m "$DRAFT" -n "$METHCSV" ${WT:+-c "$WT"} ${XG:+-g "$XG"} ) 2>&1 | tail -14 | sed 's/^/     /'
  RC=${PIPESTATUS[0]:-0}          # exit status of gapseq, NOT of tail
  [ "$RC" -eq 0 ] && ok "fill exited 0" || no "fill exited $RC"
else
  no "no usable draft -- skipping"
fi

# ---------------------------------------------------------------- E. both
say "E. assemble and gate BOTH archaea"
for d in "$OUT/$MTP4" "$OUT/$HYD"; do
  m=$(ls "$d"/*.xml 2>/dev/null | grep -v draft | head -1)
  [ -n "$m" ] && cp "$m" "$OUT/models/$(basename "$d").xml"
done
ls -l "$OUT/models" | sed 's/^/     /'
"$PY" - "$OUT/models" "$METHCSV" "$SEDCSV" <<'PYEOF'
import sys, os, csv, glob, re, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, methcsv, sedcsv = sys.argv[1:4]
def cpds(p):
    with open(p, newline='') as fh:
        rd=csv.DictReader(fh)
        idc=next((c for c in rd.fieldnames if c.lower() in ('compounds','compound','id','cpd')), rd.fieldnames[0])
        return [r[idc].strip() for r in rd]
meth, sed = cpds(methcsv), cpds(sedcsv)
SUGAR={'cpd00027','cpd00154','cpd00158','cpd00280'}
def suffix(m):
    for r in m.get_exchange_reactions():
        mt=re.match(r'R_EX_(cpd\d+)_([A-Za-z]\d*)$', r)
        if mt: return mt.group(2)
    return 'e'
fails=0; seen=0
for p in sorted(glob.glob(os.path.join(mdir,'*.xml'))):
    seen+=1
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    sfx=suffix(m); fmt=lambda x: "R_EX_{}_{}".format(x,sfx)
    mm=m.copy()
    Environment.from_compounds(meth, fmt_func=fmt).apply(mm,inplace=True,warning=False)
    s=FBA(mm); g=s.fobj if s and s.fobj else 0.0
    exch=set(mm.get_exchange_reactions())
    rows=[(r,v) for r,v in (s.values.items() if s and s.values else []) if r in exch and abs(v)>1e-6]
    ch4=[v for r,v in rows if any(re.match(r'^methane',(mm.metabolites[k].name or ''),re.I) for k in mm.reactions[r].stoichiometry)]
    sug=[r for r,v in rows if v<0 and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in mm.reactions[r].stoichiometry)]
    c1,c2,c3 = g>1e-6, (max(ch4) if ch4 else 0)>1e-6, not sug
    print("   %-22s suffix=%-3s rxn=%-5d  growth=%-9.6f methane=%-8.4f sugar=%s"
          % (n, sfx, len(m.reactions), g, max(ch4) if ch4 else 0.0, sug or 'none'))
    print("     %s grows   %s methanogenic   %s sugar-free"
          % ('OK  ' if c1 else 'FAIL', 'OK  ' if c2 else 'FAIL', 'OK  ' if c3 else 'FAIL'))
    if not (c1 and c2 and c3): fails+=1
print()
print("GATE_STATUS: %s  (%d of %d model(s) failed)" % ('PASSED' if fails==0 and seen>0 else 'FAILED', fails, seen))
sys.exit(0 if (fails==0 and seen>0) else 3)
PYEOF
GRC=$?
if [ "$GRC" -eq 0 ]; then
  ok "both archaea grow on ETO_METH, make methane, eat no sugar"
  note "43l can now assemble the eleven-member community from $OUT/models plus the"
  note "nine bacterial models in 43_smetana/models"
  echo; echo "STEP_STATUS: COMPLETE  $(date -Is)"
else
  no "not all archaea passed -- 43l should not run yet"
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"
fi
