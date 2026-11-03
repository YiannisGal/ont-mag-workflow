#!/usr/bin/env bash
# 43g_control_fixed.sh -- repair the control, diagnose mip=n/a, then run.
#
# 43f established, and these are not re-opened:
#   - smetana.legacy.Community merges correctly: 2338 reactions, 169 pool reactions,
#     objective R_Community_Growth, pools named R_EX_M_<cpd>_e_pool. Every claim I made
#     in 43c and 43d about a disconnected community was an artefact of my using
#     reframed's Community instead of SMETANA's.
#   - ETO_SED is 24 compounds by exact match. The media db is correct.
#   - MRO on the pair = 0.3. MIP = n/a, which is a separate defect, diagnosed in B here.
#
# 43f's CONTROL failed, and the fault was mine twice over:
#
#   1. The knockout deleted every reaction containing M_<cpd>_c that did not also contain
#      M_<cpd>_e. That includes reactions that CONSUME the compound -- and for any amino
#      acid the consumer is bio1. So I was deleting the biomass reaction. The knockout
#      then failed both the "dead alone" and "rescued" tests, for the wrong reason.
#   2. It reported "only 0 usable auxotroph(s)" without saying which condition each of
#      the 114 candidates failed. A control that fails uninformatively is barely better
#      than one that cannot fail. Section C prints the reason for every candidate.
#
# The repair: read bio1's OWN substrate list. Those are the compounds the organism must
# have. For each that also has an exchange and a transporter, delete only the reactions
# that PRODUCE it, protecting bio1, the transporter and the exchange.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43g"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
REL="$WS2/43_smetana/models_sed"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43g  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- A. biomass
say "A. what bio1 actually requires -- the objective is NOT the biomass reaction"
note "reframed reports the objective as R_EX_cpd11416_c, the biomass DRAIN. The"
note "reaction that assembles biomass is bio1, and its substrates are what the cell"
note "must obtain. Those are the only sensible auxotrophy targets."
"$PY" - "$REL" <<'PYEOF'
import sys, os, warnings
from reframed import load_cbmodel
rel = sys.argv[1]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_cbmodel(os.path.join(rel,'SED_das_semibin_003.xml'), flavor='fbc2')
exch = set(m.get_exchange_reactions())
obj  = [r for r,c in m.get_objective().items() if c]
print("   objective reaction        : %s  (exchange? %s)" % (obj, [o in exch for o in obj]))
bio = [r for r in m.reactions if r.lower().startswith('r_bio') and r not in exch]
print("   candidate biomass reactions: %s" % bio)
if bio:
    st = m.reactions[bio[0]].stoichiometry
    subs = [k for k,v in st.items() if v < 0]
    print("   %s has %d substrates, %d products" % (bio[0], len(subs), len(st)-len(subs)))
    have_ex = [k for k in subs if ('R_EX_%s' % k.replace('M_','').replace('_c','_e')) in exch]
    print("   substrates that also have an exchange reaction: %d" % len(have_ex))
    for k in have_ex[:12]:
        print("      %-18s %s" % (k, m.metabolites[k].name[:44]))
PYEOF

# ---------------------------------------------------------------- B. mip
say "B. why is MIP n/a while MRO is 0.3?"
note "MIP asks how many compounds the community can stop importing by trading"
note "internally. It needs a MINIMAL medium first, and returns None if the community"
note "cannot reach min_growth on one. This finds out which."
"$PY" - "$REL" "$MDB" <<'PYEOF'
import sys, os, csv, warnings, traceback
from reframed import load_cbmodel, Environment
from smetana.interface import Community
import smetana.interface as I
rel, mdb = sys.argv[1], sys.argv[2]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row.get('medium')=='ETO_SED': med.append(row['compound'].strip())
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel,n+'.xml'), flavor='fbc2')
ma, mb = load('LIQ_das_vamb_017'), load('SED_das_semibin_003')
c = Community('pc',[ma,mb])
env = Environment.from_compounds(med, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x))
from reframed import FBA
mm = c.merged.copy(); env.apply(mm, inplace=True, warning=False)
s = FBA(mm); g = s.fobj if s and s.fobj else 0.0
print("   community growth on ETO_SED            = %.6f" % g)
print("   (SMETANA's default min_growth is 0.1)")
if g < 0.1: print("   >>> below min_growth -- that alone explains mip=n/a")
else:       print("   >>> above min_growth, so the block is in the minimal-medium step")
try:
    r = I.minimal_environment(c, aerobic=False, min_mass_weight=False, min_growth=0.1,
                              max_uptake=10, validate=False, verbose=False)
    print("   minimal_environment -> %s compounds" % (len(r) if r else r))
except TypeError:
    try:
        r = I.minimal_environment(c, aerobic=False)
        print("   minimal_environment -> %s" % (len(r) if r else r))
    except Exception as e:
        print("   minimal_environment raised %s: %s" % (type(e).__name__, e))
except Exception as e:
    print("   minimal_environment raised %s: %s" % (type(e).__name__, e))
for mg in (0.1, 0.01, 0.001):
    try:
        v = I.mip_score(c, environment=env, min_growth=mg, verbose=False)
        print("   mip_score(min_growth=%-5s) -> %s" % (mg, v))
    except Exception as e:
        print("   mip_score(min_growth=%-5s) -> %s: %s" % (mg, type(e).__name__, e))
PYEOF

# ---------------------------------------------------------------- C. control
say "C. the control, rebuilt -- and every rejection explained"
"$PY" - "$REL" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, warnings, json
from reframed import load_cbmodel, save_cbmodel, Environment, FBA
rel, mdb, out = sys.argv[1], sys.argv[2], sys.argv[3]
SRC='SED_das_semibin_003'
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row.get('medium')=='ETO_SED': med.append(row['compound'].strip())

def load():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel,SRC+'.xml'), flavor='fbc2')

def grow(m, extra=()):
    mm = m.copy()
    Environment.from_compounds(list(med)+list(extra),
        fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(mm, inplace=True, warning=False)
    s = FBA(mm)
    return s.fobj if s and s.fobj else 0.0

base = load()
exch = set(base.get_exchange_reactions())
bio  = [r for r in base.reactions if r.lower().startswith('r_bio') and r not in exch]
if not bio:
    print("   no bio1-like reaction found"); print("CONTROL_STATUS: FAILED"); sys.exit(3)
BIO = bio[0]
print("   biomass reaction: %s" % BIO)
b0 = grow(base); print("   %s on ETO_SED: %.4f" % (SRC, b0))

subs = [k for k,v in base.reactions[BIO].stoichiometry.items() if v < 0]
targets=[]
for cm in subs:
    cpd = cm.replace('M_','')[:-2] if cm.endswith('_c') else None
    if not cpd: continue
    em = 'M_%s_e' % cpd
    ex = 'R_EX_%s_e' % cpd
    if cpd in med: continue                       # already supplied; no auxotrophy possible
    if em not in base.metabolites: continue       # no external form
    if ex not in exch: continue                   # cannot be fed
    trans = [r for r,rx in base.reactions.items()
             if cm in rx.stoichiometry and em in rx.stoichiometry]
    if not trans: continue                        # no way in
    targets.append((cpd, cm, em, ex, trans))
print("   biomass substrates that are feedable (exchange + transporter, not in medium): %d" % len(targets))

def producers(m, cm, protect):
    out=[]
    for r_id, rxn in m.reactions.items():
        if r_id in protect: continue
        st = rxn.stoichiometry
        if cm not in st: continue
        coef = st[cm]
        if (coef > 0 and rxn.ub > 0) or (coef < 0 and rxn.lb < 0):
            out.append(r_id)
    return out

rows=[]; found=[]
for cpd, cm, em, ex, trans in targets:
    if len(found) == 2: break
    m = load()
    protect = set(m.get_exchange_reactions()) | set(trans) | {BIO}
    kill = producers(m, cm, protect)
    if not kill:
        rows.append((cpd, base.metabolites[cm].name[:30], 0, None, None, "no producing reaction to remove")); continue
    for r in kill: m.remove_reaction(r)
    alone   = grow(m)
    rescued = grow(m, extra=[cpd])
    if alone >= 1e-6:
        why = "still grows alone (%.4f) -- alternative route remains" % alone
    elif rescued < 1e-6:
        why = "dead even when %s supplied -- knockout removed more than intended" % cpd
    else:
        why = "USABLE"
    rows.append((cpd, base.metabolites[cm].name[:30], len(kill), alone, rescued, why))
    if why == "USABLE":
        m.id = 'AUXO_%s' % cpd
        save_cbmodel(m, os.path.join(out,'AUXO_%s.xml' % cpd), flavor='fbc2')
        found.append(dict(cpd=cpd, name=base.metabolites[cm].name, removed=len(kill)))

print()
print("   %-12s %-32s %5s %9s %9s  %s" % ('compound','name','killed','alone','rescued','verdict'))
for cpd,name,k,a,r,why in rows:
    print("   %-12s %-32s %5d %9s %9s  %s" % (cpd, name, k,
          ('%.4f'%a) if a is not None else '-', ('%.4f'%r) if r is not None else '-', why))
print()
if len(found) < 2:
    print("   usable auxotrophs: %d of %d targets tried" % (len(found), len(rows)))
    print("CONTROL_STATUS: FAILED  no mutually-dependent pair -- reasons above")
    sys.exit(3)
json.dump(found, open(os.path.join(out,'control_pair.json'),'w'), indent=2)
print("   pair: AUXO_%s (%s) + AUXO_%s (%s)" % (found[0]['cpd'],found[0]['name'],
                                                found[1]['cpd'],found[1]['name']))
print("CONTROL_STATUS: BUILT")
PYEOF
CRC=$?
if [ "$CRC" -ne 0 ]; then
  bad "control still not built -- the table above says why for each candidate."
  note "sections D and E are NOT run. A zero from an unvalidated detector is not a result."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "control pair built"

say "D. does SMETANA see the cross-feeding it must see?"
KD="$OUT/control"; mkdir -p "$KD"; cp "$OUT"/AUXO_*.xml "$KD/"
( cd "$KD" && "$CLI" -o ctrl -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d AUXO_*.xml 2>&1 ) | sed 's/^/   /' | tail -12
CTRL="$KD/ctrl_detailed.tsv"; NROW=0
[ -s "$CTRL" ] && NROW=$(($(wc -l < "$CTRL") - 1))
echo "   ctrl_detailed.tsv rows: $NROW"
[ "$NROW" -gt 0 ] && { head -8 "$CTRL" | column -t -s$'\t' | sed 's/^/     /'; }
if [ "$NROW" -eq 0 ]; then
  bad "two obligate mutualists and SMETANA reports nothing. The detector is blind."
  note "every zero so far, 43e's included, is uninterpretable. Section E NOT run."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "control passed with $NROW row(s) -- from here a zero is about the organisms"

say "E. the in-situ ten on ETO_SED"
TD="$OUT/ten"; mkdir -p "$TD"; cp "$REL"/*.xml "$TD/"
note "$(ls "$TD"/*.xml | wc -l) models, global first"
( cd "$TD" && "$CLI" -o ten -m ETO_SED --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | sed 's/^/   /' | tail -10
[ -s "$TD/ten_global.tsv" ] && { echo "   --- ten_global.tsv ---"; column -t -s$'\t' "$TD/ten_global.tsv" | sed 's/^/     /'; }
note "now detailed. This is the slow one."
( cd "$TD" && "$CLI" -o ten -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | sed 's/^/   /' | tail -10
DET="$TD/ten_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET") - 1)); echo "   ten_detailed.tsv rows: $N"
  if [ "$N" -gt 0 ]; then
    echo "   --- strongest 15 ---"
    { head -1 "$DET"; tail -n +2 "$DET" | sort -t$'\t' -k9,9gr | head -15; } | column -t -s$'\t' | sed 's/^/     /'
    echo "   --- donors ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(d in c) printf "     %-26s %d\n",d,c[d]}' | sort -k2,2nr
    echo "   --- receivers ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(d in c) printf "     %-26s %d\n",d,c[d]}' | sort -k2,2nr
  else
    note "zero rows -- and now that is a statement about the community"
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
