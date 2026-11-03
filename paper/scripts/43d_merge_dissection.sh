#!/usr/bin/env bash
# 43d_merge_dissection.sh -- why does the merge produce ZERO pool reactions?
#
# 43c settled what the problem is NOT. The medium is fine, the members grow (10/10 on
# ETO_SED), and the biomass reaction is a red herring -- forcing R_bio1 changed nothing.
# The finding is one line:
#
#       merged: 2073 reactions, 0 pool reactions
#
# A community with no pool is two models in one file. Nothing can pass between them, so
# every score is zero by construction and section F's 170 opened boundary reactions had
# nothing to open onto.
#
# Live hypotheses, all about the MERGE:
#
#   H3 RELABEL.  41e's `_e0`->`_e` sed pass renamed metabolite IDs but may have left the
#                SBML compartment declared as e0. My 41e control compared reaction counts
#                and single-model FBA -- neither touches compartment metadata, so it could
#                not have seen this. (Rule 51.)
#   H4 WRONG OBJ. SMETANA uses `community.merged` (interface.py:116). I used
#                `community.merged_model`. If they differ, or if the merged model's
#                objective is never set, FBA returns zero whatever the pool looks like.
#   H5 API.      reframed needs the external compartment identified explicitly and finds
#                nothing for ModelSEED-flavoured SBML.
#
# This script MEASURES and runs one narrow experiment per hypothesis. It does not rewrite
# any model in place. Section G runs SMETANA's OWN entry points -- not my reconstruction
# of them -- because doing otherwise is the mistake that produced 41e and again 43c.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/diag_43d"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
echo "43d  $(date -Is)"; echo "out: $OUT"
PY="$WS2/envs/smetana/bin/python"
[ -x "$PY" ] || { bad "no smetana python at $PY"; echo "STEP_STATUS: FAILED"; exit 1; }
export PYTHONWARNINGS="ignore::UserWarning"      # the biomass warning is settled; stop the noise

REL="$WS2/43_smetana/models_sed"                 # relabelled + refilled on ETO_SED
[ -d "$REL" ] || { bad "no $REL"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. the merge, read
say "A. Community.__init__ in full -- what builds the pool, and from what"
"$PY" - <<'PYEOF'
import inspect
import reframed.community.model as C
src = inspect.getsource(C.Community.__init__)
for i,l in enumerate(src.split('\n')): print("   %3d | %s" % (i+1, l.rstrip()))
print()
print("   --- the merge routine ---")
fn = None
for name in ('merge_models','_merge_models','build','_build'):
    if hasattr(C.Community, name): fn = getattr(C.Community, name); print("   found Community.%s" % name); break
if fn is None:
    print("   no obvious merge method; Community attributes:")
    print("   ", [a for a in dir(C.Community) if not a.startswith('__')])
else:
    body = inspect.getsource(fn).split('\n')
    for i,l in enumerate(body):
        if any(k in l for k in ('external','pool','compartment','objective','biomass','_e_','extracellular')):
            print("   %3d | %s" % (i+1, l.rstrip()))
PYEOF

say "B. how reframed decides a compartment is external"
"$PY" - <<'PYEOF'
import inspect, re
from reframed.core.model import Compartment
print("   Compartment signature:", inspect.signature(Compartment.__init__))
import reframed.io.sbml as S
src = inspect.getsource(S).split('\n')
for i,l in enumerate(src):
    if 'external' in l.lower(): print("   %4d | %s" % (i+1, l.rstrip()))
PYEOF

# ---------------------------------------------------------------- C. compartments
say "C. H3 -- do the relabelled models still declare an external compartment?"
note "comparing a relabelled model against an ORIGINAL gapseq model. This compares"
note "STRUCTURE, not biology: the originals were filled on ETO_FILL and are not"
note "interchangeable with the refilled set for any result."
ORIG=$(find "$WS2/42_gems" -name '*.xml' ! -name '*-draft.xml' 2>/dev/null | head -1)
if [ -z "$ORIG" ]; then no "no original model found under 42_gems -- H3 cannot be tested"; fi
"$PY" - "$REL" "${ORIG:-NONE}" "$OUT" <<'PYEOF'
import sys, os, glob
from reframed import load_cbmodel
rel, orig, out = sys.argv[1], sys.argv[2], sys.argv[3]

def report(path, tag):
    print("\n   [%s] %s" % (tag, os.path.basename(path)))
    try:
        m = load_cbmodel(path, flavor='fbc2')
    except Exception as e:
        print("      load failed: %s" % e); return None
    print("      compartments:")
    for cid, c in m.compartments.items():
        ext = getattr(c, 'external', None)
        print("        %-10s name=%-22s external=%s" % (cid, getattr(c,'name','?'), ext))
    n_ext = sum(1 for c in m.compartments.values() if getattr(c,'external',False))
    print("      compartments flagged external: %d   <-- 0 means no pool can be built" % n_ext)
    # where do metabolites actually live
    from collections import Counter
    cnt = Counter(met.compartment for met in m.metabolites.values())
    print("      metabolites per compartment: %s" % dict(cnt))
    ex = m.get_exchange_reactions()
    print("      exchange reactions: %d  e.g. %s" % (len(ex), list(ex)[:3]))
    if ex:
        r = m.reactions[list(ex)[0]]
        print("      first exchange stoichiometry: %s" % r.stoichiometry)
        for mid in r.stoichiometry:
            print("        metabolite %s lives in compartment %r" % (mid, m.metabolites[mid].compartment))
    return m

a = report(os.path.join(rel, sorted(os.listdir(rel))[0]) if os.path.isdir(rel) else rel, "RELABELLED")
if orig != 'NONE':
    b = report(orig, "ORIGINAL")
PYEOF

# ---------------------------------------------------------------- D. merged obj
say "D. H4 -- is 'merged' the same object as 'merged_model', and is its objective set?"
"$PY" - "$REL" <<'PYEOF'
import sys, os, warnings
from reframed import load_cbmodel
from reframed.community.model import Community
rel = sys.argv[1]
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel, n+'.xml'), flavor='fbc2')
ma, mb = load('LIQ_das_vamb_017'), load('SED_das_semibin_003')
comm = Community('d', [ma, mb])
print("   Community attributes:", [a for a in dir(comm) if not a.startswith('_')])
mm  = getattr(comm, 'merged_model', None)
mg  = getattr(comm, 'merged', None)
print("   merged_model is None? %s   merged is None? %s" % (mm is None, mg is None))
if mm is not None and mg is not None:
    print("   merged IS merged_model? %s" % (mm is mg))
    print("   merged_model reactions=%d   merged reactions=%d" % (len(mm.reactions), len(mg.reactions)))
for tag, model in (('merged_model', mm), ('merged', mg)):
    if model is None: continue
    obj = {r:c for r,c in model.get_objective().items() if c}
    pool = [r for r in model.reactions if r.endswith('_pool')]
    cb   = [r for r in model.reactions if 'community_biomass' in r]
    print("   [%s] objective=%s  pool_rxns=%d  community_biomass_rxns=%s"
          % (tag, obj if obj else 'EMPTY  <-- FBA returns 0 with no objective', len(pool), cb))
    print("   [%s] compartments: %s" % (tag, list(model.compartments.keys())))
PYEOF

# ---------------------------------------------------------------- E. the experiment
say "E. the experiment -- flag the external compartment, then re-merge"
note "this is the H3/H5 test. If declaring 'e' external makes pool reactions appear,"
note "the merge was never given an external compartment to pool on."
note "It is run on COPIES in memory. Nothing on disk is modified."
"$PY" - "$REL" "$OUT" <<'PYEOF'
import sys, os, warnings, json
from reframed import load_cbmodel
from reframed.community.model import Community
rel, out = sys.argv[1], sys.argv[2]
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel, n+'.xml'), flavor='fbc2')

def pools(models, label):
    comm = Community('x', models)
    res = {}
    for tag in ('merged_model','merged'):
        m = getattr(comm, tag, None)
        if m is None: continue
        res[tag] = len([r for r in m.reactions if r.endswith('_pool')])
    print("   %-46s %s" % (label, res))
    return comm, res

print()
A, B = 'LIQ_das_vamb_017', 'SED_das_semibin_003'
ma, mb = load(A), load(B)
_, base = pools([ma, mb], "as loaded")

# which compartment do the exchange metabolites actually sit in?
ex = list(ma.get_exchange_reactions())
extcomp = None
if ex:
    for mid in ma.reactions[ex[0]].stoichiometry:
        extcomp = ma.metabolites[mid].compartment
print("   exchange metabolites live in compartment: %r" % extcomp)

results = {'as_loaded': base}
if extcomp is not None:
    ma2, mb2 = load(A), load(B)
    for m in (ma2, mb2):
        for cid, c in m.compartments.items():
            try: c.external = (cid == extcomp)
            except Exception as e: print("   cannot set external on %s: %s" % (cid, e))
    _, fixed = pools([ma2, mb2], "external=True on %r" % extcomp)
    results['external_flagged'] = fixed

    # and does it then grow, unconstrained?
    from reframed import FBA, Environment
    ma3, mb3 = load(A), load(B)
    for m in (ma3, mb3):
        for cid, c in m.compartments.items():
            try: c.external = (cid == extcomp)
            except Exception: pass
    comm3 = Community('z', [ma3, mb3])
    mdl = getattr(comm3, 'merged', None) or getattr(comm3, 'merged_model')
    obj = {r:c for r,c in mdl.get_objective().items() if c}
    print("   merged objective after flagging: %s" % (obj if obj else 'EMPTY'))
    if not obj:
        cb = [r for r in mdl.reactions if 'community_biomass' in r]
        if cb:
            mdl.set_objective({cb[0]: 1})
            print("   objective was empty; set to %s" % cb[0])
    env = Environment.complete(mdl, max_uptake=1000.0)
    env.apply(mdl, inplace=True, warning=False)
    s = FBA(mdl); g = s.fobj if s and s.fobj else 0.0
    print("   community growth, unconstrained, after flagging = %.6f" % g)
    results['growth_after_fix'] = g

with open(os.path.join(out,'merge_experiment.json'),'w') as fh: json.dump(results, fh, indent=2, default=str)
print()
b = max(base.values()) if base else 0
f = max(results.get('external_flagged',{}).values()) if results.get('external_flagged') else 0
if b == 0 and f > 0:
    print("VERDICT_CODE: H3H5   flagging the external compartment creates %d pool reactions." % f)
elif b > 0:
    print("VERDICT_CODE: H4     pools exist as loaded (%d); the earlier zero was read off the wrong attribute." % b)
else:
    print("VERDICT_CODE: DEEPER flagging the compartment did not create pools. The merge is")
    print("                     rejecting these models for a reason not yet identified.")
PYEOF

# ---------------------------------------------------------------- F. SMETANA itself
say "F. SMETANA's OWN entry point -- not my reconstruction of it (the 41e/43c mistake)"
note "whatever E concluded, the number that matters is the one SMETANA reports."
"$PY" - "$REL" "$WS2/41_env_checks/ETO_SED-medium.csv" <<'PYEOF'
import sys, os, csv, warnings, inspect
rel, medf = sys.argv[1], sys.argv[2]
from reframed import load_cbmodel
from reframed.community.model import Community
import smetana.interface as I
print("   smetana.interface exports:", [n for n in dir(I) if not n.startswith('_')])
for fn in ('sc_score','mu_score','mp_score','mip_score','mro_score'):
    if hasattr(I, fn):
        print("   %-10s %s" % (fn, inspect.signature(getattr(I, fn))))
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel, n+'.xml'), flavor='fbc2')
ma, mb = load('LIQ_das_vamb_017'), load('SED_das_semibin_003')
comm = Community('f', [ma, mb])
cpds=[]
with open(medf) as fh:
    for row in csv.DictReader(fh):
        for k in ('compounds','compound','id','cpd'):
            if k in row and row[k]: cpds.append(row[k].strip()); break
print()
for fn in ('mro_score','mip_score'):
    if not hasattr(I, fn): continue
    for label, kw in (('no medium', {}), ('ETO_SED', {'environment': None})):
        try:
            f = getattr(I, fn)
            params = inspect.signature(f).parameters
            args = {}
            if 'min_mol_weight' in params: args['min_mol_weight'] = False
            if label == 'ETO_SED' and 'environment' in params:
                from reframed import Environment
                args['environment'] = Environment.from_compounds(
                    cpds, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x))
            r = f(comm, **args)
            print("   %-10s %-10s -> %s" % (fn, label, r))
        except Exception as e:
            print("   %-10s %-10s -> EXCEPTION %s: %s" % (fn, label, type(e).__name__, e))
PYEOF

say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
