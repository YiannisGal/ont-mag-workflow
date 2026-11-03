#!/usr/bin/env bash
# 43c_biomass_diagnosis.sh -- why does every member grow and the community not?
#
# 43b fixed the medium: 10/10 models now grow on ETO_SED (0.66-1.17). The positive
# control -- methanogen + fermenter, which MUST exchange something -- is still zero,
# and reframed emitted:
#
#     UserWarning: Suspicious biomass identifier: R_EX_cpd11416_c
#
# Two candidate causes, and this script is built to tell them apart rather than to
# confirm the one I like:
#
#   H1 BIOMASS.  reframed has taken the biomass DRAIN (EX_cpd11416_c, the sink that
#                empties ModelSEED's biomass pseudo-metabolite) as the biomass
#                reaction instead of bio1. It is a boundary reaction, so the community
#                merge will also treat it as an exchange and rewire it. Individual FBA
#                is unharmed because that uses the SBML objective.
#
#   H2 MEDIUM.   SMETANA constrains the MERGED model's pool exchanges,
#                R_EX_M_<met>_e_pool (41d, smetana/interface.py:112). If our compound
#                IDs do not land on that naming, the community is fed nothing.
#
# The discriminator is section E: run the community with NO medium at all.
#   still zero  -> structural, H1 lives, H2 dies.
#   grows       -> the medium is not reaching the community, H2 lives, H1 dies.
#   BOTH grow   -> neither; the positive control itself is wrong. Also a real outcome.
#
# Nothing is fixed here. This script only measures. 43d applies whatever E selects.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/diag_43c"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
echo "43c  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- A. python
say "A. the interpreter, by positive identification"
PY="$WS2/envs/smetana/bin/python"
if [ ! -x "$PY" ]; then bad "no smetana python at $PY"; echo "STEP_STATUS: FAILED"; exit 1; fi
"$PY" - <<'EOF' || { echo "STEP_STATUS: FAILED"; exit 1; }
import sys
print("  python  ", sys.executable)
import reframed, smetana
print("  reframed", reframed.__version__, reframed.__file__)
print("  smetana ", getattr(smetana,'__version__','?'), smetana.__file__)
EOF
ok "interpreter and both libraries import"

# ---------------------------------------------------------------- B. sources
say "B. what the libraries actually do -- read, do not assume (41e)"
note "how reframed decides what the biomass reaction is:"
"$PY" - <<'EOF'
import inspect, re, reframed.core.cbmodel as M
src = inspect.getsource(M)
# the warning we saw, and its neighbourhood
for m in re.finditer(r'Suspicious biomass', src):
    lines = src[:m.start()].count('\n')
    body = src.split('\n')
    lo, hi = max(0, lines-22), min(len(body), lines+4)
    for i in range(lo, hi):
        print("   %4d | %s" % (i+1, body[i]))
EOF
note "how SMETANA constrains the community medium:"
"$PY" - <<'EOF'
import inspect, re, smetana.interface as I
body = inspect.getsource(I).split('\n')
for i, l in enumerate(body):
    if 'fmt_func' in l or '_pool' in l or 'Environment' in l:
        print("   %4d | %s" % (i+1, l.rstrip()))
EOF
note "what Community does with member biomass reactions:"
"$PY" - <<'EOF'
import inspect, re
try:
    import reframed.community.model as C
except Exception:
    import reframed.community as C
body = inspect.getsource(C).split('\n')
hits = [i for i,l in enumerate(body) if 'biomass' in l.lower()]
for i in hits:
    print("   %4d | %s" % (i+1, body[i].rstrip()))
print("   (%d line(s) mentioning biomass in %s)" % (len(hits), C.__file__))
EOF

# ---------------------------------------------------------------- C. models
say "C. which models -- located, then verified, not assumed"
mapfile -t DIRS < <(find "$WS2/43_smetana" -type d -name '*' -print 2>/dev/null | while read -r d; do
  n=$(find "$d" -maxdepth 1 -name '*.xml' 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] && printf '%s\t%s\t%s\n' "$n" "$(date -r "$d" +%Y-%m-%dT%H:%M)" "$d"
done | sort -k1,1nr -k2,2r)
if [ "${#DIRS[@]}" -eq 0 ]; then bad "no directory under 43_smetana holds .xml models"; echo "STEP_STATUS: FAILED"; exit 1; fi
note "directories holding models (count, mtime, path):"
for d in "${DIRS[@]}"; do printf '     %s\n' "$d"; done
MODELS=$(printf '%s\n' "${DIRS[0]}" | cut -f3)
NMOD=$(printf '%s\n' "${DIRS[0]}" | cut -f1)
note "chosen: the newest directory with the most models -> $MODELS ($NMOD)"
note "that is a heuristic. The gate below is what makes it safe: these models must"
note "each grow on ETO_SED, which only the 43b-refilled set does."

MED="$WS2/41_env_checks/ETO_SED-medium.csv"
[ -s "$MED" ] || { bad "no ETO_SED medium at $MED"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "medium: $MED ($(($(wc -l < "$MED")-1)) compounds)"

"$PY" - "$MODELS" "$MED" "$OUT" <<'EOF'
import sys, os, glob, csv, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, medf, out = sys.argv[1], sys.argv[2], sys.argv[3]

cpds = []
with open(medf) as fh:
    for row in csv.DictReader(fh):
        for k in ('compounds','compound','id','cpd'):
            if k in row and row[k]:
                cpds.append(row[k].strip()); break
cpds = [c for c in cpds if c]
print("  medium compounds parsed: %d  e.g. %s" % (len(cpds), ', '.join(cpds[:5])))

rows = []
for f in sorted(glob.glob(os.path.join(mdir, '*.xml'))):
    name = os.path.basename(f).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = load_cbmodel(f, flavor='fbc2')
    # what reframed thinks the biomass reaction is
    bio = m.biomass_reaction
    # what the SBML objective actually is
    obj = [r for r, c in m.get_objective().items() if c]
    # is bio1 present under any spelling?
    cands = [r for r in m.reactions if 'bio' in r.lower() and 'EX_' not in r]
    exch = set(m.get_exchange_reactions())
    rows.append(dict(name=name, bio=bio, obj=';'.join(obj), nrxn=len(m.reactions),
                     nexch=len(exch), bio_is_exch=bio in exch,
                     obj_is_exch=all(o in exch for o in obj) if obj else None,
                     cands=';'.join(cands[:4])))

w = max(len(r['name']) for r in rows)
print()
print("  %-*s  %-22s %-22s %5s %5s  %s" % (w,'model','biomass_reaction','SBML objective','rxn','exch','biomass is an exchange?'))
for r in rows:
    print("  %-*s  %-22s %-22s %5d %5d  %s" % (w, r['name'], r['bio'], r['obj'], r['nrxn'], r['nexch'],
          'YES  <-- ' if r['bio_is_exch'] else 'no'))
print()
nbad = sum(1 for r in rows if r['bio_is_exch'])
print("  %d/%d models have a biomass_reaction that is also an exchange reaction." % (nbad, len(rows)))
mism = [r for r in rows if r['obj'] and r['bio'] not in r['obj'].split(';')]
print("  %d/%d models where biomass_reaction != SBML objective." % (len(mism), len(rows)))
if mism:
    print("  candidate true biomass reactions present in the model, e.g.:")
    for r in mism[:3]:
        print("    %-24s objective=%-16s candidates=%s" % (r['name'], r['obj'], r['cands'] or '(none)'))
with open(os.path.join(out,'biomass_inventory.tsv'),'w') as fh:
    fh.write('model\tbiomass_reaction\tsbml_objective\tn_rxn\tn_exch\tbiomass_is_exchange\tcandidates\n')
    for r in rows:
        fh.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\n' % (r['name'],r['bio'],r['obj'],r['nrxn'],r['nexch'],r['bio_is_exch'],r['cands']))
EOF

# ---------------------------------------------------------------- D. gate
say "D. GATE -- these must be the refilled models"
"$PY" - "$MODELS" "$MED" <<'EOF'
import sys, os, glob, csv, warnings
from reframed import load_cbmodel, Environment, FBA
mdir, medf = sys.argv[1], sys.argv[2]
cpds=[]
with open(medf) as fh:
    for row in csv.DictReader(fh):
        for k in ('compounds','compound','id','cpd'):
            if k in row and row[k]: cpds.append(row[k].strip()); break
zero=0; tot=0
for f in sorted(glob.glob(os.path.join(mdir,'*.xml'))):
    tot+=1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(f, flavor='fbc2')
    env=Environment.from_compounds(cpds, fmt_func=lambda x: "R_EX_{}_e".format(x))
    env.apply(m, inplace=True, warning=False)
    s=FBA(m); g=s.fobj if s and s.fobj else 0.0
    if g<1e-6: zero+=1
    print("  %-24s growth on ETO_SED = %.4f" % (os.path.basename(f).replace('.xml',''), g))
print()
if zero:
    print("GATE_STATUS: FAILED  %d/%d models do not grow -- wrong model set" % (zero,tot))
    sys.exit(3)
print("GATE_STATUS: PASSED  %d/%d grow" % (tot-zero,tot))
EOF
GRC=$?
if [ "$GRC" -ne 0 ]; then
  bad "gate failed -- the chosen directory is not the 43b-refilled set."
  note "sections E and F are NOT run. Point WS2/43_smetana at the refilled models."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "gate passed -- proceeding with sections E and F INSIDE the gate (Rule 60)"

# ---------------------------------------------------------------- E. discriminator
say "E. THE DISCRIMINATOR -- community with medium vs community with no medium"
note "positive control: LIQ_das_vamb_017 (Methanofollis, archaeon) + SED_das_semibin_003"
note "if both conditions are zero  -> structural (H1 biomass)"
note "if no-medium grows, ETO_SED does not -> medium not reaching the pool (H2)"
note "if both grow -> neither hypothesis; the earlier control call was wrong"
"$PY" - "$MODELS" "$MED" "$OUT" <<'EOF'
import sys, os, csv, warnings, itertools
from reframed import load_cbmodel, Environment, FBA
from reframed.community.model import Community
mdir, medf, out = sys.argv[1], sys.argv[2], sys.argv[3]
A, B = 'LIQ_das_vamb_017', 'SED_das_semibin_003'

def find(n):
    p = os.path.join(mdir, n + '.xml')
    return p if os.path.exists(p) else None
pa, pb = find(A), find(B)
if not pa or not pb:
    print("  [FAIL] positive-control members not found in %s" % mdir); sys.exit(4)

cpds=[]
with open(medf) as fh:
    for row in csv.DictReader(fh):
        for k in ('compounds','compound','id','cpd'):
            if k in row and row[k]: cpds.append(row[k].strip()); break

def load(p):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(p, flavor='fbc2')

results = {}
for label, force_bio in (('as reframed loads it', False), ('biomass forced to the SBML objective', True)):
    ma, mb = load(pa), load(pb)
    detail = []
    for m in (ma, mb):
        before = m.biomass_reaction
        if force_bio:
            obj = [r for r, c in m.get_objective().items() if c]
            cand = [o for o in obj if o not in set(m.get_exchange_reactions())]
            if not cand:
                # objective is itself an exchange; fall back to a non-exchange 'bio' reaction
                cand = [r for r in m.reactions if 'bio' in r.lower() and r not in set(m.get_exchange_reactions())]
            if cand:
                m.biomass_reaction = cand[0]
        detail.append((m.id, before, m.biomass_reaction))
    for mid, b0, b1 in detail:
        print("    %-28s biomass %-22s -> %s" % (mid, b0, b1))
    comm = Community('pc', [ma, mb])
    cm = comm.merged_model
    npool = len([r for r in cm.reactions if r.endswith('_pool')])
    print("    merged: %d reactions, %d pool reactions" % (len(cm.reactions), npool))

    for cond, use_med in (('ETO_SED', True), ('NO MEDIUM (unconstrained uptake)', False)):
        mm = cm.copy()
        if use_med:
            env = Environment.from_compounds(cpds, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x))
            env.apply(mm, inplace=True, warning=False)
        s = FBA(mm)
        g = s.fobj if s and s.fobj else 0.0
        results[(label, cond)] = g
        print("    [%s] %-34s community growth = %.6f" % ('OK ' if g > 1e-6 else '-- ', cond, g))
    print()

print("  summary")
print("  %-40s %14s %14s" % ('', 'ETO_SED', 'no medium'))
for label in ('as reframed loads it', 'biomass forced to the SBML objective'):
    print("  %-40s %14.6f %14.6f" % (label,
          results[(label,'ETO_SED')], results[(label,'NO MEDIUM (unconstrained uptake)')]))
print()
asis_med   = results[('as reframed loads it','ETO_SED')] > 1e-6
asis_nomed = results[('as reframed loads it','NO MEDIUM (unconstrained uptake)')] > 1e-6
fix_med    = results[('biomass forced to the SBML objective','ETO_SED')] > 1e-6

if not asis_med and not asis_nomed and fix_med:
    print("VERDICT: H1 -- the biomass reaction. Forcing it off the exchange restores growth.")
    print("VERDICT_CODE: H1")
elif not asis_med and asis_nomed:
    print("VERDICT: H2 -- the medium is not reaching the community pool. Growth without a")
    print("         medium, none with one, and the biomass fix does not change that.")
    print("VERDICT_CODE: H2")
elif not asis_med and not asis_nomed and not fix_med:
    print("VERDICT: NEITHER -- zero in every condition including unconstrained uptake.")
    print("         An unconstrained community that cannot grow is a structural defect in")
    print("         the merge itself. Do not proceed to SMETANA; 43d must dissect Community.")
    print("VERDICT_CODE: NONE")
else:
    print("VERDICT: the community grows as loaded. The earlier zero came from something")
    print("         else in 43b's control -- most likely how the score was called, not the")
    print("         model. 43d re-runs the control through SMETANA's own entry points.")
    print("VERDICT_CODE: CONTROL")
with open(os.path.join(out,'discriminator.tsv'),'w') as fh:
    fh.write('condition\tmedium\tgrowth\n')
    for (l,c),g in results.items(): fh.write('%s\t%s\t%.6f\n' % (l,c,g))
EOF
ERC=$?
[ "$ERC" -ne 0 ] && { bad "section E did not complete (rc=$ERC)"; echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1; }

# ---------------------------------------------------------------- F. sanity
say "F. a control on the control -- can this test report a failure?"
note "an all-zero table is only meaningful if a working pair would report non-zero."
note "so: the same two models merged, unconstrained, with EVERY exchange opened wide."
"$PY" - "$MODELS" <<'EOF'
import sys, os, warnings
from reframed import load_cbmodel, FBA
from reframed.community.model import Community
mdir = sys.argv[1]
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(mdir, n + '.xml'), flavor='fbc2')
ma, mb = load('LIQ_das_vamb_017'), load('SED_das_semibin_003')
comm = Community('sanity', [ma, mb])
cm = comm.merged_model.copy()
n = 0
for r in cm.reactions:
    if r.endswith('_pool') or '_EX_' in r or r.startswith('R_EX_'):
        cm.reactions[r].lb = -1000.0; cm.reactions[r].ub = 1000.0; n += 1
s = FBA(cm); g = s.fobj if s and s.fobj else 0.0
print("  opened %d boundary reactions; community growth = %.6f" % (n, g))
if g > 1e-6:
    print("  [OK]   the test CAN report growth -- an all-zero table above is a real zero.")
else:
    print("  [FAIL] even wide open the merged community cannot grow.")
    print("         The zero is in the merge, not in the medium and not in the biology.")
EOF

say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
