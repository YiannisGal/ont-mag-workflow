#!/usr/bin/env bash
# 43h_refill_reversal.sh -- the gapfilling medium decided the answer. Re-test it.
#
# RUN THIS IN tmux. The detailed section took ~3.5 h last time and your SSH dropped twice:
#     tmux new -s 43h
#     bash scripts/43h_refill_reversal.sh 2>&1 | tee logs/43h.log
#     (detach with Ctrl-b then d; return with  tmux attach -t 43h)
#
# WHAT 43g ESTABLISHED
#   - the detector works. Two engineered auxotrophs, 59 interaction rows, lysine
#     transferred at smetana=1.0. A zero from SMETANA is now interpretable.
#   - the in-situ ten on ETO_SED: 0 detailed rows, MRO 0.79, MIP n/a.
#
# WHY I DO NOT BELIEVE THAT ZERO
#   43b refilled every model ON ETO_SED. Gapfilling adds reactions until the organism
#   grows on the medium it is given. Filling on a 24-compound minimal medium therefore
#   MAKES every model a prototroph on those 24 compounds, and a prototroph needs nothing
#   from anybody. "No cross-feeding on ETO_SED" is a restatement of the gapfilling
#   instruction, not a result about the lagoon. Rule 63.
#
#   The models filled on ETO_FILL (51 compounds, rich) keep their biosynthetic gaps.
#   Those gaps ARE the dependencies SMETANA is built to find. They are the correct input
#   for community analysis, and I reversed away from them in 41l/43b on the strength of
#   "the community does not grow on ETO_SED" -- a number produced with reframed's
#   Community instead of SMETANA's. That number is void. Rule 69: when a correction
#   reverses a decision, re-test the failure the original decision was avoiding.
#
# This script does not assume which directory holds which. Section A decides by
# measurement. Sections C-D re-take the void measurement with the right class.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43h"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
SKIP_DETAILED="${SKIP_DETAILED:-0}"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43h  $(date -Is)"; echo "out: $OUT"; echo "SKIP_DETAILED=$SKIP_DETAILED"

# ---------------------------------------------------------------- A. provenance
say "A. which directory is which -- by measurement, not by name"
note "a set filled on ETO_SED grows on ETO_SED by construction. A set filled on"
note "ETO_FILL should have members that CANNOT grow on ETO_SED. That difference"
note "identifies them without trusting a directory name."
"$PY" - "$WS2" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, glob, json, warnings
from reframed import load_cbmodel, Environment, FBA
ws2, mdb, out = sys.argv[1], sys.argv[2], sys.argv[3]
media = {}
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        media.setdefault(row['medium'], []).append(row['compound'].strip())
print("   media in db: %s" % {k: len(v) for k, v in media.items()})

def growth(path, cpds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = load_cbmodel(path, flavor='fbc2')
    Environment.from_compounds(cpds, fmt_func=lambda x: "R_EX_{}_e".format(x)
                               ).apply(m, inplace=True, warning=False)
    s = FBA(m)
    return (s.fobj if s and s.fobj else 0.0), len(m.reactions)

summary = {}
for d in sorted(glob.glob(os.path.join(ws2, '43_smetana', 'models*'))):
    xs = sorted(glob.glob(os.path.join(d, '*.xml')))
    if len(xs) < 5: continue
    dead_sed = 0; tot_rxn = 0; rows = []
    for p in xs:
        g, n = growth(p, media['ETO_SED'])
        tot_rxn += n
        if g < 1e-6: dead_sed += 1
        rows.append((os.path.basename(p).replace('.xml',''), g, n))
    summary[d] = dict(n=len(xs), dead_on_sed=dead_sed, mean_rxn=tot_rxn/len(xs), rows=rows)
    print("\n   %s" % d)
    print("     %d models, mean %.0f reactions, %d cannot grow on ETO_SED"
          % (len(xs), tot_rxn/len(xs), dead_sed))
    for n_, g, r in rows:
        print("       %-24s %6d rxn   ETO_SED growth %.4f%s" % (n_, r, g, "   <-- dead" if g < 1e-6 else ""))

json.dump({k: {kk: vv for kk, vv in v.items() if kk != 'rows'} for k, v in summary.items()},
          open(os.path.join(out, 'model_sets.json'), 'w'), indent=2)
prot = [d for d, v in summary.items() if v['dead_on_sed'] == 0]
gapped = [d for d, v in summary.items() if v['dead_on_sed'] > 0]
print()
print("   PROTOTROPHIC sets (all grow on ETO_SED, filled on it): %s" % [os.path.basename(d) for d in prot])
print("   GAP-RETAINING sets (some cannot):                      %s" % [os.path.basename(d) for d in gapped])
if not gapped:
    print("   >>> no gap-retaining set exists. Every model on disk was filled to prototrophy.")
    print("   >>> the ETO_FILL models must be rebuilt before community analysis can mean anything.")
open(os.path.join(out,'chosen.txt'),'w').write((gapped[0] if gapped else '') + "\n")
PYEOF
CHOSEN=$(cat "$OUT/chosen.txt" 2>/dev/null | head -1)
if [ -z "$CHOSEN" ]; then
  bad "no gap-retaining model set on disk."
  note "every set was filled to prototrophy on ETO_SED, so no dependency can exist."
  note "the fix is to re-run gapseq fill on ETO_FILL and keep BOTH sets. That is a"
  note "rebuild, not a diagnosis, so it gets its own script rather than happening here."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "gap-retaining set: $CHOSEN"

# ---------------------------------------------------------------- B. void measurement
say "B. the measurement that caused the reversal, re-taken with SMETANA's class"
note "43/43b concluded the ETO_FILL community could not grow on ETO_SED. That used"
note "reframed's Community. This uses smetana.legacy.Community."
"$PY" - "$CHOSEN" "$MDB" <<'PYEOF'
import sys, os, csv, glob, warnings
from reframed import load_cbmodel, Environment, FBA
from smetana.interface import Community
d, mdb = sys.argv[1], sys.argv[2]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
models=[]
for p in sorted(glob.glob(os.path.join(d,'*.xml'))):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        models.append(load_cbmodel(p, flavor='fbc2'))
c = Community('eto', models)
mm = c.merged.copy()
Environment.from_compounds(med, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x)
                           ).apply(mm, inplace=True, warning=False)
s = FBA(mm); g = s.fobj if s and s.fobj else 0.0
print("   members                         : %d" % len(models))
print("   merged reactions                : %d" % len(mm.reactions))
print("   pool reactions                  : %d" % len([r for r in mm.reactions if '_pool' in r]))
print("   COMMUNITY growth on ETO_SED     : %.6f" % g)
if g > 1e-6:
    print("   >>> the community DOES grow. The failure that justified the 43b refill was")
    print("       an artefact of the wrong class. The refill was not needed.")
else:
    print("   >>> still zero. The refill had a real cause; it must be revisited another way.")
PYEOF

# ---------------------------------------------------------------- C. mip
say "C. why MIP is n/a -- can every member grow at once?"
note "MIP needs all members to reach min_growth simultaneously. Community biomass"
note "0.45 does not imply that: the community objective can be carried by a subset."
"$PY" - "$CHOSEN" "$MDB" <<'PYEOF'
import sys, os, csv, glob, warnings
from reframed import load_cbmodel, Environment, FBA
from smetana.interface import Community
d, mdb = sys.argv[1], sys.argv[2]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
models=[]
for p in sorted(glob.glob(os.path.join(d,'*.xml'))):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        models.append(load_cbmodel(p, flavor='fbc2'))
c = Community('eto', models)
mm = c.merged.copy()
Environment.from_compounds(med, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x)
                           ).apply(mm, inplace=True, warning=False)
# find each member's growth reaction inside the merged model
growth_rxns = {}
for m in models:
    hits = [r for r in mm.reactions if m.id in r and 'bio' in r.lower()]
    if hits: growth_rxns[m.id] = hits[0]
print("   located %d member growth reactions in the merged model" % len(growth_rxns))
print()
print("   each member's maximum growth INSIDE the community, others unconstrained:")
for mid, r in growth_rxns.items():
    m2 = mm.copy(); m2.set_objective({r: 1})
    s = FBA(m2); g = s.fobj if s and s.fobj else 0.0
    print("     %-26s %.4f%s" % (mid, g, "   <-- below min_growth 0.1" if g < 0.1 else ""))
print()
print("   now all members forced to grow simultaneously at 0.1:")
m3 = mm.copy()
for mid, r in growth_rxns.items(): m3.reactions[r].lb = 0.1
s = FBA(m3)
print("     feasible? %s   objective=%s" % (s is not None and s.fobj is not None,
                                            None if not s else s.fobj))
if s is None or s.fobj is None:
    print("     >>> infeasible. THAT is why mip_score returns None. Not a bug.")
PYEOF

# ---------------------------------------------------------------- D. runs
say "D. global scores -- gap-retaining set on ETO_SED"
TD="$OUT/gapped"; mkdir -p "$TD"; cp "$CHOSEN"/*.xml "$TD/"
( cd "$TD" && "$CLI" -o gap -m ETO_SED --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | sed 's/^/   /' | tail -10
[ -s "$TD/gap_global.tsv" ] && { echo "   --- gap_global.tsv ---"; column -t -s$'\t' "$TD/gap_global.tsv" | sed 's/^/     /'; }

if [ "$SKIP_DETAILED" = "1" ]; then
  note "SKIP_DETAILED=1 -- stopping before the slow run."
  say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0
fi

say "E. detailed -- who feeds whom. THIS IS THE ~3 HOUR STEP"
note "started $(date -Is). If your SSH drops, tmux keeps it alive."
if [ -s "$TD/gap_detailed.tsv" ]; then
  note "gap_detailed.tsv already present -- not re-running. Delete it to force."
else
  ( cd "$TD" && "$CLI" -o gap -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | sed 's/^/   /' | tail -10
fi
DET="$TD/gap_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET") - 1)); echo "   gap_detailed.tsv rows: $N"
  if [ "$N" -gt 0 ]; then
    echo "   --- strongest 20 predicted interactions ---"
    { head -1 "$DET"; tail -n +2 "$DET" | sort -t$'\t' -k9,9gr | head -20; } | column -t -s$'\t' | sed 's/^/     /'
    echo "   --- donors (compounds donated with smetana>0) ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    echo "   --- receivers ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
    echo "   --- most-traded compounds ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$5]+=$9} END{for(x in c) printf "     %-20s %.2f\n",x,c[x]}' | sort -k2,2gr | head -20
    echo
    note "compare: the ETO_SED-filled set gave 0 rows. If this set gives many, the"
    note "difference is the gapfilling medium, and that belongs in the methods."
  else
    no "zero rows here too -- then the gapfilling medium is NOT the explanation,"
    no "and the ten really are self-sufficient on ETO_SED."
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
