#!/usr/bin/env bash
# 43f_control_then_run.sh -- build a control that must succeed, then believe the zeros.
#
# 43e's CLI run finished cleanly and wrote pc_detailed.tsv with a header and NO rows.
# That is SMETANA reporting no cross-feeding between the methanogen and one fermenter on
# ETO_SED. It may well be true. It is not yet interpretable, because nothing in that run
# could have told us the difference between "no interaction" and "not detecting
# interaction". Rule 2.
#
# Settled by 43e, recorded here so it is not re-litigated:
#   - smetana.interface.Community IS smetana.legacy.Community, SMETANA's own class.
#     reframed's Community is a different class. Passing reframed's was my error and is
#     the sole cause of every AttributeError and every hand-built zero in 43c/43d.
#   - versions are correct: smetana 1.2.1 requires reframed>=1.6.0, 1.6.0 is installed.
#     Nothing needs pinning. My proposed pip install in 43e was aimed at a bug that
#     does not exist -- do not run it.
#   - "47 ETO_SED rows" in 43e was a bad count: ETO_SED matches ETO_SED_NOSO4 as a
#     substring. 24 + 23 = 47. The media db is correct.
#
# ORDER MATTERS HERE. Section D builds a positive control that MUST produce cross-feeding.
# Sections E and F -- the runs we actually want -- execute only if D passes. A zero from
# an instrument that has not been shown to register a one is not a measurement.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43f"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
REL="$WS2/43_smetana/models_sed"; MDB="$WS2/41_env_checks/etoliko_media.tsv"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43f  $(date -Is)"; echo "out: $OUT"
for f in "$PY" "$CLI"; do [ -x "$f" ] || { bad "missing $f"; echo "STEP_STATUS: FAILED"; exit 1; }; done
for f in "$REL" "$MDB"; do [ -e "$f" ] || { bad "missing $f"; echo "STEP_STATUS: FAILED"; exit 1; }; done

# ---------------------------------------------------------------- A. the medium
say "A. what ETO_SED actually is, by exact match (43e counted by substring)"
awk -F'\t' 'NR==1{next} $1=="ETO_SED"{n++} END{printf "  ETO_SED  exact rows: %d\n", n+0}' "$MDB"
awk -F'\t' 'NR==1{next} $1=="ETO_SED_NOSO4"{n++} END{printf "  ETO_SED_NOSO4  exact rows: %d\n", n+0}' "$MDB"
awk -F'\t' 'NR==1{next} $1=="ETO_FILL"{n++} END{printf "  ETO_FILL exact rows: %d\n", n+0}' "$MDB"
note "ETO_SED compounds:"
awk -F'\t' 'NR==1{next} $1=="ETO_SED"{printf "     %-12s %s\n",$3,$4}' "$MDB"

# ---------------------------------------------------------------- B. right class
say "B. the community, built by SMETANA's own class this time"
"$PY" - "$REL" <<'PYEOF'
import sys, os, warnings
from reframed import load_cbmodel
from smetana.interface import Community
rel = sys.argv[1]
def load(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel, n+'.xml'), flavor='fbc2')
ma, mb = load('LIQ_das_vamb_017'), load('SED_das_semibin_003')
c = Community('pc', [ma, mb])
print("   class            : %s.%s" % (type(c).__module__, type(c).__name__))
m = c.merged
pool = [r for r in m.reactions if '_pool' in r]
obj  = {r:v for r,v in m.get_objective().items() if v}
print("   merged reactions : %d" % len(m.reactions))
print("   pool reactions   : %d   <-- 43c/43d reported 0 because I built the wrong class" % len(pool))
print("   objective        : %s" % obj)
print("   example pools    : %s" % pool[:4])
PYEOF

# ---------------------------------------------------------------- C. global
say "C. global scores on the pair (MIP/MRO) -- fast, and tells us the machinery runs"
CD="$OUT/pair"; mkdir -p "$CD"
cp "$REL/LIQ_das_vamb_017.xml" "$REL/SED_das_semibin_003.xml" "$CD/"
( cd "$CD" && "$CLI" -o pair -m ETO_SED --mediadb "$MDB" --flavor fbc2 -g \
    LIQ_das_vamb_017.xml SED_das_semibin_003.xml 2>&1 ) | sed 's/^/   /' | tail -15
for f in "$CD"/pair_global.tsv; do
  [ -s "$f" ] && { echo "   --- $(basename "$f") ---"; column -t -s$'\t' "$f" 2>/dev/null | sed 's/^/     /' || cat "$f"; }
done

# ---------------------------------------------------------------- D. the control
say "D. POSITIVE CONTROL -- two auxotrophs cut from one of your own models"
note "take SED_das_semibin_003. Find two compounds it makes, that it also has an"
note "exchange reaction for, and whose biosynthesis can be knocked out so the model"
note "cannot grow on ETO_SED alone but is rescued when the compound is supplied."
note "Pair the two knockouts. Each still makes what the other lacks, so cross-feeding"
note "is not a hypothesis here -- it is the only way either can grow."
note "If SMETANA reports nothing for THIS pair, it is not detecting cross-feeding."
"$PY" - "$REL" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, warnings, json
from reframed import load_cbmodel, save_cbmodel, Environment, FBA
rel, mdb, out = sys.argv[1], sys.argv[2], sys.argv[3]
SRC = 'SED_das_semibin_003'

med = []
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row.get('medium') == 'ETO_SED':
            med.append(row['compound'].strip())
print("   ETO_SED: %d compounds" % len(med))

def load():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(rel, SRC + '.xml'), flavor='fbc2')

def grow(m, extra=()):
    mm = m.copy()
    env = Environment.from_compounds(list(med) + list(extra),
                                     fmt_func=lambda x: "R_EX_{}_e".format(x))
    env.apply(mm, inplace=True, warning=False)
    s = FBA(mm)
    return s.fobj if s and s.fobj else 0.0

base = load()
b0 = grow(base)
print("   %s grows on ETO_SED: %.4f" % (SRC, b0))
if b0 < 1e-6:
    print("CONTROL_STATUS: FAILED  source model does not grow; cannot build auxotrophs")
    sys.exit(3)

# candidate compounds: has an exchange, is NOT already in the medium
exch = base.get_exchange_reactions()
cands = []
for r in exch:
    st = base.reactions[r].stoichiometry
    if len(st) != 1: continue
    mid = list(st)[0]
    cpd = mid.replace('M_','').rsplit('_',1)[0]
    if cpd in med: continue
    name = base.metabolites[mid].name
    cands.append((cpd, mid, name))
print("   exchangeable compounds not already in the medium: %d" % len(cands))

def knockout(m, cpd):
    """remove reactions that PRODUCE cpd in the cytosol, keeping transport so it can be imported"""
    cm = 'M_%s_c' % cpd; em = 'M_%s_e' % cpd
    killed = []
    for r_id in list(m.reactions.keys()):
        st = m.reactions[r_id].stoichiometry
        if cm not in st: continue
        if em in st: continue                      # transport: keep, that is the import route
        if r_id in m.get_exchange_reactions(): continue
        killed.append(r_id)
    for r_id in killed: m.remove_reaction(r_id)
    return killed

found = []
for cpd, mid, name in cands:
    if len(found) == 2: break
    m = load()
    if ('M_%s_c' % cpd) not in m.metabolites: continue
    killed = knockout(m, cpd)
    if not killed: continue
    alone   = grow(m)                 # must be dead
    rescued = grow(m, extra=[cpd])    # must be alive
    if alone < 1e-6 and rescued > 1e-6:
        found.append(dict(cpd=cpd, name=name, killed=len(killed), rescued=rescued))
        p = os.path.join(out, 'AUXO_%s.xml' % cpd)
        m.id = 'AUXO_%s' % cpd
        save_cbmodel(m, p, flavor='fbc2')
        print("   [OK]   auxotroph for %-10s %-34s  %d rxn removed, alone=%.4f rescued=%.4f"
              % (cpd, name[:34], len(killed), alone, rescued))

print()
if len(found) < 2:
    print("   only %d usable auxotroph(s) built." % len(found))
    print("CONTROL_STATUS: FAILED  cannot form a mutually-dependent pair")
    sys.exit(3)
with open(os.path.join(out,'control_pair.json'),'w') as fh: json.dump(found, fh, indent=2)
print("   pair: AUXO_%s (%s)  +  AUXO_%s (%s)" % (found[0]['cpd'], found[0]['name'],
                                                   found[1]['cpd'], found[1]['name']))
print("   each lacks what the other still makes. Cross-feeding is obligatory.")
print("CONTROL_STATUS: BUILT")
PYEOF
DRC=$?
if [ "$DRC" -ne 0 ]; then
  bad "could not build the control pair -- sections E and F are NOT run."
  note "without a control, a zero from SMETANA cannot be distinguished from a blind"
  note "instrument. That was the whole lesson of 43c."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi

say "D2. run SMETANA on the control pair -- this MUST report cross-feeding"
KD="$OUT/control"; mkdir -p "$KD"; cp "$OUT"/AUXO_*.xml "$KD/" 2>/dev/null
( cd "$KD" && "$CLI" -o ctrl -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d AUXO_*.xml 2>&1 ) | sed 's/^/   /' | tail -12
CTRL="$KD/ctrl_detailed.tsv"
NROW=0; [ -s "$CTRL" ] && NROW=$(($(wc -l < "$CTRL") - 1))
echo "   ctrl_detailed.tsv rows: $NROW"
[ "$NROW" -gt 0 ] && { echo "   --- first rows ---"; head -6 "$CTRL" | column -t -s$'\t' | sed 's/^/     /'; }
if [ "$NROW" -eq 0 ]; then
  bad "the control pair reports NO cross-feeding."
  note "two organisms that cannot grow without each other, and SMETANA sees nothing."
  note "the detector is blind. Every zero so far -- including 43e's -- means nothing."
  note "sections E and F are NOT run."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "control passed: $NROW interaction row(s). The detector registers cross-feeding."
note "so from here, a zero is a result about the organisms, not about the method."

# ---------------------------------------------------------------- E. the ten
say "E. the in-situ ten on ETO_SED -- global (MIP/MRO)"
TD="$OUT/ten"; mkdir -p "$TD"; cp "$REL"/*.xml "$TD/"
note "$(ls "$TD"/*.xml | wc -l) models. Global mode; this is the faster of the two."
( cd "$TD" && "$CLI" -o ten -m ETO_SED --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | sed 's/^/   /' | tail -15
[ -s "$TD/ten_global.tsv" ] && { echo "   --- ten_global.tsv ---"; column -t -s$'\t' "$TD/ten_global.tsv" | sed 's/^/     /'; }

say "F. the in-situ ten on ETO_SED -- detailed (who feeds whom)"
note "slower. If it runs long, that is expected; it solves many LPs per compound."
( cd "$TD" && "$CLI" -o ten -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d "*.xml" 2>&1 ) | sed 's/^/   /' | tail -15
DET="$TD/ten_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET") - 1)); echo "   ten_detailed.tsv rows: $N"
  if [ "$N" -gt 0 ]; then
    echo "   --- strongest 15 predicted interactions ---"
    { head -1 "$DET"; tail -n +2 "$DET" | sort -t$'\t' -k9,9gr | head -15; } | column -t -s$'\t' | sed 's/^/     /'
    echo
    echo "   --- donors, by number of compounds donated ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(d in c) printf "     %-26s %d\n", d, c[d]}' | sort -k2,2nr
    echo "   --- receivers ---"
    tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(d in c) printf "     %-26s %d\n", d, c[d]}' | sort -k2,2nr
  else
    note "zero rows -- and now that means something, because the control passed."
  fi
fi

say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
