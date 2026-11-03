#!/usr/bin/env bash
# 43e_use_smetana_properly.sh -- stop building the community myself.
#
# 43d found the cause and it is not biological:
#
#     mro_score -> AttributeError: 'Community' object has no attribute 'merged'
#     mip_score -> AttributeError: 'Community' object has no attribute 'copy'
#
# The models are fine. Section C showed both declare "Extracellular space" external,
# 163 metabolites in it, exchanges pointing at it. merge_models() runs and produces a
# 2073-reaction model with a community_growth objective. My "0 pool reactions" count was
# wrong too -- I counted names ending in _pool, a string taken from SMETANA's fmt_func,
# not from what merge_models actually names things.
#
# The real fault: I passed reframed.community.model.Community into SMETANA. SMETANA
# expects the class IT uses, which has .merged and .copy(). Either that is a different
# class, or smetana 1.2.1 is pinned to a reframed older than the 1.6.0 installed here.
#
# This script does not build a community by hand anywhere. It:
#   A  compares the Community classes by identity and by API
#   B  reads what smetana 1.2.1 DECLARES it needs, against what is installed
#   C  if a usable class exists, uses SMETANA's own loader and scoring
#   D  runs the smetana COMMAND LINE -- the tool as its author intends it
#
# It changes nothing and installs nothing. If a version pin is required, section B says
# so and section E prints the exact command, for you to run deliberately.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/diag_43e"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
echo "43e  $(date -Is)"; echo "out: $OUT"
PY="$WS2/envs/smetana/bin/python"
CLI="$WS2/envs/smetana/bin/smetana"
REL="$WS2/43_smetana/models_sed"
export PYTHONWARNINGS="ignore::UserWarning"
[ -x "$PY" ] || { bad "no python at $PY"; echo "STEP_STATUS: FAILED"; exit 1; }
[ -d "$REL" ] || { bad "no models at $REL"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- A. which class
say "A. which Community class does SMETANA actually use?"
"$PY" - <<'PYEOF'
import smetana.interface as I
import reframed
import reframed.community.model as RM
cands = {
  'smetana.interface.Community'      : getattr(I, 'Community', None),
  'reframed.Community (top level)'   : getattr(reframed, 'Community', None),
  'reframed.community.model.Community': getattr(RM, 'Community', None),
}
for name, cls in cands.items():
    if cls is None:
        print("   %-38s ABSENT" % name); continue
    print("   %-38s %s.%s  id=%s" % (name, cls.__module__, cls.__name__, hex(id(cls))))
    print("   %-38s has .merged=%-5s  has .copy=%-5s  has .merged_model=%s" % (
        '', hasattr(cls,'merged'), hasattr(cls,'copy'), hasattr(cls,'merged_model')))
ids = {n: hex(id(c)) for n,c in cands.items() if c is not None}
print()
print("   all the same class? %s" % (len(set(ids.values())) == 1))
usable = [n for n,c in cands.items() if c is not None and hasattr(c,'merged') and hasattr(c,'copy')]
print("   classes providing BOTH .merged and .copy: %s" % (usable if usable else "NONE"))
PYEOF

# ---------------------------------------------------------------- B. versions
say "B. what smetana DECLARES it needs, vs what is installed"
"$PY" - <<'PYEOF'
from importlib.metadata import distribution, version, PackageNotFoundError
for pkg in ('smetana','reframed'):
    try: print("   installed  %-10s %s" % (pkg, version(pkg)))
    except PackageNotFoundError: print("   installed  %-10s NOT FOUND" % pkg)
try:
    reqs = distribution('smetana').requires or []
    print("   smetana requires:")
    for r in reqs: print("      %s" % r)
    pin = [r for r in reqs if 'reframed' in r.lower()]
    print()
    if pin:
        print("   >>> the reframed requirement is: %s" % '; '.join(pin))
    else:
        print("   >>> smetana declares NO reframed pin -- so pip installed the newest,")
        print("       and the newest changed the Community API. That is the whole bug.")
except Exception as e:
    print("   could not read smetana metadata: %s" % e)
PYEOF

# ---------------------------------------------------------------- C. own loader
say "C. SMETANA's OWN loader and scoring -- no community built by me"
note "load_communities() is how the CLI turns model files into communities."
"$PY" - "$REL" "$WS2/41_env_checks/ETO_SED-medium.csv" <<'PYEOF'
import sys, os, csv, inspect, traceback
import smetana.interface as I
rel, medf = sys.argv[1], sys.argv[2]
A, B = 'LIQ_das_vamb_017', 'SED_das_semibin_003'
paths = [os.path.join(rel, n + '.xml') for n in (A, B)]
for p in paths:
    if not os.path.exists(p): print("   missing %s" % p); sys.exit(1)

print("   load_communities signature: %s" % inspect.signature(I.load_communities))
comm = None
try:
    comms = I.load_communities(paths, None, flavor='fbc2')
    print("   load_communities returned: %r" % type(comms))
    if isinstance(comms, dict):
        for k, v in comms.items():
            print("     community %r -> %s" % (k, type(v).__name__))
            comm = v
    else:
        comm = comms
except Exception as e:
    print("   load_communities raised %s: %s" % (type(e).__name__, e))
    traceback.print_exc()

if comm is not None:
    print("   community class: %s.%s" % (type(comm).__module__, type(comm).__name__))
    print("   has .merged=%s  has .copy=%s" % (hasattr(comm,'merged'), hasattr(comm,'copy')))
    cpds = []
    with open(medf) as fh:
        for row in csv.DictReader(fh):
            for k in ('compounds','compound','id','cpd'):
                if k in row and row[k]: cpds.append(row[k].strip()); break
    from reframed import Environment
    env = Environment.from_compounds(cpds, fmt_func=lambda x: "R_EX_M_{}_e_pool".format(x))
    for fn in ('mro_score','mip_score'):
        for label, kw in (('no medium', {}), ('ETO_SED', {'environment': env})):
            try:
                r = getattr(I, fn)(comm, verbose=False, **kw)
                print("   %-10s %-10s -> %s" % (fn, label, r))
            except Exception as e:
                print("   %-10s %-10s -> %s: %s" % (fn, label, type(e).__name__, e))
else:
    print("   no community object obtained -- section D is the remaining route.")
PYEOF

# ---------------------------------------------------------------- D. the CLI
say "D. the smetana COMMAND LINE -- the tool as its author intends it"
if [ ! -x "$CLI" ]; then
  no "no smetana executable at $CLI"
else
  ok "found $CLI"
  "$CLI" --help 2>&1 | sed 's/^/   /' | head -40
  MDB="$WS2/41_env_checks/etoliko_media.tsv"
  if [ ! -s "$MDB" ]; then
    no "no media db at $MDB -- cannot run the CLI on ETO_SED"
  else
    note "media db header and the ETO_SED rows it holds:"
    head -1 "$MDB" | sed 's/^/     /'
    awk 'NR==1{next} /ETO_SED/{n++} END{printf "     ETO_SED rows: %d\n", n+0}' "$MDB"
    CD="$OUT/cli"; mkdir -p "$CD"
    cp "$REL/LIQ_das_vamb_017.xml" "$REL/SED_das_semibin_003.xml" "$CD/" 2>/dev/null
    note "running: smetana -o pc -m ETO_SED --mediadb <db> --flavor fbc2 -d <2 models>"
    ( cd "$CD" && "$CLI" -o pc -m ETO_SED --mediadb "$MDB" --flavor fbc2 -d \
        LIQ_das_vamb_017.xml SED_das_semibin_003.xml 2>&1 ) | sed 's/^/   /' | tail -30
    echo
    note "what the CLI wrote:"
    ls -l "$CD" | sed 's/^/     /'
    for f in "$CD"/*detailed*.tsv "$CD"/*.tsv; do
      [ -s "$f" ] || continue
      echo "     --- $(basename "$f") ($(($(wc -l < "$f")-1)) rows) ---"
      head -8 "$f" | sed 's/^/       /'
      break
    done
  fi
fi

# ---------------------------------------------------------------- E. verdict
say "E. verdict"
"$PY" - <<'PYEOF'
import smetana.interface as I
from importlib.metadata import distribution, version
cls = getattr(I, 'Community', None)
usable = cls is not None and hasattr(cls,'merged') and hasattr(cls,'copy')
reqs = []
try: reqs = [r for r in (distribution('smetana').requires or []) if 'reframed' in r.lower()]
except Exception: pass
print("   smetana %s + reframed %s" % (version('smetana'), version('reframed')))
if usable:
    print("   VERDICT_CODE: CLASS   SMETANA's own Community has .merged and .copy.")
    print("     The zero came from me passing the wrong class. Sections C/D show the")
    print("     numbers the tool gives when used as intended.")
else:
    print("   VERDICT_CODE: VERSION SMETANA's Community lacks .merged/.copy, so no way of")
    print("     calling it from this environment will work. reframed 1.6.0 renamed the")
    print("     attribute to merged_model. The environment must be pinned.")
    print()
    print("     Declared requirement: %s" % (reqs or "none -- pip took the newest reframed"))
    print("     Proposed fix, to run DELIBERATELY and not from inside a diagnostic:")
    print()
    print("       $WS2/envs/smetana/bin/pip install 'reframed==1.4.0'")
    print()
    print("     That venv is ours (workspace_2/envs), not a shared conda env, so this is")
    print("     safe under the standing constraint. Re-run 41m's positive control after.")
PYEOF

say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
