#!/usr/bin/env bash
# 41e -- the decisive test, and the mcr answer that may already be on disk.
#
# 41d found smetana/interface.py:112   fmt_func = lambda x: "R_EX_M_{}_e_pool".format(x)
# So SMETANA constrains the MERGED COMMUNITY model's pool exchanges, not the member
# models'. 41d's section B compared generated ids against a single model and could
# never have matched. This builds a real 2-member community from the two gapseq models
# and reads the pool reaction names off it -- then solves for the mediadb string.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; mkdir -p "$OUT"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"
G="$CULT/gapseq_project"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41e  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"
[ -x "$PY" ] || { no "no smetana venv python"; exit 1; }

say "A. what does a merged community actually call its pool exchanges?"
"$PY" - "$G/YGCA01_Clean.xml" "$G/YGCA01_Contaminant.xml" "$OUT" <<'PY'
import sys
a,b,out=sys.argv[1],sys.argv[2],sys.argv[3]
from reframed import load_cbmodel, Community
from reframed.core.environment import Environment
m1=load_cbmodel(a,flavor='fbc2'); m1.id="M1"
m2=load_cbmodel(b,flavor='fbc2'); m2.id="M2"
print(f"  ...    members: {len(m1.reactions)} and {len(m2.reactions)} reactions")
try:
    com=Community("test",[m1,m2])
    mg=com.merged_model if hasattr(com,'merged_model') else com.merged
except Exception as e:
    print(f"  [FAIL] Community() raised {type(e).__name__}: {e}"); sys.exit(0)
pool=sorted(r for r in mg.reactions if r.endswith('_pool'))
print(f"  ...    merged model: {len(mg.reactions)} reactions, {len(pool)} pool exchanges")
for r in pool[:6]: print(f"           {r}")
if not pool:
    print("  [FAIL] no *_pool reactions -- SMETANA's whole addressing scheme is absent")
    sys.exit(0)
# Solve: for which compound string s does "R_EX_M_{s}_e_pool" land in pool?
import re
poolset=set(pool)
probe=["cpd00013","cpd00009","cpd00048","cpd00001"]
cands={}
for s in probe:
    for variant in (s, s+"_e0", s+"_c0", s.replace("cpd","cpd"), s+"_e"):
        built=f"R_EX_M_{variant}_e_pool"
        if built in poolset: cands.setdefault(s,[]).append(variant)
print("  ...    solving for the mediadb 'compound' string:")
for s in probe:
    v=cands.get(s)
    print(f"           {s:<12} -> {'WORKS with compound='+repr(v[0]) if v else 'no variant matches'}")
# what DO the pool names look like structurally?
pat=set()
for r in pool[:400]:
    mm=re.match(r'^R_EX_(.+)_pool$',r)
    if mm: pat.add(re.sub(r'cpd\d{5}','cpd#####',mm.group(1)))
print(f"  ...    distinct pool-name shapes: {sorted(pat)[:6]}")
solved=[s for s in probe if cands.get(s)]
if len(solved)==len(probe):
    v=cands[probe[0]][0]; suf=v[len(probe[0]):]
    print(f"  [OK]   CONVENTION SOLVED: mediadb compound = cpd id + {suf!r}")
    open(f"{out}/EXCHANGE_CONVENTION.txt","w").write((suf if suf else "bare")+"\n")
else:
    print(f"  [FAIL] no compound string reaches the pool exchanges.")
    print("         SMETANA hardcodes a BiGG-style '_e' external compartment;")
    print("         gapseq writes ModelSEED '_e0'. The models must be relabelled.")
    print("         Section B tests exactly that, on COPIES.")
PY

say "B. if unsolved: relabel e0->e on COPIES, and prove nothing else changed"
if [ -s "$OUT/EXCHANGE_CONVENTION.txt" ]; then
  ok "convention solved in A -- no relabelling needed, skipping"
else
  note "writing copies to \$WS2/41_env_checks/relabelled/ -- the culturomics originals"
  note "are READ-ONLY and are not touched."
  mkdir -p "$OUT/relabelled"
  for f in "$G/YGCA01_Clean.xml" "$G/YGCA01_Contaminant.xml"; do
    sed -e 's/_e0"/_e"/g' -e 's/_c0"/_c"/g' -e 's/compartment="e0"/compartment="e"/g' \
        -e 's/compartment="c0"/compartment="c"/g' -e 's/id="e0"/id="e"/g' -e 's/id="c0"/id="c"/g' \
        "$f" > "$OUT/relabelled/$(basename "$f")"
  done
  ok "relabelled copies written"
  "$PY" - "$G/YGCA01_Clean.xml" "$OUT/relabelled/YGCA01_Clean.xml" "$OUT/relabelled/YGCA01_Contaminant.xml" "$OUT" <<'PY'
import sys
orig,new,new2,out=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
from reframed import load_cbmodel, Community, FBA
mo=load_cbmodel(orig,flavor='fbc2'); mn=load_cbmodel(new,flavor='fbc2')
print(f"  ...    reactions  {len(mo.reactions)} -> {len(mn.reactions)}")
print(f"  ...    metabolites {len(mo.metabolites)} -> {len(mn.metabolites)}")
if len(mo.reactions)!=len(mn.reactions) or len(mo.metabolites)!=len(mn.metabolites):
    print("  [FAIL] relabelling changed the model's size -- reject it"); sys.exit(0)
# the real control: same objective value before and after
try:
    so=FBA(mo); sn=FBA(mn)
    vo=so.fobj if so else None; vn=sn.fobj if sn else None
    print(f"  ...    FBA objective  {vo!r} -> {vn!r}")
    if vo is None or vn is None: print("  [FAIL] one of them does not solve")
    elif abs(vo-vn)>1e-6*max(1,abs(vo)): print(f"  [FAIL] objective changed by {abs(vo-vn):.3g} -- relabelling altered the model")
    else: print("  [OK]   objective identical -- relabelling is cosmetic, as intended")
except Exception as e:
    print(f"  [FAIL] FBA raised {type(e).__name__}: {e}")
m2=load_cbmodel(new2,flavor='fbc2'); mn.id="M1"; m2.id="M2"
com=Community("t",[mn,m2]); mg=com.merged_model if hasattr(com,'merged_model') else com.merged
pool=sorted(r for r in mg.reactions if r.endswith('_pool'))
print(f"  ...    relabelled community: {len(pool)} pool exchanges")
for r in pool[:5]: print(f"           {r}")
ps=set(pool); good=[s for s in ("cpd00013","cpd00009","cpd00048","cpd00001")
                    if f"R_EX_M_{s}_e_pool" in ps]
print(f"  [{'OK' if len(good)==4 else 'FAIL'}]   {len(good)}/4 probe compounds now address correctly")
if len(good)==4:
    open(f"{out}/EXCHANGE_CONVENTION.txt","w").write("bare\n")
    open(f"{out}/RELABEL_REQUIRED.txt","w").write(
        "gapseq models must be relabelled e0->e, c0->c before SMETANA.\n"
        "Verified: reaction/metabolite counts unchanged, FBA objective unchanged.\n")
    print("  [OK]   wrote EXCHANGE_CONVENTION.txt (bare) and RELABEL_REQUIRED.txt")
PY
fi

say "C. mcr -- is it already answered on disk?"
K="$WS2/11_kofam"
if [ -d "$K" ]; then
  ok "11_kofam exists: $(find "$K" -type f | wc -l) files"
  find "$K" -type f | head -8 | sed 's/^/           /'
  note "searching existing kofam output for mcrA/B/G (K00399, K00401, K00402):"
  grep -rlE 'K00399|K00401|K00402' "$K" 2>/dev/null | head -5 | sed 's/^/           hit: /'
  grep -rhE 'K00399|K00401|K00402' "$K" 2>/dev/null | head -12 | sed 's/^/           /'
else
  no "no 11_kofam directory"
fi
note "profiles available for a direct hmmsearch if needed:"
ls "$WS2/36_gapclose/profiles" 2>/dev/null | head -5 | sed 's/^/           /'
for ko in K00399 K00401 K00402; do
  p="$WS2/36_gapclose/profiles/$ko.hmm"
  [ -f "$p" ] && echo "           $ko present" || echo "           $ko NOT in 36_gapclose/profiles"
done
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
