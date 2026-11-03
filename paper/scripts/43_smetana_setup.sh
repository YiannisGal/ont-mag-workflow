#!/usr/bin/env bash
# 43 -- prepare and PROVE the SMETANA setup before running anything that produces a
# number for the thesis.
#
# A. relabel e0->e on copies, with the control that caught nothing last time but must
#    still run: reaction/metabolite counts and FBA objective unchanged.
# B. growth of each model ALONE on ETO_SED. The methanogen is expected to FAIL -- that
#    is the dependency the chapter is about, not a bug.
# C. POSITIVE CONTROL: methanogen + fermenter, detailed mode. They MUST interact.
#    41h's control pair was two bins of one isolate, so a zero there was uninformative.
#    This pair is informative in both directions.
set -uo pipefail
WS2="${WS2:-${WS2}}"
GEMS="$WS2/42_gems"; ENV="$WS2/41_env_checks"; OUT="$WS2/43_smetana"
M="$OUT/models"; mkdir -p "$M"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"
T="${STAGE_TIMEOUT:-7200}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "43  $(date -Is)"

say "A. relabel to the compartment names SMETANA addresses"
while IFS=$'\t' read -r w dom f; do
  MX=$(find "$GEMS/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1)
  [ -z "$MX" ] && { no "$w has no model"; continue; }
  sed -e 's/_e0"/_e"/g' -e 's/_c0"/_c"/g' -e 's/compartment="e0"/compartment="e"/g' \
      -e 's/compartment="c0"/compartment="c"/g' -e 's/id="e0"/id="e"/g' -e 's/id="c0"/id="c"/g' \
      "$MX" > "$M/$w.xml"
done < "$GEMS/targets.tsv"
ok "$(ls "$M"/*.xml 2>/dev/null | wc -l) copies written to $M"
note "control: counts and FBA objective must be identical before and after"
"$PY" - "$GEMS" "$M" <<'PY' 2>&1 | grep -v Warning
import sys,os,glob,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
gems,mdir=sys.argv[1],sys.argv[2]
bad=0
for line in open(os.path.join(gems,'targets.tsv')):
    w=line.split('\t')[0]
    orig=[p for p in glob.glob(os.path.join(gems,w,'*.xml')) if not p.endswith('-draft.xml')]
    new=os.path.join(mdir,w+'.xml')
    if not orig or not os.path.exists(new): continue
    a=load_cbmodel(orig[0],flavor='fbc2'); b=load_cbmodel(new,flavor='fbc2')
    ga=getattr(FBA(a),'fobj',None); gb=getattr(FBA(b),'fobj',None)
    same=(len(a.reactions)==len(b.reactions) and len(a.metabolites)==len(b.metabolites)
          and ga is not None and gb is not None and abs(ga-gb)<=1e-9*max(1,abs(ga)))
    if not same: bad+=1
    print(f"  [{'OK' if same else 'FAIL'}]   {w:<24}{len(a.reactions)}->{len(b.reactions)} rxn  "
          f"obj {ga!r} -> {gb!r}")
print(f"  [{'OK' if bad==0 else 'FAIL'}]   relabelling is cosmetic for {'all' if bad==0 else 'NOT all'} models")
PY

say "B. can each model grow ALONE on ETO_SED?"
note "filled on ETO_FILL, analysed on ETO_SED. A member that cannot grow alone here"
note "is DEPENDENT -- that is the finding, not a failure."
"$PY" - "$M" "$ENV/ETO_SED-medium.csv" "$GEMS/targets.tsv" <<'PY' 2>&1 | grep -v Warning
import sys,os,csv,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
from reframed.core.environment import Environment
mdir,med,tg=sys.argv[1],sys.argv[2],sys.argv[3]
cpds=[r['compounds'].split('[')[0].strip() for r in csv.DictReader(open(med))]
dep=[];ind=[]
for line in open(tg):
    w,dom,_=line.rstrip('\n').split('\t')
    p=os.path.join(mdir,w+'.xml')
    if not os.path.exists(p): continue
    m=load_cbmodel(p,flavor='fbc2')
    env=Environment.from_compounds(cpds, fmt_func=lambda x: f"R_EX_{x}_e")
    hit=sum(1 for k in env if k in m.reactions)
    env.apply(m,inplace=True,warning=False)
    g=getattr(FBA(m),'fobj',None) or 0.0
    (ind if g>1e-9 else dep).append(w)
    print(f"  ...    {w:<24}{dom:<9}{hit:>3}/{len(env)} exch matched   growth={g:.4f}"
          + ("" if g>1e-9 else "   <-- DEPENDENT"))
print(f"\n  ...    {len(ind)} self-sufficient, {len(dep)} dependent on ETO_SED")
if dep: print("         dependent: "+", ".join(dep))
if not dep:
    print("  [WARN] EVERY member is self-sufficient. Cross-feeding scores will be near")
    print("         zero by construction, and that would be an artefact of the medium,")
    print("         not a finding. Check that ETO_SED is not too rich.")
else:
    print("  [OK]   at least one member requires a partner -- the question is well posed")
PY

say "C. POSITIVE CONTROL -- methanogen + fermenter must interact"
ARCH=$(awk -F'\t' '$2=="Archaea"{print $1}' "$GEMS/targets.tsv" | head -1)
FERM=$(awk -F'\t' '$2=="Bacteria"{print $1}' "$GEMS/targets.tsv" | head -1)
note "pair: $ARCH (methanogen) + $FERM (fermenter)"
note "the methanogen has no H2, formate or acetate on ETO_SED unless a partner supplies it"
P="$OUT/poscontrol"; rm -rf "$P"; mkdir -p "$P"
cp "$M/$ARCH.xml" "$M/$FERM.xml" "$P/"
cp "$ENV/etoliko_media.tsv" "$P/"
( cd "$P" && timeout "$T" "$SM" --flavor fbc2 --anaerobic --molweight -d \
    -m ETO_SED --mediadb etoliko_media.tsv -o pc "$ARCH.xml" "$FERM.xml" ) >"$P/pc.log" 2>&1
rc=$?
echo "           exit=$rc"
F="$P/pc_detailed.tsv"
if [ -s "$F" ]; then
  N=$(( $(wc -l < "$F") - 1 ))
  if [ "$N" -gt 0 ]; then
    ok "$N predicted metabolite exchange(s)"
    head -1 "$F" | sed 's/^/           /'
    sort -t$'\t' -k8 -gr "$F" 2>/dev/null | head -12 | sed 's/^/           /'
    note "any row with $ARCH as receiver is the methanogen being fed. That is the"
    note "mechanism the chapter claims, reproduced in silico."
  else
    no "ZERO exchanges between a methanogen and a fermenter on a medium with no H2,"
    no "formate or acetate. That is not biology -- something in the setup is wrong."
    tail -12 "$P/pc.log" | sed 's/^/           /'
  fi
else
  no "no detailed output"; tail -15 "$P/pc.log" | sed 's/^/           /'
fi

say "D. verdict"
if [ -s "$F" ] && [ "$(( $(wc -l < "$F") - 1 ))" -gt 0 ]; then
  ok "setup proven on a pair whose answer is known in advance. Runs 1-3 can proceed."
else
  no "do NOT run the ten-genome analysis until the positive control passes."
fi
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
