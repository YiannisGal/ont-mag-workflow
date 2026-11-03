#!/usr/bin/env bash
# 43b -- all ten models grow 0.0000 alone on ETO_SED, so the positive control's zero
# was trivial. Find out what is actually blocking growth, then refill on the analysis
# medium.
#
# 🔴 MY ERROR, stated plainly. 41k said "fill on ETO_SED, analyse on ETO_SED". 41l
# reversed it to "fill on ETO_FILL" on the argument that filling on the analysis medium
# would make every member self-sufficient and erase cross-feeding. The reversal produced
# models that cannot grow on ETO_SED AT ALL, which erases everything instead.
# The argument was wrong because gapseq's fill is MINIMAL: it adds the fewest reactions
# that permit growth, so a model filled on ETO_SED is minimally viable, not lavishly
# self-sufficient -- and MIP is defined exactly to measure what interaction saves
# relative to that.
#
# 43 section B was also a one-sided control: it warned if every member was self-
# sufficient and said nothing about every member being dead. It reported 10/10 at zero
# as "the question is well posed".
set -uo pipefail
WS2="${WS2:-${WS2}}"
GEMS="$WS2/42_gems"; ENV="$WS2/41_env_checks"; OUT="$WS2/43_smetana"
GS="${GS:-$WS2/envs/gapseq_local}"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"
T="${STAGE_TIMEOUT:-7200}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "43b  $(date -Is)"
PIXI=$(find "$HOME/Documents/UniqueELT_Solid" -maxdepth 6 -type d -name default -path '*envs*' \
       2>/dev/null | grep -v damaged | while read -r e; do [ -x "$e/bin/gapseq" ] && echo "$e"; done)
[ -n "$PIXI" ] && export PATH="$PIXI/bin:$PATH"

say "A. what exactly is blocking growth? add one group at a time"
"$PY" - "$OUT/models" "$ENV" "$GEMS/targets.tsv" <<'PY' 2>&1 | grep -v Warning
import sys,os,csv,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
from reframed.core.environment import Environment
mdir,env,tg=sys.argv[1],sys.argv[2],sys.argv[3]
def ids(p): return [r['compounds'].split('[')[0].strip() for r in csv.DictReader(open(p))]
SED=ids(os.path.join(env,'ETO_SED-medium.csv'))
FILL=ids(os.path.join(env,'ETO_FILL-medium.csv'))
AA={"cpd00023","cpd00033","cpd00035","cpd00039","cpd00041","cpd00051","cpd00053","cpd00054",
    "cpd00060","cpd00065","cpd00066","cpd00069","cpd00084","cpd00107","cpd00119","cpd00129",
    "cpd00132","cpd00156","cpd00161","cpd00322"}
CUR={"cpd00029","cpd11640","cpd00047"}; ROUTE={"cpd00116","cpd00441","cpd00098"}; B12={"cpd00166"}
tests=[("ETO_SED",SED),("SED+aminoacids",SED+[c for c in FILL if c in AA]),
       ("SED+currencies",SED+[c for c in FILL if c in CUR]),
       ("SED+B12",SED+[c for c in FILL if c in B12]),
       ("SED+route",SED+[c for c in FILL if c in ROUTE]),("ETO_FILL",FILL)]
names=[l.split('\t')[0] for l in open(tg)]
print("           "+f"{'genome':<24}"+"".join(f"{t[0]:>16}" for t in tests))
grew={t[0]:0 for t in tests}
for w in names:
    p=os.path.join(mdir,w+'.xml')
    if not os.path.exists(p): continue
    row=f"           {w:<24}"
    for lab,cp in tests:
        m=load_cbmodel(p,flavor='fbc2')
        e=Environment.from_compounds(sorted(set(cp)), fmt_func=lambda x: f"R_EX_{x}_e")
        e.apply(m,inplace=True,warning=False)
        g=getattr(FBA(m),'fobj',None) or 0.0
        if g>1e-9: grew[lab]+=1
        row+=f"{g:>16.4f}"
    print(row)
print()
for lab,_ in tests: print(f"  ...    {lab:<16} {grew[lab]}/{len(names)} genomes grow")
if grew["ETO_FILL"]==0:
    print("  [FAIL] nothing grows even on the medium it was FILLED on. That is a broken")
    print("         model set, not a medium problem. Stop and investigate.")
elif grew["ETO_SED"]==0 and grew["ETO_FILL"]>0:
    print("  [OK]   diagnosis confirmed: the models are viable only on their fill medium.")
    print("         Refilling on ETO_SED is the fix, and section B does it.")
PY

say "B. refill every draft on the ANALYSIS medium"
note "drafts already exist, so this is fill only (~12 min each, 5 at a time)"
note "the ETO_FILL models are kept at 43_smetana/models_fillrich for comparison"
mkdir -p "$OUT/models_fillrich"; cp "$OUT/models"/*.xml "$OUT/models_fillrich/" 2>/dev/null || true
cat > "$OUT/refill_one.sh" <<'INNER'
#!/usr/bin/env bash
set -uo pipefail
W="$1"; F="$2"; GS="$3"; GEMS="$4"; MED="$5"; T="$6"
D="$GEMS/$W"; cd "$D" || exit 1; B=$(basename "$F"); B="${B%.*}"
L="$D/refill.log"; OUTD="$D/sed_fill"; mkdir -p "$OUTD"
[ -s "$OUTD/$B.xml" ] && { echo "[$W] skip refill"; exit 0; }
echo "[$W] refill start $(date -Is)" | tee -a "$L"; t0=$SECONDS
timeout "$T" "$GS/gapseq" fill -m "$D/$B-draft.RDS" -n "$MED" \
    -c "$D/$B-rxnWeights.RDS" -g "$D/$B-rxnXgenes.RDS" -f "$OUTD" >>"$L" 2>&1
rc=$?; echo "[$W] refill exit=$rc $((SECONDS-t0))s" | tee -a "$L"
M=$(find "$OUTD" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1)
[ -z "$M" ] && { echo "[$W] refill wrote no model" | tee -a "$L"; exit 1; }
echo "[$W] refilled: $(basename "$M")" | tee -a "$L"
INNER
chmod +x "$OUT/refill_one.sh"
awk -F'\t' -v g="$GS" -v m="$GEMS" -v med="$ENV/ETO_SED-medium.csv" -v t="$T" \
  '{print $1"\t"$3"\t"g"\t"m"\t"med"\t"t}' "$GEMS/targets.tsv" \
  | xargs -P 5 -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -ra x <<< "{}"; '"$OUT"'/refill_one.sh "${x[@]}"'

say "C. relabel the refilled models and check growth on ETO_SED"
M2="$OUT/models_sed"; mkdir -p "$M2"
while IFS=$'\t' read -r w dom f; do
  MX=$(find "$GEMS/$w/sed_fill" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null | head -1)
  [ -z "$MX" ] && { no "$w: no refilled model"; continue; }
  sed -e 's/_e0"/_e"/g' -e 's/_c0"/_c"/g' -e 's/compartment="e0"/compartment="e"/g' \
      -e 's/compartment="c0"/compartment="c"/g' -e 's/id="e0"/id="e"/g' -e 's/id="c0"/id="c"/g' \
      "$MX" > "$M2/$w.xml"
done < "$GEMS/targets.tsv"
ok "$(ls "$M2"/*.xml 2>/dev/null | wc -l) refilled models"
"$PY" - "$M2" "$ENV/ETO_SED-medium.csv" "$GEMS/targets.tsv" "$OUT/models_fillrich" <<'PY' 2>&1 | grep -v Warning
import sys,os,csv,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
from reframed.core.environment import Environment
m2,med,tg,old=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
cpds=[r['compounds'].split('[')[0].strip() for r in csv.DictReader(open(med))]
alive=0; tot=0
for line in open(tg):
    w,dom,_=line.rstrip('\n').split('\t'); p=os.path.join(m2,w+'.xml')
    if not os.path.exists(p): continue
    tot+=1
    m=load_cbmodel(p,flavor='fbc2')
    n_new=len(m.reactions)
    op=os.path.join(old,w+'.xml'); n_old=len(load_cbmodel(op,flavor='fbc2').reactions) if os.path.exists(op) else 0
    e=Environment.from_compounds(cpds, fmt_func=lambda x: f"R_EX_{x}_e")
    e.apply(m,inplace=True,warning=False)
    g=getattr(FBA(m),'fobj',None) or 0.0
    if g>1e-9: alive+=1
    print(f"  [{'OK' if g>1e-9 else 'FAIL'}]   {w:<24}{dom:<9}{n_old:>5} -> {n_new:>5} rxn   growth on ETO_SED={g:.4f}")
print(f"\n  [{'OK' if alive==tot else 'FAIL'}]   {alive}/{tot} grow on the analysis medium")
if alive<tot:
    print("         Members that still cannot grow cannot enter a SMETANA run.")
PY

say "D. positive control, again"
ARCH=$(awk -F'\t' '$2=="Archaea"{print $1}' "$GEMS/targets.tsv" | head -1)
FERM=$(awk -F'\t' '$2=="Bacteria"{print $1}' "$GEMS/targets.tsv" | head -1)
P="$OUT/poscontrol2"; rm -rf "$P"; mkdir -p "$P"
cp "$M2/$ARCH.xml" "$M2/$FERM.xml" "$P/" 2>/dev/null
cp "$ENV/etoliko_media.tsv" "$P/"
note "$ARCH + $FERM on ETO_SED"
( cd "$P" && timeout "$T" "$SM" --flavor fbc2 --anaerobic -d -m ETO_SED \
    --mediadb etoliko_media.tsv -o pc "$ARCH.xml" "$FERM.xml" ) >"$P/pc.log" 2>&1
F="$P/pc_detailed.tsv"
if [ -s "$F" ] && [ "$(( $(wc -l < "$F") - 1 ))" -gt 0 ]; then
  ok "$(( $(wc -l < "$F") - 1 )) predicted exchange(s)"
  head -1 "$F" | sed 's/^/           /'
  sort -t$'\t' -k8 -gr "$F" | head -12 | sed 's/^/           /'
  ok "SETUP PROVEN -- runs 1-3 can proceed"
else
  no "still zero. Do not proceed."; tail -10 "$P/pc.log" | sed 's/^/           /'
fi
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
