#!/usr/bin/env bash
# 41h -- why is mip/mro 'n/a', and does the medium reach the model at all?
#
# 41g's check counted rows in result_global.tsv and passed. The row said n/a for both
# scores. Counting rows tests the shape of an answer, not the answer. Every check below
# asserts on a VALUE, and the A/B in section C is designed so that one arm must differ
# from the other -- if both arms agree, the experiment failed, not the medium.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; S="$OUT/smoke_run"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41h  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"
[ -d "$S" ] || { no "$S missing -- run 41g first"; exit 1; }

say "A. read how smetana builds the community and applies the medium"
"$PY" - <<'PY'
import inspect, smetana.interface as I
src=inspect.getsource(I)
lines=src.splitlines()
hit=[i for i,l in enumerate(lines) if 'fmt_func' in l]
for h in hit:
    a,b=max(0,h-22),min(len(lines),h+10)
    print(f"  ...    interface.py around line {h+1}:")
    for l in lines[a:b]: print("           "+l)
    break
print("  ...    public callables in smetana.interface:")
print("           "+", ".join(n for n in dir(I) if not n.startswith('_') and callable(getattr(I,n))))
PY

say "B. does the community grow AT ALL, medium aside?"
"$PY" - "$S/M1.xml" "$S/M2.xml" <<'PY' 2>&1 | grep -v UserWarning | grep -v '^  warn'
import sys, warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, Community, FBA
m1=load_cbmodel(sys.argv[1],flavor='fbc2'); m1.id='M1'
m2=load_cbmodel(sys.argv[2],flavor='fbc2'); m2.id='M2'
for n,m in (('M1',m1),('M2',m2)):
    s=FBA(m); print(f"  ...    {n} alone, unconstrained: objective = {getattr(s,'fobj',None)!r}")
c=Community('t',[m1,m2])
mg=c.merged_model if hasattr(c,'merged_model') else c.merged
s=FBA(mg); print(f"  ...    merged community, unconstrained: objective = {getattr(s,'fobj',None)!r}")
ex=[r for r in mg.reactions if 'EX_' in r]
pool=[r for r in mg.reactions if r.endswith('_pool')]
print(f"  ...    merged: {len(mg.reactions)} reactions, {len(ex)} with EX_, {len(pool)} ending _pool")
print("  ...    a sample of the merged model's exchange-like reaction names:")
for r in sorted(ex)[:8]: print(f"           {r}")
if not pool:
    print("  [FAIL] no *_pool reactions. smetana's fmt_func R_EX_M_{}_e_pool cannot match")
    print("         ANYTHING here, so the medium is empty and growth is unconstrained-or-zero.")
    print("         That is a sufficient explanation for mip=n/a, mro=n/a.")
PY

say "C. A/B: with the medium vs without it. The arms MUST differ."
cd "$S" || exit 1
runarm(){ # $1 = label, rest = smetana args ; runs ONCE, values read back from the file
  local lab="$1"; shift
  rm -f "ab_${lab}"*.tsv
  timeout 900 "$SM" --flavor fbc2 --anaerobic -o "ab_${lab}" "$@" M1.xml M2.xml \
    > "ab_${lab}.log" 2>&1
}
val(){ [ -s "ab_$1_global.tsv" ] && awk -F'\t' 'NR==2{print $4"|"$5}' "ab_$1_global.tsv" || echo "NOFILE"; }
runarm nomedium
runarm withmedium -m SMOKE --mediadb media.tsv
A=$(val nomedium); B=$(val withmedium)
printf '           %-14s mip|mro = %s\n' "no medium"    "$A"
printf '           %-14s mip|mro = %s\n' "SMOKE medium" "$B"
echo
if [ "$A" = "$B" ]; then
  no "both arms identical ($A) -- the medium changed NOTHING."
  no "Either it never reached the model, or the failure is upstream of the medium."
else
  ok "arms differ -- the medium IS reaching the model"
fi
if [ "$A" = "n/a|n/a" ]; then
  no "even WITHOUT a medium the scores are n/a -- the medium is NOT the cause."
  note "reasons reported in ab_nomedium.log:"
  grep -iE 'error|infeasib|unable|cannot|fail|empty|no growth' ab_nomedium.log 2>/dev/null \
    | head -8 | sed 's/^/           /'
  [ -s ab_nomedium.log ] || note "(log empty)"
fi
note "detailed mode reports per-compound interactions; try it for a second view:"
timeout 900 "$SM" --flavor fbc2 --anaerobic -d -o det M1.xml M2.xml > det.log 2>&1
if [ -s det_detailed.tsv ]; then
  n=$(( $(wc -l < det_detailed.tsv) - 1 ))
  ok "det_detailed.tsv: $n interaction row(s)"
  head -4 det_detailed.tsv | sed 's/^/           /'
  [ "$n" -eq 0 ] && no "zero interactions between two models of the SAME organism is expected;"
  [ "$n" -eq 0 ] && note "this pair is a poor positive control -- see the note below."
else
  no "no detailed output"; tail -6 det.log | sed 's/^/           /'
fi
cd - >/dev/null

say "D. the control problem with this pair"
note "M1 = YGCA01_Clean, M2 = YGCA01_Contaminant. These are two bins from ONE isolate."
note "Two near-identical organisms should show LOW mip and mro -- so a zero here is"
note "ambiguous: it may be correct biology, not a broken pipeline."
note "A real positive control needs two metabolically DIFFERENT organisms. Once any two"
note "MAG models exist, use a fermenter and a methanogen -- they must interact."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
