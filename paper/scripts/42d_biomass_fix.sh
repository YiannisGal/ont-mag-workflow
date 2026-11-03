#!/usr/bin/env bash
# 42d -- three MAGs have no recoverable 16S, so gapseq draft cannot infer Gram staining
# and refuses. Supply the biomass class from GTDB-Tk phylum instead, which is better
# evidence than a 16S guess and is already in MAG_TAXONOMY.tsv.
#
# The flag spelling is DISCOVERED, not assumed: 42b's detector already produced false
# negatives for draft and fill by parsing an R traceback as help text.
set -uo pipefail
WS2="${WS2:-${WS2}}"
GS="${GS:-$WS2/envs/gapseq_local}"; OUT="$WS2/42_gems"
FILLMED="$WS2/41_env_checks/ETO_FILL-medium.csv"
export TMPDIR="${TMPDIR42:-${SCRATCH_ROOT}/tmp42}"; mkdir -p "$TMPDIR"
T="${STAGE_TIMEOUT:-21600}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "42d  $(date -Is)"
PIXI=$(find "$HOME/Documents/UniqueELT_Solid" -maxdepth 6 -type d -name default -path '*envs*' \
       2>/dev/null | grep -v damaged | while read -r e; do [ -x "$e/bin/gapseq" ] && echo "$e"; done)
[ -n "$PIXI" ] && export PATH="$PIXI/bin:$PATH"
SRC="$GS/src/generate_GSdraft.R"

say "A. how does draft spell the biomass option, and what values does it take?"
[ -f "$SRC" ] || { no "no $SRC"; exit 1; }
note "option definition:"
grep -nE '\-\-biomass|"biomass"|opt\$biomass' "$SRC" | head -10 | sed 's/^/           /'
note "accepted values (strings compared against biomass in the source):"
grep -noE '"(Gram_pos|Gram_neg|gram_pos|gram_neg|archaea|Archaea|auto|pos|neg)"' "$SRC" \
  | sort -u -t: -k2 | head -12 | sed 's/^/           /'

say "B. which genomes still need a model, and what class each should get"
python3 - "$WS2" "$OUT" <<'PY'
import sys,os,csv,glob
ws2,out=sys.argv[1],sys.argv[2]
tp=min(glob.glob(os.path.join(ws2,'**','MAG_TAXONOMY.tsv'),recursive=True),key=lambda p:p.count(os.sep))
tax=list(csv.DictReader(open(tp),delimiter='\t'))
nk=next(k for k in tax[0] if k.strip().lower()=='working_name')
pk=next(k for k in tax[0] if k.strip().lower()=='phylum')
dk=next(k for k in tax[0] if k.strip().lower() in ('domain','kingdom'))
gk=next(k for k in tax[0] if k.strip().lower()=='genus')
info={r[nk]:(r[pk],r[dk],r[gk]) for r in tax}
# Gram character by phylum. Monoderm phyla are Gram-positive; the rest diderm.
POS={'bacillota','firmicutes','firmicutes_a','firmicutes_b','firmicutes_c',
     'actinobacteriota','actinobacteria','chloroflexota','deinococcota'}
rows=[]
for line in open(os.path.join(out,'targets.tsv')):
    w,d,f=line.rstrip('\n').split('\t')
    have=glob.glob(os.path.join(out,w,'*.xml'))
    have=[h for h in have if not h.endswith('-draft.xml') and os.path.getsize(h)>0]
    if have: continue
    ph,dom,gen=info.get(w,('?','?','?'))
    if 'archae' in (dom or '').lower(): cls='archaea'
    elif (ph or '').strip().lower() in POS: cls='Gram_pos'
    else: cls='Gram_neg'
    rows.append((w,f,cls,ph,gen))
    print(f"           {w:<24}{gen:<20}{ph:<20}-> {cls}")
if not rows: print("  [OK]   nothing left to build"); sys.exit(0)
print(f"  ...    {len(rows)} genome(s) to retry")
open(os.path.join(out,'retry.tsv'),'w').write(
    "".join(f"{w}\t{f}\t{c}\n" for w,f,c,_,_ in rows))
PY
[ -s "$OUT/retry.tsv" ] || { ok "nothing to do"; exit 0; }

say "C. find the working invocation on ONE genome before touching the rest"
read -r W F CLS < <(head -1 "$OUT/retry.tsv" | tr '\t' ' ')
B=$(basename "$F"); B="${B%.*}"; D="$OUT/$W"
note "test genome: $W  (class $CLS)"
WORKING=""
for spec in "--biomass $CLS" "-b $CLS"; do
  rm -f "$D/$B-draft.RDS"
  note "trying: gapseq draft $spec ..."
  ( cd "$D" && timeout 900 "$GS/gapseq" draft -r "$D/$B-all-Reactions.tbl" \
      -t "$D/$B-Transporter.tbl" -c "$F" -u 200 -l 100 $spec ) >"$D/try.log" 2>&1
  rc=$?
  if [ -s "$D/$B-draft.RDS" ]; then
    ok "SUCCESS with '$spec' (exit=$rc)"; WORKING="$spec"; break
  else
    no "'$spec' produced no RDS (exit=$rc)"; tail -4 "$D/try.log" | sed 's/^/             /'
  fi
done
if [ -z "$WORKING" ]; then
  no "no invocation produced a draft. Not proceeding to the other genomes."
  note "full log: $D/try.log"; exit 1
fi

say "D. apply to all, then fill"
cat > "$OUT/retry_one.sh" <<'INNER'
#!/usr/bin/env bash
set -uo pipefail
W="$1"; F="$2"; CLS="$3"; GS="$4"; OUT="$5"; MED="$6"; T="$7"; SPEC="$8"
D="$OUT/$W"; cd "$D" || exit 1; L="$D/gapseq.log"; B=$(basename "$F"); B="${B%.*}"
SPEC="${SPEC/\$CLS/$CLS}"
if [ ! -s "$D/$B-draft.RDS" ]; then
  echo "[$W] draft($CLS) start $(date -Is)" | tee -a "$L"; t0=$SECONDS
  # shellcheck disable=SC2086
  timeout "$T" "$GS/gapseq" draft -r "$D/$B-all-Reactions.tbl" -t "$D/$B-Transporter.tbl" \
      -c "$F" -u 200 -l 100 $SPEC >>"$L" 2>&1; rc=$?
  echo "[$W] draft exit=$rc $((SECONDS-t0))s" | tee -a "$L"
  [ -s "$D/$B-draft.RDS" ] || { echo "[$W] draft wrote no RDS" | tee -a "$L"; exit 12; }
fi
model(){ find "$D" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1; }
if [ -z "$(model)" ]; then
  echo "[$W] fill start $(date -Is)" | tee -a "$L"; t0=$SECONDS
  timeout "$T" "$GS/gapseq" fill -m "$D/$B-draft.RDS" -n "$MED" \
      -c "$D/$B-rxnWeights.RDS" -g "$D/$B-rxnXgenes.RDS" >>"$L" 2>&1; rc=$?
  echo "[$W] fill exit=$rc $((SECONDS-t0))s" | tee -a "$L"
  [ -z "$(model)" ] && { echo "[$W] fill wrote no model" | tee -a "$L"; exit 13; }
fi
echo "[$W] COMPLETE model=$(basename "$(model)")" | tee -a "$L"
INNER
chmod +x "$OUT/retry_one.sh"
# substitute the discovered flag name, keeping each genome's own class
FLAG=$(printf '%s' "$WORKING" | awk '{print $1}')
awk -F'\t' -v g="$GS" -v o="$OUT" -v m="$FILLMED" -v t="$T" -v fl="$FLAG" \
  '{print $1"\t"$2"\t"$3"\t"g"\t"o"\t"m"\t"t"\t"fl" "$3}' "$OUT/retry.tsv" \
  | xargs -P 3 -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -ra x <<< "{}"; '"$OUT"'/retry_one.sh "${x[@]}"'

say "E. all ten"
PY="$WS2/envs/smetana/bin/python"; D=0; N=0
while IFS=$'\t' read -r w dom f; do
  N=$((N+1)); MX=$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 | head -1)
  if [ -z "$MX" ]; then printf '  [FAIL] %-24s %-9s no model\n' "$w" "$dom"; continue; fi
  D=$((D+1))
  "$PY" - "$MX" "$w" "$dom" <<'PY' 2>/dev/null | grep -v Warning
import sys,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
p,w,d=sys.argv[1],sys.argv[2],sys.argv[3]
m=load_cbmodel(p,flavor='fbc2'); ex=[r for r in m.reactions if 'EX_' in r]
g=getattr(FBA(m),'fobj',None)
print(f"  [{'OK' if (g or 0)>1e-9 else 'FAIL'}]   {w:<24}{d:<9}{len(m.reactions):>5} rxn {len(ex):>4} exch  growth={g!r}")
PY
done < "$OUT/targets.tsv"
echo "           $D/$N models"
[ "$D" -eq "$N" ] && ok "the in-situ ten are modelled -- SMETANA can run" \
                  || no "still incomplete"
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
