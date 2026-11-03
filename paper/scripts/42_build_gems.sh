#!/usr/bin/env bash
# 42 -- build genome-scale models for the in-situ ten. Long-running; resumable.
#
# gapseq writes its outputs to the CURRENT DIRECTORY. Running 8 genomes concurrently in
# one directory would have them overwrite each other's -all-Reactions.tbl and produce a
# silently scrambled result. Each genome therefore gets its own working directory.
#
# -O is mandatory: it forces offline mode. Its absence (and my misreading of -U) cost a
# day and rewrote 2,435 reference files in a shared environment. -U is never passed.
#
# Every stage runs under a wall-clock budget and every stage is skipped if its output
# already exists, so the script can be re-run after an interruption without redoing work.
set -uo pipefail
WS2="${WS2:-${WS2}}"
GS="${GS:-$WS2/envs/gapseq_local}"
OUT="$WS2/42_gems"; ENVCK="$WS2/41_env_checks"
FILLMED="$ENVCK/ETO_FILL-medium.csv"
export TMPDIR="${TMPDIR42:-${SCRATCH_ROOT}/tmp42}"
STAGE_TIMEOUT="${STAGE_TIMEOUT:-21600}"     # 6 h per stage
JOBS="${JOBS:-8}"; KTHREADS="${KTHREADS:-3}"
mkdir -p "$OUT" "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "42_build_gems  $(date -Is)"

say "A. preflight -- refuse to start on a bad setup"
FATAL=0
[ -x "$GS/gapseq" ] || { no "no gapseq at $GS/gapseq"; FATAL=1; }
[ -d "$GS/dat" ]    || { no "no $GS/dat -- gapseq would fall back to a shared install"; FATAL=1; }
[ -s "$FILLMED" ]   || { no "no fill medium at $FILLMED (run 41m)"; FATAL=1; }
if grep -qE 'cpd00007|cpd00209|cpd00075' "$FILLMED" 2>/dev/null; then
  no "FILL MEDIUM CONTAINS AN OXIDANT -- refusing"; FATAL=1
else ok "fill medium is anoxic"; fi
# gapseq needs R and the aligners; those live in a pixi env. FOUR envs are named
# "default" -- 38b already learned this and the rule is to filter by bin/gapseq, which
# is the only positive evidence that an env is the gapseq one. Picking by sort order
# once selected checkm2_sandbox.
mapfile -t PIXIS < <(find "$HOME/Documents/UniqueELT_Solid" -maxdepth 6 -type d -name default \
        -path '*envs*' 2>/dev/null | grep -v damaged | sort)
note "${#PIXIS[@]} candidate pixi env(s); selecting by presence of bin/gapseq:"
SEL=()
for e in "${PIXIS[@]}"; do
  if [ -x "$e/bin/gapseq" ]; then SEL+=("$e"); printf '           KEEP   %s\n' "$e"
  else printf '           reject %s  (no bin/gapseq)\n' "$e"; fi
done
if [ "${#SEL[@]}" -eq 0 ]; then
  no "no pixi env provides bin/gapseq"; FATAL=1
elif [ "${#SEL[@]}" -gt 1 ]; then
  no "${#SEL[@]} envs provide bin/gapseq -- refusing to choose"; FATAL=1
else
  PIXI="${SEL[0]}"; ok "pixi env: $PIXI"; export PATH="$PIXI/bin:$PATH"
fi

# dependency check AFTER PATH is set, and reporting where each one resolves from.
if [ "$FATAL" -eq 0 ]; then
  ESSENTIAL="Rscript blastp exonerate bedtools hmmsearch"
  OPTIONAL="barrnap diamond"
  for b in $ESSENTIAL; do
    w=$(command -v "$b" 2>/dev/null)
    if [ -n "$w" ]; then printf '           %-12s %s\n' "$b" "$w"
    else no "missing ESSENTIAL dependency: $b"; FATAL=1; fi
  done
  for b in $OPTIONAL; do
    w=$(command -v "$b" 2>/dev/null)
    if [ -n "$w" ]; then printf '           %-12s %s\n' "$b" "$w"
    else
      note "optional dependency absent: $b"
      [ "$b" = barrnap ] && note "  barrnap gives gapseq 16S-based pathway hints; without it the" \
        && note "  reaction set is built from protein homology alone. Not fatal, but RECORD IT."
    fi
  done
fi
[ "$FATAL" -eq 0 ] && ok "all dependencies present" || { no "preflight failed"; exit 1; }
df -h "$WS2" | tail -1 | sed 's/^/           /'

say "B. resolve the genome list"
python3 - "$WS2" "$ENVCK" "$OUT" <<'PY'
import sys,os,csv
ws2,envck,out=sys.argv[1],sys.argv[2],sys.argv[3]
ins=os.path.join(envck,'in_situ_10.txt')
if not os.path.exists(ins): print("  [FAIL] no in_situ_10.txt"); sys.exit(1)
ids=[l.strip() for l in open(ins) if l.strip()]
import glob
cands=sorted(glob.glob(os.path.join(ws2,'**','MAG_CROSSWALK.tsv'),recursive=True))
if not cands:
    print("  [FAIL] MAG_CROSSWALK.tsv not found anywhere under $WS2"); sys.exit(1)
if len(cands)>1:
    print(f"  ...    {len(cands)} crosswalk files found; using the shallowest:")
    for c in cands: print(f"           {c}")
xw=min(cands,key=lambda p:p.count(os.sep))
print(f"  ...    crosswalk: {xw}")
rows=list(csv.DictReader(open(xw),delimiter='\t'))
ec=next(k for k in rows[0] if any(str(r[k]).startswith('ETO_') for r in rows))
tcands=sorted(glob.glob(os.path.join(ws2,'**','MAG_TAXONOMY.tsv'),recursive=True))
if not tcands:
    print("  [FAIL] MAG_TAXONOMY.tsv not found"); sys.exit(1)
tpath=min(tcands,key=lambda p:p.count(os.sep))
print(f"  ...    taxonomy:  {tpath}")
tax=list(csv.DictReader(open(tpath),delimiter='\t'))
nk=next(k for k in tax[0] if k.strip().lower()=='working_name')
dk=next(k for k in tax[0] if k.strip().lower() in ('domain','kingdom'))
names={r[nk] for r in tax}
wc=next(k for k in rows[0] if any(str(r[k]) in names for r in rows))
m={r[ec]:r[wc] for r in rows}
dom={r[nk]:('Archaea' if 'archae' in (r[dk] or '').lower() else 'Bacteria') for r in tax}
bins=os.path.join(ws2,'07_mags','bins')
if not os.path.isdir(bins):
    cand=[d for d in glob.glob(os.path.join(ws2,'07_mags','*')) if os.path.isdir(d)
          and len([f for f in os.listdir(d) if f.endswith(('.fa','.fna','.fasta'))])>=len(names)]
    if cand: bins=cand[0]
    print(f"  ...    genome dir: {bins}")
lines=[]; miss=[]
for i in ids:
    w=m.get(i)
    if not w: miss.append(i); continue
    fa=None
    for ext in ('.fa','.fna','.fasta'):
        p=os.path.join(bins,w+ext)
        if os.path.exists(p): fa=p; break
    if not fa: miss.append(f"{i}/{w} (no FASTA)"); continue
    lines.append(f"{w}\t{dom.get(w,'Bacteria')}\t{fa}")
print(f"  ...    {len(ids)} in-situ ids -> {len(lines)} resolved, {len(miss)} unresolved")
for l in lines:
    w,d,f=l.split('\t'); print(f"           {w:<26}{d:<10}{os.path.basename(f)}")
if miss:
    print("  [FAIL] unresolved: "+", ".join(miss)); sys.exit(1)
na=sum(1 for l in lines if l.split('\t')[1]=='Archaea')
print(f"  [OK]   {len(lines)} genomes, {na} archaea (-m Archaea), {len(lines)-na} bacteria")
open(os.path.join(out,'targets.tsv'),'w').write("\n".join(lines)+"\n")
PY
[ -s "$OUT/targets.tsv" ] || { no "no targets"; exit 1; }

say "C. build (JOBS=$JOBS, -K $KTHREADS, STAGE_TIMEOUT=${STAGE_TIMEOUT}s)"
note "each genome runs in its own directory; gapseq writes to CWD"
cat > "$OUT/one_genome.sh" <<'INNER'
#!/usr/bin/env bash
set -uo pipefail
W="$1"; DOM="$2"; FA="$3"; GS="$4"; OUT="$5"; MED="$6"; T="$7"; K="$8"
D="$OUT/$W"; mkdir -p "$D"; cd "$D" || exit 1
L="$D/gapseq.log"; : > "$L"
step(){ # name, output-to-check, command...
  local nm="$1" chk="$2"; shift 2
  if [ -s "$chk" ]; then echo "[$W] skip $nm (exists)" | tee -a "$L"; return 0; fi
  echo "[$W] $nm start $(date -Is)" | tee -a "$L"
  local t0=$SECONDS
  timeout "$T" "$@" >>"$L" 2>&1; local rc=$?
  echo "[$W] $nm exit=$rc  $((SECONDS-t0))s" | tee -a "$L"
  [ $rc -ne 0 ] && return $rc
  [ -s "$chk" ] || { echo "[$W] $nm produced no $chk" | tee -a "$L"; return 3; }
  return 0
}
B=$(basename "$FA"); B="${B%.*}"
step find      "$D/$B-all-Reactions.tbl" "$GS/gapseq" find -O -p all -b 200 -m "$DOM" -K "$K" "$FA" || exit 10
step transport "$D/$B-Transporter.tbl"   "$GS/gapseq" find-transport -O -K "$K" "$FA"              || exit 11
step draft     "$D/$B-draft.RDS"         "$GS/gapseq" draft -r "$D/$B-all-Reactions.tbl" \
      -t "$D/$B-Transporter.tbl" -c "$FA" -u 200 -l 100                                            || exit 12
# gapseq's fill output filename is NOT assumed. Detect any non-draft .xml instead.
model(){ find "$D" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null | head -1; }
if [ -n "$(model)" ]; then
  echo "[$W] skip fill (model exists: $(basename "$(model)"))" | tee -a "$L"
else
  echo "[$W] fill start $(date -Is)" | tee -a "$L"; t0=$SECONDS
  timeout "$T" "$GS/gapseq" fill -m "$D/$B-draft.RDS" -n "$MED" \
      -c "$D/$B-rxnWeights.RDS" -g "$D/$B-rxnXgenes.RDS" >>"$L" 2>&1; rc=$?
  echo "[$W] fill exit=$rc  $((SECONDS-t0))s" | tee -a "$L"
  [ $rc -ne 0 ] && exit 13
  if [ -z "$(model)" ]; then
    echo "[$W] fill exited 0 but wrote no .xml. Directory contains:" | tee -a "$L"
    ls -la "$D" | tee -a "$L"
    exit 13
  fi
fi
echo "[$W] COMPLETE $(date -Is) model=$(basename "$(model)")" | tee -a "$L"
INNER
chmod +x "$OUT/one_genome.sh"
awk -F'\t' -v g="$GS" -v o="$OUT" -v m="$FILLMED" -v t="$STAGE_TIMEOUT" -v k="$KTHREADS" \
  '{print $1"\t"$2"\t"$3"\t"g"\t"o"\t"m"\t"t"\t"k}' "$OUT/targets.tsv" \
  | xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -ra a <<< "{}"; '"$OUT"'/one_genome.sh "${a[@]}"'
note "all jobs returned"

say "D. what came out"
TOT=0; DONE=0
while IFS=$'\t' read -r w d f; do
  TOT=$((TOT+1))
  MX=$(find "$OUT/$w" -maxdepth 1 -name '*.xml' ! -name '*-draft.xml' -size +0 2>/dev/null | head -1)
  if [ -n "$MX" ]; then
    DONE=$((DONE+1)); printf '  [OK]   %-26s %s (%s)\n' "$w" "$(basename "$MX")" "$(du -h "$MX"|cut -f1)"
  else
    printf '  [FAIL] %-26s no model. last log lines:\n' "$w"
    tail -4 "$OUT/$w/gapseq.log" 2>/dev/null | sed 's/^/             /'
  fi
done < "$OUT/targets.tsv"
echo "           $DONE/$TOT models built"
[ "$DONE" -eq "$TOT" ] && ok "all genomes modelled" || no "re-run this script to resume the rest"
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
