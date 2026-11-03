#!/usr/bin/env bash
# 45_kofam_genomewide.sh -- find an existing genome-wide annotation, or make one.
#
# RUN IN tmux. Chunked and RESUMABLE: every chunk skips if its output exists, so a
# dropped connection costs nothing.
#     tmux new -s 45
#     bash scripts/45_kofam_genomewide.sh 2>&1 | tee logs/45.log
#
#   BENCH_ONLY=1   stop after the timing estimate, before committing to the full run
#   CPU=16         threads for hmmsearch (default: all cores)
#
# WHY
#   44c established that 11_kofam/hits/*.tbl is a TARGETED marker search: median 13
#   distinct KOs per genome, the three highest all archaea (45, 44, 42), which is the
#   signature of the methanogenesis marker set that 41f used for mcr 8/8.
#   Step 21's "15 B12 producers, 0 consumers" came from Bakta product strings and has
#   never been re-derived at curated thresholds. It cannot be, from that input.
#
# ORDER, as asked: look for an existing genome-wide annotation FIRST. Hours of compute
# to recreate something already on disk would be the worst outcome here.
#
# DESIGN NOTE. The obvious approach -- hmmsearch once per genome -- means 29 passes over
# the same 26,498 profiles. One pass over a COMBINED proteome is far cheaper, because
# hmmsearch's cost is dominated by walking the profile database. Bakta protein IDs
# already carry a genome prefix (ETOLIQMAG003_04260), so the result splits cleanly
# afterwards. Section D benchmarks on a profile subset and extrapolates before section E
# commits to anything.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/45_kofam_full"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp45"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
[ -d "$ENVDIR" ] && export PATH="$ENVDIR/bin:$PATH"
KOF="$WS2/11_kofam"; KOLIST="$KOF/ko_list"
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"
BENCH_ONLY="${BENCH_ONLY:-0}"
echo "45  $(date -Is)"; echo "out: $OUT"; echo "cpu: $CPU"

# ---------------------------------------------------------------- A. what exists
say "A. is a genome-wide functional annotation already on disk?"
note "any table with hundreds of KOs per genome would make section E unnecessary."
mapfile -t CAND < <(find "$WS2" -maxdepth 4 -type f \
  \( -name '*.tbl' -o -name '*kofam*' -o -name '*KEGG*' -o -name '*kegg*' \
     -o -name '*emapper*' -o -name '*eggnog*' -o -name '*annotations*' \) \
  -size +200k 2>/dev/null | grep -viE 'ko_list|profiles/' | head -40)
if [ "${#CAND[@]}" -eq 0 ]; then
  no "no large annotation tables found"
else
  printf '  %-96s %10s\n' "file" "size"
  for f in "${CAND[@]}"; do printf '  %-96s %10s\n' "${f#$WS2/}" "$(du -h "$f" | cut -f1)"; done
fi
note "and how many DISTINCT KO identifiers each contains:"
BEST=""; BESTN=0
for f in "${CAND[@]:-}"; do
  [ -s "$f" ] || continue
  n=$(grep -oE '\bK[0-9]{5}\b' "$f" 2>/dev/null | sort -u | wc -l)
  [ "$n" -gt 0 ] && printf '     %-90s %6d KOs\n' "$(basename "$f")" "$n"
  if [ "$n" -gt "$BESTN" ]; then BESTN=$n; BEST="$f"; fi
done
echo
if [ "$BESTN" -ge 1500 ]; then
  ok "richest source: $BEST ($BESTN distinct KOs) -- that looks genome-wide"
  note "section E will be SKIPPED and B12 answered from this file"
  EXISTING="$BEST"
else
  no "richest source has only $BESTN distinct KOs across everything -- not genome-wide"
  note "a genome-wide annotation of 29 MAGs should contain several thousand distinct KOs"
  EXISTING=""
fi

# ---------------------------------------------------------------- B. inputs
say "B. the inputs a full run needs"
mapfile -t FAA < <(find "$WS2" -type f -name '*.faa' 2>/dev/null | grep -viE 'test|tmp' | head -200)
echo "   proteome (.faa) files found: ${#FAA[@]}"
for f in "${FAA[@]:0:4}"; do printf '     %-80s %6s seqs\n' "${f#$WS2/}" "$(grep -c '^>' "$f")"; done
[ "${#FAA[@]}" -gt 4 ] && echo "     ... and $(( ${#FAA[@]} - 4 )) more"
PROF=$(find "$KOF" -maxdepth 2 -type d -name 'profiles' 2>/dev/null | head -1)
[ -n "$PROF" ] && ok "profiles: $PROF ($(ls "$PROF" | wc -l) entries)" || no "no profiles/ under $KOF"
[ -s "$KOLIST" ] && ok "ko_list: $KOLIST" || no "no ko_list"
for t in hmmsearch hmmpress exec_annotation; do
  command -v "$t" >/dev/null 2>&1 && printf '  [OK]   %-16s %s\n' "$t" "$(command -v "$t")" \
                                  || printf '  [--]   %-16s not on PATH\n' "$t"
done

if [ -n "$EXISTING" ]; then
  say "C. B12 from the existing annotation"
  "$PY" - "$KOLIST" "$EXISTING" "$OUT" <<'PYEOF'
import sys, re, os, json, collections
kolist, src, out = sys.argv[1:4]
defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p)>=3 and p[0].startswith('K'): defn[p[0]]=p[-1]
def kos(pat):
    rx=re.compile(pat, re.I); return {k for k,d in defn.items() if rx.search(d or '')}
PROD=kos(r'cob\(|cobalamin|cobyrinate|cobyrinic|cobinamide|precorrin|cobaltochelatase|adenosylcobinamide|cobalt-precorrin')
CONS=kos(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
print("   biosynthesis %d  B12-dependent %d KOs" % (len(PROD),len(CONS)))
per=collections.defaultdict(set)
for line in open(src, errors='replace'):
    ko=re.search(r'\bK\d{5}\b', line)
    if not ko: continue
    g=re.match(r'\s*(\S+)', line)
    gen=g.group(1) if g else '?'
    m=re.match(r'([A-Za-z0-9]+MAG\d+)', gen)
    key=m.group(1) if m else gen.split('_')[0]
    per[key].add(ko.group(0))
print("   genomes/prefixes seen: %d" % len(per))
for k,v in sorted(per.items(), key=lambda x:-len(x[1]))[:35]:
    print("     %-24s %5d KOs  biosynth=%-3d B12dep=%d" % (k, len(v), len(PROD&v), len(CONS&v)))
json.dump({k:sorted(v) for k,v in per.items()}, open(os.path.join(out,'kos_from_existing.json'),'w'))
PYEOF
  say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0
fi

# ---------------------------------------------------------------- C. combine
say "C. build one combined proteome"
COMB="$OUT/all_proteins.faa"
if [ -s "$COMB" ]; then ok "already built: $(grep -c '^>' "$COMB") sequences"; else
  : > "$COMB"
  for f in "${FAA[@]}"; do cat "$f" >> "$COMB"; done
  ok "$(grep -c '^>' "$COMB") sequences from ${#FAA[@]} file(s)"
fi
note "genome prefixes present (first field of each header, trimmed to the MAG tag):"
grep '^>' "$COMB" | sed 's/^>//' | awk '{print $1}' | sed -E 's/_[0-9]+$//' | sort | uniq -c | sort -rn | head -32 | sed 's/^/     /'
NPREF=$(grep '^>' "$COMB" | awk '{print $1}' | sed -E 's/_[0-9]+$//' | sort -u | wc -l)
echo "   distinct prefixes: $NPREF"
[ "$NPREF" -lt 20 ] && note "fewer prefixes than genomes -- the split after the run will need the crosswalk"

# ---------------------------------------------------------------- D. benchmark
say "D. BENCHMARK before committing -- 500 profiles, then extrapolate"
HMM="$OUT/kofam_all.hmm"
if [ -s "$HMM" ]; then ok "concatenated profile db already present ($(du -h "$HMM" | cut -f1))"; else
  [ -n "$PROF" ] || { bad "no profiles directory; cannot build the db"; echo "STEP_STATUS: FAILED"; exit 1; }
  note "concatenating profiles (this takes a few minutes)"
  find "$PROF" -name '*.hmm' -print0 | xargs -0 cat > "$HMM"
  ok "$(grep -c '^NAME' "$HMM") profiles, $(du -h "$HMM" | cut -f1)"
fi
NPROF=$(grep -c '^NAME' "$HMM")
BM="$OUT/bench.hmm"
if [ ! -s "$BM" ]; then
  awk '/^HMMER3/{n++} n<=500' "$HMM" > "$BM"
fi
NBM=$(grep -c '^NAME' "$BM" || echo 0)
echo "   benchmark profiles: $NBM of $NPROF"
if [ "$NBM" -gt 0 ]; then
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$OUT/bench.tbl" "$BM" "$COMB" > /dev/null 2>&1
  T1=$(date +%s); EL=$((T1-T0))
  echo "   $NBM profiles vs $(grep -c '^>' "$COMB") proteins took ${EL}s with $CPU threads"
  EST=$(( EL * NPROF / (NBM>0?NBM:1) ))
  printf '   >>> estimated full run: %d s = %.1f h\n' "$EST" "$(echo "$EST/3600" | bc -l)"
  echo "   >>> benchmark hits: $(grep -vc '^#' "$OUT/bench.tbl" 2>/dev/null || echo 0)"
else
  no "could not build a benchmark subset; skipping the estimate"
fi
if [ "$BENCH_ONLY" = "1" ]; then
  note "BENCH_ONLY=1 -- stopping here. Re-run without it to commit to the full search."
  say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0
fi

# ---------------------------------------------------------------- E. full run
say "E. the full search, in chunks so it is resumable"
CH="$OUT/chunks"; mkdir -p "$CH"
if [ -z "$(ls -A "$CH" 2>/dev/null | grep '\.hmm$')" ]; then
  note "splitting the profile db into 20 chunks"
  "$PY" - "$HMM" "$CH" <<'PYEOF'
import sys, os
hmm, ch = sys.argv[1], sys.argv[2]
blocks=[]; cur=[]
for line in open(hmm, errors='replace'):
    if line.startswith('HMMER3') and cur:
        blocks.append(cur); cur=[]
    cur.append(line)
if cur: blocks.append(cur)
n=len(blocks); k=20; size=(n+k-1)//k
print("   %d profiles -> %d chunks of <=%d" % (n,k,size))
for i in range(k):
    part=blocks[i*size:(i+1)*size]
    if not part: continue
    with open(os.path.join(ch,'c%02d.hmm'%i),'w') as fh:
        for b in part: fh.writelines(b)
PYEOF
fi
ls "$CH"/*.hmm 2>/dev/null | wc -l | sed 's/^/   chunks: /'
for c in "$CH"/*.hmm; do
  b=$(basename "$c" .hmm); o="$CH/$b.tbl"
  if [ -s "$o" ]; then printf '   [skip] %s (%s rows)\n' "$b" "$(grep -vc '^#' "$o")"; continue; fi
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$o.part" "$c" "$COMB" > /dev/null 2>&1
  RC=$?
  if [ "$RC" -eq 0 ]; then mv "$o.part" "$o"; printf '   [ok]   %s  %4ds  %s rows\n' "$b" "$(( $(date +%s)-T0 ))" "$(grep -vc '^#' "$o")"
  else rm -f "$o.part"; printf '   [FAIL] %s rc=%s\n' "$b" "$RC"; fi
done
cat "$CH"/*.tbl 2>/dev/null | grep -v '^#' > "$OUT/all_hits.tbl" || true
ok "combined: $(wc -l < "$OUT/all_hits.tbl") hit rows"

# ---------------------------------------------------------------- F. B12
say "F. B12, genome-wide, thresholded and deduplicated"
"$PY" - "$KOLIST" "$OUT/all_hits.tbl" "$OUT" <<'PYEOF'
import sys, os, re, json, collections
kolist, hits, out = sys.argv[1:4]
thr={}; defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p)<3 or not p[0].startswith('K'): continue
        try: t=float(p[1])
        except ValueError: t=None
        thr[p[0]]=t; defn[p[0]]=p[-1]
def kos(pat):
    rx=re.compile(pat, re.I); return {k for k,d in defn.items() if rx.search(d or '')}
PROD=kos(r'cob\(|cobalamin|cobyrinate|cobyrinic|cobinamide|precorrin|cobaltochelatase|adenosylcobinamide|cobalt-precorrin')
CONS=kos(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
TRAN=kos(r'\bbtu[BFCD]\b|vitamin B12 transport|cobalamin transport')
best={}
for line in open(hits, errors='replace'):
    p=line.split()
    if len(p) < 6: continue
    gene, ko = p[0], p[2]
    if not re.fullmatch(r'K\d{5}', ko): continue
    try: score=float(p[5])
    except ValueError: continue
    t=thr.get(ko)
    if t is not None and score < t: continue          # ko_list's own curated threshold
    if gene not in best or score > best[gene][1]: best[gene]=(ko,score)
per=collections.defaultdict(set)
for gene,(ko,sc) in best.items():
    m=re.match(r'([A-Za-z0-9]+?MAG\d+)', gene)
    per[m.group(1) if m else gene.rsplit('_',1)[0]].add(ko)
print("   genes with a thresholded hit: %d across %d genomes" % (len(best), len(per)))
w=max((len(k) for k in per), default=10)
print()
print("   %-*s %7s %9s %9s %9s" % (w,'genome','KOs','biosynth','B12-dep','transport'))
rows=[]
for g,ks in sorted(per.items(), key=lambda x:-len(PROD & x[1])):
    rows.append(dict(genome=g, nko=len(ks), prod=sorted(PROD&ks), cons=sorted(CONS&ks), tran=sorted(TRAN&ks)))
    print("   %-*s %7d %9d %9d %9d" % (w,g,len(ks),len(PROD&ks),len(CONS&ks),len(TRAN&ks)))
print()
for cut in (5,10,15,20):
    print("   producers at >=%2d biosynthesis KOs : %d" % (cut, sum(1 for r in rows if len(r['prod'])>=cut)))
print("   genomes with a B12-DEPENDENT enzyme : %d" % sum(1 for r in rows if r['cons']))
print("   genomes with a B12 transporter      : %d" % sum(1 for r in rows if r['tran']))
print()
print("   step 21 said 15 producers, 0 consumers, from Bakta product strings.")
print("   The cut is a CHOICE and is printed at four values above; state the one you use.")
json.dump(rows, open(os.path.join(out,'b12_genomewide.json'),'w'), indent=2)
med=sorted(r['nko'] for r in rows)[len(rows)//2] if rows else 0
print()
print("   median KOs per genome: %d" % med)
if med < 200:
    print("   >>> STILL not genome-wide. Do not report these numbers.")
else:
    print("   >>> genome-wide. These numbers are reportable, with the cut stated.")
PYEOF
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
