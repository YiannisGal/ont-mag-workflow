#!/usr/bin/env bash
# 45c_run_kofam.sh -- the genome-wide KOfam run, and B12 for all 29.
#
# RUN IN tmux. Chunked and resumable; a dropped connection costs one chunk at most.
#     tmux new -s 45c
#     bash scripts/45c_run_kofam.sh 2>&1 | tee logs/45c.log
#   BENCH_ONLY=1  stop after the timing estimate
#   CPU=28        threads (default: all cores)
#
# 45b ESTABLISHED
#   - profiles.tar.gz holds 27,864 .hmm profiles, 1.5 GB. That IS the full KOfam set.
#     My check misread it: with pipefail on, tar exits non-zero on the pipe, `|| echo 0`
#     fired, and the count became "27864\n0" which fails an integer test. Fixed below by
#     turning pipefail off for that one measurement.
#   - the 29 MAG proteomes resolve cleanly from 06_binqc/checkm2/protein_files/,
#     90,905 sequences, headers rewritten as >MAG|original.
#   - gapseq's pathway tables say ONE genome of thirteen is predicted to synthesise
#     cobalamin: LIQ_das_vamb_003 = MTP4, adenosylcobalamin biosynthesis at 66%
#     completeness, predicted=true. Everything else is predicted=false at 0-40%.
#     Corrinoids are central to methanogenesis and MTP4's methylamine methyltransferases
#     all need a corrinoid cofactor, so that is coherent -- and it makes step 21's
#     "15 producers" from Bakta product strings a ~15-fold over-count.
#
# This run tests that across all 29 genomes with curated per-KO thresholds, which is the
# only way step 21 can be closed rather than replaced by another single-source claim.
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
KOF="$WS2/11_kofam"; KOLIST="$KOF/ko_list"; TGZ="$KOF/profiles.tar.gz"
COMB="$OUT/mags_only.faa"; EXT="$OUT/profiles_full"; CH="$OUT/chunks"
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"; BENCH_ONLY="${BENCH_ONLY:-0}"; NCH="${NCH:-24}"
echo "45c  $(date -Is)"; echo "cpu: $CPU  chunks: $NCH"
[ -s "$COMB" ] || { bad "no $COMB -- run 45b first"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "proteome: $(grep -c '^>' "$COMB") sequences, $(grep '^>' "$COMB" | cut -d'|' -f1 | sort -u | wc -l) genomes"

# ---------------------------------------------------------------- A. disk
say "A. disk, before extracting 1.5 GB of profiles"
df -h ${SCRATCH_ROOT} | sed 's/^/     /'
AVAIL=$(df -Pk ${SCRATCH_ROOT} | awk 'NR==2{print $4}')
echo "   available: $(( AVAIL / 1024 / 1024 )) GB"
[ "$AVAIL" -lt 10485760 ] && { bad "under 10 GB free; extraction needs room"; exit 1; }
ok "enough room"

# ---------------------------------------------------------------- B. extract
say "B. extract and chunk the profile database"
set +o pipefail
NMEM=$(tar tzf "$TGZ" 2>/dev/null | grep -c '\.hmm$')
set -o pipefail
echo "   members in the archive: $NMEM"
[ "$NMEM" -lt 15000 ] && { bad "not the full set after all"; exit 1; }
if [ -d "$EXT" ] && [ "$(find "$EXT" -name '*.hmm' | wc -l)" -ge 15000 ]; then
  ok "already extracted: $(find "$EXT" -name '*.hmm' | wc -l) profiles"
else
  mkdir -p "$EXT"; note "extracting (a few minutes)"
  tar xzf "$TGZ" -C "$EXT"
  ok "extracted $(find "$EXT" -name '*.hmm' | wc -l) profiles"
fi
mkdir -p "$CH"
if [ "$(ls "$CH"/*.hmm 2>/dev/null | wc -l)" -lt "$NCH" ]; then
  note "building $NCH chunks by concatenating groups of profiles (no giant intermediate)"
  find "$EXT" -name '*.hmm' | sort > "$OUT/profile_list.txt"
  TOT=$(wc -l < "$OUT/profile_list.txt"); PER=$(( (TOT + NCH - 1) / NCH ))
  echo "   $TOT profiles, $PER per chunk"
  i=0
  while [ "$i" -lt "$NCH" ]; do
    C=$(printf '%s/c%02d.hmm' "$CH" "$i")
    if [ ! -s "$C" ]; then
      sed -n "$(( i*PER + 1 )),$(( (i+1)*PER ))p" "$OUT/profile_list.txt" | tr '\n' '\0' | xargs -0 cat > "$C" 2>/dev/null || true
    fi
    i=$(( i + 1 ))
  done
fi
ls "$CH"/c*.hmm 2>/dev/null | wc -l | sed 's/^/   chunks ready: /'
for c in "$CH"/c*.hmm; do printf '     %s  %s profiles\n' "$(basename "$c")" "$(grep -c '^NAME' "$c")"; done | head -4

# ---------------------------------------------------------------- C. benchmark
say "C. benchmark on one chunk -- stderr visible"
C0=$(ls "$CH"/c*.hmm | head -1); N0=$(grep -c '^NAME' "$C0")
if [ -s "$CH/$(basename "$C0" .hmm).tbl" ]; then
  ok "first chunk already done; skipping the benchmark"
  EST_NOTE="(first chunk already complete)"
else
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$CH/bench.tbl" "$C0" "$COMB" > /dev/null 2> "$OUT/bench.err"
  RC=$?; EL=$(( $(date +%s) - T0 ))
  echo "   exit=$RC  elapsed=${EL}s  profiles=$N0"
  if [ "$RC" -ne 0 ]; then
    bad "hmmsearch failed; stderr:"; head -20 "$OUT/bench.err" | sed 's/^/     /'
    echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
  fi
  echo "   hit rows: $(grep -vc '^#' "$CH/bench.tbl")"
  printf '   >>> estimated total for %s chunks: %s min (~%s h)\n' "$NCH" "$(( EL*NCH/60 ))" "$(( EL*NCH/3600 ))"
  rm -f "$CH/bench.tbl"
fi
[ "$BENCH_ONLY" = "1" ] && { note "BENCH_ONLY=1"; say DONE; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0; }

# ---------------------------------------------------------------- D. run
say "D. the full search"
note "started $(date -Is)"
for c in "$CH"/c*.hmm; do
  b=$(basename "$c" .hmm); o="$CH/$b.tbl"
  [ -s "$o" ] && { printf '   [skip] %s  %s rows\n' "$b" "$(grep -vc '^#' "$o")"; continue; }
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$o.part" "$c" "$COMB" > /dev/null 2> "$CH/$b.err"
  if [ $? -eq 0 ]; then
    mv "$o.part" "$o"
    printf '   [ok]   %s %5ds  %s rows   (%s)\n' "$b" "$(( $(date +%s)-T0 ))" "$(grep -vc '^#' "$o")" "$(date +%H:%M)"
  else
    rm -f "$o.part"; printf '   [FAIL] %s : %s\n' "$b" "$(head -2 "$CH/$b.err" | tr '\n' ' ')"
  fi
done
DONE_N=$(ls "$CH"/c*.tbl 2>/dev/null | wc -l)
echo "   chunks complete: $DONE_N of $NCH"
[ "$DONE_N" -lt "$NCH" ] && no "incomplete -- re-run this script to finish the rest"
cat "$CH"/c*.tbl 2>/dev/null | grep -v '^#' > "$OUT/all_hits.tbl" || true
ok "$(wc -l < "$OUT/all_hits.tbl") hit rows total"

# ---------------------------------------------------------------- E. B12
say "E. B12 across all 29, thresholded and deduplicated"
"$PY" - "$KOLIST" "$OUT/all_hits.tbl" "$OUT" <<'PYEOF'
import sys, os, re, json, collections, statistics
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
print("   KO sets: biosynthesis %d, B12-dependent %d, transport %d" % (len(PROD),len(CONS),len(TRAN)))
best={}
for line in open(hits, errors='replace'):
    p=line.split()
    if len(p) < 6: continue
    gene, ko = p[0], p[2]
    if not re.fullmatch(r'K\d{5}', ko): continue
    try: score=float(p[5])          # hmmsearch --tblout: 4=E-value, 5=score
    except ValueError: continue
    t=thr.get(ko)
    if t is not None and score < t: continue
    if gene not in best or score > best[gene][1]: best[gene]=(ko, score)
per=collections.defaultdict(set)
for gene,(ko,sc) in best.items():
    per[gene.split('|',1)[0]].add(ko)
print("   genes above threshold: %d across %d genomes" % (len(best), len(per)))
med=statistics.median(len(v) for v in per.values()) if per else 0
print("   median distinct KOs per genome: %.0f" % med)
if med < 200:
    print("   >>> STILL under 200 KOs per genome. Something is wrong; do not report these.")
else:
    print("   >>> genome-wide. Reportable.")
rows=[]
w=max((len(k) for k in per), default=10)
print()
print("   %-*s %7s %9s %9s %9s" % (w,'genome','KOs','biosynth','B12-dep','transport'))
for g,ks in sorted(per.items(), key=lambda x:-len(PROD & x[1])):
    r=dict(genome=g, nko=len(ks), prod=sorted(PROD&ks), cons=sorted(CONS&ks), tran=sorted(TRAN&ks))
    rows.append(r)
    print("   %-*s %7d %9d %9d %9d" % (w,g,len(ks),len(r['prod']),len(r['cons']),len(r['tran'])))
print()
for cut in (5,10,15,20,25):
    print("   producers at >=%2d biosynthesis KOs : %d" % (cut, sum(1 for r in rows if len(r['prod'])>=cut)))
print("   genomes with a B12-DEPENDENT enzyme : %d" % sum(1 for r in rows if r['cons']))
print("   genomes with a B12 transporter      : %d" % sum(1 for r in rows if r['tran']))
print()
print("   step 21 (Bakta product strings) said: 15 producers, 0 consumers.")
print("   gapseq pathways (13 genomes) said   : 1 predicted producer, LIQ_das_vamb_003 = MTP4.")
print("   The cut is a CHOICE. Five values are printed; state the one the thesis uses.")
top=sorted(rows, key=lambda r:-len(r['prod']))[:5]
print()
print("   the five with the most biosynthesis KOs:")
for r in top:
    print("     %-26s %d KOs: %s" % (r['genome'], len(r['prod']), ', '.join(r['prod'][:8])))
json.dump(rows, open(os.path.join(out,'b12_genomewide.json'),'w'), indent=2)
print()
print("   full per-genome KO sets written to %s" % os.path.join(out,'b12_genomewide.json'))
PYEOF
say "F. the by-product"
note "all_hits.tbl is now a genome-wide KO annotation for all 29 MAGs. That is the"
note "table the metabolic-potential figures should be built from, rather than Bakta"
note "product strings -- and it is what step 16b wanted."
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
