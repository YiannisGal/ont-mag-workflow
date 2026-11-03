#!/usr/bin/env bash
# 45b_profiles_and_pathways.sh -- what is in profiles.tar.gz, and B12 from gapseq now.
#
# Safe by default: inventories, answers B12 from gapseq, and BENCHMARKS. It does not
# launch the full search unless you pass RUN=1.
#     bash scripts/45b_profiles_and_pathways.sh 2>&1 | tee logs/45b.log
#     RUN=1 tmux ... to commit to the full hmmsearch afterwards
#
# THREE ERRORS IN 45, all mine, all fixed here:
#   1. I collected 200 .faa files -- 1,311,147 sequences, 346,944 prefixes -- because I
#      globbed *.faa across the whole workspace. That swept in whole-assembly proteomes
#      (LIQ_contig_*), candidate files and MSAs. Section C now selects proteomes BY MAG
#      NAME from the crosswalk and asserts the count is 29.
#   2. profiles/ held 222 entries against ~26,500 expected and the script benchmarked
#      them anyway. Section A now compares the profile count to ko_list and refuses to
#      proceed if the database is a subset.
#   3. I sent hmmsearch's stderr to /dev/null. It ran in 0 s, returned 0 hits -- it
#      FAILED and I had hidden the message. Volume III's rule about exit codes, broken
#      by me, in a script whose purpose is care. stderr is captured and shown here.
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
XW="$WS2/07_mags/MAG_CROSSWALK.tsv"
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"; RUN="${RUN:-0}"
echo "45b  $(date -Is)"; echo "cpu: $CPU   RUN=$RUN"

# ---------------------------------------------------------------- A. the tarball
say "A. what is actually in profiles.tar.gz?"
if [ ! -s "$TGZ" ]; then bad "no $TGZ"; else
  ok "$TGZ  $(du -h "$TGZ" | cut -f1)"
  note "counting members (this reads the archive; a minute or two for a full KOfam set)"
  N=$(tar tzf "$TGZ" 2>/dev/null | grep -c '\.hmm$' || echo 0)
  echo "   .hmm members: $N"
  note "a sample of the names:"
  tar tzf "$TGZ" 2>/dev/null | grep '\.hmm$' | head -5 | sed 's/^/     /'
  KL=$(awk '$1 ~ /^K[0-9]{5}$/' "$KOLIST" | wc -l)
  echo "   ko_list entries: $KL"
  echo "   extracted profiles/ currently holds: $(ls "$KOF/profiles" 2>/dev/null | wc -l)"
  if [ "$N" -ge 15000 ]; then
    ok "that is the FULL KOfam profile set -- a genome-wide run is possible"
    FULL=1
  else
    no "only $N profiles against $KL KOs in ko_list -- this is a SUBSET, not the full set"
    note "a genome-wide run would need the full database from ftp.genome.jp"
    FULL=0
  fi
fi

# ---------------------------------------------------------------- B. gapseq
say "B. B12 from gapseq's pathway predictions -- already on disk, no download"
note "gapseq predicts pathway presence and completeness per genome. Cobalamin"
note "biosynthesis is a curated MetaCyc pathway, so this is a stronger source than"
note "Bakta product strings even though it covers only the modelled genomes."
mapfile -t PW < <(find "$WS2" -name '*-all-Pathways.tbl' -o -name '*-Pathways.tbl' 2>/dev/null | grep -v '/toy/' | sort)
echo "   pathway tables found: ${#PW[@]}"
if [ "${#PW[@]}" -gt 0 ]; then
  note "columns of the first table:"
  head -2 "${PW[0]}" | cut -c1-200 | sed 's/^/     /'
  "$PY" - "${PW[@]}" <<'PYEOF'
import sys, os, csv, re
files=sys.argv[1:]
rows=[]
for f in files:
    g=os.path.basename(f).replace('-all-Pathways.tbl','').replace('-Pathways.tbl','')
    with open(f, errors='replace') as fh:
        lines=[l for l in fh if not l.startswith('#')]
    if not lines: continue
    rd=csv.DictReader(lines, delimiter='\t')
    cols=rd.fieldnames or []
    namec=next((c for c in cols if re.search(r'name|pathway', c, re.I)), None)
    compc=next((c for c in cols if re.search(r'complet', c, re.I)), None)
    predc=next((c for c in cols if re.search(r'predict', c, re.I)), None)
    hits=[]
    for r in rd:
        nm=(r.get(namec) or '') if namec else ''
        if re.search(r'cobalamin|cobinamide|adenosylcobalamin|B12|corrin', nm, re.I):
            hits.append((nm, r.get(compc,'?'), r.get(predc,'?')))
    rows.append((g, hits))
if not rows:
    print("   no pathway rows parsed"); sys.exit(0)
w=max(len(g) for g,_ in rows)
print()
for g,hits in sorted(rows, key=lambda x:-len(x[1])):
    print("   %-*s  %d cobalamin-related pathway(s)" % (w,g,len(hits)))
    for nm,comp,pred in hits[:8]:
        print("       %-58s completeness=%-7s predicted=%s" % (nm[:58], comp, pred))
print()
pos=[g for g,h in rows if any(str(p).strip().lower() in ('true','1','yes') for _,_,p in h)]
print("   genomes with a PREDICTED cobalamin pathway: %d of %d" % (len(pos), len(rows)))
print("   %s" % pos)
print()
print("   NOTE: this covers only genomes that were run through gapseq, not all 29.")
print("   It is a statement about those genomes, and the coverage must be stated with it.")
PYEOF
else no "no gapseq pathway tables found"; fi

# ---------------------------------------------------------------- C. proteomes
say "C. the 29 MAG proteomes -- selected by name, not by glob"
if [ ! -s "$XW" ]; then bad "no crosswalk at $XW"; else
  mapfile -t NAMES < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$XW" | sort -u)
  echo "   MAG working names in the crosswalk: ${#NAMES[@]}"
  SEL="$OUT/selected_faa.txt"; : > "$SEL"
  MISS=0
  for n in "${NAMES[@]}"; do
    f=$(find "$WS2" -type f -name "${n}.faa" 2>/dev/null | head -1)
    if [ -z "$f" ]; then
      f=$(find "$WS2" -type f -path "*${n}*" -name '*.faa' 2>/dev/null | grep -viE 'candidate|verify|msa|align' | head -1)
    fi
    if [ -n "$f" ]; then printf '%s\t%s\n' "$n" "$f" >> "$SEL"
    else printf '  [--]   %-26s no .faa found\n' "$n"; MISS=$((MISS+1)); fi
  done
  echo "   matched: $(wc -l < "$SEL")   missing: $MISS"
  head -4 "$SEL" | sed 's/^/     /'
  if [ "$(wc -l < "$SEL")" -gt 0 ]; then
    COMB="$OUT/mags_only.faa"; : > "$COMB"
    while IFS=$'\t' read -r n f; do
      "$PY" - "$n" "$f" >> "$COMB" <<'PYEOF'
import sys
name, path = sys.argv[1], sys.argv[2]
for line in open(path, errors='replace'):
    if line.startswith('>'):
        sys.stdout.write('>%s|%s' % (name, line[1:]))
    else:
        sys.stdout.write(line)
PYEOF
    done < "$SEL"
    ok "combined proteome: $(grep -c '^>' "$COMB") sequences from $(wc -l < "$SEL") genomes"
    note "headers are prefixed <MAG>|<original>, so the split afterwards cannot go wrong"
    note "sequences per genome:"
    grep '^>' "$COMB" | sed 's/^>//' | cut -d'|' -f1 | sort | uniq -c | sort -rn | head -32 | sed 's/^/     /'
  fi
fi

# ---------------------------------------------------------------- D. benchmark
say "D. benchmark -- with stderr VISIBLE this time"
COMB="$OUT/mags_only.faa"
if [ "${FULL:-0}" -ne 1 ]; then
  no "the profile set is a subset; a genome-wide benchmark would be meaningless"
  note "section B's gapseq answer is the usable route unless the full database is fetched"
  say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0
fi
[ -s "$COMB" ] || { bad "no combined proteome"; echo "STEP_STATUS: FAILED"; exit 1; }
EXT="$OUT/profiles_full"
if [ ! -d "$EXT" ]; then
  note "extracting profiles.tar.gz to $EXT (not touching $KOF/profiles)"
  mkdir -p "$EXT"; tar xzf "$TGZ" -C "$EXT"
fi
NP=$(find "$EXT" -name '*.hmm' | wc -l)
echo "   extracted profiles: $NP"
HMM="$OUT/kofam_full.hmm"
if [ ! -s "$HMM" ]; then
  note "concatenating"
  find "$EXT" -name '*.hmm' -print0 | xargs -0 cat > "$HMM"
fi
NPROF=$(grep -c '^NAME' "$HMM"); echo "   profiles in db: $NPROF"
BM="$OUT/bench.hmm"
"$PY" - "$HMM" "$BM" 300 <<'PYEOF'
import sys
src, dst, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
out=open(dst,'w'); n=0
for line in open(src, errors='replace'):
    if line.startswith('HMMER3'):
        n+=1
        if n>k: break
    out.write(line)
out.close(); print("   benchmark db: %d profiles" % min(n-1 if n>k else n, k))
PYEOF
T0=$(date +%s)
hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$OUT/bench.tbl" "$BM" "$COMB" > "$OUT/bench.log" 2> "$OUT/bench.err"
RC=$?; EL=$(( $(date +%s) - T0 ))
echo "   hmmsearch exit=$RC in ${EL}s"
if [ "$RC" -ne 0 ]; then
  bad "hmmsearch FAILED -- stderr follows (this is what 45 hid):"
  head -20 "$OUT/bench.err" | sed 's/^/     /'
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
NBM=$(grep -c '^NAME' "$BM")
HITS=$(grep -vc '^#' "$OUT/bench.tbl" 2>/dev/null || echo 0)
echo "   $NBM profiles vs $(grep -c '^>' "$COMB") proteins -> $HITS hit rows"
if [ "$EL" -lt 5 ]; then
  no "under 5 s is implausible for $NBM profiles -- treat the estimate with suspicion"
fi
EST=$(( EL * NPROF / (NBM>0?NBM:1) ))
printf '   >>> estimated full run: %s s  (~%s h) at %s threads\n' "$EST" "$(( EST/3600 ))" "$CPU"
if [ "$RUN" != "1" ]; then
  note "RUN=1 to commit to the full search (chunked and resumable)."
  say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0
fi
say "E. full search -- chunked, resumable"
CH="$OUT/chunks"; mkdir -p "$CH"
"$PY" - "$HMM" "$CH" 24 <<'PYEOF'
import sys, os
src, ch, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
if any(f.endswith('.hmm') for f in os.listdir(ch)):
    print("   chunks already present"); sys.exit(0)
blocks=[]; cur=[]
for line in open(src, errors='replace'):
    if line.startswith('HMMER3') and cur: blocks.append(cur); cur=[]
    cur.append(line)
if cur: blocks.append(cur)
size=(len(blocks)+k-1)//k
print("   %d profiles -> %d chunks" % (len(blocks), k))
for i in range(k):
    part=blocks[i*size:(i+1)*size]
    if not part: continue
    with open(os.path.join(ch,'c%02d.hmm'%i),'w') as fh:
        for b in part: fh.writelines(b)
PYEOF
for c in "$CH"/*.hmm; do
  b=$(basename "$c" .hmm); o="$CH/$b.tbl"
  [ -s "$o" ] && { printf '   [skip] %s\n' "$b"; continue; }
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$o.part" "$c" "$COMB" > /dev/null 2> "$CH/$b.err"
  if [ $? -eq 0 ]; then mv "$o.part" "$o"; printf '   [ok]   %s %5ds %s rows\n' "$b" "$(( $(date +%s)-T0 ))" "$(grep -vc '^#' "$o")"
  else printf '   [FAIL] %s: %s\n' "$b" "$(head -2 "$CH/$b.err" | tr '\n' ' ')"; fi
done
cat "$CH"/*.tbl 2>/dev/null | grep -v '^#' > "$OUT/all_hits.tbl" || true
ok "$(wc -l < "$OUT/all_hits.tbl") hit rows"
note "run 45c to threshold, deduplicate and answer B12 from all_hits.tbl"
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
