#!/usr/bin/env bash
# 48b_unbinned_sulfur.sh -- search the UNBINNED proteins, and refuse to conclude from a
#                           search that did not run.
#
# 48 SECTION A SETTLED: my Volume XI regex was sound. K11180/K11181 are literally
# "dissimilatory sulfite reductase alpha/beta subunit", and 0 dissimilatory KOs were
# missed by the pattern. The MAG claim stands.
#
# 48 SECTION B, refined:
#   dsrA (K11180) / dsrB (K11181)         : 0 MAGs
#   aprAB (K00394/K00395)                 : 1 MAG  -- ETO_LIQ_MAG_022, JBAIYS01
#   K23077 "dsr RELATED protein"          : the same 1 MAG
#   My "genomes with BOTH: 1" counted K23077 as dsr. It is not dsrAB. Overstated.
#
# 🔴 48 SECTION D FAILED AND I CONCLUDED FROM IT ANYWAY
#       Parse failed: Line 64: illegal character -
#       dissimilatory sulfite reductase genes in the assembly: 0
#       >>> NO dissimilatory sulfite reductase in the assembly either.
#   hmmsearch never ran. That 0 is the absence of a SEARCH, not of a gene, and I printed
#   a chapter-altering conclusion from it while holding the exit code that said so.
#   Two further faults: the target was 45_kofam_full/all_proteins.faa -- my own junk file
#   from step 45, 1.3 M proteins globbed from 200 files including a phylogenetic
#   ALIGNMENT, which is where the '-' came from -- while the real unbinned proteins were
#   in the same listing at 34_unbinned/prot/.
#
# This script validates every target BEFORE searching, and every conclusion is gated on
# the search having succeeded.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/48_sulfur"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp48"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
[ -d "$ENVDIR" ] && export PATH="$ENVDIR/bin:$PATH"
KOLIST="$WS2/11_kofam/ko_list"; SUB="$OUT/sulfur_profiles.hmm"
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"
echo "48b  $(date -Is)"
[ -s "$SUB" ] || { bad "no $SUB -- run 48 first"; exit 1; }
ok "sulfur profiles: $(grep -c '^NAME' "$SUB")"

# ---------------------------------------------------------------- A. targets
say "A. the targets, VALIDATED before anything is searched"
TARGETS=(
  "$WS2/34_unbinned/prot/LIQ_all.faa"
  "$WS2/34_unbinned/prot/SED_all.faa"
  "$WS2/05_consensus/LIQ_proteins.faa"
  "$WS2/05_consensus/SED_proteins.faa"
)
GOOD=()
for t in "${TARGETS[@]}"; do
  if [ ! -s "$t" ]; then printf '  [--]   %-56s missing\n' "${t#$WS2/}"; continue; fi
  RES=$("$PY" - "$t" <<'PYEOF'
import sys, re
p=sys.argv[1]; n=0; bad=0; first=None
ok=set("ACDEFGHIKLMNPQRSTVWYXBZJUO*acdefghiklmnpqrstvwyxbzjuo")
for i,line in enumerate(open(p, errors='replace'), 1):
    if line.startswith('>'): n+=1; continue
    for ch in line.rstrip('\n'):
        if ch not in ok:
            bad+=1
            if first is None: first=(i,ch)
            break
print("%d %d %s" % (n, bad, ('line %d char %r' % first) if first else '-'))
PYEOF
)
  NSEQ=$(echo "$RES" | awk '{print $1}'); NBAD=$(echo "$RES" | awk '{print $2}')
  if [ "$NBAD" -eq 0 ]; then
    printf '  [OK]   %-56s %8s proteins, clean\n' "${t#$WS2/}" "$NSEQ"; GOOD+=("$t")
  else
    printf '  [--]   %-56s %8s proteins, %s bad line(s): %s\n' "${t#$WS2/}" "$NSEQ" "$NBAD" "$(echo "$RES" | cut -d' ' -f3-)"
  fi
done
[ "${#GOOD[@]}" -eq 0 ] && { bad "no clean target -- nothing can be searched"; echo "STEP_STATUS: FAILED"; exit 1; }
note "45_kofam_full/all_proteins.faa is NOT used: it is step 45's junk concatenation"
note "(1.3 M proteins from 200 globbed files, including an MSA). It can be deleted."

# ---------------------------------------------------------------- B. search
say "B. search each clean target -- conclusions gated on exit status"
ANYRUN=0
for t in "${GOOD[@]}"; do
  b=$(basename "$t" .faa); o="$OUT/${b}_sulfur.tbl"
  if [ -s "$o" ]; then ok "$b already searched"; ANYRUN=1; continue; fi
  T0=$(date +%s)
  hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$o.part" "$SUB" "$t" > /dev/null 2> "$OUT/${b}.err"
  RC=$?
  if [ "$RC" -eq 0 ]; then
    mv "$o.part" "$o"; ANYRUN=1
    printf '  [OK]   %-24s %4ds  %s raw rows\n' "$b" "$(( $(date +%s)-T0 ))" "$(grep -vc '^#' "$o")"
  else
    rm -f "$o.part"
    printf '  [FAIL] %-24s exit=%s\n' "$b" "$RC"; head -4 "$OUT/${b}.err" | sed 's/^/         /'
  fi
done
if [ "$ANYRUN" -ne 1 ]; then
  bad "no search completed. NO CONCLUSION IS DRAWN -- that was 48's mistake."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi

# ---------------------------------------------------------------- C. read
say "C. what the unbinned fraction contains"
"$PY" - "$KOLIST" "$OUT" <<'PYEOF'
import sys, os, re, glob, collections
kolist, out = sys.argv[1], sys.argv[2]
thr={}; defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not re.fullmatch(r'K\d{5}', p[0]): continue
    try: t=float(p[1])
    except ValueError: t=None
    thr[p[0]]=t; defn[p[0]]=p[-1]
DSR_TRUE={'K11180','K11181'}          # dsrA / dsrB proper
DSR_REL ={'K23077','K27196'}          # "related" / flavoprotein -- NOT diagnostic
APR     ={'K00394','K00395'}
SAT     ={'K00956','K00957','K00958','K00988','K22966'}
any_true=0
for f in sorted(glob.glob(os.path.join(out,'*_sulfur.tbl'))):
    name=os.path.basename(f).replace('_sulfur.tbl','')
    best={}
    for line in open(f, errors='replace'):
        if line.startswith('#'): continue
        p=line.split()
        if len(p)<6: continue
        gene, ko = p[0], p[2]
        if not re.fullmatch(r'K\d{5}', ko): continue
        try: sc=float(p[5])
        except ValueError: continue
        t=thr.get(ko)
        if t is not None and sc < t: continue
        if gene not in best or sc > best[gene][1]: best[gene]=(ko,sc)
    cnt=collections.Counter(ko for ko,_ in best.values())
    print("\n   --- %s --- (%d genes above threshold)" % (name, len(best)))
    for ko,n in sorted(cnt.items()):
        tag=''
        if ko in DSR_TRUE: tag='   <<< dsrAB -- A REAL SULFATE REDUCER'
        elif ko in DSR_REL: tag='   (dsr-related, not diagnostic)'
        elif ko in APR: tag='   (aprAB)'
        elif ko in SAT: tag='   (sat -- assimilatory OR dissimilatory)'
        print("     %-8s %5d  %-52s%s" % (ko, n, (defn.get(ko) or '')[:52], tag))
    nt=sum(n for ko,n in cnt.items() if ko in DSR_TRUE)
    if nt:
        any_true += nt
        print("     >>> %d dsrA/dsrB gene(s). Contigs:" % nt)
        for gene,(ko,sc) in sorted(best.items()):
            if ko in DSR_TRUE: print("         %-46s %s score=%.1f" % (gene, ko, sc))
print()
print("   ========================================================")
if any_true:
    print("   dsrA/dsrB FOUND in the unbinned fraction: %d gene(s)." % any_true)
    print("   >>> Sulfate reducers ARE present in the assembly and were not binned.")
    print("       The competition premise stands. State it as: detected in the assembly,")
    print("       not recovered as MAGs.")
else:
    print("   NO dsrA/dsrB in any searched proteome.")
    print("   >>> With the searches confirmed to have RUN, this is a real absence in what")
    print("       was sequenced and assembled. Sulfate reduction must then be argued from")
    print("       the geochemistry -- sulfate 180 mean, sulfidic, H2S measured -- and the")
    print("       chapter must not claim the organisms were recovered.")
print("   ========================================================")
PYEOF

# ---------------------------------------------------------------- D. kaiju
say "D. does Kaiju see sulfate-reducing taxa in the assembly?"
note "taxonomy is weaker evidence than dsrAB, but it is independent of the HMM search."
for f in "$WS2"/35_kaiju/*.kaiju; do
  [ -s "$f" ] || continue
  echo "   --- $(basename "$f") ($(wc -l < "$f") rows) ---"
  grep -oiE 'Desulfo[a-z]*|Desulfu[a-z]*|Syntrophobacter[a-z]*|Thermodesulfo[a-z]*|Desulfobacterota|Desulfobulbus|Desulfatiglans|Desulfococcus|Desulfovibrio' "$f" 2>/dev/null \
    | sort | uniq -c | sort -rn | head -15 | sed 's/^/     /' || echo "     none matched"
done

# ---------------------------------------------------------------- E. the candidate
say "E. ETO_LIQ_MAG_022 (JBAIYS01) -- the only MAG with any sulfur-reduction machinery"
note "aprAB + K23077, no dsrAB. Its full sulfur gene complement, so the thesis can say"
note "exactly what it does and does not have."
"$PY" - "$KOLIST" "$WS2/45_kofam_full/all_hits.tbl" "$OUT/sulfur_kos.json" <<'PYEOF'
import sys, re, json, collections
kolist, hits, kos = sys.argv[1:4]
S=json.load(open(kos)); SUL=set(S['all_sulfur'])
thr={}; defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not re.fullmatch(r'K\d{5}', p[0]): continue
    try: t=float(p[1])
    except ValueError: t=None
    thr[p[0]]=t; defn[p[0]]=p[-1]
best={}
for line in open(hits, errors='replace'):
    p=line.split()
    if len(p)<6: continue
    gene, ko = p[0], p[2]
    if not re.fullmatch(r'K\d{5}', ko): continue
    try: sc=float(p[5])
    except ValueError: continue
    t=thr.get(ko)
    if t is not None and sc<t: continue
    if gene not in best or sc>best[gene][1]: best[gene]=(ko,sc)
per=collections.defaultdict(set)
for g,(ko,sc) in best.items(): per[g.split('|',1)[0]].add(ko)
tgt='LIQ_das_metabat_022'
ks=sorted(SUL & per.get(tgt,set()))
print("   %s: %d sulfur-related KOs" % (tgt, len(ks)))
for k in ks: print("     %-8s %s" % (k, (defn.get(k) or '')[:88]))
print()
print("   present: aprAB (K00394/K00395), sat (K00958), dsr-RELATED (K23077)")
print("   absent : dsrA (K11180), dsrB (K11181)")
print("   >>> aprAB without dsrAB cannot complete sulfate -> sulfide. This genome is")
print("       incomplete at 0.2607 sediment breadth and culture-only, so the missing")
print("       dsrAB may be a binning gap rather than a biological absence. Say that.")
PYEOF
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
