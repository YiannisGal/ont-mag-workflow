#!/usr/bin/env bash
# 48_sulfate_reducers.sh -- is dsrAB anywhere: in the MAGs, or in the unbinned fraction?
#
# RUN IN tmux if section D has to search a large assembly proteome; otherwise minutes.
#
# WHY SECTION A COMES FIRST -- I MAY HAVE MADE THE OPPOSITE ERROR
#   Volume XI states "dsrA and dsrB appear in none of the 29 genomes". That rests on my
#   regex  'dissimilatory sulfite reductase'.  KEGG's description for dsrA is written
#   "sulfite reductase, dissimilatory-type subunit alpha" -- DIFFERENT WORD ORDER, which
#   my pattern does not match. So the claim could be an artefact of the pattern rather
#   than a fact about the genomes.
#   Section A enumerates EVERY sulfur KO in ko_list and prints it, so the KO numbers are
#   read off the database instead of assumed. Only then is the MAG question re-asked.
#
# THE ORDER MATTERS. If A shows dsrAB has a KO my regex missed, then Volume XI's headline
# is wrong and must be withdrawn before anything is built on it. If A confirms the regex
# caught the right KOs, the absence stands and section D asks where the sulfate reducers
# actually are.
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
KOLIST="$WS2/11_kofam/ko_list"; EXT="$WS2/45_kofam_full/profiles_full"
HITS="$WS2/45_kofam_full/all_hits.tbl"; XW="$WS2/07_mags/MAG_CROSSWALK.tsv"
XAB="$WS2/26_crossmap/cross_abundance.tsv"
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"
echo "48  $(date -Is)"

# ---------------------------------------------------------------- A. the KOs
say "A. every sulfur KO in ko_list, printed -- no pattern assumed"
"$PY" - "$KOLIST" "$OUT" <<'PYEOF'
import sys, re, os, json
kolist, out = sys.argv[1], sys.argv[2]
rows=[]
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not re.fullmatch(r'K\d{5}', p[0]): continue
    d=p[-1] or ''
    if re.search(r'sulfite reductase|sulfate adenylyl|adenylylsulfate|thiosulfate|sulfur|sulfide|sulfhydrylase|APS reductase|sulfate reduct', d, re.I):
        rows.append((p[0], d))
print("   sulfur-related KOs in ko_list: %d" % len(rows))
print()
print("   --- every KO whose description mentions 'sulfite reductase' ---")
sr=[(k,d) for k,d in rows if re.search(r'sulfite reductase', d, re.I)]
for k,d in sorted(sr): print("     %-8s %s" % (k, d[:104]))
print()
print("   --- of those, which are DISSIMILATORY-type ---")
dis=[(k,d) for k,d in sr if re.search(r'dissimilat', d, re.I)]
for k,d in sorted(dis): print("     %-8s %s" % (k, d[:104]))
print("     count: %d" % len(dis))
print()
print("   --- adenylylsulfate / APS reductase (aprAB) ---")
apr=[(k,d) for k,d in rows if re.search(r'adenylylsulfate reductase|APS reductase', d, re.I)]
for k,d in sorted(apr): print("     %-8s %s" % (k, d[:104]))
print()
print("   --- sulfate adenylyltransferase (sat) -- ASSIMILATORY AND DISSIMILATORY BOTH ---")
sat=[(k,d) for k,d in rows if re.search(r'sulfate adenylyltransferase', d, re.I)]
for k,d in sorted(sat): print("     %-8s %s" % (k, d[:104]))
print()
old=re.compile(r'dissimilatory sulfite reductase', re.I)
missed=[(k,d) for k,d in dis if not old.search(d)]
print("   >>> MY VOLUME XI REGEX WAS: 'dissimilatory sulfite reductase'")
print("   >>> dissimilatory-type sulfite reductase KOs it would have MISSED: %d" % len(missed))
for k,d in missed: print("       %-8s %s" % (k, d[:100]))
if missed:
    print("   >>> Volume XI's 'no dsrAB' claim was tested with a pattern that cannot")
    print("       match these. Section B re-asks the question with the right KOs.")
else:
    print("   >>> the pattern caught every dissimilatory-type KO; the claim was sound.")
json.dump({'dsr_dissimilatory':[k for k,_ in dis], 'apr':[k for k,_ in apr],
           'sat':[k for k,_ in sat], 'all_sulfur':[k for k,_ in rows]},
          open(os.path.join(out,'sulfur_kos.json'),'w'), indent=2)
PYEOF

# ---------------------------------------------------------------- B. MAGs
say "B. re-ask the MAG question with the KO numbers from A"
"$PY" - "$OUT/sulfur_kos.json" "$HITS" "$KOLIST" "$XW" "$XAB" <<'PYEOF'
import sys, json, re, csv, collections, os
kos, hits, kolist, xw, xab = sys.argv[1:6]
S=json.load(open(kos))
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
    if t is not None and sc < t: continue
    if gene not in best or sc > best[gene][1]: best[gene]=(ko,sc)
per=collections.defaultdict(set)
for g,(ko,sc) in best.items(): per[g.split('|',1)[0]].add(ko)
bin2mag={r.get('working_name',''):r.get('mag_id','') for r in csv.DictReader(open(xw), delimiter='\t')}
meta={r['mag_id']:r for r in csv.DictReader(open(xab), delimiter='\t')}
DSR=set(S['dsr_dissimilatory']); APR=set(S['apr']); SAT=set(S['sat'])
print("   dissimilatory sulfite reductase KOs searched: %s" % sorted(DSR))
print("   aprAB KOs searched                          : %s" % sorted(APR))
print()
hit_d=[(g,sorted(DSR&s)) for g,s in per.items() if DSR&s]
hit_a=[(g,sorted(APR&s)) for g,s in per.items() if APR&s]
print("   genomes with a DISSIMILATORY sulfite reductase KO: %d" % len(hit_d))
for g,k in hit_d:
    m=bin2mag.get(g,'?'); md=meta.get(m,{})
    print("     %-24s %-16s %-18s breadth=%-7s %s" % (g,m,md.get('genus','-'),md.get('sed_breadth','-'),k))
    for kk in k: print("         %s  %s" % (kk, (defn.get(kk) or '')[:80]))
print()
print("   genomes with aprAB: %d" % len(hit_a))
for g,k in hit_a:
    m=bin2mag.get(g,'?'); md=meta.get(m,{})
    print("     %-24s %-16s %-18s breadth=%-7s %s" % (g,m,md.get('genus','-'),md.get('sed_breadth','-'),k))
print()
both=[g for g,s in per.items() if (DSR&s) and (APR&s)]
print("   >>> genomes with BOTH dsr and apr (a real sulfate reducer): %d %s" % (len(both), both))
if not both:
    print("   >>> confirmed: no MAG carries the dissimilatory pair. Sulfate reduction is")
    print("       not encoded in any recovered genome.")
PYEOF

# ---------------------------------------------------------------- C. unbinned
say "C. where is the unbinned fraction?"
note "assembly-level proteomes contain everything, binned and not. 45b saw prefixes"
note "like LIQ_contig_2376, so such files exist."
mapfile -t ASM < <(find "$WS2" -type f -name '*.faa' 2>/dev/null \
  | xargs -I{} sh -c 'n=$(grep -c "^>" "{}" 2>/dev/null || echo 0); [ "$n" -gt 20000 ] && echo "$n {}"' 2>/dev/null | sort -rn | head -10)
if [ "${#ASM[@]}" -eq 0 ]; then no "no large (>20k protein) assembly proteome found"; else
  printf '  %10s  %s\n' "proteins" "file"
  for a in "${ASM[@]}"; do printf '  %10s  %s\n' "${a%% *}" "${a#* }"; done
fi
note "and any Kaiju output from step 36b:"
find "$WS2" -maxdepth 3 -type f \( -name '*kaiju*' -o -path '*36b*' -name '*.tsv' -o -path '*36b*' -name '*.out' \) 2>/dev/null | head -10 | sed 's/^/     /'

# ---------------------------------------------------------------- D. search
say "D. search the assembly proteome for dsr and apr directly"
TARGET=$(printf '%s\n' "${ASM[@]:-}" | head -1 | sed 's/^[0-9]* //')
if [ -z "$TARGET" ] || [ ! -s "$TARGET" ]; then
  no "no assembly proteome to search -- section D cannot run"
else
  echo "   target: $TARGET ($(grep -c '^>' "$TARGET") proteins)"
  SUB="$OUT/sulfur_profiles.hmm"
  "$PY" - "$OUT/sulfur_kos.json" "$EXT" "$SUB" <<'PYEOF'
import sys, json, os, glob
kos, ext, out = sys.argv[1:4]
S=json.load(open(kos))
want=set(S['dsr_dissimilatory']) | set(S['apr']) | set(S['sat'])
found=[]
with open(out,'w') as fh:
    for k in sorted(want):
        p=os.path.join(ext,'profiles',k+'.hmm')
        if not os.path.exists(p):
            hits=glob.glob(os.path.join(ext,'**',k+'.hmm'), recursive=True)
            p=hits[0] if hits else None
        if p and os.path.exists(p):
            fh.write(open(p, errors='replace').read()); found.append(k)
print("   profiles assembled: %d of %d wanted" % (len(found), len(want)))
print("   %s" % found)
PYEOF
  if [ -s "$SUB" ]; then
    note "hmmsearch (a handful of profiles, so this is quick)"
    T0=$(date +%s)
    hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$OUT/assembly_sulfur.tbl" "$SUB" "$TARGET" > /dev/null 2> "$OUT/assembly_sulfur.err"
    RC=$?; echo "   exit=$RC in $(( $(date +%s)-T0 ))s"
    [ "$RC" -ne 0 ] && { bad "hmmsearch failed:"; head -10 "$OUT/assembly_sulfur.err" | sed 's/^/     /'; }
    "$PY" - "$KOLIST" "$OUT/assembly_sulfur.tbl" "$OUT/sulfur_kos.json" <<'PYEOF'
import sys, re, json, collections
kolist, tbl, kos = sys.argv[1:4]
S=json.load(open(kos)); DSR=set(S['dsr_dissimilatory']); APR=set(S['apr'])
thr={}; defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not re.fullmatch(r'K\d{5}', p[0]): continue
    try: t=float(p[1])
    except ValueError: t=None
    thr[p[0]]=t; defn[p[0]]=p[-1]
best={}
for line in open(tbl, errors='replace'):
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
print()
print("   thresholded hits in the assembly:")
for ko,n in cnt.most_common():
    tag=''
    if ko in DSR: tag='   <-- DISSIMILATORY SULFITE REDUCTASE'
    if ko in APR: tag='   <-- aprAB'
    print("     %-8s %4d  %s%s" % (ko, n, (defn.get(ko) or '')[:60], tag))
nd=sum(n for ko,n in cnt.items() if ko in DSR); na=sum(n for ko,n in cnt.items() if ko in APR)
print()
print("   dissimilatory sulfite reductase genes in the assembly: %d" % nd)
print("   aprAB genes in the assembly                          : %d" % na)
print()
if nd > 0:
    print("   >>> SULFATE REDUCERS ARE PRESENT IN THE ASSEMBLY but were not binned.")
    print("       The chapter's competition premise stands; state that SRB were detected")
    print("       in the assembly and not recovered as MAGs.")
    print("   >>> the contigs carrying them:")
    for gene,(ko,sc) in sorted(best.items()):
        if ko in DSR: print("       %-40s %s  score=%.1f" % (gene, ko, sc))
else:
    print("   >>> NO dissimilatory sulfite reductase in the assembly either.")
    print("       Then the premise cannot rest on this dataset at all, and the chapter")
    print("       must argue sulfate reduction from the geochemistry alone.")
PYEOF
  fi
fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
