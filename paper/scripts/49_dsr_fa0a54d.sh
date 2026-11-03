#!/usr/bin/env bash
# 49_dsr_contigs.sh -- whose are the dsrAB contigs, and is SED's absence real?
#
# Read-only. Minutes.
#
# 48b FOUND: two dsrAB operons in the LIQ assembly, unbinned.
#     LIQ_contig_12880_2  K11181 dsrB  574.3    LIQ_contig_12880_3  K11180 dsrA  592.6
#     LIQ_contig_19299_13 K11181 dsrB  599.3    LIQ_contig_19299_14 K11180 dsrA  331.7
#   and ETO_LIQ_MAG_022 (JBAIYS01) carries sat + aprAB + dsrMK (K27187/K27188) + DsrC-
#   associated genes -- the whole dissimilatory pathway EXCEPT dsrAB.
#
# THREE QUESTIONS, each stated so it can come back "no":
#   A  What organism are those contigs from? Kaiju assigned them a TAXON ID; 48's grep
#      for genus NAMES in a file of numeric IDs was guaranteed to match nothing, and I
#      reported that nothing as if it meant something. Resolved through names.dmp here.
#   B  Do the contigs plausibly belong to ETO_LIQ_MAG_022? Compared on GC, length and
#      coverage against that bin's own contigs and against every other bin. If they look
#      like several bins equally, the answer is "cannot tell" and that is a real answer.
#   C  Is dsrAB genuinely absent from SED, or just under-sampled? LIQ has 270,775
#      proteins to SED's 170,327. Subsampling LIQ to SED's size and re-searching says
#      whether 4 genes would still be detectable at that depth.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/49_dsr"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp49"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
[ -d "$ENVDIR" ] && export PATH="$ENVDIR/bin:$PATH"
SUB="$WS2/48_sulfur/sulfur_profiles.hmm"; KOLIST="$WS2/11_kofam/ko_list"
CONTIGS=(LIQ_contig_12880 LIQ_contig_19299)
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"
echo "49  $(date -Is)"; echo "contigs of interest: ${CONTIGS[*]}"

# ---------------------------------------------------------------- A. taxonomy
say "A. what does Kaiju call these contigs? (taxon IDs resolved to names)"
note "48's grep looked for 'Desulfo' in a file of NUMERIC taxon ids. That could only"
note "ever return nothing, and I reported the nothing. This resolves the ids properly."
KJ=$(ls "$WS2"/35_kaiju/LIQ.kaiju 2>/dev/null | head -1)
NAMES=$(find "$WS2" "$HOME" -maxdepth 6 -name 'names.dmp' 2>/dev/null | head -1)
echo "   kaiju: ${KJ:-none}"
echo "   names.dmp: ${NAMES:-not found}"
if [ -n "${KJ:-}" ]; then
  note "format of the kaiju file (first 3 rows):"
  head -3 "$KJ" | cut -c1-140 | sed 's/^/     /'
  for c in "${CONTIGS[@]}"; do
    echo "   --- $c ---"
    grep -m 20 -E "(^|[[:space:]])${c}(_[0-9]+)?([[:space:]]|$)" "$KJ" 2>/dev/null | head -8 | cut -c1-160 | sed 's/^/     /' \
      || echo "     no rows for $c"
  done
  if [ -n "${NAMES:-}" ]; then
    note "resolving every taxon id seen for those contigs:"
    "$PY" - "$KJ" "$NAMES" "${CONTIGS[@]}" <<'PYEOF'
import sys, re, collections
kj, names = sys.argv[1], sys.argv[2]
contigs = sys.argv[3:]
want=set()
rx=re.compile(r'(%s)_\d+' % '|'.join(re.escape(c) for c in contigs))
rows=collections.defaultdict(collections.Counter)
for line in open(kj, errors='replace'):
    p=line.rstrip('\n').split('\t')
    if len(p) < 3: continue
    m=rx.search(p[1] if len(p)>1 else '')
    if not m: continue
    tid=p[2].strip()
    if tid.isdigit(): rows[m.group(1)][tid]+=1; want.add(tid)
if not want:
    print("     no taxon ids found for those contigs in this file")
    sys.exit(0)
nm={}
for line in open(names, errors='replace'):
    p=[x.strip() for x in line.split('|')]
    if len(p)>=4 and p[0] in want and p[3]=='scientific name': nm[p[0]]=p[1]
for c, cnt in rows.items():
    print("     %s:" % c)
    for tid,n in cnt.most_common(6):
        print("       taxid %-9s n=%-4d %s" % (tid, n, nm.get(tid,'(name not found)')))
PYEOF
  else
    no "names.dmp not found -- taxon ids cannot be resolved to names here"
    note "the ids above are still the answer; they just need a lookup"
  fi
fi

# ---------------------------------------------------------------- B. binning
say "B. do the dsrAB contigs look like ETO_LIQ_MAG_022?"
note "GC and length against every bin. This can legitimately return 'cannot tell'."
BINS=$(find "$WS2/07_mags" -maxdepth 2 -type d -name 'bins' | head -1)
ASM=$(find "$WS2" -maxdepth 3 -type f \( -name 'LIQ*final*.fa*' -o -name 'LIQ*contigs*.fa*' -o -name 'LIQ*assembly*.fa*' \) 2>/dev/null | head -3)
echo "   bins dir: ${BINS:-none}"
echo "   candidate LIQ assemblies:"; printf '     %s\n' ${ASM:-none}
"$PY" - "$OUT" "${BINS:-none}" "$WS2" "${CONTIGS[@]}" <<'PYEOF'
import sys, os, glob, re
out, bins, ws2 = sys.argv[1], sys.argv[2], sys.argv[3]
contigs=set(sys.argv[4:])
def gc(seq):
    s=seq.upper(); g=s.count('G')+s.count('C'); n=g+s.count('A')+s.count('T')
    return 100.0*g/n if n else 0.0
# find the contigs in whatever assembly holds them
found={}
cands=[]
for pat in ('*.fa','*.fna','*.fasta'):
    cands += glob.glob(os.path.join(ws2,'*',pat)) + glob.glob(os.path.join(ws2,'*','*',pat))
cands=[c for c in cands if os.path.getsize(c) > 5_000_000][:40]
print("   scanning %d large fasta file(s) for the contigs" % len(cands))
for f in cands:
    name=None; buf=[]
    try: fh=open(f, errors='replace')
    except Exception: continue
    for line in fh:
        if line.startswith('>'):
            if name in contigs and name not in found:
                found[name]=(f, ''.join(buf))
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if name in contigs and name not in found: found[name]=(f, ''.join(buf))
    fh.close()
    if len(found)==len(contigs): break
if not found:
    print("   contigs not located in any scanned assembly -- cannot compare composition")
    sys.exit(0)
print()
for c,(f,seq) in found.items():
    print("   %-22s len=%-9d GC=%.1f%%   in %s" % (c, len(seq), gc(seq), os.path.basename(f)))
# bin GC profiles
if bins != 'none' and os.path.isdir(bins):
    print()
    print("   %-26s %8s %8s %7s" % ('bin','n_contig','bp','GC%'))
    prof=[]
    for b in sorted(glob.glob(os.path.join(bins,'*.fa'))):
        seqs=[]; buf=[]
        for line in open(b, errors='replace'):
            if line.startswith('>'):
                if buf: seqs.append(''.join(buf)); buf=[]
            else: buf.append(line.strip())
        if buf: seqs.append(''.join(buf))
        tot=''.join(seqs)
        if not tot: continue
        prof.append((os.path.basename(b).replace('.fa',''), len(seqs), len(tot), gc(tot)))
    for n,nc,bp,g in sorted(prof, key=lambda x:-x[2])[:40]:
        print("   %-26s %8d %8d %7.1f" % (n,nc,bp,g))
    print()
    for c,(f,seq) in found.items():
        cg=gc(seq)
        near=sorted(prof, key=lambda x: abs(x[3]-cg))[:5]
        print("   %s (GC %.1f%%) is closest in GC to:" % (c, cg))
        for n,nc,bp,g in near:
            print("       %-26s GC %.1f%%  (delta %.1f)%s" % (n,g,abs(g-cg),
                  '   <-- the aprAB/dsrMK genome' if 'metabat_022' in n else ''))
        print("   >>> GC alone is weak evidence. If metabat_022 is not clearly nearest,")
        print("       say the contigs could not be assigned and leave it there.")
PYEOF

# ---------------------------------------------------------------- C. depth
say "C. is SED's lack of dsrAB real, or a depth artefact?"
note "subsample LIQ proteins down to SED's count and re-search. If 4 dsrAB genes are"
note "still found at SED's depth, then SED's absence is meaningful. If they vanish,"
note "the comparison cannot be made and that must be said."
LIQ="$WS2/34_unbinned/prot/LIQ_all.faa"; SED="$WS2/34_unbinned/prot/SED_all.faa"
if [ -s "$LIQ" ] && [ -s "$SED" ] && [ -s "$SUB" ]; then
  NS=$(grep -c '^>' "$SED"); NL=$(grep -c '^>' "$LIQ")
  echo "   LIQ $NL proteins, SED $NS proteins  (ratio $(echo "scale=2; $NL/$NS" | bc))"
  for REP in 1 2 3; do
    SUBF="$OUT/liq_sub_${REP}.faa"
    [ -s "$SUBF" ] || "$PY" - "$LIQ" "$NS" "$REP" > "$SUBF" <<'PYEOF'
import sys, random
src, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
recs=[]; name=None; buf=[]
for line in open(src, errors='replace'):
    if line.startswith('>'):
        if name: recs.append((name, ''.join(buf)))
        name=line; buf=[]
    else: buf.append(line)
if name: recs.append((name, ''.join(buf)))
random.seed(seed); random.shuffle(recs)
for h,s in recs[:n]: sys.stdout.write(h); sys.stdout.write(s)
PYEOF
    O="$OUT/liq_sub_${REP}.tbl"
    if [ ! -s "$O" ]; then
      hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$O.part" "$SUB" "$SUBF" >/dev/null 2>"$OUT/sub${REP}.err" \
        && mv "$O.part" "$O" || { no "replicate $REP failed"; head -3 "$OUT/sub${REP}.err" | sed 's/^/       /'; continue; }
    fi
    N=$("$PY" - "$KOLIST" "$O" <<'PYEOF'
import sys, re
kolist, tbl = sys.argv[1], sys.argv[2]
thr={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not re.fullmatch(r'K\d{5}', p[0]): continue
    try: thr[p[0]]=float(p[1])
    except ValueError: thr[p[0]]=None
best={}
for line in open(tbl, errors='replace'):
    if line.startswith('#'): continue
    p=line.split()
    if len(p)<6: continue
    g,k=p[0],p[2]
    if not re.fullmatch(r'K\d{5}',k): continue
    try: s=float(p[5])
    except ValueError: continue
    t=thr.get(k)
    if t is not None and s<t: continue
    if g not in best or s>best[g][1]: best[g]=(k,s)
print(sum(1 for k,_ in best.values() if k in ('K11180','K11181')))
PYEOF
)
    printf '   replicate %s at SED depth: %s dsrAB gene(s)\n' "$REP" "$N"
  done
  note "full-depth LIQ had 4. If the subsamples keep finding dsrAB, SED's zero is a"
  note "real difference; if they lose it, the two assemblies cannot be compared."
else no "inputs missing -- section C cannot run"; fi
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
