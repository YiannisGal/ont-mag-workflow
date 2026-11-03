#!/usr/bin/env bash
# 49b_read_the_contigs.sh -- identify the dsrAB contigs by reading what else is on them.
#
# Read-only. Minutes.
#
# 49 SETTLED SECTION C, and it is a real result:
#   LIQ subsampled to SED's exact protein count (170,327), three independent replicates,
#   found 3 dsrAB genes every time against 4 at full depth. SED's zero is NOT a depth
#   artefact -- dsrAB is genuinely absent from the sediment assembly and present in the
#   liquid one. Worth stating plainly, because the enrichment medium contains NO SULFATE
#   (Volume VI): sulfate reducers that grow fermentatively or syntrophically can persist
#   without it, and that is the reading the chapter needs.
#
# 49 SECTIONS A AND B FAILED, both on my file assumptions -- but they SAID SO rather than
# concluding, which is the fix from 48D working as intended:
#   A  35_kaiju/LIQ.kaiju is keyed on READ UUIDs (f09cae34-3e62-...), not contig names.
#      It is a read-level run. Contig taxonomy cannot come from it at all.
#   B  my glob capped at the first 40 large fastas and missed 04_polish/LIQ_asm_final.fasta
#      -- which the script had itself listed as a candidate two lines earlier.
#
# A BETTER APPROACH THAN COMPOSITION. Two contigs carry a few dozen proteins. Search THOSE
# against the full KOfam set -- seconds -- and read what sits beside dsrAB. A ribosomal
# protein operon on the same contig names the organism far more reliably than GC does.
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
ASM="$WS2/04_polish/LIQ_asm_final.fasta"
PROT="$WS2/34_unbinned/prot/LIQ_all.faa"
KOLIST="$WS2/11_kofam/ko_list"; CH="$WS2/45_kofam_full/chunks"
BINS="$WS2/07_mags/bins"
C1=LIQ_contig_12880; C2=LIQ_contig_19299
CPU="${CPU:-$(nproc 2>/dev/null || echo 8)}"
echo "49b  $(date -Is)"

# ---------------------------------------------------------------- A. contigs
say "A. the contigs themselves, from a PINNED path"
[ -s "$ASM" ] || { bad "no $ASM"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "assembly: $ASM  ($(grep -c '^>' "$ASM") contigs)"
"$PY" - "$ASM" "$C1" "$C2" "$OUT" <<'PYEOF'
import sys, os
asm, c1, c2, out = sys.argv[1:5]
want={c1,c2}; got={}; hdr={}
name=None; buf=[]
for line in open(asm, errors='replace'):
    if line.startswith('>'):
        if name in want: got[name]=''.join(buf)
        full=line[1:].rstrip('\n'); name=full.split()[0]; hdr[name]=full; buf=[]
    else: buf.append(line.strip())
if name in want: got[name]=''.join(buf)
if not got:
    print("   NEITHER contig found in this assembly.")
    print("   headers look like: ")
    n=0
    for line in open(asm, errors='replace'):
        if line.startswith('>'):
            print("     %s" % line[1:].rstrip()[:100]); n+=1
            if n>=5: break
    sys.exit(3)
with open(os.path.join(out,'dsr_contigs.fasta'),'w') as fh:
    for k,s in got.items(): fh.write('>%s\n%s\n' % (hdr.get(k,k), s))
for k,s in got.items():
    u=s.upper(); g=u.count('G')+u.count('C'); n=g+u.count('A')+u.count('T')
    print("   %-22s length=%-9d GC=%.2f%%" % (k, len(s), 100.0*g/n if n else 0))
    print("       header: %s" % hdr.get(k,'')[:120])
print("   wrote dsr_contigs.fasta")
PYEOF
[ $? -eq 0 ] || { bad "contigs not recoverable; stopping"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- B. proteins
say "B. every protein encoded on those two contigs"
"$PY" - "$PROT" "$C1" "$C2" "$OUT" <<'PYEOF'
import sys, os, re
prot, c1, c2, out = sys.argv[1:5]
rx=re.compile(r'^(%s|%s)_(\d+)$' % (re.escape(c1), re.escape(c2)))
keep=[]; name=None; buf=[]
def flush():
    if name:
        m=rx.match(name.split()[0])
        if m: keep.append((name, ''.join(buf)))
for line in open(prot, errors='replace'):
    if line.startswith('>'):
        flush(); name=line[1:].rstrip('\n'); buf=[]
    else: buf.append(line)
flush()
p=os.path.join(out,'dsr_contig_proteins.faa')
with open(p,'w') as fh:
    for h,s in keep: fh.write('>%s\n%s' % (h,s))
from collections import Counter
c=Counter(rx.match(h.split()[0]).group(1) for h,_ in keep)
print("   proteins extracted: %d" % len(keep))
for k,v in c.items(): print("     %-22s %d proteins" % (k,v))
print("   wrote %s" % os.path.basename(p))
PYEOF

# ---------------------------------------------------------------- C. annotate
say "C. annotate them against the FULL KOfam set"
note "a few dozen proteins against 27,864 profiles -- fast, because the cost scales"
note "with sequence count once the profiles are read."
FAA="$OUT/dsr_contig_proteins.faa"
if [ ! -s "$FAA" ]; then no "no proteins extracted"; else
  OUT_TBL="$OUT/dsr_contig_kofam.tbl"
  if [ -s "$OUT_TBL" ]; then ok "already annotated"; else
    : > "$OUT_TBL.part"; RC=0
    for c in "$CH"/c*.hmm; do
      hmmsearch --cpu "$CPU" -E 1e-5 --tblout "$OUT/_t.tbl" "$c" "$FAA" > /dev/null 2>"$OUT/_t.err" || { RC=1; break; }
      grep -v '^#' "$OUT/_t.tbl" >> "$OUT_TBL.part" || true
    done
    if [ "$RC" -eq 0 ]; then mv "$OUT_TBL.part" "$OUT_TBL"; ok "annotated"; else
      rm -f "$OUT_TBL.part"; bad "a chunk failed:"; head -3 "$OUT/_t.err" | sed 's/^/     /'
    fi
    rm -f "$OUT/_t.tbl" "$OUT/_t.err"
  fi
  [ -s "$OUT_TBL" ] && "$PY" - "$KOLIST" "$OUT_TBL" <<'PYEOF'
import sys, re, collections
kolist, tbl = sys.argv[1], sys.argv[2]
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
    g,k=p[0],p[2]
    if not re.fullmatch(r'K\d{5}',k): continue
    try: s=float(p[5])
    except ValueError: continue
    t=thr.get(k)
    if t is not None and s<t: continue
    if g not in best or s>best[g][1]: best[g]=(k,s)
print("   proteins with a thresholded KO: %d" % len(best))
print()
print("   %-26s %-8s %7s  %s" % ('protein','KO','score','definition'))
for g,(k,s) in sorted(best.items(), key=lambda x:(x[0].rsplit('_',1)[0], int(x[0].rsplit('_',1)[1]))):
    mark=''
    d=defn.get(k,'')
    if k in ('K11180','K11181'): mark=' <<< dsrAB'
    elif re.search(r'ribosomal protein', d, re.I): mark=' <<< RIBOSOMAL -- taxonomic marker'
    print("   %-26s %-8s %7.1f  %s%s" % (g,k,s,d[:60],mark))
rib=[(g,k) for g,(k,s) in best.items() if re.search(r'ribosomal protein', defn.get(k,''), re.I)]
print()
print("   ribosomal proteins on these contigs: %d" % len(rib))
if rib:
    print("   >>> a ribosomal operon beside dsrAB is the strongest identification available")
    print("       without a database search. Their sequences are in dsr_contig_proteins.faa")
    print("       and can be placed against the 12_phylogeny reference set.")
else:
    print("   >>> no ribosomal proteins on these contigs, so they cannot be placed")
    print("       taxonomically from content alone. Say the source organism is undetermined.")
PYEOF
fi

# ---------------------------------------------------------------- D. composition
say "D. GC of the contigs against every bin"
"$PY" - "$OUT/dsr_contigs.fasta" "$BINS" <<'PYEOF'
import sys, os, glob
fa, bins = sys.argv[1], sys.argv[2]
def gc(s):
    u=s.upper(); g=u.count('G')+u.count('C'); n=g+u.count('A')+u.count('T')
    return 100.0*g/n if n else 0.0
seqs={}; name=None; buf=[]
for line in open(fa, errors='replace'):
    if line.startswith('>'):
        if name: seqs[name]=''.join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
if name: seqs[name]=''.join(buf)
prof=[]
for b in sorted(glob.glob(os.path.join(bins,'*.fa'))):
    tot=[]
    for line in open(b, errors='replace'):
        if not line.startswith('>'): tot.append(line.strip())
    t=''.join(tot)
    if t: prof.append((os.path.basename(b).replace('.fa',''), len(t), gc(t)))
for c,s in seqs.items():
    cg=gc(s)
    print("\n   %s  GC %.2f%%  length %d" % (c, cg, len(s)))
    for n,bp,g in sorted(prof, key=lambda x: abs(x[2]-cg))[:6]:
        star='   <-- the aprAB/dsrMK genome' if 'metabat_022' in n else ''
        print("     %-26s GC %5.2f%%  delta %4.2f  %8d bp%s" % (n,g,abs(g-cg),bp,star))
print()
print("   GC alone cannot assign a contig to a bin. If metabat_022 is not clearly")
print("   nearest, the honest statement is that the contigs were not assigned.")
PYEOF
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
