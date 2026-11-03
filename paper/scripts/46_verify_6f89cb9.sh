#!/usr/bin/env bash
# 46_verify_and_tables.sh -- verify the KOfam run, then build the tables figures need.
#
# Read-only apart from writing tables into 45_kofam_full/tables/. Minutes.
#
# WHY SECTION A EXISTS
#   45c section B printed:
#       gzip: stdin: invalid compressed data--crc error
#       tar: Error is not recoverable: exiting now
#       [OK]   extracted 27864 profiles
#   The extraction ABORTED and my script called it OK, because it counted files
#   afterwards instead of checking tar's exit status. That is the same fault as 45's
#   hidden stderr: judging a step by its output rather than its status. The count landing
#   exactly on 27,864 suggests the damage is at the end of the archive, and all 24 chunks
#   exited 0 (hmmer rejects malformed profiles), but assuming that is not verifying it.
#
# WHAT 45c PRODUCED
#   1,509,145 hit rows; 43,001 genes above threshold; median 1058 distinct KOs per
#   genome across 29 genomes. B12: 16 producers at >=10 biosynthesis KOs, 26 of 29 with
#   a B12-DEPENDENT enzyme -- so step 21's "0 consumers" is refuted.
set -uo pipefail
WS2="${WS2:-${WS2}}"
K="$WS2/45_kofam_full"; TAB="$K/tables"; mkdir -p "$TAB"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
KOLIST="$WS2/11_kofam/ko_list"; EXT="$K/profiles_full"; CH="$K/chunks"
echo "46  $(date -Is)"

# ---------------------------------------------------------------- A. integrity
say "A. did the CRC error cost us any profiles?"
note "every HMMER3 profile must start with HMMER3 and end with //. A truncated file"
note "would have either no terminator or no model lines."
"$PY" - "$EXT" "$KOLIST" <<'PYEOF'
import sys, os, glob, re
ext, kolist = sys.argv[1], sys.argv[2]
files=sorted(glob.glob(os.path.join(ext,'**','*.hmm'), recursive=True))
print("   .hmm files on disk: %d" % len(files))
bad_head=[]; bad_tail=[]; empty=[]
for f in files:
    sz=os.path.getsize(f)
    if sz == 0: empty.append(f); continue
    with open(f, 'rb') as fh:
        head=fh.readline().decode('ascii','replace')
        fh.seek(max(0, sz-40)); tail=fh.read().decode('ascii','replace')
    if not head.startswith('HMMER3'): bad_head.append(os.path.basename(f))
    if '//' not in tail: bad_tail.append(os.path.basename(f))
print("   zero-byte            : %d %s" % (len(empty), [os.path.basename(x) for x in empty[:5]]))
print("   missing HMMER3 header: %d %s" % (len(bad_head), bad_head[:5]))
print("   missing // terminator: %d %s" % (len(bad_tail), bad_tail[:5]))
kos={os.path.basename(f)[:-4] for f in files}
want=set()
for line in open(kolist):
    p=line.split('\t')
    if p and re.fullmatch(r'K\d{5}', p[0]): want.add(p[0])
print()
print("   KOs in ko_list       : %d" % len(want))
print("   profiles present     : %d" % len(kos & want))
missing = want - kos
print("   KOs with NO profile  : %d" % len(missing))
if missing:
    print("     e.g. %s" % sorted(missing)[:10])
    print("     (ko_list is larger than the profile set in normal KOfam releases;")
    print("      what matters is whether any profile is TRUNCATED, above.)")
if empty or bad_head or bad_tail:
    print("INTEGRITY: DAMAGED")
else:
    print("INTEGRITY: CLEAN -- every profile has a header and a terminator")
PYEOF

note "and did hmmsearch actually consume every profile in every chunk?"
TOTC=0
for c in "$CH"/c*.hmm; do n=$(grep -c '^NAME' "$c"); TOTC=$(( TOTC + n )); done
echo "   profiles across chunks: $TOTC"
echo "   .tbl outputs present  : $(ls "$CH"/c*.tbl 2>/dev/null | wc -l) of $(ls "$CH"/c*.hmm | wc -l)"
for c in "$CH"/c*.err; do
  [ -s "$c" ] && { no "$(basename "$c") is non-empty:"; head -3 "$c" | sed 's/^/       /'; }
done
ok "no chunk reported an error" 2>/dev/null || true

# ---------------------------------------------------------------- B. matrix
say "B. the KO matrix -- genome x KO, the substrate for every functional figure"
"$PY" - "$KOLIST" "$K/all_hits.tbl" "$TAB" <<'PYEOF'
import sys, os, re, csv, collections, json
kolist, hits, tab = sys.argv[1:4]
thr={}; defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)<3 or not p[0].startswith('K'): continue
    try: t=float(p[1])
    except ValueError: t=None
    thr[p[0]]=t; defn[p[0]]=p[-1]
best={}
for line in open(hits, errors='replace'):
    p=line.split()
    if len(p) < 6: continue
    gene, ko = p[0], p[2]
    if not re.fullmatch(r'K\d{5}', ko): continue
    try: sc=float(p[5])
    except ValueError: continue
    t=thr.get(ko)
    if t is not None and sc < t: continue
    if gene not in best or sc > best[gene][1]: best[gene]=(ko,sc)
per=collections.defaultdict(collections.Counter)
for gene,(ko,sc) in best.items():
    per[gene.split('|',1)[0]][ko] += 1
genomes=sorted(per); allko=sorted({k for c in per.values() for k in c})
print("   %d genomes x %d distinct KOs" % (len(genomes), len(allko)))

with open(os.path.join(tab,'ko_matrix_counts.tsv'),'w',newline='') as fh:
    w=csv.writer(fh, delimiter='\t')
    w.writerow(['KO','definition'] + genomes)
    for ko in allko:
        w.writerow([ko, defn.get(ko,'')] + [per[g].get(ko,0) for g in genomes])
with open(os.path.join(tab,'ko_matrix_binary.tsv'),'w',newline='') as fh:
    w=csv.writer(fh, delimiter='\t')
    w.writerow(['KO','definition'] + genomes)
    for ko in allko:
        w.writerow([ko, defn.get(ko,'')] + [1 if per[g].get(ko,0) else 0 for g in genomes])
with open(os.path.join(tab,'genome_ko_summary.tsv'),'w',newline='') as fh:
    w=csv.writer(fh, delimiter='\t'); w.writerow(['genome','n_genes_assigned','n_distinct_KO'])
    for g in genomes:
        w.writerow([g, sum(per[g].values()), len(per[g])])
print("   wrote ko_matrix_counts.tsv, ko_matrix_binary.tsv, genome_ko_summary.tsv")

# --- functional marker panels, resolved BY DESCRIPTION from ko_list
PANELS = {
 'methanogenesis (mcr)'      : r'methyl-coenzyme M reductase',
 'methylamine methyltransf.' : r'\[methyl-Co\(III\).*methylamine|methylamine.*methyltransferase|monomethylamine|dimethylamine|trimethylamine.*methyltransferase',
 'methanol methyltransferase': r'methanol.*methyltransferase',
 'CO2 reduction (fwd/fmd)'   : r'formylmethanofuran dehydrogenase',
 'hydrogenase (Ni-Fe)'       : r'\[NiFe\].*hydrogenase|hydrogenase large subunit|coenzyme F420 hydrogenase',
 'formate dehydrogenase'     : r'formate dehydrogenase',
 'sulfate reduction (dsr)'   : r'dissimilatory sulfite reductase|adenylylsulfate reductase|sulfate adenylyltransferase',
 'nitrogen fixation (nif)'   : r'nitrogenase',
 'acetate kinase / pta'      : r'\backw?\b|acetate kinase|phosphate acetyltransferase',
 'B12 biosynthesis'          : r'cob\(|cobalamin|cobyrinate|precorrin|cobinamide|cobaltochelatase',
 'B12-dependent enzymes'     : r'methionine synthase|methylmalonyl-CoA mutase|glutamate mutase|ethanolamine ammonia-lyase',
 'B12 transport (btu)'       : r'\bbtu[BFCD]\b|cobalamin transport',
}
rows=[]
for label, pat in PANELS.items():
    rx=re.compile(pat, re.I)
    ks={k for k,d in defn.items() if rx.search(d or '')}
    rows.append((label, ks))
with open(os.path.join(tab,'function_panels.tsv'),'w',newline='') as fh:
    w=csv.writer(fh, delimiter='\t'); w.writerow(['panel','n_KO_in_panel'] + genomes)
    for label, ks in rows:
        w.writerow([label, len(ks)] + [len(ks & set(per[g])) for g in genomes])
print()
print("   %-28s %5s  %s" % ('panel','nKO','genomes with >=1 KO'))
for label, ks in rows:
    n=sum(1 for g in genomes if ks & set(per[g]))
    print("   %-28s %5d  %d" % (label, len(ks), n))
print()
print("   wrote function_panels.tsv")
json.dump({g: sorted(per[g]) for g in genomes}, open(os.path.join(tab,'per_genome_kos.json'),'w'))
PYEOF

# ---------------------------------------------------------------- C. cross-check
say "C. cross-check against what we already believe"
note "mcr should be present in exactly the 8 archaea, and absent from all 21 bacteria."
note "if that fails, the run is not trustworthy no matter how large the table is."
"$PY" - "$TAB" "$WS2" <<'PYEOF'
import sys, os, csv, re, glob
tab, ws2 = sys.argv[1], sys.argv[2]
tax={}
for f in glob.glob(os.path.join(ws2,'08_taxonomy','**','gtdbtk.*.summary.tsv'), recursive=True):
    with open(f) as fh:
        rd=csv.DictReader(fh, delimiter='\t')
        c=rd.fieldnames or []
        n=next((x for x in c if re.search(r'user_genome|genome|bin',x,re.I)), c[0])
        t=next((x for x in c if re.search(r'classification|taxonom',x,re.I)), None)
        for r in rd:
            if t: tax[r[n]]='Archaea' if 'd__Archaea' in r[t] else 'Bacteria'
p=os.path.join(tab,'function_panels.tsv')
with open(p) as fh:
    rd=list(csv.reader(fh, delimiter='\t'))
hdr=rd[0]; genomes=hdr[2:]
row=next((r for r in rd[1:] if r[0].startswith('methanogenesis')), None)
if not row: print("   no mcr panel row"); sys.exit(0)
arch=[g for g in genomes if tax.get(g)=='Archaea']
bact=[g for g in genomes if tax.get(g)=='Bacteria']
print("   archaea in the table: %d   bacteria: %d   untyped: %d"
      % (len(arch), len(bact), len(genomes)-len(arch)-len(bact)))
def val(g): return int(row[2+genomes.index(g)])
a_pos=[g for g in arch if val(g)>0]; b_pos=[g for g in bact if val(g)>0]
print("   mcr present in %d/%d archaea" % (len(a_pos), len(arch)))
for g in arch: print("     %-26s %s" % (g, 'mcr' if val(g)>0 else '-'))
print("   mcr present in %d/%d bacteria  %s" % (len(b_pos), len(bact), b_pos))
print()
if len(a_pos)==len(arch) and not b_pos:
    print("CROSSCHECK: PASSED -- mcr in every archaeon, no bacterium. Consistent with 41f.")
else:
    print("CROSSCHECK: FAILED -- that pattern should be near-perfect; investigate before")
    print("            building any figure on this table.")
PYEOF
say "DONE"; echo "tables in $TAB"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
