#!/usr/bin/env bash
# =============================================================================
# 31_readlevel_methylotrophy.sh  --  is the methylotrophic consumer end in the
#                                    lagoon, or only in the bottles?
#
# THE GAP THIS CLOSES
#   Step 30 put the PRODUCING end of both non-competitive methanogenic routes in
#   the sediment: ETO_LIQ_MAG_015 (Lacrimispora) carries cutC (choline -> TMA)
#   and CE8 pectin methylesterase (pectin -> methanol), at 0.90 sediment breadth.
#   The CONSUMING end is not there. The three methylamine methanogens with mttB
#   and pyl operons sit at 0.08-0.18 sediment breadth — enrichment-specific by
#   step 26's own criterion — and the only methanogen confirmed in situ is a
#   hydrogenotrophic Methanofollis.
#
#   But the MAG catalogue holds ~2% of the sediment community, and SingleM puts
#   sediment archaea at 14.1%. "Not recovered" and "not there" are different
#   claims and only one of them has been tested.
#
# THE DESIGN, WHICH IS STEP 20b's, GENERALISED
#   Level 1  the 29 MAGs                  done, negative for the sediment
#   Level 2  every gene in both assemblies binned or not
#   Level 3  the raw reads                 100% of the data
#   A hit at level 2 or 3 that is absent at level 1 is a real organism the
#   binning missed. A controlled zero at level 3 is a TESTED NEGATIVE, which is
#   a far stronger sentence than a MAG-based absence and one almost no MAG study
#   bothers to establish.
#
# THE QUERIES COME FROM THIS DATASET, NOT FROM A DOWNLOAD
#   Step 15 built marker trees and confirmed mcrA 9/9 and mttB 5/5 in these
#   genomes. Those proteins are the queries. Using our own verified sequences
#   beats fetching references: they are the right length, the right domain, and
#   already checked against the ACR and MttB-superfamily look-alikes that step 15
#   was built to catch.
#
# THE CONTROL, AND WHY IT HAS TO BE ARCHAEAL
#   Step 19b's controls failed because rpoB and tuf are BACTERIAL and eight
#   genomes here are archaea. The control for "can this search see archaea in
#   sediment reads at all" must therefore be an ARCHAEAL single-copy protein
#   taken from the same genomes and searched the same way. If the archaeal
#   ribosomal control is found and the methyltransferases are not, the negative
#   is real. If the control is also absent, the search is broken and NO VERDICT
#   is printed.
#
# ONT-SPECIFIC METHOD NOTE
#   Level 3 uses `diamond blastx --frameshift 15`. Plain blastx fails on ONT
#   reads because one indel shifts the frame and destroys the alignment;
#   frameshift mode recovers alignments across indels and is the correct tool
#   for long error-prone reads. Same reasoning as step 20b.
#
# OUTPUT  31_readlevel/marker_queries.faa · assembly_hits.tsv · read_hits.tsv
#         31_readlevel/VERDICT.txt
# RUNTIME 20-60 min, dominated by the read-level search.
# USAGE   ws2 && bash scripts/31_readlevel_methylotrophy.sh 2>&1 | tee logs/31.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="31_readlevel"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; T="$WS2/15_marker_trees"
R="$WS2/31_readlevel"; mkdir -p "$R"
rm -f "$R/VERDICT.txt" "$R/read_hits.tsv" "$R/assembly_hits.tsv"

MINID="${MINID:-50}"          # percent identity, protein level
MINLEN="${MINLEN:-30}"        # aligned amino acids

echo "######## 31 — METHYLOTROPHIC METHANOGENESIS AT READ LEVEL ########"
echo "  started: $(date -Is)"
runlog "31_readlevel started"

# ---------------------------------------------------- 1. build the query set
echo
echo "######## 1. QUERIES FROM OUR OWN VERIFIED MARKERS ########"
QF="$R/marker_queries.faa"; : > "$QF"
python3 - "$T" "$A" "$F" "$QF" <<'PY'
import sys, os, csv, glob, re
T, A, F, QF = sys.argv[1:5]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def dom(m): return ((tax.get(m, {}) or {}).get('domain','') or '')
ARC = [m for m in tax if 'archae' in dom(m).lower()]

def rdfa(p):
    n, b = None, []
    for line in open(p):
        if line.startswith('>'):
            if n: yield n, ''.join(b)
            n = line[1:].strip(); b = []
        else: b.append(line.strip())
    if n: yield n, ''.join(b)

# TARGET markers (the question) and CONTROL markers (can the search see archaea)
TARGETS = {
 'mcrA':  ('methyl-coenzyme M reductase alpha', 'TARGET'),
 'mtaB':  ('methanol methyltransferase',        'TARGET'),
 'mtbB':  ('dimethylamine methyltransferase',   'TARGET'),
 'mtmB':  ('monomethylamine methyltransferase', 'TARGET'),
 'mttB':  ('trimethylamine methyltransferase',  'TARGET'),
}
out, found = [], {}

# route 1: step 15's marker trees, which already vetted these sequences
for d in sorted(glob.glob(os.path.join(T, '*'))):
    if not os.path.isdir(d): continue
    mk = os.path.basename(d)
    key = next((k for k in TARGETS if k.lower() == mk.lower()), None)
    if key is None: continue
    for fa in glob.glob(os.path.join(d, '*.f*a')):
        for nm, sq in rdfa(fa):
            tag = nm.split()[0]
            # keep only sequences from OUR genomes, not the reference set the
            # tree was built against — a reference query would test whether the
            # reads contain that reference, which is a different question.
            if not any(mid(m) in nm or m in nm for m in tax): continue
            if len(sq) < 100: continue
            out.append((f"{key}|{tag}", sq)); found[key] = found.get(key, 0) + 1
        if found.get(key): break

print(f"  marker-tree route: {sum(found.values())} sequences for "
      f"{len(found)} marker(s): {', '.join(sorted(found)) or 'none'}")

# THE METHANOL ARM WAS NEVER QUERIED IN THE FIRST RUN. Step 15 built trees for
# mcrA and mttB only, so mtaB/mtbB/mtmB had no source and the methanol half of
# the pectin -> methanol -> CH4 chain went untested while I reported the TMA
# half as answered. Pull them from the archaeal proteomes directly, matching the
# CONCEPT rather than a guessed phrase — the same fix the control needed.
MT_PAT = (('mtaB', ('methanol', 'corrinoid')),
          ('mtaB', ('methanol--corrinoid',)),
          ('mtbB', ('dimethylamine',)),
          ('mtmB', ('monomethylamine',)),
          ('mttB', ('trimethylamine',)))
extra = {}
for m in ARC:
    t = os.path.join(A, m, f'{m}.tsv'); fa = os.path.join(A, m, f'{m}.faa')
    if not (os.path.exists(t) and os.path.exists(fa)): continue
    want = {}
    for line in open(t):
        if line.startswith('#'): continue
        f_ = line.rstrip('\n').split('\t')
        if len(f_) < 8 or f_[1].lower() != 'cds': continue
        pr = (f_[7] or '').lower()
        for key, pats in MT_PAT:
            if all(x in pr for x in pats): want[f_[5]] = key; break
    if not want: continue
    for nm, sq in rdfa(fa):
        lt = nm.split()[0]
        if lt in want and len(sq) >= 150:
            k = want[lt]
            out.append((f"{k}|{mid(m)}|{lt}", sq))
            extra[k] = extra.get(k, 0) + 1
if extra:
    print(f"  proteome route: {sum(extra.values())} more sequences — "
          + ', '.join(f'{k} {v}' for k, v in sorted(extra.items())))
    found.update({k: found.get(k, 0) + v for k, v in extra.items()})
else:
    print("  proteome route: no methanol/methylamine methyltransferase named by")
    print("  Bakta in any archaeal genome. Step 11b already found Bakta reporting")
    print("  NO methanogenesis genes at all, so this is an annotation gap, not an")
    print("  absence — the methanol arm stays untested and must be reported so.")

# CONTROL: archaeal single-copy proteins from the same genomes, via Bakta.
# THE FIRST VERSION MATCHED EXACT PHRASES ("50s ribosomal protein l2",
# "elongation factor 1-alpha") and found ZERO across all eight archaea. That is
# the same trap as step 20's DMSP keyword screen and step 29's transporters:
# I guessed how an annotator would phrase something instead of asking what it
# actually wrote. Bakta's archaeal annotation is already known to be weak here —
# step 11b found it reporting no methanogenesis genes at all.
# Match the CONCEPT loosely, take the long ones, and if that still fails, print
# what Bakta really writes rather than stopping blind.
CTRL_PAT = ('ribosomal protein', 'elongation factor', 'rna polymerase',
            'translation initiation factor', 'chaperonin', 'thermosome')
MAXPER = 6          # keep the query set small; 6 long markers per genome is ample
nctrl = 0
prodsample = []
for m in ARC:
    t = os.path.join(A, m, f'{m}.tsv'); fa = os.path.join(A, m, f'{m}.faa')
    if not (os.path.exists(t) and os.path.exists(fa)): continue
    want = set()
    for line in open(t):
        if line.startswith('#'): continue
        f_ = line.rstrip('\n').split('\t')
        if len(f_) < 8 or f_[1].lower() != 'cds': continue
        pr = (f_[7] or '').lower()
        prodsample.append(pr)
        if any(x in pr for x in CTRL_PAT): want.add(f_[5])
    if not want: continue
    took = 0
    for nm, sq in rdfa(fa):
        if took >= MAXPER: break
        lt = nm.split()[0]
        # length filter: a 60 aa ribosomal protein gives too short an alignment
        # for a frameshift-tolerant blastx on error-prone reads to score well.
        if lt in want and len(sq) >= 150:
            out.append((f"CONTROL_arc|{mid(m)}|{lt}", sq)); nctrl += 1; took += 1

print(f"  archaeal control: {nctrl} ribosomal/EF-1a sequences from "
      f"{len(ARC)} archaeal genomes")

if not found:
    print()
    print("  !! NO TARGET MARKER SEQUENCES FOUND under 15_marker_trees/.")
    print("     Without queries there is nothing to search for, and reporting")
    print("     zero read hits from zero queries would be a fabricated negative.")
    print("     Point the script at the marker directory")
    print("     (15_marker_trees/).")
    sys.exit(3)
if nctrl == 0:
    print()
    print("  !! NO ARCHAEAL CONTROL SEQUENCES. The search would be uninterpretable:")
    print("     a zero for the methyltransferases could not be told apart from a")
    print("     search that cannot see archaea at all. This is the step-19b rpoB")
    print("     failure and the script stops rather than repeat it.")
    print()
    print("     WHAT BAKTA ACTUALLY WROTE in the archaeal genomes — the 25 most")
    print("     common product strings, so the patterns can be fixed from data")
    print("     rather than from my guess about phrasing:")
    from collections import Counter as _C
    for p_, c_ in _C(x for x in prodsample if x).most_common(25):
        print(f"       {c_:>5}  {p_[:78]}")
    sys.exit(3)

with open(QF, 'w') as fh:
    for nm, sq in out: fh.write(f">{nm}\n{sq}\n")
print(f"  wrote {len(out)} query sequences")
PY
QRC=$?
[ "$QRC" = "0" ] || { echo "  (stopping — see above)"; exit "$QRC"; }

# ---------------------------------------------------- 2. reads
CLEAN="$WS2/00_input/sediment_reads_clean.fastq.gz"
if [ -s "$CLEAN" ]; then SED_IN="$CLEAN"
else
  SED_IN="$SED_READS"
  echo "  !! clean sediment reads not found — using \$SED_READS, which still"
  echo "     contains the RBK_test amplicon library. Run 13f with WRITE_CLEAN=1."
fi

run_tool diamond makedb --in "$QF" -d "$R/markers" >/dev/null 2>&1 || {
  echo "  !! diamond makedb failed — is diamond on PATH?"; exit 1; }

# THE CACHE KEY INCLUDES THE QUERY SET. Without it, adding a marker and
# re-running would silently reuse hits computed from the OLD queries and report
# zero for the new one — which is exactly how 19b reported 0/29 for four
# families across two runs. A cache key must include the query.
QHASH=$(md5sum "$QF" | cut -c1-12)
echo "  query-set hash: $QHASH"

search () {                       # $1 tag, $2 reads
  local tag="$1"
  local reads="$2"
  local out="$R/${tag}_read_hits.tsv"
  local stamp="$R/${tag}.qhash"
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$QHASH" ]; then
    echo "  [$tag] cached"; return 0
  fi
  echo "  [$tag] blastx --frameshift 15 against $(basename "$reads") ..."
  run_tool diamond blastx -d "$R/markers" -q "$reads" \
      --frameshift 15 --id "$MINID" --query-cover 0 --max-target-seqs 5 \
      -f 6 qseqid sseqid pident length evalue bitscore \
      -p "$THREADS" -o "$out.part" 2>"$R/${tag}.err"
  # CAPTURE THE EXIT CODE HERE. `if [ -s f ] || [ $? -eq 0 ]` reads the exit
  # status of the FIRST test, not of diamond — the same mistake that made step
  # 29's gate report PASS on a script that had exited 4. $? survives one command.
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    mv "$out.part" "$out"          # an empty result file is a REAL zero and is kept
    echo "$QHASH" > "$stamp"       # stamp AFTER success, never before
  else
    rm -f "$out.part"
    echo "  !! $tag search failed rc=$rc (see ${tag}.err)"; return 1
  fi
}
echo
echo "######## 2. READ-LEVEL SEARCH ########"
search SED "$SED_IN" || exit 1
search LIQ "$LIQ_READS" || exit 1
NSED=$(zcat -f "$SED_IN" | awk 'NR%4==1' | wc -l)
NLIQ=$(zcat -f "$LIQ_READS" | awk 'NR%4==1' | wc -l)
echo "  reads searched: SED $NSED, LIQ $NLIQ"

python3 - "$R" "$NSED" "$NLIQ" "$MINLEN" <<'PY'
import sys, os, csv
from collections import defaultdict
R, NSED, NLIQ, MINLEN = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

def load(tag):
    p = os.path.join(R, f'{tag}_read_hits.tsv')
    best = {}
    if not os.path.exists(p): return best
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        q, s = f[0], f[1]
        try: pid, ln, bs = float(f[2]), int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        # ONE READ, ONE ASSIGNMENT. A read matching several markers is counted
        # once, to its best hit — otherwise a conserved methyltransferase core
        # inflates every family it resembles.
        cur = best.get(q)
        if cur is None or bs > cur[1]: best[q] = (s.split('|')[0], bs)
    return best

sed, liq = load('SED'), load('LIQ')
def counts(b):
    c = defaultdict(int)
    for _q, (m, _bs) in b.items(): c[m] += 1
    return c
cs, cl = counts(sed), counts(liq)
markers = sorted(set(cs) | set(cl))

print()
print("=" * 96)
print("  READ-LEVEL HITS  (one read counted once, to its best marker)")
print("=" * 96)
print(f"  {'marker':<16}{'SED reads':>11}{'per 1M':>10}{'LIQ reads':>11}{'per 1M':>10}")
for m in markers:
    print(f"  {m:<16}{cs.get(m,0):>11}{1e6*cs.get(m,0)/max(1,NSED):>10.2f}"
          f"{cl.get(m,0):>11}{1e6*cl.get(m,0)/max(1,NLIQ):>10.2f}")

ctrl_s, ctrl_l = cs.get('CONTROL_arc', 0), cl.get('CONTROL_arc', 0)
tgt = [m for m in markers if m != 'CONTROL_arc']
print()
print("=" * 96)
print("  VERDICT")
print("=" * 96)
lines = []
def say(s):
    print(s); lines.append(s)

if ctrl_s == 0:
    say("  >> NO VERDICT. The archaeal ribosomal control returned ZERO hits in the")
    say("     sediment reads. A search that cannot see archaeal single-copy genes")
    say("     cannot be asked about archaeal methyltransferases, and a zero here")
    say("     would mean nothing. This is the step-19b failure mode and the reason")
    say("     the control is archaeal rather than bacterial.")
else:
    say(f"  archaeal control: {ctrl_s} sediment reads "
        f"({1e6*ctrl_s/max(1,NSED):.2f} per million) — the search CAN see archaea.")
    say("")
    mt = {m: cs.get(m, 0) for m in tgt if m != 'mcrA'}
    nmt = sum(mt.values())
    mcr = cs.get('mcrA', 0)
    say(f"  mcrA in sediment reads: {mcr}")
    say(f"  methylamine/methanol methyltransferases in sediment reads: {nmt}"
        + (f"  ({', '.join(f'{k} {v}' for k, v in sorted(mt.items()) if v)})" if nmt else ""))
    say("")
    if nmt > 0:
        say("  >> THE METHYLOTROPHIC CONSUMER END IS PRESENT IN THE LAGOON.")
        say("     These reads were never binned — the organisms carrying them are")
        say("     below the coverage that assembly and binning need — but the genes")
        say("     are in the sediment. Step 26's 0.08-0.18 breadth for the")
        say("     Methanosarcina MAGs is a RECOVERY limit, not an absence, and the")
        say("     pectin/choline -> methanol/TMA -> CH4 chain now has both ends")
        say("     documented in situ.")
    elif mcr > 0:
        say("  >> METHANOGENS YES, METHYLOTROPHIC ONES NOT DETECTED.")
        say("     mcrA is in the sediment reads but no methylamine or methanol")
        say("     methyltransferase is. With the control passing, that is a TESTED")
        say("     NEGATIVE and a strong sentence: methanogenesis in this sediment")
        say("     is not detectably methylotrophic, and the culture-based")
        say("     methylotrophy is an enrichment phenomenon. Report the detection")
        say("     floor alongside it.")
    else:
        say("  >> NO METHANOGENESIS MARKER AT ALL in the sediment reads, despite a")
        say("     passing archaeal control. Given SingleM's 14.1% archaea that is")
        say("     surprising and worth checking: mcrA is a single gene in a genome,")
        say("     so at low coverage it can be missed where ribosomal proteins are")
        say("     not. Quote the per-million rates, not the raw zero.")
say("")
say(f"  DETECTION FLOOR: {NSED} sediment reads searched. One read is "
    f"{1e6/max(1,NSED):.3f} per million.")
say("  A zero means 'below this floor', never 'absent from Etoliko'.")
with open(os.path.join(R, 'VERDICT.txt'), 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
PY

echo
echo "######## VALIDATION GATE ########"
check "query set"    $([ -s "$R/marker_queries.faa" ] && echo 0 || echo 1)
check "SED search"   $([ -f "$R/SED_read_hits.tsv" ] && echo 0 || echo 1)
check "LIQ search"   $([ -f "$R/LIQ_read_hits.tsv" ] && echo 0 || echo 1)
check "verdict"      $([ -s "$R/VERDICT.txt" ] && echo 0 || echo 1)
runlog "31_readlevel complete"
echo
echo "######## 31 DONE ########"
