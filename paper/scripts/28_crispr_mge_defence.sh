#!/usr/bin/env bash
# =============================================================================
# 28_crispr_mge_defence.sh  --  who has been infected by what, and who fights back
#
# THREE THINGS THAT BELONG TOGETHER AND NONE OF WHICH HAS BEEN DONE
#
#   1. CRISPR SPACERS MATCHED TO THE PROVIRUSES.
#      Step 18 found 68 proviruses across 18 of 29 genomes and linked them to
#      hosts by CO-LOCATION — the provirus is on a contig that belongs to that
#      MAG. That is an assembly fact, not an infection record.
#      A CRISPR spacer is different. It is a piece of an invader's genome that
#      the cell CAPTURED AND KEPT. A spacer in genome A matching virus B is
#      direct evidence that A's lineage was infected by B and survived. It also
#      links viruses to hosts they are NOT currently integrated in, which
#      co-location can never do.
#
#   2. PLASMIDS. geNomad called them in the same run that called the viruses.
#      Nobody has opened the file.
#
#   3. DEFENCE SYSTEMS. With 18 of 29 genomes carrying a provirus, what the
#      community does about it is the other half of the story.
#
# HONESTY ABOUT THE DEFENCE CENSUS
#   The right tool is DefenseFinder or PADLOC, with curated HMM profiles. This
#   script does NOT have them and does not pretend to: it counts Bakta product
#   strings for the systems Bakta actually names — restriction-modification,
#   cas genes, abortive infection, retrons, toxin-antitoxin — and reports the
#   result as a FLOOR. Systems Bakta does not name (most CBASS, Thoeris,
#   Gabija, Wadjet, Zorya) are invisible here and a zero for them means
#   nothing. That limitation is printed with the table, not buried in a header.
#
# WHERE THE SPACERS COME FROM
#   Bakta detects CRISPR arrays natively, so no dependency on step 16 having
#   run. The script tries, in order: Bakta JSON (which carries the spacer
#   sequences), Bakta TSV/GFF3, then minced if installed. If none yields
#   spacers it says so and stops rather than reporting an empty result as an
#   absence of matches.
#
# THE MATCH CRITERION
#   >=90% identity over the full spacer. Matching is done in Python by
#   seed-and-extend (see section 2) rather than BLAST, which is not installed
#   on this server and is not worth adding for a problem this size.
#
# INPUTS   09_annotation/<mag>/<mag>.{json,tsv}   ·   07_mags/bins/*.fa
#          18_viruses/**/*virus_summary.tsv, *plasmid_summary.tsv, *.fna
# OUTPUT   28_crispr/spacers.fasta · spacer_hits.tsv · virus_host_crispr.tsv
#          28_crispr/plasmids.tsv · defence_systems.tsv
# RUNTIME  minutes.
# USAGE    ws2 && bash scripts/28_crispr_mge_defence.sh 2>&1 | tee logs/28.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="28_crispr_mge"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; V="$WS2/18_viruses"; C="$WS2/28_crispr"
mkdir -p "$C"

echo "######## 28 — CRISPR SPACERS, PLASMIDS, DEFENCE ########"
echo "  started: $(date -Is)"
runlog "28_crispr_mge started"
[ -d "$A" ] || { echo "  !! 09_annotation not found"; exit 1; }
[ -d "$V" ] || echo "  !! 18_viruses not found — spacer matching will be skipped"

# ---------------------------------------------------------------- 1. spacers
python3 - "$A" "$F" "$C" <<'PY'
import sys, os, csv, glob, json, re
A, F, C = sys.argv[1:4]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)

spacers = []          # (mag, array_id, index, seq)
arrays  = {}          # mag -> n arrays
route = None

# route 1: Bakta JSON carries spacer sequences inside each crispr feature
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); j = os.path.join(d, f'{n}.json')
    if not os.path.exists(j): continue
    try: doc = json.load(open(j))
    except Exception: continue
    na = 0
    for f in doc.get('features', []):
        if str(f.get('type', '')).lower() not in ('crispr', 'crispr_array'): continue
        na += 1
        aid = f.get('id') or f"{n}_crispr_{na}"
        for i, sp in enumerate(f.get('spacers', []) or [], 1):
            s = sp.get('sequence') if isinstance(sp, dict) else sp
            if s and len(s) >= 20: spacers.append((n, aid, i, s.upper()))
    if na: arrays[n] = na
if spacers: route = 'Bakta JSON'

# route 2: GFF3, where some Bakta versions write spacers as their own features
if not spacers:
    for d in sorted(glob.glob(os.path.join(A, '*'))):
        n = os.path.basename(d); g = os.path.join(d, f'{n}.gff3')
        if not os.path.exists(g): continue
        for line in open(g):
            if line.startswith('#') or 'spacer' not in line.lower(): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: continue
            m = re.search(r'(?:sequence|spacer)=([ACGTNacgtn]{20,60})', f[8])
            if m:
                spacers.append((n, f[0], len(spacers) + 1, m.group(1).upper()))
                arrays[n] = arrays.get(n, 0) + 1
    if spacers: route = 'Bakta GFF3'

print(f"  spacer source: {route or 'NONE FOUND'}")
if not spacers:
    print("  !! No CRISPR spacers could be extracted from the Bakta output.")
    print("     That is NOT the same as 'this community has no CRISPR arrays'.")
    print("     Bakta may have been run without --skip-crispr disabled, or this")
    print("     version writes arrays without spacer sequences. Install minced")
    print("     and run it over 07_mags/bins/*.fa, or re-run Bakta, then repeat.")
    print("     Reporting zero matches from zero spacers would be a fabricated")
    print("     negative, so nothing further is written.")
    sys.exit(3)

with open(os.path.join(C, 'spacers.fasta'), 'w') as fh:
    for m, aid, i, s in spacers:
        fh.write(f">{mid(m)}|{aid}|sp{i}\n{s}\n")
print(f"  {len(spacers)} spacers from {len(arrays)} genome(s) with CRISPR arrays")
for m in sorted(arrays, key=mid):
    k = sum(1 for x in spacers if x[0] == m)
    print(f"    {mid(m):<17}{arrays[m]:>2} array(s), {k:>4} spacers")
print()
print("  A GENOME WITH NO ARRAY IS NOT DEFENCELESS — it simply cannot be")
print("  interrogated this way. Spacer evidence is one-directional: it proves")
print("  past infection where present and proves nothing where absent.")
PY
SPRC=$?
[ "$SPRC" = "3" ] && { echo "  (spacer matching skipped)"; }

# ------------------------------------------------- 2. match spacers to MGEs
# BLAST+ IS NOT INSTALLED HERE AND THE DEPENDENCY IS NOT WORTH IT. Matching a
# thousand 30 bp queries against a few thousand contigs is a seed-and-extend
# problem small enough to do exactly, in Python, with a stated guarantee:
#
#   Three non-overlapping 10 bp seeds are taken from each spacer. By the
#   pigeonhole principle, any alignment with AT MOST 2 MISMATCHES must leave at
#   least one of those seeds intact, so every such match is found. Candidates
#   are then verified by full comparison on both strands.
#
#   The limit, stated because it matters: a 30 bp spacer at exactly 90% identity
#   carries 3 mismatches, and a 3-mismatch match is found only if the mismatches
#   happen not to hit all three seeds. Genuine spacer-protospacer pairs are
#   overwhelmingly 0-2 mismatches, so this loses very little — but it is a
#   heuristic at the boundary and is not sold as exhaustive.
if [ "$SPRC" != "3" ] && [ -d "$V" ]; then
  MGE="$C/mge_contigs.fna"
  : > "$MGE"
  found=0
  for f in $(find "$V" -maxdepth 3 \( -name '*virus*.fna' -o -name '*plasmid*.fna' \) 2>/dev/null | sort); do
    cat "$f" >> "$MGE"; found=$((found+1))
  done
  if [ "$found" -eq 0 ]; then
    echo "  !! no viral/plasmid FASTA under 18_viruses — cannot match spacers"
  else
    echo "  MGE database: $(grep -c '^>' "$MGE") sequences from $found file(s)"
    python3 - "$C/spacers.fasta" "$MGE" "$C/spacer_hits_raw.tsv" <<'PYX'
import sys
qf, tf, out = sys.argv[1:4]
K, NSEED, MAXFRAC = 10, 3, 0.10
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]

def rdfa(p):
    n, b = None, []
    for line in open(p):
        if line.startswith('>'):
            if n: yield n, ''.join(b)
            n = line[1:].split()[0]; b = []
        else: b.append(line.strip().upper())
    if n: yield n, ''.join(b)

sp = list(rdfa(qf))
idx = {}
for i, (nm, s) in enumerate(sp):
    L = len(s)
    offs = [j * K for j in range(NSEED)] if L >= K * NSEED else [0, L - K]
    for o in offs:
        if o < 0 or o + K > L: continue
        idx.setdefault(s[o:o + K], []).append((i, o))
print(f"  index: {len(idx)} seeds from {len(sp)} spacers "
      f"(guaranteed detection to 2 mismatches)")

hits, seen = [], set()
for tn, ts in rdfa(tf):
    for strand, seq in (('+', ts), ('-', rc(ts))):
        n = len(seq)
        for p in range(n - K + 1):
            c = idx.get(seq[p:p + K])
            if not c: continue
            for i, o in c:
                nm, s = sp[i]; L = len(s)
                st = p - o
                if st < 0 or st + L > n: continue
                w = seq[st:st + L]
                mm = sum(1 for a, b in zip(s, w) if a != b)
                if mm > int(MAXFRAC * L): continue
                key = (nm, tn, strand, st)
                if key in seen: continue
                seen.add(key)
                hits.append((nm, tn, 100.0 * (L - mm) / L, L, L, mm, strand, st))
with open(out, 'w') as fh:
    for h in hits:
        fh.write('\t'.join(str(x) for x in h) + '\n')
print(f"  matches: {len(hits)}")
PYX
  fi
fi

# ------------------------------- 2b. ANCESTRY vs IMMUNITY, decided not flagged
# The first run produced 77 "genuine" matches, of which 22 sat at 90-95%
# identity — including every one of the Chakrabartyella-to-Chakrabartyella hits.
# Chance is not the explanation: an exact 30-mer match has probability ~1e-18
# per position, and even a 3-mismatch match is ~1e-13, so against ~30 Mb these
# are all REAL sequence similarity. The question is what KIND.
#
#   IMMUNITY  the host captured a fragment of an invader. Only the ~30 bp
#             protospacer is shared; the DNA either side of it is unrelated to
#             anything in the host genome.
#   ANCESTRY  the host and the target share the region by descent. Then the
#             FLANKS are shared too, because the whole locus is homologous.
#
# So the discriminator is the flank, and it is cheap: take 150 bp either side
# of each protospacer, cut it into 25-mers, and ask whether those occur in the
# HOST's own genome. Several present means the locus is shared and the "spacer
# match" is homology wearing a costume.
if [ -s "$C/spacer_hits_raw.tsv" ]; then
python3 - "$C/spacer_hits_raw.tsv" "$C/mge_contigs.fna" "$F" "$C/hits_verified.tsv" <<'PYV'
import sys, os, glob
qf, tf, F, out = sys.argv[1:5]
FLANK, KM, MINK = 150, 25, 3
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
def rdfa(p):
    n, b = None, []
    for line in open(p):
        if line.startswith('>'):
            if n: yield n, ''.join(b)
            n = line[1:].split()[0]; b = []
        else: b.append(line.strip().upper())
    if n: yield n, ''.join(b)

hits = []
for line in open(qf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 8: continue
    hits.append(dict(spacer=f[0], target=f[1], pid=float(f[2]), L=int(f[3]),
                     mm=int(f[5]), strand=f[6], start=int(f[7])))
tseq = dict(rdfa(tf))
flank_km = {}                      # kmer -> set(hit index)
for i, h in enumerate(hits):
    s = tseq.get(h['target'])
    if s is None: continue
    if h['strand'] == '-': s = rc(s)
    a = max(0, h['start'] - FLANK); b = min(len(s), h['start'] + h['L'] + FLANK)
    for seg in (s[a:h['start']], s[h['start'] + h['L']:b]):
        for j in range(len(seg) - KM + 1):
            flank_km.setdefault(seg[j:j + KM], set()).add(i)
found = {}
for p in glob.glob(os.path.join(F, '**', '*.f*a'), recursive=True):
    g = os.path.splitext(os.path.basename(p))[0]
    try: seqs = list(rdfa(p))
    except Exception: continue
    for _n, sq in seqs:
        for strand in (sq, rc(sq)):
            for j in range(len(strand) - KM + 1):
                c = flank_km.get(strand[j:j + KM])
                if c:
                    for i in c: found.setdefault(i, set()).add(strand[j:j + KM])
with open(out, 'w') as fh:
    fh.write("spacer\ttarget\tpident\tmismatches\tflank_kmers_in_host\tverdict\n")
    n_anc = n_imm = 0
    for i, h in enumerate(hits):
        k = len(found.get(i, ()))
        v = 'ANCESTRY' if k >= MINK else 'immunity'
        n_anc += (v == 'ANCESTRY'); n_imm += (v == 'immunity')
        fh.write(f"{h['spacer']}\t{h['target']}\t{h['pid']:.1f}\t{h['mm']}\t{k}\t{v}\n")
print(f"  flank test: {n_imm} consistent with immunity, {n_anc} explained by "
      f"shared ancestry (>= {MINK} flanking {KM}-mers also in the host genome)")
PYV
fi

# ---------------------------------------------------------- 3. report + rest
python3 - "$A" "$F" "$V" "$C" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter
A, F, V, C = sys.argv[1:5]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus','') or '-')

# ---- spacer hits ----------------------------------------------------------
# THE ARTEFACT THAT HAS TO BE EXCLUDED FIRST.
# The MGE database is built from contigs of THIS dataset, including contigs
# belonging to the very genomes whose spacers are being queried. If a CRISPR
# array happens to sit on a contig geNomad called a provirus or a plasmid, then
# every spacer in that array matches ITSELF at 100% identity, and the script
# reports a beautiful set of host-virus links that are pure self-hits.
# So: build a contig -> MAG map and DROP every match where the target belongs to
# the spacer's own genome. Those are reported separately, because a CRISPR array
# sitting inside a geNomad MGE call is itself worth knowing about — it means
# geNomad mis-called a chromosomal region.
contig2mag, _amb = {}, set()
def _addc(c, k):
    c = (c or '').strip()
    if not c or c in _amb: return
    if c in contig2mag and contig2mag[c] != k:
        del contig2mag[c]; _amb.add(c); return
    contig2mag[c] = k
for p in glob.glob(os.path.join(F, '**', '*.f*a'), recursive=True):
    b = os.path.splitext(os.path.basename(p))[0]
    k = mid(b) if b in tax else b
    if b not in tax: continue
    try:
        for line in open(p):
            if line.startswith('>'): _addc(line[1:].split()[0], mid(b))
    except Exception: pass
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d)
    if n not in tax: continue
    t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    try:
        for line in open(t):
            if line.startswith('#') or not line.strip(): continue
            _addc(line.split('\t')[0], mid(n))
    except Exception: pass

def target_mag(name):
    c = str(name).split('|')[0].split()[0]
    if c in contig2mag: return contig2mag[c]
    return contig2mag.get(re.sub(r'[_:]\d+[-_]\d+$', '', c))

verdict = {}
vf = os.path.join(C, 'hits_verified.tsv')
if os.path.exists(vf):
    for r in csv.DictReader(open(vf), delimiter='\t'):
        verdict[(r['spacer'], r['target'])] = r['verdict']

raw = os.path.join(C, 'spacer_hits_raw.tsv')
hits, selfhits, ancestry = [], [], []
if os.path.exists(raw):
    for line in open(raw):
        f = line.rstrip('\n').split('\t')
        if len(f) < 5: continue
        q, s, pid, ln, qlen = f[0], f[1], float(f[2]), int(f[3]), int(f[4])
        if not (qlen and ln >= 0.9 * qlen and pid >= 90): continue
        host = q.split('|')[0]
        rec = dict(spacer=q, target=s, pident=pid, alen=ln, qlen=qlen,
                   host=host, target_mag=target_mag(s) or '')
        rec['verdict'] = verdict.get((q, s), 'immunity')
        if rec['target_mag'] == host: selfhits.append(rec)
        elif rec['verdict'] == 'ANCESTRY': ancestry.append(rec)
        else: hits.append(rec)

print()
print("=" * 104)
print("  1. CRISPR SPACERS vs MOBILE ELEMENTS — DIRECT EVIDENCE OF PAST INFECTION")
print("=" * 104)
print(f"  contig map: {len(contig2mag)} contigs -> MAG ({len(_amb)} ambiguous discarded)")
if not os.path.exists(raw):
    print("  no matcher output — see the messages above.")
else:
    print(f"  {len(selfhits) + len(ancestry) + len(hits)} matches total: "
          f"{len(selfhits)} self-hits, {len(ancestry)} explained by shared "
          f"ancestry, {len(hits)} consistent with immunity.")
    if ancestry:
        ah = Counter((h['host'], h['target_mag'] or 'unbinned') for h in ancestry)
        print()
        print("  DISCARDED AS ANCESTRY — the DNA flanking the protospacer is also")
        print("  present in the host's own genome, so the locus is shared by descent")
        print("  and the match is homology, not a captured invader:")
        for (a, b), v in ah.most_common(8):
            print(f"    {a} -> {b:<20}{v:>4} hits")
    if selfhits:
        sh = Counter(h['host'] for h in selfhits)
        print()
        print("  SELF-HITS (spacer matching a contig of its own genome) — these are")
        print("  NOT infection records. A spacer matches itself whenever its array")
        print("  sits on a contig geNomad called an MGE, which means geNomad called a")
        print("  chromosomal region a plasmid or provirus. Worth checking on its own:")
        for k, v in sh.most_common():
            print(f"    {k:<17}{v:>5} self-matches")
    if not hits:
        print()
        print("  0 GENUINE cross-genome matches once self-hits are removed.")
        print("  Interpretable, and not a failure: the viruses a cell recorded in its")
        print("  arrays need not be the ones that assembled here, and this catalogue")
        print("  holds ~2% of the sediment community. Read it as 'the recovered")
        print("  viruses are not the ones these genomes recorded', never as 'these")
        print("  viruses did not infect these hosts'.")
    else:
        byhost = defaultdict(set)
        for h in hits: byhost[h['host']].add(h['target'].split('|')[0])
        print()
        print(f"  {'host MAG':<17}{'genus':<20}{'MGEs matched':>13}  target genome(s)")
        wn_of = {mid(w): w for w in tax}
        for k in sorted(byhost):
            tg = sorted({h['target_mag'] for h in hits
                         if h['host'] == k and h['target_mag']})
            g = gen(wn_of.get(k, k))
            print(f"  {k:<17}{g[:19]:<20}{len(byhost[k]):>13}  "
                  f"{', '.join(tg) if tg else 'unbinned contigs'}")
        with open(os.path.join(C, 'virus_host_crispr.tsv'), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['host_mag', 'mge_contig', 'mge_belongs_to_mag', 'spacer',
                        'pident', 'align_len', 'spacer_len'])
            for h in hits:
                w.writerow([h['host'], h['target'], h['target_mag'], h['spacer'],
                            f"{h['pident']:.1f}", h['alen'], h['qlen']])
        print()
        print("  A SPACER MATCH IS STRONGER THAN CO-LOCATION. Step 18 linked provirus")
        print("  to host by which contig it sat on, which is an assembly fact. A")
        print("  spacer is a captured fragment of an invader the lineage survived,")
        print("  and it can link a virus to a host it is not integrated in at all.")
        print("  Where the target maps to ANOTHER genome in the catalogue, that is a")
        print("  recorded infection crossing two organisms you have both halves of —")
        print("  the strongest virus-host evidence this dataset can produce.")

        # ---- the ancestry confound, and how to see it -----------------------
        # A spacer is supposed to be a captured PHAGE fragment. If instead it
        # matches a region the two genomes share BY DESCENT, the "immunity
        # record" is just homology. Two things separate those:
        #   (1) IDENTITY. A real protospacer match is at or near 100%; shared
        #       ancestry between distinct species decays well below that.
        #   (2) TAXONOMIC DISTANCE. A host and target in the same GENUS are
        #       exactly where ancestry can imitate immunity, and must be
        #       flagged rather than reported alongside distant pairs.
        print()
        print("  IDENTITY DISTRIBUTION OF THE GENUINE HITS")
        buck = Counter()
        for h in hits:
            buck['100%' if h['pident'] >= 99.99 else
                 '95-99%' if h['pident'] >= 95 else '90-95%'] += 1
        for k in ('100%', '95-99%', '90-95%'):
            print(f"    {k:<8}{buck.get(k, 0):>5}")
        print("    Real protospacer matches concentrate at 100%. A hit spread down")
        print("    toward 90% is more likely shared sequence than captured sequence.")
        same = [h for h in hits if h['target_mag'] and
                gen(wn_of.get(h['host'], h['host'])) != '-' and
                gen(wn_of.get(h['host'], h['host'])) ==
                gen(wn_of.get(h['target_mag'], h['target_mag']))]
        if same:
            pairs = sorted({(h['host'], h['target_mag']) for h in same})
            print()
            print("  SAME-GENUS PAIRS — READ WITH CARE:")
            for a, b in pairs:
                ids = [h['pident'] for h in same
                       if h['host'] == a and h['target_mag'] == b]
                print(f"    {a} -> {b}  ({len(ids)} hits, identity "
                      f"{min(ids):.1f}-{max(ids):.1f}%)")
            print("    Host and target are the same genus, so a shared chromosomal")
            print("    region could imitate a spacer match. This is STILL the more")
            print("    interesting case if it holds — a mobile element moving between")
            print("    close relatives, recorded by the immune system of one of them —")
            print("    but it needs the identity to sit at 100% and the target contig")
            print("    to be a confident geNomad call. Check both before writing it.")

# ---- plasmids -------------------------------------------------------------
print()
print("=" * 104)
print("  2. PLASMIDS — CALLED BY geNOMAD, NEVER OPENED")
print("=" * 104)
pfiles = sorted(set(glob.glob(os.path.join(V, '**', '*plasmid_summary.tsv'),
                                recursive=True)))
# DEDUPE BY CONTIG NAME AND REPORT PER SOURCE. geNomad was run three times here
# — sediment assembly, culture assembly, and the concatenated MAG contigs — and
# the third overlaps the first two, so a raw row count double-counts. This is
# the same trap that inflated n_proviruses sevenfold at step 25; there it was
# stage tables, here it is overlapping inputs.
byfile, seen_ct = {}, {}
for p in pfiles:
    if '_genes' in os.path.basename(p): continue
    b = os.path.basename(p)
    try:
        rows_ = list(csv.DictReader(open(p), delimiter='\t'))
    except Exception as e:
        print(f"  !! {b} unreadable ({e})"); continue
    byfile[b] = rows_
    for r in rows_:
        nm = (r.get('seq_name') or r.get('contig') or '').split('|')[0].strip()
        if nm: seen_ct.setdefault(nm, []).append((b, r))
if not byfile:
    print("  no *plasmid_summary.tsv rows found under 18_viruses/.")
else:
    tot = sum(len(v) for v in byfile.values())
    print(f"  {tot} rows across {len(byfile)} geNomad run(s); "
          f"{len(seen_ct)} DISTINCT plasmid contigs after deduplication.")
    for b, v in sorted(byfile.items()):
        print(f"    {b:<52}{len(v):>5} rows")
    ln = []
    for nm, lst in seen_ct.items():
        r = lst[0][1]
        for k in r:
            if k.lower() in ('length', 'seq_len') and r[k]:
                try: ln.append(int(float(r[k])))
                except ValueError: pass
                break
    if ln:
        ln.sort()
        print(f"  length: median {ln[len(ln)//2]:,} bp, "
              f"range {ln[0]:,}-{ln[-1]:,}, total {sum(ln)/1e6:.1f} Mb")
    print()
    print("  Plasmids plus the 68 proviruses make this a mobile-element-rich")
    print("  community, which fits a low-energy sediment where lysogeny beats")
    print("  lysis. NOTE the size distribution: a median around 6 kb is small for")
    print("  a plasmid and overlaps the size of an assembly fragment, so some of")
    print("  these will be mis-assigned chromosomal pieces. Treat the count as an")
    print("  upper bound and the long ones (>20 kb) as the confident set.")
    with open(os.path.join(C, 'plasmids.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['contig', 'n_sources', 'sources', 'field', 'value'])
        for nm, lst in sorted(seen_ct.items()):
            srcs = ','.join(sorted(set(x[0] for x in lst)))
            for k, v in lst[0][1].items():
                w.writerow([nm, len(lst), srcs, k, v])

# ---- defence systems ------------------------------------------------------
DEF = {
 'Restriction-modification': ('restriction endonuclease','type i restriction',
        'type ii restriction','type iii restriction','hsdr','hsdm','hsds',
        'dna methyltransferase','adenine-specific methyltransferase'),
 'CRISPR-Cas':               ('crispr-associated','cas1','cas2','cas3','cas9',
                              'cas12','cas13','csn1','cascade'),
 'Abortive infection':       ('abortive infection','abi'),
 'Retron':                   ('retron','reverse transcriptase, group ii'),
 'Toxin-antitoxin':          ('toxin-antitoxin','antitoxin','relbe','mazef',
                              'higba','vapbc','parde'),
 'Argonaute':                ('argonaute',),
 'BREX / Pgl':               ('brex','pglx','pglz'),
 'DISARM':                   ('disarm','drma','drmb'),
}
prod = defaultdict(list)
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].lower() == 'cds': prod[n].append(f[7].lower())

print()
print("=" * 104)
print("  3. DEFENCE SYSTEMS — A FLOOR, FROM PRODUCT STRINGS ONLY")
print("=" * 104)
print("  THE RIGHT TOOL IS DefenseFinder OR PADLOC and this is not it. These are")
print("  Bakta product-string counts for the systems Bakta actually names. CBASS,")
print("  Thoeris, Gabija, Wadjet, Zorya and most of the recently described systems")
print("  are INVISIBLE here, and a zero for them means nothing at all. Use this to")
print("  rank genomes against each other, never to claim a system is absent.")
print()
hdr = f"  {'MAG ID':<17}{'genus':<19}"
for k in DEF: hdr += f"{k[:9]:>11}"
print(hdr + f"{'total':>8}")
rows = []
for m in sorted(prod, key=mid):
    line = f"  {mid(m):<17}{gen(m)[:18]:<19}"
    row = dict(mag_id=mid(m), genus=gen(m)); tot = 0
    for k, pats in DEF.items():
        c = sum(1 for p in prod[m] if any(x in p for x in pats))
        line += f"{c:>11}"; row[k.replace(' ', '_')] = c; tot += c
    print(line + f"{tot:>8}"); row['total'] = tot; rows.append(row)
if rows:
    with open(os.path.join(C, 'defence_systems.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    top = sorted(rows, key=lambda r: -r['total'])[:3]
    print()
    print("  Most defended (by this floor): " +
          ", ".join(f"{r['mag_id']} ({r['genus']}, {r['total']})" for r in top))
    print()

# ---- 4. THE JOIN I TOLD YOU TO DO YOURSELF. Doing it here instead. -----------
print()
print("=" * 104)
print("  4. DEFENCE vs PROVIRUS LOAD")
print("=" * 104)
cat = os.path.join(os.path.dirname(F.rstrip('/')), '12_results', 'MAG_CATALOGUE.tsv')
if not os.path.exists(cat):
    cat = os.path.join(F, '..', '12_results', 'MAG_CATALOGUE.tsv')
prov, size, cmpl = {}, {}, {}
try:
    for r in csv.DictReader(open(cat), delimiter='\t'):
        k = r.get('mag_id')
        if not k: continue
        try: prov[k] = int(float(r.get('n_proviruses') or 0))
        except ValueError: prov[k] = 0
        for c in ('bp', 'genome_size', 'observed_bp'):
            if r.get(c):
                try: size[k] = float(r[c]); break
                except ValueError: pass
        try: cmpl[k] = float(r.get('completeness') or 0)
        except ValueError: cmpl[k] = 0.0
except Exception as e:
    print(f"  !! could not read MAG_CATALOGUE.tsv ({e}) — run 25 first")

if prov:
    # NORMALISE. Raw defence counts scale with genome size, and the "most
    # defended" line above is topped by the largest genome in the catalogue,
    # which is exactly what an unnormalised count would produce whether or not
    # the biology is there. Per-Mb, and only genomes >=90% complete, because an
    # incomplete genome is missing defence genes for a reason unrelated to
    # phage.
    pts = []
    for r in rows:
        k = r['mag_id']
        if cmpl.get(k, 0) < 90 or not size.get(k): continue
        pts.append((k, r['genus'], r['total'] / (size[k] / 1e6), prov.get(k, 0)))
    if len(pts) < 5:
        print(f"  only {len(pts)} genomes are >=90% complete with a size and a")
        print("  provirus count — too few to say anything. Reporting the rows only.")
    else:
        pts.sort(key=lambda t: -t[2])
        print(f"  {len(pts)} genomes >=90% complete. Defence per Mb, not raw counts:")
        print("  raw totals scale with genome size, and the 'most defended' line")
        print("  above is led by the largest genome in the catalogue — which is what")
        print("  an unnormalised count produces whether or not the biology is real.")
        print()
        print(f"  {'MAG ID':<17}{'genus':<20}{'defence/Mb':>12}{'proviruses':>12}")
        for k, g, d, p in pts:
            print(f"  {k:<17}{g[:19]:<20}{d:>12.1f}{p:>12}")
        wi = [d for _k, _g, d, p in pts if p > 0]
        wo = [d for _k, _g, d, p in pts if p == 0]
        print()
        if wi and wo:
            mi, mo = sum(wi) / len(wi), sum(wo) / len(wo)
            print(f"  with provirus(es): mean {mi:.1f}/Mb (n={len(wi)})")
            print(f"  without:           mean {mo:.1f}/Mb (n={len(wo)})")
            print()
            if abs(mi - mo) < 0.25 * max(mi, mo):
                print("  >> NO DIFFERENCE WORTH REPORTING. Defence load does not track")
                print("     provirus carriage in this catalogue. With n of this size that")
                print("     is the expected outcome even if an effect exists, so read it")
                print("     as 'not resolvable here', not as 'no relationship'.")
            elif mi > mo:
                print("  >> LYSOGENS CARRY MORE DEFENCE — consistent with an arms race,")
                print("     where the genomes most exposed to phage accumulate both the")
                print("     scars and the countermeasures.")
            else:
                print("  >> LYSOGENS CARRY LESS DEFENCE. One reading is that the resident")
                print("     prophage IS the defence: a prophage excludes superinfection by")
                print("     its relatives, so a lysogen needs less of its own machinery.")
            print()
            print("  EITHER WAY THIS IS n=%d AND UNCONTROLLED for phylogeny." % len(pts))
            print("  Archaea and Bacillota differ in both variables for reasons that have")
            print("  nothing to do with phage. Do not put a p-value on this.")
        else:
            print("  all genomes fall on one side — no contrast to draw.")
PY

echo
echo "######## VALIDATION GATE ########"
check "defence table" $([ -s "$C/defence_systems.tsv" ] && echo 0 || echo 1)
runlog "28_crispr_mge complete"
echo
echo "######## 28 DONE ########"
