#!/usr/bin/env bash
# =============================================================================
# 23b_pul_analysis.sh  --  is Gaoshiqia a public good or a private digester?
#
# WHY THIS REPLACES THE SECRETOME ATTEMPT
#   DeepSig is unusable: it pins tensorflow 2.2.0, which no longer exists for
#   linux-64 on any channel. The package is unmaintained and no solver flag will
#   fix it. Rather than chase a signal-peptide predictor, ask the question the
#   secretome was only a proxy for — and answer it from the annotation already
#   on disk.
#
# THE QUESTION
#   Step 20 gave ETO_LIQ_MAG_002 (Gaoshiqia) 221 glycoside hydrolases, 12
#   polysaccharide lyases and 120 SusC/SusD proteins — an order of magnitude
#   more CAZymes than anything else in the catalogue. The chapter wants to call
#   it the base of the food chain: it breaks down polysaccharide, everyone else
#   eats the products, the fermenters feed the methanogens.
#
#   THAT MAY BE BACKWARDS. SusC is a TonB-dependent transporter and SusD a
#   surface lipoprotein. Together they form a POLYSACCHARIDE UTILISATION LOCUS
#   (PUL) — the Bacteroidota strategy of binding polymer at the cell surface,
#   cleaving it there, and importing the oligosaccharides into the periplasm
#   before anything else can reach them. It is a SELFISH uptake system. The
#   literature calls it exactly that, and it is the reason Bacteroidota dominate
#   polysaccharide-rich environments without feeding their neighbours.
#
#   So the discriminating question is not "how many CAZymes" but WHERE THEY SIT:
#     GH inside a PUL (next to susC/susD)  -> selfish, periplasmic, private
#     GH standing alone, far from any PUL  -> candidate secreted public good
#
# HOW THIS DECIDES IT WITHOUT A NEW TOOL
#   PULs are defined by GENOMIC CONTEXT: a susC/susD pair with CAZymes beside
#   them. Bakta already gives locus tags in genomic order with products, so the
#   loci can be reconstructed directly. This is how PULs are annotated in the
#   literature (PULDB/dbCAN-PUL use the same adjacency logic).
#
#   The output is the fraction of each genome's CAZymes that sit inside a PUL.
#   High fraction = private digester. Low fraction = the enzymes are dispersed
#   and more likely exported.
#
# WHAT THIS CANNOT DO
#   It cannot prove secretion. A GH outside a PUL is a CANDIDATE public good,
#   not a demonstrated one — that needs a signal-peptide prediction (SignalP 6.0
#   is free for academic use but requires a manual download from DTU) or
#   proteomics. Say "PUL-associated" and "non-PUL", never "secreted".
#
# OUTPUT   23_extras/pul/
#     pul_summary.tsv     per MAG: CAZymes, PULs, % of CAZymes inside a PUL
#     pul_loci.tsv        every PUL, its span, and the CAZymes it contains
#     nonPUL_cazymes.tsv  the candidate public goods
#
# RUNTIME  2-5 min — it only reads Bakta tables.
# USAGE    ws2 && bash scripts/23b_pul_analysis.sh 2>&1 | tee logs/23b_pul.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="23b_pul"
A="$WS2/09_annotation"; F="$WS2/07_mags"; P="$WS2/23_extras/pul"; mkdir -p "$P"

echo "######## 23b — POLYSACCHARIDE UTILISATION LOCI ########"
echo "  started: $(date -Is)"
runlog "23b_pul_analysis started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }

python3 - "$A" "$F" "$P" <<'PY'
import sys, os, csv, glob
from collections import Counter, defaultdict

A, F, P = sys.argv[1], sys.argv[2], sys.argv[3]
WINDOW = 10          # genes either side of a susC/susD pair
MIN_PAIR_GAP = 5     # susC and susD must be within this many genes of each other

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return (tax.get(m, {}) or {}).get('genus', '') or '-'

CAZY = ('glycoside hydrolase','glycosyl hydrolase','polysaccharide lyase',
        'carbohydrate esterase','glucanase','xylanase','chitinase','amylase',
        'pectinesterase','pectate lyase','cellulase','agarase','alginate lyase',
        'galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
        'arabinofuranosidase','levanase','pullulanase','chondroitin')
SUSC = ('susc','tonb-dependent receptor','tonb-dependent transporter')
SUSD = ('susd','ragb','sushi')

def classify(prod):
    p = (prod or '').lower()
    if any(x in p for x in SUSC): return 'susC'
    if any(x in p for x in SUSD): return 'susD'
    if any(x in p for x in CAZY): return 'CAZy'
    return ''

summary, loci, nonpul = [], [], []
print(f"\n  {'MAG ID':<17}{'genus':<20}{'CAZy':>6}{'susC':>6}{'susD':>6}"
      f"{'PULs':>6}{'inPUL':>9}{'%priv':>8}{'%priv-hi':>9}")
print("  %priv   = strict, needs a susC/susD PAIR")
print("  %priv-hi= upper bound, susC proximity alone (susD is poorly conserved")
print("            and often unannotated, so the strict rule undercounts)")
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    m = os.path.basename(d); t = os.path.join(d, f'{m}.tsv')
    if not os.path.exists(t): continue

    # Bakta writes features in genomic order per contig; keep that order, and
    # keep contigs separate so adjacency never spans a contig break.
    order = defaultdict(list)
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8 or f[1].lower() != 'cds': continue
        order[f[0]].append((f[5], f[7]))      # (locus_tag, product)

    kinds = {}
    for ctg, genes in order.items():
        for i, (lt, pr) in enumerate(genes):
            k = classify(pr)
            if k: kinds[(ctg, i)] = (k, lt, pr)

    n_caz = sum(1 for v in kinds.values() if v[0] == 'CAZy')
    n_c   = sum(1 for v in kinds.values() if v[0] == 'susC')
    n_d   = sum(1 for v in kinds.values() if v[0] == 'susD')

    # find susC/susD pairs, then the CAZymes within WINDOW of the pair
    in_pul, npul = set(), 0
    for ctg, genes in order.items():
        idx_c = [i for i in range(len(genes)) if kinds.get((ctg, i), ('',))[0] == 'susC']
        idx_d = [i for i in range(len(genes)) if kinds.get((ctg, i), ('',))[0] == 'susD']
        used_d = set()
        for i in idx_c:
            near = [j for j in idx_d if j not in used_d and abs(j - i) <= MIN_PAIR_GAP]
            if not near: continue
            j = min(near, key=lambda x: abs(x - i)); used_d.add(j)
            npul += 1
            lo, hi = max(0, min(i, j) - WINDOW), min(len(genes), max(i, j) + WINDOW + 1)
            cz = [genes[k][0] for k in range(lo, hi)
                  if kinds.get((ctg, k), ('',))[0] == 'CAZy']
            for k in range(lo, hi):
                if kinds.get((ctg, k), ('',))[0] == 'CAZy':
                    in_pul.add((ctg, k))
            loci.append(dict(mag_id=mid(m), genus=gen(m), contig=ctg,
                             susC=genes[i][0], susD=genes[j][0],
                             span_genes=hi - lo, n_cazymes=len(cz),
                             cazymes=';'.join(cz)))
    # SENSITIVITY BOUND. Requiring a susC/susD PAIR undercounts PULs whenever
    # susD is not annotated — susD is far less conserved than susC and Bakta
    # often leaves it as a generic product. ETO_LIQ_MAG_013 shows the problem
    # directly: 21 susC, ZERO susD, so it scores 0% by the strict rule despite
    # 73 CAZymes. Compute an UPPER bound from susC proximity alone.
    near_c = set()
    for ctg, genes in order.items():
        for i in range(len(genes)):
            if kinds.get((ctg, i), ('',))[0] != 'susC': continue
            for k in range(max(0, i - WINDOW), min(len(genes), i + WINDOW + 1)):
                if kinds.get((ctg, k), ('',))[0] == 'CAZy':
                    near_c.add((ctg, k))
    priv = 100 * len(in_pul) / n_caz if n_caz else 0.0
    priv_hi = 100 * len(near_c) / n_caz if n_caz else 0.0
    for (ctg, i), (k, lt, pr) in kinds.items():
        if k == 'CAZy' and (ctg, i) not in in_pul:
            nonpul.append(dict(mag_id=mid(m), genus=gen(m), locus=lt, product=pr))
    if n_caz or n_c:
        print(f"  {mid(m):<17}{gen(m)[:19]:<20}{n_caz:>6}{n_c:>6}{n_d:>6}"
              f"{npul:>6}{len(in_pul):>9}{priv:>8.1f}%{priv_hi:>9.1f}%")
    summary.append(dict(mag_id=mid(m), working_name=m, genus=gen(m),
                        n_cazymes=n_caz, n_susC=n_c, n_susD=n_d, n_PULs=npul,
                        cazymes_in_PUL=len(in_pul), pct_in_PUL=f"{priv:.1f}",
                        cazymes_near_susC=len(near_c), pct_in_PUL_upper=f"{priv_hi:.1f}"))

for name, rr in (('pul_summary.tsv', summary), ('pul_loci.tsv', loci),
                 ('nonPUL_cazymes.tsv', nonpul)):
    with open(os.path.join(P, name), 'w', newline='') as fh:
        if rr:
            w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rr[0]))
            w.writeheader(); w.writerows(rr)

# ---------------------------------------------------------------- the verdict
g = next((r for r in summary if 'gaoshiqia' in (r['genus'] or '').lower()), None)
print()
print("=" * 100)
print("  VERDICT — PUBLIC GOOD OR PRIVATE DIGESTER?")
print("=" * 100)
if not g:
    print("  Gaoshiqia not found in the annotation set.")
else:
    print(f"  {g['mag_id']} ({g['genus']}): {g['n_cazymes']} CAZymes, "
          f"{g['n_susC']} susC, {g['n_susD']} susD, {g['n_PULs']} PULs")
    print(f"  {g['pct_in_PUL']}% of its CAZymes sit inside a polysaccharide "
          f"utilisation locus (upper bound {g['pct_in_PUL_upper']}%).")
    print()
    pct = float(g['pct_in_PUL'])
    if pct >= 50:
        print("  >> PRIVATE DIGESTER. The majority of its carbohydrate-active enzymes")
        print("     are PUL-associated: bound at the cell surface, cleaved there, and")
        print("     imported as oligosaccharides. This is selfish uptake, and it means")
        print("     Gaoshiqia is NOT the community's polysaccharide provider.")
        print()
        print("     CONSEQUENCE FOR THE CHAPTER. Do not write that it 'opens complex")
        print("     organic matter for everyone else' — that was my phrasing at step 20")
        print("     and the gene arrangement contradicts it. The carbon entering the")
        print("     fermenter-methanogen chain has to come from somewhere else, which")
        print("     makes the pectin -> methanol route MORE important, not less:")
        print("     pectinesterase releases free methanol into the medium rather than")
        print("     importing it.")
    elif pct >= 20:
        print("  >> MIXED. A substantial minority of CAZymes stand outside any PUL and")
        print("     are candidates for export. Check nonPUL_cazymes.tsv for what they")
        print("     act on — if they target the abundant polymers, a partial public")
        print("     good is defensible.")
    else:
        print("  >> DISPERSED. Few CAZymes are PUL-associated, so the selfish-uptake")
        print("     reading does not apply and the public-good interpretation stands.")

print()
print("  ACROSS THE CATALOGUE")
withpul = [r for r in summary if int(r['n_PULs']) > 0]
print(f"    genomes with at least one PUL: {len(withpul)}/{len(summary)}")
for r in sorted(withpul, key=lambda x: -int(x['n_PULs']))[:8]:
    print(f"      {r['mag_id']:<17}{r['genus'][:20]:<21}{r['n_PULs']:>3} PULs, "
          f"{r['pct_in_PUL']}% of CAZymes private")
orphan = [r for r in summary if int(r['n_susC']) >= 5 and int(r['n_susD']) == 0]
if orphan:
    print()
    print("  susC WITHOUT susD — PULs my strict rule cannot see:")
    for r in orphan:
        print(f"      {r['mag_id']:<17}{r['genus'][:20]:<21}{r['n_susC']:>3} susC, "
              f"0 susD, {r['n_cazymes']} CAZymes, upper-bound "
              f"{r['pct_in_PUL_upper']}% PUL-associated")
    print("    susD is poorly conserved and frequently left as a generic product,")
    print("    so these are very likely real PULs. Use the upper bound for them.")

print()
print("  LIMIT, TO STATE PLAINLY: this is genomic CONTEXT, not localisation. A")
print("  non-PUL CAZyme is a candidate for export, not a demonstrated secreted")
print("  protein. Confirming that needs SignalP 6.0 (free for academic use, manual")
print("  download from DTU) or proteomics. Write 'PUL-associated' and 'non-PUL'.")
PY

echo
echo "######## VALIDATION GATE ########"
check "PUL summary" $([ -s "$P/pul_summary.tsv" ] && echo 0 || echo 1)
runlog "23b_pul_analysis complete"
echo
echo "######## 23b DONE ########"
