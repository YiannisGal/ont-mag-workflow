#!/usr/bin/env bash
# =============================================================================
# 23c_secreted_vs_pul.sh  --  do the two secretome results actually disagree?
#
# THE APPARENT CONTRADICTION
#   23b said Gaoshiqia is a PRIVATE digester: 54.2% of its 253 CAZymes sit inside
#   a polysaccharide utilisation locus, which is surface-bound selfish uptake.
#   23 then found it has 800 signal-peptide proteins (15.3%, by far the highest
#   in the catalogue) including 114 SECRETED CAZymes — five times more than any
#   other genome. Read naively, the second result says it exports its enzymes
#   after all.
#
# WHY THAT READING IS PROBABLY WRONG
#   Gaoshiqia is a Bacteroidota — a GRAM-NEGATIVE. A signal peptide in a
#   Gram-negative gets a protein across the INNER membrane, into the periplasm.
#   It does not get it out of the cell. Reaching the medium needs a second,
#   separate export step across the outer membrane. So SignalP-positive means
#   "not cytoplasmic", NOT "released into the environment".
#
#   And PUL enzymes are precisely non-cytoplasmic: they work at the cell surface
#   and in the periplasm, which is why they carry signal peptides. The two
#   results are then not in conflict at all — they are measuring the same
#   proteins from two directions.
#
#   The arithmetic already hints at it: 137 PUL-associated + 114 secreted = 251,
#   against 253 CAZymes total. That is suspiciously close to complete overlap.
#
# THE TEST
#   Intersect the two sets directly, by locus tag.
#     HIGH overlap  -> the same proteins. Non-cytoplasmic AND PUL-bound, i.e.
#                      surface-attached selfish uptake. 23b stands.
#     LOW overlap   -> two distinct populations: a private PUL system PLUS a
#                      genuinely exported enzyme set. Gaoshiqia would then be
#                      doing both, and the public-good reading partly revives.
#
#   This is the discriminating measurement, and neither script could make it
#   alone.
#
# OUTPUT   23_extras/pul/secreted_vs_pul.tsv
#          per genome: secreted CAZymes, PUL CAZymes, overlap, and the verdict
#
# RUNTIME  seconds.
# USAGE    ws2 && bash scripts/23c_secreted_vs_pul.sh 2>&1 | tee logs/23c.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="23c_sec_vs_pul"
A="$WS2/09_annotation"; F="$WS2/07_mags"; E="$WS2/23_extras"; P="$E/pul"

echo "######## 23c — SECRETED vs PUL-ASSOCIATED ########"
echo "  started: $(date -Is)"
runlog "23c_secreted_vs_pul started"
[ -s "$E/secretome.tsv" ] || { echo "  !! run 23_genome_extras.sh (DeepSig) first"; exit 1; }
[ -s "$P/pul_loci.tsv" ]  || { echo "  !! run 23b_pul_analysis.sh first"; exit 1; }

python3 - "$A" "$F" "$E" "$P" <<'PY'
import sys, os, csv, glob
from collections import defaultdict

A, F, E, P = sys.argv[1:5]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return (tax.get(m, {}) or {}).get('genus', '') or '-'
def phy(m): return ((tax.get(m, {}) or {}).get('phylum', '') or '')

# ---- locus tags inside a PUL, per genome
pul_loci = defaultdict(set)
id2wn = {mid(w): w for w in tax}
for r in csv.DictReader(open(os.path.join(P, 'pul_loci.tsv')), delimiter='\t'):
    wn = id2wn.get(r.get('mag_id',''), r.get('mag_id',''))
    for lt in (r.get('cazymes') or '').split(';'):
        if lt.strip(): pul_loci[wn].add(lt.strip())

# ---- CAZyme locus tags with a signal peptide, per genome
CAZY = ('glycoside hydrolase','glycosyl hydrolase','polysaccharide lyase',
        'carbohydrate esterase','glucanase','xylanase','chitinase','amylase',
        'pectinesterase','pectate lyase','cellulase','agarase','alginate lyase',
        'galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
        'arabinofuranosidase','levanase','pullulanase','chondroitin')

prod = defaultdict(dict)
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    for line in open(t):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].lower() == 'cds':
            prod[n][f[5]] = f[7]

sec_caz = defaultdict(set)
for p in sorted(glob.glob(os.path.join(E, 'deepsig', '*.tsv'))):
    m = os.path.basename(p)[:-4]
    for line in open(p):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 3: continue
        if 'signal' not in ' '.join(f[1:]).lower(): continue
        pr = (prod.get(m, {}).get(f[0]) or '').lower()
        if any(x in pr for x in CAZY): sec_caz[m].add(f[0])

print()
print("=" * 104)
print("  DO THE TWO MEASURES DESCRIBE THE SAME PROTEINS?")
print("=" * 104)
print(f"  {'MAG ID':<17}{'genus':<20}{'gram':<7}{'secCAZy':>8}{'pulCAZy':>9}"
      f"{'both':>7}{'sec&notPUL':>12}{'% of sec':>10}")
rows = []
for m in sorted(set(sec_caz) | set(pul_loci), key=mid):
    s, u = sec_caz.get(m, set()), pul_loci.get(m, set())
    both = s & u
    only_s = s - u
    gram = 'neg' if any(x in phy(m).lower() for x in
                        ('bacteroidota','pseudomonadota','desulfobacterota',
                         'verrucomicrobiota','planctomycetota')) else 'pos/arc'
    pct = 100 * len(only_s) / len(s) if s else 0.0
    print(f"  {mid(m):<17}{gen(m)[:19]:<20}{gram:<7}{len(s):>8}{len(u):>9}"
          f"{len(both):>7}{len(only_s):>12}{pct:>9.1f}%")
    rows.append(dict(mag_id=mid(m), genus=gen(m), gram=gram,
                     secreted_cazymes=len(s), pul_cazymes=len(u),
                     overlap=len(both), secreted_not_pul=len(only_s),
                     pct_secreted_outside_pul=f"{pct:.1f}"))
with open(os.path.join(P, 'secreted_vs_pul.tsv'), 'w', newline='') as fh:
    if rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

g = next((r for r in rows if 'gaoshiqia' in r['genus'].lower()), None)
print()
print("=" * 104)
print("  VERDICT")
print("=" * 104)
if not g:
    print("  Gaoshiqia not found.")
else:
    pct = float(g['pct_secreted_outside_pul'])
    print(f"  {g['mag_id']} ({g['genus']}, Gram-{g['gram']}): "
          f"{g['secreted_cazymes']} signal-peptide CAZymes, "
          f"{g['pul_cazymes']} PUL-associated, {g['overlap']} are BOTH.")
    print(f"  {g['secreted_not_pul']} carry a signal peptide but sit outside any PUL "
          f"({pct:.1f}% of the secreted set).")
    print()
    if pct < 40:
        print("  >> THE TWO MEASURES AGREE. Most signal-peptide CAZymes are the PUL")
        print("     enzymes themselves. In a Gram-negative a signal peptide reaches the")
        print("     PERIPLASM, not the medium — so 'secreted' here means surface- and")
        print("     periplasm-facing, which is exactly what a PUL is. 23b stands:")
        print("     Gaoshiqia is a private digester, and the high signal-peptide count")
        print("     is a consequence of its PUL machinery rather than evidence against.")
    else:
        print("  >> TWO DISTINCT SETS. A substantial fraction of its signal-peptide")
        print("     CAZymes are NOT PUL-associated, so Gaoshiqia appears to run a")
        print("     private PUL system AND export a separate enzyme set. The")
        print("     public-good reading partly revives — but note that even these")
        print("     need a second export step across the outer membrane, which a")
        print("     signal peptide alone does not demonstrate.")
    print()
    print("  EITHER WAY, THE GRAM-NEGATIVE CAVEAT IS THE POINT. Do not write")
    print("  'secreted' for a Bacteroidota on the strength of a signal peptide.")
    print("  Write 'predicted non-cytoplasmic' and say which compartment the")
    print("  genomic context supports.")

print()
print("  THE CONTRAST WORTH REPORTING")
print("  Bacteroidota carry high signal-peptide fractions because Gram-negatives")
print("  route far more proteins to the periplasm. Compare within a Gram type, not")
print("  across: Gaoshiqia 15.3%, Lentimicrobium 13.5% and GCA-2746605 11.2% are")
print("  the three Bacteroidota; almost every Bacillota sits near 5%. That is cell")
print("  envelope architecture, not ecological strategy, and reporting it as the")
print("  latter would be a mistake.")
PY

echo
echo "######## VALIDATION GATE ########"
check "secreted vs PUL table" $([ -s "$P/secreted_vs_pul.tsv" ] && echo 0 || echo 1)
runlog "23c_secreted_vs_pul complete"
echo
echo "######## 23c DONE ########"
