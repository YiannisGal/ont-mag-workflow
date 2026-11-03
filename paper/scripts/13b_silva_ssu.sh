#!/usr/bin/env bash
# =============================================================================
# 13b_silva_ssu.sh  --  SSU rRNA profiling: the rare and uncultured majority
#
# WHY THIS AND NOT JUST KRAKEN2
#   Kraken2 can only report lineages that have SEQUENCED GENOMES. Etoliko
#   sediment is dominated by uncultured taxa — 22 of 29 MAGs have no
#   species-level match in GTDB, and MAGs cover only ~2% of the sediment
#   community. A genome-database classifier is structurally blind to exactly
#   the organisms this chapter is about.
#
#   SILVA is built from rRNA, not genomes, so it contains tens of thousands of
#   environmental clades that were never cultured or assembled. For rare taxa
#   in an anoxic sediment it is far more sensitive than any Kraken2 build.
#
#   SECOND PURPOSE: this produces a profile in the SAME UNITS as the 16S
#   amplicon chapter — SSU rRNA, SILVA taxonomy. The metagenome and amplicon
#   datasets become directly comparable, from the same DNA, which is the bridge
#   the integration phase needs. Arc_Otu22 can be searched against these reads
#   directly.
#
# METHOD, AND WHY IT SUITS ONT
#   ONT reads are long, so a read overlapping an rRNA operon often spans most
#   of the SSU — far more taxonomic signal than a 150 bp fragment. Reads are
#   mapped to SILVA SSU NR99 with minimap2 (map-ont), filtered on identity and
#   alignment length, and the best hit's SILVA lineage is taken.
#
# THE LIMITATION, STATED PLAINLY
#   rRNA read counts are NOT cell counts: rRNA operon copy number varies 1-15
#   between taxa, so abundances are biased toward high-copy lineages. The same
#   bias affects amplicon data, which is precisely why the two are comparable.
#   Report as relative rRNA abundance, never as cell abundance.
#
# OUTPUT   13_profile/silva/
#     {SED,LIQ}_ssu_profile.tsv     ranked lineages with read counts
#     lineages_not_in_MAGs.tsv      THE KEY TABLE: present in the community,
#                                   absent from the 29-MAG catalogue
#
# RUNTIME  1-3 h (mapping the full read sets), download ~1 GB.
# USAGE    ws2 && bash scripts/13b_silva_ssu.sh 2>&1 | tee logs/13b_silva.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13b_silva"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
# Output directory is overridable so the same PAF can be re-filtered at a
# different threshold WITHOUT destroying the previous profile. The first
# version hardcoded this, so a relaxed re-run silently replaced the strict one
# and the comparison then read the same numbers twice.
S="${SILVA_OUT:-$WS2/13_profile/silva}"; F="$WS2/07_mags"; mkdir -p "$S"
echo "  [out] $S   (override with SILVA_OUT=)"

MINID="${MINID:-0.80}"      # ONT error rate makes a strict identity cutoff wrong
MINLEN="${MINLEN:-400}"     # alignment length, not read length

echo "######## 13b — SILVA SSU rRNA PROFILING ########"
echo "  started: $(date -Is)  min identity $MINID  min aln length $MINLEN"
runlog "13b_silva_ssu started"

# ------------------------------------------------------------- DATABASE
SILVA="${SILVA_SSU:-${DB_DIR}/silva/SILVA_SSURef_NR99.fasta}"
mkdir -p "$(dirname "$SILVA")"
if [ ! -s "$SILVA" ]; then
  echo "  SILVA SSU NR99 not found — downloading (~200 MB compressed)"
  BASEU="https://ftp.arb-silva.de/current/Exports"
  # the release number changes; discover the current SSURef_NR99 file
  FN=$(wget -qO- "$BASEU/" 2>/dev/null \
       | grep -oE 'SILVA_[0-9._]+_SSURef_NR99_tax_silva\.fasta\.gz' | head -1)
  if [ -z "$FN" ]; then
    echo "  !! could not auto-detect the SILVA filename."
    echo "     Browse $BASEU/ , pick SILVA_<ver>_SSURef_NR99_tax_silva.fasta.gz, then:"
    echo "       wget -c -O ${SILVA}.gz <url> && gunzip ${SILVA}.gz"
    echo "     or set SILVA_SSU=/path/to/existing.fasta and re-run."
    exit 1
  fi
  echo "  found: $FN"
  wget -c -O "${SILVA}.gz" "$BASEU/$FN" || { echo "  !! download failed"; exit 1; }
  gunzip -f "${SILVA}.gz"
fi
echo "  [db] $SILVA  ($(grep -c '^>' "$SILVA" 2>/dev/null || echo '?') sequences)"

MMI="$S/silva.mmi"
[ -s "$MMI" ] || { echo "  indexing SILVA ..."; run_tool minimap2 -x map-ont -d "$MMI" "$SILVA" 2>/dev/null; }

# ------------------------------------------------------------- MAP
map_sample() {  # map_sample <CODE> <reads>
  local code="$1" reads="$2"
  local paf="$S/${code}_vs_silva.paf.gz"
  [ -s "$paf" ] && { echo "  $code: already mapped"; return 0; }
  [ -s "$reads" ] || { echo "  !! missing $reads"; return 1; }
  echo "  $code: mapping to SILVA ($(date +%H:%M)) ..."
  # PAF is enough and far smaller than SAM; secondary alignments suppressed so
  # each read contributes its single best lineage assignment
  run_tool minimap2 -x map-ont -t "$THREADS" --secondary=no -c "$MMI" "$reads" 2>/dev/null \
    | gzip > "$paf"
  echo "    alignments: $(zcat "$paf" | wc -l)"
}
map_sample SED "$SED_READS"
map_sample LIQ "$LIQ_READS"

# ------------------------------------------------------------- PROFILE
python3 - "$S" "$SILVA" "$F" "$MINID" "$MINLEN" "$WS2" <<'PY'
import sys, os, gzip, csv
from collections import Counter, defaultdict

S, SILVA, F, MINID, MINLEN, WS2 = sys.argv[1:7]
MINID, MINLEN = float(MINID), int(MINLEN)

# SILVA headers: >ACC.start.end Domain;Phylum;Class;Order;Family;Genus;Species
lineage = {}
with open(SILVA) as fh:
    for line in fh:
        if line.startswith('>'):
            h = line[1:].rstrip('\n')
            acc, _, tax = h.partition(' ')
            lineage[acc] = tax
print(f"  SILVA lineages loaded: {len(lineage)}")

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus']

def split_lin(t):
    p = [x.strip() for x in t.split(';') if x.strip()]
    d = {}
    for i, r in enumerate(RANKS):
        d[r] = p[i] if i < len(p) else ''
    return d

profiles = {}
for code in ('SED', 'LIQ'):
    paf = os.path.join(S, f'{code}_vs_silva.paf.gz')
    if not os.path.exists(paf): continue
    best = {}
    n_aln = 0
    with gzip.open(paf, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12: continue
            n_aln += 1
            q, tname = f[0], f[5]
            try:
                nmatch, alen, mapq = int(f[9]), int(f[10]), int(f[11])
            except ValueError:
                continue
            if alen < MINLEN: continue
            ident = nmatch / alen if alen else 0
            if ident < MINID: continue
            score = nmatch
            if q not in best or score > best[q][0]:
                best[q] = (score, tname, ident, alen)
    print(f"\n  {code}: {n_aln} alignments -> {len(best)} reads assigned "
          f"(identity >= {MINID:.0%}, alignment >= {MINLEN} bp)")

    counts = {r: Counter() for r in RANKS}
    full = Counter()
    for q, (sc, tname, ident, alen) in best.items():
        tax = lineage.get(tname) or lineage.get(tname.split()[0]) or ''
        if not tax: continue
        d = split_lin(tax)
        for r in RANKS:
            if d[r]: counts[r][d[r]] += 1
        key = ';'.join(d[r] for r in RANKS if d[r])
        if key: full[key] += 1
    profiles[code] = (counts, full, len(best))

    tot = sum(counts['domain'].values()) or 1
    print(f"    domain composition: " + ', '.join(
        f"{k} {100*v/tot:.1f}%" for k, v in counts['domain'].most_common(4)))
    print(f"\n    top 20 genera (% of SSU reads)")
    gtot = sum(counts['genus'].values()) or 1
    for g, c in counts['genus'].most_common(20):
        print(f"      {c:>7}  {100*c/gtot:>6.2f}%  {g}")

    with open(os.path.join(S, f'{code}_ssu_profile.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['rank', 'taxon', 'reads', 'pct_of_assigned'])
        for r in RANKS:
            t = sum(counts[r].values()) or 1
            for k, c in counts[r].most_common():
                w.writerow([r, k, c, f"{100*c/t:.4f}"])

# --------------------------------------------- WHAT THE MAGs MISS
print()
print("=" * 100)
print("  LINEAGES IN THE COMMUNITY BUT ABSENT FROM THE 29-MAG CATALOGUE")
print("=" * 100)
mag_taxa = set()
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        for k in ('phylum', 'class', 'order', 'family', 'genus'):
            v = (r.get(k) or '').strip()
            if v: mag_taxa.add(v.lower())
print(f"  taxa represented in the MAG catalogue: {len(mag_taxa)}")

rows = []
for code, (counts, full, nassigned) in profiles.items():
    gtot = sum(counts['genus'].values()) or 1
    missing = [(g, c) for g, c in counts['genus'].most_common()
               if g.lower() not in mag_taxa]
    covered = sum(c for g, c in counts['genus'].items() if g.lower() in mag_taxa)
    print(f"\n  {code}: {100*covered/gtot:.1f}% of genus-assigned SSU reads belong to "
          f"genera present in the MAG set")
    print(f"  top 25 genera NOT represented by any MAG:")
    for g, c in missing[:25]:
        print(f"      {c:>7}  {100*c/gtot:>6.2f}%  {g}")
        rows.append(dict(sample=code, genus=g, reads=c,
                         pct_of_ssu=f"{100*c/gtot:.4f}"))

if rows:
    with open(os.path.join(S, 'lineages_not_in_MAGs.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(rows)

print()
print("  This is the honest answer to 'what is in this lagoon that your")
print("  catalogue does not describe'. A high-abundance genus with no MAG is a")
print("  specific, nameable gap — far stronger than saying 2% was recovered.")
print()
print("  CAVEAT: rRNA read counts are not cell counts. Operon copy number varies")
print("  1-15 across taxa, so abundances are biased toward high-copy lineages.")
print("  The SAME bias affects 16S amplicon data, which is exactly why these two")
print("  datasets are comparable to each other. Report as relative rRNA abundance.")
PY

echo
echo "######## VALIDATION GATE ########"
for c in SED LIQ; do
  check "$c mapped to SILVA" $([ -s "$S/${c}_vs_silva.paf.gz" ] && echo 0 || echo 1)
  check "$c SSU profile"     $([ -s "$S/${c}_ssu_profile.tsv" ] && echo 0 || echo 1)
done
check "missing-lineage table" $([ -s "$S/lineages_not_in_MAGs.tsv" ] && echo 0 || echo 1)

runlog "13b_silva_ssu complete"
echo
echo "######## 13b DONE ########"
