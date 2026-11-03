#!/usr/bin/env bash
# =============================================================================
# 17_comparative_genomics.sh  --  what differs between close relatives?
#
# THREE QUESTIONS THE MAG CATALOGUE RAISES BUT CANNOT ANSWER
#
#   1. THE Chakrabartyella PAIR.
#      ETO_SED_MAG_002 (sediment) and ETO_LIQ_MAG_012 (enrichment) are 99.80%
#      ANI over 90% of the genome — the same organism, sampled once from the
#      lagoon and once after it grew under defined amendments. Any gene present
#      in one and absent from the other is a candidate for what enabled
#      enrichment. This is a natural experiment that the dataset handed us.
#
#      CAVEAT BUILT INTO THE OUTPUT: the two MAGs are 88.9% and 89.6% complete,
#      so a gene "absent" from one may simply be in the missing 10%. Differences
#      are reported WITH the completeness caveat attached, and only differences
#      involving genes present in the more complete genome are worth pursuing.
#
#   2. CO-OCCURRING CONGENERS.
#      Three Chakrabartyella species, three Methanofollis, two Methanosarcina,
#      two Sedimentibacter, two Lacrimispora — all in one small lagoon system.
#      Co-occurring congeners must differ somewhere or one would exclude the
#      other. Pangenome analysis asks where.
#
#   3. DID THE FAST GROWERS SIMPLY WIN?
#      Enrichment selects for growth rate. Codon usage bias in ribosomal
#      proteins relative to the whole genome is a established proxy for maximal
#      growth rate (the basis of gRodon's models). Comparing culture MAGs with
#      sediment MAGs asks whether enrichment simply selected the fast growers.
#      Reported as a RELATIVE RANKING, never as a doubling time — the absolute
#      calibration requires a trained model this script does not implement.
#
# INPUT    09_annotation/<MAG>/<MAG>.faa and .ffn, 07_mags/MAG_TAXONOMY.tsv
# OUTPUT   17_compgen/
#
# RUNTIME  10-30 min (MMseqs2 clustering dominates).
# USAGE    ws2 && bash scripts/17_comparative_genomics.sh 2>&1 | tee logs/17_compgen.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="17_compgen"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; G="$WS2/17_compgen"
mkdir -p "$G/clusters" "$G/pairs"

echo "######## 17 — COMPARATIVE GENOMICS ########"
echo "  started: $(date -Is)"
runlog "17_comparative_genomics started"
[ -d "$A" ] || { echo "  !! $A missing — run 10_annotate.sh"; exit 1; }
record_version "mmseqs" "${ENV_mmseqs:-env_mmseqs}" mmseqs version

# =====================================================================
# 1. ORTHOLOGUE CLUSTERING ACROSS ALL 29 MAGs
# =====================================================================
echo
echo "######## PROTEIN CLUSTERING (MMseqs2) ########"
ALLFAA="$G/all_proteins.faa"
if [ ! -s "$ALLFAA" ]; then
  : > "$ALLFAA"
  for d in "$A"/*/; do
    n=$(basename "$d"); faa="$d/${n}.faa"
    [ -s "$faa" ] || continue
    awk -v m="$n" '/^>/{sub(/^>/,""); split($0,a," "); print ">"m"@"a[1]; next}{print}' "$faa" >> "$ALLFAA"
  done
fi
echo "  proteins: $(grep -c '^>' "$ALLFAA")"

CLUST="$G/clusters/all_cluster.tsv"
if [ ! -s "$CLUST" ]; then
  echo "  clustering at 50% identity / 80% coverage (orthologue-scale) ..."
  run_tool mmseqs easy-cluster "$ALLFAA" "$G/clusters/all" "$TMPDIR/mmseqs_$$" \
      --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads "$THREADS" \
      > "$G/clusters/mmseqs.log" 2>&1 \
    || echo "  !! mmseqs failed — see $G/clusters/mmseqs.log"
fi
[ -s "$CLUST" ] || { echo "  !! no cluster table produced"; exit 1; }
echo "  clusters: $(cut -f1 "$CLUST" | sort -u | wc -l)"

# =====================================================================
# 2 + 3. PANGENOMES, THE PAIR, AND GROWTH-RATE PROXY
# =====================================================================
python3 - "$G" "$F" "$A" <<'PY'
import sys, os, csv, glob, math
from collections import defaultdict, Counter

G, F, A = sys.argv[1], sys.argv[2], sys.argv[3]

tax = {}
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    tax = {r['working_name']: r for r in csv.DictReader(open(p), delimiter='\t')}
mid = {k: (v.get('mag_id') or k) for k, v in tax.items()}

# cluster table: rep <tab> member ; ids are  <mag>@<locus>
members = defaultdict(list)
for line in open(os.path.join(G, 'clusters', 'all_cluster.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 2: continue
    members[f[0]].append(f[1])

def mag_of(x): return x.split('@')[0]

# per-MAG cluster membership
mag_clusters = defaultdict(set)
for rep, mem in members.items():
    for m in mem:
        mag_clusters[mag_of(m)].add(rep)
mags = sorted(mag_clusters)
print(f"  MAGs with clustered proteins: {len(mags)}")

# ---------------------------------------------------------------- PANGENOMES
GROUPS = defaultdict(list)
for m in mags:
    g = (tax.get(m, {}) or {}).get('genus', '')
    if g: GROUPS[g].append(m)
GROUPS = {g: v for g, v in GROUPS.items() if len(v) >= 2}

print()
print("=" * 100)
print("  PANGENOMES OF CO-OCCURRING CONGENERS")
print("=" * 100)
if not GROUPS:
    print("  no genus has 2+ MAGs")
rows = []
for g, ms in sorted(GROUPS.items(), key=lambda kv: -len(kv[1])):
    sets = {m: mag_clusters[m] for m in ms}
    core = set.intersection(*sets.values())
    pan = set.union(*sets.values())
    print(f"\n  {g}  ({len(ms)} genomes)")
    print(f"    pangenome {len(pan)} clusters | core {len(core)} "
          f"({100*len(core)/len(pan):.0f}%) | accessory {len(pan)-len(core)}")
    for m in sorted(ms, key=lambda m: mid.get(m, m)):
        uniq = sets[m] - set.union(*(sets[o] for o in ms if o != m))
        t = tax.get(m, {})
        print(f"      {mid.get(m,m):<17} {t.get('sample','?'):<4} "
              f"{len(sets[m]):>5} clusters, {len(uniq):>4} unique to it "
              f"(comp {t.get('completeness','?')}%)")
        rows.append(dict(genus=g, mag_id=mid.get(m, m), working_name=m,
                         sample=t.get('sample', ''), n_clusters=len(sets[m]),
                         n_unique=len(uniq), core=len(core), pan=len(pan),
                         completeness=t.get('completeness', '')))
if rows:
    with open(os.path.join(G, 'pangenomes.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(rows)

# ------------------------------------------------- THE Chakrabartyella PAIR
PAIR = [m for m in mags if mid.get(m) in ('ETO_SED_MAG_002', 'ETO_LIQ_MAG_012')]
print()
print("=" * 100)
print("  THE 99.8% ANI PAIR — sediment vs enrichment, same organism")
print("=" * 100)
if len(PAIR) != 2:
    print("  pair not found among clustered MAGs")
else:
    a, b = sorted(PAIR, key=lambda m: tax.get(m, {}).get('sample', ''))   # SED first
    ta, tb = tax.get(a, {}), tax.get(b, {})
    print(f"  {mid[a]} ({ta.get('sample')}, {ta.get('completeness')}% complete)  vs  "
          f"{mid[b]} ({tb.get('sample')}, {tb.get('completeness')}% complete)")
    only_a = mag_clusters[a] - mag_clusters[b]
    only_b = mag_clusters[b] - mag_clusters[a]
    shared = mag_clusters[a] & mag_clusters[b]
    print(f"  shared {len(shared)} clusters | only in {ta.get('sample')} {len(only_a)} "
          f"| only in {tb.get('sample')} {len(only_b)}")

    # what ARE those genes? read products from the Bakta tables
    prod = {}
    for m in (a, b):
        tsv = os.path.join(A, m, f'{m}.tsv')
        if not os.path.exists(tsv): continue
        for line in open(tsv):
            if line.startswith('#'): continue
            f_ = line.rstrip('\n').split('\t')
            if len(f_) > 7 and f_[1] in ('cds', 'CDS'):
                prod[f"{m}@{f_[5]}"] = (f_[6], f_[7])

    def describe(clusters, owner, label, limit=25):
        print(f"\n  --- present only in {label} ({len(clusters)} clusters) ---")
        seen = Counter()
        for c in clusters:
            for x in members[c]:
                if mag_of(x) == owner and x in prod:
                    gene, pr = prod[x]
                    if 'hypothetical' in pr.lower(): pr = 'hypothetical protein'
                    seen[(gene or '-', pr)] += 1
                    break
        hyp = sum(v for (g_, p_), v in seen.items() if p_ == 'hypothetical protein')
        named = [(k, v) for k, v in seen.items() if k[1] != 'hypothetical protein']
        print(f"      {hyp} hypothetical, {len(named)} with a product name")
        for (gene, pr), n in sorted(named, key=lambda kv: -kv[1])[:limit]:
            print(f"      {n:>3}x {gene:<12} {pr[:70]}")

    describe(only_a, a, f"{ta.get('sample')} ({mid[a]})")
    describe(only_b, b, f"{tb.get('sample')} ({mid[b]})")

    with open(os.path.join(G, 'pairs', 'chakrabartyella_pair.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['cluster', 'present_in', 'gene', 'product'])
        for clusters, owner, lab in ((only_a, a, mid[a]), (only_b, b, mid[b])):
            for c in clusters:
                for x in members[c]:
                    if mag_of(x) == owner and x in prod:
                        w.writerow([c, lab, prod[x][0], prod[x][1]]); break

    print()
    print("  CAVEAT: both genomes are ~89% complete, so a cluster 'absent' from one")
    print("  may lie in its missing 10%. Treat these as CANDIDATES; confirm any gene")
    print("  that matters by mapping reads to it directly, not by absence from a MAG.")

# ------------------------------------------------- GROWTH RATE PROXY
print()
print("=" * 100)
print("  CODON USAGE BIAS — a relative proxy for maximal growth rate")
print("=" * 100)

def codon_counts(seq):
    c = Counter()
    for i in range(0, len(seq) - 2, 3):
        cod = seq[i:i+3].upper()
        if len(cod) == 3 and set(cod) <= set('ACGT'):
            c[cod] += 1
    return c

def read_ffn(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name: yield name, ''.join(buf)
                name, buf = line[1:], []
            else: buf.append(line)
    if name: yield name, ''.join(buf)

results = []
for m in mags:
    ffn = os.path.join(A, m, f'{m}.ffn')
    tsv = os.path.join(A, m, f'{m}.tsv')
    if not (os.path.exists(ffn) and os.path.exists(tsv)): continue
    ribo = set()
    for line in open(tsv):
        if line.startswith('#'): continue
        f_ = line.rstrip('\n').split('\t')
        if len(f_) > 7 and f_[1] in ('cds', 'CDS') and 'ribosomal protein' in f_[7].lower():
            ribo.add(f_[5])
    if len(ribo) < 15: continue

    hi, all_ = Counter(), Counter()
    for hdr, seq in read_ffn(ffn):
        tag = hdr.split()[0]
        cc = codon_counts(seq)
        all_ += cc
        if tag in ribo: hi += cc

    # per-amino-acid relative synonymous codon usage distance between the
    # highly expressed set and the genome as a whole
    CODON_AA = {}
    bases = 'TCAG'
    aas = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
    i = 0
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                CODON_AA[b1+b2+b3] = aas[i]; i += 1
    by_aa = defaultdict(list)
    for cod, aa in CODON_AA.items():
        if aa != '*': by_aa[aa].append(cod)

    d, n = 0.0, 0
    for aa, cods in by_aa.items():
        if len(cods) < 2: continue
        th, ta = sum(hi[c] for c in cods), sum(all_[c] for c in cods)
        if th < 20 or ta < 20: continue
        for c in cods:
            d += abs(hi[c]/th - all_[c]/ta)
        n += 1
    if n == 0: continue
    cub = d / n
    t = tax.get(m, {})
    results.append((cub, mid.get(m, m), t.get('sample', ''), t.get('genus', ''),
                    len(ribo)))

if results:
    results.sort(reverse=True)
    print(f"  {'rank':<5}{'MAG ID':<17}{'smp':<5}{'genus':<20}{'CUB':>8}{'ribo genes':>12}")
    for i, (cub, mag, smp, gen, nr) in enumerate(results, 1):
        print(f"  {i:<5}{mag:<17}{smp:<5}{(gen or '-')[:20]:<20}{cub:>8.4f}{nr:>12}")
    sed = [c for c, _, s, _, _ in results if s == 'SED']
    liq = [c for c, _, s, _, _ in results if s == 'LIQ']
    if sed and liq:
        ms = sum(sed)/len(sed); ml = sum(liq)/len(liq)
        print()
        print(f"  mean CUB — sediment MAGs {ms:.4f} (n={len(sed)}) | "
              f"culture MAGs {ml:.4f} (n={len(liq)})")
        print("  Higher bias = stronger selection for translational efficiency =")
        print("  faster potential growth. If culture MAGs are consistently higher,")
        print("  enrichment selected fast growers and the cultures are NOT a")
        print("  random sample of the sediment community.")
    with open(os.path.join(G, 'codon_usage_bias.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['rank','mag_id','sample','genus','codon_usage_bias','n_ribosomal_genes'])
        for i, (cub, mag, smp, gen, nr) in enumerate(results, 1):
            w.writerow([i, mag, smp, gen, f"{cub:.4f}", nr])
else:
    print("  not enough ribosomal-protein annotations to compute the proxy")

print()
print("  CUB here is a RELATIVE index, not a doubling time. It ranks genomes")
print("  against each other within this dataset; converting to hours requires a")
print("  calibrated model (gRodon) that this script does not implement.")
PY

echo
echo "######## VALIDATION GATE ########"
check "cluster table"     $([ -s "$G/clusters/all_cluster.tsv" ] && echo 0 || echo 1)
check "pangenome table"   $([ -s "$G/pangenomes.tsv" ] && echo 0 || echo 1)
check "codon usage table" $([ -s "$G/codon_usage_bias.tsv" ] && echo 0 || echo 1)
check "pair comparison"   $([ -s "$G/pairs/chakrabartyella_pair.tsv" ] && echo 0 || echo 1)

runlog "17_comparative_genomics complete"
echo
echo "######## 17 DONE ########"
