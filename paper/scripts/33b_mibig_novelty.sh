#!/usr/bin/env bash
# =============================================================================
# 33b_mibig_novelty.sh  --  do any of the 141 clusters match a KNOWN compound?
#
# WHY THIS COSTS NOTHING
#   Step 33 ran antiSMASH with `--cb-knownclusters`, so every region has already
#   been compared against MIBiG. The results are on disk in the per-genome output
#   and have never been read. This step parses them; it runs no new search.
#
# THE QUESTION, AND WHY THE ANSWER MATTERS EITHER WAY
#   22 of 29 MAGs have no GTDB species-level match. 81% of sediment genera have
#   no cultured representative (step 24). If that novelty extends to secondary
#   metabolism, most of these 141 regions should match nothing in MIBiG.
#
#   A high proportion with NO known-cluster hit is therefore a RESULT, not an
#   absence: uncharacterised biosynthetic potential in an uncharacterised
#   community. And any region that DOES match a characterised BGC is worth
#   naming, because it is the one place where a compound can be guessed.
#
# 🔴 WHAT "SIMILARITY" MEANS HERE, AND WHAT IT DOES NOT
#   knownclusterblast scores the number of GENES in the query region with BLAST
#   hits to genes in a MIBiG cluster. It is a measure of shared gene content, not
#   of product identity. Two clusters can share most of their genes and make
#   different molecules — a single tailoring enzyme or one adenylation-domain
#   specificity changes the product. So:
#
#     "N of M genes have homologues in BGC0000xxx"      <- what this measures
#     "this organism makes <compound>"                   <- NOT what this measures
#
#   Every number below is reported in the first phrasing. The second must not
#   appear in the thesis on this evidence.
#
# HOW THE DATA IS FOUND
#   antiSMASH's output layout differs between releases, so nothing about it is
#   asserted here. The script looks for the knownclusterblast text files first
#   (stable, human-readable), falls back to walking the JSON for any node whose
#   key mentions knowncluster, and REPORTS which route it used plus one raw
#   example, so the parse can be audited rather than trusted.
#
# INPUT   33_antismash/as/<genome>/  (knownclusterblast/ and/or <genome>.json)
#         33_antismash/BGC_SUMMARY.tsv   for in-situ status
# OUTPUT  33_antismash/mibig_hits.tsv · NOVELTY.txt
# RUNTIME under a minute.
# USAGE   ws2 && bash scripts/33b_mibig_novelty.sh 2>&1 | tee logs/33b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="33b_mibig"
B="$WS2/33_antismash"

MINSIM="${MINSIM:-0.50}"     # fraction of query-cluster genes with a homologue

echo "######## 33b — MIBiG COMPARISON, ALREADY COMPUTED ########"
echo "  started: $(date -Is)"
runlog "33b_mibig started"

[ -d "$B/as" ] || { echo "  !! no $B/as — run 33 first"; exit 1; }
[ -s "$B/BGC_SUMMARY.tsv" ] || { echo "  !! no BGC_SUMMARY.tsv — run 33 first"; exit 1; }

echo "  scanning $(ls -1d "$B"/as/*/ 2>/dev/null | wc -l) genome output directories"
NKCB=$(find "$B/as" -maxdepth 2 -type d -name 'knownclusterblast' 2>/dev/null | wc -l)
NJSON=$(find "$B/as" -maxdepth 2 -name '*.json' 2>/dev/null | wc -l)
echo "  knownclusterblast directories: $NKCB"
echo "  JSON files: $NJSON"
if [ "$NKCB" -eq 0 ] && [ "$NJSON" -eq 0 ]; then
  echo "  !! neither route is available. --cb-knownclusters may not have run."
  echo "     Check one log:  grep -i clusterblast $B/as_logs/*.log | head"
  runlog "**33b ABORTED: no knownclusterblast output**"
  exit 2
fi

python3 - "$B" "$MINSIM" <<'PY'
import sys, os, re, csv, glob, json
from collections import Counter, defaultdict

B, MINSIM = sys.argv[1], float(sys.argv[2])
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t: say('  ' + t); say('=' * W)

# ------------------------------------------------- in-situ status from step 33
meta = {}
for r in csv.DictReader(open(os.path.join(B, 'BGC_SUMMARY.tsv')), delimiter='\t'):
    meta[r['genome']] = r

# ------------------------------------------------- regions, from the GBK names
# A region GBK is named <record>.region00N.gbk, which gives contig and region
# number without needing to trust any other layout. Products come from the
# region feature only (a CDS /product is not a cluster type).
FEATKEY = re.compile(r'^ {5}(\S+)\s')
def region_products(path):
    prods, inreg = [], False
    for line in open(path, errors='replace'):
        m = FEATKEY.match(line)
        if m:
            inreg = m.group(1) in ('region', 'protocluster', 'cand_cluster')
            continue
        s = line.strip()
        if inreg and s.startswith('/product='):
            prods.append(s.split('=', 1)[1].strip('"'))
    return sorted(set(prods))

regions = {}                       # (genome, record, regnum) -> products
for d in sorted(glob.glob(os.path.join(B, 'as', '*'))):
    if not os.path.isdir(d): continue
    g = os.path.basename(d)
    for gbk in sorted(glob.glob(os.path.join(d, '*.region*.gbk'))):
        bn = os.path.basename(gbk)
        m = re.match(r'(.+)\.region0*(\d+)\.gbk$', bn)
        if not m: continue
        regions[(g, m.group(1), int(m.group(2)))] = region_products(gbk)
print(f"  regions indexed from GBK filenames: {len(regions)}")

# ------------------------------------------------------------------ ROUTE A
# knownclusterblast/<record>_c<N>.txt. The format has been stable for years:
#   a "Table of genes ..." block listing query genes, then "Significant hits:"
#   as a numbered list, then "Details:" with one ">>" block per hit carrying
#   "Number of proteins with BLAST hits to this cluster:".
hits = defaultdict(list)           # (genome, record, regnum) -> [dict]
route = None
example = None

kcb_files = glob.glob(os.path.join(B, 'as', '*', 'knownclusterblast', '*.txt'))
if kcb_files:
    route = 'knownclusterblast text'
    for f in kcb_files:
        g = os.path.basename(os.path.dirname(os.path.dirname(f)))
        bn = os.path.basename(f)
        m = re.match(r'(.+)_c(\d+)\.txt$', bn)
        if not m: continue
        rec, num = m.group(1), int(m.group(2))
        txt = open(f, errors='replace').read()
        if example is None: example = (f, txt[:1200])

        # query gene count: lines in the table before "Significant hits"
        nq = 0
        mtab = re.search(r'Table of genes.*?:\s*\n(.*?)\n\s*\n', txt, re.S)
        if mtab:
            nq = len([l for l in mtab.group(1).splitlines() if l.strip()])

        # one block per hit
        for blk in txt.split('>>')[1:]:
            acc = re.search(r'^\s*\d+\.\s*(\S+)', blk, re.M)
            src = re.search(r'^Source:\s*(.+)$', blk, re.M)
            typ = re.search(r'^Type:\s*(.+)$', blk, re.M)
            nb  = re.search(r'Number of proteins with BLAST hits to this cluster:\s*(\d+)',
                            blk)
            sc  = re.search(r'Cumulative BLAST score:\s*(\d+)', blk)
            if not acc: continue
            n = int(nb.group(1)) if nb else 0
            hits[(g, rec, num)].append(dict(
                accession=acc.group(1),
                compound=(src.group(1).strip() if src else ''),
                mibig_type=(typ.group(1).strip() if typ else ''),
                genes_matched=n,
                query_genes=nq,
                frac=(n / nq if nq else 0.0),
                score=int(sc.group(1)) if sc else 0))

# ------------------------------------------------------------------ ROUTE B
if not hits and glob.glob(os.path.join(B, 'as', '*', '*.json')):
    route = 'JSON walk'
    def walk(o, depth=0):
        """Yield any dict living under a key that mentions knowncluster."""
        if depth > 8: return
        if isinstance(o, dict):
            for k, v in o.items():
                if 'knowncluster' in str(k).lower():
                    yield v
                else:
                    yield from walk(v, depth + 1)
        elif isinstance(o, list):
            for v in o[:200]:
                yield from walk(v, depth + 1)

    for jf in sorted(glob.glob(os.path.join(B, 'as', '*', '*.json'))):
        g = os.path.basename(os.path.dirname(jf))
        try: data = json.load(open(jf))
        except Exception as e:
            print(f"  [warn] {g}: JSON unreadable ({e})"); continue
        for node in walk(data):
            if example is None:
                example = (jf, json.dumps(node)[:1200])
            res = node.get('results') if isinstance(node, dict) else None
            if not isinstance(res, list): continue
            for r in res:
                if not isinstance(r, dict): continue
                num = r.get('region_number') or r.get('region')
                rank = r.get('ranking') or []
                rec = r.get('record_id') or r.get('seq_record') or '?'
                for entry in rank:
                    if not (isinstance(entry, list) and len(entry) >= 2): continue
                    ref, sco = entry[0], entry[1]
                    if not isinstance(ref, dict): continue
                    n = 0
                    if isinstance(sco, dict):
                        n = sco.get('hits') or sco.get('core_gene_hits') or 0
                    try: num = int(num)
                    except (TypeError, ValueError): continue
                    hits[(g, rec, num)].append(dict(
                        accession=ref.get('accession', '?'),
                        compound=ref.get('description', ''),
                        mibig_type=ref.get('cluster_type', ''),
                        genes_matched=int(n), query_genes=0,
                        frac=0.0, score=0))

rule('0. WHERE THE NUMBERS CAME FROM')
if route is None:
    say("  !! no knownclusterblast results could be parsed by either route.")
    say("     Nothing is written. A table of zeros here would read as 'nothing")
    say("     matches MIBiG', which is a claim this run has not tested.")
    open(os.path.join(B, 'NOVELTY.txt'), 'w').write('\n'.join(rep) + '\n')
    sys.exit(3)
say(f"  route: {route}")
say(f"  regions with at least one MIBiG hit parsed: {len(hits)} of {len(regions)}")
if example:
    say()
    say(f"  raw example from {os.path.relpath(example[0], B)} — check the parse against it:")
    for line in example[1].splitlines()[:14]:
        say('    | ' + line[:92])

# ------------------------------------------------------------- per region
rows = []
for key, prods in sorted(regions.items()):
    g, rec, num = key
    hs = sorted(hits.get(key, []), key=lambda h: (-h['genes_matched'], -h['score']))
    m = meta.get(g, {})
    best = hs[0] if hs else None
    rows.append(dict(
        mag_id=m.get('mag_id', g), genome=g, contig=rec, region=num,
        products=';'.join(prods) or 'unknown',
        in_situ=m.get('in_situ', ''), sed_breadth=m.get('sed_breadth', ''),
        n_mibig_hits=len(hs),
        best_accession=best['accession'] if best else '',
        best_compound=best['compound'] if best else '',
        best_type=best['mibig_type'] if best else '',
        genes_matched=best['genes_matched'] if best else 0,
        query_genes=best['query_genes'] if best else 0,
        frac_matched=round(best['frac'], 3) if best else 0.0))

with open(os.path.join(B, 'mibig_hits.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ summary
rule('1. HOW MUCH OF THIS IS ALREADY DESCRIBED?')
nreg = len(rows)
nohit = [r for r in rows if r['n_mibig_hits'] == 0]
withhit = [r for r in rows if r['n_mibig_hits'] > 0]
strong = [r for r in withhit if r['frac_matched'] >= MINSIM]
say(f"  regions total                          {nreg:>5}")
say(f"  no MIBiG hit at all                    {len(nohit):>5}  ({100*len(nohit)/nreg:.0f}%)")
say(f"  some gene-level similarity             {len(withhit):>5}  ({100*len(withhit)/nreg:.0f}%)")
if any(r['query_genes'] for r in withhit):
    say(f"  >= {MINSIM:.0%} of query genes matched         {len(strong):>5}  "
        f"({100*len(strong)/nreg:.0f}%)")
else:
    say(f"  (query-gene counts unavailable on this route — the >= {MINSIM:.0%} line")
    say("   cannot be computed, and only raw gene-hit counts are meaningful)")

rule('2. THE REGIONS THAT DO RESEMBLE SOMETHING CHARACTERISED')
top = sorted(withhit, key=lambda r: (-r['frac_matched'], -r['genes_matched']))[:15]
if top:
    say(f"  {'mag_id':<20}{'reg':>4}{'match':>7}{'genes':>7}  {'accession':<12}{'products':<26}compound")
    say('-' * W)
    for r in top:
        fr = f"{100*r['frac_matched']:.0f}%" if r['query_genes'] else '-'
        gm = f"{r['genes_matched']}/{r['query_genes']}" if r['query_genes'] else str(r['genes_matched'])
        say(f"  {r['mag_id']:<20}{r['region']:>4}{fr:>7}{gm:>7}  "
            f"{r['best_accession'][:11]:<12}{r['products'][:25]:<26}{r['best_compound'][:28]}")
    say()
    say("  READ THIS AS SHARED GENE CONTENT, NOT AS A COMPOUND IDENTIFICATION.")
    say("  A cluster can share most of its genes with a characterised BGC and make")
    say("  a different molecule — one tailoring enzyme, or one adenylation-domain")
    say("  specificity, is enough. The defensible sentence is 'N of M genes have")
    say("  homologues in <accession>', never 'this organism makes <compound>'.")
else:
    say("  None. No region resembles a characterised cluster at gene level.")

rule('3. NOVELTY BY PRODUCT CLASS')
byclass_tot = Counter(); byclass_hit = Counter()
for r in rows:
    for p in r['products'].split(';'):
        if not p: continue
        byclass_tot[p] += 1
        if r['n_mibig_hits'] > 0: byclass_hit[p] += 1
say(f"  {'class':<32}{'regions':>8}{'with hit':>10}{'% described':>13}")
say('-' * W)
for p, n in byclass_tot.most_common():
    h = byclass_hit[p]
    say(f"  {p:<32}{n:>8}{h:>10}{100*h/n:>12.0f}%")

rule('4. NOVELTY BY IN-SITU STATUS')
for lab, want in (('confirmed in sediment', 'yes'), ('enrichment-only', 'no')):
    grp = [r for r in rows if r['in_situ'] == want]
    if not grp: continue
    nh = sum(1 for r in grp if r['n_mibig_hits'] == 0)
    say(f"  {lab:<26}{len(grp):>5} regions   {nh:>5} with no MIBiG hit "
        f"({100*nh/len(grp):.0f}%)")

rule('5. WHAT TO REPORT')
frac_novel = 100 * len(nohit) / nreg
say(f"  {len(nohit)} of {nreg} regions ({frac_novel:.0f}%) have no gene-level resemblance to any")
say("  characterised cluster in MIBiG.")
say()
if frac_novel >= 60:
    say("  That is the result, and it lines up with everything else this project has")
    say("  measured about novelty: 22 of 29 MAGs with no GTDB species match, 81% of")
    say("  sediment genera with no cultured representative (step 24), 0.47% of reads")
    say("  classifiable against RefSeq (step 13). The secondary metabolism is as")
    say("  uncharacterised as the organisms.")
    say()
    say("  It is a POSITIVE claim about uncharacterised potential, not an absence:")
    say("  antiSMASH found the clusters, and MIBiG simply has nothing like them.")
else:
    say("  A majority of regions resemble something characterised, which is worth")
    say("  checking against the taxonomy — widely-distributed classes (terpenes,")
    say("  quorum-sensing autoinducers) are well represented in MIBiG and will")
    say("  match across large phylogenetic distances.")
say()
say("  LIMITS, all of which belong in the text:")
say("   1. Gene-content similarity, not compound identity (see section 2).")
say("   2. MIBiG is biased toward organisms people culture — Actinobacteria,")
say("      Bacillus, Pseudomonas, fungi. A low hit rate for anaerobic sediment")
say("      Clostridia partly reflects that bias, not only genuine novelty.")
say("   3. antiSMASH regions are predictions of potential. Nothing here shows a")
say("      compound is made, let alone under lagoon conditions.")
say("   4. Archaeal regions were called under bacterial rules and their terpene")
say("      calls overlap ether-lipid biosynthesis (step 33).")

open(os.path.join(B, 'NOVELTY.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(B, 'mibig_hits.tsv')}")
print(f"  written: {os.path.join(B, 'NOVELTY.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "novelty analysis completed" "$RC" || OK=1
check "mibig_hits.tsv written" $([ -s "$B/mibig_hits.tsv" ] && echo 0 || echo 1) || OK=1
NR=$(($(wc -l < "$B/mibig_hits.tsv" 2>/dev/null || echo 1) - 1))
NG=$(ls -1 "$B"/as/*/*.region*.gbk 2>/dev/null | wc -l)
check "one row per region ($NR/$NG)" $([ "$NR" -eq "$NG" ] && echo 0 || echo 1) || OK=1
[ "$OK" -eq 0 ] && runlog "33b_mibig complete" || runlog "**33b_mibig FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
