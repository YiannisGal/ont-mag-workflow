#!/usr/bin/env bash
# =============================================================================
# 25_final_catalogue.sh  --  one row per genome, everything this project knows
#
# WHAT THIS IS
#   The deliverable. Twenty-four analysis steps produced tables scattered across
#   as many directories; this joins them on the MAG identifier and emits the
#   table that goes into the thesis, plus a per-genome narrative paragraph.
#
# DESIGN RULE: JOIN ON THE PUBLISHED ID, NEVER ON FILENAME OR DIRECTORY ORDER.
#   This is the specific failure that produced workspace_1's three conflicting
#   MAG sets (38 / 22 / 17) with no way to reconcile them — audit finding F1.
#   Every input here is matched by mag_id or working_name looked up through
#   MAG_TAXONOMY.tsv, and any input that cannot be matched is REPORTED as
#   unmatched rather than silently dropped.
#
# EVERY COLUMN CARRIES ITS EVIDENCE LEVEL
#   A catalogue that mixes a CheckM2 completeness with a keyword-matched CAZyme
#   count, and prints them in adjacent columns at the same weight, invites the
#   reader to trust both equally. Columns are therefore grouped and the header
#   block states, for each group, what produced it and how far it can be taken.
#
# WHAT IT WILL NOT DO
#   It will not invent a value. A missing input becomes an empty cell and a line
#   in the MISSING report — never a zero, because a zero is a measurement.
#
# INPUT   whatever exists of:
#   07_mags/MAG_TAXONOMY.tsv, MAG_CROSSWALK.tsv, coverage/
#   11_kofam/ko_modules.tsv          12_phylogeny/derep_clusters.tsv
#   15_marker_trees/*/               18_viruses/
#   19b_auxotrophy/                  20_substrates/substrate_matrix.tsv
#   21_network/nodes.tsv             23_extras/{rrna_copy_number,genome_stats}.tsv
#   23_extras/pul/pul_summary.tsv    22_singlem/
#
# OUTPUT  12_results/
#   MAG_CATALOGUE.tsv    the full wide table
#   MAG_CATALOGUE.md     thesis-formatted, grouped, with the evidence key
#   MAG_PROFILES.md      a paragraph per genome
#   MISSING.txt          every input that was absent or unmatched
#
# RUNTIME  under a minute.
# USAGE    ws2 && bash scripts/25_final_catalogue.sh 2>&1 | tee logs/25_catalogue.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="25_catalogue"
R="$WS2/12_results"; mkdir -p "$R"

echo "######## 25 — FINAL MAG CATALOGUE ########"
echo "  started: $(date -Is)"
runlog "25_final_catalogue started"

python3 - "$WS2" "$R" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter

WS2, R = sys.argv[1], sys.argv[2]
missing = []

def readtsv(path, label):
    p = os.path.join(WS2, path)
    if not os.path.exists(p):
        missing.append(f"{label}: NOT FOUND at {path}")
        return []
    try:
        rows = list(csv.DictReader(open(p), delimiter='\t'))
        print(f"  [ok]   {label:<26} {len(rows):>5} rows   {path}")
        return rows
    except Exception as e:
        missing.append(f"{label}: unreadable ({e})")
        return []

print("\n  INPUTS")
tax_rows = readtsv('07_mags/MAG_TAXONOMY.tsv', 'taxonomy')
if not tax_rows:
    raise SystemExit("  !! MAG_TAXONOMY.tsv is the join key for everything. Cannot continue.")

# ---- the spine: one record per MAG, keyed BOTH ways
mags, by_wn, by_id = {}, {}, {}
for r in tax_rows:
    wn = (r.get('working_name') or '').strip()
    mi = (r.get('mag_id') or wn).strip()
    if not mi: continue
    rec = dict(mag_id=mi, working_name=wn)
    # Column names vary between the tables this project produced, so alias them
    # instead of assuming. n_trna / rrna / n_contigs / msa_pct all came back
    # 0/29 in the first build purely because the header did not match.
    ALIAS = {
      'sample':('sample','Sample'), 'domain':('domain','Domain'),
      'phylum':('phylum',), 'class':('class',), 'order':('order',),
      'family':('family',), 'genus':('genus',), 'species':('species',),
      'completeness':('completeness','Completeness','compl'),
      'contamination':('contamination','Contamination','contam'),
      'mimag':('mimag','MIMAG','mimag_quality','quality','tier'),
      'gunc_css':('gunc_css','CSS','clade_separation_score','css'),
      # real header names, read from MAG_TAXONOMY.tsv rather than guessed
      'n_trna':('n_trna_types','n_trna','tRNA','trna','n_tRNA','tRNAs'),
      'has_5S':('has_5S',), 'has_16S':('has_16S',), 'has_23S':('has_23S',),
      'n_contigs':('n_contigs','contigs','num_contigs','n_scaffolds'),
      'bp':('bp','total_bp','length','genome_size','size_bp'),
      'classification_method':('classification_method','method'),
      'msa_pct':('msa_percent','msa_pct','msa','aln_pct'),
      'red_value':('red_value',), 'novelty':('novelty',),
      'deepest_rank':('deepest_rank',), 'gtdb_note':('gtdb_note',),
      'closest_ref':('closest_ref','fastani_reference','closest_reference',
                     'closest_placement_reference'),
      'ani':('ani','fastani_ani','closest_placement_ani','skani_ani'),
    }
    lower = {k.lower(): v for k, v in r.items()}
    for dst, srcs in ALIAS.items():
        for c in srcs:
            v = r.get(c)
            if v in (None, ''): v = lower.get(c.lower())
            if v not in (None, ''): rec[dst] = v; break
    # MIMAG asks for 5S+16S+23S, so collapse the three flags into one column
    def _yes(v): return str(v).strip().lower() in ('1','y','yes','true','t')
    present = [n for n, k in (('5S','has_5S'), ('16S','has_16S'), ('23S','has_23S'))
               if _yes(rec.get(k))]
    if any(k in rec for k in ('has_5S','has_16S','has_23S')):
        rec['rrna'] = '+'.join(present) if present else 'none'
    mags[mi] = rec; by_wn[wn] = mi; by_id[mi] = mi

def _norm(x):
    return re.sub(r'\.(fa|fasta|fna|gff|tsv)$', '', (x or '').strip(), flags=re.I)

# case-insensitive lookup, so CoverM's "Genome" column and a lowercase
# working_name both resolve. The first version matched case-sensitively and
# silently failed to join 31 coverage rows.
by_id_l = {k.lower(): v for k, v in by_id.items()}
by_wn_l = {k.lower(): v for k, v in by_wn.items()}

def key(x):
    x = _norm(x)
    if not x: return None
    return (by_id.get(x) or by_wn.get(x)
            or by_id_l.get(x.lower()) or by_wn_l.get(x.lower()))

def key_anycol(row):
    """Last resort: scan EVERY field for something that resolves to a MAG.
    Guessing column names is what lost the coverage tables; this does not
    guess."""
    for v in row.values():
        k = key(v)
        if k: return k
    return None

def merge(rows, label, fields,
          idcols=('mag_id','working_name','genome','Genome','bin','Bin','MAG',
                  'genome_name','Genome name','name','sample','Sample')):
    """Join by identifier. Anything that will not match is REPORTED."""
    hit = miss = 0
    for r in rows:
        k = None
        for c in idcols:
            if r.get(c):
                k = key(r[c])
                if k: break
        if not k: k = key_anycol(r)
        if not k:
            miss += 1; continue
        hit += 1
        for src, dst in fields.items():
            v = r.get(src)
            if v not in (None, ''): mags[k][dst] = v
    if miss:
        missing.append(f"{label}: {miss} row(s) could not be matched to a MAG id")
    print(f"         {label:<26} matched {hit}, unmatched {miss}")

# ---- crosswalk (binner, provenance, size)
merge(readtsv('07_mags/MAG_CROSSWALK.tsv', 'crosswalk'), 'crosswalk',
      {'binner':'binner','das_refined':'das_refined','contigs':'n_contigs',
       'bp':'bp','total_bp':'bp','avg_gene_length':'avg_gene_len',
       'coding_density':'coding_density','raw_bin_file':'source_bin'})

# ---- coverage
# CoverM names its columns after the BAM/sample, e.g. "SED_das_vamb_004.fa Mean"
# or "reads.fastq.gz Covered Fraction", so a fixed field map cannot work either.
# Find the columns by what they CONTAIN, and say which ones were used.
for code in ('SED','LIQ'):
    rows = readtsv(f'07_mags/coverage/{code}_coverage.tsv', f'coverage {code}')
    if not rows:
        continue
    hdr = list(rows[0].keys())
    def pick(*words):
        for h in hdr:
            hl = h.lower()
            if all(w in hl for w in words): return h
        return None
    c_mean = pick('mean') or pick('coverage') or pick('depth')
    c_brdt = pick('covered', 'fraction') or pick('breadth')
    print(f"         coverage {code}: depth column = {c_mean!r}, "
          f"breadth column = {c_brdt!r}")
    if not c_mean:
        missing.append(f"coverage {code}: no depth-like column found in "
                       f"{hdr[:6]}")
    fmap = {}
    if c_mean: fmap[c_mean] = f'cov_{code}'
    if c_brdt: fmap[c_brdt] = f'breadth_{code}'
    if fmap:
        # CoverM emits an "unmapped" pseudo-genome row. It is expected and is
        # not a join failure, so drop it before merging rather than reporting
        # it as a problem every run.
        real = [r for r in rows
                if 'unmapped' not in ' '.join(str(v) for v in r.values()).lower()]
        merge(real, f'coverage {code}', fmap)

# ---- network degree
merge(readtsv('21_network/nodes.tsv', 'network degree'), 'network degree',
      {'degree_crossfeed':'deg_crossfeed','degree_competition':'deg_competition',
       'degree_dependency':'deg_dependency'})

# ---- rRNA operons, genome size, PULs
merge(readtsv('23_extras/rrna_copy_number.tsv', 'rRNA operons'), 'rRNA operons',
      {'n_16S':'n_16S','n_23S':'n_23S','n_5S':'n_5S','operon_floor':'rrna_operons'})
merge(readtsv('23_extras/genome_stats.tsv', 'genome size'), 'genome size',
      {'estimated_bp':'est_genome_bp'})
merge(readtsv('23_extras/pul/pul_summary.tsv', 'PULs'), 'PULs',
      {'n_cazymes':'n_cazymes','n_PULs':'n_PULs','pct_in_PUL':'pct_cazy_in_PUL',
       'pct_in_PUL_upper':'pct_cazy_in_PUL_hi'})
merge(readtsv('23_extras/secretome.tsv', 'secretome'), 'secretome',
      {'pct_secreted':'pct_secreted','secreted_cazymes':'secreted_cazymes'})

# ---- auxotrophy (19b preferred, 19 fallback)
aux_src = '19b_auxotrophy' if os.path.exists(
    os.path.join(WS2, '19b_auxotrophy/pathway_completeness.tsv')) else '19_auxotrophy'
if aux_src == '19_auxotrophy':
    missing.append("auxotrophy: using step 19 (route-unaware). Run 19b — its Ala/B6/"
                   "Tyr calls were marker-set failures.")
auxrows = readtsv(f'{aux_src}/auxotrophies.tsv', f'auxotrophy ({aux_src})')
byaux = defaultdict(list)
for r in auxrows:
    k = key(r.get('mag_id') or r.get('working_name'))
    if k: byaux[k].append(r.get('pathway',''))
for k, v in byaux.items():
    mags[k]['auxotrophies'] = ';'.join(sorted(set(x for x in v if x)))
    mags[k]['n_auxotrophies'] = str(len(set(v)))

# ---- substrate matrix -> a few headline counts
sub = readtsv('20_substrates/substrate_matrix.tsv', 'substrates')
for r in sub:
    k = key(r.get('mag_id') or r.get('working_name'))
    if not k: continue
    for src, dst in (('KW_CAZY_GH','n_GH'), ('KW_CAZY_PL','n_PL'),
                     ('KW_SUS','n_susCD'), ('KW_PEPT_exo','n_peptidase_exo')):
        if r.get(src) not in (None, ''): mags[k][dst] = r[src]

# ---- methanogenesis modules
for r in readtsv('11_kofam/ko_modules.tsv', 'KEGG modules'):
    k = key(r.get('mag_id') or r.get('working_name') or r.get('genome'))
    if not k: continue
    mod = (r.get('module') or '').strip()
    val = r.get('completeness') or r.get('fraction') or r.get('present')
    if mod and val not in (None, ''): mags[k][f'mod_{mod}'] = val

# ---- proviruses
# TWO BUGS in the previous version, both silent, and the second is the reason
# n_proviruses stayed empty even though 18_viruses/provirus_host_links.tsv
# exists and the glob was finding it:
#   (1) the counter was keyed by CONTIG name, which never matches a MAG key;
#   (2) the counter was never joined into `mags` AT ALL. n_proviruses was
#       declared as a column and printed in the profiles, but no line of code
#       ever assigned it. A column that is only ever read is a column that is
#       always empty, and the completeness block dutifully reported the GAP
#       without being able to say it was self-inflicted.
#
# Fixing (1) needs a contig -> MAG map. Two independent sources are used and
# any contig name that resolves to more than one MAG is DISCARDED rather than
# arbitrated — Bakta renames sequences unless told otherwise, so `contig_1`
# can legitimately appear in every genome, and silently picking one would
# invent provirus assignments.
prov = Counter()
contig2mag, _amb = {}, set()

def _addc(c, k):
    c = (c or '').strip()
    if not c or c in _amb: return
    if c in contig2mag and contig2mag[c] != k:
        del contig2mag[c]; _amb.add(c); return
    contig2mag[c] = k

# primary source: the MAG fastas — their headers are the true contig names
for p in glob.glob(os.path.join(WS2, '07_mags', '**', '*.f*a'), recursive=True):
    k = key(os.path.splitext(os.path.basename(p))[0])
    if not k: continue
    try:
        for line in open(p):
            if line.startswith('>'): _addc(line[1:].split()[0], k)
    except Exception: pass
_from_fasta = len(contig2mag)

# secondary source: Bakta's per-genome .tsv, field 0 = sequence id
for d in sorted(glob.glob(os.path.join(WS2, '09_annotation', '*'))):
    if not os.path.isdir(d): continue
    n = os.path.basename(d); k = key(n)
    if not k: continue
    t = os.path.join(d, f'{n}.tsv')
    if not os.path.exists(t): continue
    try:
        for line in open(t):
            if line.startswith('#') or not line.strip(): continue
            _addc(line.split('\t')[0], k)
    except Exception: pass

print(f"  contig map: {len(contig2mag)} contigs -> MAG "
      f"({_from_fasta} from fastas, {len(_amb)} ambiguous names discarded)")

def _tomag(name):
    """geNomad writes 'contig|provirus_123_4567' and sometimes a bare
    coordinate suffix. Strip progressively and test against the map."""
    if not name: return None
    c = str(name).split('|')[0].split()[0]
    if c in contig2mag: return contig2mag[c]
    return contig2mag.get(re.sub(r'[_:]\d+[-_]\d+$', '', c))

# ---- BUG (3), found by reading the previous run's own output.
# The glob `*provirus*.tsv` matched EVERY geNomad stage table, and geNomad
# emits one row per provirus in each of them:
#     LIQ_asm_final_provirus.tsv                          76 rows
#     LIQ_asm_final_provirus_aggregated_classification    76 rows
#     LIQ_asm_final_provirus_taxonomy                     76 rows
#     LIQ_asm_final_provirus_features                     76 rows
#     LIQ_asm_final_provirus_marker_classification        76 rows
#     LIQ_asm_final_provirus_nn_classification            76 rows
#     LIQ_asm_final_virus_summary                         76 rows
# The SAME 76 proviruses, seven times, all summed into one Counter. The set of
# MAGs carrying a provirus came out right (18/29, matching step 18) because
# multiplying every count by seven does not change which counts are non-zero —
# which is exactly why the error was invisible in the completeness block. The
# per-genome NUMBERS were inflated ~7x.
#
# Two rules from this:
#   - A glob is a claim about which files exist. `*provirus*` claimed one file
#     per run and matched eight.
#   - Counting rows across files needs an identity for the thing being counted.
#     A Counter over rows counts rows; a set over provirus identifiers counts
#     proviruses.
#
# Preference order: step 18's own product first, geNomad only as a fallback,
# and never both.
prov_ids = defaultdict(set)          # mag -> {provirus identity}
link = os.path.join(WS2, '18_viruses', 'provirus_host_links.tsv')
src = None
if os.path.exists(link) and os.path.getsize(link) > 0:
    src = 'provirus_host_links.tsv'
    nrow = 0
    try:
        for r in csv.DictReader(open(link), delimiter='\t'):
            nrow += 1
            k = None
            for c in ('mag_id', 'mag', 'genome', 'bin', 'host', 'host_mag'):
                if r.get(c):
                    k = key(r[c]) or _tomag(r[c])
                    if k: break
            if not k:
                for v in r.values():
                    k = key(v) or _tomag(v)
                    if k: break
            if not k: continue
            ident = None
            for c in ('provirus', 'provirus_id', 'virus_id', 'seq_name',
                      'contig', 'source_seq', 'scaffold', 'coordinates'):
                if r.get(c): ident = f"{c}={r[c]}"; break
            prov_ids[k].add(ident or tuple(sorted(
                (a, b) for a, b in r.items() if isinstance(b, str))))
        ndd = sum(len(v) for v in prov_ids.values())
        print(f"  [ok]   proviruses {src:<34} {nrow} rows -> {ndd} distinct "
              f"across {len(prov_ids)} MAGs")
        if nrow and ndd < nrow:
            print(f"         {nrow - ndd} duplicate row(s) collapsed — this file "
                  f"is not one row per provirus")
    except Exception as e:
        missing.append(f"provirus_host_links.tsv unreadable ({e})")
        prov_ids.clear(); src = None

if src is None:
    # FALLBACK: geNomad summaries only. Explicitly NOT *provirus*.tsv — that
    # glob is what caused the 7x overcount — and explicitly not *_genes.tsv,
    # which is gene-level and whose ragged rows also broke csv.DictReader.
    vfiles = sorted(p for p in glob.glob(
        os.path.join(WS2, '18_viruses', '**', '*_virus_summary.tsv'),
        recursive=True) if '_genes' not in os.path.basename(p))
    for p in vfiles:
        hit = tot = 0
        try:
            for r in csv.DictReader(open(p), delimiter='\t'):
                if 'provirus' not in (r.get('topology') or
                                      r.get('Topology') or '').lower():
                    continue
                tot += 1
                name = (r.get('seq_name') or r.get('contig') or
                        r.get('source_seq') or '')
                k = key(name) or _tomag(name)
                if not k:
                    for v in r.values():
                        k = key(v) or _tomag(v)
                        if k: break
                if k:
                    prov_ids[k].add(f"{os.path.basename(p)}::{name}")
                    hit += 1
        except Exception as e:
            missing.append(f"provirus table {os.path.basename(p)}: unreadable ({e})")
            continue
        print(f"  [ok]   proviruses {os.path.basename(p):<34} "
              f"{tot} provirus rows, {hit} assigned to a MAG")
    if not vfiles:
        missing.append("proviruses: no *_virus_summary.tsv under 18_viruses/ "
                       "and no provirus_host_links.tsv — run 18_viruses_mges.sh")
    src = 'geNomad *_virus_summary.tsv'

prov = Counter({k: len(v) for k, v in prov_ids.items()})
if prov:
    print(f"  proviruses: source = {src}; "
          f"{sum(prov.values())} proviruses across {len(prov)} MAGs "
          f"(max {max(prov.values())} in one genome)")
else:
    missing.append("proviruses: source found but no row resolved to a MAG. May "
                   "be REAL — proviruses on contigs that were never binned "
                   "belong to no MAG — but check the contig map line above.")

# ---- and the fix for bug (2): actually write the column
for k, n in prov.items():
    mags[k]['n_proviruses'] = n
# ---- marker-tree verdicts
for p in glob.glob(os.path.join(WS2, '15_marker_trees', '*', '*verdict*.tsv')):
    marker = os.path.basename(os.path.dirname(p))
    for r in csv.DictReader(open(p), delimiter='\t'):
        k = None
        for c in ('mag_id','working_name','genome','query'):
            if r.get(c):
                k = key(re.split(r'[|_]', r[c])[0]) or key(r[c])
                if k: break
        v = r.get('verdict') or r.get('call')
        if k and v: mags[k][f'tree_{marker}'] = v.split('!')[0].strip()

# ---- dereplication cluster
for r in readtsv('12_phylogeny/derep_clusters.tsv', 'derep clusters'):
    k = key(r.get('mag_id') or r.get('working_name') or r.get('genome'))
    if k and r.get('cluster'): mags[k]['derep_cluster'] = r['cluster']

# =====================================================================
# WRITE
# =====================================================================
GROUPS = [
 ("IDENTITY",   ['mag_id','working_name','sample','binner','das_refined','derep_cluster'],
  "assigned at step 08d; binner provenance from the DAS_Tool crosswalk"),
 ("QUALITY",    ['completeness','contamination','mimag','gunc_css','n_trna','rrna',
                 'n_contigs','bp','est_genome_bp','avg_gene_len','coding_density'],
  "CheckM2 1.1.0 completeness/contamination; GUNC 1.0.3 chimerism (CSS); MIMAG "
  "tier from barrnap + tRNAscan-SE. est_genome_bp = observed / completeness."),
 ("TAXONOMY",   ['domain','phylum','class','order','family','genus','species',
                 'classification_method','msa_pct','red_value','closest_ref','ani',
                 'novelty','deepest_rank'],
  "GTDB-Tk 2.7.2 against GTDB R232 (skani, not FastANI). Placements below ~50% "
  "MSA are provisional."),
 ("ABUNDANCE",  ['cov_SED','breadth_SED','cov_LIQ','breadth_LIQ','rrna_operons',
                 'n_16S','n_23S','n_5S'],
  "CoverM 0.7.0 within each sample only. rRNA operon counts are FLOORS — operons "
  "collapse during assembly — and are trustworthy only above 90% completeness."),
 ("METABOLISM", [],   # filled below with mod_* columns
  "KOfam HMMs with per-family adaptive thresholds (step 11b). Bakta alone was "
  "shown to be unreliable here: it found zero methanogenesis genes in eight "
  "confirmed methanogens and four false dsrA."),
 ("VERIFIED",   [],   # filled below with tree_* columns
  "PHYLOGENETIC verdicts, not HMM hits. Patristic distance to curated reference "
  "sets. These are the only marker calls that survived independent checking."),
 ("SUBSTRATE",  ['n_GH','n_PL','n_susCD','n_cazymes','n_PULs','pct_cazy_in_PUL',
                 'pct_cazy_in_PUL_hi','n_peptidase_exo','pct_secreted','secreted_cazymes'],
  "CAZyme and peptidase counts are Bakta PRODUCT-STRING matches — weaker than KO "
  "evidence, useful for comparison between genomes, not as absolute counts. "
  "pct_cazy_in_PUL is genomic context (susC/susD adjacency), not localisation."),
 ("INTERACTION",['n_auxotrophies','auxotrophies','deg_crossfeed','deg_competition',
                 'deg_dependency','n_proviruses'],
  "Auxotrophy from step 19b: scored per alternative route, >=3 markers, only in "
  "genomes >=90% complete. Network degrees are POTENTIAL interactions from gene "
  "content, never measured fluxes."),
]
allkeys = set()
for m in mags.values(): allkeys |= set(m)
mod_cols  = sorted(k for k in allkeys if k.startswith('mod_'))
tree_cols = sorted(k for k in allkeys if k.startswith('tree_'))
GROUPS[4] = (GROUPS[4][0], mod_cols,  GROUPS[4][2])
GROUPS[5] = (GROUPS[5][0], tree_cols, GROUPS[5][2])

cols = []
for _, cs, _ in GROUPS:
    for c in cs:
        if c not in cols: cols.append(c)

order = sorted(mags.values(), key=lambda r: (r.get('sample','') != 'SED', r['mag_id']))
with open(os.path.join(R, 'MAG_CATALOGUE.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    for r in order: w.writerow(r)
print(f"\n  written: {os.path.join(R,'MAG_CATALOGUE.tsv')}  "
      f"({len(order)} genomes x {len(cols)} columns)")

# ---- markdown, grouped, with the evidence key
with open(os.path.join(R, 'MAG_CATALOGUE.md'), 'w') as fh:
    fh.write("# Etoliko Lagoon MAG catalogue\n\n")
    fh.write(f"{len(order)} metagenome-assembled genomes "
             f"({sum(1 for r in order if r.get('sample')=='SED')} sediment, "
             f"{sum(1 for r in order if r.get('sample')=='LIQ')} enrichment culture).\n\n")
    fh.write("## How to read this table\n\n")
    fh.write("Columns are grouped by what produced them. **The groups are not "
             "equally strong** and should not be quoted as though they were.\n\n")
    for name, cs, note in GROUPS:
        if not cs: continue
        fh.write(f"- **{name}** — {note}\n")
    fh.write("\n")
    for name, cs, note in GROUPS:
        cs = [c for c in cs if any(r.get(c) for r in order)]
        if not cs: continue
        fh.write(f"\n## {name}\n\n{note}\n\n")
        fh.write('| MAG | ' + ' | '.join(cs) + ' |\n')
        fh.write('|' + '---|' * (len(cs) + 1) + '\n')
        for r in order:
            fh.write(f"| {r['mag_id']} | " +
                     ' | '.join(str(r.get(c, '')) for c in cs) + ' |\n')
print(f"  written: {os.path.join(R,'MAG_CATALOGUE.md')}")

# ---- one paragraph per genome
def para(r):
    bits = []
    t = ' '.join(x for x in (r.get('genus'), r.get('species')) if x) or \
        (r.get('family') or r.get('phylum') or 'unclassified')
    src = 'the sediment' if r.get('sample') == 'SED' else 'an enrichment culture'
    bits.append(f"**{r['mag_id']}** ({t}) was recovered from {src} by "
                f"{r.get('binner','an unrecorded binner')}")
    q = []
    if r.get('completeness'): q.append(f"{float(r['completeness']):.1f}% complete")
    if r.get('contamination'): q.append(f"{float(r['contamination']):.1f}% contaminated")
    if r.get('mimag'): q.append(f"MIMAG {r['mimag']}")
    if q: bits.append(', '.join(q))
    if r.get('est_genome_bp'):
        bits.append(f"estimated genome size {int(r['est_genome_bp'])/1e6:.2f} Mb")
    if r.get('rrna_operons'):
        bits.append(f"at least {r['rrna_operons']} rRNA operon(s)")
    s = '; '.join(bits) + '.'
    ex = []
    for c in tree_cols:
        if str(r.get(c, '')).upper().startswith('CONFIRM'):
            ex.append(c.replace('tree_', '') + ' phylogenetically confirmed')
    if r.get('n_PULs') and r['n_PULs'] not in ('0', 0):
        ex.append(f"{r['n_PULs']} polysaccharide utilisation loci "
                  f"({r.get('pct_cazy_in_PUL','?')}% of its CAZymes)")
    if r.get('n_proviruses'): ex.append(f"{r['n_proviruses']} predicted provirus(es)")
    if r.get('n_auxotrophies'):
        ex.append(f"{r['n_auxotrophies']} firm auxotrophy(ies): "
                  f"{r.get('auxotrophies','')}")
    if ex: s += ' ' + '. '.join(x[0].upper() + x[1:] for x in ex) + '.'
    return s

with open(os.path.join(R, 'MAG_PROFILES.md'), 'w') as fh:
    fh.write("# Per-genome profiles\n\n")
    fh.write("Generated from the catalogue. Every statement traces to a column; "
             "nothing here is inferred beyond what the tables contain.\n\n")
    for r in order:
        fh.write(para(r) + "\n\n")
print(f"  written: {os.path.join(R,'MAG_PROFILES.md')}")

# ---- what is missing
with open(os.path.join(R, 'MISSING.txt'), 'w') as fh:
    fh.write("Inputs absent or unmatched when the catalogue was built.\n")
    fh.write("An empty cell in the catalogue means NOT MEASURED, never zero.\n\n")
    for m in missing: fh.write(m + "\n")
print()
print("=" * 96)
print("  COMPLETENESS OF THE CATALOGUE ITSELF")
print("=" * 96)
EXPECTED_PARTIAL = {
  'cov_SED':'only the 4 sediment MAGs have sediment coverage',
  'breadth_SED':'as cov_SED',
  'cov_LIQ':'only the 25 culture MAGs have culture coverage',
  'breadth_LIQ':'as cov_LIQ',
  'species':'only genomes with a GTDB species-level match — 22/29 are novel',
  'genus':'two MAGs are novel genera and have none',
  'n_auxotrophies':'19b calls auxotrophies only in genomes >=90% complete',
  'auxotrophies':'as n_auxotrophies',
  'n_proviruses':'step 18 found proviruses in 18 of 29 MAGs; the other 11 '
                 'genuinely carry none. A blank here is a result, not a gap',
}
for c in cols:
    n = sum(1 for r in order if r.get(c) not in (None, ''))
    if n < len(order):
        note = EXPECTED_PARTIAL.get(c)
        tag = f"   expected — {note}" if note else "   <-- GAP"
        print(f"    {c:<24}{n:>3}/{len(order)} populated{tag}")
if missing:
    print(f"\n  {len(missing)} input issue(s) — see MISSING.txt:")
    for m in missing: print("    - " + m)
else:
    print("\n  every expected input was found and matched.")
PY

echo
echo "######## VALIDATION GATE ########"
check "catalogue TSV"  $([ -s "$R/MAG_CATALOGUE.tsv" ] && echo 0 || echo 1)
check "catalogue MD"   $([ -s "$R/MAG_CATALOGUE.md" ] && echo 0 || echo 1)
check "profiles"       $([ -s "$R/MAG_PROFILES.md" ] && echo 0 || echo 1)
NR=$(( $(wc -l < "$R/MAG_CATALOGUE.tsv" 2>/dev/null || echo 1) - 1 ))
check "29 genomes"     $([ "$NR" -eq 29 ] && echo 0 || echo 1)

runlog "25_final_catalogue complete"
echo
echo "######## 25 DONE ########"
