"""SP01 step 1 -- which genomes were Combined_analysis's "31 MAGs" / "34 genomes", and how it called markers.
Read-only on Combined_analysis and workspace_1. Writes provenance_31.tsv and combined_analysis_marker_method.md."""
import os, glob, hashlib, csv, re, subprocess

OUT, WS2, B = os.environ['OUT'], os.environ['WS2'], os.environ['B']
CA = f'{B}/Combined_analysis/workspace'
WS1 = f'{B}/Metagenome_analysis_PhD/workspace'

def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''): h.update(blk)
    return h.hexdigest()

# ---- reference sets
xw = list(csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t'))
final = {}                                               # md5 -> ETO id
for r in xw:
    p = f"{WS2}/07_mags/bins/{r['working_name']}.fa"
    final[md5(p)] = r['mag_id']
assert len(final) == 29, f'expected 29 final MAG md5s, got {len(final)}'
qc = {r['bin']: r for r in csv.DictReader(open(f'{WS2}/06_binqc/MAG_QC.tsv'), delimiter='\t')}
union = {}                                               # md5 -> union bin name
for p in glob.glob(f'{WS2}/06_binqc/union_bins/*.fa'):
    union[md5(p)] = os.path.basename(p)[:-3]
REJECTED = {'SED_das_metabat_003', 'SED_das_unknown_000a'}
ws1 = {}
for ext in ('*.fa', '*.fasta', '*.fna'):
    for p in glob.glob(f'{WS1}/**/{ext}', recursive=True):
        try:
            if os.path.getsize(p) < 200 * 1024 * 1024: ws1.setdefault(md5(p), os.path.relpath(p, WS1))
        except OSError: pass

# ---- the files to classify
files = sorted(glob.glob(f'{CA}/00_inputs/metagenomes/mags/*')) + sorted(glob.glob(f'{CA}/05_all_31_mags_analysis/*'))
rows = []
for p in files:
    m = md5(p); rel = os.path.relpath(p, CA)
    is_fasta = p.endswith(('.fa', '.fasta', '.fna'))
    if m in final:
        cat, hit = 'a_final_MAG', final[m]
    elif m in union and union[m] in REJECTED:
        cat, hit = 'b_GUNC_rejected_bin', union[m]
    elif m in union:
        cat, hit = 'union_bin_not_final', union[m]
    elif m in ws1:
        cat, hit = 'c_workspace1_bin', ws1[m]
    elif not is_fasta:
        cat, hit = 'not_a_genome(table)', '-'
    else:
        cat, hit = 'd_unknown', '-'
    rows.append((rel, m, cat, hit))

# which 29 final IDs did Combined_analysis's FASTAs NOT cover?
covered = {r[3] for r in rows if r[2] == 'a_final_MAG'}
with open(f'{OUT}/provenance_31.tsv', 'w') as fh:
    fh.write('file\tmd5\tmatch_category\tmatched_final_ID\n')
    for r in rows: fh.write('\t'.join(r) + '\n')
    for mid in sorted(set(final.values()) - covered):
        fh.write(f'(absent from Combined_analysis)\t-\tfinal_MAG_not_in_Combined\t{mid}\n')

# ---- row counts of the Combined tables, to show what "31"/"34"/"36" count
def nrows(p):
    with open(p) as fh: return sum(1 for _ in fh) - 1
def ncols(p):
    with open(p) as fh: return len(fh.readline().rstrip('\n').split('\t')) - 1
cnt = {os.path.relpath(p, CA): (nrows(p), ncols(p)) for p in files if p.endswith('.tsv')}
nfa = sum(1 for r in rows if r[0].endswith('.fa'))

# ---- method: verbatim excerpts
def lines(path, a, b):
    with open(path) as fh: L = fh.readlines()
    return ''.join(L[a - 1:b])
def find(path, pat):
    with open(path) as fh:
        for i, l in enumerate(fh, 1):
            if re.search(pat, l): return i
    raise SystemExit(f'pattern {pat!r} not found in {path}')
s02 = f'{CA}/scripts/02_run_read_recruitment.sh'
s04 = f'{CA}/scripts/04_build_metabolic_matrix.sh'
s06 = f'{CA}/scripts/06_expand_to_all_31_mags.sh'
s11 = f'{WS2}/scripts/11b_kofam_archaea.sh'
a02 = find(s02, 'MAG_Lysinibacillus_sphaericus')
a04, b04 = find(s04, r'^MARKERS = \{'), find(s04, r'^df = pd.DataFrame')
a06, b06 = find(s06, r'^MARKERS = \{'), find(s06, r'^df = pd.DataFrame')
a11 = find(s11, r'run_tool hmmsearch')
KOLIST = os.environ['KOLIST']
defs = {}
for l in open(KOLIST):
    f = l.rstrip('\n').split('\t')
    if len(f) > 11: defs[f[0]] = f[11]
used = {'mcrA':'K00399','mcrB':'K00401','mcrG':'K00402','mtrA':'K00577','mer':'K00320','mtd':'K00319',
        'dsrA':'K11180','dsrB':'K11181','aprA':'K00394','aprB':'K00395','sat':'K00958','fdhA':'K00123',
        'hydA':'K00436','echA':'K14083','ureA':'K01428','ureB':'K01429','ureC':'K01430','nifH':'K02588',
        'nifD':'K02586','nifK':'K02591','fadA':'K00632','fadB':'K01782','fadD':'K01897'}
md = [
 '# How Combined_analysis called its marker matrix (and what "31" / "34" / "36" are)', '',
 'Source tree (read-only): `' + CA + '`', '',
 '## What the numbers count', '',
 f'- `00_inputs/metagenomes/mags/` holds **{nfa} MAG FASTAs**. Classification is in `provenance_31.tsv`.',
]
for k, (r, c) in sorted(cnt.items()):
    md.append(f'- `{k}`: {r} data rows x {c} data columns')
md += ['',
 '"31" appears only in the step name (`06_expand_to_all_31_mags.sh`, `05_all_31_mags_analysis/`), in RUN_05 and in Chapter 6 prose. '
 'No file, table or loop in Combined_analysis enumerates 31 genomes: the matrix is built from `glob(mags/*.fa)` (29) plus 5 hard-coded isolates = 34 columns, '
 'written to a file named `all_36_genomes_...`. The MAG_QC.tsv/master catalog rows are the workspace_2 union-bin QC table, not a MAG set.', '',
 '## Source of the species-style MAG names', '',
 f'`scripts/02_run_read_recruitment.sh` lines {a02}-{a02+3} (verbatim):', '', '```bash', lines(s02, a02, a02 + 3).rstrip(), '```', '',
 'The same mapping is re-typed in `04_build_metabolic_matrix.sh` (`mag_file_map`). These are hand labels on four workspace_2 bins; '
 'they are not GTDB-Tk R232 names from 08d.', '',
 '## Marker calling in 06_expand_to_all_31_mags.sh (the 34-genome matrix)', '',
 f'Verbatim, lines {a06}-{b06}:', '', '```python', lines(s06, a06, b06).rstrip(), '```', '',
 f'04_build_metabolic_matrix.sh (the 9-genome matrix) uses the same MARKERS dict, lines {a04}-{b04}.', '',
 '## The hit tables it read', '',
 f'MAG hits came from `workspace_2/11_kofam/hits/<bin>.tbl`, produced by `scripts/11b_kofam_archaea.sh` line {a11}:', '',
 '```bash', lines(s11, a11, a11 + 1).rstrip(), '```', '',
 '## Method, stated plainly', '',
 '1. **MAGs:** every row of an `hmmsearch -E 1e-5 --tblout` table whose KO is in MARKERS adds 1 to the cell. '
 '**No KOfam curated threshold** and no score-type or length filter. Cells are raw row counts, so any homolog above E=1e-5 counts.',
 '2. **Isolates:** substring match of `gene=<marker>`, `product=<marker>` or ` <marker> ` against each Prokka GFF line. '
 'This is name matching, not KO assignment. It is a different method from the MAGs, so the columns are not comparable.',
 '3. **KO labels are wrong in several rows** (checked against `11_kofam/ko_list`):', '',
 '| label used | KO used | ko_list definition of that KO |', '|---|---|---|',
]
for g, k in used.items():
    md.append(f'| {g} | {k} | {defs.get(k, "?")} |')
md += ['',
 'Label errors visible above: ureA/ureC are swapped (K01428 is ureC, K01430 is ureA); `hydA`=K00436 is the NiFe hoxH, not [FeFe] HydA; '
 '`fadB`=K01782 is the multifunctional fadJ/fadB of E. coli-type bacteria; `echA`=K14083 is the trimethylamine methyltransferase MttB, not the Ech hydrogenase (K14086); '
 '`fdhA`=K00123 is the Fdh alpha-subunit, and without a threshold it also collects NADH:quinone oxidoreductase and other Complex-I/molybdopterin paralogs.', '',
 'Conclusion: the Combined_analysis matrix is presence of any E<1e-5 homolog for MAGs and a name search for isolates. '
 'It is not a curated-threshold KO call and must not be cited as one. SP01 replaces it (`markers_long.tsv`, `markers_wide.tsv`).']
open(f'{OUT}/combined_analysis_marker_method.md', 'w').write('\n'.join(md) + '\n')

from collections import Counter
print('provenance categories:', dict(Counter(r[2] for r in rows)))
print('final MAGs absent from Combined_analysis:', sorted(set(final.values()) - covered))
