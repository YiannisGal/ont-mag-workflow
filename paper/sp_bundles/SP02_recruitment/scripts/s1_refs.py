"""SP02 step 1 -- two competitive references with genome~contig names; per-genome lengths; md5 of every input.
Index A: 29 MAGs + YGCA01_Clean + YGCA01_Contaminant + YGCA10 + YGCC01 + background.
Index B: identical, with YGCA07 instead of YGCA10.
Isolates: Manuscript_Analysis/Final_Genomes/*.fasta; every copy on the server is md5-identical to
Manuscript_Analysis_MASTER/important_files/genomes_final (checked here). The mixed YGCA01 file is not used.
Background (optional in the spec): the 11 GCA accessions in all_12_genomes.txt (the 12th entry is YGCC01 itself).
That list covers the YGCC01 comparison set only; no Clostridium_J/Terrisporobacter/Paraclostridium backgrounds exist in it."""
import os, csv, hashlib
B, WS2, OUT = os.environ['B'], os.environ['WS2'], os.environ['OUT']
MA = f'{B}/UniqueELT_Solid/Manuscript_Analysis'
GF = f'{B}/UniqueELT_Solid/Manuscript_Analysis_MASTER/important_files/genomes_final'
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
os.makedirs(f'{OUT}/refs', exist_ok=True)
G = []
for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t'):
    G.append((r['mag_id'], 'MAG', f"{WS2}/07_mags/bins/{r['working_name']}.fa"))
iso = {}
for g in ('YGCA01_Clean', 'YGCA01_Contaminant', 'YGCA10', 'YGCC01', 'YGCA07'):
    p = f'{MA}/Final_Genomes/{g}.fasta'
    assert md5(p) == md5(f'{GF}/{g}.fasta'), f'{g} differs from genomes_final'
    iso[g] = p
bg = []
for l in open(f'{B}/UniqueELT_Solid/Manuscript_Analysis_MASTER/tables/species_comparison_ANI/all_12_genomes.txt'):
    p = l.strip()
    if 'GCA_' in os.path.basename(p): bg.append(('BG_' + '_'.join(os.path.basename(p).split('_')[:2]), 'background', p))
IX = {'A': G + [(g, 'isolate', iso[g]) for g in ('YGCA01_Clean', 'YGCA01_Contaminant', 'YGCA10', 'YGCC01')] + bg,
      'B': G + [(g, 'isolate', iso[g]) for g in ('YGCA01_Clean', 'YGCA01_Contaminant', 'YGCA07', 'YGCC01')] + bg}
with open(f'{OUT}/refs/genomes.tsv', 'w') as t:
    t.write('index\tgenome\trole\tfasta\tmd5\tlength_bp\tn_contigs\n')
    for ix, gl in IX.items():
        with open(f'{OUT}/refs/{ix}.fa.part', 'w') as o:
            for g, role, p in gl:
                L = n = 0
                for l in open(p):
                    if l.startswith('>'):
                        n += 1; o.write(f'>{g}~{l[1:].split()[0]}\n')
                    else:
                        o.write(l); L += len(l.strip())
                t.write(f'{ix}\t{g}\t{role}\t{p}\t{md5(p)}\t{L}\t{n}\n')
        os.rename(f'{OUT}/refs/{ix}.fa.part', f'{OUT}/refs/{ix}.fa')
        print(ix, len(gl), 'genomes', md5(f'{OUT}/refs/{ix}.fa'))
with open(f'{OUT}/refs/reads_md5.tsv', 'w') as t:
    for p in (os.environ['SED'], os.environ['LIQ']): t.write(f'{p}\t{md5(p)}\n')
