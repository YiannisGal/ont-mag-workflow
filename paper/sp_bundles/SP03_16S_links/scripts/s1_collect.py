"""SP03 step 1 -- every 16S copy of the 29 MAGs and 5 isolates via barrnap (kingdom by GTDB domain), plus inputs md5.
Header: >GENOME|cN|contig:start-end(strand)|len|partial_flag
The supplied inputs/all_strains_16S.fasta contains YGCC01 plus Sporosarcina REFERENCE genomes (GCA_*, S_ureae, S_psychrophila),
not the other four isolates, so isolates are taken from barrnap on their own assemblies; the supplied YGCC01 copy is cross-checked."""
import os, csv, subprocess, hashlib, glob
B, WS2, OUT = os.environ['B'], os.environ['WS2'], os.environ['OUT']
MA = f'{B}/UniqueELT_Solid/Manuscript_Analysis'
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
tax = {r['mag_id']: r for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_TAXONOMY.tsv'), delimiter='\t')}
G = []
for m, r in sorted(tax.items()):
    G.append((m, 'MAG', f"{WS2}/07_mags/bins/{r['working_name']}.fa", 'arc' if r['domain'] == 'Archaea' else 'bac', r['genus'] or r['family']))
for g, d in (('YGCA01_C', 'YGCA01_C'), ('YGCA01_Tp', 'YGCA01_Contam'), ('YGCA07', 'YGCA07'), ('YGCA10', 'YGCA10'), ('YGCC01', 'YGCC01')):
    G.append((g, 'isolate', f'{MA}/Functional_{d}/{d}.fna', 'bac', ''))
def rc(s):
    return s[::-1].translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))
seqs = {}
with open(f'{OUT}/rrna/genome_16S.fa', 'w') as o, open(f'{OUT}/rrna/genomes.tsv', 'w') as t, open(f'{OUT}/inputs_md5.tsv', 'w') as mt:
    t.write('genome\tsource\tdomain_kingdom\ttaxon\tfasta\tmd5\tn_16S\tn_16S_full\tcopies\n')
    mt.write('file\tmd5\n')
    for g, src, fa, k, taxon in G:
        gff = f'{OUT}/rrna/{g}.gff'; fout = f'{OUT}/rrna/{g}.rrna.fa'
        cmd = ['barrnap', '--kingdom', k, '--threads', os.environ['CPU'], '--outseq', fout, fa]
        print('CMD', ' '.join(cmd), flush=True)
        with open(gff, 'w') as gh: subprocess.run(cmd, check=True, stdout=gh, stderr=subprocess.DEVNULL)
        n = full = 0; copies = []; name = None; cur = {}
        for l in open(fout):
            if l.startswith('>'): name = l[1:].strip(); cur[name] = []
            else: cur[name].append(l.strip())
        # barrnap 1.10 (Infernal) headers are contig:start0-end(strand) with no gene name: take 16S from the GFF
        for gl in open(gff):
            if gl.startswith('#'): continue
            f = gl.rstrip('\n').split('\t')
            if len(f) < 9 or 'Name=16S_rRNA' not in f[8]: continue
            key = f'{f[0]}:{int(f[3]) - 1}-{f[4]}({f[6]})'
            s = ''.join(cur[key]); n += 1; partial = 'partial' in f[8].lower()
            full += not partial
            hid = f'{g}|c{n}|{key}|{len(s)}|{"partial" if partial else "full"}'
            o.write(f'>{hid}\n{s}\n'); seqs[hid] = s; copies.append(f'c{n}:{len(s)}{"p" if partial else ""}')
        t.write(f"{g}\t{src}\t{k}\t{taxon}\t{fa}\t{md5(fa)}\t{n}\t{full}\t{','.join(copies) or '-'}\n")
        mt.write(f'{fa}\t{md5(fa)}\n')
    for p in glob.glob(f'{OUT}/../inputs/*.fasta'): mt.write(f'{os.path.realpath(p)}\t{md5(p)}\n')
    for p in [f'{B}/Etoliko_138_2/Bacteria/out_3/otus.fa', f'{B}/Etoliko_138_2/Bacteria/out_3/otutab_norm.txt',
              f'{B}/Etoliko_138_2/Bacteria/out_3/taxonomy_138_2/taxonomy.tsv',
              f'{B}/Etoliko_138_2/Bacteria/out_3/qplot_input_fixed/Bacteria Metadata.txt',
              f'{B}/Etoliko_138_2/Archaea/out/otus.fa', f'{B}/Etoliko_138_2/Archaea/out/otutab_norm.txt',
              f'{B}/Etoliko_138_2/Archaea/out/taxonomy_138_2/taxonomy.tsv',
              f'{B}/Etoliko_138_2/ONT/final_output/cluster_seq_trimmed.fa', f'{B}/Etoliko_138_2/ONT/final_output/cluster_tab_trimmed.txt',
              f'{B}/Etoliko_138_2/ONT/final_output/taxonomy.tsv', f'{B}/Etoliko_138_2/ONT/final_output/Metadata_minion.txt']:
        mt.write(f'{p}\t{md5(p)}\n')
# positive control
open(f'{OUT}/rrna/control_Arc_Otu22.fa', 'w').write(open(f'{OUT}/../inputs/Arc_Otu22_fulllength_16S.fasta').read())
# cross-check the supplied YGCC01 16S against barrnap's copies
sup = ''.join(l.strip() for l in open(f'{OUT}/../inputs/YGCC01_16S.fasta') if not l.startswith('>')).upper()
hit = [h for h, s in seqs.items() if h.startswith('YGCC01|') and (s.upper() == sup or rc(s.upper()) == sup)]
print('supplied YGCC01_16S.fasta identical to barrnap copy:', hit or 'NONE (reported, not used)')
# table03 claim: 7 MAGs lack 16S
t3 = {r['MAG']: r['16S'] for r in csv.DictReader(open(f'{WS2}/tables/table03_mag_catalogue.tsv'), delimiter='\t')}
ours = {l.split('\t')[0]: int(l.split('\t')[6]) for l in list(open(f'{OUT}/rrna/genomes.tsv'))[1:]}
mism = [(m, t3.get(m), ours[m]) for m in t3 if (t3[m] == 'yes') != (ours.get(m, 0) > 0)]
print('table03 "16S" column vs barrnap now: mismatches', mism or 'none')
print('MAGs without 16S:', sorted(m for m in ours if m.startswith('ETO') and ours[m] == 0))
