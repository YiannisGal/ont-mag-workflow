"""SP05 step 1 -- inputs with md5: the two dsr contigs, every protein called on them (34_unbinned Prodigal gene calls, the set the
findings report counts as 26), step-58 dsrAB references, and ETO_LIQ_MAG_017's mcrA fragments joined as in step 54
(adjacent same-strand ORFs of one locus concatenated in gene order)."""
import os, hashlib, csv
WS2, OUT = os.environ['WS2'], os.environ['OUT']
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
def rf(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith('>'): n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
C = ('LIQ_contig_12880', 'LIQ_contig_19299')
ctg = rf(f'{WS2}/55_dsr_id/dsr_contigs.fasta')
with open(f'{OUT}/work/contigs.fa', 'w') as o:
    for k, v in ctg.items():
        if k.split()[0] in C: o.write(f'>{k.split()[0]}\n{v}\n')
prot = rf(f'{WS2}/34_unbinned/prot/LIQ_all.faa')
cp = {k.split()[0]: (k, v.rstrip('*')) for k, v in prot.items() if k.split()[0].rsplit('_', 1)[0] in C}
with open(f'{OUT}/work/contig_proteins.faa', 'w') as o:
    for k, (h, v) in sorted(cp.items(), key=lambda x: (x[0].rsplit('_', 1)[0], int(x[0].rsplit('_', 1)[1]))):
        o.write(f'>{k}\n{v}\n')
refs = rf(f'{WS2}/58_dsr_finish/refs.faa')
open(f'{OUT}/work/dsr_refs.faa', 'w').write(''.join(f'>{k.split()[0]}\n{v}\n' for k, v in refs.items()))
mags = rf(f'{WS2}/45_kofam_full/mags_only.faa')
mg = {k.split()[0].split('|')[1]: v.rstrip('*') for k, v in mags.items() if k.startswith('LIQ_das_vamb_004|')}
J = {'MAG017_mcrA_joined_3836_37+38': ('LIQ_contig_3836_37', 'LIQ_contig_3836_38'),
     'MAG017_mcrA_joined_13851_5+6': ('LIQ_contig_13851_5', 'LIQ_contig_13851_6')}
with open(f'{OUT}/work/mag017_mcrA.faa', 'w') as o:
    for name, parts in J.items():
        o.write(f'>{name}\n' + ''.join(mg[p] for p in parts) + '\n')
    for p in sorted({x for ps in J.values() for x in ps}): o.write(f'>MAG017_fragment_{p}\n{mg[p]}\n')
with open(f'{OUT}/work/kofam_query.faa', 'w') as o:
    for f in ('contig_proteins.faa', 'dsr_refs.faa', 'mag017_mcrA.faa'): o.write(open(f'{OUT}/work/{f}').read())
with open(f'{OUT}/inputs_md5.tsv', 'w') as o:
    for p in (f'{WS2}/55_dsr_id/dsr_contigs.fasta', f'{WS2}/34_unbinned/prot/LIQ_all.faa', f'{WS2}/58_dsr_finish/refs.faa',
              f'{WS2}/45_kofam_full/mags_only.faa', f'{WS2}/15_marker_trees/mcrA/mcrA_all.faa', f'{WS2}/15_marker_trees/mcrA/mcrA_verdicts.tsv',
              f'{WS2}/54_mcr_fragmentation/hits_with_coords.tsv', f'{WS2}/07_mags/bins/LIQ_das_vamb_004.fa', f'{WS2}/07_mags/bins/LIQ_das_metabat_022.fa',
              f'{WS2}/58_dsr_finish/coverage.tsv', os.environ['KOLIST']):
        o.write(f'{p}\t{md5(p)}\n')
print('contigs', len([k for k in ctg if k.split()[0] in C]), 'contig proteins', len(cp), 'dsr refs', len(refs))
