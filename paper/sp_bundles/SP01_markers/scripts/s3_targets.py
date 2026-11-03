"""SP01 step 3a -- target KO profiles + one proteome of the 34 genomes (29 MAGs + 5 isolates).
Headers: >GENOME|protein  (GENOME = ETO_* id for MAGs, isolate name for isolates)."""
import os, csv, hashlib
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
PROF = os.environ['PROFILES']

# the spec's KO list; K00577 (the true mtrA) is ADDED because the spec's "mtrA K00584" is mtrH in ko_list
SPEC = [('mcrA','K00399'),('mcrB','K00401'),('mcrG','K00402'),('mtrH(spec:mtrA)','K00584'),('mtrA','K00577'),
        ('mer','K00320'),('mtd','K00319'),('dsrA','K11180'),('dsrB','K11181'),('aprA','K00394'),('aprB','K00395'),
        ('sat','K00958'),('fdhA','K05299'),('fdhA','K22015'),('fdhF','K00123'),('hydA','K00533'),('hydA','K17997'),
        ('echA','K14086'),('nifH','K02588'),('nifD','K02586'),('nifK','K02591'),('ureA','K01430'),('ureB','K01429'),
        ('ureC','K01428'),('ureD','K03190'),('ureE','K03187'),('ureF','K03188'),('ureG','K03189'),('fadA','K00632'),
        ('fadB','K01825'),('fadD','K01897'),('fadE','K00249'),('narG','K00370'),('merA','K00520'),('arcA','K01478'),
        ('arcB','K00611'),('arcC','K00926')]
# M00087 beta-oxidation (KEGG REST, fetched 2026-09-29): all KOs of its definition
M00087 = ['K00232','K00249','K00255','K06445','K09479','K01692','K07511','K13767','K00022','K07516',
          'K01825','K01782','K07514','K07515','K10527','K00632','K07508','K07509','K07513']
kos = []
for g, k in SPEC + [('M00087', k) for k in M00087]:
    if k not in kos: kos.append(k)
os.makedirs(f'{OUT}/markers', exist_ok=True)
with open(f'{OUT}/markers/targets.tsv', 'w') as t:
    t.write('KO\tgene\n')
    for g, k in SPEC: t.write(f'{k}\t{g}\n')
    for k in M00087:
        if k not in dict((kk, 1) for _, kk in SPEC): t.write(f'{k}\tM00087_only\n')
with open(f'{OUT}/markers/targets.hmm', 'w') as o:
    for k in kos:
        o.write(open(f'{PROF}/{k}.hmm').read())
print('target profiles:', len(kos))

xw = list(csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t'))
name2id = {}
for r in xw:
    name2id[r['working_name']] = r['mag_id']; name2id[r['qc_bin_name']] = r['mag_id']
seen = set(); n = 0
with open(f'{OUT}/markers/all34.faa', 'w') as o, open(f'{OUT}/markers/inputs.tsv', 'w') as tab:
    tab.write('genome\tsource\tfile\tmd5\n')
    for line in open(f'{WS2}/45_kofam_full/mags_only.faa'):
        if line.startswith('>'):
            b, p = line[1:].split()[0].split('|', 1)
            mid = name2id[b]; seen.add(mid); n += 1
            o.write(f'>{mid}|{p}\n')
        else: o.write(line)
    # the proteomes step 45 used, by md5
    for l in open(f'{WS2}/45_kofam_full/selected_faa.txt'):
        b, p = l.split()
        tab.write(f'{name2id[b]}\tMAG\t{p}\t{hashlib.md5(open(p,"rb").read()).hexdigest()}\n')
    for l in list(open(f'{OUT}/isolates/inputs.tsv'))[1:]:
        g, p, m, _ = l.rstrip('\n').split('\t'); tab.write(f'{g}\tisolate\t{p}\t{m}\n')
    o.write(open(f'{OUT}/isolates/isolates.faa').read())
assert len(seen) == 29, seen
print('MAG proteins', n, 'genomes', len(seen))
