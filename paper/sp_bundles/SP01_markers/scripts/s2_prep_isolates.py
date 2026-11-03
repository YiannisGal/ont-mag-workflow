"""SP01 step 2a -- concatenate the five Prokka proteomes as >ISOLATE|locus_tag; record md5s.
One Prokka .faa per isolate (Functional_<id>/). The Manuscript_Analysis_MASTER copies are byte-identical;
every other .faa under UniqueELT_Solid is a Prodigal/BUSCO/GTDB-Tk/CheckM2 call, not Prokka."""
import os, hashlib
OUT, B = os.environ['OUT'], os.environ['B']
MA = f'{B}/UniqueELT_Solid/Manuscript_Analysis'
ISO = {  # SP01 name -> Prokka faa
    'YGCA01_C':  f'{MA}/Functional_YGCA01_C/YGCA01_C.faa',
    'YGCA01_Tp': f'{MA}/Functional_YGCA01_Contam/YGCA01_Contam.faa',
    'YGCA07':    f'{MA}/Functional_YGCA07/YGCA07.faa',
    'YGCA10':    f'{MA}/Functional_YGCA10/YGCA10.faa',
    'YGCC01':    f'{MA}/Functional_YGCC01/YGCC01.faa',
}
with open(f'{OUT}/isolates/isolates.faa', 'w') as out, open(f'{OUT}/isolates/inputs.tsv', 'w') as tab:
    tab.write('genome\tfaa\tmd5\tn_proteins\n')
    for g, p in ISO.items():
        n = 0
        for line in open(p):
            if line.startswith('>'):
                n += 1; out.write(f'>{g}|{line[1:].split()[0]}\n')
            else: out.write(line)
        tab.write(f'{g}\t{p}\t{hashlib.md5(open(p,"rb").read()).hexdigest()}\t{n}\n')
        print(g, n)
