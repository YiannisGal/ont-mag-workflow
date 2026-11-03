# Build CGView (classic 1.0) inputs and Proksee upload tables from Bakta GBFF + eggNOG COG.
import sys,os
from Bio import SeqIO
m,gbff,egg,od,pd=sys.argv[1:6]
cog={}
for l in open(egg):
    if l.startswith("#"): continue
    p=l.rstrip("\n").split("\t")
    if len(p)==2 and p[1]: cog[p[0]]=p[1][0]
    elif len(p)>6 and p[6] not in("","-"): cog[p[0]]=p[6][0]
recs=sorted(SeqIO.parse(gbff,"genbank"),key=lambda r:-len(r.seq))
off=0; fa=open(f"{od}/{m}.fasta","w"); seq=[]
cds=open(f"{od}/{m}_cds_cog.tab","w"); rna=open(f"{od}/{m}_rna.tab","w"); ctg=open(f"{od}/{m}_contigs.tab","w")
H="seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n"
for f in(cds,rna,ctg): f.write(H)
px=open(f"{pd}/{m}_COG_features.tsv","w"); px.write("name\ttype\tcontig\tstart\tstop\tstrand\tlegend\n")
for i,r in enumerate(recs):
    seq.append(str(r.seq)); L=len(r.seq)
    ctg.write(f"{m}\tcontig\tother\t{off+1}\t{off+L}\t.\t{'+' if i%2==0 else '-'}\t.\n")
    for f in r.features:
        s,e=int(f.location.start)+1,int(f.location.end); st="+" if f.location.strand==1 else "-"
        if f.type=="CDS":
            lt=f.qualifiers.get("locus_tag",[""])[0]; c=cog.get(lt)
            cds.write(f"{m}\tbakta\t{c if c else 'CDS'}\t{off+s}\t{off+e}\t.\t{st}\t.\n")
            px.write(f"{lt}\tCDS\t{r.id}\t{s}\t{e}\t{1 if st=='+' else -1}\t{('COG '+c) if c else 'No COG'}\n")
        elif f.type in("rRNA","tRNA"):
            rna.write(f"{m}\tbakta\t{f.type}\t{off+s}\t{off+e}\t.\t{st}\t.\n")
    off+=L
fa.write(f">{m}\n"+"".join("".join(seq)[i:i+80]+"\n" for i in range(0,off,80)))
print(m,"contigs",len(recs),"bp",off,"COG-assigned",sum(1 for k in cog))
