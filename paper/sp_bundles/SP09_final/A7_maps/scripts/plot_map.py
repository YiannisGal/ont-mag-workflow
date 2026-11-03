# Circular genome map (CGView/Proksee style) from a Bakta GBFF with pyCirclize.
import sys,re,collections
import numpy as np
from Bio import SeqIO
from pycirclize import Circos
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
m,gbff,outp,domain,circ=sys.argv[1:6]
EGG={}
if len(sys.argv)>6:
    for l in open(sys.argv[6]):
        if l.startswith('#'): continue
        p=l.rstrip('\n').split('\t')
        if len(p)==2 and p[1]: EGG[p[0]]=p[1][0]
        elif len(p)>6 and p[6] not in ('','-'): EGG[p[0]]=p[6][0]
recs=sorted(SeqIO.parse(gbff,"genbank"),key=lambda r:-len(r.seq))
COG={"J":"#ff9e9e","A":"#ffb6c1","K":"#ff6f6f","L":"#c03d3d","B":"#f4c2c2","D":"#ffd27f","Y":"#f5deb3","V":"#ffb347","T":"#ffa500","M":"#c9a227","N":"#e6c35c","Z":"#d2b48c","W":"#bfa98a","U":"#a0d468","O":"#4caf50","X":"#7cb342","C":"#7ec8e3","G":"#4fa3d1","E":"#1f78b4","F":"#6a8fd8","H":"#9c88d9","I":"#b39ddb","P":"#7b68ee","Q":"#5e35b1","R":"#bdbdbd","S":"#9e9e9e"}
COGN={"A":"RNA processing","B":"Chromatin","W":"Extracellular structures","Y":"Nuclear structure","Z":"Cytoskeleton","J":"Translation","K":"Transcription","L":"Replication/repair","D":"Cell cycle","V":"Defense","T":"Signal transduction","M":"Cell envelope","N":"Motility","U":"Secretion","O":"PTM/chaperones","X":"Mobilome","C":"Energy","G":"Carbohydrate","E":"Amino acid","F":"Nucleotide","H":"Coenzyme","I":"Lipid","P":"Inorganic ion","Q":"Secondary metab.","R":"General function","S":"Unknown function"}
NOCOG="#e0e0e0"
CNAME={}
import os as _os
_cf=_os.path.join(_os.path.dirname(sys.argv[6]) if len(sys.argv)>6 else "",m,"cog_classify.tsv")
if _os.path.exists(_cf):
    for _l in open(_cf):
        _p=_l.rstrip("\n").split("\t")
        if len(_p)>5 and _p[5]:
            gn=_p[5]
            if gn.startswith('Mcr'): continue  # McrA/Mcr1 = restriction enzyme / flavin reductase, not methyl-CoM reductase (Mcor*)
            gn={"McorA":"mcrA","McorB":"mcrB","McorG":"mcrG"}.get(gn,gn)
            CNAME[_p[0]]=gn[0].lower()+gn[1:]
def cog(f):
    lt=f.qualifiers.get('locus_tag',[''])[0]
    if lt in EGG: return EGG[lt]
    for x in f.qualifiers.get("db_xref",[]):
        mm=re.match(r"COG:([A-Z])$",x)
        if mm: return mm.group(1)
    return None
ARCH={"ETO_LIQ_MAG_003","ETO_LIQ_MAG_005","ETO_LIQ_MAG_007","ETO_LIQ_MAG_010","ETO_LIQ_MAG_011"}
key=re.compile(r"^(mcr[ABG]|mtr[A-H]|mta[ABC]|mtb[ABC]|mtt[BC])$") if m in ARCH else re.compile(r"^(spo0A|spoIIIE|spoVT|sigF|sigG|sigK|katA|sodA|cspA|sspA|ssp[A-Z])$")
stats=collections.Counter(); cogc=collections.Counter(); keyf=[]
circos=Circos({r.id:len(r.seq) for r in recs},space=0 if len(recs)==1 else (0.3 if len(recs)>30 else 1))
for i,r in enumerate(recs):
    sec=circos.get_sector(r.id)
    # ring 1: contigs (alternating shade)
    sec.add_track((94,99)).axis(fc="#5b6c8f" if i%2==0 else "#a6b4cf",ec="none")
    fw=sec.add_track((86,93)); rv=sec.add_track((79,86))
    rr=sec.add_track((73,78)); gc=sec.add_track((58,72)); sk=sec.add_track((42,57))
    for f in r.features:
        s,e=int(f.location.start),int(f.location.end)
        if f.type=="CDS":
            c=cog(f); col=COG.get(c,NOCOG) if c else NOCOG; stats["CDS"]+=1; cogc[c or "-"]+=1
            (fw if f.location.strand==1 else rv).rect(s,e,fc=col,ec="none")
            lt0=f.qualifiers.get("locus_tag",[""])[0]
            g=CNAME.get(lt0,"") if m in ARCH else (f.qualifiers.get("gene",[""])[0] or CNAME.get(lt0,""))
            if key.match(g): keyf.append((r.id,(s+e)/2,g))
        elif f.type=="rRNA":
            stats["rRNA"]+=1; rr.rect(s,e,fc="#d62728",ec="none")
            p=f.qualifiers.get("product",[""])[0]
            if "16S" in p: keyf.append((r.id,(s+e)/2,"16S rRNA"))
        elif f.type=="tRNA": stats["tRNA"]+=1; rr.rect(s,e,fc="#2ca02c",ec="none")
    seq=str(r.seq).upper(); L=len(seq); win=max(500,min(10000,L//200 or 500)); step=win//2 or 1
    pos=np.arange(0,max(L-win,1),step)
    g=np.array([(seq[p:p+win].count("G")+seq[p:p+win].count("C"))/max(1,win) for p in pos]) if len(pos) else np.array([])
    k=np.array([(seq[p:p+win].count("G")-seq[p:p+win].count("C"))/max(1,seq[p:p+win].count("G")+seq[p:p+win].count("C")) for p in pos]) if len(pos) else np.array([])
    if len(pos)>1:
        x=pos+win/2; gm=(sum(seq.count(b) for b in "GC")/L)
        gc.fill_between(x,np.clip(g-gm,-0.15,0.15),np.zeros_like(g),vmin=-0.15,vmax=0.15,fc="#333333",ec="none")
        k=np.clip(k,-0.3,0.3); kp=np.where(k>0,k,0); kn=np.where(k<0,k,0)
        sk.fill_between(x,kp,np.zeros_like(kp),vmin=-0.3,vmax=0.3,fc="#2a9d8f",ec="none")
        sk.fill_between(x,kn,np.zeros_like(kn),vmin=-0.3,vmax=0.3,fc="#e76f51",ec="none")
    if i==0:
        tr=sec.add_track((99,100))
        if L>200000: tr.xticks_by_interval(500000 if L>2e6 else 100000,label_formatter=lambda v:f"{v/1e6:.1f} Mb",label_size=6,outer=True)
# label key genes / 16S (dedupe nearby)
lab=collections.defaultdict(list)
for sid,x,gname in keyf: lab[sid].append((x,gname))
for sid,v in lab.items():
    sec=circos.get_sector(sid); t=sec.add_track((99,100))
    v.sort(); groups=[]
    for x,gname in v:
        if groups and x-groups[-1][-1][0]<20000 and ("rRNA" in gname)==("rRNA" in groups[-1][-1][1]): groups[-1].append((x,gname))
        else: groups.append([(x,gname)])
    last=-1e9
    for gr in groups:
        xs=[a for a,_ in gr]; x=sum(xs)/len(xs)
        if "rRNA" in gr[0][1]: lab="16S rRNA" if len(gr)==1 else f"16S rRNA ×{len(gr)}"
        else:
            pre=collections.OrderedDict()
            for _,gn in gr: pre.setdefault(gn[:3],[]); (gn[3:] not in pre[gn[:3]]) and pre[gn[:3]].append(gn[3:])
            lab="–".join(k+"".join(s) for k,s in pre.items())
        if x-last<sec.size*0.012: continue
        t.annotate(x,lab,label_size=5.5,text_kws=dict(style="normal" if "rRNA" in lab else "italic"),line_kws=dict(lw=0.4))
        last=x
tot=sum(len(r.seq) for r in recs); gcp=sum(str(r.seq).upper().count(b) for r in recs for b in "GC")/tot*100
circos.text(f"{m}\n{tot/1e6:.2f} Mb · {len(recs)} contig{'s' if len(recs)>1 else ''}{' (circular)' if circ=='Y' else ''}\nGC {gcp:.1f}% · {stats['CDS']} CDS · {stats['rRNA']} rRNA · {stats['tRNA']} tRNA",size=8)
fig=circos.plotfig(figsize=(9,9),dpi=600)
used=[c for c in COG if cogc[c]>0]
h=[Patch(fc=COG[c],label=f"{c} {COGN.get(c,'')}") for c in used]+[Patch(fc=NOCOG,label="No COG assigned")]
h2=[Patch(fc="#5b6c8f",label="Contigs (alternating, by length)"),Patch(fc="#d62728",label="rRNA"),Patch(fc="#2ca02c",label="tRNA"),Patch(fc="#333333",label="GC content (dev. from mean)"),Patch(fc="#2a9d8f",label="GC skew +"),Patch(fc="#e76f51",label="GC skew −")]
fig.legend(handles=h2,loc="upper left",fontsize=6,frameon=False,title="Rings (outside→inside: contigs, CDS +, CDS −, RNA, GC, skew)",title_fontsize=6)
fig.legend(handles=h,loc="lower left",fontsize=5.5,frameon=False,ncol=2,title="CDS colour: COG category",title_fontsize=6)
fig.patch.set_facecolor("white")
fig.savefig(f"{outp}.png",dpi=600,facecolor="white"); fig.savefig(f"{outp}.svg",facecolor="white")
print(m,dict(stats),"COG-assigned:",stats["CDS"]-cogc["-"],"key labels:",len(keyf))
