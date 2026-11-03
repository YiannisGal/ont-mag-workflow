# Data for the Figure-5-style panels: per candidate, references with ANI, AAI, dGC, dDDH (pending).
import csv,os,gzip,sys
csv.field_size_limit(sys.maxsize)
S="${WS2}/99_thesis_audit/SP09_final"; P=S+"/A11_species_panels"
def gc(p):
    g=c=n=0
    for l in open(p):
        if l[0]!=">": u=l.strip().upper(); g+=u.count("G"); c+=u.count("C"); n+=len(u)-u.count("N")
    return (g+c)/n*100
def gpath(x):
    for p in(f"{S}/A2_refs/genomes/{x}.fna",f"{S}/A1_gtdbtk/genomes/{x}.fa",f"{S}/A13_extra_candidates/genomes/{x}.fa"):
        if os.path.exists(p): return p
aai={}
for r in csv.DictReader(open(S+"/A4_aai/all_pairs_aai.tsv"),delimiter="\t"): aai[(r["genome1"],r["genome2"])]=aai[(r["genome2"],r["genome1"])]=float(r["AAI"])
for f in os.listdir(S+"/A13_extra_candidates/aai/pairs"):
    a,b=f[:-4].split("__"); aai[(a,b)]=aai[(b,a)]=float(open(S+"/A13_extra_candidates/aai/pairs/"+f).read().splitlines()[-1].split("\t")[4])
meta={}
for d in("bac120","ar53"):
    for r in csv.DictReader(gzip.open(f"{S}/A2_refs/gtdb_metadata/{d}_metadata_r232.tsv.gz","rt"),delimiter="\t"):
        meta[r["accession"][3:]]=(r["gtdb_taxonomy"].split("s__")[-1],"type strain" in r["gtdb_type_designation_ncbi_taxa"],r["ncbi_strain_identifiers"].split(",")[0].strip())
def anis(m):
    d=f"{S}/A3_ani" if os.path.exists(f"{S}/A3_ani/{m}_skani.tsv") else f"{S}/A13_extra_candidates/ani"
    a={}
    for r in csv.DictReader(open(f"{d}/{m}_skani.tsv"),delimiter="\t"):
        k=os.path.basename(r["Ref_file"]).rsplit(".",1)[0]; a[k]=max(a.get(k,0),float(r["ANI"]))
    for l in open(f"{d}/{m}_fastani.tsv"):
        p=l.split(); k=os.path.basename(p[1]).rsplit(".",1)[0]; a[k]=max(a.get(k,0),float(p[2]))
    return a
pm={}
for r in csv.DictReader(open(P+"/panel_members.tsv"),delimiter="\t"): pm.setdefault(r["candidate"],[]).append(r["accession"])
MF=["ETO_LIQ_MAG_005","ETO_LIQ_MAG_007","ETO_LIQ_MAG_011"]
tree={m:f"{S}/A6_trees/v2/{m}/{m}.treefile" for m in pm}
for m in MF: tree[m]=f"{S}/A13_extra_candidates/trees/methanofollis/methanofollis.treefile"
tree["ETO_LIQ_MAG_008"]=f"{S}/A13_extra_candidates/trees/digj01/digj01.treefile"
for m in["ETO_LIQ_MAG_005","ETO_LIQ_MAG_011","ETO_LIQ_MAG_008"]:
    pm[m]=[l.split("\t")[0] for l in open(f"{S}/A2_refs/{m}_refs.tsv").read().splitlines()[1:]]
for m in MF: pm[m]=[x for x in pm.get(m,[])]+[o for o in MF if o!=m]
# Methanofollis panel 007 also uses the shared tree: references = its A2 set
pm["ETO_LIQ_MAG_007"]=[l.split("\t")[0] for l in open(f"{S}/A2_refs/ETO_LIQ_MAG_007_refs.tsv").read().splitlines()[1:]]+["ETO_LIQ_MAG_005","ETO_LIQ_MAG_011"]
for m,refs in pm.items():
    A=anis(m); g0=gc(gpath(m)); tips=[]
    # Etoliko-vs-Etoliko ANI from skani triangle
    tri={}
    for r in csv.DictReader(open(S+"/A12_MQ_ani_screen/etoliko_all_vs_all.tsv"),delimiter="\t"):
        a=os.path.basename(r["Ref_file"]).rsplit(".",1)[0]; b=os.path.basename(r["Query_file"]).rsplit(".",1)[0]; tri[(a,b)]=tri[(b,a)]=float(r["ANI"])
    with open(f"{P}/data/{m}.tsv","w") as f:
        f.write("tip\tlabel\ttype\tetoliko\tANI\tAAI\tdDDH\tdGC\n")
        f.write(f"{m}\t{m}\tFALSE\tTRUE\t100\t100\tNA\t0\n")
        for a in refs:
            if a.startswith("ETO_"):
                f.write(f"{a}\t{a} (Etoliko MAG)\tFALSE\tTRUE\t{tri.get((m,a),'NA')}\t{aai.get((m,a),'NA')}\tNA\t{gc(gpath(a))-g0:.2f}\n"); continue
            sp,t,st=meta.get(a,("?",False,"")); lab=f"{sp}{' '+st if st else ''}{'ᵀ' if t else ''}"
            f.write(f"ref_{a}\t{lab}\t{str(t).upper()}\tFALSE\t{A.get(a,'NA')}\t{aai.get((m,a),'NA')}\tNA\t{gc(gpath(a))-g0:.2f}\n")
    open(f"{P}/data/{m}.tree","w").write(tree[m])
    print(m,len(refs),"refs; tree:",os.path.basename(tree[m]))
