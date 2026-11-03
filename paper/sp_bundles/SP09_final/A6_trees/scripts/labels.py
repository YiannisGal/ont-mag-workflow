# build per-set label tables: tip id -> display label, type flag, Etoliko flag
import gzip,csv,sys,os
csv.field_size_limit(sys.maxsize)
S=sys.argv[1]; G=f"{S}/A2_refs/gtdb_metadata"
need=set()
for m in os.listdir(f"{S}/A6_trees/v2"):
    if m.startswith("ETO_") and os.path.isdir(f"{S}/A6_trees/v2/{m}"): need|=set(open(f"{S}/A6_trees/v2/{m}/genomes.txt").read().split())
info={}
for d in("bac120","ar53"):
    for r in csv.DictReader(gzip.open(f"{G}/{d}_metadata_r232.tsv.gz","rt"),delimiter="\t"):
        a=r["accession"][3:]
        if a in need:
            sp=r["gtdb_taxonomy"].split("s__")[-1]; t="type strain" in r["gtdb_type_designation_ncbi_taxa"]
            info[a]=(sp,t,r["ncbi_strain_identifiers"].split(",")[0].strip() if r["ncbi_strain_identifiers"] not in("","none") else "")
for m in os.listdir(f"{S}/A6_trees/v2"):
    p=f"{S}/A6_trees/v2/{m}"
    if not (m.startswith("ETO_") and os.path.isdir(p)): continue
    og=set(open(f"{p}/outgroup.txt").read().split())
    with open(f"{p}/labels.tsv","w") as f:
        f.write("id\tlabel\ttype\tetoliko\toutgroup\n"); f.write(f"{m}\t{m}\tFALSE\tTRUE\tFALSE\n")
        for a in open(f"{p}/genomes.txt").read().split():
            sp,t,st=info.get(a,("?",False,""))
            lab=f"{sp}{' '+st if st else ''}{'ᵀ' if t else ''} ({a})"
            f.write(f"ref_{a}\t{lab}\t{str(t).upper()}\tFALSE\t{str(a in og).upper()}\n")
print("ok")
