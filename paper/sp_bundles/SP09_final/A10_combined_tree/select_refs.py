# Pick GTDB r232 species representatives as context for the combined Etoliko tree.
import gzip,csv,sys,re,collections
csv.field_size_limit(sys.maxsize)
W="${WS2}"; S=W+"/99_thesis_audit/SP09_final"
eto=[r["gtdb_classification"] for r in csv.DictReader(open(W+"/07_mags/MAG_TAXONOMY.tsv"),delimiter="\t")]
eto+=["d__Bacteria;p__Bacillota;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium_J;s__",
      "d__Bacteria;p__Bacillota;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Terrisporobacter;s__",
      "d__Bacteria;p__Bacillota;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Paraclostridium;s__",
      "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales_A;f__Planococcaceae;g__Sporosarcina;s__"]
orders={t.split(";")[3] for t in eto}; genera={t.split(";")[5] for t in eto if t.split(";")[5]!='g__'}; phyla={t.split(";")[1] for t in eto}
named=lambda x:not re.search(r"\d",x.split("__")[1])
best={}
for d in("bac120","ar53"):
    for r in csv.DictReader(gzip.open(f"{S}/A2_refs/gtdb_metadata/{d}_metadata_r232.tsv.gz","rt"),delimiter="\t"):
        if r["gtdb_representative"]!="t": continue
        t=r["gtdb_taxonomy"].split(";"); dom,p,c,o,f=t[:5]
        g=t[5]
        if g in genera and r["gtdb_type_species_of_genus"]=="t": key=("genus",g)
        elif o in orders: key=("family",f)
        elif p in phyla: key=("class",c)
        elif named(p): key=("phylum",p)
        else: continue
        q=("type strain" in r["gtdb_type_designation_ncbi_taxa"], float(r["checkm2_completeness"] or 0)-5*float(r["checkm2_contamination"] or 0))
        if key[0]=="genus": best[key]=((True,999),r["accession"][3:],r["gtdb_taxonomy"],d); continue
        if key not in best or q>best[key][0]: best[key]=(q,r["accession"][3:],r["gtdb_taxonomy"],d)
# cap families per order at 8 (type material, then quality); keep families that hold Etoliko genomes
efam={t.split(";")[4] for t in eto}
byo=collections.defaultdict(list)
for k,v in best.items():
    if k[0]=="family": byo[v[2].split(";")[3]].append((k[1] not in efam,not v[0][0],-v[0][1],k))
for o,l in byo.items():
    for x in sorted(l)[8:]:
        if x[3][1] not in efam: del best[x[3]]
with open(f"{S}/A10_combined_tree/refs.tsv","w") as o:
    o.write("accession\tlevel\tgroup\tgtdb_taxonomy\tdomain\ttype_strain\n")
    for (lv,g),(q,a,t,d) in sorted(best.items()): o.write(f"{a}\t{lv}\t{g}\t{t}\t{d}\t{q[0]}\n")
c=collections.Counter((v[3],k[0]) for k,v in best.items()); print(c, "total",len(best))
