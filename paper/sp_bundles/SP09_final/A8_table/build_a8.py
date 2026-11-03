# Final A8 protologue table: 9 candidates + 2 controls, all evidence in one row each.
import csv,os,re,collections
W="${WS2}"; S=W+"/99_thesis_audit/SP09_final"
tx={r["mag_id"]:r for r in csv.DictReader(open(W+"/07_mags/MAG_TAXONOMY.tsv"),delimiter="\t")}
br={r["mag_id"]:r for r in csv.DictReader(open(W+"/26_crossmap/cross_abundance.tsv"),delimiter="\t")}
# GTDB-Tk RED: SP09 rerun for the 8 HQ, workspace_2 run (same version/release) for the others
red={}
for f in [S+"/A1_gtdbtk/out/gtdbtk.bac120.summary.tsv",S+"/A1_gtdbtk/out/gtdbtk.ar53.summary.tsv"]:
    for r in csv.DictReader(open(f),delimiter="\t"): red[r["user_genome"]]=r["red_value"]
wn={m:tx[m]["working_name"] for m in tx}
for f in [W+"/08_taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv",W+"/08_taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"]:
    for r in csv.DictReader(open(f),delimiter="\t"):
        for m,w in wn.items():
            if w==r["user_genome"] and m not in red: red[m]=r["red_value"]
def fasta(p):
    s=[];c=[]
    for l in open(p):
        if l[0]==">":
            if c: s.append("".join(c))
            c=[]
        else: c.append(l.strip().upper())
    s.append("".join(c)); return s
circ={}
for smp,f in(("SED","sediment"),("LIQ","liquid_cultures")):
    for l in open(f"{W}/03_assembly/{f}_assembly_info.txt"):
        if l[0]!="#": p=l.split("\t"); circ[f"{smp}_{p[0]}"]=p[3]
# ANI (max of skani/FastANI) + AF
def ani(m):
    d=f"{S}/A3_ani" if os.path.exists(f"{S}/A3_ani/{m}_skani.tsv") else f"{S}/A13_extra_candidates/ani"
    best=(0,"","")
    for r in csv.DictReader(open(f"{d}/{m}_skani.tsv"),delimiter="\t"):
        if float(r["ANI"])>best[0]: best=(float(r["ANI"]),os.path.basename(r["Ref_file"]).rsplit(".",1)[0],r["Align_fraction_query"])
    for l in open(f"{d}/{m}_fastani.tsv"):
        p=l.split()
        if float(p[2])>best[0]: best=(float(p[2]),os.path.basename(p[1]).rsplit(".",1)[0],f"{int(p[3])/int(p[4])*100:.1f}")
    return best
sp={}
for f in os.listdir(S+"/A2_refs"):
    if f.endswith("_refs.tsv") and f.startswith("ETO_"):
        for r in csv.DictReader(open(f"{S}/A2_refs/{f}"),delimiter="\t"): sp[r["accession"]]=r["species"].replace("s__","")
for r in csv.DictReader(open(S+"/A2_refs/controls_refs.tsv",encoding="utf-8"),delimiter="\t",fieldnames=["a","t"]): sp[r["a"]]=r["t"].split("s__")[-1]
# AAI
aai={}
for f in [S+"/A4_aai/all_pairs_aai.tsv"]:
    for r in csv.DictReader(open(f),delimiter="\t"): aai[(r["genome1"],r["genome2"])]=aai[(r["genome2"],r["genome1"])]=float(r["AAI"])
for f in os.listdir(S+"/A13_extra_candidates/aai/pairs"):
    a,b=f[:-4].split("__"); l=open(S+"/A13_extra_candidates/aai/pairs/"+f).read().splitlines()[-1].split("\t"); aai[(a,b)]=aai[(b,a)]=float(l[4])
def aai_row(m):
    rows=list(csv.DictReader(open(f"{S}/A2_refs/{m}_refs.tsv"),delimiter="\t"))
    by=collections.defaultdict(list); [by[r["gtdb_taxonomy_genus"]].append(r["accession"]) for r in rows]
    out=[]
    for g,mem in by.items():
        mx=max(aai.get((m,a),0) for a in mem); intra=[aai[(x,y)] for i,x in enumerate(mem) for y in mem[i+1:] if (x,y) in aai]
        out.append((mx,g.replace("g__",""),min(intra) if intra else None))
    out.sort(reverse=True); return out
pocp={}
for r in csv.DictReader(open(S+"/A5_pocp/pocp.tsv"),delimiter="\t"):
    pocp.setdefault(r["candidate"],[]).append((float(r["POCPu"]),r["genus"].replace("g__","")))
for l in open(S+"/A13_extra_candidates/pocp/pocp.tsv").read().splitlines()[1:]:
    q,s,v=l.split("\t"); pocp.setdefault(q,[]).append((float(v),""))
cog={}
for m in os.listdir(S+"/A7_maps/cog"):
    p=f"{S}/A7_maps/cog/{m}/cog_classify.tsv"
    if os.path.exists(p): cog[m]=sum(1 for _ in open(p))-1
trna={}
for m in os.listdir(S+"/A7_maps/bakta"):
    p=f"{S}/A7_maps/bakta/{m}/{m}.tsv"
    if os.path.exists(p):
        aas=set(); cds=0
        for l in open(p):
            q=l.split("\t")
            if len(q)>7 and q[1]=="tRNA": aas.add(q[7].split("-")[1] if "-" in q[7] else q[7])
            if len(q)>1 and q[1]=="cds": cds+=1
        trna[m]=(len(aas),cds)
ddh={}
for r in csv.DictReader(open(S+"/dDDH_results/tygs_dDDH_all_MAG_pairs.tsv"),delimiter="\t"):
    if r["mag"] not in ddh or float(r["d4"])>ddh[r["mag"]][0]: ddh[r["mag"]]=(float(r["d4"]),r["ci_d4"],r["subject"])
verdict={
"ETO_SED_MAG_001":"new species (Sedimentibacter)",
"ETO_LIQ_MAG_001":"new species; novel genus-level lineage (RED + tree), but AAI/POCP genus criteria not met",
"ETO_LIQ_MAG_003":"new species (MTP4)",
"ETO_LIQ_MAG_007":"new species (Methanofollis)",
"ETO_LIQ_MAG_009":"new species (Chakrabartyella)",
"ETO_LIQ_MAG_010":"new species (Methanosarcina)",
"ETO_LIQ_MAG_005":"new species (Methanofollis); MIMAG medium",
"ETO_LIQ_MAG_011":"new species (Methanofollis); MIMAG medium",
"ETO_LIQ_MAG_008":"new species (DIGJ01); MIMAG medium",
"ETO_LIQ_MAG_002":"control: Gaoshiqia sediminis (confirmed)",
"ETO_LIQ_MAG_004":"control: Lysinibacillus sphaericus_D (confirmed)"}
tree={"ETO_SED_MAG_001":"sister to S. sp032006665 (97/100)","ETO_LIQ_MAG_001":"sister to CALKIU01 (100/100), next to Robertmurraya (100/100)",
"ETO_LIQ_MAG_003":"basal to MTP4 (65.7/73); MTP4+MAG 100/100","ETO_LIQ_MAG_007":"own branch in Methanofollis (100/100 clade)",
"ETO_LIQ_MAG_009":"within Chakrabartyella (88.6/85)","ETO_LIQ_MAG_010":"within Methanosarcina near sp052057765 (62.6/93)",
"ETO_LIQ_MAG_005":"sister to 011 (100/100), near M. formosanus","ETO_LIQ_MAG_011":"sister to 005 (100/100), near M. formosanus",
"ETO_LIQ_MAG_008":"within DIGJ01","ETO_LIQ_MAG_002":"sister to G. sediminis A06ᵀ (88.3/90)","ETO_LIQ_MAG_004":"within L. sphaericus complex (96.3/100)"}
order=["ETO_SED_MAG_001","ETO_LIQ_MAG_001","ETO_LIQ_MAG_003","ETO_LIQ_MAG_007","ETO_LIQ_MAG_009","ETO_LIQ_MAG_010","ETO_LIQ_MAG_005","ETO_LIQ_MAG_011","ETO_LIQ_MAG_008","ETO_LIQ_MAG_002","ETO_LIQ_MAG_004"]
cols=["mag_id","set","GTDB_classification","size_Mb","contigs","N50","GC_pct","circular","CheckM2_completeness","CheckM2_contamination","GUNC_CSS","MIMAG","5S","16S","23S","tRNA_types","CDS","COG_assigned_pct","RED","nearest_ANI","nearest_reference","AF_query_pct","max_AAI_own_genus","intra_genus_AAI_min","POCPu_vs_type_species","tree_placement","dDDH_max_d4","dDDH_CI","dDDH_nearest","SED_breadth","SED_relabund","SeqCode_type_genome","supported_rank"]
rows=[]
for m in order:
    t=tx[m]; s=fasta(W+"/07_mags/bins/"+t["working_name"]+".fa"); L=sorted(map(len,s),reverse=True); T=sum(L); a=0
    for n50 in L:
        a+=n50
        if a>=T/2: break
    gc=sum(x.count("G")+x.count("C") for x in s)/sum(len(x)-x.count("N") for x in s)*100
    hdr=[l[1:].split()[0] for l in open(W+"/07_mags/bins/"+t["working_name"]+".fa") if l[0]==">"]
    ci="yes" if len(s)==1 and circ.get(hdr[0])=="Y" else "no"
    an=ani(m); ar=aai_row(m) if os.path.exists(f"{S}/A2_refs/{m}_refs.tsv") and m not in("ETO_LIQ_MAG_002","ETO_LIQ_MAG_004") else []
    own=t["genus"]
    if m=="ETO_LIQ_MAG_001":
        mx=f"{ar[0][0]:.1f} (max over 25 genera: {ar[0][1]})"; im="see A4 (per genus)"
    elif ar:
        mx=f"{ar[0][0]:.1f}"; im=f"{ar[0][2]:.1f}" if ar[0][2] else "NA"
    else: mx=im="n/a"
    pc=pocp.get(m,[])
    pcs=(f"{min(pc)[0]:.1f}–{max(pc)[0]:.1f} (25 genera)" if m=="ETO_LIQ_MAG_001" else f"{pc[0][0]:.1f}") if pc else "n/a"
    comp,cont=float(t["completeness"]),float(t["contamination"]); sq=[]
    if comp<=90: sq.append("completeness ≤90")
    if cont>=5: sq.append("contamination ≥5")
    if t["has_16S"]!="yes": sq.append("no 16S")
    if len(s)>100: sq.append(">100 contigs")
    tr=trna.get(m,("",""))
    rows.append([m,"control" if m in("ETO_LIQ_MAG_002","ETO_LIQ_MAG_004") else ("HQ candidate" if t["mimag"]=="high" else "MQ SeqCode-eligible"),
        t["gtdb_classification"].split(";",2)[-1].replace("o__","").replace("f__","").replace("g__","").replace("s__",""),
        f"{T/1e6:.2f}",len(s),n50,f"{gc:.1f}",ci,t["completeness"],t["contamination"],t["gunc_css"],t["mimag"],t["has_5S"],t["has_16S"],t["has_23S"],
        tr[0] or t["n_trna_types"],tr[1],f"{cog.get(m,0)/tr[1]*100:.0f}" if tr[1] and m in cog else "",
        red.get(m,"N/A"),f"{an[0]:.2f}" if an[0]>=80 else "<80",f"{sp.get(an[1],'')} ({an[1]})" if an[1] else "none ≥80%",an[2],mx,im,pcs,tree[m],f"{ddh[m][0]:.1f}" if m in ddh else "n.d.",ddh[m][1] if m in ddh else "",ddh[m][2] if m in ddh else "",
        br[m]["sed_breadth"],br[m].get("sed_relabund",""),"yes" if not sq else "no: "+"; ".join(sq),verdict[m]])
with open(S+"/A8_table/novel_taxa_summary.tsv","w") as f:
    f.write("\t".join(cols)+"\n"); [f.write("\t".join(map(str,r))+"\n") for r in rows]
with open(S+"/A8_table/novel_taxa_summary.md","w") as f:
    f.write("| "+" | ".join(cols)+" |\n|"+"---|"*len(cols)+"\n"); [f.write("| "+" | ".join(str(x).replace("|","/") for x in r)+" |\n") for r in rows]
for r in rows: print(r[0],r[3],r[4],r[7],r[18],r[19],r[22],r[24],r[29])
