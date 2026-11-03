# Figure-6-style circular tree: Etoliko MAGs + culturomics isolates + GTDB r232 context
suppressMessages({library(ggtree);library(treeio);library(ape);library(ggplot2);library(ggnewscale);library(dplyr)})
setwd(commandArgs(TRUE)[1])
tr<-read.tree("combined_both_domains.nwk")
an<-read.delim("tip_annotation.tsv",stringsAsFactors=FALSE,colClasses="character")
an$etoliko<-an$etoliko=="TRUE"; names(an)[names(an)=="label"]<-"disp"; rownames(an)<-an$tip
n<-Ntip(tr); eto<-an$tip[an$etoliko]
# maximal reference-only clades -> collapse; label with deepest shared GTDB rank
ranks<-c("phylum","class","order","genus")
desc<-lapply((n+1):(n+tr$Nnode),function(k) extract.clade(tr,k)$tip.label)
names(desc)<-(n+1):(n+tr$Nnode)
refonly<-sapply(desc,function(x) !any(x %in% eto) & length(x)>1 & length(unique(an[x,"phylum"]))==1)
par<-tr$edge[match((n+1):(n+tr$Nnode),tr$edge[,2]),1]; names(par)<-(n+1):(n+tr$Nnode)
maximal<-names(refonly)[refonly & !(as.character(par[names(refonly)]) %in% names(refonly)[refonly])]
clab<-sapply(maximal,function(k){x<-an[desc[[k]],]; for(r in rev(ranks)) if(length(unique(x[[r]]))==1) return(paste0(unique(x[[r]]))); unique(x$domain)[1]})
incl<-unique(unlist(desc[maximal])); an$single<-!an$etoliko & !(an$tip %in% incl)
p<-ggtree(tr,layout="fan",open.angle=8,linewidth=0.25) %<+% an
# domain shading
dom_node<-c(Bacteria=getMRCA(tr,an$tip[an$domain=="Bacteria"]),Archaea=getMRCA(tr,an$tip[an$domain=="Archaea"]))
p<-p+geom_hilight(data=data.frame(node=dom_node,dom=names(dom_node)),aes(node=node,fill=dom),alpha=0.10,extend=1.75)+
  scale_fill_manual(values=c(Bacteria="#c0504d",Archaea="#7f7f7f"),name="Domain")+new_scale_fill()
# Etoliko-containing phyla
ph<-sort(unique(an$phylum[an$etoliko]))
phn<-sapply(ph,function(x){tp<-an$tip[an$phylum==x]; if(length(tp)>1) getMRCA(tr,tp) else which(tr$tip.label==tp)})
phcol<-setNames(c("#8c6d31","#e7ba52","#5b8a72","#843c39","#3182bd","#9c9ede")[seq_along(ph)],ph)
p<-p+geom_hilight(data=data.frame(node=phn,ph=ph),aes(node=node,fill=ph),alpha=0.28,extend=1.1)+scale_fill_manual(values=phcol,name="Phyla with Etoliko genomes")+new_scale_fill()
for(k in maximal) p<-collapse(p,as.numeric(k),mode="max",fill="grey45",colour="grey30",alpha=0.85)
cl<-data.frame(node=as.numeric(maximal),lab=clab)
p<-p+geom_cladelab(data=cl,mapping=aes(node=node,label=lab),fontsize=2.1,barsize=0,offset=0.06,angle="auto",horizontal=FALSE,textcolour="grey20")
p<-p+geom_tiplab(aes(label=ifelse(single,reflab,NA)),colour="grey35",size=1.6,offset=0.05)
p<-p+geom_tiplab(aes(label=ifelse(etoliko,disp,NA)),colour="#C00000",fontface="bold",size=1.9,offset=0.05,align=FALSE)
# rings: source, novelty, sediment rel. abundance
ring<-an[,c("source","novelty")]; ring[ring==""]<-NA; rownames(ring)<-an$tip
p2<-gheatmap(p,ring[,"source",drop=FALSE],offset=1.05,width=0.05,colnames=FALSE,color=NA)+
  scale_fill_manual(values=c("SED MAG"="#7b4f2c","LIQ MAG"="#d9a441","Isolate"="#2c7fb8"),na.value=NA,name="Source",na.translate=FALSE)+new_scale_fill()
p3<-gheatmap(p2,ring[,"novelty",drop=FALSE],offset=1.17,width=0.05,colnames=FALSE,color=NA)+
  scale_fill_manual(values=c("new species (HQ)"="#08306b","new species (MQ, SeqCode-eligible)"="#4292c6","new species-level lineage (ANI screen)"="#c6dbef"),na.value=NA,name="Novelty (SP09)",na.translate=FALSE)+new_scale_fill()
ra<-data.frame(sed=as.numeric(an$sed_relabund),row.names=an$tip)
p4<-gheatmap(p3,ra,offset=1.29,width=0.05,colnames=FALSE,color=NA)+scale_fill_gradient(low="#e0ecf4",high="#0c2c84",na.value=NA,name="Sediment rel.\nabundance (%)",trans="sqrt")
p4<-p4+theme(legend.position="right",legend.text=element_text(size=7),legend.title=element_text(size=8),plot.background=element_rect(fill="white",colour=NA))+
  labs(caption="GTDB-Tk 2.7.2 (r232) bac120/ar53 markers, IQ-TREE 3.1.1 (LG+F+G4). Bacterial and archaeal trees midpoint-rooted and joined at the root for display only.\nGrey triangles: collapsed GTDB reference clades without Etoliko genomes. Red: Etoliko MAGs (ETO_*) and culturomics isolates (YG*).")
for(ext in c("pdf","svg")) ggsave(paste0("Etoliko_combined_tree.",ext),p4,width=14,height=13,device=if(ext=="pdf") cairo_pdf else grDevices::svg)
ggsave("Etoliko_combined_tree.png",p4,width=14,height=13,dpi=400,bg="white")
cat("collapsed",length(maximal),"clades\n")
