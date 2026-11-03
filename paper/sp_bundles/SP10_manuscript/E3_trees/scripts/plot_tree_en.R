# SP10 E3: re-draw of SP09 A6/A13 plot_tree.R in the manuscript style.
# Same treefile + labels.tsv; changes: no title/subtitle, sans-serif >= 7 pt, 17 cm wide,
# Etoliko tips bold in the domain colour (Bacteria #1F5FA6, Archaea #D9622B).
suppressMessages({library(ggtree);library(ape);library(ggplot2)})
a<-commandArgs(TRUE); m<-a[1]; d<-a[2]; out<-a[3]
tr<-read.tree(file.path(d,paste0(m,".treefile")))
lab<-read.delim(file.path(d,"labels.tsv"),stringsAsFactors=FALSE); names(lab)[names(lab)=="label"]<-"disp"
arch<-length(Sys.glob(file.path(d,"*.ar53.msa.faa")))>0
ecol<-if(arch) "#D9622B" else "#1F5FA6"
nl<-tr$node.label; if(is.null(nl)) nl<-rep("",tr$Nnode); nl[is.na(nl)]<-""
sp<-do.call(rbind,lapply(strsplit(ifelse(nl=="","NA/NA",nl),"/"),function(x) as.numeric(x[1:2])))
tr$node.label<-ifelse(!is.na(sp[,1]) & sp[,1]>=80 & sp[,2]>=95, paste0(sp[,1],"/",sp[,2]),"")
n<-Ntip(tr); h_cm<-max(8, n*0.36+1.5)
p<-ggtree(tr,linewidth=0.3) %<+% lab +
  geom_tiplab(aes(label=disp,color=etoliko,fontface=ifelse(etoliko,"bold","plain")),size=2.5,family="sans") +
  geom_nodelab(size=2.2,hjust=1.05,vjust=-0.5,color="grey30",family="sans") +
  scale_color_manual(values=c("TRUE"=ecol,"FALSE"="black"),guide="none") +
  geom_treescale(width=0.05,fontsize=2.5,family="sans",y=-1.2,offset=0.4) + scale_y_continuous(expand=expansion(add=c(1.8,0.6))) + xlim(0,max(node.depth.edgelength(tr))*2.1) +
  theme(text=element_text(family="sans"),plot.background=element_rect(fill="white",colour=NA))
for(ext in c("pdf","png")) ggsave(file.path(out,paste0(m,"_tree.",ext)),p,width=17,height=h_cm,units="cm",limitsize=FALSE,
  dpi=600,bg="white",device=if(ext=="pdf") cairo_pdf else "png")
cat("plotted",m,n,"tips",round(h_cm,1),"cm\n")
