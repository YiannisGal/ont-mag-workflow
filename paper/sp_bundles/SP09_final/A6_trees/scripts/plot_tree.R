suppressMessages({library(ggtree);library(ape);library(ggplot2)})
a<-commandArgs(TRUE); m<-a[1]; d<-a[2]
tr<-read.tree(file.path(d,paste0(m,".treefile")))
lab<-read.delim(file.path(d,"labels.tsv"),stringsAsFactors=FALSE); names(lab)[names(lab)=="label"]<-"disp"
# support: IQ-TREE node labels "SH/UFB"
nl<-tr$node.label; if(is.null(nl)) nl<-rep("",tr$Nnode); nl[is.na(nl)]<-""; sp<-do.call(rbind,lapply(strsplit(ifelse(nl=="","NA/NA",nl),"/"),function(x) as.numeric(x[1:2])))
tr$node.label<-ifelse(!is.na(sp[,1]) & sp[,1]>=80 & sp[,2]>=95, paste0(sp[,1],"/",sp[,2]),"")
n<-Ntip(tr); h<-max(6, n*0.16)
p<-ggtree(tr,linewidth=0.3) %<+% lab +
  geom_tiplab(aes(label=disp,color=etoliko,fontface=ifelse(etoliko,"bold","plain")),size=2.3,align=FALSE) +
  geom_nodelab(size=1.6,hjust=1.15,vjust=-0.4,color="grey30") +
  scale_color_manual(values=c("TRUE"="#C00000","FALSE"="black"),guide="none") +
  geom_treescale(width=0.05,fontsize=2) + xlim(0,max(node.depth.edgelength(tr))*1.9) +
  ggtitle(paste0(m," — GTDB-Tk ",ifelse(length(Sys.glob(file.path(d,"*.ar53.msa.faa"))),"ar53","bac120")," markers, IQ-TREE 3.1.1"),subtitle="Node labels: SH-aLRT/UFBoot, shown where ≥80/≥95. ᵀ = type strain. Red bold = Etoliko MAG.") +
  theme(plot.title=element_text(size=9),plot.subtitle=element_text(size=7),plot.background=element_rect(fill="white",colour=NA))
for(ext in c("pdf","svg")) ggsave(file.path(d,"..","..",paste0(m,"_tree.",ext)),p,width=8,height=h,limitsize=FALSE,device=if(ext=="pdf") cairo_pdf else grDevices::svg)
cat("plotted",m,"\n")
