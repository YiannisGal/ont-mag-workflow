# S_mcrA_tree and S_dsr_trees from the SP05 IQ-TREE runs (no recomputation).
# mcrA: restriction-enzyme "McrA" pruned; rooted on the ACR alpha clade when monophyletic, else drawn unrooted as written.
# dsr: one panel per subunit, drawn as written by IQ-TREE (no root is implied), UFBoot shown >= 70.
suppressMessages({library(ggtree);library(ape);library(ggplot2);library(patchwork)})
a<-commandArgs(TRUE); TD<-a[1]; LD<-a[2]; OUT<-a[3]
pal<-c("query"="#D9622B","MCR alpha"="#7a3b1b","ACR alpha"="#1F5FA6","ACR other subunits"="#8fb3dd",
       "ACR component A2"="#5a5a5a","reductive"="#1b7f5f","oxidative (rDsr)"="#b8860b")
ufb<-function(tr){ nl<-tr$node.label; nl[is.na(nl)]<-""; v<-suppressWarnings(as.numeric(sub(".*/","",nl))); ifelse(!is.na(v)&v>=70,as.character(round(v)),"") }
draw<-function(tr,lab,legend_title="",show_legend=TRUE){
  tr$node.label<-ufb(tr); names(lab)[names(lab)=="label"]<-"disp"
  p<-ggtree(tr,linewidth=0.25) %<+% lab +
    geom_tiplab(aes(label=disp,colour=cat,fontface=ifelse(cat=="query","bold","plain")),size=2.5,family="sans",show.legend=FALSE)+
    geom_tippoint(aes(colour=cat),size=1.1)+
    geom_nodelab(size=2.2,hjust=1.15,vjust=-0.35,colour="grey35",family="sans")+
    scale_colour_manual(values=pal,name=legend_title,breaks=intersect(names(pal),unique(lab$cat)),
      labels=sub("^query$","this study",intersect(names(pal),unique(lab$cat))))+
    geom_treescale(width=0.2,fontsize=2.5,family="sans",y=-2,offset=0.6)+
    scale_y_continuous(expand=expansion(add=c(3,1)))+
    xlim(0,max(node.depth.edgelength(tr))*1.7)+
    theme(text=element_text(family="sans",size=8),legend.position=if(show_legend) "bottom" else "none",
          legend.text=element_text(size=7),legend.title=element_blank(),plot.background=element_rect(fill="white",colour=NA))+
    guides(colour=guide_legend(override.aes=list(size=2.5),nrow=2))
  p
}
# ---- mcrA
tr<-read.tree(file.path(TD,"mcrA.treefile")); lab<-read.delim(file.path(LD,"mcrA_tips.tsv"),stringsAsFactors=FALSE)
tr<-drop.tip(tr,lab$tip[lab$cat=="DROP"]); lab<-lab[lab$cat!="DROP",]
acr<-lab$tip[lab$cat=="ACR alpha"]
if(length(acr)>1 && is.monophyletic(tr,acr,reroot=TRUE)){ tr<-root(tr,outgroup=acr,resolve.root=TRUE); rooted<-"ACR alpha" } else { rooted<-"as written by IQ-TREE (unrooted)" }
p<-draw(tr,lab)
h<-max(10,Ntip(tr)*0.3+3)
for(ext in c("pdf","png")) ggsave(file.path(OUT,paste0("S_mcrA_tree.",ext)),p,width=17,height=h,units="cm",limitsize=FALSE,dpi=600,bg="white",device=if(ext=="pdf") cairo_pdf else "png")
cat("mcrA:",Ntip(tr),"tips, rooted:",rooted,"\n")
# ---- dsrA / dsrB
ps<-list()
for(g in c("dsrA","dsrB")){
  tr<-read.tree(file.path(TD,paste0(g,".contree"))); lab<-read.delim(file.path(LD,paste0(g,"_tips.tsv")),stringsAsFactors=FALSE)
  ps[[g]]<-draw(tr,lab,show_legend=(g=="dsrB"))+labs(tag=if(g=="dsrA") "a" else "b")+theme(plot.tag=element_text(face="bold",size=9))
  cat(g,Ntip(tr),"tips\n")
}
out<-ps$dsrA/ps$dsrB
for(ext in c("pdf","png")) ggsave(file.path(OUT,paste0("S_dsr_trees.",ext)),out,width=17,height=34,units="cm",dpi=600,bg="white",device=if(ext=="pdf") cairo_pdf else "png")
