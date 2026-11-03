# SP10 E3: re-draw of SP09 A11 plot_panel.R. Same data; no caption, sans >= 7 pt, 17 cm, Etoliko bold in domain colour.
suppressMessages({library(ggtree);library(ape);library(ggplot2);library(aplot)})
a<-commandArgs(TRUE); m<-a[1]; P<-a[2]; OUT<-a[3]; ecol<-a[4]
d<-read.delim(file.path(P,"data",paste0(m,".tsv")),stringsAsFactors=FALSE,na.strings="NA")
tr<-read.tree(readLines(file.path(P,"data",paste0(m,".tree")),warn=FALSE)[1])
og<-tr$tip.label[!(tr$tip.label %in% d$tip) & grepl("^ref_",tr$tip.label)]
labs<-unlist(lapply(Sys.glob(file.path(dirname(readLines(file.path(P,"data",paste0(m,".tree")),warn=FALSE)[1]),"labels.tsv")),function(f){x<-read.delim(f,stringsAsFactors=FALSE);setNames(x$label,x$id)}))
# keep panel tips + outgroup tips (outgroup = tips listed in outgroup.txt)
ogf<-file.path(dirname(readLines(file.path(P,"data",paste0(m,".tree")),warn=FALSE)[1]),"outgroup.txt")
ogs<-if(file.exists(ogf)) paste0("ref_",readLines(ogf)) else character(0)
keep<-intersect(tr$tip.label,c(d$tip,ogs)); tr<-keep.tip(tr,keep)
for(o in intersect(ogs,tr$tip.label)) d<-rbind(d,data.frame(tip=o,label=sub(" \\(G.*","",ifelse(is.na(labs[o]),o,labs[o])),type=FALSE,etoliko=FALSE,ANI=NA,AAI=NA,dDDH=NA,dGC=NA))
d$label<-sub(" \\((GC[AF]_[0-9.]+)\\)$","",d$label)
p<-ggtree(tr,linewidth=0.35)
td<-p$data[p$data$isTip,]; ord<-td$label[order(-td$y)]; d<-d[match(ord,d$tip),]; d$y<-rev(seq_along(ord))
lab<-ggplot(d,aes(x=0,y=factor(tip,levels=rev(ord))))+geom_text(aes(label=label,colour=etoliko,fontface=ifelse(etoliko,"bold.italic","italic")),hjust=0,size=2.5,family="sans")+
  scale_colour_manual(values=c("TRUE"=ecol,"FALSE"="black"),guide="none")+xlim(0,1)+theme_void()
long<-rbind(data.frame(tip=d$tip,col="ANI",v=d$ANI),data.frame(tip=d$tip,col="AAI",v=d$AAI),data.frame(tip=d$tip,col="dDDH",v=d$dDDH))
og_tip<-d$tip %in% ogs; long$txt<-ifelse(is.na(long$v),ifelse(long$col=="dDDH" | (long$col=="AAI" & rep(og_tip,3)),"n.d.","<80"),formatC(round(long$v,1),format="fg"))
long$tip<-factor(long$tip,levels=rev(ord))
long$col<-factor(long$col,levels=c("ANI","AAI","dDDH"))
rel<-ggplot(long,aes(col,tip))+geom_tile(aes(fill=v),colour="white",linewidth=0.6)+geom_text(aes(label=txt,colour=!is.na(v)&v>60),size=2.5,family="sans")+
  scale_fill_gradientn(colours=c("#def5e5","#5ec9b5","#3b7fa6","#3b3a6b","#0b0405"),limits=c(20,100),na.value="grey88",name="Relatedness (%)")+
  scale_colour_manual(values=c("TRUE"="white","FALSE"="black"),guide="none")+scale_x_discrete(position="top")+
  theme_minimal()+theme(axis.title=element_blank(),axis.text.y=element_blank(),panel.grid=element_blank(),axis.text.x=element_text(face="bold",size=8))
gcd<-data.frame(tip=factor(d$tip,levels=rev(ord)),col="Δ G+C",v=d$dGC)
mx<-max(2,max(abs(gcd$v),na.rm=TRUE))
gcp<-ggplot(gcd,aes(col,tip))+geom_tile(aes(fill=v),colour="white",linewidth=0.6)+geom_text(aes(label=ifelse(is.na(v),"",formatC(round(v,2),format="fg"))),size=2.5,family="sans")+
  scale_fill_gradient2(low="#3b5ba5",mid="white",high="#b2182b",midpoint=0,limits=c(-mx,mx),na.value="grey88",name="Δ G+C (%)")+
  scale_x_discrete(position="top")+theme_minimal()+theme(axis.title=element_blank(),axis.text.y=element_blank(),panel.grid=element_blank(),axis.text.x=element_text(face="bold",size=8))
n<-length(ord); yl<-c(0.5,n+0.5)
pt<-p+scale_y_continuous(limits=yl,expand=c(0,0))+theme(plot.margin=margin(28,0,5,5))
fix<-function(g) g+scale_y_discrete(expand=c(0,0))+coord_cartesian(ylim=yl,clip="off")
library(patchwork)
out<-pt+fix(rel)+fix(gcp)+fix(lab)+plot_layout(widths=c(2.6,1.5,0.5,3.4),guides="collect") & theme(text=element_text(family="sans",size=8),legend.position="bottom",legend.title=element_text(size=7),legend.text=element_text(size=7),plot.background=element_rect(fill="white",colour=NA))
h<-max(9,0.5*length(ord)+3)
for(ext in c("pdf","png")) ggsave(file.path(OUT,paste0(m,"_species_panel.",ext)),out,width=17,height=h,units="cm",dpi=600,bg="white",device=if(ext=="pdf") cairo_pdf else "png")
cat("panel",m,length(ord),"tips\n")
