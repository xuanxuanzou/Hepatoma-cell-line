### Drug_targets heatmap
library(dplyr)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)


snRNA <- readRDS(infile)
Drug_targets <- read.table("./Drug_targets.tsv",header = F)

plot_data1 <- snRNA@assays$RNA@scale.data[as.vector(Drug_targets$V1[Drug_targets$V1%in%rownames(snRNA@assays$RNA@scale.data)]),]
dim(plot_data1)

### cell line and seurat cluster

setwd(outdir)
SKHep_meta<-read.table("SKHep_cell_cluster.txt",sep ="\t",header=T,row.names = 1,check.names=F)
MHCC97H_meta<-read.table("MHCC97H_cell_cluster.txt",sep ="\t",header=T,row.names = 1,check.names=F)
MHCC97L_meta<-read.table("MHCC97L_cell_cluster.txt",sep ="\t",header=T,row.names = 1,check.names=F)
Huh7_meta<-read.table("Huh7_cell_cluster.txt",sep ="\t",header=T,row.names = 1,check.names=F)
HepG2_meta<-read.table("HepG2_cell_cluster.txt",sep ="\t",header=T,row.names = 1,check.names=F)
HepG2_meta$ID<-paste0("HepG2_",rownames(HepG2_meta))
colnames(HepG2_meta)[1]<-"Idents(snRNA)"

meta <- rbind(SKHep_meta,MHCC97H_meta,MHCC97L_meta,Huh7_meta,HepG2_meta)
rownames(meta) <- meta$ID
colnames(meta)[1] <- "Cluster"
meta$Cellline <- unlist(lapply(X=as.character(rownames(meta)), FUN = function(x){return(strsplit(x,split="_")[[1]][1])}))
#meta$Cluster <- paste(meta$Cellline,meta$Cluster,sep = "_")

#meta$Cellline<-gsub("SKHep","SK-HEP-1",meta$Cellline)
rownames(meta) <- gsub("SKHep","SK-Hep",rownames(meta))
rownames(meta) <- gsub("Huh7","Huh-7",rownames(meta))
rownames(meta) <- gsub(".Huh-7",".Huh7",rownames(meta))

meta <- meta[colnames(plot_data1),]
all(rownames(meta)==colnames(plot_data1))
meta$Cluster=paste0("C",meta$Cluster)
head(meta)

### annotation_col_order
annotation_col <- meta[,c(1,3)]
head(annotation_col)
annotation_col$Cellline <- factor(annotation_col$Cellline,level=c("Huh7","HepG2","MHCC97L","MHCC97H","SKHep"),ordered=TRUE)
annotation_col$potential <- factor(annotation_col$Cellline,level=c("Huh7","HepG2","MHCC97L","MHCC97H","SKHep"),labels=1:5,ordered=TRUE)
#annotation_col$cluster <- factor(annotation_col$cluster,level=c(0:10))
annotation_col_order <- annotation_col[order(annotation_col$potential,annotation_col$Cluster),]
annotation_col_order <- annotation_col_order[,-3]

#annotation_col_order$Cellline<-factor(annotation_col_order$cellline)
annotation_col_order$Cluster <- factor(annotation_col_order$Cluster)
head(annotation_col_order)
table(annotation_col_order$Cluster)
table(annotation_col_order$Cellline)

plot_data2 <- plot_data1[,rownames(annotation_col_order)]
all(rownames(annotation_col_order)==colnames(plot_data2))

### ann_colors
#c(Huh7 = "#E41A1C", HepG2 = "#377EB8", MHCC97L = "#4DAF4A", MHCC97H = "#984EA3", SKHEP1 = "#FF7F00"), #
ann_colors = list(
  Cellline = c(Huh7 = "#E41A1C", HepG2 = "#377EB8", MHCC97L = "#4DAF4A", MHCC97H = "#984EA3", SKHep = "#FF7F00"),#brewer.pal(5,"Set1"),
  Cluster = c(C0="#8DD3C7", C1="#FFFFB3", C2="#BEBADA", C3="#FB8072", C4="#80B1D3", C5="#FDB462", C6="#B3DE69")#brewer.pal(7,"Set3")
)



### gsva feature, E_EMT, M_EMT, Hypoxia, ki-67
setwd("../00.genesets")
Angiogenesis <- read.table("Angiogenesis.tsv",header=F)
Check_point <- read.table("Check_point.tsv",header=F)
EMT <- read.table("EMT.tsv",header=F)
Immune_escape <-read.table("Immune_escape.tsv",header=F)
Immune_surveilance <-read.table("Immune_surveilance.tsv",header=F)
stem <- read.table("STEMNESS.tsv",header=F)

Angiogenesis_data <- snRNA@assays$RNA@data[as.vector(Angiogenesis$V1[Angiogenesis$V1%in%rownames(snRNA@assays$RNA@data)]),]
Angiogenesis_data_ave <- colMeans(as.matrix(Angiogenesis_data))

Check_point_data <- snRNA@assays$RNA@data[as.vector(Check_point$V1[Check_point$V1%in%rownames(snRNA@assays$RNA@data)]),]
Check_point_data_ave <- colMeans(as.matrix(Check_point_data))

emt_data <- snRNA@assays$RNA@data[as.vector(EMT$V1[EMT$V1%in%rownames(snRNA@assays$RNA@data)]),]
emt_data_ave <- colMeans(as.matrix(emt_data))

Immune_escape_data <- snRNA@assays$RNA@data[as.vector(Immune_escape$V1[Immune_escape$V1%in%rownames(snRNA@assays$RNA@data)]),]
Immune_escape_data_ave <- colMeans(as.matrix(Immune_escape_data))

Immune_surveilance_data <- snRNA@assays$RNA@data[as.vector(Immune_surveilance$V1[Immune_surveilance$V1%in%rownames(snRNA@assays$RNA@data)]),]
Immune_surveilance_data_ave <- colMeans(as.matrix(Immune_surveilance_data))

stem_data <- snRNA@assays$RNA@data[as.vector(stem$V1[stem$V1%in%rownames(snRNA@assays$RNA@data)]),]
stem_data_ave <- colMeans(as.matrix(stem_data))

E_EMT <- read.table("E_EMT.tsv",sep = "\t")
E_EMT_data <- snRNA@assays$RNA@data[as.vector(E_EMT$V1[E_EMT$V1%in%rownames(snRNA@assays$RNA@data)]),]
E_EMT_data_ave <- colMeans(as.matrix(E_EMT_data))

M_EMT <- read.table("M_EMT.tsv",sep = "\t")
M_EMT_data <- snRNA@assays$RNA@data[as.vector(M_EMT$V1[M_EMT$V1%in%rownames(snRNA@assays$RNA@data)]),]
M_EMT_data_ave <- colMeans(as.matrix(M_EMT_data))

Hypoxia <- read.table("Huh7.C97L.C97H_Hypoxia_2.txt",sep = "\t")
Hypoxia_data <- snRNA@assays$RNA@data[as.vector(Hypoxia$V1[Hypoxia$V1%in%rownames(snRNA@assays$RNA@data)]),]
Hypoxia_data_ave <- colMeans(as.matrix(Hypoxia_data))

genes <- read.table("metabolism_pathway_gene_list.txt",sep="\t",header=T,stringsAsFactors=F)
geneGly <- genes[genes$Type %in% "Glycolysis","Gene"]
geneTCA <- genes[genes$Type %in% "TCA_cycle","Gene"]
Gly_data <- snRNA@assays$RNA@data[as.vector(geneGly[geneGly %in% rownames(snRNA@assays$RNA@data)]),]
Gly_data_ave <- colMeans(as.matrix(Gly_data))
TCA_data <- snRNA@assays$RNA@data[as.vector(geneTCA[geneTCA %in% rownames(snRNA@assays$RNA@data)]),]
TCA_data_ave <- colMeans(as.matrix(TCA_data))

celllines <- sapply(strsplit(colnames(snRNA@assays$RNA@data),"_"),"[",1)
vlndata <- cbind(Gly_data_ave,TCA_data_ave,celllines)
vlndata <- vlndata[vlndata[,"celllines"]=="Huh-7" | vlndata[,"celllines"]=="MHCC97L" | vlndata[,"celllines"]=="MHCC97H",]
vlndata[vlndata[,"celllines"]=="Huh-7","celllines"] <- "Huh7"
vlndata <- as.data.frame(vlndata)
vlndata$celllines <- factor(vlndata$celllines,levels=c("Huh7","MHCC97L","MHCC97H"))


plot_data3 <- rbind(Angiogenesis_data_ave,Check_point_data_ave,emt_data_ave,Immune_escape_data_ave,Immune_surveilance_data_ave,stem_data_ave,E_EMT_data_ave,M_EMT_data_ave,Hypoxia_data_ave,Gly_data_ave,TCA_data_ave,snRNA@assays$RNA@data["MKI67",])
rownames(plot_data3) <- c("Angiogenesis","Check_point","EMT","Immune_escape","Immune_surveilance","stem","E_EMT","M_EMT","Hypoxia","Gly","TCA","MKI67")
dim(plot_data3)

plot_data4 <- plot_data3[,rownames(annotation_col_order)]
all(rownames(annotation_col_order)==colnames(plot_data4))

#plot_data4[plot_data4 < -1] = -1 
#plot_data4[plot_data4 > 1] = 1 

bk <- unique(c(seq(-4,4, length=100)))

plot_data4 <- apply(plot_data4,2,as.numeric)
rownames(plot_data4) <- rownames(plot_data3)
pdf("./cellline_cluster_feature12_ClusterCell.pdf",12,6)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks=bk, color=plasma(100),fontsize_row=10,fontsize_col=15, cluster_row=F, cluster_col=F,scale="row", show_colnames=F)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks=bk, color=plasma(100),fontsize_row=10, fontsize_col=15, cluster_row=F, cluster_col=T, scale="row", show_colnames=F)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks=bk, color=colorRampPalette(colors = c("blue","white","red"))(100),fontsize_row=10,fontsize_col=15, cluster_row=F, cluster_col=F,scale="row", show_colnames=F)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks=bk, color=colorRampPalette(colors = c("blue","white","red"))(100),fontsize_row=10, fontsize_col=15, cluster_row=F, cluster_col=T, scale="row", show_colnames=F)
dev.off()

png("./cellline_cluster_feature12_ClusterCell_bluered.png",res=100,width=800, height=400)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors=ann_colors,color=colorRampPalette(colors = c("blue","white","red"))(100),fontsize_row=10,fontsize_col=15, cluster_row=F, cluster_col=F,scale="row", show_colnames=F)
dev.off()


png("./cellline_cluster_feature12_ClusterCell_plasma.png",res=100,width=800, height=400)
pheatmap(plot_data4, annotation_col = annotation_col_order,annotation_colors=ann_colors,color=plasma(100),fontsize_row=10,fontsize_col=15, cluster_row=F, cluster_col=F,scale="none", show_colnames=F)
dev.off()





all(colnames(plot_data4)==rownames(annotation_col_order))
cellline_violin_data <-  t(plot_data4)
cellline_violin_data <- as.data.frame(cellline_violin_data)
cellline_violin_data$Cellline <- annotation_col_order$Cellline

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(9)

for (i in colnames(cellline_violin_data)[1:9]){
  p1<-ggplot(cellline_violin_data,aes(x=factor(Cellline),y=cellline_violin_data[,i],fill=Cellline,color=Cellline))+
    geom_violin(fill=NA,scale="width")+
    geom_boxplot(fill=NA,width=0.1)+
    #geom_jitter(shape=21)+
    labs(x="",y=i)+
    theme_classic()+
    #scale_y_continuous(limits = c(2,4))+
    scale_fill_manual(values = a)+
    scale_color_manual(values = a)+
    theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
    guides(fill=FALSE)
  
  pdf(paste0("cellline_",i,"_violin.pdf"),10,6)
  print(p1)
  dev.off()
} 




### E,E-M,M
### annotation_col_order
EM_info <- read.table("./5CellLines_E-M.stat.txt",sep = "\t",header=T)
EM_info$sample<-paste(EM_info$cell,EM_info$sample,sep="_")
EM_info$sample<-gsub("SK-HEP-1_","SK-Hep_",EM_info$sample)
EM_info$sample<-gsub("Huh7_","Huh-7_",EM_info$sample)
rownames(EM_info)<-EM_info$sample

annotation_col<-EM_info[,c(2,3)]
annotation_col$cell<-gsub("SK-HEP-1","SKHep",annotation_col$cell)
head(annotation_col)
annotation_col$cell<-factor(annotation_col$cell,level=c("Huh7","HepG2","MHCC97L","MHCC97H","SKHep"),ordered=TRUE)
annotation_col$potential<-factor(annotation_col$cell,level=c("Huh7","HepG2","MHCC97L","MHCC97H","SKHep"),labels=1:5,ordered=TRUE)
#annotation_col$cluster<-factor(annotation_col$cluster,level=c(0:10))
annotation_col_order<-annotation_col[order(annotation_col$stat,annotation_col$potential),]
annotation_col_order<-annotation_col_order[,-3]

plot_data5 <- plot_data3[,rownames(annotation_col_order)]
all(rownames(annotation_col_order)==colnames(plot_data5))

ann_colors = list(
  cell = c(Huh7 = "#E41A1C", HepG2 = "#377EB8", MHCC97L = "#4DAF4A", MHCC97H = "#984EA3", SKHep = "#FF7F00"),#brewer.pal(5,"Set1"),
  stat = c(E="#8DD3C7", EM="#FFFFB3", M="#BEADA")#brewer.pal(7,"Set3")
)

bk <- unique(c(seq(-4,4, length=100)))

pdf("EM_feature9.pdf",12,6)
#par(mai=c(0.5,2,2,2))
#colorRampPalette(colors=c("dark blue","blue","royal blue","sky blue","white","salmon","red","maroon","black"))(100)
#colorRampPalette(colors = c("blue","white","red"))(100)
pheatmap(plot_data5, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks = bk, color = plasma(100),fontsize_row=10,fontsize_col = 15, cluster_row = F, cluster_col = F,scale = "row", show_colnames = F)
pheatmap(plot_data5, annotation_col = annotation_col_order,annotation_colors = ann_colors, breaks = bk, color = colorRampPalette(colors = c("blue","white","red"))(100),fontsize_row=10,fontsize_col = 15, cluster_row = F, cluster_col = F,scale = "row", show_colnames = F)
dev.off()

all(colnames(plot_data5)==rownames(annotation_col_order))
EM_violin_data <-  t(plot_data5)
EM_violin_data <- as.data.frame(EM_violin_data)
EM_violin_data$stat <- annotation_col_order$stat

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a = getPalette(9)

for (i in colnames(EM_violin_data)[1:9]){
  p1<-ggplot(EM_violin_data,aes(x=factor(stat),y=EM_violin_data[,i],fill=stat,color=stat))+
    geom_violin(fill=NA,scale="width")+
    geom_boxplot(fill=NA,width=0.1)+
    #geom_jitter(shape=21)+
    labs(x="",y=i)+
    theme_classic()+
    #scale_y_continuous(limits = c(2,4))+
    scale_fill_manual(values = a)+
    scale_color_manual(values = a)+
    theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
    guides(fill=FALSE)
  
  pdf(paste0("EM_",i,"_violin.pdf"),8,6)
  print(p1)
  dev.off()
} 

