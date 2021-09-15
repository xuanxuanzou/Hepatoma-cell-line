
library(pheatmap)
library(data.table)

ATAC <- read.table("scATAC5_CellLines_GeneActivityScore.txt.gz",sep="\t",stringsAsFactors=F,comment.char="")

ann_colors = list(
  Class = c(Huh7=hue_pal()(5)[1], HepG2=hue_pal()(5)[2], MHCC97L=hue_pal()(5)[3], MHCC97H=hue_pal()(5)[4], SK.Hep=hue_pal()(5)[5]))


################# 1. sample30 #################
set.seed(6)
ATACda <- ATAC[,c(sample(grep("Huh7",colnames(ATAC)),30),sample(grep("HepG2",colnames(ATAC)),30), sample(grep("MHCC97L",colnames(ATAC)),30),sample(grep("MHCC97H",colnames(ATAC)),30),sample(grep("SK.Hep",colnames(ATAC)),30))]
corDa <- cor(ATACda)

annota = data.frame(Class=factor(c(rep("Huh7",30),rep("HepG2",30),rep("MHCC97L",30),rep("MHCC97H",30),rep("SK.Hep",30)),levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK.Hep")))
rownames(annota) = colnames(ATACda)


pdf("ATAC_cor.30cells_Cor.pdf",width=13,height=13)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, color=colorRampPalette(colors=c("dark blue","blue", "royal blue","sky blue","white","salmon","red","maroon","black"))(100), fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
dev.off()



################# 2. sample50 #################
set.seed(6)
ATACda <- ATAC[,c(sample(grep("Huh7",colnames(ATAC)),50),sample(grep("HepG2",colnames(ATAC)),50), sample(grep("MHCC97L",colnames(ATAC)),50),sample(grep("MHCC97H",colnames(ATAC)),50),sample(grep("SK.Hep",colnames(ATAC)),50))]
corDa <- cor(ATACda)

annota = data.frame(Class=factor(c(rep("Huh7",50),rep("HepG2",50),rep("MHCC97L",50),rep("MHCC97H",50),rep("SK.Hep",50)),levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK.Hep")))
rownames(annota) = colnames(ATACda)


pdf("ATAC_cor.50cells_Cor.pdf",width=25,height=25)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, color=colorRampPalette(colors=c("dark blue","blue", "royal blue","sky blue","white","salmon","red","maroon","black"))(100), fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
dev.off()




################# 5. sample20 #################
set.seed(6)
ATACda <- ATAC[,c(sample(grep("Huh7",colnames(ATAC)),20),sample(grep("HepG2",colnames(ATAC)),20), sample(grep("MHCC97L",colnames(ATAC)),20),sample(grep("MHCC97H",colnames(ATAC)),20),sample(grep("SK.Hep",colnames(ATAC)),20))]
corDa <- cor(ATACda)

annota = data.frame(Class=factor(c(rep("Huh7",20),rep("HepG2",20),rep("MHCC97L",20),rep("MHCC97H",20),rep("SK.Hep",20)),levels=c("Huh7","HepG2","MHCC97L","MHCC97H","SK.Hep")))
rownames(annota) = colnames(ATACda)


pdf("ATAC_cor.20cells_Cor.pdf",width=8,height=8)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, color=colorRampPalette(colors=c("dark blue","blue", "royal blue","sky blue","white","salmon","red","maroon","black"))(100), fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
pheatmap(corDa,scale="none", cluster_rows=T,cluster_cols=T, border_color=NA, fontsize_row=6, fontsize_col=6, angle_col=90, annotation_col=annota, annotation_row=annota, annotation_colors=ann_colors)
dev.off()
