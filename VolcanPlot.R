library(Seurat)
library(ggplot2)
library(ggrepel)


### 1.Huh7
Huh7 <- readRDS("./Huh7.rds")

cluster.markers <- FindMarkers(Huh7, ident.1 = 3,logfc.threshold=0 , min.pct = 0)
cluster.markers$log2FoldChange <- log2(exp(cluster.markers$avg_logFC))
#cluster.markers[which(cluster.markers$p_val_adj < 1e-50),"p_val_adj"] <- 1e-50
cluster.markers[which(cluster.markers$p_val_adj == 1),"p_val_adj"] <- as.numeric(substr(paste0("1e-",runif(length(which(cluster.markers$p_val_adj == 0)),min=0.1,max=0.5)),1,6))
cluster.markers$color <- ifelse(cluster.markers$p_val_adj<0.05 & abs(cluster.markers$log2FoldChange)>= 1,ifelse(cluster.markers$log2FoldChange > 1,'red','blue'),'grey')
color <- c(red = "red",blue="blue",grey="grey")
genes <- c("CA9","ENO2","SLC6A8","BNIP3","FAM162A","BNIP3L","INSIG2","NDRG1","LDHA","PLOD2","ALDOC","ANGPTL4","ZNF395","HILPDA")
cluster.markers$label = ""
cluster.markers[rownames(cluster.markers) %in% genes,"label"] <- rownames(cluster.markers)[rownames(cluster.markers) %in% genes]
p <- ggplot(cluster.markers, aes(log2FoldChange,-log10(p_val_adj), col = color)) +geom_point() +theme_bw() +
  scale_color_manual(values = color) +labs(x="log2 fold change",y="-log10 (adjusted P)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="black",lwd=0.6) +
  theme(legend.position="none",panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=18))
p1 <- p+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=5,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),segment.color = "black",show.legend = FALSE)
pdf("./Huh7.Volcan.pdf")
print(p1)
dev.off()
png("./Huh7.Volcan.png",type = "cairo",res=500,,width = 3000, height = 3500)

print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank()))
dev.off()
png("./Huh7.Volcan1.png",type = "cairo",res=500,,width = 3000, height = 3500)
print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank())+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=0,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),segment.color = "black",show.legend = FALSE))
dev.off()
write.table(cluster.markers,"./Huh7.Volcan.txt",col.names=T,row.names=T,sep="\t",quote=F)



############ 2.97L
C97L <- readRDS("./97L.rds")

cluster.markers <- FindMarkers(C97L, ident.1 = 6,logfc.threshold=0 , min.pct = 0)
cluster.markers$log2FoldChange <- log2(exp(cluster.markers$avg_logFC))
#cluster.markers[which(cluster.markers$p_val_adj < 1e-50),"p_val_adj"] <- 1e-50
cluster.markers[which(cluster.markers$p_val_adj == 1),"p_val_adj"] <- as.numeric(substr(paste0("1e-",runif(length(which(cluster.markers$p_val_adj == 0)),min=0.1,max=0.5)),1,6))
cluster.markers$color <- ifelse(cluster.markers$p_val_adj<0.05 & abs(cluster.markers$log2FoldChange)>= 1,ifelse(cluster.markers$log2FoldChange > 1,'red','grey'),'grey')
color <- c(red = "red",grey="grey",grey="grey")
genes <- c("CA9","ENO2","SLC6A8","BNIP3","FAM162A","BNIP3L","INSIG2","NDRG1","LDHA","PLOD2","ALDOC","ANGPTL4","ZNF395","HILPDA")
cluster.markers$label = ""
cluster.markers[rownames(cluster.markers) %in% genes,"label"] <- rownames(cluster.markers)[rownames(cluster.markers) %in% genes]
p <- ggplot(cluster.markers, aes(log2FoldChange,-log10(p_val_adj), col = color)) +geom_point() +theme_bw() +
  scale_color_manual(values = color) +labs(x="log2 fold change",y="-log10 (adjusted P)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="black",lwd=0.6) +
  theme(legend.position="none",panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=18))
p1 <- p+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=5, box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),segment.color = "black",show.legend = FALSE)
pdf("./C97L.Volcan.pdf")
print(p1)
dev.off()
png("./C97L.Volcan.png",type = "cairo",res=500,width = 3000, height = 3500)
print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank()))
dev.off()

png("./C97L.Volcan1.png",type = "cairo",res=500,width = 3000, height = 3500)
print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank())+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=0,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),segment.color = "black",show.legend = FALSE))
dev.off()
write.table(cluster.markers,"./C97L.Volcan.txt",col.names=T,row.names=T,sep="\t",quote=F)



############ 3.97H
C97H <- readRDS("./97H.rds")

cluster.markers <- FindMarkers(C97H, ident.1=5, logfc.threshold=0 , min.pct = 0)
cluster.markers$log2FoldChange <- log2(exp(cluster.markers$avg_logFC))
#cluster.markers[which(cluster.markers$p_val_adj < 1e-50),"p_val_adj"] <- 1e-50
cluster.markers[which(cluster.markers$p_val_adj == 1),"p_val_adj"] <- as.numeric(substr(paste0("1e-",runif(length(which(cluster.markers$p_val_adj == 0)),min=0.1,max=0.5)),1,6))
cluster.markers$color <- ifelse(cluster.markers$p_val_adj<0.05 & abs(cluster.markers$log2FoldChange)>= 1,ifelse(cluster.markers$log2FoldChange > 1,'red','grey'),'grey')
color <- c(red = "red",grey="grey",grey="grey")
genes <- c("CA9","ENO2","SLC6A8","BNIP3","FAM162A","BNIP3L","INSIG2","NDRG1","LDHA","PLOD2","ALDOC","ANGPTL4","ZNF395","HILPDA")
cluster.markers$label = ""
cluster.markers[rownames(cluster.markers) %in% genes,"label"] <- rownames(cluster.markers)[rownames(cluster.markers) %in% genes]
p <- ggplot(cluster.markers, aes(log2FoldChange,-log10(p_val_adj), col = color)) +geom_point() +theme_bw() +
  scale_color_manual(values = color) +labs(x="log2 fold change",y="-log10 (adjusted P)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="black",lwd=0.6) +
  theme(legend.position="none",panel.grid=element_blank(),axis.title=element_text(size=20),axis.text=element_text(size=18))
p1 <- p+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=5,box.padding = unit(0.5, "lines"),point.padding = unit(0.5, "lines"),segment.color = "black",show.legend = FALSE)
pdf("./C97H.Volcan.pdf")
print(p1)
dev.off()
png("./C97H.Volcan.png",type = "cairo",res=500,width = 3000, height = 3500)
print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank()))
dev.off()
png("./C97H.Volcan1.png",type = "cairo",res=500,width = 3000, height = 3500)
print(p+theme(axis.text = element_blank())+theme(axis.title = element_blank())+geom_text_repel(data = cluster.markers, aes(x=cluster.markers$log2FoldChange,y = -log10(cluster.markers$p_val_adj),label = label),size=0,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"),segment.color = "black",show.legend = FALSE))
dev.off()
write.table(cluster.markers,"./C97H.Volcan.txt",col.names=T,row.names=T,sep="\t",quote=F)

