subclustering<-function(combined.c1){
  #DefaultAssay(combined.c1) <- "RNA"
  #Going through the guided-clustering tutorial's workflow
  #combined.c1 <- NormalizeData(combined.c1)
  #combined.c1 <- FindVariableFeatures(combined.c1)
  DefaultAssay(combined.c1)<-"integrated"
  combined.c1 <- ScaleData(combined.c1)
  combined.c1 <- RunPCA(combined.c1)
  # Elbow plot to identify dimensionality
  ElbowPlot(combined.c1)
  # Use the elbow from above to sub-cluster the cells
  combined.c1 <- FindNeighbors(combined.c1, dims = 1:14)
  combined.c1 <- FindClusters(combined.c1, resolution = 0.1)
  
  return(combined.c1 <- RunUMAP(combined.c1, dims = 1:10))
}
inhibitory_subcluster<-subclustering(Inhibitory_seurat_object)

DimPlot(inhibitory_subcluster, reduction = "umap", split.by = "condition",label = T, pt.size=0.5)+ theme(text = element_text(size=25, face="bold"))
#Plot heatmap with Common marker for inhibitory neurons subclass

DoHeatmap(inhibitory_subcluster,features=c("Lamp5","Sv2c","Vip","Sst","Pvalb","Calb1","Cck","Calb2","Sv2a","Sv2b","Inpp4b"), size=5, angle=90, label=T,slot="scale.data", assay="integrated")+theme(text=element_text(size=15, face="bold"))

markers_inhibitory<- FindAllMarkers(inhibitory_subcluster,thresh.use = 0.25,only.pos = TRUE)
inhibitory_subcluster<-RenameIdents(inhibitory_subcluster, "0"="Sv2b and Cck",
                                    "1"="Pvalb",
                                    "2"="Sst",
                                    "3"="Vip",
                                    "4"="Lamp5 and Sv2c",
                                    "5"="Others",
                                    "6"="Others",
                                    "7"="Inpp4b",
                                    "8"="Sst")


color1<-c("Pvalb"="orange","Sst"="green", "Vip"="blue","Lamp5 and Sv2c"="magenta","Sv2b and Cck"="red","Others"="grey","Inpp4b"= "cyan")


order1 = c("Pvalb","Sst","Vip","Lamp5 and Sv2c","Sv2b and Cck","Inpp4","Others")

heatmap_inhibit<-DoHeatmap(inhibitory_subcluster,features=c("Lamp5","Sv2c","Vip","Sst","Pvalb","Calb1","Cck","Calb2","Sv2a","Sv2b","Inpp4b"), size=10, angle=90, label=T,slot="scale.data", assay="integrated", group.colors = color1)+theme(text=element_text(size=20, face="bold"))

inhibit<-DimPlot(inhibitory_subcluster, reduction = "umap",label = T,label.size = 8, pt.size=1, cols=color1, order=order1)+theme(text = element_text(size=25, face="bold"))
ggsave("inhibitory_dim.tiff", ggplotify::as.ggplot(inhibit), dpi=1200, width=40, height=30, unit="cm")
ggsave("heatmap_inhibit.tiff", ggplotify::as.ggplot(heatmap_inhibit), dpi=1200, width=40, height=32, unit="cm")

##Proportion of cells of inhibitory neurons
table(inhibitory_subcluster@meta.data$condition)
prop.table(table(Idents(inhibitory_subcluster)))
table(Idents(inhibitory_subcluster), inhibitory_subcluster@meta.data$condition)
tab<-prop.table(table(Idents(inhibitory_subcluster), inhibitory_subcluster@meta.data$condition), margin = 2)
tab<-as.data.frame(tab)
names(tab)<-c("cell_type","condition","prop")
tab$prop<-signif(tab$prop,3)
tab

