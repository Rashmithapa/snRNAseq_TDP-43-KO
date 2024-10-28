##Secreted inflammatory molecules
secreted<-c("Il1a","Il1b","C1qa","C1qb","Tnf","Tnfa","C1q","Csf2","Il6","Il1b","Cspg4" , "Cxcl1","Ccl2","Ccl5","Cxcl8","Il17","Il17a","Cxcl10","Il10","Tgfa","Tgfb1","C1ra","C2","C3","C3a","C3b","C4","C4a","C4b","C5","C5a","C5b","C6","C7","C8", "Nos2","Mmp12","Infg","Il2","Il4","Vegfb")
#cre-negative versus cre-positive in KO
divide<-Astrocytes_seurat_object
subset0<-subset(divide,subset= condition%in%"2m_KO")
subset1<-subset(divide,subset= condition%in%"2m_control") #|condition%in%"FAD-positive"|condition%in%"FAD-negative")
DimPlot(subset0)
FeaturePlot(subset0, features = "cre-transgene")

cre_expression = GetAssayData(object = subset0, assay = "RNA", slot = "data")["cre-transgene",]

pos_cre = names(which(cre_expression>0))
neg_cre = names(which(cre_expression<=0))

#Subsetting the cre positive and negative cells
cre_pos_cells = subset(subset0,cells=pos_cre)
cre_neg_cells = subset(subset0,cells=neg_cre)

Idents(cre_pos_cells) <- 'cre-positive'
Idents(cre_neg_cells)<- 'cre-negative'
#merge the two seurat object
cre_merged=merge(x = cre_pos_cells, y = c(cre_neg_cells, subset1))

cre_merged$cre<-Idents(cre_merged)

secreted<-c("Il1a","Il1b","C1qa","C1qb","Tnf","Tnfa","C1q","Csf2","Il6","Il1b","Cspg4" , "Cxcl1","Ccl2","Ccl5","Cxcl8","Il17","Il17a","Cxcl10","Il10","Tgfa","Tgfb1","C1ra","C2","C3","C3a","C3b","C4","C4a","C4b","C5","C5a","C5b","C6","C7","C8", "Nos2","Mmp12","Infg","Il2","Il4","Vegfb")
#synaptogenic<-c("Gpc4", "Gpc6","Sparc","Sparcl1", "Thbs1", "Thbs2")

avgexp_sec=AverageExpression(cre_merged, assays = "RNA",slot = "data", group.by= c("cre"),features = secreted)

avgexp_sec=as.matrix(avgexp_sec$RNA)

aa<-pheatmap(
  avgexp_sec,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",  # Set border_color to the desired color
  fontsize = 20,
  scale="row")
ggsave("avgsec_astro.pdf", ggplotify::as.ggplot(aa), dpi=1200, width=10, height=25, unit="cm")

##Disease associated astrocytes genes in astrocytes within the 6 samples
astro_genes<-c("Nrxn1","Nrg3","Gpc5","Erbb4","Ntm", "Dclk1", "Slc1a3","Gria2","Luzp2", "Slc7a10", "Mfge8","Id3", "Aqp4", "Myoc", "Id1", "Fabp7", "Ctsb", "Vim", "Osmr", "Serpina3n", "Gsn", "Ggta1","Cst3", "Gfap", "Vim", "Clu")


avgexp_astro=AverageExpression(Astrocytes_seurat_object, assays = "RNA",slot = "data", group.by= c("condition"),features = astro_genes)

avgexp_astro=as.matrix(avgexp_astro$RNA)
#avgexp_daa<- scale(avgexp_daa)

aa<-pheatmap(
  avgexp_astro,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  fontsize = 14,scale="row")


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
  combined.c1 <- FindClusters(combined.c1, resolution = 0.3)
  
  return(combined.c1 <- RunUMAP(combined.c1, dims = 1:10))
}

astrocytes_subcluster<-subclustering(Astrocytes_seurat_object)

astro_subcluster<-DimPlot(astrocytes_subcluster, reduction = "umap", split.by = "condition",label = T, pt.size=1,label.size = 4)
ggsave("astro_subcluster.tiff", ggplotify::as.ggplot(astro_subcluster), dpi=1200, width=50, height=15, unit="cm")



pan_marker<-c("Lnc2", "Steap4", "S1pr3","Timp1","Hspb1","Cxcl10", "Cd44", "Osmr", "Cp", "Serpina3n", "Aspg", "Vim", "Gfap")
a1_marker<- c("H2-T23", "Serping1", "H2-D1", "Ggta1", "Iigp1" , "Gbp2", "Fkbp5", "Ugt1a", "Psmb8", "Srgn", "Amigo2", "Fbln5")
a2_marker<-c("Clcf1", "Tgm1", "Ptx3", "S100a10", "Sphk1", "Cd109", "Ptgs2", "Emp1", "Slc10a6", "Tm4sf1", "B3gnt5", "Cd14")

Astrocytes_seurat_object<- AddModuleScore(object = Astrocytes_seurat_object, features = list(a1_marker), name = "A1_score", assay="RNA")
Astrocytes_seurat_object<- AddModuleScore(object = Astrocytes_seurat_object, features = list(a2_marker), name = "A2_score", assay="RNA")
Astrocytes_seurat_object<- AddModuleScore(object = Astrocytes_seurat_object, features = list(pan_marker), name = "PAN_score", assay="RNA")
Astrocytes_seurat_object<- AddModuleScore(object = Astrocytes_seurat_object, features = list(hom_astro), name = "hom_astro", assay="RNA")

astrocytes_subcluster<- AddModuleScore(object = astrocytes_subcluster, features = list(a1_marker), name = "A1_score", assay="RNA")
astrocytes_subcluster<- AddModuleScore(object = astrocytes_subcluster, features = list(a2_marker), name = "A2_score", assay="RNA")
astrocytes_subcluster<- AddModuleScore(object = astrocytes_subcluster, features = list(pan_marker), name = "PAN_score", assay="RNA")


pan<-FeaturePlot(astrocytes_subcluster, features = "PAN_score1", split.by="condition",pt.size=2, cols= rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(legend.position = c(0.93,0.95))+
  patchwork::plot_layout(ncol = 3, nrow = 2)
ggsave("pan.png", ggplotify::as.ggplot(pan), dpi=300, width=50, height=30, unit="cm")

a1<-FeaturePlot(astrocytes_subcluster, features = "A1_score1", split.by="condition",pt.size = 2, cols= rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(legend.position = c(0.93,0.95))+
  patchwork::plot_layout(ncol = 3, nrow = 2)
ggsave("a1.png", ggplotify::as.ggplot(a1), dpi=300, width=50, height=30, unit="cm")


a2<-FeaturePlot(astrocytes_subcluster, features = "A2_score1", split.by="condition",pt.size=2, cols= rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(legend.position = c(0.93,0.95))+
  patchwork::plot_layout(ncol = 3, nrow = 2)

ggsave("a2.png", ggplotify::as.ggplot(a2), dpi=300, width=50, height=30, unit="cm")hom_astro<-c("Nrxn1","Nrg3","Gpc5","Erbb4","Ntm", "Dclk1", "Slc1a3","Gria2","Luzp2", "Slc7a10", "Mfge8")

