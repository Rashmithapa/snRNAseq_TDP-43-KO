m1_micro<-c("Irg1","Ccl2","Csf2","Marcksl1","Il23a","Cmpk2","Rtp4","Slc7a11","Maff","Gpr85","Pde4b","Trex1","Ier3","Ifi47","Nfkbie","Oasl2","Gdf15","Cdkn1a","Slc7a2","Gem","F10","Oas2","Il15","Rnd3","Gbp2","Lif","Irf1","Fabp3","Src","Slc39a14")

m2_micro<-c("Fbxo32","Cebpa","Gadd45g","N4bp2l1","Tns1","Lpin1","Mnt","Rgs2","Gprc5b","Rassf3","Rasa3","Ypel3","Mknk2","Pparg","Tpcn1","Pdcd4","Arhgap18","Dnmt1","Gmnn","Rapgef6","Idh1","Klc4","Mxi1","Dennd4c","H2afv","St6gal1","Ldlrap1","Mef2d","Sgsh","Pcmtd2")


proliferate_micro<-c("Cd86", "Egr1", "Pdgfb", "Fos", "Ptgs2", "Cxcl2", "F3", "Nfkb1", "Socs3", "Bcl6", "Il1b", "Ccl4", "Il12b", "Tnfsf9", "Junb","Il1a")

phagocytic_micro<-c("C1qa", "C1qb", "C1qc", "Tyrobp", "Trem2", "Aif1", "B2m", "Prdx5", "Fcer1g", "Cstb", "Ctsz", "Cd63",  "Cd68","Fau","Ftl1")

ifn1_micro<-c("Trim30a", "Oasl2", "Cxcl10", "Iftm3", "Bst2", "Rsad2", "Isg15", "Ift1", "Gbp2", "Ift3", "Ift2", "Cxcl10")

apc_micro<-c("Cd74", "H2-Aa", "Cd52" , "Ccl6", "Cd74"," H2-Eb1", "H2-Ab1", "H2-K1",  "H2-D1")
DAM_markers<-c("Apoe", "Ctsd", "Lpl", "Tyrobp", "Trem2")
hom_genes<-c("P2ry12", "P2ry13", "Cx3cr1", "CD33","Tmem119","Nav2")

Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(m1_micro), name = "m1_microglia", assay="RNA")
Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(m2_micro), name = "m2_microglia", assay="RNA")
Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(DAM_markers), name = "DAM_markers", assay="RNA")
Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(hom_genes), name = "HOM_genes", assay="RNA")
Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(proliferate_micro), name = "proliferate_micro", assay="RNA")

Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(phagocytic_micro), name = "phagocytic_micro", assay="RNA")

Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(ifn1_micro), name = "ifn1_micro", assay="RNA")

Microglia_seurat_object<- AddModuleScore(object = Microglia_seurat_object, features = list(apc_micro), name = "apc_micro", assay="RNA")

##Heatmap for disease-associated and homeostatic microglia
dam<-c("Ifi204","H2-K1","H2-D1","Ttr","Oasl2","B2m","Lyz2","Ifi204","H2-K1","Ncf1","Etl4","Apoe","Itgax", "Tyrobp","Ctsb","Cstd","Cst7","C1qa","Gbp2","Clec7a","Trem2","Axl","Lpl","Spp1","Csf1","Cd9","Ccl6")

avgexp_dam = AverageExpression(Microglia_seurat_object, assays = "RNA",slot = "data", group.by= c("condition"),features = dam)
avgexp_dam=as.matrix(avgexp_dam$RNA)

#avgexp_dam<- scale(avgexp_dam)

avgexp_dam<-Heatmap(avgexp_dam,cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, width = unit(120, "mm"),height = unit(250, "mm"), # Use col instead of lgd for specifying the color function 
                    name = "Expression",  # Optional: Set the legend title
                    #row_title = "Genes",
                    border = TRUE, col =  rev(brewer.pal(11, "RdYlBu")))

#ggsave("avgexp_dam.pdf", ggplotify::as.ggplot(avgexp_dam), dpi=1200, width=20, height=30, unit="cm",limitsize = FALSE)

avgexp_dam<-pheatmap(
  avgexp_dam,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",  # Set border_color to the desired color,
  fontsize = 20,
  scale="row")
ggsave("avgexp_dam.pdf", ggplotify::as.ggplot(avgexp_dam), dpi=1200, width=12, height=30, unit="cm")


hom_micro<-c("P2ry12" , "P2ry13" , "Cx3cr1" , "CD33","Tmem119","Pmepa1","Selplg","Olfml3","Txnip","Maf","Serinc3","Sall1","Tgfbr2","Tgfbr2","Itga6", "Atp8a2","Abcc3","Cmklr1","Gpr34")

avgexp_hom=AverageExpression(Microglia_seurat_object, assays = "RNA",slot = "data", group.by= c("condition"),features = hom_micro)

avgexp_hom=as.matrix(avgexp_hom$RNA)

aa<-pheatmap(
  avgexp_hom,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",  # Set border_color to the desired color
  fontsize = 20,
  scale="row")
ggsave("avgexp_hom.pdf", ggplotify::as.ggplot(aa), dpi=1200, width=15, height=30, unit="cm")
