library("Seurat")
library("dplyr")
library("patchwork")
library("cowplot")
library("AnnotationHub")
library("ensembldb")
library("purrr")
library("ggplot2")
library("clusterProfiler")
library("biomaRt")#For pathway analysis
library("org.Mm.eg.db")
library("tibble")
library("ggupset")
library("enrichplot")
library("DOSE")
library("reshape")
library("tibble")
library("ggpubr")

##Importing data
KO.6mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/filtered_feature_bc_matrix_TDP43_KO_6m/")
ctl.6mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/filtered_feature_bc_matrix_TDP43_control_6m/")
posfad.data<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/filtered_feature_bc_matrix_5Xfad_pos/")
negfad.data<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/filtered_feature_bc_matrix_5Xfad_neg/")
KO1.2mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/KO1_filtered_feature_bc_matrix_2m/")
KO2.2mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/KO2_filtered_feature_bc_matrix_2m/")
TDP.2mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/TDP_filtered_feature_bc_matrix_2m/")
WT.2mdata<-Read10X(data.dir = "/Users/rashmithapa/Library/CloudStorage/OneDrive-UniversityofWyoming/Desktop/single nucleus RNA seq analysis/TDP43_all_samples/all_filtered_feature_matrix/WT_filtered_feature_bc_matrix_2m/")

##Creating Seurat Object and filter to remove gene that appears in less than 1% of cells
KO1_2m<-CreateSeuratObject(counts=KO1.2mdata, project="KO1_2m",min.cells=62, min.features=1)
KO2_2m<-CreateSeuratObject(counts=KO2.2mdata, project="KO2_2m",min.cells=60, min.features=1)
TDP_2m<-CreateSeuratObject(counts=TDP.2mdata, project="TDP_2m",min.cells=55, min.features=1)
WT_2m<-CreateSeuratObject(counts=WT.2mdata, project="WT_2m",min.cells=61, min.features=1)
KO_6m<-CreateSeuratObject(counts=KO.6mdata, project="KO_6m",min.cells=54, min.features=1)
ctl_6m<-CreateSeuratObject(counts=ctl.6mdata, project="ctl_6m",min.cells=51, min.features=1)
posfad<-CreateSeuratObject(counts=posfad.data, project="posfad",min.cells=41, min.features=1)
negfad<-CreateSeuratObject(counts=negfad.data, project="negfad",min.cells=42, min.features=1)

## Add a column in the the seurat object metadata with condition
KO1_2m@meta.data$condition<-"2m_KO"

KO2_2m@meta.data$condition<-"2m_KO"

TDP_2m@meta.data$condition<-"2m_control"

WT_2m@meta.data$condition<-"2m_control"

KO_6m@meta.data$condition<-"6m_KO"

ctl_6m@meta.data$condition<-"6m_control"

posfad@meta.data$condition<-"FAD-positive"

negfad@meta.data$condition<-"FAD-negative"

##Quality control
KO1_2m[["percent.mt"]]<-PercentageFeatureSet(KO1_2m, pattern ="^mt-")
VlnPlot(KO1_2m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

KO2_2m[["percent.mt"]]<-PercentageFeatureSet(KO2_2m, pattern ="^mt-")
VlnPlot(KO2_2m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

WT_2m[["percent.mt"]]<-PercentageFeatureSet(WT_2m, pattern ="^mt-")
VlnPlot(WT_2m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

TDP_2m[["percent.mt"]]<-PercentageFeatureSet(TDP_2m, pattern ="^mt-")
VlnPlot(TDP_2m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

KO_6m[["percent.mt"]]<-PercentageFeatureSet(KO_6m, pattern ="^mt-")
VlnPlot(KO_6m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

ctl_6m[["percent.mt"]]<-PercentageFeatureSet(ctl_6m, pattern ="^mt-")
VlnPlot(ctl_6m, features= c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

posfad[["percent.mt"]]<-PercentageFeatureSet(posfad, pattern ="^mt-")
VlnPlot(posfad, features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

negfad[["percent.mt"]]<-PercentageFeatureSet(negfad, pattern ="^mt-")
VlnPlot(negfad, features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

##Subsetting to filter out the cells with lowest 10% and highest 1% nFeature_RNA and expressing mitochondrial gene percentage less than 10%
KO1_2m<- subset(KO1_2m, subset= nFeature_RNA>562 & nFeature_RNA< 6182 & percent.mt<=10)
KO2_2m<- subset(KO2_2m, subset= nFeature_RNA>683 & nFeature_RNA< 6211 & percent.mt<=10)
WT_2m<- subset(WT_2m, subset= nFeature_RNA>692 & nFeature_RNA< 6357 & percent.mt<=10)
TDP_2m<- subset(TDP_2m, subset= nFeature_RNA>729 & nFeature_RNA< 5903 & percent.mt<=10)
KO_6m<- subset(KO_6m, subset= nFeature_RNA>507 & nFeature_RNA< 5633 & percent.mt<=10)
ctl_6m<- subset(ctl_6m, subset= nFeature_RNA>504 & nFeature_RNA<5641  & percent.mt<=10)
posfad<- subset(posfad, subset= nFeature_RNA>746 & nFeature_RNA< 4073 & percent.mt<=10)
negfad<- subset(negfad, subset= nFeature_RNA>513 & nFeature_RNA< 5901 & percent.mt<=10)

##Removing mt gene
KO1_2m <- KO1_2m[!grepl("^mt-", rownames(KO1_2m)), ]
KO2_2m <- KO2_2m[!grepl("^mt-", rownames(KO2_2m)), ]
WT_2m <- WT_2m [!grepl("^mt-", rownames(WT_2m)), ]
TDP_2m <-TDP_2m[!grepl("^mt-", rownames(TDP_2m)), ]
KO_6m <- KO_6m[!grepl("^mt-", rownames(KO_6m)), ]
ctl_6m <- ctl_6m[!grepl("^mt-", rownames(ctl_6m)), ]
posfad <- posfad[!grepl("^mt-", rownames(posfad)), ]
negfad <- negfad[!grepl("^mt-", rownames(negfad)), ]

##Normalize
KO_6m<-NormalizeData(KO_6m)
ctl_6m<-NormalizeData(ctl_6m)
posfad<-NormalizeData(posfad)
negfad<-NormalizeData(negfad)
KO1_2m<-NormalizeData(KO1_2m)
KO2_2m<-NormalizeData(KO2_2m)
WT_2m<-NormalizeData(WT_2m) 
TDP_2m<-NormalizeData(TDP_2m) 

##Identify top variable feature genes
KO1_2m<-FindVariableFeatures(KO1_2m, selection.method = "vst", nfeatures= 2000)
KO2_2m<-FindVariableFeatures(KO2_2m, selection.method = "vst", nfeatures= 2000)
WT_2m<-FindVariableFeatures(WT_2m, selection.method = "vst", nfeatures= 2000)
TDP_2m<-FindVariableFeatures(TDP_2m, selection.method = "vst", nfeatures= 2000)
KO_6m<-FindVariableFeatures(KO_6m, selection.method = "vst", nfeatures= 2000)
ctl_6m<-FindVariableFeatures(ctl_6m, selection.method = "vst", nfeatures= 2000)
posfad<-FindVariableFeatures(posfad, selection.method = "vst", nfeatures= 2000)
negfad<-FindVariableFeatures(negfad, selection.method = "vst", nfeatures= 2000)

##Integrate and merge samples
merge.list<-list(KO1_2m, KO2_2m, WT_2m, TDP_2m, KO_6m, ctl_6m, posfad, negfad)
features <- SelectIntegrationFeatures(object.list = merge.list)
anchors<- FindIntegrationAnchors(object.list = merge.list, dims=1:30, anchor.features = features)
merge<- IntegrateData(anchorset = anchors, dims=1:30) 
DefaultAssay(merge)<-"integrated"
##Scaling
merge<-ScaleData(merge)

##Dimensionality reduction
merge<- RunPCA(merge, npcs=30)
merge<-RunUMAP(merge, reduction="pca", dims=1:20)
ElbowPlot(object = merge, ndims = 40)
merge<- FindNeighbors(merge, reduction= "pca", dims=1:30)
merge<- FindClusters(object= merge, resolution=1.6 , random.seed = 0)
DimPlot(merge, reduction="umap", group.by= "condition")
DimPlot(merge, reduction= "umap", label=TRUE)
DimPlot(merge, reduction= "umap", split.by = "condition", label=T)

DefaultAssay(merge)<-"RNA"

#Find all markers for each cluster
library("multtest")
library("metap")
all_markers<- FindAllMarkers(merge,thresh.use = 0.25,only.pos = TRUE)

##Cell-type identification
merge<-RenameIdents(object = merge, 
                    "0"=	"Astrocytes",
                    "1"=	"Excitatory neurons",
                    "2"=	"Oligodendrocytes",
                    "3"=	"Excitatory neurons",
                    "4"=	"Inhibitory neurons",
                    "5"=	"Excitatory neurons",
                    "6"= "Excitatory neurons",
                    "7"=	"Excitatory neurons",
                    "8"=	"Excitatory neurons",
                    "9"=	"Excitatory neurons",
                    "10"=	"Excitatory neurons",
                    "11"=	"Microglia",
                    "12"=	"Inhibitory neurons",
                    "13"=	"Excitatory neurons",
                    "14"=	"Excitatory neurons",
                    "15"=	"Excitatory neurons",
                    "16"=	"Excitatory neurons",
                    "17"=	"Excitatory neurons",
                    "18"=	"Excitatory neurons",
                    "19"=	"OPC",
                    "20"=	"Excitatory neurons",
                    "21"=	"Excitatory neurons",
                    "22"=	"Astrocytes",
                    "23"=	"Inhibitory neurons",
                    "24"=	"Excitatory neurons",
                    "25"=	"Inhibitory neurons",
                    "26"=	"Inhibitory neurons",
                    "27"=	"Unknown",
                    "28"=	"Excitatory neurons",
                    "29"=	"Inhibitory neurons",
                    "30"=	"Excitatory neurons",
                    "31"=	"Excitatory neurons",
                    "32"=	"Unknown",#Mixed
                    "33"=	"Inhibitory neurons",
                    "34"=	"Astrocytes",
                    "35"=	"Oligodendrocytes",
                    "36"=	"Endothelial cells",
                    "37"=	"Inhibitory neurons",
                    "38"=	"Unknown",#Mixed
                    "39"=	"Excitatory neurons",
                    "40"=	"Inhibitory neurons",
                    "41"=	"Unknown",#Mixed
                    "42"=	"Excitatory neurons",
                    "43"=	"Excitatory neurons",
                    "44"=	"Excitatory neurons",
                    "45"=	"Unknown",#Mixed
                    "46"=	"Inhibitory neurons")



##Cre-feature plot
cre_feature<-FeaturePlot(copy, features = "cre-transgene", pt.size=0.5, slot="data", cols = alpha( c("lightblue","grey", "red"), 1), order=T)+ theme(legend.position = c(0.93, 0.95), axis.title = element_text(size=20), axis.text =element_text(size=15) , plot.caption = element_text(size=20)) 
ggsave("cre_feature.tiff", ggplotify::as.ggplot(cre_feature), dpi=1200, width=30, height=20, unit="cm")

###HEatmap fpr the marker genes
markers<-c("Aqp4", "Slc1a2","Gja1","Gfap", "Gad1", "Gad2","Sst","Npy","Vip","Pvalb","Lamp5","Cx3cr1","Csf1r","Hexb","Havcr2","Cldn11","Plp1","Mbp","Mog","Sept4","Pdgfra","Vcan","Cspg4","Olig1","Cldn5","Rgs5","Vtn","Flt1","Slc17a7","Slc17a6","Satb2","Syt1","Snap25","Nrgn","Grin1","Fezf2","Bcl11b","Foxp2","Rorb")
downsample<- subset(merge, downsample = 300)
my_level<-c("Excitatory neurons","Inhibitory neurons", "Astrocytes","Microglia","Oligodendrocytes","OPC","Endothelial cells","Unknown")

downsample@active.ident <- factor(x = downsample@active.ident, levels = my_level)

DoHeatmap(downsample, features=markers, slot="scale.data", assay="integrated",size=5, angle= 40,group.colors=color)+ theme(axis.text=element_text(size=10, face="bold"))


###Differentially expressed genes
##Excitatory neurons
Excitatory_seurat_object<-subset(merge, idents= "Excitatory neurons")
Excitatory_seurat_object$celltype_subject<-paste(Idents(Excitatory_seurat_object), Excitatory_seurat_object$condition, sep= "_")

Excitatory_seurat_object@meta.data
Excitatory_seurat_object$celltype<-Idents(Excitatory_seurat_object)
Idents(Excitatory_seurat_object)<- "celltype_subject"
DE_excitatory_2m_all<- Excitatory_seurat_object%>%
  FindMarkers(ident.1= "Excitatory neurons_2m_KO", ident.2= "Excitatory neurons_2m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))


DE_excitatory_6m_all<- Excitatory_seurat_object%>%
  FindMarkers(ident.1= "Excitatory neurons_6m_KO", ident.2= "Excitatory neurons_6m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))


DE_excitatory_fad_all<- Excitatory_seurat_object%>%
  FindMarkers(ident.1= "Excitatory neurons_FAD-positive", ident.2= "Excitatory neurons_FAD-negative",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

##Inhibitory neurons
Inhibitory_seurat_object<-subset(merge, idents= "Inhibitory neurons")
Inhibitory_seurat_object$celltype_subject<-paste(Idents(Inhibitory_seurat_object), Inhibitory_seurat_object$condition, sep= "_")

Inhibitory_seurat_object@meta.data
Inhibitory_seurat_object$celltype<-Idents(Inhibitory_seurat_object)
Idents(Inhibitory_seurat_object)<- "celltype_subject"
DE_inhibitory_2m_all<- Inhibitory_seurat_object%>%
  FindMarkers(ident.1= "Inhibitory neurons_2m_KO", ident.2= "Inhibitory neurons_2m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_inhibitory_6m_all<- Inhibitory_seurat_object%>%
  FindMarkers(ident.1= "Inhibitory neurons_6m_KO", ident.2= "Inhibitory neurons_6m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_inhibitory_fad_all<- Inhibitory_seurat_object%>%
  FindMarkers(ident.1= "Inhibitory neurons_FAD-positive", ident.2= "Inhibitory neurons_FAD-negative",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

#Astrocytes
Astrocytes_seurat_object<-subset(merge, idents= "Astrocytes")
Astrocytes_seurat_object$celltype_subject<-paste(Idents(Astrocytes_seurat_object), Astrocytes_seurat_object$condition, sep= "_")

Astrocytes_seurat_object@meta.data
Astrocytes_seurat_object$celltype<-Idents(Astrocytes_seurat_object)
Idents(Astrocytes_seurat_object)<- "celltype_subject"

DE_Astrocytes_2m_all<- Astrocytes_seurat_object%>%
  FindMarkers(ident.1= "Astrocytes_2m_KO", ident.2= "Astrocytes_2m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Astrocytes_6m_all<- Astrocytes_seurat_object%>%
  FindMarkers(ident.1= "Astrocytes_6m_KO", ident.2= "Astrocytes_6m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Astrocytes_fad_all<- Astrocytes_seurat_object%>%
  FindMarkers(ident.1= "Astrocytes_FAD-positive", ident.2= "Astrocytes_FAD-negative",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

##Microglia
Microglia_seurat_object$celltype_subject<-paste(Idents(Microglia_seurat_object), Microglia_seurat_object$condition, sep= "_")

Microglia_seurat_object@meta.data
Microglia_seurat_object$celltype<-Idents(Microglia_seurat_object)
Idents(Microglia_seurat_object)<- "celltype_subject"

DE_Microglia_2m_all<- Microglia_seurat_object%>%
  FindMarkers(ident.1= "Microglia_2m_KO", ident.2= "Microglia_2m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Microglia_6m_all<- Microglia_seurat_object%>%
  FindMarkers(ident.1= "Microglia_6m_KO", ident.2= "Microglia_6m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Microglia_fad_all<- Microglia_seurat_object%>%
  FindMarkers(ident.1= "Microglia_FAD-positive", ident.2= "Microglia_FAD-negative",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

##Oligodendrocytes
Oligodendrocytes_seurat_object<-subset(merge, idents= "Oligodendrocytes")
Oligodendrocytes_seurat_object$celltype_subject<-paste(Idents(Oligodendrocytes_seurat_object), Oligodendrocytes_seurat_object$condition, sep= "_")

Oligodendrocytes_seurat_object@meta.data
Oligodendrocytes_seurat_object$celltype<-Idents(Oligodendrocytes_seurat_object)
Idents(Oligodendrocytes_seurat_object)<- "celltype_subject"

DE_Oligodendrocytes_2m_all<- Oligodendrocytes_seurat_object%>%
  FindMarkers(ident.1= "Oligodendrocytes_2m_KO", ident.2= "Oligodendrocytes_2m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Oligodendrocytes_6m_all<- Oligodendrocytes_seurat_object%>%
  FindMarkers(ident.1= "Oligodendrocytes_6m_KO", ident.2= "Oligodendrocytes_6m_control",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

DE_Oligodendrocytes_fad_all<- Oligodendrocytes_seurat_object%>%
  FindMarkers(ident.1= "Oligodendrocytes_FAD-positive", ident.2= "Oligodendrocytes_FAD-negative",logfc.threshold = 0)%>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

##Cell proportion 
table(merge@meta.data$condition)
prop.table(table(Idents(merge)))
table(Idents(merge), merge@meta.data$condition)
tab<-prop.table(table(Idents(merge), merge@meta.data$condition), margin = 2)
tab<-as.data.frame(tab)
names(tab)<-c("cell_type","condition","prop")
tab$prop<-signif(tab$prop,3)
tab
#ggplot(tab,aes(fill=cell_type,y= prop, x= condition)) + geom_bar(position="stack",stat = "identity", width= 0.5)+ ylab("Fraction of different cell types")+xlab("")+geom_text(aes(label = prop), size = 3, position = position_stack(vjust=0.5))

cell_prop<-ggplot(tab, aes(fill = cell_type, y = prop, x = "")) +
  geom_col(color = "white") +
  coord_polar(theta="y") +
  facet_wrap(~ condition, ncol = 2) +theme_void()+ scale_fill_manual(values=rev(color), limits = rev(order))+theme(text=element_text(size=25, face="bold"))

##Heatmap for Calcium, oxidative phosphorylation and synapse-related genes
#Excitatory neurons heatmap
calcium<-c("Atp2a2","Atp2b3","Atp1b1","Atp2b2","Atp2b4","Calm1","Calm2", "Calm3","Camk2a","Itpk1","Itpka","Slc8a2","Cacnb3","Cacna1a","Cacna1h","Cacng2","P2ry14","Ryr1","Gnb5","Gng3", "Gng2", "Gnb2", "Gnao1")

synapse_20<-c("Syp","Syt1","Syt2","Syt7","Syt13","Grin1","Grin2a","Grin2b","Grik3","Homer1","Slc17a7","Gria1","Gria2","Syn1","Syn2","Dlg4","Snap25","Sv2a","Sv2b","Vamp2","Rab3a")

etc_20<-c("Cox7c","Cox8a","Ndufa4","Cox4i1","Cox5a","Cox16","Atp5g1","Atp5a1","Ndufs7","Ndufb8","Cycs","Atp5b","Ndufc2","Atp5h","Uqcrh","Ndufb10","Cox6b1","Coq10b","Ndufv1","Ldha")

#E_sub<-subset(Inhibitory_seurat_object,subset= condition%in%"2m_control"| condition%in%"2m_KO" | condition%in%"FAD-negative"| condition%in%"FAD-positive")

avgexp_etc = AverageExpression(Inhibitory_seurat_object, assays = "RNA",slot = "data", group.by= c("condition"),features = etc_20)
avgexp_etc=as.matrix(avgexp_etc$RNA)

aa<-pheatmap(
  avgexp_etc,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",  # Set border_color to the desired color,
  fontsize = 20,
  scale="row")
ggsave("etc_I_20.pdf", ggplotify::as.ggplot(aa), dpi=300, width=12, height=30, unit="cm")


##Volcanoplot
top_up<-DE_Oligodendrocytes_6m%>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(n=10, wt= avg_log2FC) %>%
  pull(gene)
bottom_down<-DE_Oligodendrocytes_6m %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(n=-10, wt= avg_log2FC) %>%
  pull(gene)

volcanoplot<-EnhancedVolcano(DE_Oligodendrocytes_6m,
                             lab = DE_Oligodendrocytes_6m[,'gene'],
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             selectLab = c(top_up,bottom_down, "B2m", "H2-K1", "H2-D1"),
                             title = " ",
                             subtitle = bquote(italic( )),
                             axisLabSize = 20,
                             pCutoff = 0.05,
                             FCcutoff = 0.25,
                             pointSize = 1.5,
                             labSize = 5,
                             labCol = 'black',
                             labFace = 'bold',
                             boxedLabels = FALSE,
                             colAlpha = 1,
                             legendPosition = 'right',
                             legendLabSize = 15,
                             legendIconSize = 4.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             colConnectors = 'black',max.overlaps = 20,col = c('grey', 'forestgreen', 'royalblue', 'red2'))
