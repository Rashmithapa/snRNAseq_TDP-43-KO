# snRNAseq_TDP-43-KO
Title: Single-nucleus RNA Sequencing Reveals GABAergic Vulnerability and Reactive Gliosis Driven by Loss of TDP-43
Link to the data in GEO database (GEO: GSE309401): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309401

R version 4.2.3 (2023-03-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 14.5

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape_0.8.9           DOSE_3.22.1             enrichplot_1.16.2      
 [4] ggupset_0.4.0           tibble_3.2.1            org.Mm.eg.db_3.15.0    
 [7] biomaRt_2.52.0          clusterProfiler_4.4.4   ggplot2_3.5.1          
[10] purrr_1.0.2             ensembldb_2.20.2        AnnotationFilter_1.20.0
[13] GenomicFeatures_1.48.4  AnnotationDbi_1.58.0    Biobase_2.56.0         
[16] GenomicRanges_1.48.0    GenomeInfoDb_1.32.4     IRanges_2.30.1         
[19] S4Vectors_0.34.0        AnnotationHub_3.4.0     BiocFileCache_2.4.0    
[22] dbplyr_2.5.0            BiocGenerics_0.42.0     cowplot_1.1.3          
[25] patchwork_1.2.0         dplyr_1.1.4             SeuratObject_4.1.3     
[28] Seurat_4.3.0           

loaded via a namespace (and not attached):
  [1] utf8_1.2.4                    R.utils_2.12.3               
  [3] spatstat.explore_3.3-1        reticulate_1.38.0            
  [5] tidyselect_1.2.1              RSQLite_2.3.7                
  [7] htmlwidgets_1.6.4             grid_4.2.3                   
  [9] BiocParallel_1.30.4           Rtsne_0.17                   
 [11] scatterpie_0.2.3              munsell_0.5.1                
 [13] codetools_0.2-19              ica_1.0-3                    
 [15] future_1.33.2                 miniUI_0.1.1.1               
 [17] withr_3.0.0                   spatstat.random_3.3-1        
 [19] colorspace_2.1-0              GOSemSim_2.22.0              
 [21] progressr_0.14.0              filelock_1.0.3               
 [23] knitr_1.48                    rstudioapi_0.16.0            
 [25] ROCR_1.0-11                   tensor_1.5                   
 [27] listenv_0.9.1                 MatrixGenerics_1.8.1         
 [29] GenomeInfoDbData_1.2.8        polyclip_1.10-6              
 [31] farver_2.1.2                  bit64_4.0.5                  
 [33] downloader_0.4                treeio_1.20.2                
 [35] parallelly_1.37.1             vctrs_0.6.5                  
 [37] generics_0.1.3                xfun_0.46                    
 [39] R6_2.5.1                      graphlayouts_1.1.1           
 [41] spatstat.univar_3.0-0         bitops_1.0-7                 
 [43] spatstat.utils_3.0-5          cachem_1.1.0                 
 [45] fgsea_1.22.0                  gridGraphics_0.5-1           
 [47] DelayedArray_0.22.0           promises_1.3.0               
 [49] BiocIO_1.6.0                  scales_1.3.0                 
 [51] ggraph_2.1.0                  gtable_0.3.5                 
 [53] globals_0.16.3                goftest_1.2-3                
 [55] tidygraph_1.3.0               rlang_1.1.4                  
 [57] splines_4.2.3                 rtracklayer_1.56.1           
 [59] lazyeval_0.2.2                spatstat.geom_3.3-2          
 [61] BiocManager_1.30.23           yaml_2.3.9                   
 [63] reshape2_1.4.4                abind_1.4-5                  
 [65] httpuv_1.6.15                 qvalue_2.28.0                
 [67] tools_4.2.3                   ggplotify_0.1.2              
 [69] RColorBrewer_1.1-3            ggridges_0.5.6               
 [71] Rcpp_1.0.13                   plyr_1.8.9                   
 [73] progress_1.2.3                zlibbioc_1.42.0              
 [75] RCurl_1.98-1.13               prettyunits_1.2.0            
 [77] deldir_2.0-2                  viridis_0.6.5                
 [79] pbapply_1.7-2                 zoo_1.8-12                   
 [81] SummarizedExperiment_1.26.1   ggrepel_0.9.5                
 [83] cluster_2.1.4                 fs_1.6.4                     
 [85] magrittr_2.0.3                data.table_1.15.4            
 [87] scattermore_1.2               DO.db_2.9                    
 [89] lmtest_0.9-40                 RANN_2.6.1                   
 [91] ProtGenerics_1.28.0           fitdistrplus_1.2-1           
 [93] matrixStats_1.3.0             hms_1.1.3                    
 [95] mime_0.12                     evaluate_0.24.0              
 [97] xtable_1.8-4                  XML_3.99-0.16                
 [99] gridExtra_2.3                 compiler_4.2.3               
[101] shadowtext_0.1.4              KernSmooth_2.23-20           
[103] crayon_1.5.3                  R.oo_1.26.0                  
[105] htmltools_0.5.8.1             ggfun_0.1.5                  
[107] later_1.3.2                   aplot_0.2.3                  
[109] tidyr_1.3.1                   DBI_1.2.3                    
[111] tweenr_2.0.3                  MASS_7.3-58.2                
[113] rappdirs_0.3.3                Matrix_1.5-3                 
[115] cli_3.6.3                     R.methodsS3_1.8.2            
[117] parallel_4.2.3                igraph_2.0.3                 
[119] pkgconfig_2.0.3               GenomicAlignments_1.32.1     
[121] sp_2.1-4                      plotly_4.10.4                
[123] spatstat.sparse_3.1-0         xml2_1.3.6                   
[125] ggtree_3.4.4                  XVector_0.36.0               
[127] yulab.utils_0.1.4             stringr_1.5.1                
[129] digest_0.6.36                 sctransform_0.4.1            
[131] RcppAnnoy_0.0.22              spatstat.data_3.1-2          
[133] Biostrings_2.64.1             rmarkdown_2.27               
[135] leiden_0.4.3.1                fastmatch_1.1-4              
[137] tidytree_0.4.6                uwot_0.2.2                   
[139] restfulr_0.0.15               curl_5.2.1                   
[141] shiny_1.8.1.1                 Rsamtools_2.12.0             
[143] rjson_0.2.21                  lifecycle_1.0.4              
[145] nlme_3.1-162                  jsonlite_1.8.8               
[147] viridisLite_0.4.2             fansi_1.0.6                  
[149] pillar_1.9.0                  lattice_0.20-45              
[151] KEGGREST_1.36.3               fastmap_1.2.0                
[153] httr_1.4.7                    survival_3.5-3               
[155] GO.db_3.15.0                  interactiveDisplayBase_1.34.0
[157] glue_1.7.0                    png_0.1-8                    
[159] BiocVersion_3.15.2            bit_4.0.5                    
[161] ggforce_0.4.2                 stringi_1.8.4                
[163] blob_1.2.4                    memoise_2.0.1                
[165] ape_5.7-1                     irlba_2.3.5.1                
[167] future.apply_1.11.2
![image](https://github.com/user-attachments/assets/2c3bff27-b6d5-4d24-9259-80473e717131)

