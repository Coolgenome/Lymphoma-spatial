### Data reading in, preprocessing, cleaning, and cell type identification

### load essential packages ###
library(Seurat)
library(Harmony)
library(dplyr)

### read in raw data and construct Seurat object ###
# For TMA slide F1
# read in data
F1_count <- read.csv("F1_6_exprMat_file.csv")
F1_count <- filter(F1_count, F1_count$cell_ID != 0) # All transcripts not assigned to a cell are show with a “cell_ID” value of 0. These are filtered out.
rownames(F1_count) <- paste0("F1_",F1_count$cell_ID,"_",F1_count$fov)
F1_count <- as.matrix(t(F1_count[, 3:1002])) # get count matrix with gene in the row names and cell ID in the column names
F1_meta <- read.csv("F1_6_metadata_file.csv")
rownames(F1_meta) <- paste0("F1_",F1_meta$cell_ID,"_",F1_meta$fov)
F1_meta$Slide <- "F1"

# For TMA slide F2
# read in data
F2_count <- read.csv("F2_1_exprMat_file.csv")
F2_count <- filter(F2_count, F2_count$cell_ID != 0)
rownames(F2_count) <- paste0("F2_",F2_count$cell)
F2_count <- as.matrix(t(F2_count[, 4:1003])) # get count matrix with gene in the row names and cell ID in the column names
F2_meta <- read.csv("F2_1_metadata_file.csv")
rownames(F2_meta) <- paste0("F2_",F2_meta$cell)
F2_meta$Slide <- "F2"

# It should be noted that due to the version differences of TMA slides F1,5,6 and slides F2,3,4, there are slight differences in the gene panels and meta data terms.
# Here we take the common genes and common meta data terms from both datasets
common_genes <- intersect(rownames(F1_count), rownames(F2_count))
common_terms <- intersect(colnames(F1_meta), colnames(F2_meta))
F1_count <- F1_count[common_genes,]; F1_meta <- F1_meta[,common_terms]
F2_count <- F2_count[common_genes,]; F2_meta <- F2_meta[,common_terms]

# Create Seurat objects
F1_obj <- CreateSeuratObject(counts = F1_count)
F1_obj <- AddMetaData(F1_obj, F1_meta)
F2_obj <- CreateSeuratObject(counts = F2_count)
F2_obj <- AddMetaData(F2_obj, F2_meta)
#The Seurat objects for other TMA slides are constructed in a similar way.

#Merge all objects and initial quality filtering
Lymphoma_data <- merge(F1_obj, y = c(F2_obj, F3_obj, F4_obj, F5_obj, F6_obj))
Lymphoma_data <- JoinLayers(Lymphoma_data)
Lymphoma_data <- subset(Lymphoma_data, nCount_RNA >= 20 & nFeature_RNA >= 10) # Cells with low nCount and low nFeature were filtered out.

### Data normalization, dimension reduction, unsupervised clustering following standardized steps in Seurat. ###
Lymphoma_data <- SCTransform(Lymphoma_data, assay = "RNA") #use scTransform for data normalization
DefaultAssay(Lymphoma_data) <- "SCT"
Lymphoma_data <- RunPCA(Lymphoma_data,features=VariableFeatures(object=Lymphoma_data)) #linear dimension reduction with PCA
ElbowPlot(Lymphoma_data, ndims = 50, reduction = "pca") # to determine the number of dimensions used for downstream analysis, here we select the top 40 PCs.
Lymphoma_data <- RunHarmony(Lymphoma_data, "Slide", plot_convergence=T, max.iter.harmony = 50) #remove TMA slide-specific batch effects, integrate data in the PCA space.
Lymphoma_data <- FindNeighbors(Lymphoma_data, reduction="harmony", dims=1:40)
Lymphoma_data <- FindClusters(Lymphoma_data,resolution=0.1) 
Lymphoma_data <- RunUMAP(Lymphoma_data, reduction = "harmony", dims = 1:40)

# check the UMAP embeddings and clustering of cells
DimPlot(Lymphoma_data, reduction = "umap",label = T, raster = T)

# check the expression of key lineage markers in each cluster
T_lineage_markers <- c("CD3D","CD3E","CD2","CD4","CD40LG","CD8A","CD8B","GZMK","GZMB","GZMH","NKG7")
B_lineage_markers <- c("MS4A1","CD19","CD79A","TCL1A","CD38","JCHAIN","MZB1","XBP1","IGHG1","IGHG2","IGHA1","IGHM")
Myeloid_lineage_markers <- c("CD14","CD68","CD163","LYZ","S100A8","S100A9","APOE","C1QB","C1QC","HLA-DQA1")
Stromal_lineage_markers <- c("VIM","COL1A1","DCN","TIMP1","FN1","ACTA2","VWF","PECAM1","MGP","CLU","VCAM1")
Epithelial_lineage_markers <- c("EPCAM","CDH1","KRT7","KRT8","KRT6A/B/C")
RBC_markers <- c("HBB","HBA1/2")

# To visualization gene expressions, we here use DotPlot function.
DotPlot(Lymphoma_data, features = T_lineage_markers, scale.by = "size")+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))
# Similarly we can visualize the expression of other lineage markers.

# Identify DEGs expressed in each cluster
Lymphoma_markers <- FindAllMarkers(Lymphoma_data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
Lymphoma_markers_top30 <- Lymphoma_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n=30)

# By examining the expression of lineage markers and the DEGs in each cluster, we can identify some cell clusters with the expression of multiple markers from different cell lineages,
# and also some clusters with very low expression level of any lineage markers. These findings may be due to multiple technical reasons, such as inaccurate cell segmentations. Since
# these cells without clear lineage identities are not helpful for downstream analysis, we will remove these clusters from the object and re-run the above steps (data normalization,
# dimension reduction, clustering). Then, we will re-examine the expression of lineage markers and DEGs of each new cluster, and determine if any additional clusters need to be removed.
# In such a way, multiple rounds of filtering and clustering were performed until we get clear clustering results, based on which major cell types and detailed cell states were defined.