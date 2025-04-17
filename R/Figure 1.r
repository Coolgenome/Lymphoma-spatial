### Figure 1 ###
### The codes are separated by Figures, and a detailed description of the analysis flow are provided to ensure readers' understanding and the transparency and reproducibility of the results.

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

### Figure 1a was created with Biorender ###

### Figure 1b ###
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

### Here, for convenience, we provided the demo data count matrix and detailed metadata ###
Lymphoma_count <- readRDS("./demo_data/Lymphoma_count.rds")
Lymphoma.meta <- readRDS("./demo_data/Lymphoma.meta.rds")

### Explanation of metadata ###
# CenterX_local_px, CenterY_local_px: x and y coordinates in each FOV
# CenterX_global_px, CenterY_global_px: x and y coordinates in the whole image
# UMAP1, UMAP2: UMAP embeddings
# Slide: TMA slide ID
# FOV: FOV ID
# sampl_ID: sample ID (detailed sample information can be found in Table S1)
# major_lineage: major cell types identified in the dataset, as shown in Figure 1b.
# cell_state: detailed cell states identified in the dataset, as shown in Figure 2a.

Lymphoma_data <- CreateSeuratObject(counts=Lymphoma_count)
Lymphoma_data <- AddMetaData(Lymphoma_data, Lymphoma.meta)
Lymphoma_data <- SCTransform(Lymphoma_data, assay = "RNA")

Lymphoma_data <- saveRDS("./demo_data/Lymphoma_data.rds")


### UMAP for major cell types ###
ggplot(Lymphoma.meta, aes(x=UMAP1, y=UMAP2, color=major_lineage))+
      geom_point(size=0.1, alpha=0.3)+
      scale_color_manual(values=c("#D10000","#96F148","#0000FF","#ffed6f","#fb9a99","#006837","#6a3d9a"))+
      theme_classic()


### Figure 1c ###
major_cell_type_marker <- c("CD3D","CD3E","CD3G","CD2","NKG7","CCL5","MS4A1","CD79A","TCL1A","CD37","MZB1",
                            "JCHAIN","CD68","CD163","APOC1","APOE","C1QB","C1QC","COL1A1","COL1A2","COL3A1",
                            "DCN","FN1","TIMP1","ACTA2","HSPA1A/B","HSPB1","KRT5","KRT6A/B/C","KRT17","KRT16",
                            "HBB","HBA1/2")

DotPlot(Lymphoma_data,features = major_cell_type_marker, scale.by = "size", group.by = "major_lineage", scale.max = 70, scale = 10)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 1d and Extended Data Figure 1 ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

# Barplot showing major lineage compositions for each sample in CosMx and CODEX datasets (Figure 1d bottom)
Major_lineage_prop <- readRDS('./Major_lineage_prop.rds')
for (i in 1:nrow(Major_lineage_prop)) {
      sample_name <- Major_lineage_prop$sample_ID[i]
      sample <- data.frame(cell_type = c("T","B","Myeloid","Stromal","T", "B","Myeloid","Stromal"),
                           cell_prop = c(as.numeric(Major_lineage_prop[i,c(2:5)]), as.numeric(Major_lineage_prop[i,c(6:9)])),
                           modality = c(rep('CODEX',4), rep('CosMx',4)))
      sample$cell_type <- factor(sample$cell_type, levels=rev(c("T","B","Myeloid","Stromal")))
      
      ggplot(sample, aes(cell_prop, modality, fill=cell_type))+
            geom_bar(stat = "identity",position = "fill")+
            ggtitle(sample_name)+
            theme_classic()+
            scale_fill_manual(values=c("#ffed6f","#0000FF","#96F148","#D10000"))+
            theme(axis.ticks.length = unit(0.2,'cm'))+
            guides(fill=guide_legend(title = NULL))
      ggsave(paste0("./sample_level_cell_prop/",sample_name, ".pdf"), width=6,height=1.5)
}