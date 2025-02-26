### Figure 1 ###

### Figure 1a was created with Biorender ###

### Figure 1b ###
library(Seurat)

### data normalization, dimension reduction, unsupervised clustering following standardized steps in Seurat ###
Lymphoma_data <- readRDS("./Lymphoma_data.rds") 

DimPlot(Lymphoma_data, reduction = "umap",label = F, raster = T, group.by="major_lineage", 
        cols=c("#D10000","#96F148","#0000FF","#ffed6f","#fb9a99","#006837","#6a3d9a"))

### Figure 1c ###
major_cell_type_marker <- c("CD3D","CD3E","CD3G","CD2","NKG7","CCL5","MS4A1","CD79A","TCL1A","CD37","MZB1",
                            "JCHAIN","CD68","CD163","APOC1","APOE","C1QB","C1QC","COL1A1","COL1A2","COL3A1",
                            "DCN","FN1","TIMP1","ACTA2","HSPA1A/B","HSPB1","KRT5","KRT6A/B/C","KRT17","KRT16",
                            "HBB","HBA1/2")

DotPlot(Lymphoma_data,features = major_cell_type_marker, scale.by = "size", group.by = "major_lineage", scale.max = 70, scale = 10)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 1d ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.