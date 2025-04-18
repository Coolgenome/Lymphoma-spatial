### Figure 1 ###
### The codes are separated by Figures, and a detailed description of the analysis flow are provided to ensure readers' understanding and the transparency and reproducibility of the results.

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

### Figure 1a was created with Biorender ###

### Figure 1b ###
### Data reading in, preprocessing, cleaning, and cell type and state identification are described in the separate script Preprocessing.r
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
