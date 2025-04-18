### Figure 2 ###

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

### Data reading in, preprocessing, cleaning, and cell type and state identification are described in the separate script Preprocessing.r
### Here for demonstrating the workflow, we directly provide the demo data, including count matrix and metadata. The processing of demo data is described in Figure 1.r
### load data object ###
Lymphoma_data <- readRDS("./demo_data/Lymphoma_data.rds") ### This is saved from the step of Figure 1b.

### Figure 2a ###
### UMAP for cell states ###
ggplot(Lymphoma_data@meta.data, aes(x=UMAP1, y=UMAP2, color=cell_state))+
      geom_point(size=0.1, alpha=0.3)+
      scale_color_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5",
                                  "#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837",
                                  "#6a3d9a"))+
      theme_classic()

### Figure 2b ###
cell_state_marker  <-c("MS4A1","CD79A","TCL1A","CD37","IGHM","TUBB","PCNA","TYMS","STMN1","MZB1","XBP1","IGHG1","IGHG2","IGHA1","CD44",
                       "JCHAIN","CD40","CD38","TNFRSF17","CD3D","CD3E","CD3G","CD2","NKG7","CCL5","CD68","CD163","APOC1","APOE",
                       "C1QC","C1QB","C1QA","SPP1","MMP12","DUSP1","FOS","JUN","CXCL8","IL1B","CXCL1/2/3","LGALS1","COL1A1","COL1A2",
                       "COL3A1","DCN","FN1","TIMP1","ACTA2","CCL21","CCL19","VWF","PECAM1","MGP","CLU","VCAM1","HSPA1A/B","HSPB1",
                       "KRT5","KRT6A/B/C","KRT17","KRT16","HBB","HBA1/2")

DotPlot(Lymphoma_data,features = cell_state_marker, scale.by = "size", group.by = "cell_state", scale.max = 60, scale = 10)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#225ea8","#4575b4","#6baed6","#abd9e9","#e0f3f8","#ffffbf","#d73027","#800026"))

### Figure 2c was created with BioRender ###

### Figure 2d ###
### neighborhood analysis ###
# spatial neighborhood was calculated based on the distance between cell centroids, as detailed in the manuscript Method. Here we use 200px as the neighborhood searching radius. #
Lymphoma.meta <- Lymphoma_data@meta.data
major.ord = unique(Lymphoma.meta$cell_state)
my_neighbor_list <- list()
for (clstype in major.ord){
    my.neighbor = c()
    print (clstype)
    
    for(sam in unique(Lymphoma.meta$FOV)) {
    print (sam)
    tmp.meta = Lymphoma.meta[Lymphoma.meta$cell_state == clstype & Lymphoma.meta$FOV == sam, ]
    tmp.meta1 = Lymphoma.meta[Lymphoma.meta$FOV == sam, ]
    
    if(clstype %in% unique(tmp.meta1$cell_state)){

        tmp.mtx = matrix(0,nrow = length(major.ord),ncol = nrow(tmp.meta))
        rownames(tmp.mtx) = major.ord
        colnames(tmp.mtx) = rownames(tmp.meta)
        
        for(i in 1:nrow(tmp.meta)) {
        xloc = tmp.meta[i, 'CenterX_local_px']
        yloc = tmp.meta[i, 'CenterY_local_px']
        cutoff = 200
        idx1 = tmp.meta1[, 'CenterX_local_px'] <=  (xloc + cutoff) & tmp.meta1[,'CenterX_local_px'] >=  (xloc - cutoff)
        idx2 = tmp.meta1[, 'CenterY_local_px'] <=  (yloc + cutoff) & tmp.meta1[,'CenterY_local_px'] >=  (yloc - cutoff)
        
        square = tmp.meta1[idx1 & idx2,]
        dis = apply(square, 1, function(x) {sqrt((as.numeric(x['CenterX_local_px']) - as.numeric(xloc))^2 + (as.numeric(x['CenterY_local_px']) - as.numeric(yloc))^2)} )
        
        tmp.freq = table(tmp.meta1[match(names(dis[dis <= cutoff]), rownames(tmp.meta1)), 'cell_state'])
        tmp.mtx[,i] = tmp.freq[match(major.ord,names(tmp.freq))]
        }
        
        tmp.mtx[which(is.na(tmp.mtx))] = 0
        
        my.neighbor = cbind(my.neighbor, tmp.mtx)
    }
    }
    
    my_neighbor_list[[clstype]] <- my.neighbor
    
    write.csv(t(my.neighbor), paste0('./neighborhood/',clstype,'.neighbor.sub.csv'),quote = F)
}

saveRDS(my_neighbor_list,"./demo_data/my_neighbor_list.rds")

# Following the described steps, we obtianed the nerighborhood composition matrix for each cell state #
# merge neighborhood composition count matrix for all cells #
cell_neighbor_df <- c()
for(i in 1:length(my_neighbor_list)){
    df <- as.data.frame(t(my_neighbor_list[[i]]))
    df$center_cell_state <- names(my_neighbor_list)[i]
    cell_neighbor_df <- rbind(cell_neighbor_df, df)
}

# aggregrate neighborhood cell counts by center cell type #
cell_state <- unique(Lymphoma.meta$cell_state)
Neighbor_sum <- matrix(0, nrow=19, ncol=19)
rownames(Neighbor_sum) = cell_state
colnames(Neighbor_sum) = cell_state

for(i in cell_state){
    df_filter <- filter(cell_neighbor_df, cell_neighbor_df$center_cell_state == i)
    df_filter <- df_filter[,cell_state]
    df_filter_sum = colSums(df_filter)
    Neighbor_sum[i,] = df_filter_sum
}

Neighbor_sum <- as.data.frame(Neighbor_sum)
Neighbor_sum$center_cell_state <- rownames(Neighbor_sum)

# calculate cell proportion in each neighborhood #
cell.prop <- as.data.frame(prop.table(as.matrix(Neighbor_sum[,0:19]), margin=1))
cell.prop$center_cell_state <- rownames(cell.prop)

cell.prop <- reshape2::melt(cell.prop, id.vars=c("center_cell_state"),
                            measure.vars= cell_state,
                            variable.name = "neighbor_cell_state", value.name="neighbor_cell_prop")

cell.prop$center_cell_state <- factor(cell.prop$center_cell_state, levels=c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                                                            "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                                                            "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                                                            "C16_Stressed","C17_Epithelial","C18_RBC"))              
cell.prop$neighbor_cell_state <- factor(cell.prop$neighbor_cell_state, levels=c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                                                                "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                                                                "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                                                                "C16_Stressed","C17_Epithelial","C18_RBC"))              

# plot #
ggplot(cell.prop,aes(center_cell_state, neighbor_cell_prop, fill=neighbor_cell_state))+
    geom_bar(stat = "identity",position = "fill")+
    ggtitle("Neighboring cell proportion for each cell state")+
    theme_classic()+
    theme(axis.ticks.length = unit(0.5,'cm'))+
    theme(axis.text.x = element_text(angle=60,hjust = 1))+
    guides(fill=guide_legend(title = "Cell type"))+
    scale_fill_manual(values = c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5","#bc80bd",
                                 "#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))

### Figure 2e-g ###
# k-means clustering based on neighborhood matrix to obtain 7 unique spatial niches. Here we directly provide the results with CN allocation of cells.
spatial_niche <- readRDS("./demo_data/spatial_niche.rds")

# project CN clusters in PCA space (Figure 2e) #
cell_neighbor_pca <- prcomp(cell_neighbor_df[,1:19])
cell_neighbor_pca_coord <- as.data.frame(cell_neighbor_pca$x)
cell_neighbor_pca_coord <- rownames_to_column(cell_neighbor_pca_coord, var="Barcode")
cell_neighbor_pca_coord <- left_join(cell_neighbor_pca_coord, spatial_niche, by="Barcode")

library(plotly)
plot_ly(cell_neighbor_pca_coord, x = ~PC1, y = ~PC2, z = ~PC3, color = ~CN_cluster,
        colors=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"),
        alpha=0.7, sizes= c(50,50))

# plot cell composition in each CN (Figure 2f) #
Lymphoma.meta <- rownames_to_column(Lymphoma.meta, var="Barcode")
Lymphoma.meta <- left_join(Lymphoma.meta, spatial_niche, by="Barcode")

cellnum <- table(Lymphoma.meta$cell_state, Lymphoma.meta$CN_cluster)
cellnum
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("cell_state","CN_cluster","Proportion")
cell.prop$cell_state <- factor(cell.prop$cell_state, levels=c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                                              "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                                              "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                                              "C16_Stressed","C17_Epithelial","C18_RBC"))

ggplot(cell.prop,aes(CN_cluster,Proportion,fill=cell_state))+
    geom_bar(stat = "identity",position = "fill")+
    ggtitle("Cell composition in each spatial niche")+
    theme_bw()+
    scale_fill_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5","#bc80bd",
                               "#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    guides(fill=guide_legend(title = NULL))

# DEGs among spatial niches (Figure 2g) #
spatial_niche <- column_to_rownames(spatial_niche, var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, spatial_niche)

CN_marker <- c("CD3D","CD3E","CD2","CD8A","GZMB","GZMK","NKG7","PRF1","GZMH","CXCL9","CXCL10","CCL5","CCL2","MZB1","IGHG1",
               "IGHG2","IGKC","CD68","C1QA","C1QB","C1QC","LYZ","APOE","APOC1","S100A9","COL1A1","COL1A2","COL3A1","FN1","DCN",
               "ACTA2","TIMP1","COL6A1","COL6A2","CD79A","MS4A1","CD19","TCL1A","IGHM","CD37","STMN1","TYMS","TOP2A","BCL2",
               "BTK","SYK","CD24","CD52")

DotPlot(Lymphoma_data,features = CN_marker, scale.by = "size", group.by = "CN_cluster", scale.max = 50, scale = 10, col.max=1.5)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 2h ###
### FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).
