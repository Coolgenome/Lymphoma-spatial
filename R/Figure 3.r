### Figure 3 ###

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)


### load data objects ###
Lymphoma_data <- readRDS("./demo_data/Lymphoma_data.rds") ### This is saved from Figure 1b.
spatial_niche <- readRDS("./demo_data/spatial_niche.rds")
spatial_niche <- column_to_rownames(spatial_niche, var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, spatial_niche)

### Figure 3a ###
### for visualization of the expression of naive and memory markers ###
T_cell_naive_mem_markers <- c("SELL","TCF7","CCR7","IL7R","ANXA1","GPR183")

DotPlot(subset(Lymphoma_data, cell_state == "C5_T"),features = T_cell_naive_mem_markers, group.by="CN_cluster", scale.by = "size", scale=10, col.max = 1.5)+
        RotatedAxis()+
        scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))
# similar for cytotoxic markers and exhaustion markers

### Figure 3b ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 3c-3d were generated in similar ways as Figure 3a ###

### Figure 3e was generated with BioRender ###

### Figure 3f ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 3g ###
# cell-cell communication was performed using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5). After running the program, a matrix was obtianed, with each row representing a cell pair and each column representing a ligand-receptor pair #
# We then add metadata for cell 1 and cell 2 in each interacting cell pair (each row), including cell state, CN cluster, and sample ID #
cell_cell_communication_res <- readRDS("./demo_data/cell_cell_communication_res.rds")
Lymphoma.meta <- Lymphoma_data@meta.data
meta_for_cell1 <- select(Lymphoma.meta, c("cell_state","CN_cluster","sample_ID"))
colnames(meta_for_cell1)[1:2] <- c("cell_state1","CN_cluster1")
meta_for_cell1 <- rownames_to_column(meta_for_cell1, var="cell1")
cell_cell_communication_res <- left_join(cell_cell_communication_res, meta_for_cell1, by="cell1")

meta_for_cell2 <- select(Lymphoma.meta, c("cell_state","CN_cluster"))
colnames(meta_for_cell2)[1:2] <- c("cell_state2","CN_cluster2")
meta_for_cell2 <- rownames_to_column(meta_for_cell2, var="cell2")
cell_cell_communication_res <- left_join(cell_cell_communication_res, meta_for_cell2, by="cell2")

cell_cell_communication_res$cell_pair <- paste0(cell_cell_communication_res$cell1,":",cell_cell_communication_res$cell2)

# define functions for calculating aggregrated communication cell pair prop. and communication intensity for each CN in each sample
#a is the cluster id of sender cell, b is the cluster id of receiver cell, mol_a is the ligand from cell a, mol_b is the receptor from cell b
aggretrate_sum_sample_CN <- function(a,b,mol_a,mol_b){
    #a is the cluster# of sender cell, b is the cluster# of receiver cell
    cell_cell_communication_res$specimen_CN1_CN2 <- paste0(cell_cell_communication_res$sample_ID,"_", cell_cell_communication_res$CN_cluster1, "_", cell_cell_communication_res$CN_cluster2)
    data_1 <- filter(cell_cell_communication_res, cell_cell_communication_res$cell_state1 %in% a & cell_cell_communication_res$cell_state2 %in% b) 
    data_1 <- data_1[,c(paste0(mol_a,".",mol_b,"_a2b"),"cell_pair","cell1","cell2","cell_state1","CN_cluster1","cell_state2","CN_cluster2","sample_ID","specimen_CN1_CN2")]
    data_2 <- filter(cell_cell_communication_res, cell_cell_communication_res$cell_state1 %in% b & cell_cell_communication_res$cell_state2 %in% a)
    data_2 <- data_2[,c(paste0(mol_a,".",mol_b,"_b2a"),"cell_pair","cell1","cell2","cell_state1","CN_cluster1","cell_state2","CN_cluster2","sample_ID","specimen_CN1_CN2")]
    colnames(data_1)[1] <- "intensity"
    colnames(data_2)[1] <- "intensity"
    data <- rbind(data_1, data_2)
    data_filter <- filter(data, data$CN_cluster1==data$CN_cluster2)
    data_filter <- filter(data_filter, !is.na(data_filter$sample_ID))
    data_filter$LR <- paste0(mol_a,"-",mol_b)
    data_filter$count_ind <- 0
    data_filter$count_ind[data_filter$intensity > 0] <- 1
    group_sum <- data_filter %>% group_by(specimen_CN1_CN2) %>% summarise(sum = sum(intensity))
    group_sum <- as.data.frame(group_sum)
    group_count <- data_filter %>% group_by(specimen_CN1_CN2) %>% summarise(count_value = sum(count_ind))
    group_count <- as.data.frame(group_count)
    group_count_sum <- as.data.frame(table(data_filter$specimen_CN1_CN2))
    colnames(group_count_sum) <- c("specimen_CN1_CN2","group_count_sum")
    group_prop_df <- left_join(group_count_sum, group_count, by="specimen_CN1_CN2")
    group_prop_df$group_prop <- group_prop_df$count_value/group_prop_df$group_count_sum
    data_filter$group_sum <- NA
    data_filter$group_prop <- NA
    group <- group_count$specimen_CN1_CN2
    for (i in 1:length(group)) {
        data_filter$group_sum[data_filter$specimen_CN1_CN2 == group[i]] <- group_sum$sum[which(group_sum$specimen_CN1_CN2 == group[i])]
        data_filter$group_count[data_filter$specimen_CN1_CN2 == group[i]] <- group_prop_df$count_value[which(group_prop_df$specimen_CN1_CN2 == group[i])]
        data_filter$group_count_sum[data_filter$specimen_CN1_CN2 == group[i]] <- group_prop_df$group_count_sum[which(group_prop_df$specimen_CN1_CN2 == group[i])]
        data_filter$group_prop[data_filter$specimen_CN1_CN2 == group[i]] <- group_prop_df$group_prop[which(group_prop_df$specimen_CN1_CN2 == group[i])]*100
    }
    
    return(data_filter)
}

# calculate the PD-L1:PD-1 interaction between C6_TAM_APOE_C1Q and C5_T
CD274_PDCD1_mye_T_sample <- aggretrate_sum_sample_CN("C6_TAM_APOE_C1Q","C5_T","CD274","PDCD1")
CD274_PDCD1_mye_T_sample_for_plot <- unique(select(CD274_PDCD1_mye_T_sample, c("CN_cluster1","sample_ID","specimen_CN1_CN2","LR", "group_sum","group_prop","group_count","group_count_sum")))

### calculate total T cell count in each spatial niche ###
Lymphoma.meta$specimen_CN1_CN2 <- paste0(Lymphoma.meta$sample_ID,"_",Lymphoma.meta$CN_cluster,"_",Lymphoma.meta$CN_cluster)
sample_CN_cell_count_table <- as.data.frame(table(Lymphoma.meta$specimen_CN1_CN2, Lymphoma.meta$cell_state))
colnames(sample_CN_cell_count_table) <- c("specimen_CN1_CN2","cell_state","Count")
sample_CN_cell_count_table_T_cell <- filter(sample_CN_cell_count_table, sample_CN_cell_count_table$cell_state=="C5_T")
colnames(sample_CN_cell_count_table_T_cell)[3] <- "Total_T_count"

CD274_PDCD1_mye_T_sample_for_plot <- left_join(CD274_PDCD1_mye_T_sample_for_plot, sample_CN_cell_count_table_T_cell, by="specimen_CN1_CN2")

### calculate normalized interaction proportion, which equals the count of interaction cell pairs involving T cells divided by total number of T cells in the corresponding sample
CD274_PDCD1_mye_T_sample_for_plot$interac_prop <- CD274_PDCD1_mye_T_sample_for_plot$group_count / CD274_PDCD1_mye_T_sample_for_plot$Total_T_count

### calculate normalized interaction intensity, which equals the intensity of interactions involving T cells divided by total number of T cells in the corresponding sample
CD274_PDCD1_mye_T_sample_for_plot$group_sum_by_cell_count <- CD274_PDCD1_mye_T_sample_for_plot$group_sum / CD274_PDCD1_mye_T_sample_for_plot$Total_T_count

### visualizing normalized interaction proportion in each spatial niche
ggplot(CD274_PDCD1_mye_T_sample_for_plot, aes(x=CN_cluster1, y=interac_prop))+
    geom_boxplot(color="#969696",fill="white", outlier.alpha = 0, width=0.8)+
    geom_point(aes(color=CN_cluster1), position="jitter",size=2, alpha=0.5)+
    scale_fill_manual(values = c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
    scale_color_manual(values = c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
    theme_classic()

# similar for calculating the PD-L1:PD-1 interaction between C0_Tumor-B and C5_T

### Figure 3h ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).