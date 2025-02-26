### Figure 3 ###

### Figure 3a ###
library(Seurat)

Lymphoma_data <- readRDS("./Lymphoma_data.rds") 

Lymphoma_data_T_cell <- subset(Lymphoma_data, major_lineage == "T")
Lymphoma_data_T_cell <- SCTransform(Lymphoma_data_T_cell, assay = "RNA")

saveRDS(Lymphoma_data_T_cell,"Lymphoma_data_T_cell.rds")

T_cell_naive_mem_markers <- c("SELL","TCF7","CCR7","IL7R","ANXA1","GPR183")
T_cell_cytotoxicity_markers <- c("GZMK","GZMB","PRF1","NKG7","GMZH","GNLY","XCL1/2")
T_cell_exhaustion_markers <- c("LAG3","HAVCR2","TNFRSF9","TIGIT","PDCD1","ENTPD1","TOX")

DotPlot(Lymphoma_data_T_cell,features = T_cell_naive_mem_markers, group.by="spatial_niche", scale.by = "size", col.max = 1.5)+
        RotatedAxis()+
        scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

DotPlot(Lymphoma_data_T_cell,features = T_cell_cytotoxicity_markers, group.by="spatial_niche", scale.by = "size", col.max = 1.5)+
        RotatedAxis()+
        scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

DotPlot(Lymphoma_data_T_cell,features = T_cell_exhaustion_markers, group.by="spatial_niche", scale.by = "size", col.max = 1.5)+
        RotatedAxis()+
        scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 3b ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 3c-3d were generated in similar ways as Figure 3a ###

### Figure 3e was generated with BioRender ###

### Figure 3f ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 3g ###
# cell-cell communication was performed using Spyrrow in python.
cell_cell_communication_res <- readRDS("cell_cell_communication_res.rds")     # Commun_merge

# define functions for calculating aggregrated communication cell pair prop. and communication intensity
#a is the cluster id of sender cell, b is the cluster id of receiver cell, mol_a is the ligand from cell a, mol_b is the receptor from cell b
aggretrate_sum_sample_CN <- function(a,b,mol_a,mol_b){
      #a is the cluster# of sender cell, b is the cluster# of receiver cell
      data_1 <- filter(Commun_merge, Commun_merge$cell_type1 %in% a & Commun_merge$cell_type2 %in% b) 
      data_1 <- data_1[,c(paste0(mol_a,".",mol_b,"_a2b"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide","tumor_site","Specimen.name.2","specimen_CN1_CN2")]
      data_2 <- filter(Commun_merge, Commun_merge$cell_type1 %in% b & Commun_merge$cell_type2 %in% a)
      data_2 <- data_2[,c(paste0(mol_a,".",mol_b,"_b2a"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide","tumor_site","Specimen.name.2","specimen_CN1_CN2")]
      colnames(data_1)[1] <- "intensity"
      colnames(data_2)[1] <- "intensity"
      data <- rbind(data_1, data_2)
      data_filter <- filter(data, data$CN_cluster1==data$CN_cluster2)
      data_filter <- filter(data_filter, !is.na(data_filter$Specimen.name.2))
      data_filter$LR <- paste0(mol_a,"-",mol_b)
      data_filter$count_ind <- 0
      data_filter$count_ind[data_filter$intensity > 0] <- 1
      group_sum <- data_filter %>% group_by(specimen_CN1_CN2) %>% summarise(sum = sum(intensity))
      group_sum <- as.data.frame(group_sum)
      #print("group_sum:"); print(group_sum)
      group_count <- data_filter %>% group_by(specimen_CN1_CN2) %>% summarise(count_value = sum(count_ind))
      group_count <- as.data.frame(group_count)
      group_count_sum <- as.data.frame(table(data_filter$specimen_CN1_CN2))
      colnames(group_count_sum) <- c("specimen_CN1_CN2","group_count_sum")
      group_prop_df <- left_join(group_count_sum, group_count, by="specimen_CN1_CN2")
      group_prop_df$group_prop <- group_prop_df$count_value/group_prop_df$group_count_sum
      #print("group_prop:"); print(group_prop)
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

CD274_PDCD1_mye_T_sample <- aggretrate_sum_sample_CN(2,1,"CD274","PDCD1")
CD274_PDCD1_mye_T_sample_for_plot <- unique(select(CD274_PDCD1_mye_T_sample, c("CN_cluster1","Specimen.name.2","specimen_CN1_CN2","LR","group_sum","group_prop","group_count","group_count_sum")))

### calculate T cell count in each spatial niche ###
spatial_niche$specimen_CN1_CN2 <- paste0(spatial_niche$sample_ID,"_",spatial_niche$CN_cluster_new,"_",spatial_niche$CN_cluster_new)
sample_CN_cell_count_table <- as.data.frame(table(spatial_niche$specimen_CN1_CN2, spatial_niche$cell_state))
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


### Figure 3h ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.