### Figure 6 ###

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(FSA)  

### load data objects ###
Lymphoma_data <- readRDS("./demo_data/Lymphoma_data.rds") ### This is saved from the step of Figure 1b
spatial_niche <- readRDS("./demo_data/spatial_niche.rds")
spatial_niche <- column_to_rownames(spatial_niche, var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, spatial_niche)
Lymphoma.meta <- Lymphoma_data@meta.data
Lymphoma.meta$Barcode <- rownames(Lymphoma.meta)
tumor_sample_comparison <- readRDS("./demo_data/tumor_sample_comparison.rds")
Lymphoma.meta <- left_join(Lymphoma.meta, tumor_sample_comparison, by="sample_ID")
my_neighbor_list <- readRDS("./demo_data/my_neighbor_list.rds") ### This is saved from the step of Figure 2d
cell_cell_communication_res <- readRDS("./demo_data/cell_cell_communication_res.rds")

### Figure 6a  compare cell state proportion ###
Tumor_site_comparison <- filter(Lymphoma.meta, Lymphoma.meta$Comparison2_tumor_anatomical_site %in% c("Nodal","IPS","EN-O"))
Tumor_site_comparison_FOV_meta <- unique(select(Tumor_site_comparison, c("FOV","Comparison2_tumor_anatomical_site")))

Tumor_site_FOV_table <- table(Tumor_site_comparison$FOV, Tumor_site_comparison$cell_state)
Tumor_site_FOV_table_prop <- as.data.frame(prop.table(Tumor_site_FOV_table, margin=1))
colnames(Tumor_site_FOV_table_prop)<-c("FOV","cell_state","Freq")
Tumor_site_FOV_table_prop <- left_join(Tumor_site_FOV_table_prop, Tumor_site_comparison_FOV_meta, by="FOV") # add back metadata

# compare C0_Tumor-B cell prop. in the FOVs from tumors in different sites
ggplot(filter(Tumor_site_FOV_table_prop, Tumor_site_FOV_table_prop$cell_state=="C0_Tumor-B"),
    aes(Comparison2_tumor_anatomical_site,Freq, fill=Comparison2_tumor_anatomical_site))+
    geom_boxplot(outliers=F)+
    geom_point(position="jitter", aes(color=Comparison2_tumor_anatomical_site))+
    ggtitle("Cell proportion of C0_Tumor-B in each FOV")+
    scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
    scale_color_manual(values=c("#fff7bc","#fe9929","#993404"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    theme(legend.position = "none", plot.title = element_text(size=10), 
          axis.title = element_text(size=8), 
          axis.text.x = element_text(angle = 45, hjust=1))+
    guides(fill=guide_legend(title = NULL))+
    xlab("Tumor site")+
    ylab("Prop. in all cells of an FOV")

# statistical test
kruskal.test(Freq ~ Comparison2_tumor_anatomical_site, data = filter(Tumor_site_FOV_table_prop, Tumor_site_FOV_table_prop$cell_state=="C0_Tumor-B"))
# For pairwaise comparisons, Dunn's post-hoc test was used. 
dunnTest(Freq ~ Comparison2_tumor_anatomical_site, data = filter(Tumor_site_FOV_table_prop, Tumor_site_FOV_table_prop$cell_state=="C0_Tumor-B"))

# similar for comparison of other cell prop. in the FOVs from nodal, IPS, and EN-O.

### Figure 6b  CN compositions tumors from different anatomical sites ###
Tumor_site_comparison_nodal <- filter(Tumor_site_comparison, Tumor_site_comparison$Comparison2_tumor_anatomical_site == "Nodal")
Tumor_site_comparison_nodal_table <- table(Tumor_site_comparison_nodal$CN_cluster)
Tumor_site_comparison_nodal_prop <- as.data.frame(prop.table(Tumor_site_comparison_nodal_table))

# show CN composition of nodal cases #
ggplot(Tumor_site_comparison_nodal_prop, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1)+
    coord_polar("y", start=0)+
    theme_void()+
    scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

# similar for IPS and EN-O cases

### Figure 6c-6d  comparison of cell state prop. in the neighborhood of C0_Tumor-B ###
C0_neighbor_log <- t(my_neighbor_list[["C0_Tumor-B"]])
C0_neighbor_log <- as.data.frame(log(C0_neighbor_log+1, base=2))
C0_neighbor_log <- rownames_to_column(C0_neighbor_log, var="Barcode")
C0_neighbor_log <- left_join(C0_neighbor_log, Lymphoma.meta, by="Barcode") # add metadata

C0_neighbor_log_comparison_tumor_site <- C0_neighbor_log[,c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                                            "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                                            "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                                            "C16_Stressed","C17_Epithelial","C18_RBC","Barcode","Comparison2_tumor_anatomical_site")]
C0_neighbor_log_comparison_tumor_site <- na.omit(C0_neighbor_log_comparison_tumor_site)

# compare the number of C0_Tumor-B in the neighborhood of C0 #
ggplot(C0_neighbor_log_comparison_tumor_site, aes(x=Comparison2_tumor_anatomical_site, y=`C0_Tumor-B`, fill=Comparison2_tumor_anatomical_site))+
    geom_boxplot(outlier.alpha = 0.5, color="#525252",width=0.5)+
    ggtitle("C0_Tumor-B in the neighborhood of C0_Tumor-B")+
    theme_classic()+
    scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    theme(legend.position = "none", plot.title = element_text(size=14), 
            axis.title = element_text(size=14), 
            axis.title.x = element_text(size=14),
            axis.text.x = element_text(angle = 45, hjust=1, size=14,),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=14))+
    guides(fill=guide_legend(title = NULL))+
    xlab("Tumor site")+
    ylab("log2(Number of neighboring cells + 1)")

# statistical test
kruskal.test(`C0_Tumor-B` ~ Comparison2_tumor_anatomical_site, data = C0_neighbor_log_comparison_tumor_site)
# For pairwaise comparisons, Dunn's post-hoc test was used.
dunnTest(`C0_Tumor-B` ~ Comparison2_tumor_anatomical_site, data = C0_neighbor_log_comparison_tumor_site)

# similar for comparison of the number of C5_T in the neighborhood of C0 #

### Figure 6e  CN allocation of T cells ###
T_cell_metadata <- filter(Lymphoma.meta, Lymphoma.meta$cell_state=="C5_T")
T_cell_niche_comparison_tumor_site <- table(T_cell_metadata$Comparison2_tumor_anatomical_site, T_cell_metadata$CN_cluster)
T_cell_niche_comparison_tumor_site_prop <- as.data.frame(prop.table(T_cell_niche_comparison_tumor_site, margin=1))
colnames(T_cell_niche_comparison_tumor_site_prop) <- c("Tumor_site","CN_cluster","Proportion")

# show CN allocations of T cells #
ggplot(T_cell_niche_comparison_tumor_site_prop,aes(Tumor_site,Proportion,fill=CN_cluster))+
    geom_bar(stat = "identity",position = "fill")+
    ggtitle("CN composition of T cells")+
    theme_bw()+
    scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    guides(fill=guide_legend(title = NULL))

### Figure 6f ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 6g-6h, 6l ###
cell_level_sample_meta_tumor_site <- select(Lymphoma.meta, c("Barcode","Comparison2_tumor_anatomical_site")) %>% column_to_rownames(var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, cell_level_sample_meta_tumor_site)

T_cell_comparison_tumor_site <- subset(Lymphoma_data, subset = cell_state=="C5_T" & Comparison2_tumor_anatomical_site %in% c("Nodal", "IPS", "EN-O"))
T_cell_comparison_tumor_site <- SCTransform(T_cell_comparison_tumor_site, assay = "RNA")

T_cytotoxicity_signature <- list(c("GZMK", "GPR183", "GNLY", "GZMH", "CX3CR1", "NKG7", "TBX21", "PRF1", "CD300A", "EOMES"))
T_exhaustion_signature <- list(c("GZMB", "HAVCR2", "ENTPD1", "VCAM1", "LAG3", "TNFRSF9", "CTLA4", "TIGIT", "TNFRSF18", "PDCD1", "IFNG", "CXCR6", "GNLY", "DUSP4"))
cell_cycle_signature <- list("G1S_sig" = c("PCNA","TYMS","HMGB2","DNMT1"),
                             "G2M_sig" = c("TOP2A","UBE2C","HMGB2","NUSAP1","CENPF","BIRC5","PTTG1","MKI67","CDKN3","STMN1"))

T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=T_cytotoxicity_signature, name = "T_cytotoxicity_signature",nbin=10,ctrl=20)
T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=T_exhaustion_signature, name = "T_exhaustion_signature",nbin=10,ctrl=20)
T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=cell_cycle_signature, name = "cell_cycle_signature",nbin=10,ctrl=20)

# plot #
ggplot(T_cell_comparison_tumor_site@meta.data, aes(x=Comparison2_tumor_anatomical_site, y=T_cytotoxicity_signature1, fill=Comparison2_tumor_anatomical_site))+
    geom_violin(scale = "width",)+
    geom_boxplot(width=0.15, outliers=F, fill="white")+
    scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
    theme_classic()+
    theme(legend.position = "none", plot.title = element_text(size=14), 
        axis.title = element_text(size=14), 
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust=1, size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14))+
    ggtitle('T_cytotoxicity_signature')+
    guides(fill=guide_legend(title = NULL))+
    xlab("Tumor site")+
    ylab("Signature score")

# statistical test
pairwise.t.test(T_cell_comparison_tumor_site@meta.data$T_cytotoxicity_signature1, T_cell_comparison_tumor_site@meta.data$Comparison2_tumor_anatomical_site, p.adjust.method = "BH")

# similar for plotting and comparing the G2/M signature score and exhuastion signature score

### Figure 6i ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 6j-6k was from CODEX images ###

### Figure 6m ###
Exhaustion_markers <- c("LAG3","HAVCR2","PDCD1","ENTPD1")

DotPlot(T_cell_comparison_tumor_site,features = Exhaustion_markers, group.by = "Comparison2_tumor_anatomical_site", scale=10, scale.by = "size")+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#253494","#4575b4","#abd9e9","#e0f3f8","#ffffbf","#d73027","#800026"))

### Figure 6n ###
# cell-cell communication was performed using Spyrrow in python, as detailed in descriptions for Figure 3g.
# define functions for calculating aggregrated communication cell pair prop. and communication intensity in different tumor sites
#a is the cluster id of sender cell, b is the cluster id of receiver cell, mol_a is the ligand from cell a, mol_b is the receptor from cell b
#in the original analysis, due to the high cell pair number involved in this step, we introduce a random step in the function to downsample the whole dataset to improve computational efficiency. Here since we are working on downsampled data, we omitted this step.

meta_for_cell1 <- select(Lymphoma.meta, c("cell_state","CN_cluster","sample_ID"))
colnames(meta_for_cell1)[1:2] <- c("cell_state1","CN_cluster1")
meta_for_cell1 <- rownames_to_column(meta_for_cell1, var="cell1")
cell_cell_communication_res <- left_join(cell_cell_communication_res, meta_for_cell1, by="cell1")

meta_for_cell2 <- select(Lymphoma.meta, c("cell_state","CN_cluster"))
colnames(meta_for_cell2)[1:2] <- c("cell_state2","CN_cluster2")
meta_for_cell2 <- rownames_to_column(meta_for_cell2, var="cell2")
cell_cell_communication_res <- left_join(cell_cell_communication_res, meta_for_cell2, by="cell2")

cell_cell_communication_res$cell_pair <- paste0(cell_cell_communication_res$cell1,":",cell_cell_communication_res$cell2)

tumor_sample_comparison_for_commun <- select(tumor_sample_comparison, c("sample_ID","Comparison2_tumor_anatomical_site"))
colnames(tumor_sample_comparison_for_commun)[2] <- "tumor_site"
cell_cell_communication_res <- left_join(cell_cell_communication_res, tumor_sample_comparison_for_commun, by="sample_ID")

aggretrate_tumor_site <- function(a,b,mol_a,mol_b){
    #cell_cell_communication_res_downsample <- cell_cell_communication_res[sample(nrow(cell_cell_communication_res), 500000), ]
    data_1 <- filter(cell_cell_communication_res, cell_cell_communication_res$cell_state1 %in% a & cell_cell_communication_res$cell_state2 %in% b) 
    data_1 <- data_1[,c(paste0(mol_a,".",mol_b,"_a2b"),"cell_pair","cell1","cell2","cell_state1","CN_cluster1","cell_state2","CN_cluster2","tumor_site")]
    data_2 <- filter(cell_cell_communication_res, cell_cell_communication_res$cell_state1 %in% b & cell_cell_communication_res$cell_state2 %in% a)
    data_2 <- data_2[,c(paste0(mol_a,".",mol_b,"_b2a"),"cell_pair","cell1","cell2","cell_state1","CN_cluster1","cell_state2","CN_cluster2","tumor_site")]
    colnames(data_1)[1] <- "intensity"
    colnames(data_2)[1] <- "intensity"
    data <- rbind(data_1, data_2)
    data_filter <- filter(data, !is.na(data$tumor_site))
    data_filter$LR <- paste0(mol_a,"-",mol_b)
    data_filter$count_ind <- 0
    data_filter$count_ind[data_filter$intensity > 0] <- 1
    group_mean <- data_filter %>% group_by(tumor_site) %>% summarise(mean_value = mean(intensity))
    group_mean <- as.numeric(group_mean$mean_value)
    group_count <- data_filter %>% group_by(tumor_site) %>% summarise(count_value = sum(count_ind))
    group_count_sum <- table(data_filter$tumor_site)
    group_prop <- as.numeric(group_count$count_value/group_count_sum)
    data_filter$group_mean <- NA
    data_filter$group_prop <- NA
    group <- unique(data_filter$tumor_site)[order(unique(data_filter$tumor_site))]
    for (i in 1:length(group)) {
        data_filter$group_mean[data_filter$tumor_site== group[i]] <- group_mean[i]
        data_filter$group_count_sum[data_filter$tumor_site== group[i]] <- group_count_sum[i]
        data_filter$group_prop[data_filter$tumor_site== group[i]] <- group_prop[i]*100
    }
    
    return(data_filter)
}

cell_state_other <- c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                          "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU","C16_Stressed","C17_Epithelial","C18_RBC")

CD274_PDCD1_other_T_sum <- aggretrate_tumor_site(cell_state_other,"C5_T","CD274","PDCD1") 
CD274_PDCD1_other_T_sum_for_plot <- unique(select(CD274_PDCD1_other_T_sum, c("tumor_site","LR","group_mean","group_count_sum","group_prop")))

ggplot(CD274_PDCD1_other_T_sum_for_plot, aes(x = LR, y = tumor_site, color = group_mean, size = group_prop)) + 
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))+
    geom_point(stat = 'identity') + 
    xlab("") + ylab("") + 
    ggtitle("Interaction intensity") + 
    coord_flip()+
    theme_classic()
  