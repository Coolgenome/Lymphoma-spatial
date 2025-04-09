### Figure 5 ###

### load essential packages ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)


### load data objects ###
Lymphoma_data <- readRDS("./Lymphoma_data.rds")
spatial_niche <- readRDS("./spatial_niche.rds")
spatial_niche <- column_to_rownames(spatial_niche, var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, spatial_niche)
Lymphoma.meta <- Lymphoma_data@meta.data
Lymphoma.meta$Barcode <- rownames(Lymphoma.meta)
tumor_sample_comparison <- readRDS("./tumor_sample_comparison.rds")
Lymphoma.meta <- left_join(Lymphoma.meta, tumor_sample_comparison, by="sample_ID")
my_neighbor_list <- readRDS("./my_neighbor_list.rds") ### This is saved from the step of Figure 2d

### Figure 5a  compare cell state proportion ###
EBV_comparison <- filter(Lymphoma.meta, Lymphoma.meta$Comparison1_nodal_EBV_status %in% c("Negative","Positive"))
EBV_comparison_FOV_meta <- unique(select(EBV_comparison, c("FOV","Comparison1_nodal_EBV_status")))

EBV_FOV_table <- table(EBV_comparison$FOV, EBV_comparison$cell_state)
EBV_FOV_table_prop <- as.data.frame(prop.table(EBV_FOV_table, margin=1))
colnames(EBV_FOV_table_prop)<-c("FOV","cell_state","Freq")
EBV_FOV_table_prop <- left_join(EBV_FOV_table_prop, EBV_comparison_FOV_meta, by="FOV") # add back metadata

# compare C0_Tumor-B cell prop. in the FOVs from EBV+ and EBV- cases
ggplot(filter(EBV_FOV_table_prop, EBV_FOV_table_prop$cell_state=="C0_Tumor-B"),
    aes(Comparison1_nodal_EBV_status, Freq, fill=Comparison1_nodal_EBV_status))+
    geom_boxplot(outliers=F)+
    geom_point(position="jitter", aes(color=Comparison1_nodal_EBV_status))+
    ggtitle("Cell proportion of C0_Tumor-B in each FOV")+
    scale_fill_manual(values=c("#a6cee3","#fdbf6f"))+
    scale_color_manual(values=c("#a6cee3","#fdbf6f"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    theme(legend.position = "none", plot.title = element_text(size=10), 
          axis.title = element_text(size=8), 
          axis.text.x = element_text(angle = 45, hjust=1))+
    guides(fill=guide_legend(title = NULL))+
    xlab("EBV status")+
    ylab("Prop. in all cells of an FOV")

# similar for comparison of C5_T cell prop. in the FOVs from EBV+ and EBV- cases

### Figure 5b  CN compositions in EBV+/- tumors ###
EBV_comparison_pos <- filter(EBV_comparison, EBV_comparison$Comparison1_nodal_EBV_status == "Positive")
EBV_comparison_pos_table <- table(EBV_comparison_pos$CN_cluster)
EBV_comparison_pos_prop <- as.data.frame(prop.table(EBV_comparison_pos_table))

# show CN composition of EBV positive cases #
ggplot(EBV_comparison_pos_prop, aes(x="", y=Freq, fill=Var1))+
    geom_bar(stat="identity", width=1)+
    coord_polar("y", start=0)+
    theme_void()+
    scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

# similar for EBV negative cases

### Figure 5c  CN allocation of T cells ###
T_cell_metadata <- filter(Lymphoma.meta, Lymphoma.meta$cell_state=="C5_T")
T_cell_niche_comparison_EBV <- table(T_cell_metadata$Comparison1_nodal_EBV_status, T_cell_metadata$CN_cluster)
T_cell_niche_comparison_EBV_prop <- as.data.frame(prop.table(T_cell_niche_comparison_EBV, margin=1))
colnames(T_cell_niche_comparison_EBV_prop) <- c("EBV_status","CN_cluster","Proportion")

# show CN allocation of T cells in EBV positive and negative cases #
ggplot(T_cell_niche_comparison_EBV_prop, aes(EBV_status,Proportion,fill=CN_cluster))+
    geom_bar(stat = "identity",position = "fill")+
    ggtitle("CN composition of T cells")+
    theme_classic()+
    scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    guides(fill=guide_legend(title = NULL))

### Figure 5d ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 5e  obtian the neighborhood composition matrix of C0_Tumor-B cells ###
C0_neighbor <- as.data.frame(t(my_neighbor_list[["C0_Tumor-B"]]))
C0_neighbor <- rownames_to_column(C0_neighbor, var="Barcode")
C0_neighbor <- left_join(C0_neighbor, Lymphoma.meta, by="Barcode")

C0_neighbor_comparison_EBV <- C0_neighbor[,c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                             "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                             "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                             "C16_Stressed","C17_Epithelial","C18_RBC","Barcode","Comparison1_nodal_EBV_status")]
C0_neighbor_comparison_EBV <- na.omit(C0_neighbor_comparison_EBV)

# calculate neighborhood compositions of C0_Tumor-B cells in EBV negative cases
C0_neighbor_comparison_EBV_neg <- filter(C0_neighbor_comparison_EBV, C0_neighbor_comparison_EBV$Comparison1_nodal_EBV_status=="Negative")
C0_neighbor_comparison_EBV_neg$Comparison1_nodal_EBV_status <- NULL
C0_neighbor_comparison_EBV_neg <- column_to_rownames(C0_neighbor_comparison_EBV_neg, var="Barcode")

C0_neighbor_comparison_EBV_neg_sum <- apply(C0_neighbor_comparison_EBV_neg, 2, sum)
C0_neighbor_comparison_EBV_neg_sum <- data.frame(cell_state = names(C0_neighbor_comparison_EBV_neg_sum),
                                                 neighbor_cell_count = C0_neighbor_comparison_EBV_neg_sum)
C0_neighbor_comparison_EBV_neg_sum$cell_state <- factor(C0_neighbor_comparison_EBV_neg_sum$cell_state, 
                                                        levels = c("C0_Tumor-B","C1_PC_IgG","C2_PC_IgA","C3_Resting-B","C4_PC_IgM","C5_T",
                                                                    "C6_TAM_APOE_C1Q","C7_TAM_SPP1","C8_Mac_DUSP1","C9_Mac_CXCL8","C10_Mac_MT2A",
                                                                    "C11_FRC","C12_HEV","C13_Endothelial_VWF","C14_VSMC","C15_Stromal_CLU",
                                                                    "C16_Stressed","C17_Epithelial","C18_RBC"))
C0_neighbor_comparison_EBV_neg_sum$prop <- C0_neighbor_comparison_EBV_neg_sum$neighbor_cell_count/sum(C0_neighbor_comparison_EBV_neg_sum$neighbor_cell_count)

# plot #
ggplot(C0_neighbor_comparison_EBV_neg_sum, aes(x="", y=prop, fill=cell_state)) +
    geom_bar(stat="identity", width=1)+
    coord_polar("y", start=0)+
    theme_void()+
    scale_fill_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF",
                                "#fff7fb","#fccde5","#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a",
                                "#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))

# similar for neighborhood compositions of C0_Tumor-B cells in EBV positive cases

### Figure 5f-5i ###
cell_level_sample_meta_EBV <- select(Lymphoma.meta, c("Barcode","Comparison1_nodal_EBV_status")) %>% column_to_rownames(var="Barcode")
Lymphoma_data <- AddMetaData(Lymphoma_data, cell_level_sample_meta_EBV)

T_cell_comparison_EBV <- subset(Lymphoma_data, cell_state=="C5_T" & Comparison1_nodal_EBV_status %in% c("Negative","Positive"))
T_cell_comparison_EBV <- SCTransform(T_cell_comparison_EBV, assay = "RNA")

T_cell_comparison_EBV_CN1 <- subset(T_cell_comparison_EBV, subset = CN_cluster_new == "CN1")
T_cell_comparison_EBV_CN1 <- SCTransform(T_cell_comparison_EBV_CN1, assay = "RNA")

T_cytotoxicity_signature <- list(c("GZMK", "GPR183", "GNLY", "GZMH", "CX3CR1", "NKG7", "TBX21", "PRF1", "CD300A", "EOMES"))
T_exhaustion_signature <- list(c("GZMB", "HAVCR2", "ENTPD1", "VCAM1", "LAG3", "TNFRSF9", "CTLA4", "TIGIT", "TNFRSF18", "PDCD1", "IFNG", "CXCR6", "GNLY", "DUSP4"))

T_cell_comparison_EBV_CN1 <- AddModuleScore(T_cell_comparison_EBV_CN1, features=T_cytotoxicity_signature, name = "T_cytotoxicity_signature",nbin=10,ctrl=20)
T_cell_comparison_EBV_CN1 <- AddModuleScore(T_cell_comparison_EBV_CN1, features=T_exhaustion_signature, name = "T_exhaustion_signature",nbin=10,ctrl=20)

T_cell_comparison_EBV_CN1_metadata <- T_cell_comparison_EBV_CN1@meta.data
T_cell_comparison_EBV_CN1_metadata$barcode <- rownames(T_cell_comparison_EBV_CN1_metadata)

T_cell_comparison_EBV_CN1_exp <- T_cell_comparison_EBV_CN1@assays$SCT$scale.data[c("HAVCR2","GZMB"),]
T_cell_comparison_EBV_CN1_exp <- as.data.frame(t(T_cell_comparison_EBV_CN1_exp))
T_cell_comparison_EBV_CN1_exp$barcode <- rownames(T_cell_comparison_EBV_CN1_exp)
T_cell_comparison_EBV_CN1_exp <- left_join(T_cell_comparison_EBV_CN1_exp, T_cell_comparison_EBV_CN1_metadata, by="barcode") 

# define functions for half violin plot #
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                        draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                grp <- data[1, "group"]
                                newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                
                                if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                            1))
                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                aesthetics$alpha <- rep(1, nrow(quantiles))
                                both <- cbind(quantiles, aesthetics)
                                quantile_grob <- GeomPath$draw_panel(both, ...)
                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                }
                                else {
                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                }
                            }
                        )
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                            draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                            show.legend = NA, inherit.aes = TRUE) {
                                                                    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
                                                                    position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                                                                    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
                                                                    }

# plot #
ggplot(T_cell_comparison_EBV_CN1_exp, aes(x=Comparison1_nodal_EBV_status, y=T_cytotoxicity_signature1, fill=Comparison1_nodal_EBV_status))+
    geom_split_violin(scale = "width")+
    geom_boxplot(width=0.15, outliers=F, fill="white")+
    scale_fill_manual(values=c("#a6cee3","#fdbf6f"))+
    theme_classic()+
    theme(legend.position = "none", plot.title = element_text(size=14), 
        axis.title = element_text(size=14), 
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust=1, size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14))+
    ggtitle('T_cytotoxicity_signature')+
    guides(fill=guide_legend(title = NULL))+
    xlab("EBV_status")+
    ylab("Signature score")
# similar for plotting gene expression and exhuastion signatrue score

### Figure 5j ###
# FOV images from CosMx SMI were plotted using Spyrrow (https://github.com/liuyunho/Spyrrow) in python (version 3.10.5).

### Figure 5k was from CODEX images ###