### Figure 5 ###
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)

### Figure 5a ###
spatial_niche <- readRDS("spatial_niche.rds") ### refer to Figure 2 code
EBV_comparison <- filter(spatial_niche, spatial_niche$EBV_status %in% c("Negative","Positive"))

ggplot(filter(spatial_niche, spatial_niche$cell_state=="C0_Tumor-B"),
       aes(EBV_status, Freq, fill=EBV_status))+
       geom_boxplot()+
       ggtitle("Cell proportion of C0_Tumor-B in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#a6cee3","#fdbf6f"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("EBV status")+
       ylab("Prop. in all cells of an FOV")

ggplot(filter(spatial_niche, spatial_niche$cell_state=="C5_T"),
       aes(EBV_status, Freq, fill=EBV_status))+
       geom_boxplot()+
       ggtitle("Cell proportion of C5_T in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#a6cee3","#fdbf6f"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("EBV status")+
       ylab("Prop. in all cells of an FOV")
 
### Figure 5b ###
EBV_comparison_pos <- filter(EBV_comparison, EBV_comparison$EBV_status == "Positive")
EBV_comparison_pos_table <- table(EBV_comparison_pos$spatial_niche)
EBV_comparison_pos_prop <- as.data.frame(prop.table(EBV_comparison_pos_table))

ggplot(EBV_comparison_pos_prop, aes(x="", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

EBV_comparison_neg <- filter(spatial_niche, spatial_niche$EBV_status == "Negative")
EBV_comparison_neg_table <- table(EBV_comparison_neg$spatial_niche)
EBV_comparison_neg_prop <- as.data.frame(prop.table(EBV_comparison_neg_table))

ggplot(EBV_comparison_neg_prop, aes(x="", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

### Figure 5c ###
T_cell_metadata <- filter(spatial_niche, spatial_niche$cell_state=="C5_T")
T_cell_niche_comparison_EBV <- table(T_cell_metadata$EBV_status, T_cell_metadata$spatial_niche)
cell.prop1<-as.data.frame(prop.table(T_cell_niche_comparison_EBV))
colnames(cell.prop1)<-c("EBV_status","spatial_niche","Proportion")

ggplot(cell.prop1,aes(Proportion,EBV_status,fill=spatial_niche))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("CN composition of T cells")+
  theme_bw()+
  scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
  theme(axis.ticks.length = unit(0.5,'cm'))+
  guides(fill=guide_legend(title = NULL))

### Figure 5d ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 5e ###
Lymphoma.meta <- readRDS("./Lymphoma.meta.for_neighborhood.rds")
my_neighbor_list <- readRDS("./my_neighbor_list.rds")

C0_neighbor <- as.data.frame(t(my_neighbor_list$C0))
C0_neighbor <- rownames_to_column(C0_neighbor, var="cell_barcode")
C0_neighbor <- left_join(C0_neighbor, Lymphoma.meta, by="cell_barcode")

C0_neighbor_comparison2 <- C0_neighbor[,c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18",
                                          "C3","C4","C7","C13","C15","C6","C12","C17","cell_state",
                                          "EBV_status")]
C0_neighbor_comparison2 <- na.omit(C0_neighbor_comparison2)

# calculate neighborhood compositions of C0_Tumor-B cells in EBV negative cases
C0_neighbor_comparison2_EBV_neg <- filter(C0_neighbor_comparison2, C0_neighbor_comparison2$EBV_status=="Negative")
C0_neighbor_comparison2_EBV_neg$EBV_status <- NULL
C0_neighbor_comparison2_EBV_neg <- column_to_rownames(C0_neighbor_comparison2_EBV_neg, var="cell_barcode")

C0_neighbor_comparison2_EBV_neg_sum <- apply(C0_neighbor_comparison2_EBV_neg, 2, sum)
C0_neighbor_comparison2_EBV_neg_sum <- data.frame(cell_state = names(C0_neighbor_comparison2_EBV_neg_sum),
                                                  neighbor_cell_count = C0_neighbor_comparison2_EBV_neg_sum)
C0_neighbor_comparison2_EBV_neg_sum$cell_state <- factor(C0_neighbor_comparison2_EBV_neg_sum$cell_state, 
                                                         levels = c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18",
                                                                    "C3","C4","C7","C13","C15","C6","C12","C17"))
C0_neighbor_comparison2_EBV_neg_sum$prop <- C0_neighbor_comparison2_EBV_neg_sum$neighbor_cell_count/sum(C0_neighbor_comparison2_EBV_neg_sum$neighbor_cell_count)

ggplot(C0_neighbor_comparison2_EBV_neg_sum, aes(x="", y=prop, fill=cell_state)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF",
                             "#fff7fb","#fccde5","#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a",
                             "#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))

# calculate neighborhood compositions of C0_Tumor-B cells in EBV positive cases
C0_neighbor_comparison2_EBV_pos <- filter(C0_neighbor_comparison2, C0_neighbor_comparison2$EBV_status=="Positive")
C0_neighbor_comparison2_EBV_pos$EBV_status <- NULL
C0_neighbor_comparison2_EBV_pos <- column_to_rownames(C0_neighbor_comparison2_EBV_pos, var="cell_barcode")

C0_neighbor_comparison2_EBV_pos_sum <- apply(C0_neighbor_comparison2_EBV_pos, 2, sum)
C0_neighbor_comparison2_EBV_pos_sum <- data.frame(cell_state = names(C0_neighbor_comparison2_EBV_pos_sum),
                                                  neighbor_cell_count = C0_neighbor_comparison2_EBV_pos_sum)
C0_neighbor_comparison2_EBV_pos_sum$cell_state <- factor(C0_neighbor_comparison2_EBV_pos_sum$cell_state, 
                                                         levels = c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18",
                                                                    "C3","C4","C7","C13","C15","C6","C12","C17"))
C0_neighbor_comparison2_EBV_pos_sum$prop <- C0_neighbor_comparison2_EBV_pos_sum$neighbor_cell_count/sum(C0_neighbor_comparison2_EBV_pos_sum$neighbor_cell_count)

ggplot(C0_neighbor_comparison2_EBV_pos_sum, aes(x="", y=prop, fill=cell_state)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF",
                             "#fff7fb","#fccde5","#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a",
                             "#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))

### Figure 5f-5i ###
Lymphoma_data_T_cell <- readRDS("./Lymphoma_data_T_cell.rds")  # refer to Figure 3a

T_cell_comparison_EBV <- subset(Lymphoma_data_T_cell, EBV_status %in% c("Negative","Positive"))
T_cell_comparison_EBV <- SCTransform(T_cell_comparison_EBV, assay = "RNA")

T_cell_comparison_EBV_CN1 <- subset(T_cell_comparison_EBV, CN_cluster_new == "CN1")
T_cell_comparison_EBV_CN1 <- SCTransform(T_cell_comparison_EBV_CN1, assay = "RNA")

T_cytotoxicity_signature <- list(c("GZMK", "GPR183", "FGFBP2", "GNLY", "GZMH", "CX3CR1", "FCGR3A", "PLEK", "KLRD1", "HOPX", "KLRG1", 
                                   "NKG7", "ZEB2", "TBX21", "PRF1", "CD300A", "CD244", "EOMES", "GZMM", "SEMA4A"))
T_exhaustion_signature <- list(c("CXCL13", "GZMB", "HAVCR2", "CCL3", "ENTPD1", "VCAM1", "LAG3", "TNFRSF9", "CTLA4", "CCL4L1", 
                                 "TIGIT", "TNFRSF18", "PDCD1", "IFNG", "ACP5", "LAYN", "CXCR6", "GNLY", "PHLDA1", "DUSP4"))

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
ggplot(T_cell_comparison_EBV_CN1_exp, aes(x=EBV_status, y=T_cytotoxicity_signature1, fill=EBV_status))+
    geom_split_violin(scale = "width")+
    geom_boxplot(width=0.15, outlier.size = 0.5, outlier.alpha = 0, outlier.color = "#f0f0f0")+
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
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 5k was from CODEX images ###