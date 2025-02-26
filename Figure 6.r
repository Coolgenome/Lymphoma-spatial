### Figure 6 ###
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)

### Figure 6a ###
spatial_niche <- readRDS("spatial_niche.rds") ### refer to Figure 2 code
Tumor_site_comparison <- filter(spatial_niche, spatial_niche$Tumor_site %in% c("Nodal","IPS","EN-O"))

ggplot(filter(Tumor_site_comparison, Tumor_site_comparison$cell_state=="C0_Tumor-B"),
       aes(Tumor_site,Freq, fill=Tumor_site))+
       geom_boxplot()+
       ggtitle("Cell proportion of C0_Tumor-B in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("Tumor site")+
       ylab("Prop. in all cells of an FOV")

ggplot(filter(Tumor_site_comparison, Tumor_site_comparison$Cluster_by_cell_type=="C5_T"),
       aes(Tumor_site, Freq, fill=Tumor_site))+
       geom_boxplot()+
       ggtitle("Cell proportion of C5_T in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("Tumor site")+
       ylab("Prop. in all cells of an FOV")

ggplot(filter(Tumor_site_comparison, Tumor_site_comparison$Cluster_by_cell_type=="C6_TAM_APOE_C1Q"),
       aes(Tumor_site, Freq, fill=Tumor_site))+
       geom_boxplot()+
       ggtitle("Cell proportion of C6_TAM_APOE_C1Q in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("Tumor site")+
       ylab("Prop. in all cells of an FOV")

ggplot(filter(Tumor_site_comparison, Tumor_site_comparison$Cluster_by_cell_type=="C11_FRC"),
       aes(Tumor_site, Freq, fill=Tumor_site))+
       geom_boxplot()+
       ggtitle("Cell proportion of C11_FRC in each FOV")+
       theme_classic()+
       scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(legend.position = "none", plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text.x = element_text(angle = 45, hjust=1))+
       guides(fill=guide_legend(title = NULL))+
       xlab("Tumor site")+
       ylab("Prop. in all cells of an FOV")

### Figure 6b ###
Tumor_site_comparison_nodal <- filter(Tumor_site_comparison, Tumor_site_comparison$Tumor_site == "Nodal")
Tumor_site_comparison_nodal_table <- table(Tumor_site_comparison_nodal$spatial_niche)
Tumor_site_comparison_nodal_prop <- as.data.frame(prop.table(Tumor_site_comparison_nodal_table))

ggplot(Tumor_site_comparison_nodal_prop, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#ed2224","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

Tumor_site_comparison_IPS <- filter(Tumor_site_comparison, Tumor_site_comparison$Tumor_site == "IPS")
Tumor_site_comparison_IPS_table <- table(Tumor_site_comparison_IPS$spatial_niche)
Tumor_site_comparison_IPS_prop <- as.data.frame(prop.table(Tumor_site_comparison_IPS_table))

ggplot(Tumor_site_comparison_IPS_prop, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#ed2224","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

Tumor_site_comparison_EN_O <- filter(Tumor_site_comparison, Tumor_site_comparison$Tumor_site == "EN-O")
Tumor_site_comparison_EN_O_table <- table(Tumor_site_comparison_EN_O$spatial_niche)
Tumor_site_comparison_EN_O_prop <- as.data.frame(prop.table(Tumor_site_comparison_EN_O_table))
   
ggplot(Tumor_site_comparison_EN_O_prop, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values=c("#ed2224","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))

### Figure 6c-6d ###
Lymphoma.meta <- readRDS("./Lymphoma.meta.for_neighborhood.rds")
my_neighbor_list <- readRDS("./my_neighbor_list.rds")

C0_neighbor_log <- my_neighbor_list$C0
C0_neighbor_log <- log(C0_neighbor_log+1, base=2)
C0_neighbor_log <- as.data.frame(t(C0_neighbor_log))
C0_neighbor_log <- rownames_to_column(C0_neighbor_log, var="cell_barcode")
C0_neighbor_log <- left_join(C0_neighbor_log, Lymphoma.meta, by="cell_barcode")

C0_neighbor_log_comparison_tumor_site <- C0_neighbor_log[,c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18",
                                                            "C3","C4","C7","C13","C15","C6","C12","C17","cell_state",
                                                            "Tumor_site")]
C0_neighbor_log_comparison_tumor_site <- na.omit(C0_neighbor_log_comparison_tumor_site)

ggplot(C0_neighbor_log_comparison_tumor_site, aes(x=Tumor_site, y=C0, fill=Tumor_site))+
  geom_boxplot(outlier.alpha = 0.5, color="#525252",width=0.5)+
  ggtitle("C0_Tumor-B in the neighborhood of C0")+
  theme_classic()+
  scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
  theme(axis.ticks.length = unit(0.5,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=14, family="Arial"), 
        axis.title = element_text(size=14, family="Arial"), 
        axis.title.x = element_text(size=14, family="Arial"),
        axis.text.x = element_text(angle = 45, hjust=1, size=14, family="Arial"),
        axis.title.y = element_text(size=14, family="Arial"),
        axis.text.y = element_text(size=14, family="Arial"))+
  guides(fill=guide_legend(title = NULL))+
  xlab("Tumor site")+
  ylab("log2(Number of neighboring cells + 1)")

ggplot(C0_neighbor_log_comparison_tumor_site, aes(x=Tumor_site, y=C1, fill=Tumor_site))+
  geom_boxplot(outlier.alpha = 0.5, color="#525252",width=0.5)+
  ggtitle("C5_T in the neighborhood of C0")+
  theme_classic()+
  scale_fill_manual(values=c("#fff7bc","#fe9929","#993404"))+
  theme(axis.ticks.length = unit(0.5,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=14, family="Arial"), 
        axis.title = element_text(size=14, family="Arial"), 
        axis.title.x = element_text(size=14, family="Arial"),
        axis.text.x = element_text(angle = 45, hjust=1, size=14, family="Arial"),
        axis.title.y = element_text(size=14, family="Arial"),
        axis.text.y = element_text(size=14, family="Arial"))+
  guides(fill=guide_legend(title = NULL))+
  xlab("Tumor site")+
  ylab("log2(Number of neighboring cells + 1)")

### Figure 6e ###
T_cell_metadata <- filter(spatial_niche, spatial_niche$cell_state=="C5_T")
T_cell_niche_comparison_tumor_site <- table(T_cell_metadata$Tumor_site, T_cell_metadata$spatial_niche)
cell.prop1<-as.data.frame(prop.table(T_cell_niche_comparison_tumor_site))
colnames(cell.prop1)<-c("Tumor_site","spatial_niche","Proportion")

ggplot(cell.prop1,aes(Proportion,Tumor_site,fill=spatial_niche))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("CN composition of T cells")+
  theme_bw()+
  scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
  theme(axis.ticks.length = unit(0.5,'cm'))+
  guides(fill=guide_legend(title = NULL))

### Figure 6f ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 6g-6h, 6l ###
Lymphoma_data_T_cell <- readRDS("./Lymphoma_data_T_cell.rds")  # refer to Figure 3a

T_cell_comparison_tumor_site <- subset(Lymphoma_data_T_cell, Tumor_site %in% c("Nodal", "IPS", "EN-O"))
T_cell_comparison_tumor_site <- SCTransform(T_cell_comparison_tumor_site, assay = "RNA")

T_cytotoxicity_signature <- list(c("GZMK", "GPR183", "FGFBP2", "GNLY", "GZMH", "CX3CR1", "FCGR3A", "PLEK", "KLRD1", "HOPX", "KLRG1", 
                                   "NKG7", "ZEB2", "TBX21", "PRF1", "CD300A", "CD244", "EOMES", "GZMM", "SEMA4A"))
T_exhaustion_signature <- list(c("CXCL13", "GZMB", "HAVCR2", "CCL3", "ENTPD1", "VCAM1", "LAG3", "TNFRSF9", "CTLA4", "CCL4L1", 
                                 "TIGIT", "TNFRSF18", "PDCD1", "IFNG", "ACP5", "LAYN", "CXCR6", "GNLY", "PHLDA1", "DUSP4"))
cell_cycle_signature <- list("G1S_sig" = c("PCNA","TYMS","HMGB2","DNMT1"),
                             "G2M_sig" = c("TOP2A","UBE2C","HMGB2","NUSAP1","CENPF","BIRC5","PTTG1","MKI67","CDKN3","STMN1"))

T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=T_cytotoxicity_signature, name = "T_cytotoxicity_signature",nbin=10,ctrl=20)
T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=T_exhaustion_signature, name = "T_cytotoxic_signature",nbin=10,ctrl=20)
T_cell_comparison_tumor_site <- AddModuleScore(T_cell_comparison_tumor_site, features=cell_cycle_signature, name = "cell_cycle_signature",nbin=10,ctrl=20)

# plot #
ggplot(T_cell_comparison_tumor_site@meta.data, aes(x=Tumor_site, y=T_cytotoxicity_signature1, fill=Tumor_site))+
    geom_violin(scale = "width")+
    geom_boxplot(width=0.15, outlier.size = 0.5, outlier.alpha = 0, outlier.color = "#f0f0f0")+
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

# similar for plotting the G2/M signature score and exhuastion signature score

### Figure 6i ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 6j-6k was from CODEX images ###

### Figure 6m ###
Exhaustion_markers <- c("LAG3","HAVCR2","PDCD1","ENTPD1")

DotPlot(T_cell_comparison_tumor_site,features = Exhaustion_markers, group.by = "Tumor_site", scale=10, scale.by = "size")+
  RotatedAxis()+
  scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))
  
### Figure 6n ###
# cell-cell communication was performed using Spyrrow in python, see the document "Cell cell communication" for detailed information.
cell_cell_communication_res <- readRDS("cell_cell_communication_res.rds")     # Commun_merge

# define functions for calculating aggregrated communication cell pair prop. and communication intensity in different tumor sites
#a is the cluster id of sender cell, b is the cluster id of receiver cell, mol_a is the ligand from cell a, mol_b is the receptor from cell b
aggretrate_sum_tumor_site <- function(a,b,mol_a,mol_b){
      data_1 <- filter(Commun_merge, Commun_merge$cell_type1 %in% a & Commun_merge$cell_type2 %in% b) 
      data_1 <- data_1[,c(paste0(mol_a,".",mol_b,"_a2b"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide","tumor_site")]
      data_2 <- filter(Commun_merge, Commun_merge$cell_type1 %in% b & Commun_merge$cell_type2 %in% a)
      data_2 <- data_2[,c(paste0(mol_a,".",mol_b,"_b2a"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide","tumor_site")]
      colnames(data_1)[1] <- "intensity"
      colnames(data_2)[1] <- "intensity"
      data <- rbind(data_1, data_2)
      data_filter <- filter(data, !is.na(data$tumor_site))
      data_filter$LR <- paste0(mol_a,"-",mol_b)
      data_filter$count_ind <- 0
      data_filter$count_ind[data_filter$intensity > 0] <- 1
      group_sum <- data_filter %>% group_by(tumor_site) %>% summarise(sum = sum(intensity))
      group_sum <- as.numeric(group_sum$sum)
      print("group_sum:"); print(group_sum)
      group_count <- data_filter %>% group_by(tumor_site) %>% summarise(count_value = sum(count_ind))
      group_count_sum <- table(data_filter$tumor_site)
      group_prop <- as.numeric(group_count$count_value/group_count_sum)
      print("group_prop:"); print(group_prop)
      data_filter$group_sum <- NA
      data_filter$group_prop <- NA
      group <- unique(data_filter$tumor_site)[order(unique(data_filter$tumor_site))]
      for (i in 1:length(group)) {
        data_filter$group_sum[data_filter$tumor_site== group[i]] <- group_sum[i]
        data_filter$group_count_sum[data_filter$tumor_site== group[i]] <- group_count_sum[i]
        data_filter$group_prop[data_filter$tumor_site== group[i]] <- group_prop[i]*100
      }
      
      return(data_filter)
}

CD274_PDCD1_other_T_sum <- aggretrate_sum_tumor_site(c(0,2:18),1,"CD274","PDCD1") 
CD274_PDCD1_other_T_sum_for_plot <- unique(select(CD274_PDCD1_other_T_sum, c("tumor_site","LR","group_sum","group_count_sum","group_prop")))
    
ggplot(CD274_PDCD1_T_filter_for_plot, aes(x = LR, y = tumor_site, color = group_sum, size = group_prop)) + 
  scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))+
  geom_point(stat = 'identity') + 
  xlab("") + ylab("") + 
  ggtitle("Interaction intensity") + 
  coord_flip()+
  theme_classic()