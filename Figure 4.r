### Figure 4 ###

### Figure 4a ###
library(Seurat)

Lymphoma_data <- readRDS("./Lymphoma_data.rds") 

Lymphoma_data_B_cell <- subset(Lymphoma_data, major_lineage == "B")
Lymphoma_data_B_cell <- SCTransform(Lymphoma_data_B_cell, assay = "RNA")

B_proliferation_marker <- c("TYMS","STMN1","TOP2A","PCNA","TUBB","BCL2")

DotPlot(Lymphoma_data_B_cell,features = B_proliferation_marker, group.by="spatial_niche", scale.by = "size",col.max = 1.5)+RotatedAxis()+
        scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

cell_cycle_signature <- list("G1S_sig" = c("PCNA","TYMS","HMGB2","DNMT1"),
                             "G2M_sig" = c("TOP2A","UBE2C","HMGB2","NUSAP1","CENPF","BIRC5","PTTG1","MKI67","CDKN3","STMN1"))

Lymphoma_data_B_cell <- AddModuleScore(Lymphoma_data_B_cell, features=cell_cycle_signature, name = "cell_cycle_signature",nbin=10,ctrl=20)

ggplot(Lymphoma_data_B_cell@meta.data, aes(x=spatial_niche, y=cell_cycle_signature2, fill=spatial_niche))+
       geom_violin(scale='width')+
       geom_boxplot(width=0.3, outlier.size = 0, outlier.alpha = 0, outlier.color = "#f0f0f0")+
       scale_fill_manual(values=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"))+
       theme(legend.position = "none", plot.title = element_text(size=14), 
             axis.title = element_text(size=14), 
             axis.title.x = element_text(size=14),
             axis.text.x = element_text(angle = 45, hjust=1, size=14),
             axis.title.y = element_text(size=14),
             axis.text.y = element_text(size=14))+
       theme_classic()+
       ggtitle('G2M signature')+
       xlab("spatial_niche")+
       ylab("Signature score")


### Figure 4b-4c ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.

### Figure 4d ###
# cell-cell communication was performed using Spyrrow in python, see the document "Cell cell communication" for detailed information.
cell_cell_communication_res <- readRDS("cell_cell_communication_res.rds")     # Commun_merge

# define functions for calculating aggregrated communication cell pair prop. and communication intensity
#a is the cluster id of sender cell, b is the cluster id of receiver cell, mol_a is the ligand from cell a, mol_b is the receptor from cell b
aggretrate <- function(a,b,mol_a,mol_b){
      #a is the cluster# of sender cell, b is the cluster# of receiver cell
      data_1 <- filter(Commun_merge, Commun_merge$cell_type1==a & Commun_merge$cell_type2 == b) 
      data_1 <- data_1[,c(paste0(mol_a,".",mol_b,"_a2b"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide")]
      data_2 <- filter(Commun_merge, Commun_merge$cell_type1==b & Commun_merge$cell_type2 == a)
      data_2 <- data_2[,c(paste0(mol_a,".",mol_b,"_b2a"),"cell_pair","cell1","cell2","cell_type1","CN_cluster1",
                          "cell_type2","CN_cluster2","Slide")]
      colnames(data_1)[1] <- "intensity"
      colnames(data_2)[1] <- "intensity"
      data <- rbind(data_1, data_2)
      data_filter <- filter(data, data$CN_cluster1==data$CN_cluster2)
      data_filter$LR <- paste0(mol_a,"-",mol_b)
      data_filter$count_ind <- 0
      data_filter$count_ind[data_filter$intensity > 0] <- 1
      group_mean <- data_filter %>% group_by(CN_cluster1) %>% summarise(mean_value = mean(intensity))
      group_mean <- as.numeric(group_mean$mean_value)
      print("group_mean:"); print(group_mean)
      group_count <- data_filter %>% group_by(CN_cluster1) %>% summarise(count_value = sum(count_ind))
      group_count_sum <- table(data_filter$CN_cluster1)
      group_prop <- as.numeric(group_count$count_value/group_count_sum)
      print("group_prop:"); print(group_prop)
      data_filter$group_mean <- NA
      data_filter$group_prop <- NA
      group <- unique(data_filter$CN_cluster1)[order(unique(data_filter$CN_cluster1))]
      for (i in 1:length(group)) {
        data_filter$group_mean[data_filter$CN_cluster1== group[i]] <- group_mean[i]
        data_filter$group_prop[data_filter$CN_cluster1== group[i]] <- group_prop[i]*100
      }
      
      return(data_filter)
}
# aggregrate results for multiple sender cell types/ligands to one receiver cell type/receptor
# b is the receiver cell type, mol_b is the receptor on b, mol_a is the vector of ligand, this function calculates signals from all senders
agg_table <- function(b, mol_b, mol_a){
      table <- data.frame()
      for(i in mol_a){
        for(j in 0:18){
          df <- aggretrate(j,b,i,mol_b)
          df_unique <- unique(select(df, c("CN_cluster1","LR","group_mean","group_prop"))) 
          df_unique$sender <- paste0("C",j)
          table <- rbind(table, df_unique)
        }
      }
      return(table)
}

C0_CXCR4_agg_table <- agg_table(0, "CXCR4", "CXCL12")
C0_CXCR4_agg_table$sender <- factor(C0_CXCR4_agg_table$sender, levels = c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18",
                                                                              "C3","C4","C7","C13","C15","C6","C12","C17"))

# visualization of cell-cell communication result, here we choose to use bubble plot, with the bubble size indicating the proportion of cell pairs
# involved in interaction, and bubble color  indicates the overall interaction intensity normalized by the total number of cell pairs formed between
# two cell types, so that each bubble can reveal how a cell type as a whole was influenced by a certain interaction with another cell type within a 
# spatial niche.
ggplot(filter(C0_CXCR4_agg_table, C0_CXCR4_agg_table$CN_cluster1==c("CN5","CN6")), 
           aes(x = CN_cluster1, y = sender, color = group_mean, size = group_prop)) + 
      scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))+
      geom_point(stat = 'identity') + 
      xlab("") + ylab("") + 
      ggtitle("CXCL12-CXCR4") + 
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Figure 4e was generated in a similar way to Figure 4d ###

### Figure 4f ###
# FOV images from CosMx SMI were plotted using Spyrrow in python.