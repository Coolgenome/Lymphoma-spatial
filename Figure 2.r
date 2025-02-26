### Figure 2 ###

### Figure 2a ###
library(Seurat)

Lymphoma_data <- readRDS("./Lymphoma_data.rds") 

DimPlot(Lymphoma_data, reduction = "umap",label = F, raster = T, group.by="cell_state", 
        cols=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5",
               "#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837",
               "#6a3d9a"))

### Figure 2b ###

cell_state_marker  <-c("MS4A1","CD79A","TCL1A","CD37","IGHM","TUBB","TUBB4B","TYMS","STMN1","MZB1","XBP1",
                       "IGHG1","IGHG2","IGHA1","CD44","JCHAIN","CD40","CD38","TNFRSF17","CD3D","CD3E","CD3G",
                       "CD2","NKG7","CCL5","CD68","CD163","APOC1","APOE","C1QC","C1QB","C1QA","SPP1","MMP12",
                       "DUSP1","FOS","JUN","CXCL8","IL1B","CXCL1/2/3","LGALS1","COL1A1","COL1A2","COL3A1",
                       "DCN","FN1","TIMP1","ACTA2","CCL21","CCL19","VWF","PECAM1","MGP","CLU","VCAM1",
                       "HSPA1A/B","HSPB1","KRT5","KRT6A/B/C","KRT17","KRT16","HBB","HBA1/2")

DotPlot(Lymphoma_data,features = cell_state_marker, scale.by = "size", group.by = "cell_state", scale.max = 60, scale = 10)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 2c was created with BioRender ###

### Figure 2d ###
### neighborhood analysis ###
# calculating spatial neighborhood #
Lymphoma.meta <- readRDS("./Lymphoma.meta.for_neighborhood.rds") ### with columns cell_barcode, EBV_status, Tumor_site

major.ord = unique(Lymphoma.meta$Sub_Cls)
neighborhood_radius = 200
my_neighbor_list <- list()
for (clstype in major.ord){
      my.neighbor = c()
      print (clstype)
      
      for(sam in unique(Lymphoma.meta$Slide_FOV)) {
        print (sam)
        tmp.meta = Lymphoma.meta[Lymphoma.meta$Sub_Cls == clstype & Lymphoma.meta$Slide_FOV == sam, ]
        tmp.meta1 = Lymphoma.meta[Lymphoma.meta$Slide_FOV == sam, ]
        
        if(clstype %in% unique(tmp.meta1$Sub_Cls)){
          
          tmp.mtx = matrix(0,nrow = length(major.ord),ncol = nrow(tmp.meta))
          rownames(tmp.mtx) = major.ord
          colnames(tmp.mtx) = rownames(tmp.meta)
          
          for(i in 1:nrow(tmp.meta)) {
            xloc = tmp.meta[i, 'CenterX_local_px']
            yloc = tmp.meta[i, 'CenterY_local_px']
            cutoff = neighborhood_radius
            idx1 = tmp.meta1[, 'CenterX_local_px'] <=  (xloc + cutoff) & tmp.meta1[,'CenterX_local_px'] >=  (xloc - cutoff)
            idx2 = tmp.meta1[, 'CenterY_local_px'] <=  (yloc + cutoff) & tmp.meta1[,'CenterY_local_px'] >=  (yloc - cutoff)
          
            square = tmp.meta1[idx1 & idx2,]
            dis = apply(square, 1, function(x) {sqrt((as.numeric(x['CenterX_local_px']) - as.numeric(xloc))^2 + (as.numeric(x['CenterY_local_px']) - as.numeric(yloc))^2)} )
            
            tmp.freq = table(tmp.meta1[match(names(dis[dis <= cutoff]), rownames(tmp.meta1)), 'Sub_Cls'])
            tmp.mtx[,i] = tmp.freq[match(major.ord,names(tmp.freq))]
          }
          
          tmp.mtx[which(is.na(tmp.mtx))] = 0
          
          my.neighbor = cbind(my.neighbor, tmp.mtx)
        }
      }
      
      my_neighbor_list[[clstype]] <- my.neighbor
      
      write.csv(t(my.neighbor), paste0(clstype,'.neighbor.sub.csv'),quote = F)
}

saveRDS(my_neighbor_list,"my_neighbor_list.rds")

# statistics of spatial neighborhood #
# aggregrate neighborhood cell counts for all cells #
cell_neighbor_df <- c()
for(i in 1:length(my_neighbor_list)){
   cell_neighbor_df <- cbind(cell_neighbor_df, my_neighbor_list[[i]])
}

# aggregrate neighborhood cell counts by center cell type #
Neighbor_sum <- matrix(0, nrow=19, ncol=19)
rownames(Neighbor_sum) = major.ord
colnames(Neighbor_sum) = major.ord

for(i in major.ord){
    df_filter <- filter(cell_neighbor_df, cell_neighbor_df$Sub_Cls == i)
    df_filter <- df_filter[,major.ord]
    df_filter_sum = colSums(df_filter)
    Neighbor_sum[i,] = df_filter_sum
}
    
Neighbor_sum <- as.data.frame(Neighbor_sum)
Neighbor_sum$center_cell_type <- rownames(Neighbor_sum)

# calculate cell prop. in each neighborhood #
cell.prop1 <- as.data.frame(prop.table(as.matrix(Neighbor_sum[,0:19]), margin=1))
cell.prop1$center_cell_type <- rownames(cell.prop1)

cell.prop1 <- reshape2::melt(cell.prop1, id.vars=c("center_cell_type"),
                             measure.vars=c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18","C3","C4","C7","C13","C15","C6","C12","C17"),
                             variable.name = "neighbor_cell_type", value.name="neighbor_cell_prop")

cell.prop1$center_cell_type <- factor(cell.prop1$center_cell_type, levels=c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18","C3","C4","C7","C13","C15","C6","C12","C17"))              
cell.prop1$neighbor_cell_type <- factor(cell.prop1$neighbor_cell_type, levels=c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18","C3","C4","C7","C13","C15","C6","C12","C17"))              

# plot #
ggplot(cell.prop1,aes(center_cell_type,neighbor_cell_prop,fill=neighbor_cell_type))+
       geom_bar(stat = "identity",position = "fill")+
       ggtitle("Neighboring cell proportion for each cell type")+
       theme_classic()+
       theme(axis.ticks.length = unit(0.5,'cm'))+
       theme(axis.text.x = element_text(angle=60,hjust = 1))+
       guides(fill=guide_legend(title = "Cell type"))+
       scale_fill_manual(values = c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5","#bc80bd",
                                    "#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))


### Figure 2e-g ###
### k-means clustering based on neighborhood matrix ###
spatial_niche <- readRDS("spatial_niche.rds")  ### provide spatial_niche.rds (cell level table, with columns cell_barcode, cell_state, spatial_niche, Sample_ID, EBV_status, Tumor_site)

# project CN clusters in PCA space (Figure 2e) #
cell_neighbor_df_scaled <- scale(t(cell_neighbor_df)) ### need to check matrix if transpose is needed
cell_neighbor_pca <- prcomp(cell_neighbor_df_scaled)
cell_neighbor_pca_coord <- as.data.frame(cell_neighbor_pca$x)
cell_neighbor_pca_coord <- rownames_to_column(cell_neighbor_pca_coord, var="center_cell_barcode")
cell_neighbor_pca_coord <- left_join(cell_neighbor_pca_coord, spatial_niche, by="center_cell_barcode")

library(plotly)
plot_ly(cell_neighbor_pca_coord, x = ~PC1, y = ~PC2, z = ~PC3, color = ~CN_cluster_new,
        colors=c("#ed2224","#fbb14d","#3B50A3","#eee817","#79c479","#55C7F3","#dad9d9"),
        alpha=0.3, sizes= c(100,100))


# plot cell composition in each CN (Figure 2f) #
cellnum1 <- table(spatial_niche$cell_state, spatial_niche$spatial_niche)
cellnum1
cell.prop1<-as.data.frame(prop.table(cellnum1))
colnames(cell.prop1)<-c("cell_state","spatial_niche","Proportion")
cell.prop1$Cell_type <- factor(cell.prop1$Cell_type, levels=c("C0","C5","C9","C11","C14","C1","C2","C8","C10","C16","C18","C3","C4","C7","C13","C15","C6","C12","C17"))

ggplot(cell.prop1,aes(Proportion,spatial_niche,fill=cell_state))+
    geom_bar(stat = "identity",position = "fill")+
    ggtitle("Cell composition in each spatial niche")+
    theme_bw()+
    scale_fill_manual(values=c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5","#bc80bd",
                               "#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837","#6a3d9a"))+
    theme(axis.ticks.length = unit(0.5,'cm'))+
    guides(fill=guide_legend(title = NULL))

# DEGs among spatial niches #
Lymphoma_data <- AddMetaData(Lymphoma_data, select(cell_neighbor_km_cluster,"spatial_niche"))

CN_marker <- c("CD3D","CD3E","CD2","CD8A","GZMB","GZMK","NKG7","PRF1","GZMH","GXCL9","CXCL10","CCL5","CCL2","MZB1","IGHG1",
               "IGHG2","IGKC","C1QA","C1QB","C1QC","LYZ","APOE","APOC1","S100A9","COL1A1","COL1A2","COL3A1","FN1","DCN",
               "ACTA2","TIMP1","COL6A1","COL6A2","CD79A","MS4A1","CD19","TCL1A","IGHM","CD37","STMN1","TYMS","TOP2A","BCL2",
               "BTK","SYK","CD24","CD52") 

DotPlot(Lymphoma_data,features = CN_marker, scale.by = "size", group.by = "spatial_niche", scale.max = 60, scale = 10)+
    RotatedAxis()+
    scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))

### Figure 2h ###
### FOV images from CosMx SMI were plotted using Spyrrow in python, see the document "CosMx data visualization" for detailed information.