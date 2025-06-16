### Statistical tests used in analysis ###

### load essential packages ###
library(Seurat)
library(FSA)

# Identify DEGs expressed in each cell cluster (Prerpocessing)
Lymphoma_markers <- FindAllMarkers(Lymphoma_data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
Lymphoma_markers_top30 <- Lymphoma_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n=30)

# Comparing the normalized PD-L1 : PD1 interaction proportion in each spatial niche (Figure 3)
kruskal.test(interac_prop ~ CN_cluster1, data = CD274_PDCD1_mye_T_sample_for_plot)

# Comparing C0_Tumor-B cell prop. in the FOVs from EBV+ and EBV- cases (Figure 5)
wilcox.test(Freq ~ Comparison1_nodal_EBV_status, data = filter(EBV_FOV_table_prop, EBV_FOV_table_prop$cell_state=="C0_Tumor-B"))

# Comparing cytoxicity signature score in T cells from CN1 of EBV+ and EBV- cases (Figure 5)
t.test(T_cytotoxicity_signature1 ~ Comparison1_nodal_EBV_status, data = T_cell_comparison_EBV_CN1_exp)

# DEG analysis between T cells from CN1 of EBV+ and EBV- groups (Figure 5)
Idents(T_cell_comparison_EBV_CN1) <- T_cell_comparison_EBV_CN1$Comparison1_nodal_EBV_status
T_cell_comparison_EBV_CN1_markers <- FindAllMarkers(T_cell_comparison_EBV_CN1,logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)

# Comparing C0_Tumor-B cell prop. in the FOVs from tumors in different anatomical sites (Figure 6)
kruskal.test(Freq ~ Comparison2_tumor_anatomical_site, data = filter(Tumor_site_FOV_table_prop, Tumor_site_FOV_table_prop$cell_state=="C0_Tumor-B"))
# For pairwaise comparisons, Dunn's post-hoc test was used. 
dunnTest(Freq ~ Comparison2_tumor_anatomical_site, data = filter(Tumor_site_FOV_table_prop, Tumor_site_FOV_table_prop$cell_state=="C0_Tumor-B"))

# Comparing the number of C0_Tumor-B in the neighborhoods centered by C0 (Figure 6)
kruskal.test(`C0_Tumor-B` ~ Comparison2_tumor_anatomical_site, data = C0_neighbor_log_comparison_tumor_site)
# For pairwaise comparisons, Dunn's post-hoc test was used.
dunnTest(`C0_Tumor-B` ~ Comparison2_tumor_anatomical_site, data = C0_neighbor_log_comparison_tumor_site)

# Comparing cytoxicity signature score in T cells from tumors in different anatomical sites (Figure 7)
pairwise.t.test(T_cell_comparison_tumor_site@meta.data$T_cytotoxicity_signature1, T_cell_comparison_tumor_site@meta.data$Comparison2_tumor_anatomical_site, p.adjust.method = "BH")
