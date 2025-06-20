# Lymphoma-spatial

This repository includes the codes used for generating the main results presented in the paper Multi-modal spatial characterization of tumor-immune microenvironments identifies targetable inflammatory niches in diffuse large B-cell lymphoma.

The code was originally developed in R (version 4.2.0), with softwares as detailed in the manuscript. To be compatible with updated versions, the code has also been tested in R (version 4.3.1), with the following required packages:
Seurat (version 5.1.0);
SeuratObject (version 5.0.2);
harmony (version 1.2.0);
ggplot2 (version 3.5.1);
dplyr (version 1.1.4);
tidyverse (version 2.0.0);
reshape2 (version 1.4.4);
plotly (version 4.10.4);
FSA (version 0.9.5).

Please refer to the folder R/ for related code. A detailed description of the function of each part of the code is provided in the related scripts.
We also provided demo data including around 80k cells and the cell- and sample-level metadata. It is reasonable to expect that the results generated based on demo data can be different from what was presented in the paper, while we hope the demo data can be used as input files for demonstrating the workflow of analysis.

Here we provide the outline of analysis covered by each script:  

Preprocessing.r covers the following analysis:  
> •	Data reading in  
> •	Data preprocessing  
> •	Data cleaning  
> •	DEG analysis (FindAllMarkers function in Seurat, with Wilcoxon rank sum test)  
> •	Cell type and state identification  

Figure 1.r covers the following analysis:  
> •	Dimension reduction (UMAP) plot showing major cell types  
> •	DEGs for major cell types  
> •	Comparative analysis of major cell type compositions across CosMx SMI and CODEX datasets  

Figure 2.r covers the following analysis:  
> • Dimension reduction (UMAP) plot showing cell states  
> • DEGs for cell states  
> •	Identification of cellular neighborhoods  
> •	Comparative analysis of neighborhood compositions of different cell states  
> •	Identification of spatial niches  
> •	Comparative analysis of cell state compositions across different spatial niches  
> •	Identifying DEGs across spatial niches  

Figure 3.r covers the following analysis:  
> •	Niche-specific cell function state analysis related to T cell chemotaxis, activation, and exhaustion  
> •	Neighborhood-based cell-cell communication analysis  
> •	Comparative analysis of PD-L1:PD-1 interactions across different spatial niches (Kruskal-Wallis test)  

Figure 4.r covers the following analysis:  
> •	Niche-specific cell function state analysis for C0_Tumor-B cells  
> •	Comparative analysis of CXCL12:CXCR4 interactions across different spatial niches  

Figure 5.r covers the following analysis:  
> •	Comparative analysis of cellular and spatial niche compositions between EBV+ and EBV- nodal lesions (Wilcoxon rank sum test)  
> •	Comparative analysis of neighborhood compositions of C0_Tumor-B cells between EBV+ and EBV- nodal lesions  
> •	Comparative analysis of T cell spatial locations and phenotypes between EBV+ and EBV- nodal lesions (Student's t-test for signature score; FindAllMarkers function in Seurat for DEGs)  

Figure 6-7.r covers the following analysis:  
> •	Comparative analysis of cellular and spatial niche compositions across different tumor anatomical sites (Kruskal-Wallis test with Dunn's post-hoc test for pairwise comparisons)  
> •	Comparative analysis of neighborhood compositions of C0_Tumor-B cells in different anatomical sites (Kruskal-Wallis test with Dunn's post-hoc test for pairwise comparisons)  
> •	Comparative analysis of T cell spatial locations, phenotypes, and immune checkpoint interactions across different tumor anatomical sites (Pairwise t-test)  


A brief workflow of the analysis is shown below:  
![Workflow](https://github.com/user-attachments/assets/a92c89ef-8781-4791-ac06-54a5f72da5e8)




For any questions, please leave your comment in GitHub or contact Yibo Dai (ydai6@mdanderson.org). We will help address the issues as soon as possible.
