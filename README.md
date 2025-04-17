# Lymphoma-spatial

This repository includes the codes used for generating the main results presented in the paper Multi-modal spatial characterization of tumor-immune microenvironments identifies targetable inflammatory niches in diffuse large B-cell lymphoma.

The code was originally developed in R (version 4.2.0), with softwares as detailed in the manuscript. To be compatible with updated versions, the code has also been tested in R (version 4.3.1), with the following required packages:
Seurat (version 5.1.0);
SeuratObject (version 5.0.2);
ggplot2 (version 3.5.1);
dplyr (version 1.1.4);
tidyverse (version 2.0.0);
reshape2 (version 1.4.4);
plotly (version 4.10.4);
FSA (version 0.9.5).

Please refer to the folder R/ for related code. A detailed description of the function of each part of the code is provided in the related scripts.
We also provided demo data including around 80k cells and the cell- and sample-level metadata. It is reasonable to expect that the results generated based on demo data can be different from what was presented in the paper, while we hope the demo data can be used as input files for demonstrating the workflow of analysis.

For any questions, please leave your comment in GitHub or contact Yibo Dai (ydai6@mdanderson.org). We will help address the issues as soon as possible.
