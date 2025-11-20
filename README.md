# Gene-Expression-Analysis-Pipeline
Project Overview:
This project analyzes genomic microarray data to identify key molecular drivers of Leukemia. By comparing gene expression profiles between Leukemia patients and Normal healthy controls, I identified significant biomarkers, distinct molecular subtypes, and potential therapeutic drug targets.

Data Source
Dataset: Gene Expression Omnibus (GEO)

Accession: [GSE48558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48558)

Platform: Affymetrix Human Genome U133 Plus 2.0 Array

Sample Size: ~170 samples (Leukemia vs. Normal)

Languagues:

R: Data preprocessing, Normalization (RMA), Quality Control.

Python: Statistical analysis, Visualization.

Libraries: pandas, seaborn, matplotlib, gseapy (Enrichr), scikit-learn, lifelines.


![PCA Plot](Gene%20Expression%20Analysis%20Pipeline/figures/Three%20Plots%20of%20PC1,%20PC2,%20PC3.png)

![Volcano Plot](Gene%20Expression%20Analysis%20Pipeline/figures/Volcano%20Plot%20-%20Leukemia%20vs%20Normal.png)

![Heatmap](Gene%20Expression%20Analysis%20Pipeline/figures/Heatmap%20-%20Top%20100%20Significant%20Genes.png)

![GO_Down](Gene%20Expression%20Analysis%20Pipeline/figures/GO%20Biological%20Processes%20(Down-Regulated%20in%20Leukemia).png)

![GO_Up](Gene%20Expression%20Analysis%20Pipeline/figures/GO%20Biological%20Processes%20(Up-Regulated%20in%20Leukemia).png)

