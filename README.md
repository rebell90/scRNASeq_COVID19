# scRNASeq_COVID19
scRNA-Seq steps 1-4

Single Cell Analysis of the GEO experiment GSE161918 "Time-resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19". In this step, only samples from initial timepoints are used in scRNAseq via the attached supplementary file "GSE161918_AllBatches_SeuratObj.rds.gz".
Data is loaded and converted to a Single Cell Experiment object and subsequently subsetted to timepoint T0 and Batch1.
The code for the project is contained in subfolders with the following topics:

1) libraries/configure/data-load
2) quality control
3) normalization
4) gene variance modeling
5) dimensionality reduction
6) clustering
7) marker gene detection




In following projects, I will aggregate the data for both bulk and "pseudo-bulk" analysis, including Differential Expression,
along with Functional Enrichment and Pathway Analysis.

This project was inspired by the eBook "Orchestrating single-cell analysis with Bioconductor" (see citation below):
Amezquita, R.A., Lun, A.T.L., Becht, E. et al. Orchestrating single-cell analysis with Bioconductor.
Nat Methods 17, 137–145 (2020). https://doi.org/10.1038/s41592-019-0654-x
