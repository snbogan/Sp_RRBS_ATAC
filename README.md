# Sp_RRBS_ATAC
This repository is work in process, and will contain all code and data used by Sam Bogan, Marie Strader, and Gretchen Hofmann for analyzing the effect of differential methylation on expression and splicing in the purple urchin during transgenerational plasticity. Analyses also include multifactorial models of transcription as a function of methylation and chromatin state (i.e., ATAC-seq data).

In the table of contents below, each lettered entry represents a folder in this repo with an R markdown file, its knitted report, and the data that it imports.

Table of Contents:

A. A_Differential_Exp_Splicing_Meth (estimations of differential expression, differential splicing, and differential methylation)
        
        1. A1_edgeR_DE_Sp_RRBS.Rmd/.html(differential expression)
                i. Input data: matrix of RNA-seq read counts and an R markdown (edgeR_DE_Sp_RRBS_ATAC) for modeling and reporting differential expression           
                of these reads.
                ii. Code: R markdown documenting edgeR differential expression analysis
                iii. html knit of A1_edgeR_DE_Sp_RRBS.Rmd
                iv. Outut data: two csvs for differential expression data for contrasts between maternal and developmental treatment groups
