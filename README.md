# Sp_RRBS_ATAC
This repository is a work in process, and will contain all code and data used by Sam Bogan, Marie Strader, and Gretchen Hofmann for analyzing the effect of differential methylation on expression and splicing in the purple urchin during transgenerational plasticity. Analyses also include multifactorial models of transcription as a function of methylation and chromatin state (i.e., ATAC-seq data).

In the table of contents below, each lettered entry represents a folder in this repo with an R markdown file, its knitted report, and the data that it imports and exports.

Table of Contents:

A. A_Differential_Exp_Splicing_Meth (estimations of differential expression, differential splicing, and differential methylation)
        
        1. A1_edgeR_DE_Sp_RRBS.Rmd/.md(differential expression)
                i. Input data: matrix of RNA-seq read counts
                ii. Code: R markdown documenting edgeR differential expression analysis (A1_edgeR_DE_Sp_RRBS.Rmd)
                iii. Markdown: knit of A1_edgeR_DE_Sp_RRBS.Rmd
                iv. Outut data: two CSVs of diff expression from contrasts between maternal and developmental treatment groups
                
        2. A2_edgeR_DS_Sp_RRBS.Rmd/.md(differential exon use)
                i. Input data: matrix of RNA-seq exon read counts
                ii. Code: R markdown documenting edgeR differential exon use analysis (A2_edgeR_DS_Sp_RRBS.Rmd)
                iii. Markdown: knit of A2_edgeR_DS_Sp_RRBS.Rmd
                iv. Outut data: two CSVs of diff exon use from contrasts between maternal and developmental treatment groups, two
                    dataframes containing DEU ~ exon number linear model coefficients for genes that do (non_p_coef_filt.csv) and   
                    don't (non_spur_models) exhibit spurious transciprtion and two dfs of -log p-vals with geneids for GO-term 
                    Mann Whitney U tests
