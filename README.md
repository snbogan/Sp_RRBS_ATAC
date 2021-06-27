# Sp_RRBS_ATAC

## Authors
* Samuel Bogan, University of California, Santa Barbara, Dept. of Ecology, Evolution, and Marine Biology (UCSB)
* Marie Strader, Auburn University, Dept. of Biological Sciences

## Description
This repository is a work in process, and will contain all code and data used by Sam Bogan, Marie Strader, and Gretchen Hofmann for analyzing the effect of differential methylation on expression and splicing in the purple urchin during transgenerational plasticity. Analyses also include multifactorial models of transcription as a function of methylation and chromatin state (i.e., ATAC-seq data).

In the table of contents below, each lettered entry represents a folder in this repo with an R markdown file, its knitted report, and the data that it imports and exports.

All R scripts were run using R version 3.6.1. This research was funded by a United States National Science Foundation award (IOS-1656262) to Gretchen Hofmann, UCSB. In addition, diving and boating resources were provided by Santa Barbara Coastal Long-Term Ecological Research program (NSF award OCE-1232779; Director: Dr. Daniel Reed).

## Table of Contents:

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
                    
        3. A3_edgeR_DM_Sp_RRBS.Rmd/.md(differential methylation)
                i. Input data: Bismark .cov files from each library and df's of CpG % methylation data grouped by feature type
                ii. Code: R markdown documenting edgeR differential methylation analyses of CpGs and genomic features and 
                binomial general linearized models testing effect of feature position on differential methylation (A3_edgeR_DM_Sp_RRBS.Rmd)
                iii. Markdown: knit of A3_edgeR_DM_Sp_RRBS.Rmd
                iv. Output data: twelve CSVs of differential methylation coefficients across CpGs and genomic feature type in response
                to maternal and developmental treatments
