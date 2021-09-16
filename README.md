# Sp_RRBS_ATAC

## Authors
* Samuel Bogan, University of California, Santa Barbara, Dept. of Ecology, Evolution, and Marine Biology (UCSB)
* Marie Strader, Auburn University, Dept. of Biological Sciences

## Description
This repository is a work in process, and will contain all code and data used by Sam Bogan, Marie Strader, and Gretchen Hofmann for analyzing the effect of differential methylation on expression and splicing in the purple urchin during transgenerational plasticity. Analyses also include multifactorial models of transcription as a function of methylation and chromatin state (i.e., ATAC-seq data).

All R scripts were run using R version 3.6.1. Input files for Section A (see table of contents) were produced using scripts available at: https://github.com/mariestrader/S.purp_RRBS_RNAseq_2019. Raw RRBS and RNA-seq read are available through the NCBI Short Read Archive under the accession PRJNA548926 and ATAC-seq .bed files are available under GEO experssion omnibus Bioproject PRJNA377768. ATAC-seq .bed files are large and must be imported to Sp_RRBS_ATAC/B_Integr_ATAC_RRBS_RNA/Input_data/ if Section B is to be run in full.

This research was funded by a United States National Science Foundation award (IOS-1656262) to Gretchen Hofmann, UCSB. This research was also supported by the Santa Barbara Channel LTER (OCE-1831937).

## Table of Contents:

    A. A_Differential_Exp_Splicing_Meth (estimations of differential expression, differential splicing, and differential methylation)
        
        1. A1_edgeR_DE_Sp_RRBS.Rmd/.md (differential expression)
                i. Input data: matrix of RNA-seq read counts
                ii. Code: R markdown documenting edgeR differential expression analysis (A1_edgeR_DE_Sp_RRBS.Rmd)
                iii. Markdown: knit of A1_edgeR_DE_Sp_RRBS.Rmd
                iv. Outut data: CSVs of diff expression from contrasts between maternal and developmental treatment groups
                
        2. A2_edgeR_DS_Sp_RRBS.Rmd/.md (differential exon use)
                i. Input data: matrix of RNA-seq exon read counts
                ii. Code: R markdown documenting edgeR differential exon use analysis (A2_edgeR_DS_Sp_RRBS.Rmd)
                iii. Markdown: knit of A2_edgeR_DS_Sp_RRBS.Rmd
                iv. Outut data: (a) CSVs of diff exon use from contrasts between maternal and developmental treatment groups, and (b)
                    dataframes containing DEU ~ exon number linear model coefficients for genes that do (non_p_coef_filt.csv) and   
                    don't (non_spur_models) exhibit spurious transciprtion and two dfs of -log p-vals with geneids for GO-term 
                    Mann Whitney U tests
                    
        3. A3_edgeR_DM_Sp_RRBS.Rmd/.md (differential methylation)
                i. Input data: Bismark .cov files from each library and df's of CpG % methylation data grouped by feature type
                ii. Code: R markdown documenting edgeR differential methylation analyses of CpGs and genomic features and 
                binomial general linearized models testing effect of feature position on differential methylation (A3_edgeR_DM_Sp_RRBS.Rmd)
                iii. Markdown: knit of A3_edgeR_DM_Sp_RRBS.Rmd
                iv. Output data: CSVs of differential methylation coefficients across CpGs and genomic feature type in response
                to maternal and developmental treatments
                
    B. B_Integr_ATAC_RRBS_RNA (integration of ATAC-seq, RRBS, and RNA-seq data)
        
        1. B1_ATAC_Summary_Sp_RRBS.Rmd (counting/summarizing ATAC-seq data from 39 hpf larvae)
                i. Input data: IMPORTED BY USER - 3 replicate .bed files of ATAC-seq transposase inserts from 39 hpf *S. purpuratus*
                ii. Code: RUN ON HIGH PERFORMANCE SYSTEM OR CLUSTER - B1_ATAC_summary_Sp_RRBS.Rmd R and bash scripts for summarizing 
                ATAC-seq data
                iii. Output data: 
                        a. .Rdata file containing intermediate files that is perdiodically saved while running .R script
                        b. ATAC_summary_peakSummaryALL.Rdata - summarized ATAC-seq insert densities per gene by genomic feature
                        c. ATAC_summary_annots_prism_cons.Rdata - summarized ATAC-seq data for consensus inserts only
                        d. .Rdata file containing annotations of ATAC-seq inserts
                        e. EXPORTED BY USER DURING EXECUTION - combined .bed file of all ATAC-seq replicates
                        
        2. B2_Integr_Sp_RRBS.Rmd (combinatorial integrations of ATAC-seq, RRBS, and RNA-seq data)
              i. Input data:
                        a. ATAC_summary_peakSummaryALL.Rdata
                        b. Strongylocentrotus_purpuratus.Spur_3.1.42.gff3.gz
                        c. sp3_1_GCF.gff3.gz
                        d. GenePageGeneralInfo_AllGenes.csv
                        e. (From section A output) A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_GE_table_filt.csv
                        f. (From section A output) A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_GE_table_filt.csv
                        g. (From section A output) A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_1kb_promoter_DM_df.csv
                        h. (From section A output) A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_intron_gene_DM_df.csv
                        i. (From section A output) A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_exon_gene_DM_df.csv
              iv. Output data: dataframes containing nearly all combinations of baseline and differential expression, splicing + baseline and differential measures of   
                               methylation across genic features, chromatin accessibility across genic features, and metrics for genic architecture.
                               
    C. C_Models_and_Figures (modeling and plotting gene regulation as a function of DNA methylation)
    
        1. C1_Models_Figures_Sp_RRBS.Rmd (Model fitting and selection with raw data and fitted data plots)
              i. Input data: this folder remains empty because input data are read in from B_Integr_ATAC_RRBS_RNA/Output_data
              ii. Code: MUST BE RUN WITH AT LEAST 3 CORES - C1_Models_Figures_Sp_RRBS.Rmd
              iii. GO_MWU-master: master directory for GO enrichment script by Wright, Aglyamova, Meyer, & Matz, 2015 (https://github.com/z0on/GO_MWU)
              iii. Output_data: Results of Fisher's exact and Mann Whitney U tests of GO-term enrichment among gene groups
              iv. Output_figures: Main-text and supplemental figures

