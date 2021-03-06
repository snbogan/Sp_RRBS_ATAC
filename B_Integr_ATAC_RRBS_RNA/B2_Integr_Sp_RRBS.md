B2\_Integr\_Sp\_RRBS
================
Sam Bogan
6/27/2021

``` r
# Load required packages
library( tidyverse )
```

    ## Warning: package 'tidyverse' was built under R version 3.6.2

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.2     ✓ dplyr   1.0.6
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 3.6.2

    ## Warning: package 'tibble' was built under R version 3.6.2

    ## Warning: package 'tidyr' was built under R version 3.6.2

    ## Warning: package 'readr' was built under R version 3.6.2

    ## Warning: package 'purrr' was built under R version 3.6.2

    ## Warning: package 'dplyr' was built under R version 3.6.2

    ## Warning: package 'forcats' was built under R version 3.6.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library( ape )
```

    ## Warning: package 'ape' was built under R version 3.6.2

# Read and wrangle chromatin accessibility data

``` r
# Read in chromatin accesibility df containg genewise Tn5 insert densities per genomic feature
load( "Output_data/ATAC_summary_peakSummaryALL.Rdata" )

# Convert NA counts to 0's for each genomic feature type
peaks_summary_all$peak_density_Introns[ is.na( peaks_summary_all$peak_density_Introns ) ] = 0
peaks_summary_all$peak_density_Intron1[ is.na( peaks_summary_all$peak_density_Intron1 ) ] = 0
peaks_summary_all$peak_density_TSS[ is.na( peaks_summary_all$peak_density_TSS ) ] = 0
peaks_summary_all$peak_density_Exons[ is.na( peaks_summary_all$peak_density_Exons ) ] = 0

# Divide all counts by 3 to create average across three replicates
peaks_summary_all$peak_density_Introns <- ( peaks_summary_all$peak_density_Introns) / 3
peaks_summary_all$peak_density_Intron1 <- ( peaks_summary_all$peak_density_Intron1 ) / 3
peaks_summary_all$peak_density_TSS <- ( peaks_summary_all$peak_density_TSS ) / 3
peaks_summary_all$peak_density_Exons <- ( peaks_summary_all$peak_density_Exons ) / 3

# Convert 'spu_id' column to 'geneid'
names( peaks_summary_all )[ names( peaks_summary_all ) == "spu_id" ] <- "geneid"

### Create an index of exon lengths and intron lengths for normalization of Tn5 insert counts
## Exon = length of summed exon features
## Intron length = gene length - summed exon lengths

# Read in gff file as df
Spur_3.1_gff <- read.gff( "Input_data/Strongylocentrotus_purpuratus.Spur_3.1.42.gff3" )

# Calc gene lengths
Spur_3.1_gff_genes <- filter( Spur_3.1_gff, type == "gene" )
Spur_3.1_gff_genes$geneid <- gsub( ".*ID=gene:", "", Spur_3.1_gff_genes$attributes )
Spur_3.1_gff_genes$geneid <- gsub( ";.*", "", Spur_3.1_gff_genes$geneid )
Spur_3.1_gff_genes$Len <- Spur_3.1_gff_genes$end - Spur_3.1_gff_genes$start

## Calc exon lenghts
# First, output length of each exon
Spur_3.1_gff_exons <- filter( Spur_3.1_gff, type == "exon" )
Spur_3.1_gff_exons$geneid <- gsub( ".*Parent=transcript:", "", Spur_3.1_gff_exons$attributes )
Spur_3.1_gff_exons$geneid <- gsub( "-tr.*", "", Spur_3.1_gff_exons$geneid )
Spur_3.1_gff_exons$Len <- Spur_3.1_gff_exons$end - Spur_3.1_gff_exons$start

## Next, sum all exon lengths per gene
Spur_3.1_gff_exons_sum <- aggregate( Spur_3.1_gff_exons$Len, 
                                     by = list( Category = Spur_3.1_gff_exons$geneid ), 
                                     FUN = sum )

# Adjust column names after aggregating exon legnths per gene
names( Spur_3.1_gff_exons_sum ) = c( "geneid", "Len" )

# Create df of geneids, gene lengths, and summed exon lengths
Spur_3.1_gff_gene_exon_len <- merge( data.frame( gene_len = Spur_3.1_gff_genes$Len,
                                                 geneid = Spur_3.1_gff_genes$geneid ),
                                     data.frame( geneid = Spur_3.1_gff_exons_sum$geneid,
                                                 exon_len = Spur_3.1_gff_exons_sum$Len ), 
                                      by = "geneid" )

# Calculate intron length
Spur_3.1_gff_gene_exon_len$intron_len <- Spur_3.1_gff_gene_exon_len$gene_len - 
  Spur_3.1_gff_gene_exon_len$exon_len

# Merge gene/exon/intron length info with peaks_summary_all
peaks_summary_all_corr <- merge( peaks_summary_all, Spur_3.1_gff_gene_exon_len, by = "geneid" )

# Normalize TSS, intron, and exon accessibility by length
peaks_summary_all_corr$TSS_dens_norm <- peaks_summary_all_corr$peak_density_TSS / 1000

peaks_summary_all_corr$intron_dens_norm <- as.numeric( 
  ifelse( peaks_summary_all_corr$intron_len == "0", "0",
          peaks_summary_all_corr$peak_density_Introns /
            peaks_summary_all_corr$intron_len ) )

peaks_summary_all_corr$exon_dens_norm <- peaks_summary_all_corr$peak_density_Exons /
  peaks_summary_all_corr$exon_len
```

# Combine chromatin accessibility and gene expression data

``` r
## Maternal correlations between promoter DM and transcript DE
# Read in maternal differential expression data
mat_edgeR_GE <- read.csv( "~/Documents/GitHub/Sp_RRBS_ATAC/A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_GE_table_filt.csv" )

# Change transcript ID to geneid
mat_edgeR_GE$X <- gsub( "transcript:", "", mat_edgeR_GE$X )
mat_edgeR_GE$X <- gsub( "-tr", "", mat_edgeR_GE$X )

# Subset maternal GE df down to geneid, logFC, and logCPM
keep_vars <- c( "X", "logFC", "logCPM" )
mat_edgeR_GE <- mat_edgeR_GE[ keep_vars ]

# Change 'X' column to 'geneid'
names( mat_edgeR_GE )[ names( mat_edgeR_GE ) == "X" ] <- "geneid"

# Merge maternal GE and chromatin access df's
ATAC_dens_by_GE <- merge( mat_edgeR_GE, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
ATAC_dens_by_GE[ is.na( ATAC_dens_by_GE ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Output maternal GE by chromatin access df
write.csv( ATAC_dens_by_GE, "Output_data/ATAC_dens_by_GE.csv", row.names = FALSE )
```

# Combine chromatin accessibility and transcript variant/splicing data

``` r
## Retrieve data on transcript variant number per gene from Spur_3.1 assembly and annotation
# Read in sp3_1 annotation for Spur_3.1
sp3_gff <- read.gff( "Input_data/sp3_1_GCF.gff3" )

# Filter annotation for transcripts only
sp3_transcripts <- filter( sp3_gff, type == "mRNA" )

# Create transcript id and geneid columns
sp3_transcripts$product <- gsub( ".*product=", "", sp3_transcripts$attributes )
sp3_transcripts$Locus <- gsub( ".*gene=", "", sp3_transcripts$attributes )
sp3_transcripts$Locus <- gsub( ";.*", "", sp3_transcripts$Locus )

# Measure number of transcript variants by counting rows per geneid
sp3_transcripts_freq <- data.frame( table( sp3_transcripts$Locus ) )

# Read in Echinobase geneinfo sheet
geneinfo <- read.csv( "Input_data/GenePageGeneralInfo_AllGenes.csv" )

# Filter for rows containing SPU_IDs -> SPU_clean
geneinfo_spu <- filter( geneinfo, grepl( "SPU_", SPU ) )
geneinfo_spu$ID <- gsub( ".*SPU", "", geneinfo_spu$SPU )
geneinfo_spu$ID <- substr( geneinfo_spu$ID, 1, 7 )

# Create geneid column w/ SPU_ID
geneinfo_spu$geneid <- paste( "SPU", geneinfo_spu$ID, sep = "" )

# Add SPU_ID/geneid to transcript variant frequency df
names( sp3_transcripts_freq )[ names( sp3_transcripts_freq ) == "Var1"] <- "Locus"
sp3_transcripts_freq <- merge( sp3_transcripts_freq,
                                    geneinfo_spu,
                                    by = "Locus" )

# Create df of geneid and # of transcript variants
spu_mRNA_var_df <- data.frame( geneid = sp3_transcripts_freq$geneid,
                               mRNA_var = sp3_transcripts_freq$Freq )

# Create binary variable for presence of transcript variants
spu_mRNA_var_df$splice <- ifelse( spu_mRNA_var_df$mRNA_var > 1, TRUE, FALSE )

# Merge transcript spu_mRNA_var_df with chromatin acessibility data
ATAC_dens_by_splicing <- merge( spu_mRNA_var_df, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
ATAC_dens_by_splicing[ is.na( ATAC_dens_by_splicing ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Export chromatin accesibility by constitutive splicing data
write.csv( ATAC_dens_by_splicing, "Output_data/ATAC_dens_by_splicing.csv", row.names = FALSE )
```

# Combine diff methylation, diff expression, and chromatin accessibility data

``` r
# Change 'X' column to 'geneid' and other ammendments necessary before merging two edgeR tables
names( mat_edgeR_GE )[ names( mat_edgeR_GE ) == "X" ] <- "geneid"
names( mat_edgeR_GE )[ names( mat_edgeR_GE ) == "logFC" ] <- "GE_logFC"
names( mat_edgeR_GE )[ names( mat_edgeR_GE ) == "logCPM" ] <- "GE_logCPM"

# Read in maternal promoter diff meth data
dm_promoter_df <- 
  read.csv( "~/Documents/GitHub/Sp_RRBS_ATAC/A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_1kb_promoter_DM_df.csv" )

# Subset maternal GE df down to geneid, logFC, and logCPM
dm_promoter_df <- dm_promoter_df[ keep_vars ]

# Change 'X' column to 'geneid'
names( dm_promoter_df )[ names( dm_promoter_df ) == "X" ] <- "geneid"
names( dm_promoter_df )[ names( dm_promoter_df ) == "logFC" ] <- "pr_meth_logFC"
names( dm_promoter_df )[ names( dm_promoter_df ) == "logCPM" ] <- "pr_meth_logCPM"

# Combine DM promoter and DE data (x = methylation; y = gene expression)
dm_promoter_by_DE <- merge( dm_promoter_df, mat_edgeR_GE, by = "geneid" ) # 1,114 rows

# Add chromatin accessibility data
dm_promoter_by_DE_ATAC <- merge( dm_promoter_by_DE, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_promoter_by_DE_ATAC[ is.na( dm_promoter_by_DE_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
## Maternal correlations between intron DM and transcript DE

# Read in maternal intron diff meth data
dm_intron_df <- 
  read.csv( "~/Documents/GitHub/Sp_RRBS_ATAC/A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_intron_gene_DM_df.csv" )

# Fix geneids
dm_intron_df$X <- gsub( "-I.*", "", dm_intron_df$X )

# Change 'X' column to 'geneid'
names( dm_intron_df )[ names( dm_intron_df ) == "X" ] <- "geneid"
names( dm_intron_df )[ names( dm_intron_df ) == "logFC" ] <- "int_meth_logFC"
names( dm_intron_df )[ names( dm_intron_df ) == "logCPM" ] <- "int_meth_logCPM"

# Combine DM promoter and DE data (x = methylation; y = gene expression)
dm_intron_by_DE <- merge( dm_intron_df, mat_edgeR_GE, by = "geneid" ) # 4,028 rows

# Add chromatin accessibility data
dm_intron_by_DE_ATAC <- merge( dm_intron_by_DE, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_intron_by_DE_ATAC[ is.na( dm_intron_by_DE_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Read in maternal exon diff meth data
dm_exon_df <- 
  read.csv( "~/Documents/GitHub/Sp_RRBS_ATAC/A_Differential_Exp_Splicing_Meth/Output_data/mat_edgeR_exon_gene_DM_df.csv" )

# Fix geneids
dm_exon_df$X <- gsub( "-I.*", "", dm_exon_df$X )

# Change 'X' column to 'geneid'
names( dm_exon_df )[ names( dm_exon_df ) == "X" ] <- "geneid"
names( dm_exon_df )[ names( dm_exon_df ) == "logFC" ] <- "ex_meth_logFC"
names( dm_exon_df )[ names( dm_exon_df ) == "logCPM" ] <- "ex_meth_logCPM"

# Combine DM promoter and DE data (x = methylation; y = gene expression)
dm_exon_by_DE <- merge( dm_exon_df, mat_edgeR_GE, by = "geneid" ) # 3,229 rows

# Add chromatin accessibility data
dm_exon_by_DE_ATAC <- merge( dm_exon_by_DE, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_exon_by_DE_ATAC[ is.na( dm_exon_by_DE_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
## Merge exon DM, intron DM, and maternal DE dfs
# Subset intron and exon DM df's down to geneid, logFC, and logCPM
keep_vars_dm <- c( 1, 2, 3 )
dm_intron_df <- dm_intron_df[ keep_vars_dm ]
dm_exon_df <- dm_exon_df[ keep_vars_dm ]

# Merge intron and exon DM df's
dm_intron_exon_by_DE <- merge( merge( dm_intron_df, dm_exon_df, by = "geneid" ),
                            mat_edgeR_GE, by = "geneid" )

# Add chromatin accessibility data
dm_intron_exon_by_DE_ATAC <- merge( dm_intron_exon_by_DE, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_intron_exon_by_DE_ATAC[ is.na( dm_intron_exon_by_DE_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Exports CSV's of diff meth by GE and chromatin access for each genomic feature type or group
write.csv( dm_promoter_by_DE_ATAC, "Output_data/dm_promoter_by_DE_ATAC.csv", row.names = FALSE )
write.csv( dm_intron_by_DE_ATAC, "Output_data/dm_intron_by_DE_ATAC.csv", row.names = FALSE )
write.csv( dm_exon_by_DE_ATAC, "Output_data/dm_exon_by_DE_ATAC.csv", row.names = FALSE )
write.csv( dm_intron_exon_by_DE_ATAC, "Output_data/dm_intron_exon_by_DE_ATAC.csv", row.names = FALSE )
```

# Combine chromatin accessibility and methylation data

Does not include gene expression data

``` r
# Merge promoter meth and chromatin access
# Add chromatin accessibility data
dm_promoter_by_ATAC <- merge( dm_promoter_df, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_promoter_by_ATAC[ is.na( dm_promoter_by_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Now, introns
dm_intron_by_ATAC <- merge( dm_intron_df, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_intron_by_ATAC[ is.na( dm_intron_by_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Now, exons
dm_exon_by_ATAC <- merge( dm_exon_df, peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_exon_by_ATAC[ is.na( dm_exon_by_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Lastly, introns and exons
dm_intron_exon_by_ATAC <- merge( merge( dm_intron_df, dm_exon_df, by = "geneid" ),
                            peaks_summary_all_corr, by = "geneid", all.x = TRUE )
dm_intron_exon_by_ATAC[ is.na( dm_intron_exon_by_ATAC ) ] <- 0
```

    ## Warning in `[<-.factor`(`*tmp*`, thisvar, value = 0): invalid factor level, NA
    ## generated

``` r
# Exports CSV's of diff meth by chromatin access, without GE, for each genomic feature type or group
write.csv( dm_promoter_by_ATAC, "Output_data/dm_promoter_by_ATAC.csv", row.names = FALSE )
write.csv( dm_intron_by_ATAC, "Output_data/dm_intron_by_ATAC.csv", row.names = FALSE )
write.csv( dm_exon_by_ATAC, "Output_data/dm_exon_by_ATAC.csv", row.names = FALSE )
write.csv( dm_intron_exon_by_DE_ATAC, "Output_data/dm_intron_exon_by_ATAC.csv", row.names = FALSE )
```
