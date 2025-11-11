# Functional insight into CDK7 in Cryptococcus neoformans

This study uncovers the role of the cyclin-dependent kinase **CDK7** in
the fungal pathogen *Cryptococcus neoformans*. Using genetic tagging,
phosphoproteomics, and transcriptomics, our study shows that fungal CDK7
forms an active CAK complex regulating transcription, splicing, and
cell-cycle progression. Human CDK7 inhibitors, especially **SY-1365
(Mevociclib)**, **BS-181** effectively block fungal growth and act
synergistically with existing antifungals, highlighting CDK7 as a
promising **new antifungal target**.

### CnCDK7 is more similar to humans

Alignment of amino acids sequences shows that CnCDK7 is similar to human
CDK7

![](README_files/figure-commonmark/percent_similar-1.png)

## RNA-seq analysis of ***C. neoformans*** under CDK7 inhibition

### Data processing and QC

### Alignment statistics

``` r
library(magrittr)

star_dir <- "/Users/pooja/Documents/CDK7_project/RNASeq/RNASe1_1h_Mev_BS181/BAM/"

star_align_log_files <- fs::dir_ls(star_dir, 
                                   glob = "*Log.final.out" ,
                                   recurse = T,type = "file")

names(star_align_log_files) <- gsub(x=basename(star_align_log_files), pattern ="*_star_alignLog.final.out", replacement = "")

parcutils::get_star_align_log_summary_plot(x = star_align_log_files,
                                col_total_reads = "red3", 
                                col_mapped_reads  = "blue4") 
```

![](README_files/figure-commonmark/align_stats-1.png)

#### Generate count matrix

``` r
source("run_rsubread.R")
library(magrittr)

gff_file <- "/Users/pooja/Documents/CDK7_project/RNASeq/RNASeq_1h_2022/GENOME/FungiDB-59_CneoformansH99.gff"


dir = list.files(star_dir, pattern = "*.bam$", full.names = TRUE)

cc <- run_rsubread(dir = dir, gff_file = gff_file)

readr::write_delim(cc, file = "/Users/pooja/Documents/CDK7_project/RNASeq/RNASe1_1h_Mev_BS181/Mev_BS181_1h_read_count.txt", delim = "\t")
```

#### DESeq2-based differential expression and normalization

``` r
library(magrittr)

cc <- readr::read_delim( "/Users/pooja/Documents/CDK7_project/RNASeq/RNASe1_1h_Mev_BS181/Mev_BS181_1h_read_count.txt", delim="\t", col_names = TRUE) %>% 
  dplyr::select(
    gene_name,dplyr::matches("^WT_DMSO|^WT_BS|^WT_Mev"))

sample_info <- cc %>% colnames() %>% .[-1]  %>%
              tibble::tibble(samples = . , groups = rep(c( "WT_BS181_1h", "WT_DMSO", "WT_Mev") , each = 3)) 


knitr::kable(sample_info)
```

| samples                    | groups      |
|:---------------------------|:------------|
| WT_BS181_1h_1              | WT_BS181_1h |
| WT_BS181_1h_2              | WT_BS181_1h |
| WT_BS181_1h_3              | WT_BS181_1h |
| WT_DMSO_Set1               | WT_DMSO     |
| WT_DMSO_Set2_sampled_reads | WT_DMSO     |
| WT_DMSO_Set3               | WT_DMSO     |
| WT_Mevo_Set1               | WT_Mev      |
| WT_Mevo_Set2               | WT_Mev      |
| WT_Mevo_Set3               | WT_Mev      |

``` r
control=c("WT_DMSO")
treatment = c("WT_BS181_1h", "WT_Mev")

res <- parcutils::run_deseq_analysis(counts = cc ,
                         sample_info = sample_info,
                         column_geneid = "gene_name" ,
                         cutoff_lfc = 1,
                         cutoff_pval = 0.05,
                         min_counts = 10,
                         group_numerator = treatment,
                         group_denominator = control,
                         regul_based_upon = 1)
```

#### Normalized count matrix

``` r
norm_mat <- parcutils::get_normalised_expression_matrix(x = res, 
                                            samples = NULL,
                                            genes = NULL,
                                            summarise_replicates = FALSE, 
                                            )
```

### Data quality assessment

#### Correlation between replicates

``` r
parcutils::get_corr_heatbox(x = res, show_corr_values = T, cluster_samples = F, plot_type = "lower")
```

![](README_files/figure-commonmark/QC-1.png)

#### PCA plot

``` r
parcutils::get_pca_plot(x = res, label_replicates = TRUE)
```

![](README_files/figure-commonmark/pca-1.png)

## **Functional insights into CnCDK7 using CDK7 inhibitors: SY-1365 & BS-181**

The following section highlights functional insights into the role of
CDK7 in *Cryptococcus neoformans* (Cn) using CDK7 inhibitors as
molecular tools. The data integrate multiple microbiology, molecular and
multi-omics approaches to elucidate CDK7 function in this WHO-priority
fungal pathogen, employing two human CDK7 inhibitors with demonstrated
antifungal potential.

1.  [SY-1365](https://github.com/sethiyap/CnCDK7_analysis/blob/main/Paper_SY-1365.md)
2.  [BS-181](https://github.com/sethiyap/CnCDK7_analysis/blob/main/Paper_BS181.md)
