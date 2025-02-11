
# spoon

<!-- badges: start -->

[![R-CMD-check](https://github.com/kinnaryshah/spoon/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kinnaryshah/spoon/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<img src="vignettes/hex.png" width="50"/>

## Introduction

`spoon` is a method to address the mean-variance relationship in
spatially resolved transcriptomics data. Current approaches rank
spatially variable genes based on either p-values or some effect size,
such as the proportion of spatially variable genes. However, previous
work in RNA-sequencing has shown that a technical bias, referred to as
the “mean-variance relationship”, exists in these data in that the
gene-level variance is correlated with mean RNA expression. We found
that there is a “mean-variance relationship” in spatial transcriptomics
data, and so we propose `spoon`, a statistical framework to prioritize
spatially variable genes that is not confounded by this relationship. We
fit a spline curve to estimate the mean-variance relationship. Then,
similar to using weights in a weighted least squares model, we used
weights that we plugged into a Gaussian Process Regression model fit
with a nearest-neighbor Gaussian process model to the preprocessed
expression measurements for each gene, i.e. one model per gene. `spoon`
removes the bias and leads to a more informative set of spatially
variable genes.

The `generate_weights()` function calculates individual observation
weights, where an individual observation is a UMI (unique molecular
identifier) count value for a specific gene and sample. If the desired
SVG detection method accepts weights, then the individual observation
weights can be used as inputs. If the desired SVG detection method does
not accept weights, then the Delta method is leveraged to rescale the
data and covariates by the weights. These scaled data and covariates are
used as inputs into the desired SVG detection function.

Bioconductor houses the infrastructure to store and analyze spatially
resolved transcriptomics data for R users, including many SVG detection
methods. This method addresses the mean-variance relationship
confounding SVG detection, which is related to these other Bioconductor
packages. Additionally, `spoon` is inspired by `limma::voom()` , which
is a popular Bioconductor package.

## Installation

The following code will install the latest release version of the
`spoon` package from Bioconductor. Additional details are shown on the
[Bioconductor](https://bioconductor.org/packages/spoon) page.

``` r
install.packages("BiocManager")
BiocManager::install("spoon")
```

The latest development version can also be installed from the `devel`
version of Bioconductor or from
[GitHub](https://github.com/kinnaryshah/spoon).

## Input data format

We recommend the input data be provided as a
[SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment)
Bioconductor object. The outputs are stored in the `rowData` of the
`SpatialExperiment` object. The examples below use this input data
format.

The inputs can also be provided as a numeric matrix of raw counts and a
numeric matrix of spatial coordinates.

## Tutorial

**Load packages and data**

``` r
library(nnSVG)
library(STexampleData)
```

    ## Loading required package: ExperimentHub

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: AnnotationHub

    ## Loading required package: BiocFileCache

    ## Loading required package: dbplyr

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:ExperimentHub':
    ## 
    ##     cache

    ## The following object is masked from 'package:AnnotationHub':
    ## 
    ##     cache

    ## Loading required package: SpatialExperiment

``` r
library(SpatialExperiment)
library(BRISC)
```

    ## Loading required package: RANN

    ## Loading required package: parallel

    ## Loading required package: rdist

    ## Loading required package: pbapply

    ## The ordering of inputs x (covariates) and y (response) in BRISC_estimation has been changed BRISC 1.0.0 onwards.
    ##   Please check the new documentation with ?BRISC_estimation.

``` r
library(BiocParallel)
library(scuttle)
library(Matrix)
```

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(spoon)

spe <- Visium_mouseCoronal()
```

    ## see ?STexampleData and browseVignettes('STexampleData') for documentation

    ## loading from cache

**Preprocessing**

``` r
# keep spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

# filter out low quality genes
spe <- filter_genes(spe)
```

    ## Gene filtering: removing mitochondrial genes

    ## removed 13 mitochondrial genes

    ## Gene filtering: retaining genes with at least 3 counts in at least 0.5% (n = 14) of spatial locations

    ## removed 21883 out of 32272 genes due to low expression

``` r
# calculate logcounts (log-transformed normalized counts) using scran package
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

# choose a small number of genes for this example to run quickly
set.seed(3)
ix_random <- sample(seq_len(nrow(spe)), 10)
spe <- spe[ix_random, ]

# remove spots with zero counts
spe <- spe[, colSums(logcounts(spe)) > 0]
```

**Step 1: generate weights**

``` r
weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 1,
                                                     RNGseed = 4))
```

    ## 21.8962962962963% of observations had their weight stabilized

**Step 2: weighted SVG detection**

``` r
spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 1, RNGseed = 5))
```

    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.
    ## Warning in nnSVG(spe[i, ], X = matrix(w[, i]), assay_name =
    ## "weighted_logcounts"): Rows (genes) and/or columns (spots) containing all zero
    ## counts have been found. Please see examples in tutorial for code to filter out
    ## zeros and/or low-expressed genes to avoid errors.

**Show results**

``` r
# display results
rowData(spe)
```

    ## DataFrame with 10 rows and 11 columns
    ##                               gene_id   gene_name    feature_type weighted_mean
    ##                           <character> <character>     <character>     <numeric>
    ## ENSMUSG00000030282 ENSMUSG00000030282        Cmas Gene Expression      2.648393
    ## ENSMUSG00000022601 ENSMUSG00000022601      Zbtb11 Gene Expression      0.950334
    ## ENSMUSG00000040220 ENSMUSG00000040220        Gas8 Gene Expression      0.403708
    ## ENSMUSG00000020704 ENSMUSG00000020704       Asic2 Gene Expression      1.320341
    ## ENSMUSG00000019173 ENSMUSG00000019173       Rab5c Gene Expression      2.830449
    ## ENSMUSG00000042156 ENSMUSG00000042156       Dzip1 Gene Expression      0.919617
    ## ENSMUSG00000050856 ENSMUSG00000050856       Atp5k Gene Expression     40.248319
    ## ENSMUSG00000033594 ENSMUSG00000033594     Spata2l Gene Expression      0.769744
    ## ENSMUSG00000037003 ENSMUSG00000037003        Tns2 Gene Expression      0.362365
    ## ENSMUSG00000026817 ENSMUSG00000026817         Ak1 Gene Expression      2.459033
    ##                    weighted_LR_stat weighted_sigma.sq weighted_tau.sq
    ##                           <numeric>         <numeric>       <numeric>
    ## ENSMUSG00000030282        502.84436         1.2603218        1.513975
    ## ENSMUSG00000022601        179.37299         0.7362656        1.079308
    ## ENSMUSG00000040220          7.10898         0.0161402        0.763033
    ## ENSMUSG00000020704        378.39817         1.2209177        1.177761
    ## ENSMUSG00000019173         38.96786         0.1380796        1.397702
    ## ENSMUSG00000042156         87.34277         0.1366747        1.315811
    ## ENSMUSG00000050856        141.25775         3.4727291       11.651505
    ## ENSMUSG00000033594        308.28764         0.6077159        1.038926
    ## ENSMUSG00000037003        228.10584         0.8156734        0.582560
    ## ENSMUSG00000026817        190.23953         0.3057840        1.490159
    ##                    weighted_prop_sv weighted_phi weighted_padj weighted_rank
    ##                           <numeric>    <numeric>     <numeric>     <numeric>
    ## ENSMUSG00000030282        0.4542851     6.153181   0.00000e+00             1
    ## ENSMUSG00000022601        0.4055278     2.217423   0.00000e+00             6
    ## ENSMUSG00000040220        0.0207145     2.117208   2.85960e-02            10
    ## ENSMUSG00000020704        0.5089959     1.872951   0.00000e+00             2
    ## ENSMUSG00000019173        0.0899084    16.532180   3.45332e-09             9
    ## ENSMUSG00000042156        0.0940971     4.405543   0.00000e+00             8
    ## ENSMUSG00000050856        0.2296135    27.166381   0.00000e+00             7
    ## ENSMUSG00000033594        0.3690638     2.204040   0.00000e+00             3
    ## ENSMUSG00000037003        0.5833599     0.179875   0.00000e+00             4
    ## ENSMUSG00000026817        0.1702637    11.461125   0.00000e+00             5

**Other SVG detection tools**

We provided a function to use the weights with
[nnSVG](https://www.nature.com/articles/s41467-023-39748-z) for more
accurate spatially variable gene detection. The weights can also be used
with other spatially variable gene detection tools using the following
procedure:

``` r
assay_name <- "logcounts"
weighted_logcounts <- t(weights)*assays(spe)[[assay_name]]
assay(spe, "weighted_logcounts") <- weighted_logcounts
```

`weighted_logcounts` can be accessed from
`assay(spe, "weighted_logcounts")`. Then, `weighted_logcounts` should be
used as the input counts matrix and `weights` as the input covariate
matrix in a spatially variable detection tool.

## Session information

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Ventura 13.6.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] spoon_1.1.3                 Matrix_1.7-0               
    ##  [3] scuttle_1.15.4              BiocParallel_1.39.0        
    ##  [5] BRISC_1.0.6                 pbapply_1.7-2              
    ##  [7] rdist_0.0.5                 RANN_2.6.2                 
    ##  [9] STexampleData_1.13.3        SpatialExperiment_1.15.1   
    ## [11] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.3
    ## [13] Biobase_2.65.1              GenomicRanges_1.57.1       
    ## [15] GenomeInfoDb_1.41.2         IRanges_2.39.2             
    ## [17] S4Vectors_0.43.2            MatrixGenerics_1.17.0      
    ## [19] matrixStats_1.4.1           ExperimentHub_2.13.1       
    ## [21] AnnotationHub_3.13.3        BiocFileCache_2.13.0       
    ## [23] dbplyr_2.5.0                BiocGenerics_0.51.3        
    ## [25] nnSVG_1.9.0                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1        dplyr_1.1.4             blob_1.2.4             
    ##  [4] filelock_1.0.3          Biostrings_2.73.2       fastmap_1.2.0          
    ##  [7] digest_0.6.37           mime_0.12               lifecycle_1.0.4        
    ## [10] KEGGREST_1.45.1         RSQLite_2.3.7           magrittr_2.0.3         
    ## [13] compiler_4.4.1          rlang_1.1.4             tools_4.4.1            
    ## [16] utf8_1.2.4              yaml_2.3.10             knitr_1.48             
    ## [19] S4Arrays_1.5.10         bit_4.5.0               curl_5.2.3             
    ## [22] DelayedArray_0.31.14    abind_1.4-8             withr_3.0.1            
    ## [25] purrr_1.0.2             grid_4.4.1              fansi_1.0.6            
    ## [28] beachmat_2.21.6         cli_3.6.3               rmarkdown_2.28         
    ## [31] crayon_1.5.3            generics_0.1.3          rstudioapi_0.16.0      
    ## [34] httr_1.4.7              rjson_0.2.23            DBI_1.2.3              
    ## [37] cachem_1.1.0            zlibbioc_1.51.1         AnnotationDbi_1.67.0   
    ## [40] BiocManager_1.30.25     XVector_0.45.0          vctrs_0.6.5            
    ## [43] jsonlite_1.8.9          bit64_4.5.2             magick_2.8.5           
    ## [46] glue_1.8.0              codetools_0.2-20        BiocVersion_3.20.0     
    ## [49] UCSC.utils_1.1.0        tibble_3.2.1            pillar_1.9.0           
    ## [52] rappdirs_0.3.3          htmltools_0.5.8.1       GenomeInfoDbData_1.2.13
    ## [55] R6_2.5.1                evaluate_1.0.0          lattice_0.22-6         
    ## [58] png_0.1-8               memoise_2.0.1           Rcpp_1.0.13            
    ## [61] SparseArray_1.5.44      xfun_0.48               pkgconfig_2.0.3
