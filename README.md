---
editor_options: 
  markdown: 
    wrap: 72
---

# spoon ![](vignettes/hex.png){width="50"}

## Introduction

`spoon` is a method to address the mean-variance relationship in
spatially resolved transcriptomics data. Current approaches rank
spatially variable genes based on either p-values or some effect size,
such as the proportion of spatially variable genes. However, previous
work in RNA-sequencing has shown that a technical bias, referred to as
the "mean-variance relationship", exists in these data in that the
gene-level variance is correlated with mean RNA expression. We found
that there is a "mean-variance relationship" in spatial transcriptomics
data, and so we propose `spoon`, a statistical framework to prioritize
spatially variable genes that is not confounded by this relationship. We
fit a spline curve to estimate the mean-variance relationship. Then,
similar to using weights in a weighted least squares model, we used
weights that we plugged into a Gaussian Process Regression model fit
with a nearest-neighbor Gaussian process model to the preprocessed
expression measurements for each gene, i.e. one model per gene. `spoon`
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

```{r, eval=FALSE}
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

```{r}
library(nnSVG)
library(STexampleData)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

spe <- Visium_mouseCoronal()
```

**Preprocessing**

```{r}

# keep spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

# filter out low quality genes
spe <- filter_genes(spe)

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

```{r}
weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 1,
                                                     RNGseed = 4))
```

**Step 2: weighted SVG detection**

```{r}
spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 1, RNGseed = 5))
```

**Show results**

```{r}
# display results
rowData(spe)
```

**Other SVG detection tools**

We provided a function to use the weights with
[nnSVG](https://www.nature.com/articles/s41467-023-39748-z) for more
accurate spatially variable gene detection. The weights can also be used
with other spatially variable gene detection tools using the following
procedure:

```{r, eval=FALSE}
assay_name <- "logcounts"
weighted_logcounts <- t(weights)*assays(spe)[[assay_name]]
assay(spe, "weighted_logcounts") <- weighted_logcounts
```

`weighted_logcounts` can be accessed from
`assay(spe, "weighted_logcounts")`. Then, `weighted_logcounts` should be
used as the input counts matrix and `weights` as the input covariate
matrix in a spatially variable detection tool.

## Session information

```{r}
sessionInfo()
```
