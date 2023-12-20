#build small spe and matrix to quickly test out changes to functions
#following nnSVG vignette
library(nnSVG)
library(STexampleData)
library(ggplot2)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(purrr)

# load example dataset from STexampleData package
  spe <- Visium_humanDLPFC()
  dim(spe)

  # keep spots over tissue
  spe <- spe[, colData(spe)$in_tissue == 1]
  dim(spe)

  # filter low-expressed and mitochondrial genes
  # using function from nnSVG package with default filtering parameters
  spe <- filter_genes(spe)

  # calculate logcounts (log-transformed normalized counts) using scran package
  # using library size factors
  spe <- computeLibraryFactors(spe)
  spe <- logNormCounts(spe)
  assayNames(spe)

  known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
  ix_known <- which(rowData(spe)$gene_name %in% known_genes)
  ix <- c(ix_known)

  spe <- spe[ix, ]
  dim(spe)

  spe <- spe[, colSums(logcounts(spe)) > 0]

# run nnSVG using a single thread for this example workflow
#set.seed(123)
#spe <- nnSVG(spe, n_threads = 1)

# show results
#rowData(spe)

#make sure running weighted after and before nnSVG returns the same output!!
#currently something is hardcoded from nnSVG output

set.seed(1)
weights <- generate_weights(input = spe, stabilize = TRUE, n_threads = 1, BPPARAM = NULL)
spe <- weighted_nnSVG(spe, w=weights)
rowData(spe)

#verify that using counts matrix and coords matrix input gives same output as spe object

counts_mat <- counts(spe)
logcounts_mat <- logcounts(spe)
coords_mat <- spatialCoords(spe)

set.seed(1)
weights_2 <- generate_weights(input = counts_mat, spatial_coords = coords_mat, n_threads = 1, BPPARAM = NULL)
results <- weighted_nnSVG(logcounts_mat, coords_mat, w = weights_2)

