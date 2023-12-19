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

# select small set of random genes and several known SVGs for faster runtime in this example workflow
set.seed(123)
#ix_random <- sample(seq_len(nrow(spe)), 10)
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(spe)$gene_name %in% known_genes)
ix <- c(ix_known)

spe <- spe[ix, ]
dim(spe)


spe <- spe[, colSums(logcounts(spe)) > 0]

# run nnSVG using a single thread for this example workflow
set.seed(123)
spe <- nnSVG(spe, n_threads = 1)

# show results
rowData(spe)

#make sure running weighted after and before nnSVG returns the same output!!
#currently something is hardcoded from nnSVG output

weights <- generate_weights(spe, 1, NULL)
weights_new <- stabilize_weights(weights)
spe <- weighted_nnSVG(spe, "logcounts", weights_new) #this runs but gives warnings for each gene, zero row/col sums
