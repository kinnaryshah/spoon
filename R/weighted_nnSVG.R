#' Weighted nnSVG
#'
#' Run nnSVG for SVG detection using the weights
#'
#' @details This function incorporates weights for each observation to run nnSVG

#' @param input either a SpatialExperiment object which contains a logcounts
#' matrix, or a logcounts matrix
#' @param spatial_coords matrix containing columns of spatial coordinates,
#' needed if input is a matrix
#' @param assay_name if using a SpatialExperiment object, name of the assay in
#' which the logcounts matrix is stored
#' @param w weights matrix
#' @param n_threads default = 1, number of threads for parallelization
#' @param BPPARAM optional additional argument for parallelization to use
#' BiocParallel

#' @return either spe with weighted nnSVG statistics, or matrix with weighted
#' nnSVG statistics
#'
#' @import SpatialExperiment
#' @import nnSVG
#' @import BiocParallel
#' @import Matrix

#' @export
#'
#' @examples
#' library(nnSVG)
#' library(STexampleData)
#' library(SpatialExperiment)
#' library(BRISC)
#' library(BiocParallel)
#' library(scuttle)
#' library(Matrix)
#'
#' spe <- Visium_humanDLPFC()
#'
#' # keep spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#'
#' # filter low-expressed and mitochondrial genes
#' spe <- filter_genes(spe)
#'
#' # calculate logcounts (log-transformed normalized counts) using scran package
#' spe <- computeLibraryFactors(spe)
#' spe <- logNormCounts(spe)
#'
#' known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
#' ix_known <- which(rowData(spe)$gene_name %in% known_genes)
#' ix <- c(ix_known)
#'
#' spe <- spe[ix, ]
#'
#' spe <- spe[, colSums(logcounts(spe)) > 0]
#'
#' #EXAMPLE 1 USING SPATIAL EXPERIMENT
#'
#' set.seed(1)
#' weights_1 <- generate_weights(input = spe,
#'                               stabilize = TRUE)
#' spe_results <- weighted_nnSVG(input = spe,
#'                               w = weights_1,
#'                               BPPARAM = MulticoreParam(workers = 1,
#'                                                        RNGseed = 4))
#'
#' # display results
#' rowData(spe_results)
#'
#'
#' #EXAMPLE 2 USING MATRIX
#'
#' counts_mat <- counts(spe)
#' logcounts_mat <- logcounts(spe)
#' coords_mat <- spatialCoords(spe)
#'
#' set.seed(1)
#' weights_2 <- generate_weights(input = counts_mat,
#'                               spatial_coords = coords_mat,
#'                               stabilize = TRUE)
#' results <- weighted_nnSVG(input = logcounts_mat,
#'                           spatial_coords = coords_mat,
#'                           w = weights_2,
#'                           BPPARAM = MulticoreParam(workers = 1, RNGseed = 4))
#'
#' # display results
#' print(results)
#'
weighted_nnSVG <- function(input, spatial_coords = NULL,
                           assay_name = "logcounts", w,
                           n_threads = 1,
                           BPPARAM = MulticoreParam(workers = 1)){

  # Make sure nnSVG fixed the interceptless model
  stopifnot(
    "Please update your nnSVG to minimum v1.5.3 to have the correct result" =
      packageVersion("nnSVG")>='1.5.3'
  )

  if (n_threads > 1) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }

  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }

  if(is(input, "SpatialExperiment")) {
    stopifnot(dim(t(w)) == dim(assays(spe)[[assay_name]]))
    weighted_logcounts <- t(w)*assays(spe)[[assay_name]]
    weighted_mean <- Matrix::rowMeans(weighted_logcounts)
    assay(spe, "weighted_logcounts") <- weighted_logcounts
  }
  else{
    stopifnot(dim(t(w)) == dim(input))
    weighted_logcounts_mat <- t(w)*input
    weighted_mean <- rowMeans(weighted_logcounts_mat)
  }

  #compute weighted logcounts matrix and run nnSVG with covariate
  #runs using parallelization
  ix <- seq_len(dim(spe)[1])
  output <- bplapply(ix, function(i) {

    if(is(input, "SpatialExperiment")) {
      weighted_nnSVG_i <- weighted_nnSVG_calc_spe(spe, w, i)
    }
    else{
      stopifnot(dim(spatial_coords)[1] == dim(input)[2])
      weighted_nnSVG_i <- weighted_nnSVG_calc_mat(weighted_logcounts_mat,
                                                  spatial_coords, w, i)
    }

    weighted_nnSVG_i

  }, BPPARAM = BPPARAM)

  # collapse output list into matrix
  weighted_nnSVG_output <- do.call("rbind", output)

  #return spe with weighted nnSVG params added to it
  res <- cbind(
    weighted_mean,
    weighted_LR_stat = unlist(weighted_nnSVG_output[,"weighted_LR_stat"]),
    weighted_sigma.sq = unlist(weighted_nnSVG_output[,"weighted_sigma.sq"]),
    weighted_tau.sq = unlist(weighted_nnSVG_output[,"weighted_tau.sq"]),
    weighted_prop_sv = unlist(weighted_nnSVG_output[,"weighted_prop_sv"]),
    weighted_phi = unlist(weighted_nnSVG_output[,"weighted_phi"]),
    weighted_padj = unlist(weighted_nnSVG_output[,"weighted_padj"])
  )

  res <- cbind(res,
               weighted_rank = rank(-1*res[,"weighted_LR_stat"]))

  if (is(input, "SpatialExperiment")) {
    # return in rowData of spe object
    stopifnot(nrow(spe) == nrow(res))
    rowData(spe) <- cbind(rowData(spe), res)
    spe
  } else {
    # return as numeric matrix
    stopifnot(nrow(input) == nrow(res))
    res
  }
}
