#' Weighted nnSVG
#'
#' Run nnSVG for SVG detection using the weights
#'
#' @details This function incorporates the weights for each observation to run nnSVG

#' @param input either a SpatialExperiment object which contains a logcounts matrix, or a logcounts matrix
#' @param spatial_coords matrix containing columns of spatial coordinates, needed if input is a matrix
#' @param assay_name if using a SpatialExperiment object, name of the assay in which the logcounts matrix is stored
#' @param w weights matrix
#' @param BPPARAM optional additional argument for parallelization to use BiocParallel

#' @return either spe with weighted nnSVG statistics, or matrix with weighted nnSVG statistics
#'
#' @import SpatialExperiment
#' @import nnSVG
#' @import purrr
#' @import scuttle
#' @import BiocParallel

#' @export
#'
#' @examples
#' weighted_nnSVG(spe, w=weights)
#'
#' weighted_nnSVG(logcounts_mat, coords_mat, w=weights)
#'
weighted_nnSVG <- function(input, spatial_coords = NULL,
                           assay_name = "logcounts", w,
                           BPPARAM = MulticoreParam(workers = 1)){

  # Make sure nnSVG fixed the interceptless model
  stopifnot(
    "Please update your nnSVG to minimum v1.5.3 to have the correct result" =
      packageVersion("nnSVG")>='1.5.3'
  )

  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }

  if(is(input, "SpatialExperiment")) {
    weighted_logcounts <- t(w)*logcounts(spe)
    assay(spe, "weighted_logcounts") <- weighted_logcounts
    weighted_mean <- rowMeans(weighted_logcounts)
  }
  else{
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
      weighted_nnSVG_i <- weighted_nnSVG_calc_mat(weighted_logcounts_mat, spatial_coords, w, i)
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
    weighted_prop_sv = unlist(weighted_nnSVG_output[,"weighted_prop_sv"])
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
