#' Weighted nnSVG
#'
#' Run nnSVG for SVG detection using the weights
#'
#' @details This function incorporates the weights for each observation to run nnSVG

#' @param input either a SpatialExperiment object which contains a logcounts matrix, or a logcounts matrix
#' @param assay_name if using a SpatialExperiment object, name of the assay in which the logcounts matrix is stored
#' @param w weights matrix
#'
#' @return either spe with weighted nnSVG statistics, or matrix with weighted nnSVG statistics
#'
#' @import SpatialExperiment
#' @import nnSVG
#' @import purrr
#' @import scuttle
#' @export
#'
#' @examples
#' weighted_nnSVG(spe, "logcounts", w)
#'
weighted_nnSVG <- function(input, spatial_coords = NULL,
                           assay_name = "logcounts", w){

  # Make sure nnSVG fixed the interceptless model
  stopifnot(
    "Please update your nnSVG to minimum v1.5.3 to have the correct result" =
      packageVersion("nnSVG")>='1.5.3'
  )

  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }

  #compute weighted logcounts matrix and run nnSVG with covariate
  if(is(input, "SpatialExperiment")) {
    weighted_logcounts <- t(w)*logcounts(spe)
    assay(spe, "weighted_logcounts") <- weighted_logcounts
    weighted_mean <- rowMeans(weighted_logcounts)
    G <- dim(spe)[1]
    weighted_nnSVG_list <- map(.x=c(1:G), .f=~weighted_nnSVG_calc_spe(spe, w, .x))
  }
  else{
    weighted_logcounts_mat <- t(w)*input
    weighted_mean <- rowMeans(weighted_logcounts_mat)
    G <- dim(input)[1]
    weighted_nnSVG_list <- map(.x=c(1:G), .f=~weighted_nnSVG_calc_mat(weighted_logcounts_mat, spatial_coords, w, .x))
  }

  #return spe with weighted nnSVG params added to it
  res <- cbind(
    weighted_mean,
    weighted_LR_stat = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_LR_stat')])),
    weighted_sigma.sq = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_sigma.sq')])),
    weighted_tau.sq = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_tau.sq')])),
    weighted_prop_sv = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_prop_sv')]))
  )

  res <- cbind(res,
               weighted_rank = rank(-1*results[,"weighted_LR_stat"]))

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
