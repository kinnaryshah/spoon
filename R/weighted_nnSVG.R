#' Weighted nnSVG
#'
#' Run nnSVG for SVG detection using the weights
#'
#' @details This function incorporates the weights for each observation to run nnSVG

#' @param input either a SpatialExperiment object which contains a logcounts matrix, or a matrix
#' @param assay_name name of the assay in which the logcounts matrix is stored
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
weighted_nnSVG <- function(input, assay_name = "logcounts", w){
  print(assay_name)
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }
  print("before weighting")
  weighted_logcounts <- t(w)*logcounts(spe)
  assay(spe, "weighted_logcounts") <- weighted_logcounts # assign a new entry to assays slot, nnSVG will use "logcounts" by default
  print(assayNames(spe))
  # Make sure nnSVG fixed the interceptless model
  stopifnot(
    "Please update your nnSVG to minimum v1.5.3 to have the correct result" =
      packageVersion("nnSVG")>='1.5.3'
  )
  print("before running")
  #run nnSVG with covariate
  weighted_nnSVG_list <- map(.x=c(1:dim(spe)[1]), .f=~weighted_nnSVG_calc(spe, w, .x))

  #need to manually calculate mean for weighted nnSVG
  weighted_mean <- rowMeans(weighted_logcounts)

  #return spe with weighted nnSVG params added to it
  #can add more nnSVG output if desired
  res <- cbind(
    weighted_mean,
    weighted_rank = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_rank')])),
    weighted_LR_stat = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_LR_stat')])),
    weighted_sigma.sq = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_sigma.sq')])),
    weighted_tau.sq = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_tau.sq')])),
    weighted_prop_sv = unlist(lapply(weighted_nnSVG_list, function (x) x[c('weighted_prop_sv')]))
  )

  if (is(input, "SpatialExperiment")) {
    # return in rowData of spe object
    stopifnot(nrow(spe) == nrow(res))
    rowData(spe) <- cbind(rowData(spe), res)
    spe
  } else {
    # return as numeric matrix
    stopifnot(nrow(input) == nrow(res))
    rownames(res) <- row_names
    res
  }
}
