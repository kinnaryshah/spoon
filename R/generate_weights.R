#' Generate weights
#'
#' Generate weights on the observation level for each gene
#'
#' @details This function generates weights for each observation, which are used as input to scale the data and covariates

#' @param spe SpatialExperiment object, contains a raw counts matrix to generate weights from
#' @param n_threads default = 1, number of threads for parallelization
#' @param BPPARAM optional additional argument for parallelization to use BiocParallel
#'
#' @return list of s_g, r_tilda, lambda_hat, and weights matrix
#'
#' @import SpatialExperiment
#' @import nnSVG
#' @import BRISC
#' @import BiocParallel
#' @import scuttle
#' @export
#'
#' @examples
#' generate_weights(spe)
#'
generate_weights <- function(spe, n_threads = 1, BPPARAM = NULL){

  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }

  spe <- spe[, colSums(counts(spe)) > 0]
  dim(spe)

  spe <- logNormCounts(spe)

  #calculate weights

  # Count Matrix, transpose so each row is a spot, and each column is a gene
  r <- t(as.matrix(counts(spe)))

  n <- dim(spe)[2] # Number of Cells
  G <- dim(spe)[1] # Number of Genes

  # Sample-specific Library Size
  R <- rowSums(r)
  stopifnot(length(R)==n)

  # Temporary Matrix, replicate library size for each row
  tmp_R_mat <- matrix(
    rep(R, each = G),
    byrow = TRUE, nrow = n, ncol = G
  )

  # logCPM
  y <- log2(r+0.5) - log2(tmp_R_mat+1) + log2(10^6)

  coords <- spatialCoords(spe)

  # scale coordinates proportionally
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)

  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = "AMMD", verbose = F)

  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = 10, n_omp = 1,
                             search.type = "tree", ordering = order_brisc,
                             verbose = F)

  # run BRISC using parallelization
  # run BRISC by column of y so BRISC is run per gene
  ix <- seq_len(ncol(y))
  #ix <- c(1,2,3)
  out_brisc <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[ ,i]
    suppressWarnings({
      runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = NULL,
                                  cov.model = "exponential",
                                  ordering = order_brisc, neighbor = nn_brisc,
                                  verbose = F)
      })
    })


    pred_i <- BRISC_prediction(out_i, coords = coords, X.0 = NULL, verbose = F)
    residual_i <- y_i - pred_i$prediction

    return(list(pred_i$prediction, residual_i))
  }, BPPARAM = BPPARAM)

  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)

  # *Voom Variance Modelling -------------------------------------------------

  mu_hat <- unname(as.matrix(as.data.frame(mat_brisc[,1])))
  stopifnot(dim(mu_hat) == c(n, G))

  s_g <- unname(as.data.frame(mat_brisc[,2])) |>
    apply(MARGIN = 2,  # Column wise
          FUN = sd)
  stopifnot(length(s_g) == G)

  y_bar <- colMeans(mu_hat)
  stopifnot(length(y_bar) == G)

  # Geometric Mean
  R_tilda <- exp(mean(log(R)))
  # The reason of calculating log is to avoid integer overflow

  # Log2 Counts
  # Note: slight notation abuse. Prev r denotes read counts
  r_tilda <- y_bar + log2(R_tilda) - log2(10^6)
  stopifnot(length(r_tilda)==G)

  # *PREDICT MODEL -----------------------------------------------------------------
  stopifnot(dim(mu_hat)==dim(tmp_R_mat))
  lambda_hat <- mu_hat + log2(tmp_R_mat+1) - log2(10^6)

  #gives percentage of lambda_hat values out of range
  sum(lambda_hat < range(r_tilda)[1] | lambda_hat > range(r_tilda)[2]) / (dim(spe)[1]*dim(spe)[2])
  sum(lambda_hat < range(r_tilda)[1]) / (dim(spe)[1]*dim(spe)[2])
  sum(lambda_hat > range(r_tilda)[2]) / (dim(spe)[1]*dim(spe)[2])

  spline_fit <- smooth.spline(y=sqrt(s_g), x=r_tilda)

  # NOTE: It is possible that lambda is out of range of r_tilda
  # which will produce NA predicted values due to extrapolation
  tmp_pred_sqrt_sg <- predict(
    spline_fit,
    x = c(lambda_hat)
  )$y |>
    matrix(
      nrow = n, ncol = G
    )

  w <- tmp_pred_sqrt_sg^(-4)

  return(list(s_g, r_tilda, lambda_hat, w))
}
