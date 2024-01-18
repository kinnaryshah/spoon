#' Generate weights
#'
#' Generate weights on the observation level for each gene
#'
#' @details This function generates weights for each observation, which are used as input to scale the data and covariates

#' @param spe SpatialExperiment object, contains a raw counts matrix to generate weights from
#' @param spatial_coords matrix containing columns of spatial coordinates, needed if input is a matrix
#' @param assay_name if using a SpatialExperiment object, name of the assay in which the counts matrix is stored
#' @param stabilize when TRUE, stabilize weights to avoid extrapolation (highly recommended)
#' @param n_threads default = 1, number of threads for parallelization
#' @param BPPARAM optional additional argument for parallelization to use BiocParallel
#'
#' @return weights matrix
#'
#' @import SpatialExperiment
#' @import nnSVG
#' @import BRISC
#' @import BiocParallel
#' @import scuttle
#' @export
#'
#' @examples
#' library(nnSVG)
#' library(STexampleData)
#' library(ggplot2)
#' library(SpatialExperiment)
#' library(BRISC)
#' library(BiocParallel)
#' library(scuttle)
#'
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
#' weights_1 <- generate_weights(input = spe, stabilize = TRUE)
#'
#' #EXAMPLE 2 USING MATRIX
#'
#' counts_mat <- counts(spe)
#' logcounts_mat <- logcounts(spe)
#' coords_mat <- spatialCoords(spe)
#'
#' set.seed(1)
#' weights_2 <- generate_weights(input = counts_mat, spatial_coords = coords_mat, stabilize = TRUE)
#'
generate_weights <- function(input, spatial_coords = NULL,
                             assay_name = "counts",
                             stabilize = TRUE,
                             n_threads = 1, BPPARAM = NULL){

  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }

  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }

  #calculate weights

  # Count Matrix, transpose so each row is a spot, and each column is a gene
  if (is(input, "SpatialExperiment")) {
    r <- t(as.matrix(assays(spe)[[assay_name]]))
    coords <- spatialCoords(spe)

  } else {
    r <- t(as.matrix(input))
    coords <- spatial_coords
    row_names <- rownames(input)
  }

  n <- dim(r)[1] # Number of Cells
  G <- dim(r)[2] # Number of Genes

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

  if(stabilize){
    y_bar <- r_tilda

    max_ybar <- max(y_bar)
    min_ybar <- min(y_bar)

    s_g_max_ybar <- predict(spline_fit, x=max_ybar)$y
    s_g_min_ybar <- predict(spline_fit, x=min_ybar)$y

    #this matrix has same dimensions of lambda_hat
    tmp_pred_sqrt_sg <- predict(
      spline_fit,
      x = c(lambda_hat)
    )$y |>
      matrix(
        nrow = n, ncol = G
      )

    #constrain individual observation weights that have lambda hat more extreme than range of r_tilda
    count_changes <- 0
    for (i in 1:nrow(lambda_hat)) {
      for (j in 1:ncol(lambda_hat)) {
        #if this observation is greater than the max_ybar, change the weight matrix
        if(lambda_hat[i,j] > max_ybar){
          count_changes <- count_changes + 1
          tmp_pred_sqrt_sg[i,j] <- s_g_max_ybar
        }
        #if this observation is less than the min_ybar, change the weight matrix
        if(lambda_hat[i,j] < min_ybar){
          count_changes <- count_changes + 1
          tmp_pred_sqrt_sg[i,j] <- s_g_min_ybar
        }
      }
    }

    print(paste0(count_changes/(n*G)*100, "% of observations had their weight stabilized"))

    w <- tmp_pred_sqrt_sg^(-4)

  }

  return(w)
}
