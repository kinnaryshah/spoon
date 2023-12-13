#' Stabilize weights
#'
#' Stabilize weights to mitigate extrapolation
#'
#' @details This function stabilizes the weights that are out of bounds

#' @param weights_list, output from generate_weights, contains list of s_g, r_tilda, lambda_hat, and weights matrix
#'
#' @return stabilized weights matrix
#'
#' @export
#'
#' @examples
#' stabilize_weights(weights_list)
#'
stabilize_weights <- function(weights_list){

  s_g <- weights_list[[1]]
  r_tilda <- weights_list[[2]]
  lambda_hat <- weights_list[[3]]
  w <- weights_list[[4]]

  spline_fit <- smooth.spline(y=sqrt(s_g), x=r_tilda)

  #get dimensions of spe
  n <- dim(w)[1] # Number of Cells
  G <- dim(w)[2] # Number of Genes

  #find min and max
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
    print(i)
    for (j in 1:ncol(lambda_hat)) {
      print(j)
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

  print(count_changes/(n*G))

  w_new <- tmp_pred_sqrt_sg^(-4)

  return(w_new)
}
