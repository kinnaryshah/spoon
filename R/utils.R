#sub function for weighted_nnSVG()

weighted_nnSVG_calc_spe <- function(spe, w, i){
  res <- tryCatch({
    weight_output_i <- nnSVG(spe[i,],
                             X=matrix(w[,i]),
                             assay_name = "weighted_logcounts")
    list(weighted_LR_stat = rowData(weight_output_i)$LR_stat,
         weighted_sigma.sq = rowData(weight_output_i)$sigma.sq,
         weighted_tau.sq = rowData(weight_output_i)$tau.sq,
         weighted_prop_sv = rowData(weight_output_i)$prop_sv,
         weighted_phi = rowData(weight_output_i)$phi,
         weighted_padj = rowData(weight_output_i)$padj)
  }, error=function(e){
    message("ERROR :",conditionMessage(e), "\n")
    list(weighted_LR_stat = NA,
         weighted_sigma.sq = NA,
         weighted_tau.sq = NA,
         weighted_prop_sv = NA,
         weighted_phi = NA,
         weighted_padj = NA)
  })
  return(res)
}

weighted_nnSVG_calc_mat <- function(w_logcounts_mat, coords, w, i){
  res <- tryCatch({
    weight_output_i <- nnSVG(input = t(as.matrix(w_logcounts_mat[i,])),
                             spatial_coords = coords,
                             X=matrix(w[,i]))

    list(weighted_LR_stat = weight_output_i[11],
         weighted_sigma.sq = weight_output_i[1],
         weighted_tau.sq = weight_output_i[2],
         weighted_prop_sv = weight_output_i[9],
         weighted_phi = weight_output_i[3],
         weighted_padj = weight_output_i[14])
  }, error=function(e){
    message("ERROR :",conditionMessage(e), "\n")
    list(weighted_LR_stat = NA,
         weighted_sigma.sq = NA,
         weighted_tau.sq = NA,
         weighted_prop_sv = NA,
         weighted_phi = NA,
         weighted_padj = NA)
  })
  return(res)
}
