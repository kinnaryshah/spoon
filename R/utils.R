#sub function for weighted_nnSVG()

weighted_nnSVG_calc <- function(spe, i){
  res = tryCatch({
    print(i)
    weight_output_i <- nnSVG(spe[i,], X=matrix(w[,i]), assay_name = "weighted_logcounts")
    list(weighted_rank = rowData(weight_output_i)$rank,
         weighted_LR_stat = rowData(weight_output_i)$LR_stat,
         weighted_sigma.sq = rowData(weight_output_i)$sigma.sq,
         weighted_tau.sq = rowData(weight_output_i)$tau.sq,
         weighted_prop_sv = rowData(weight_output_i)$prop_sv)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    list(weighted_rank = NA,
         weighted_LR_stat = NA,
         weighted_sigma.sq = NA,
         weighted_tau.sq = NA,
         weighted_prop_sv = NA)
  })
  return(res)
}
