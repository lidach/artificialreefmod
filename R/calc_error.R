#' @title calc_error
#' 
#' @description function to calculate median absolute relative error, median relative error, and errors within 50% and 90%
#' 
#' @param true values from model that is being compared to, i.e. "true" estimates of metric
#' @param est values being compared to "true" values, "estimated" 
#' @param nsites number of sites in model

calc_error <- function(true, est, nsites){
  list <- list(abs_error = NA, rel_error = NA, MRE = NA, MARE = NA, ninety_lwr = NA, ninety_upr = NA, fifty_lwr = NA, fifty_upr = NA)
  # relative error
  for(i in 1:nsites){
    list$rel_error[i] <- (est[i] - true[i]) / true[i]
  }

  # median absolute relative error
  list$MARE <- median(abs(list$rel_error), na.rm =TRUE)
  # median relative error 
  list$MRE <- median(list$rel_error, na.rm =TRUE)

  # 50%, 90%
  list$ninety_lwr <- quantile(list$rel_error, c(0.05))
  list$ninety_upr <- quantile(list$rel_error, c(0.95))
  list$fifty_lwr <- quantile(list$rel_error, c(0.25))
  list$fifty_upr <- quantile(list$rel_error, c(0.75))
  return(list)
}
