#' DIC
#' 
#' @importFrom extraDistr dgev
#' 
#' @description 
#' This function calculates the DIC of the GEV \code{model} 
#' 
#' @param model object from \code{\link{GEVmodel}}
#' @return DIC
#' @export 
DIC <- function(model) {
  
  Y <- model$y
  X <- model$x
  
  n.sims <- nrow(model$params)
  p      <- ncol(X)
  
  mcmc_mean <- colMeans(model$params)
  
  deviance <- rep(0, n.sims) 
  for (n in 1:n.sims) {
    deviance[n] <- - 2 * sum(extraDistr::dgev(Y, mu = model$params[n, 1] + X %*% model$params[n, 3 + 1:p], sigma = model$params[n, 2], xi = model$params[n, 3], log = TRUE))
  }
  
  aux <- - 2 * sum(extraDistr::dgev(Y, mu = mcmc_mean[1] + X %*% mcmc_mean[3 + 1:p], sigma = mcmc_mean[2], xi = mcmc_mean[3], log = TRUE))
  
  pDIC <- - aux + mean(deviance)
  
  DIC <- aux + 2 * pDIC
  
  return(DIC)
}
