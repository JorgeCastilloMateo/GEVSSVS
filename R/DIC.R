#' DIC
#' 
#' @importFrom extraDistr dgev
#' 
#' @description 
#' This function calculates the DIC of the models retuned by the 
#' \code{\link{GEVmodel}} function.
#' 
#' @param model Object from \code{\link{GEVmodel}}
#' @return DIC value
#' @export 
DIC <- function(model) {
  
  Y <- model$y
  X <- model$x
  
  n.sims <- nrow(model$params)
  p      <- ncol(X)
  
  mcmc_mean <- colMeans(model$params)
  
  deviance <- rep(0, n.sims) 
  
  if (p == 1 && all(X == 0)) {
    for (n in 1:n.sims) {
      deviance[n] <- - 2 * sum(extraDistr::dgev(Y, 
        mu = model$params[n, 1], 
        sigma = model$params[n, 2], 
        xi = model$params[n, 3], log = TRUE))
    }
    
    DhatTheta <- - 2 * sum(extraDistr::dgev(Y, 
      mu = mcmc_mean[1], 
      sigma = mcmc_mean[2], 
      xi = mcmc_mean[3], log = TRUE))
  } else {
    for (n in 1:n.sims) {
      deviance[n] <- - 2 * sum(extraDistr::dgev(Y, 
        mu = model$params[n, 1] + X %*% model$params[n, 3 + 1:p], 
        sigma = model$params[n, 2], 
        xi = model$params[n, 3], log = TRUE))
    }
    
    DhatTheta <- - 2 * sum(extraDistr::dgev(Y, 
        mu = mcmc_mean[1] + X %*% mcmc_mean[3 + 1:p], 
        sigma = mcmc_mean[2], 
        xi = mcmc_mean[3], log = TRUE))
  }
  
  pDIC <- - DhatTheta + mean(deviance)
  
  DIC <- DhatTheta + 2 * pDIC
  
  return(DIC)
}
