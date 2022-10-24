#' GEV model with SSVS
#' 
#' @importFrom stats cov
#' @importFrom stats dnorm
#' @importFrom stats rbinom
#' 
#' @description 
#' This function fits the GEV model 
#' 
#' @param Y VECTOR of data
#' @param X MATRIX of covariates (design matrix without intercept)
#' @param inits Initial values (\eqn{\sigma} transformed scale \eqn{(-\infty,\infty)})
#' @param prior Vector. Variance of the zero-mean normal prior for 
#'   (1) \eqn{\beta_0}, (2) \eqn{\log \sigma}, (3) \eqn{\xi}, and value of
#'   (4) \eqn{\tau_i^2}, (5) \eqn{c_i^2 \tau_i^2}, and (6) \eqn{p_i}
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations. (ii) 
#'   Thinning rate. (iii) Number of iterations discarded at the beginning. (iv)
#'   Report the number of iterations rate.
#' @return A \code{"GEVmodel"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{\beta_0,\sigma,\xi,\beta_1,\ldots,\beta_p,\gamma_1,\ldots,\gamma_p}}
#'   \item{\code{y}}{Data fitted}
#'   \item{\code{x}}{Covariates}
#' @export 
GEVmodel <- function(Y, 
                     X = NULL, 
                     inits = NULL, 
                     prior = c(100^2, 100^2, 3, 1 / 100^2, 100^2, 0.5),
                     n.sims = 100000,
                     n.thin = 1,
                     n.burnin = 10000,
                     n.report = 10000
                     ) {
  
  T <- length(Y)
  p <- ifelse(is.null(X), 0, ncol(X))
  d <- 3 + 2 * p
  epsilon <- 0.000001
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d-p)
  keep   <- matrix(nrow = n.sims / n.thin, ncol = d)
  if (p == 0) {
    colnames(keep) <- c("mu", "sigma", "xi")
    X <- matrix(0, nrow = T, ncol = 1)
  } else {
    colnames(keep) <- c("mu", "sigma", "xi", paste0("beta", 1:p), paste0("gamma", 1:p))
  }
  if (is.null(inits)) {
    params <- rep(0, d)
    if (p != 0) params[(4 + p):(3 + 2 * p)][((1:p) + 1) %% 2 == 0] <- 1
  } else if (d == length(inits)) {
    params <- inits
  } else {
    stop("'inits' does not have the proper length")
  }
  if (params[3] == 0) params[3] <- 0.001
  params <- c(params, logfall(T, p, params, Y, X, prior))
  I <- epsilon * diag(d - p)
  pi <- prior[6]
  taui <- sqrt(prior[4])
  ci_x_taui <- sqrt(prior[5])

  # for 1
  for (b in 1:25) {
    params <- rwBmetropolis(params, I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burnin
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[1:(3+p)]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[10000 + (1:n.burnin), ]),
    stats::cov(keepBurnin[10000 + (1:n.burnin), ])
  )
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolis(params, Sigma[[2]] + I, T, p, Y, X, prior)
    if (p != 0) {
      for (i in 1:p) {
        a <- pi * stats::dnorm(params[3 + i], mean = 0, sd = ci_x_taui)
        u <- (1 - pi) * stats::dnorm(params[3 + i], mean = 0, sd = taui)
        params[3 + p + i] <- stats::rbinom(1, size = 1, prob = a / (a + u))
      }
    }
    
    Sigma <- MuSigmaUpdate(params[1:(3+p)], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-(d+1)]
    }
    
  }
  
  keep[,2] <- exp(keep[,2])
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y
  keep.list$x <- X
  
  class(keep.list) <- "GEVmodel"
  
  return(keep.list)
}
