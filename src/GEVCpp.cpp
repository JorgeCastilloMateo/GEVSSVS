// Functions for GEV SSVS

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

/////////
// GEV //
/////////

// [[Rcpp::export]]
double logdgev(double x, double mu, double sigma, double xi) {
  
  double aux = 1 + xi * ((x - mu) / sigma);
  if (aux < 0) {
    return -INFINITY;
  } else {
    double density = - log(sigma) - (1 + 1 / xi) * log(aux) - pow(aux, - 1 / xi);
    return density;
  }
}

/////////////////
//    OTHER    //
/////////////////

// dnorm ~ log mu = 0
//
// Proportional in x to normal log-density with mu = 0
//
// @param x parameter
// @param var variance
// @return log-density
// [[Rcpp::export]]
double logdnorm(double x, double var) {
  return (- pow(x, 2) / (2 * var));
}

// Online Update of Mean and Variance
//
// Update mu and Sigma for tunning the Metropolis step
//
// @param x Vector New value
// @param mu Vector Old mean
// @param Sigma Matrix Old variance
// @param n Number of data with x
// @return List with updated mu and Sigma
// [[Rcpp::export]]
Rcpp::List MuSigmaUpdate(arma::vec x, arma::vec mu, arma::mat Sigma, int n) {
  arma::vec u = (x - mu) / n;
  mu += u;
  Sigma = (n - 2) * Sigma / (n - 1) + u * u.t() * n;
  return Rcpp::List::create(mu, Sigma);
}

/////////////////
//  GEV MODEL  //
/////////////////

// [[Rcpp::export]]
double logfall(int T, int p, arma::vec x, arma::vec Y, arma::mat X, arma::vec prior) {
  
  double logLikelihood = logdnorm(x(0), prior(0)) + 
                         logdnorm(x(1), prior(1)) +
                         logdnorm(x(2), prior(2));
  
  if (p == 0) {
    for (int t = 0; t < T; ++t) {
      logLikelihood += logdgev(Y(t), x(0), exp(x(1)), x(2));
    }
  } else {
    for (int t = 0; t < T; ++t) {
      logLikelihood += logdgev(Y(t), x(0) + (X.row(t) * x(arma::span(3, 2 + p))).eval()(0), exp(x(1)), x(2));
    }
    for (int i = 0; i < p; ++i) {
      logLikelihood += logdnorm(x(3 + i), (1 - x(3 + p + i)) * prior(3) + x(3 + p + i) * prior(4));
    }
  }
  
  return logLikelihood;
}

// [[Rcpp::export]]
arma::vec rwBmetropolis(arma::vec x, arma::mat sd, int T, int p, arma::vec Y, arma::mat X, arma::vec prior) {
  arma::vec y = x; 
  y(arma::span(0, 2 + p)) = x(arma::span(0, 2 + p)) + arma::chol(sd, "lower") * arma::randn(3 + p);
  y(3 + 2 * p) = logfall(T, p, y(arma::span(0, 2 + 2 * p)), Y, X, prior);
  double A = y(3 + 2 * p) - x(3 + 2 * p);
  if (log(R::runif(0, 1)) <= A) x = y;
  return x;
}
