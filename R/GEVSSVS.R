#' GEV-SSVS regression
#' 
#' Generalized Extreme Value (GEV) regression with Stochastic Search Variable 
#' Selection (SSVS).
#' 
#' @references 
#' Castillo-Mateo J, Asín J, Cebrián AC, Mateo-Lázaro J, Abaurrea J (2023). 
#' “Bayesian variable selection in generalized extreme value regression: Modeling annual maximum temperature.” 
#' \emph{Mathematics}, \strong{11}(3), 759. \doi{10.3390/math11030759}.
#' 
#' @docType package
#' @author Jorge Castillo-mateo <jorgecastillomateo@gmail.com>
#' @import extraDistr Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib GEVSSVS
#' @name GEVSSVS
NULL  