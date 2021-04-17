#' 2D or 3D visualization for the Dunn-Sidak power
#'
#' @param Q1 a vector of numeric value for the range of Q12 probability
#' @param Q2 a vector of numeric value for the range of Q12 probability
#' @param Q3 a vector of numeric value for the range of Q12 probability
#'
#' @return a 2D contour plot or 3D plot for the Dunn-Sidak power
#'
#' @importFrom stats  dnorm pnorm
#'
#' @examples
#'
#' Q = c(0.5, 0.5, 0.5)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#'
#' DSpower(Q, pi1, pi2, P, n)
#'
#' @export
