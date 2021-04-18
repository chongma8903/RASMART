#' 2D or 3D visualization for the Dunn-Sidak power
#'
#' @param Q1 a vector of numeric value for the range of Q12 probability
#' @param Q2 a vector of numeric value for the range of Q21 probability
#' @param Q3 a vector of numeric value for the range of Q31 probability
#' @param pi1 the respond rate of treatment 1, 2, 3 in stage I in SMART
#' @param pi2 the respond rate for patients who do not respond to a treatment in stage I, but
#'            respond to another treatment in stage II
#' @param P Randomization probability in stage I
#' @param n Sample size in SMART
#' @param ... other arguments fits to DSpower and plot3D
#'
#' @return a 2D contour plot or 3D plot for the Dunn-Sidak power
#'
#' @details One of Q1, Q2, Q3 must be a single value, and the other two should be a range of
#'          numerical value between 0 and 1.
#'
#' @importFrom stats  dnorm pnorm
#'
#' @examples
#'
#' Q1 = 0.3
#' Q2 = seq(0.1, 0.9, by = 0.005)
#' Q3 = seq(0.1, 0.9, by = 0.005)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#'
#' plot_DSpower(Q1, Q2, Q3, pi1, pi2, P, n, col = "black", nlevels = 30)
#'
#' @export
plot_DSpower <- function(Q1, Q2, Q3, pi1, pi2, P, n, ...) {

  len = c(length(Q1), length(Q2), length(Q3))
  if(sum(len == 1) != 1) stop("Only one of Q1, Q2, Q3 be a length of one numeric value!")

  if(length(Q1) == 1){
    len1 = length(Q2)
    len2 = length(Q3)
    z = matrix(NA, nrow = len1, ncol = len2)

    for(i in 1:len1){
      for(j in 1:len2){
        z[i, j] = DSpower(Q = c(Q1, Q2[i], Q3[j]), pi1, pi2, P, n)
      }
    }

   contour2D(z, x = Q2, y = Q3, ...)
  }

  if(length(Q2) == 1){
    len1 = length(Q1)
    len2 = length(Q3)
    z = matrix(NA, nrow = len1, ncol = len2)

    for(i in 1:len1){
      for(j in 1:len2){
        z[i, j] = DSpower(Q = c(Q1[i], Q2, Q3[j]), pi1, pi2, P, n)
      }
    }

    contour2D(z, x = Q1, y = Q3, ...)

  }

  if(length(Q3) == 1){
    len1 = length(Q1)
    len2 = length(Q2)
    z = matrix(NA, nrow = len1, ncol = len2)

    for(i in 1:len1){
      for(j in 1:len2){
        z[i, j] = DSpower(Q = c(Q1[i], Q2[i], Q3), pi1, pi2, P, n)
      }
    }

    contour2D(z, x = Q1, y = Q2, ...)

  }

}
























