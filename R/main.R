#' expected number of failure in SMART
#'
#' @param Q a vector of (Q12, Q21, Q31) for second stage randomization probability in SMART
#' @param pi1 the respond rate of treatment 1, 2, 3 in stage I in SMART
#' @param pi2 the respond rate for patients who do not respond to a treatment in stage I, but
#'            respond to another treatment in stage II
#' @param P Randomization probability in stage I
#' @param n Sample size in SMART
#'
#' @examples
#'
#' Q = c(0.5, 0.5, 0.5)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#'
#' enf(Q, pi1, pi2, P, n)
#'
#' @export
enf <- function(Q,
                pi1 = c(0.4, 0.35, 0.25),
                pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4),
                P = c(1/3, 1/3, 1/3),
                n = 365){

  c0 = sum(P*(1-pi1)*(1-pi2[c(2,4,6)]))
  c1 = sum(Q*P*(1-pi1)*(pi2[c(1,3,5)] - pi2[c(2,4,6)]))

  n*(c0 - c1)
}


#' Dunn-Sidak power for selecting the best regime in SMART
#'
#' @inheritParams enf
#'
#' @return Dunn-Sidak approximated power in SMART
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
DSpower <- function(Q,
                    pi1 = c(0.4, 0.35, 0.25),
                    pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4),
                    P = c(1/3, 1/3, 1/3),
                    n = 365) {

  # Calculate the mean response rate of each DTR: DTRmean
  pi12 <- pi1[[1]] + (1-pi1[[1]])*pi2[[1]]
  pi13 <- pi1[[1]] + (1-pi1[[1]])*pi2[[2]]
  pi21 <- pi1[[2]] + (1-pi1[[2]])*pi2[[3]]
  pi23 <- pi1[[2]] + (1-pi1[[2]])*pi2[[4]]
  pi31 <- pi1[[3]] + (1-pi1[[3]])*pi2[[5]]
  pi32 <- pi1[[3]] + (1-pi1[[3]])*pi2[[6]]

  DTRmean <- cbind(pi12,pi13,pi21,pi23,pi31,pi32)
  colnames(DTRmean) <- c("pi12","pi13","pi21","pi23","pi31","pi32")

  #=================================================#
  # Calculate the variance-covariance matrix: Sigma
  #=================================================#

  # Variances
  sigma12 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[1]])^2/P[[1]] +
    (1-pi1[[1]])*pi2[[1]]*(1-pi2[[1]])/(P[[1]]*Q[[1]])

  sigma13 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[2]])^2/P[[1]] +
    (1-pi1[[1]])*pi2[[2]]*(1-pi2[[2]])/(P[[1]]*(1-Q[[1]]))

  sigma21 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[3]])^2/P[[2]] +
    (1-pi1[[2]])*pi2[[3]]*(1-pi2[[3]])/(P[[2]]*Q[[2]])

  sigma23 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[4]])^2/P[[2]] +
    (1-pi1[[2]])*pi2[[4]]*(1-pi2[[4]])/(P[[2]]*(1-Q[[2]]))

  sigma31 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[5]])^2/P[[3]] +
    (1-pi1[[3]])*pi2[[5]]*(1-pi2[[5]])/(P[[3]]*Q[[3]])

  sigma32 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[6]])^2/P[[3]] +
    (1-pi1[[3]])*pi2[[6]]*(1-pi2[[6]])/(P[[3]]*(1-Q[[3]]))

  # Covariances
  sigma1213 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[1]])*(1-pi2[[2]])/P[[1]]
  sigma2123 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[3]])*(1-pi2[[4]])/P[[2]]
  sigma3132 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[5]])*(1-pi2[[6]])/P[[3]]

  # Output covariance matrix
  Sigma <- matrix(c(sigma12, sigma1213, 0, 0, 0, 0,
                    sigma1213, sigma13, 0, 0, 0, 0,
                    0, 0, sigma21, sigma2123, 0, 0,
                    0, 0, sigma2123, sigma23, 0, 0,
                    0, 0, 0, 0, sigma31, sigma3132,
                    0, 0, 0, 0, sigma3132, sigma32),
                  nrow = 6, ncol = 6, byrow = TRUE)


  # Difference of max DTRmean and other DTRmeans
  D = max(DTRmean) - DTRmean[DTRmean != max(DTRmean)]
  iD = which(DTRmean != max(DTRmean))  ## non-maximum DTRmean id
  iMax <- which.max(DTRmean)           ## maximum DTRmean id

  # Calculate the variance of D
  covD <- sapply(iD, function(x) {
    sum(Sigma[c(iMax, x), c(iMax, x)]*matrix(c(1,-1,-1,1),2))
  })

  # Calculate critical values: cv
  cv <- sqrt(n) * D / sqrt(covD)

  # return the approximate power using Dunn-Sidak inequality
  prod(sapply(cv, function(x) pnorm(x, 0, 1)))

}



#' MC power for selecting the best regime in SMART
#'
#' @inheritParams enf
#' @param m number of Monte Carlo repetitions
#' @param seed seed for generate MC repetitions
#'
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#'
#' seed = 2021
#' Q = c(0.5, 0.5, 0.5)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#' m = 10000
#'
#' MCpower(Q, pi1, pi2, P, n, m, seed)
#'
#' @export
MCpower <- function(Q = c(0.5, 0.5, 0.5),
                    pi1 = c(0.4, 0.35, 0.25),
                    pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4),
                    P = c(1/3, 1/3, 1/3),
                    n = 365,
                    m = 10000,
                    seed = 2021){

  # Calculate the mean response rate of each DTR: DTRmean
  pi12 <- pi1[[1]] + (1-pi1[[1]])*pi2[[1]]
  pi13 <- pi1[[1]] + (1-pi1[[1]])*pi2[[2]]
  pi21 <- pi1[[2]] + (1-pi1[[2]])*pi2[[3]]
  pi23 <- pi1[[2]] + (1-pi1[[2]])*pi2[[4]]
  pi31 <- pi1[[3]] + (1-pi1[[3]])*pi2[[5]]
  pi32 <- pi1[[3]] + (1-pi1[[3]])*pi2[[6]]

  DTRmean <- cbind(pi12,pi13,pi21,pi23,pi31,pi32)
  colnames(DTRmean) <- c("pi12","pi13","pi21","pi23","pi31","pi32")

  #=================================================#
  # Calculate the variance-covariance matrix: Sigma
  #=================================================#

  # Variances
  sigma12 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[1]])^2/P[[1]] +
    (1-pi1[[1]])*pi2[[1]]*(1-pi2[[1]])/(P[[1]]*Q[[1]])

  sigma13 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[2]])^2/P[[1]] +
    (1-pi1[[1]])*pi2[[2]]*(1-pi2[[2]])/(P[[1]]*(1-Q[[1]]))

  sigma21 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[3]])^2/P[[2]] +
    (1-pi1[[2]])*pi2[[3]]*(1-pi2[[3]])/(P[[2]]*Q[[2]])

  sigma23 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[4]])^2/P[[2]] +
    (1-pi1[[2]])*pi2[[4]]*(1-pi2[[4]])/(P[[2]]*(1-Q[[2]]))

  sigma31 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[5]])^2/P[[3]] +
    (1-pi1[[3]])*pi2[[5]]*(1-pi2[[5]])/(P[[3]]*Q[[3]])

  sigma32 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[6]])^2/P[[3]] +
    (1-pi1[[3]])*pi2[[6]]*(1-pi2[[6]])/(P[[3]]*(1-Q[[3]]))

  # Covariances
  sigma1213 <- pi1[[1]]*(1-pi1[[1]])*(1-pi2[[1]])*(1-pi2[[2]])/P[[1]]
  sigma2123 <- pi1[[2]]*(1-pi1[[2]])*(1-pi2[[3]])*(1-pi2[[4]])/P[[2]]
  sigma3132 <- pi1[[3]]*(1-pi1[[3]])*(1-pi2[[5]])*(1-pi2[[6]])/P[[3]]

  # Output covariance matrix
  Sigma <- matrix(c(sigma12, sigma1213, 0, 0, 0, 0,
                    sigma1213, sigma13, 0, 0, 0, 0,
                    0, 0, sigma21, sigma2123, 0, 0,
                    0, 0, sigma2123, sigma23, 0, 0,
                    0, 0, 0, 0, sigma31, sigma3132,
                    0, 0, 0, 0, sigma3132, sigma32),
                  nrow = 6, ncol = 6, byrow = TRUE)



  # Difference of max DTRmean and other DTRmeans
  D = max(DTRmean) - DTRmean[DTRmean != max(DTRmean)]
  iD = which(DTRmean != max(DTRmean))  ## non-maximum DTRmean id
  iMax <- which.max(DTRmean)           ## maximum DTRmean id

  # Calculate the MVN covariance matrix
  M <- length(iD)
  sim_Sigma <- matrix(nrow = M, ncol = M)
  cvcov <- c()

  for(i in iD){

    #Save covariance of iMax and i to calculate critical value
    covMi <- matrix(c(1, -1), nrow=1) %*% Sigma[c(iMax, i), c(iMax, i)] %*% matrix(c(1, -1), ncol=1)
    cvcov <- c(cvcov, covMi)

    for(j in iD){
      covD <- matrix(c(1, -1), nrow=1) %*% Sigma[c(iMax, i), c(iMax, j)] %*% matrix(c(1, -1), ncol=1)
      covMj <- matrix(c(1, -1), nrow=1) %*% Sigma[c(iMax, j), c(iMax, j)] %*% matrix(c(1, -1), ncol=1)
      sim_Sigma[which(iD == i), which(iD == j)] <- covD / (sqrt(covMi) * sqrt(covMj))
    }
  }


  # Generate m Monte Carlo replicates from MVN distribution
  set.seed(seed)
  simMVN <- rmvnorm(m, rep(0, M), sim_Sigma)


  # Calculate critical values
  cv <- D * sqrt(n) / sqrt(cvcov)


  # Get empirical power
  sum(apply(simMVN, 1, function(x) sum(x < cv)==M))/m

}



#' Gradient of log Dunn-Sidak power w.r.t
#'
#' @inheritParams enf
#'
#' @return a vector of gradients w.r.t the second stage randomization probability for the
#'         Dunn-Sidak approximated power in RASMART
#'
#' @importFrom Matrix bdiag
#'
#' @examples
#'
#' Q = c(0.5, 0.5, 0.5)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#'
#' DSpowerl_grad(Q, pi1, pi2, P, n)
#'
#' @export
DSpowerl_grad <- function(Q,
                          pi1 = c(0.4, 0.35, 0.25),
                          pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4),
                          P = c(1/3, 1/3, 1/3),
                          n = 365){

  # Calculate the mean response rate of each DTR: DTRmean
  pi12 <- pi1[[1]] + (1-pi1[[1]])*pi2[[1]]
  pi13 <- pi1[[1]] + (1-pi1[[1]])*pi2[[2]]
  pi21 <- pi1[[2]] + (1-pi1[[2]])*pi2[[3]]
  pi23 <- pi1[[2]] + (1-pi1[[2]])*pi2[[4]]
  pi31 <- pi1[[3]] + (1-pi1[[3]])*pi2[[5]]
  pi32 <- pi1[[3]] + (1-pi1[[3]])*pi2[[6]]

  DTRmean <- cbind(pi12,pi13,pi21,pi23,pi31,pi32)
  colnames(DTRmean) <- c("pi12","pi13","pi21","pi23","pi31","pi32")


  #================================================#
  # REORDER THE PARAMETERS IN THE FOLLOWING ORDER: #
  #       MAX DTR, MAX COMPLEMENT DTR, OTHERS      #
  # E.G.: IF A2A3 IS THE MAX, oP=c(P2, P1, P3)     #
  #       oQ=c(Q21, Q12, Q31)                      #
  #       opi1=c(pi1_A2, pi1_A1, pi1_A3)           #
  #       opi2=c(pi2_A2A3, pi2_A2A1, pi2_A1A2, pi2_A1A3, pi2_A3A1, pi2_A3A2) #
  #       oMean=c(pi23, pi21, pi12, pi13, pi31, pi32) #
  #================================================#

  iMax <- which.max(DTRmean)  ## maximum DTRmean id
  iMaxc <- iMax-1 + 2*(iMax %% 2) ## maximum DTRmean complement DTR id
  oMax3 <- ceiling(iMax/2) # Order of max regime in sequence of three (P, Q, and Response rate at stage I)
  oMax3c <- seq(3)[seq(3) != oMax3] # Orders of regimes other than max in sequence of three
  oMax6c <- seq(6)[!seq(6) %in% c(iMax, iMaxc)] # Orders of regimes other than max and max complementary in sequence of six

  # Reorder the randomization probabilities at stage I: oP
  oP <- P[c(oMax3, oMax3c)]
  names(P) <- c("P1", "P2", "P3")
  names(oP) <- names(P)[c(oMax3, oMax3c)]

  # Reorder the randomization probabilities at stage II: oQ
  QMax <- ifelse(iMax %% 2 == 0, 1-Q[oMax3], Q[oMax3]) # max regime corresponding Q
  oQ <- c(QMax, Q[oMax3c])
  names(Q) <- c("Q12", "Q21", "Q31")
  QMaxnm <- ifelse(iMax %% 2 == 0, paste0("1-", names(Q[oMax3])), names(Q[oMax3]))
  names(oQ) <- c(QMaxnm, names(Q[oMax3c]))

  # Reorder the response rate at the end of stage I: opi1
  opi1 <- pi1[c(oMax3, oMax3c)]
  names(opi1) <- dimnames(pi1)[[2]][c(oMax3, oMax3c)]

  # Reorder the subgroup means at stage II: opi2
  opi2 <- pi2[c(iMax, iMaxc, oMax6c)]
  names(opi2) <- dimnames(pi2)[[2]][c(iMax, iMaxc, oMax6c)]

  # Reorder mean vector
  oMean <- DTRmean[c(iMax, iMaxc, oMax6c)]
  names(oMean) <- dimnames(DTRmean)[[2]][c(c(iMax, iMaxc, oMax6c))]


  #========================================================================#
  # Calculate the variance-covariance matrix for reordered vectors: oSigma #
  #========================================================================#

  # Variances
  osigma12 <- opi1[[1]]*(1-opi1[[1]])*(1-opi2[[1]])^2/oP[[1]] +
    (1-opi1[[1]])*opi2[[1]]*(1-opi2[[1]])/(oP[[1]]*oQ[[1]])

  osigma13 <- opi1[[1]]*(1-opi1[[1]])*(1-opi2[[2]])^2/oP[[1]] +
    (1-opi1[[1]])*opi2[[2]]*(1-opi2[[2]])/(oP[[1]]*(1-oQ[[1]]))

  osigma21 <- opi1[[2]]*(1-opi1[[2]])*(1-opi2[[3]])^2/oP[[2]] +
    (1-opi1[[2]])*opi2[[3]]*(1-opi2[[3]])/(oP[[2]]*oQ[[2]])

  osigma23 <- opi1[[2]]*(1-opi1[[2]])*(1-opi2[[4]])^2/oP[[2]] +
    (1-opi1[[2]])*opi2[[4]]*(1-opi2[[4]])/(oP[[2]]*(1-oQ[[2]]))

  osigma31 <- opi1[[3]]*(1-opi1[[3]])*(1-opi2[[5]])^2/oP[[3]] +
    (1-opi1[[3]])*opi2[[5]]*(1-opi2[[5]])/(oP[[3]]*oQ[[3]])

  osigma32 <- opi1[[3]]*(1-opi1[[3]])*(1-opi2[[6]])^2/oP[[3]] +
    (1-opi1[[3]])*opi2[[6]]*(1-opi2[[6]])/(oP[[3]]*(1-oQ[[3]]))

  # Covariances
  osigma1213 <- opi1[[1]]*(1-opi1[[1]])*(1-opi2[[1]])*(1-opi2[[2]])/oP[[1]]
  osigma2123 <- opi1[[2]]*(1-opi1[[2]])*(1-opi2[[3]])*(1-opi2[[4]])/oP[[2]]
  osigma3132 <- opi1[[3]]*(1-opi1[[3]])*(1-opi2[[5]])*(1-opi2[[6]])/oP[[3]]

  # Output covariance matrix
  oSigma = as.matrix(bdiag(list(matrix(c(osigma12, osigma1213, osigma1213, osigma13),2),
                                matrix(c(osigma21, osigma2123, osigma2123, osigma23),2),
                                matrix(c(osigma31, osigma3132, osigma3132, osigma32),2))
  ))


  # Calculate the critival values of the DTRmean difference from the optimal DTR
  D <- oMean[1] - oMean[-1] ## Difference of max DTRmean and other DTRmeans
  covD <- sapply(1:length(D), function(x) { ## Variance of D
    t(c(1, -1)) %*% oSigma[c(1, x+1), c(1, x+1)] %*% c(1, -1)
  })
  cv <- sqrt(n) * D / sqrt(covD) ## Critical values


  # Calculate the constants C_sigma^2 before Qs: cS2
  cSA <- matrix(c(1,0,0,
                  1,0,0,
                  0,1,0,
                  0,1,0,
                  0,0,1,
                  0,0,1),
                nrow = 6, ncol = 3, byrow = TRUE)

  cS2 <- cSA %*% ((1/oP)*(1-opi1))*opi2*(1-opi2)
  rownames(cS2) <- c("cS2_12", "cS2_13", "cS2_21", "cS2_23", "cS2_31", "cS2_32")


  #==================================================#
  # Calculate terms in derivatives of power w.r.t Qs #
  #==================================================#

  # Term 1: CDF * PDF at critical values
  cdfs <- as.vector(sapply(cv, function(x) pnorm(x, 0, 1)))
  pdfs <- as.vector(sapply(cv, function(x) dnorm(x, 0, 1)))
  term1 <- 1/cdfs * pdfs

  # Term 2: sqrt(1/covD)
  term2 <- 1/sqrt(covD)

  # Term 3: (Critical value * 1/covD)/2
  term3 <- (cv * 1/covD)/2

  # Term 4: constant (cS2) * Qs
  exQA <- matrix(c(1,0,0,
                   -1,0,0,
                   0,1,0,
                   0,-1,0,
                   0,0,1,
                   0,0,-1),
                 nrow = 6, ncol = 3, byrow = TRUE)
  exQs <- (c(0,1,0,1,0,1) + exQA %*% oQ)^(-2)
  term4 <- as.vector(cS2) * as.vector(exQs)

  # Derivatives of power w.r.t. individual Q in the order of (Q for iMax, rest two Qs)
  dPQ1A <- matrix(c(1,-1,0,0,0,0,
                    1,0,0,0,0,0,
                    1,0,0,0,0,0,
                    1,0,0,0,0,0,
                    1,0,0,0,0,0),
                  nrow = 5, ncol = 6, byrow = TRUE)
  dPQ1 <- term1 %*% (term2 + term3 * (dPQ1A %*% term4))
  dPQ2 <- term1[c(2,3)] %*% (term2[c(2,3)] + term3[c(2,3)] * (c(1,-1) * term4[c(3, 4)]))
  dPQ3 <- term1[c(4,5)] %*% (term2[c(4,5)] + term3[c(4,5)] * (c(1,-1) * term4[c(5, 6)]))


  # Return gradients of power w.r.t Qs
  GradQ <- c(dPQ1, dPQ2, dPQ3)
  names(GradQ) <- names(oQ)

  GradQ
}


#' minimum failure count in RASMART subject to a lower bound of the power
#'
#' @param x initial second stage randomization probability in SMART
#' @param pi1 the respond rate of treatment 1, 2, 3 in stage I in SMART
#' @param pi2 the respond rate for patients who do not respond to a treatment in stage I, but
#'            respond to another treatment in stage II
#' @param P Randomization probability in stage I
#' @param n Sample size in SMART
#' @param power the minimum controlled power for RASMART
#' @param alpha.step The step size to search a small alpha for satisfying the Lipschitz constant
#'            such that f(x) <= f(x_0) + df(x_0)(x-x_0) + 1/(2alpha)||x-x_0||^2. Default is 0.8.
#' @param lambda_seq The penalty parameter to minimize the lagrangian function.
#' @param ... other arguments related to mfc
#'
#' @import dplyr
#'
#' @examples
#'
#' x = c(0.3, 0.3, 0.3)
#' pi1 = c(0.4, 0.35, 0.25)
#' pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4)
#' P = c(1/3, 1/3, 1/3)
#' n = 365
#' power = 0.8
#' alpha.step = 0.8
#' lambda_seq = seq(1,50, by = 1)
#'
#' mfc(x, pi1, pi2, P, n, power, alpha.step, lambda_seq)
#'
#' @export
mfc.default <- function(x = c(0.5, 0.5, 0.5),
                        pi1 = c(0.4, 0.35, 0.25),
                        pi2 = c(0.35, 0.32, 0.36, 0.56, 0.37, 0.4),
                        P = c(1/3, 1/3, 1/3),
                        n = 365,
                        power = 0.8,
                        alpha.step = 0.8,
                        lambda_seq = seq(1, 100, by = 0.5),
                        ...){
  ## initial Q
  Q = x
  optimal_seq = Lambda_seq = NULL
  # i = 0
  for (lambda in lambda_seq) {
    ## initialize Q
    old_Q = Q
    outer_itr = 0
    path = NULL ## store Q, lambda, y
    # alpha.step = 0.1
    # print(paste("i = ", i, " is running!"))
    # i = i + 1

    repeat{

      old_L = enf(old_Q, pi1, pi2, P, n) + lambda*(log(power) - log(DSpower(old_Q, pi1, pi2, P, n)))
      grad_Q = -n*P*(1-pi1)*(pi2[c(1,3,5)] - pi2[c(2,4,6)]) -
        lambda*DSpowerl_grad(old_Q, pi1, pi2, P, n)

      ## update Q and lambda
      ## use Majorization-minimization approximation to backtrack line search alpha
      inner_itr = 0
      alpha = 1
      new_Q = new_L = NULL
      check_region = NULL ## check if Q is within [0,1]

      repeat{

        new_Q = old_Q - alpha*grad_Q

        ## check if new_Q is within (0,1)
        check_region = all(new_Q > 0 & new_Q < 1)
        if(!check_region){
          alpha = alpha.step*alpha
          inner_itr = inner_itr + 1
          next
        }

        new_L = enf(new_Q, pi1, pi2, P, n) + lambda*(log(power) - log(DSpower(new_Q, pi1, pi2, P, n)))
        L_approx = old_L + grad_Q%*%c(new_Q - old_Q) + 0.5/alpha*sum(c(new_Q - old_Q)^2)

        if( (new_L <= L_approx) | inner_itr>50){
          break
        } else{
          alpha = alpha.step*alpha
          inner_itr = inner_itr + 1
        }
      }

      eps = norm(new_Q - old_Q, type = "2")/sqrt(3)

      if(eps <= 1e-6 | outer_itr > 200){
        break
      } else{
        old_Q = new_Q
        outer_itr = outer_itr + 1
        path = rbind(path, c(new_Q, enf(new_Q, pi1, pi2, P, n), new_L))
      }

    }

    optimal_seq = rbind(optimal_seq, old_Q)
    Lambda_seq = c(Lambda_seq, lambda)

    if( mean(abs(old_Q - Q)) < 1e-9) {
      break
    }
  }

  opt_power = apply(optimal_seq, 1, function(t) DSpower(t, pi1, pi2, P, n))
  ENF = apply(optimal_seq, 1, function(t) enf(t, pi1, pi2, P, n))

  AllTab = data.frame(cbind(optimal_seq, Lambda_seq, opt_power, ENF))
  colnames(AllTab) = c(names(grad_Q), "Lambda", "Power", "ENF")

  ActiveTab = AllTab %>%
    dplyr::filter(Power >= power) %>%
    dplyr::arrange(desc(Power))

  optimalQ = ActiveTab %>%
    dplyr::slice(which.min(ENF)) %>%
    dplyr::select(contains("Q"))

  minenf = ActiveTab %>%
    dplyr::slice(which.min(ENF)) %>%
    dplyr::pull(ENF)

  DSpower = ActiveTab %>%
    dplyr::slice(which.min(ENF)) %>%
    dplyr::pull(Power)

  result = list(optimalQ = as.numeric(optimalQ),
                minENF = as.numeric(minenf),
                DSpower = as.numeric(DSpower),
                ActiveTab = ActiveTab,
                AllTab = AllTab)
  names(result$optimalQ) = names(grad_Q)

  result$call <- match.call()
  class(result) <- "mfc"

  result
}


#' mfc generic
#'
#' @param x an object of class from ''mfc''.
#' @param ... further arguments passed to or from other methods.
#'
#' @keywords internal
#'
#' @seealso \code{\link{mfc.default}}
#'
#'
#' @export
mfc <- function(x, ...) UseMethod("mfc")


#' @rdname mfc
#'
#' @keywords internal
#'
#' @export
print.mfc <- function(x, ...){
  cat("Call:\n")
  print(x$call)

  cat("\nOptimal Q:\n ")
  print(x$optimalQ)
  #cat(x$optimalQ, fill = 2, labels = names(x$optimalQ))

  cat("\nMinimum ENF: ")
  cat(x$minENF)

  cat("\nDunn-Sidak power: ")
  cat(x$DSpower)

}


