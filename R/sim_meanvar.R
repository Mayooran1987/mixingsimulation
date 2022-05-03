##' This function provides the mean and variance of the expected number of CFUs in the single mixing stage.
##' @title Mean and variance of the expected number of CFUs in the single mixing stage.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return Mean and variance changes in the single mixing stage.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}.
##' This function produces a graphical display of the mean and variance changes at each mixing stage. It is helpful to identify the optimal number of revolutions of the mixture, which is a point of mixing that initiates Poisson-like homogeneity.
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha <- 0.1
##' k <- 30
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 2000
##' sim_meanvar(mu, sigma , alpha , k, distribution, n_sim)
##' @export
sim_meanvar <- function(mu, sigma , alpha , k, distribution, n_sim){
  # f_spri <- function(mu, k, alpha, distribution) {
  #   sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
  # }
  i <- NULL
  sim_meanvar_1 <- function(mu, sigma , alpha , k, distribution){
    # set.seed(1, kind = "L'Ecuyer-CMRG")
    rpoislog <- function(S, mu, sig, nu = 1, condS = FALSE, keep0 = FALSE){
      sim <- function(nr) {
        lamx <- rnorm(nr)
        x <- rpois(nr, exp(sig * lamx + mu + log(nu)))
        if (!keep0)
          x <- x[x > 0]
        return(x)
      }
      if (S < 1)
        stop("S is not positive")
      if (!is.finite(S))
        stop("S is not finite")
      if ((S/trunc(S)) != 1)
        stop("S is not an integer")
      if (sig < 0)
        stop("sig is not positive")
      if (nu < 0)
        stop("nu is not positive")
      if (condS) {
        simVec <- vector("numeric", 0)
        fac <- 2
        nr <- S
        while (length(simVec) < S) {
          simvals <- sim(nr * fac)
          simVec <- c(simVec, simvals)
          fac <- (1/(length(simvals)/(nr * fac))) * 2
          fac <- ifelse(is.finite(fac), fac, 1000)
          nr <- S - length(simvals)
        }
        simVec <- simVec[1:S]
      }
      else simVec <- sim(S)
      return(simVec)
    }
    x <-  matrix(NA, nrow = 1, ncol = k) # If we want to apply a beta algorithm to generate Dirichlet distribution's random numbers.
    for (j in 1:k) {
      x[,j] <- stats::rbeta(1,alpha, alpha*(k - j))
    }
    w <-  matrix(NA, nrow = 1, ncol = k)
    for (j in 2:k) {
      w[,1] <- x[,1]
      w[,j] <- x[,j] %*% prod(1 - x[,1:(j - 1)])
    }
    if (distribution == "Poisson-Type A") {
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[,j] <- stats::rpois(1, mu/k)
      }
    } else if (distribution == "Poisson-Type B") {
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[,j] <- stats::rpois(1, mu*w[,j])
      }
    } else if (distribution == "Lognormal-Type A") {
      M <- matrix(round(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = k, nrow = 1)
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[,j] <- stats::rbinom(1, M[,j], 1/k)
      }
    } else if (distribution == "Lognormal-Type B") {
      M <- matrix(round(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = k, nrow = 1)
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[i,j] <- stats::rbinom(1, M[,j], w[,j])
      }
    } else if (distribution == "Poisson lognormal-Type A") {
      M <- matrix(rpoislog( k, log(mu), sigma, keep0 = TRUE), ncol = k, nrow = 1)
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[,j] <- stats::rbinom(1, M[,j], 1/k)
      }
    } else if (distribution == "Poisson lognormal-Type B") {
      M <- matrix(rpoislog( k, log(mu), sigma, keep0 = TRUE), ncol = k, nrow = 1)
      sim <-  matrix(NA, nrow = 1, ncol = k)
      for (j in 1:k) {
        sim[,j] <- stats::rbinom(1, M[,j], w[,j])
      }
    } else {
      print("please choose the one of the given distribution type with case sensitive such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'")
    }
    mean_N <- apply(sim, 1, mean)
    Var_N <- apply(sim, 1, var)
    # result <- data.frame(mean_N,Var_N)
    result <- c(mean_N,Var_N)
    # cat("Calculation took", proc.time()[1], "seconds.\n")
    return(result)
  }
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  result1 <- sapply(1:n_sim, function(i){sim_meanvar_1(mu, sigma , alpha , k, distribution)})
  mean_N <- mean(result1[1,])
  Var_N <- mean(result1[2,])
  results <- as.matrix(c(mean_N,Var_N),nrow = 1,ncol = 2)
  return(results)
}



