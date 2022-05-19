##' This function calculates the resulting generated number of colony forming units in the mixed sample in the single mixing plan with single stage of the mixing.
##' @title The generated number of colony-forming units in the mixed sample by the simulation results in the single mixing plan with a single stage of the mixing.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param summary if we need to get all simulated \eqn{N'}, use \code{summary = 2}; otherwise, if we use \code{summary = 1}, the function provides the mean value of the simulated \eqn{N'}.
##' @param n_sim number of simulations
##' @return total number of colony forming units in the single mixing plan
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}. To more details, please refer the details section of  \link{compare_mixing_3}.
##' @seealso  \link{compare_mixing_3}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha <- 0.1
##' k <- 30
##' n_sim <- 20000
##' sim_single(mu, sigma, alpha, k, distribution = "Poisson lognormal-Type B", summary = 1, n_sim)
##' @export
sim_single <- function(mu, sigma , alpha , k, distribution, summary, n_sim){
  i <- NULL
  sim_single_1 <- function(mu, sigma , alpha , k, distribution){
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
      sim1 <-  matrix(NA, nrow = 2000, ncol = k)
      # to get a precise estimation, we used 2000 generated binomial random variables.
      for (j in 1:k) {
        sim1[,j] <- stats::rbinom(2000, M[,j], w[,j])
      }
      sim <- round(apply(sim1, 2, mean))
    } else {
      print("please choose the one of the given distribution type with case sensitive such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'")
    }
    result <- sum(sim)

    return(result)
  }
  if (summary == 1) {
    results <- round(mean(as.integer(lapply(1:n_sim, function(i){sim_single_1(mu, sigma , alpha , k, distribution)}))))
} else if (summary == 2) {
  results <- as.integer(lapply(1:n_sim, function(i){sim_single_1(mu, sigma , alpha , k, distribution)}))
} else {
  print("please include the summary value, which depends on your expected output")
}
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(results)
}

