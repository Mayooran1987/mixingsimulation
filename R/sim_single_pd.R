##' This function gives a probability of detection at a single stage (or revolution) of the mixing process.
##' @title The estimated value of detection probability at a single stage (or revolution) of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return The probability of detection at each stage of the mixing process.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and let \eqn{l} be the number of stages in the mixing process.
##' This function provides the probability of detection at each stage of the mixing process. The probability of detection can be determined by how many primary samples contain CFUs greater than UDL out of the number of primary samples engaged at each mixing stage.
##'
##' Therefore, the probability of detection (\eqn{p_d}) can be estimated from following formula,
##'
##' \deqn{p_d = \frac{\textnormal{Number of primary samples which contain CFUs greater than UDL}}{\textnormal{Number of primary samples}} ;}
##'
##' where the upper decision limit (UDL) depends on microorganisms and testing regulations. For example, UDL should be equal to 0 for testing Salmonella in milk powder sample if we consider 25g primary sample.
##'
##' @seealso \link{sim_single_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha <- 0.1
##' k <- 30
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 0
##' n_sim <- 2000
##' sim_single_pd(mu, sigma , alpha , k, distribution, UDL, n_sim)
##' @export
sim_single_pd <- function(mu, sigma , alpha , k, distribution, UDL, n_sim){
  i <- NULL
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
  }
  sim_single_pd <- function(mu, sigma , alpha , k, distribution){
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
    result <- length(which(sim > UDL))/k
    # cat("Calculation took", proc.time()[1], "seconds.\n")
    return(result)
  }
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  results <- mean(as.numeric(lapply(1:n_sim, function(i){sim_single_pd(mu, sigma , alpha , k, distribution)})))
  # colnames(results) <- f_spri(mu, k, alpha, distribution)
  return(results)
}




