##' This function calculates the resulting generated number of colony forming units in the mixed sample in the single mixing plan with single stage of the mixing.
##' @title The generated number of colony-forming units in the mixed sample by the simulation results in the single mixing plan with a single stage of the mixing.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param summary need to select one number from the list 1 to 5 which depends on what simulated observations are needed (default 1). If \code{summary = 1}, the function provides the expected total number of the simulated \eqn{N'} after \code{n_sim} simulations are run.
##' If \code{summary = 2}, the function provides generated CFUs in each primary sample after \code{n_sim} simulations are run. If \code{summary = 3}, the function provides the expected total number of the simulated \eqn{N'} at each simulation run.
##' If \code{summary = 4}, the function provides the mean and variance values of the simulated samples. If \code{summary = 5}, the function provides all simulated observations in the \code{n_sim} * \code{k} matrix.
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
##' sim_single(mu, sigma, alpha, k, distribution = "Poisson lognormal-Type B", n_sim, summary = 2)
##' @export
sim_single <- function(mu, sigma , alpha , k, distribution, n_sim, summary = 1){
  # #If we want to apply a gamma algorithm to generate Dirichlet distribution's random numbers.
  # #x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
  # #sm <- x%*%rep(1, k)
  # #w <- x/as.vector(sm)
  # set.seed(seed, kind = "L'Ecuyer-CMRG")
  x <-  matrix(NA, nrow = 1, ncol = k) # If we want to apply a beta algorithm to generate Dirichlet distribution's random numbers.
  for (j in 1:k) {
    x[,j] <- stats::rbeta(1,alpha, alpha*(k - j))
  }
  w <-  matrix(NA, nrow = 1, ncol = k)
  for (j in 2:k) {
    w[,1] <- x[,1]
    w[,j] <- x[,j] %*% prod(1 - x[,1:(j - 1)])
  }
  #x
  # sum(w)
  # stick_breaking_process <-  function(num_weights, alpha) {
  #   betas <-  stats::rbeta(num_weights, 1, alpha)
  #   remaining_stick_lengths <-  c(1, cumprod(1 - betas))[1:num_weights]
  #   weights <-  remaining_stick_lengths * betas
  #   return(weights)
  # }
  # num_weights <- 10
  # alpha <- 0.1
  # w <-  as.matrix(stick_breaking_process(k, alpha))
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
  if (distribution == "Poisson-Type A") {
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rpois(1, mu/k)
      }
    }
  } else if (distribution == "Poisson-Type B") {
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rpois(1, mu*w[,j])
      }
    }
  } else if (distribution == "Lognormal-Type A") {
    M <- matrix(NA, ncol = k, nrow = 1)
    for (j in 1:k) {
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      set.seed(1, kind = "L'Ecuyer-CMRG")
      M[,j] <- as.integer(stats::rlnorm(1, meanlog = log(mu), sdlog = sigma)) # lognormal distribution
    }
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rbinom(1, M[,j], 1/k)
      }
    }
  } else if (distribution == "Lognormal-Type B") {
    M <- matrix(NA, ncol = k, nrow = 1)
    for (j in 1:k) {
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      set.seed(1, kind = "L'Ecuyer-CMRG")
      M[,j] <- as.integer(stats::rlnorm(1, meanlog = log(mu), sdlog = sigma)) # lognormal distribution
    }
    # M <- matrix(as.integer(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rbinom(1, M[,j], w[,j])
      }
    }

  } else if (distribution == "Poisson lognormal-Type A") {
    # code "rpoislog" used from R package "poilog"
    # M <- matrix(NA, ncol = k, nrow = 1)
    # for (j in 1:k) {
    #   # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
    #   M[,j] <- as.integer(VGAM::rpolono(1, meanlog = log(mu), sdlog = sigma)) # poisson lognormal distribution
    #
    # }
    set.seed(1, kind = "L'Ecuyer-CMRG")
    M <- matrix(rpoislog( k, log(mu), sigma, keep0 = TRUE), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rbinom(1, M[,j], 1/k)
      }
    }
  } else if (distribution == "Poisson lognormal-Type B") {
    # M <- matrix(NA, ncol = k, nrow = 1)
    # for (j in 1:k) {
    #   # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
    #   M[,j] <- as.integer(VGAM::rpolono(1, meanlog = log(mu), sdlog = sigma)) # poisson lognormal distribution
    # }
    # M <- matrix(as.integer(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = 1, nrow = k)
    # M <- matrix(VGAM::rpolono(k, meanlog = log(mu), sdlog = sigma), ncol = k, nrow = 1)
    # M <- matrix(sads::rpoilog( k, log(mu), sig = sigma), ncol = k, nrow = 1)
    #set.seed(1, kind = "L'Ecuyer-CMRG")
    M <- matrix(rpoislog( k, log(mu), sigma, keep0 = TRUE), ncol = k, nrow = 1)
    # set.seed(1, kind = "L'Ecuyer-CMRG")
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim) {
      for (j in 1:k) {
        sim[i,j] <- stats::rbinom(1, M[,j], w[,j])
      }
    }
  } else {
    print("please choose the one of the given distribution type with case sensitive such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'")
  }
  if (summary == 1) {
    # result <- sim
    result <- round(sum(apply(sim, 2, mean)))
    # result <- mean(apply(sim, 2, mean))
    #   St.D_N <- sqrt(sum(apply(sim, 2, var)))
    #   result <- data.frame(mean_N,St.D_N)
  } else if (summary == 2) {
    # result<- rowSums(sim)
    # result <-  apply(sim, 1, sum)
    result <-  round(apply(sim, 2, mean))
  } else if (summary == 3) {
    # result<- rowSums(sim)
    # result <-  apply(sim, 1, sum)
    # result <-  sim
    result <-  round(apply(sim, 1, sum))
  } else if (summary == 4) {
    mean_N <- mean(apply(sim, 2, mean))
    Var_N <- mean(apply(sim, 2, var))
    result <- data.frame(mean_N,Var_N)
  } else if (summary == 5) {
    result <- sim
  } else {
    print("please include the summary value, which depends on your expected output")
  }
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}
