##' This function calculates the resulting total number of colony forming units in the mixed sample in the single mixing plan with single stage of the mixing.
##' @title The total number of colony-forming units in the mixed sample by the simulation results in the single mixing plan with a single stage of the mixing.
##' @param mu the average number of colony-forming units in the mixed sample, which is in logarithmic scale if we use a Lognormal/Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param summary if we need to get all simulated \eqn{N'}, use \code{summary = FALSE} otherwise function provides mean value of the simulated \eqn{N'} ( default \code{summary = TRUE}).
##' @param n_sim number of simulations
##' @return total number of colony forming units in the single mixing plan
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}.To more details, please refer the details section of  \link{compare_mixing_stages}. (to be finished later on)
##' @seealso  \link{compare_mixing_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha <- 0.1
##' k <- 30
##' n_sim <- 2000
##' sim_single(mu, sigma, alpha, k, distribution = "Poisson lognormal-Type B", n_sim)
##' @export
sim_single <- function(mu, sigma , alpha , k, distribution, n_sim, summary = TRUE){
  if (distribution == "Poisson-Type A") {
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for (i in 1:n_sim){
    for(j in 1:k){
      sim[,j] <- stats::rpois(n_sim, mu/k)
    }
    # }
  } else if (distribution == "Poisson-Type B") {
    # x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
    # sm <- x%*%rep(1, k)
    # w <- x/as.vector(sm)
    x <-  matrix(NA, nrow = 1, ncol = k) # to apply a beta algorithm to generate Dirichlet distribution's random numbers.
    for (j in 1:k){
      x[,j] <- stats::rbeta(1,alpha, alpha*(k-j))
    }
    w <-  matrix(NA, nrow = 1, ncol = k)
    for (j in 2:k){
      w[,1] <- x[,1]
      w[,j] <- x[j] %*% prod(1 - x[1:(j-1)])
    }
    # sum(w)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for (i in 1:n_sim){
    for(j in 1:k){
      sim[,j] <- stats::rpois(n_sim, mu*w[,j])
    }
    # }
  } else if (distribution == "Lognormal-Type A") {
    M <- matrix(NA, ncol = k, nrow = 1)
    for(j in 1:k){
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      M[,j] <- as.integer(stats::rlnorm(1, meanlog = log(mu), sdlog = sigma)) # lognormal distribution
    }
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for(i in 1:n_sim){
    for (j in 1:k){
      sim[,j] <- stats::rbinom(n_sim, M[,j], 1/k)
    }
    # }
  } else if (distribution == "Lognormal-Type B") {
    # If we want to apply a gamma algorithm to generate Dirichlet distribution's random numbers.
    # x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
    # sm <- x%*%rep(1, k)
    # w <- x/as.vector(sm)
    x <-  matrix(NA, nrow = 1, ncol = k) # If we want to apply a beta algorithm to generate Dirichlet distribution's random numbers.
    for (j in 1:k){
      x[,j] <- stats::rbeta(1,alpha, alpha*(k-j))
    }
    w <-  matrix(NA, nrow = 1, ncol = k)
    for (j in 2:k){
      w[,1] <- x[,1]
      w[,j] <- x[j] %*% prod(1 - x[1:(j-1)])
    }
    # sum(w)
    M <- matrix(NA, ncol = k, nrow = 1)
    for(j in 1:k){
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      M[,j] <- as.integer(stats::rlnorm(1, meanlog = log(mu), sdlog = sigma)) # lognormal distribution
    }
    # M <- matrix(as.integer(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for (i in 1:n_sim){
    for(j in 1:k){
      sim[,j] <- stats::rbinom(n_sim, M[,j], w[,j])
    }
    # }

  } else if (distribution == "Poisson lognormal-Type A") {
    M <- matrix(NA, ncol = k, nrow = 1)
    for(j in 1:k){
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      M[,j] <- as.integer(VGAM::rpolono(1, meanlog = log(mu), sdlog =sigma)) # poisson lognormal distribution
    }
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for(i in 1:n_sim){
    for (j in 1:k){
      sim[,j] <- stats::rbinom(n_sim, M[,j], 1/k)
    }
    # }
  } else if (distribution == "Poisson lognormal-Type B") {
    # If we want to apply a gamma algorithm to generate Dirichlet distribution's random numbers.
    # x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
    # sm <- x%*%rep(1, k)
    # w <- x/as.vector(sm)
    x <-  matrix(NA, nrow = 1, ncol = k) # If we want to apply a beta algorithm to generate Dirichlet distribution's random numbers.
    for (j in 1:k){
      x[,j] <- stats::rbeta(1,alpha, alpha*(k-j))
    }
    w <-  matrix(NA, nrow = 1, ncol = k)
    for (j in 2:k){
      w[,1] <- x[,1]
      w[,j] <- x[j] %*% prod(1 - x[1:(j-1)])
    }
    # sum(w)
    M <- matrix(NA, ncol = k, nrow = 1)
    for(j in 1:k){
      # M[,j] <- as.integer(stats::rnorm(1,  mu, sigma)) # normal distribution
      M[,j] <- as.integer(VGAM::rpolono(1, meanlog = log(mu), sdlog =sigma)) # poisson lognormal distribution
    }
    # M <- matrix(as.integer(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma)), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    # for (i in 1:n_sim){
    for(j in 1:k){
      sim[,j] <- stats::rbinom(n_sim, M[,j], w[,j])
    }
    # }



  } else {
    print("please choose the one of the given distribution type with case sensitive such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'")
  }
  if (summary == TRUE){
    result <- sum(apply(sim, 2, mean))
    #   St.D_N <- sqrt(sum(apply(sim, 2, var)))
    #   result <- data.frame(mean_N,St.D_N)
  } else {
    # result<- rowSums(sim)
    result <-  apply(sim, 1, sum)
  }
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}




