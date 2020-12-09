##' This function calculates the resulting total number of colony forming units in the mixed sample in the single mixing plan. (to be finished later on)
##' @title The total number of colony-forming units in the mixed sample by simulation result (in the single mixing plan).
##' @param mu the average number of colony-forming units in the mixed sample, which is in logarithmic scale if we use a lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the log10 scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"}
##' @param summary if we need to get the mean and standard deviation of simulated \eqn{N'}, use summary = TRUE ( default summary =FALSE).
##' @param n_sim number of simulations
##' @return total number of colony forming units in the single mixing plan
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}.To more details, please refer the details section of  \link{compare_mixing}. (to be finished later on)
##' @seealso \link{sim_multiple}, \link{compare_mixing}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha <- 0.1
##' k <- 10
##' n_sim <- 20000
##' sim_single(mu, sigma, alpha, k, distribution = "Poisson-Type A", n_sim, summary = TRUE)
##' @export
sim_single <- function(mu, sigma , alpha , k, distribution, n_sim, summary = FALSE){
  if (distribution == "Poisson-Type A") {
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (j in 1:n_sim){
      sim[j,] <- stats::rpois(k, mu/k)
    }
  } else if (distribution == "Poisson-Type B") {
    x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
    sm <- x%*%rep(1, k)
    w <- x/as.vector(sm)
    # y<-matrix(stats::rgamma(k,alpha,1),ncol=k, nrow=1)
    # w <- y[,i]/sum(y[,i])
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim){
      for(j in 1:k){
        sim[i,j] <- stats::rpois(1, mu*w[,j])
      }
    }
  } else if (distribution == "Lognormal-Type A") {
    M <- matrix(stats::rlnorm(k, meanlog = log(mu), sdlog = sigma), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim){
      for(j in 1:k){
        sim[i,j] <- stats::rbinom(1, floor(M[,j]), 1/k)
      }
    }
  } else if (distribution == "Lognormal-Type B") {
    x <- matrix(stats::rgamma(k,alpha), ncol = k, nrow = 1)
    sm <- x%*%rep(1, k)
    w <- x/as.vector(sm)
    M <- matrix(rlnorm(k, meanlog = log(mu), sdlog = sigma), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_sim, ncol = k)
    for (i in 1:n_sim){
      for(j in 1:k){
        sim[i,j] <- stats::rbinom(1, floor(M[,j]), w[,j])
      }
    }
  } else {
    print("please choose the one of the given distribution type with case sensitive such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'")
  }
  if (summary == TRUE){
    mean_N <- sum(apply(sim, 2, mean))
    St.D_N <- sqrt(sum(apply(sim, 2, var)))
    result <- data.frame(mean_N,St.D_N)
  } else {
    result<-  apply(sim, 1, sum)
  }
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}





