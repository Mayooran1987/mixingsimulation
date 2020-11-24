##' This function calculates the resulting total number of colony forming units in the mixed sample in the single mixing plan. (to be finished later on)
##' @title The total number of colony-forming units in the mixed sample by simulation result (in the single mixing plan).
##' @param n_iter number of iterations
##' @param mu average number of colony-forming units in a primary sample which is in logarithmic scale if we use a lognormal distribution
##' @param sigma log standard deviation of the colony-forming units in a primary sample
##' @param b concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as 'Poisson-fair' or 'Poisson-beta' or 'Lognormal-fair' or 'Lognormal-beta'
##' @param summary if we need to get the mean and standard deviation of simulated \eqn{N'}, use summary = TRUE ( default summary =FALSE).
##' @return total number of colony forming units in the single mixing plan
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary sample and \eqn{N' = \sum(N_i)} (to be finished later on)
##' @seealso \link{sim_multiple}, \link{compare_mixing}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' n_iter <- 200000
##' mu <- c(20,30,50,60)
##' sigma <- c(0.8,0.8,0.8,0.8)
##' b <- c(0.1,0.1,0.1,0.1)
##' k <- c(10,10,10,10)
##'distribution <-  c("Poisson-fair","Poisson-beta","Lognormal-fair","Lognormal-beta")
##' sim_single( n_iter, mu[2], sigma[2], b[2], k[2], distribution[2], summary = TRUE)
##' @export
sim_single <- function(n_iter, mu, sigma, b, k, distribution, summary = FALSE){
  if (distribution == "Poisson-fair") {
    sim <-  matrix(NA, nrow = n_iter, ncol = k)
    for(j in 1:k){
      sim[,j] <-  stats::rpois(n_iter, mu)
    }
  } else if (distribution == "Poisson-beta") {
    for (i in 1:length(k)) {
      # alpha <- rep(b,k)
      y<-matrix(stats::rgamma(k,b,1),ncol=k, nrow=1)
      w <- y[,i]/sum(y[,i])
      sim <-  matrix(NA, nrow = n_iter, ncol = k)
      for(j in 1:k){
        sim[,j] <-  rpois(n_iter, mu*w)
      }
    }
  } else if (distribution == "Lognormal-fair") {
    M <- matrix(stats::rlnorm(k, meanlog = mu, sdlog = sigma), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_iter, ncol = k)
    for(j in 1:k){
      sim[,j] <-  rbinom(n_iter, floor(M[,j]), 1/k)
    }
  } else if (distribution == "Lognormal-beta") {
    y <- matrix(stats::rgamma(k,b), ncol = k, nrow = 1)
    sum <- apply(y, 1, sum)
    w <- matrix(NA, ncol=k, nrow=1)
    for( j in 1:k) {
      w[ ,j] <-  y[,j]/sum
    }
    M <- matrix(rlnorm(k, meanlog = mu, sdlog = sigma), ncol = k, nrow = 1)
    sim <-  matrix(NA, nrow = n_iter, ncol = k)
    for(j in 1:k){
      sim[,j] <-  stats::rbinom(n_iter, floor(M[,j]), w[,j])
    }

  } else {
    print("please choose the one of the given distribution with case sensitive such as 'Poisson-fair' or 'Poisson-beta' or 'Lognormal-fair' or 'Lognormal-beta' ")
  }
  if (summary == TRUE){
    mean_N <- sum(apply(sim, 2, mean))
    St.D_N <- sqrt(sum(apply(sim, 2, var)))
    result <- data.frame(mean_N,St.D_N)
  } else {
    result<-  apply(sim, 1, sum)
  }
  cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}



