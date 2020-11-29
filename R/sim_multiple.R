##' This function calculates the resulting total number of colony forming units in the mixed sample in the multiple mixing plans. (to be finished later on)
##' @title The total number of colony-forming units in the mixed sample by simulation result (in the multiple mixing plan).
##' @param n_iter number of iterations
##' @param mu average number of colony-forming units in mixed sample which is in logarithmic scale if we use a lognormal distribution
##' @param sigma log standard deviation of the colony-forming units in the mixed sample
##' @param b concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"}
##' @return total number of colony forming units in the multiple mixing scheme
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}. To more details, please refer the details section of  \link{compare_mixing}. (to be finished later on)
##' @seealso \link{sim_single}, \link{compare_mixing}
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
##'distribution <-  c("Poisson-Type A","Poisson-Type A","Poisson-Type A","Poisson-Type A")
##' head(sim_multiple( n_iter, mu, sigma, b, k, distribution))
##' @export
sim_multiple <- function(n_iter, mu, sigma, b, k, distribution){
  f_spri <- function(mu, k, b, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, b = %.1f, %s)", mu, k, b, distribution)
  }
  #f_spri(mu, k, b, distribution)
  if (length(mu)!=length(k)) {
    warning("length of mu and length of k are must be equal")
  } else {
  # if (i in 1:length(k)){
  sim.sum1 <- matrix(NA, nrow = n_iter, ncol = length(k))
  for(j in 1:length(k)){
    sim.sum1[,j] <-  sim_single(n_iter, mu[j], sigma[j], b[j], k[j], distribution[j])
  }
   }
  result <- data.frame(sim.sum1)
  colnames(result) <- f_spri(mu, k, b, distribution)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}

