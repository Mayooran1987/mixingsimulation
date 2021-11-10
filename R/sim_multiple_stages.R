##' This function calculates the resulting expected total number of CFUs in the mixed sample in the multiple mixing plans at each stage of the mixing process.
##' @title The expected total number of CFUs in the mixed sample in the multiple mixing schemes at each stage of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return The expected total number of CFUs in each revolution/stage.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the mixing of \eqn{k}  primary samples and \eqn{N' = \sum N_i} and let \eqn{N_i} be the number of CFUs.
##'For this package development, we have employed the notations 'Type-A' and 'Type-B' to indicate the type of distributions, which are applied in the previous literature as 'fair' and 'beta', respectively; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##'
##' This package will consider stage-by-stage the mixing process and assumes systematically breaking clusters at every stage of the mixing. Therefore, it can be assumed the concentration parameter also systematically changes with the concentration of the contribution.
##' @seealso \link{sim_single}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' rate <- 0.01
##' l <- 25000
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(30,75)
##' distribution <- c("Poisson lognormal-Type B","Poisson lognormal-Type B")
##' n_sim <- 20000
##' colMeans(sim_multiple_stages(mu, sigma, alpha_in, k, l, rate, distribution, n_sim))
##' @export
sim_multiple_stages <- function(mu, sigma, alpha_in, k, l, rate, distribution, n_sim){
  f_spri <- function(mu, k, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, %s)", mu, k, distribution)
    }
  # stages <- 1:l
# sim.sum <- matrix(NA, nrow = n_sim, ncol = length(k))
# for(i in 1:length(k)){
  sim.sum3 <- matrix(NA, nrow = l, ncol = length(k))
  for(j in 1:length(k)){
    sim.sum3[,j] <-  sim_single_stages(mu, sigma, alpha_in, k[j], l, rate, distribution[j], n_sim)
  }
  result <- data.frame(sim.sum3)
  colnames(result) <- f_spri(mu, k, distribution)
  # result <- data.frame(stages, mean(sim.sum3))
  # colnames(result) <- c("stages", f_spri(mu, k, distribution))
  return(result)
}



