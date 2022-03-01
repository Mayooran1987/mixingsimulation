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
##' \deqn{p_d = \frac{\textnormal{Number of simulated primary samples which contain CFUs greater than UDL}}{\textnormal{Number of primary samples}} ;}
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
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
  }
  p_d <- length(which(sim_single(mu, sigma , alpha , k, distribution, n_sim, summary = 2) > UDL))/k
  # sim_single(mu, sigma , alpha , k, distribution, n_sim, summary = 5)
  # p_d <-length(which( X > UDL))/n_sim
  results <- p_d
  # results <- data.frame(sim.sum1)
  # colnames(results) <- f_spri(mu, k, alpha, distribution)
  return(results)
}




