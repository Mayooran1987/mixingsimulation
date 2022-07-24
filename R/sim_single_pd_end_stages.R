##' This function gives a probability of detection at end of stage (or revolution) of the mixing process.
##' @title The estimated value of detection probability at end of stage (or revolution) of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions /stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return The probability of detection at end of stage of the mixing process.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and let \eqn{l} be the number of stages in the mixing process.
##'This function provides probability of detection at each stage of the mixing process. At each stage (or revolution), the probability of detection (\eqn{p_d}) can be estimated by using function \link{sim_single_pd}.
##' @seealso \link{sim_single_pd_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- 25000
##' r <- 0.01
##' distribution <- "Poisson lognormal-Type B"
##' UDL <- 0
##' n_sim <- 2000
##' sim_single_pd_end_stages(mu, sigma , alpha_in, k, l, r, distribution, UDL, n_sim)
##' @export
sim_single_pd_end_stages <- function(mu, sigma , alpha_in, k, l, r, distribution, UDL, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
  }
  alpha <- matrix(NA, nrow = l , ncol = 1)
  for (j in 1:l) {
    if (j == 1) {
      alpha[j,] <- alpha_in
    } else {
      alpha[j,] <- alpha[j - 1,] + r
    }
  }
  results <- sim_single_pd(mu, sigma , alpha[l,] , k, distribution, UDL, n_sim)
  return(results)
}
