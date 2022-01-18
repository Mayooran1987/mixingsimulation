##' This function estimates the detectability in the mixing process after a specific number of revolutions for different mixing schemes.
##' @title The estimated average detectability value in the multiple mixing plans.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit of the expected total CFUs, which can be found from stabilising point when the mean is about standard deviation
##' @param n_sim number of simulations
##' @return The detectability in the mixing process after a specific number of revolutions.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and let \eqn{l} be the number of stages in the mixing process. This function estimates the detectability value after a specific number of revolutions in each mixing scheme. However, we need to apply function \link{sim_single_detectability}  if we want to estimate individual detectability values at each stage of the mixing process.
##' @seealso \link{sim_single_detectability}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' \item McCallum DA (2005) A conceptual guide to detection probability for point counts and other count-based survey methods. USDA Forest Service General Technical Report PSW-GTR-191, \href{https://www.fs.usda.gov/treesearch/pubs/32063}{pp 754-761}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(30,75)
##' l <- 25000
##' rate <- 0.01
##' distribution <- c("Poisson lognormal-Type B","Poisson lognormal-Type B")
##' UDL <- 138
##' n_sim <- 20000
##' sim_multiple_detectability(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim)
##' @export
sim_multiple_detectability <- function(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim){
  # Total_CFU <- NULL
  # mixing_scheme <- NULL
  f_spri <- function(mu, k, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, %s)", mu, k, distribution)
  }
  # f_spr <- function(l) {
  #   sprintf("Simulation results (no.revolutions = %.0f)", l)
  # }
  # stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = 1, ncol = length(distribution))
  for (j in 1:length(distribution)) {
    sim.sum3[,j] <-  sim_single_detectability(mu, sigma , alpha_in, k[j], l, rate, distribution[j], UDL, n_sim)
  }
  result <-  data.frame(sim.sum3)
  colnames(result) <- c(f_spri(mu, k, distribution))
  return(result)
}

