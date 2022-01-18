##' This function provides an estimated detectability value in the mixing process of the single mixing scheme.
##' @title The estimated detectability in the single mixing scheme.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit of the expected total CFUs, which can be found from stabilising point when the mean is about standard deviation
##' @param n_sim number of simulations
##' @return the estimated detectability at each stage of the mixing process.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and \eqn{l} be the number of stages in the mixing process.
##' This function provides an estimated detectability in the mixing process of the single mixing scheme. To more details, please refer the details section of  \link{compare_mixing_detectability}.
##' @seealso \link{sim_single}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' \item McCallum DA (2005) A conceptual guide to detection probability for point counts and other count-based survey methods. USDA Forest Service General Technical Report PSW-GTR-191, \href{https://www.fs.usda.gov/treesearch/pubs/32063}{pp 754-761}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- 25000
##' rate <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 138
##' n_sim <- 20000
##' sim_single_detectability(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim)
##' @export
sim_single_detectability <- function(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim){
  # f_spri <- function(mu, k, alpha, distribution) {
  #   sprintf("mixing plan (mu = %.1f, k = %.0f, l = %.1f, %s)", mu, k, l, distribution)
  # }
  pd <- sim_single_pd_stages(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim)
  detectability <- 1 - (1 - mean(pd))^k
  results <- detectability
  # colnames(results) <- f_spri(mu, k, l, distribution)
  return(results)
}






