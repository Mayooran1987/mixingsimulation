##' This function gives a probability of detection at each stage (or revolution) of the mixing process.
##' @title The estimated value of detection probability at each stage (or revolution) of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions /stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return The probability of detection at each stage of the mixing process.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and let \eqn{l} be the number of stages in the mixing process.
##'This function provides probability of detection at each stage of the mixing process. At each stage (or revolution), the probability of detection (\eqn{p_d}) can be estimated by using function \link{sim_single_pd}.
##' @seealso \link{sim_single_stages}
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
##' stages <- c(1:l)
##' Prob_df <-
##' data.frame(stages,sim_single_pd_stages(mu,sigma,alpha_in,k,l,r,distribution,UDL,n_sim))
##' colnames(Prob_df) <- c("no.revolutions","prob.detection")
##' plot_example <- ggplot2::ggplot(Prob_df) +
##' ggplot2::geom_line(ggplot2::aes(x = stages, y = prob.detection)) +
##' #ggplot2::stat_smooth(geom = "smooth", method = "gam", mapping = ggplot2::aes(x = no.revolutions,
##' #y = prob.detection), se = FALSE, n = 1000) +
##' ggplot2::ylab(expression("Prob. detection" ~~ (P[d[l]]))) +
##' ggplot2::theme_classic() + ggplot2::xlab(expression("No. of revolutions")) +
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) +
##' #ggplot2::ggtitle(label = f_spr(n_sim))+
##' ggthemes::scale_colour_colorblind()
##' print(plot_example)
##' @export
sim_single_pd_stages <- function(mu, sigma , alpha_in, k, l, r, distribution, UDL, n_sim){
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
  sim.sum1 <- matrix(NA, nrow = l, ncol = 1)
  for (i in 1:l) {
    sim.sum1[i,] <- sim_single_pd(mu, sigma , alpha[i,] , k, distribution, UDL, n_sim)
  }
  results <- sim.sum1
  # results <- data.frame(sim.sum1)
  # colnames(results) <- f_spri(mu, k, alpha, distribution)
  message("\033[1;31m","This simulation takes a few hours to produce the output! Thanks for your patience.")
  return(results)
}



