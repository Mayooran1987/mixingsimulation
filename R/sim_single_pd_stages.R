##' This function gives a probability of detection at each stage of the mixing process.
##' @title The estimated value of detection probability at each stage of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions /stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit of the expected total CFUs, which can be found from stabilising point when the mean is about standard deviation
##' @param n_sim number of simulations
##' @return The probability of detection at each stage of the mixing process.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and let \eqn{l} be the number of stages in the mixing process.
##'This function provides probability of detection at each stage of the mixing process. It can be determined by how many of the simulated samples exceed the upper limit set for the CFUs at each stage of the mixing process.
##'
##' The probability of detection (\eqn{p_d}) can be estimated from following formula,
##'
##' \deqn{p_d = \frac{\textnormal{number of simulated samples which are greater than UDL}}{\textnormal{number of simulations}} ;}
##'
##' where the upper decision limit (UDL) can be found from a stabilising point in the graphical display of the expected total number of CFUs versus the number of revolutions when the mean is about standard deviation.
##'
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
##' rate <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 138
##' n_sim <- 20000
##' no.revolutions <-c(1:l)
##' Prob_df <- data.frame(no.revolutions,
##' sim_single_pd_stages(mu,sigma,alpha_in,k,l,rate,distribution,UDL,n_sim))
##' colnames(Prob_df) <- c("no.revolutions","P_d")
##' cummean <- function(x){cumsum(x)/seq_along(x)}
##' cum_mean <- cummean(Prob_df[,2])
##' plot_example <- ggplot2::ggplot(Prob_df) +
##'   ggplot2::geom_line(ggplot2::aes(x = no.revolutions, y = P_d))+
##'   # ggplot2::geom_line(ggplot2::aes(x = log10(mean), y = P_d))+
##'   ggplot2::geom_line( ggplot2::aes(x = no.revolutions,y = cum_mean),color = "red",size = .75)+
##'   ggplot2::xlab(expression("Number of revolutions"))+
##'   ggplot2::ylab(expression("Probability of detection"))+
##'   ggplot2::theme_classic()+
##'   ggplot2::ggtitle(label = "Probability of detection versus number of revolutions")+
##'   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
##'   ggthemes::scale_colour_colorblind()
##' print(plot_example)
##' @export
sim_single_pd_stages <- function(mu, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, l = %.1f, %s)", mu, k, l, distribution)
  }
  X <- sim_single_stages(mu, sigma , alpha_in, k, l, rate, distribution, n_sim, summary = 3)
  sim.sum1 <- matrix(NA, nrow = l, ncol = 1)
  for (j in 1:l) {
    # sim.sum1[j,] <- length(which(sim_single(mu, sigma, alpha[1,j], k, distribution, n_sim, summary = FALSE) > UDL))/n_sim
    sim.sum1[j,] <- length(which(X[j,] > UDL))/n_sim
  }
  results <- sim.sum1
  # results <- data.frame(sim.sum1)
  colnames(results) <- f_spri(mu, k, l, distribution)
  return(results)
}
