##' This function gives a simulated number of CFUs after each stage of the mixing process.
##' @title The total number of colony-forming units in the mixed sample by the simulation results in the single mixing plan with \eqn{l}  number of stages.
##' @param l the maximum number of stages in the mixing process
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return average number of colony forming units in the single mixing plan with \eqn{l} number of stages.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and \eqn{l} be the number of stages in the mixing process.
##' This function provides simulated number of CFUs after each stages of the mixing process. To more details, please refer the details section of  \link{compare_mixing_3}.
##' @seealso \link{sim_single}
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
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 2000
##' stages <-c(1:l)
##' Prob_df <-
##' data.frame(stages,sim_single_stages(mu,sigma,alpha_in,k,l,r,distribution,n_sim))
##' colnames(Prob_df) <- c("no.revolutions","CFU")
##' plot_sim_single_stages <- ggplot2::ggplot(Prob_df,ggplot2::aes(x = no.revolutions, y = CFU)) +
##' ggplot2::geom_line() +
##' # tidyquant::geom_ma(ma_fun = SMA, n = 50, linetype = 1, colour = "red") +
##' ggplot2::xlab(expression("Number of revolutions")) +
##' ggplot2::ylab(expression("Expected total number of CFUs")) +
##' ggplot2::theme_classic() +
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
##' ggthemes::scale_colour_colorblind()
##' print(plot_sim_single_stages)
##' @export
sim_single_stages <- function(mu, sigma , alpha_in, k, l, r, distribution, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, l = %.1f, %s)", mu, k, l, distribution)
  }
  # stages <- 1:l
  alpha <- matrix(NA, nrow = 1 , ncol = l)
  for (j in 1:l) {
    if (j == 1) {
      alpha[,j] <- alpha_in
    } else {
      alpha[,j] <- alpha[,j - 1] + r
    }
  }
    sim.sum1 <- matrix(NA, nrow = l, ncol = 1)
    for (j in 1:length(alpha)) {
      sim.sum1[j,] <- sim_single(mu, sigma, alpha[,j], k, distribution, n_sim, summary = 1)
    }
    results <- sim.sum1
    colnames(results) <- f_spri(mu, k, l, distribution)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  #  message("Total number of colony-forming units in the mixed sample")
  return(results)

}

