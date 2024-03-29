##' This function provides a graphical display to compare mixing plans using a different number of revolutions based on the estimated probability of detection at the end of mixing.
##' @title Graphical comparison of mixing plans based on estimated probability of detection at the end of the mixing with different revolutions
##' @param mulower the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muupper the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sigma the standard deviation of the colony-forming units (CFUs) in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return graphical display compares mixing plans using a different number of revolutions based on the estimated probability of detection at the end of the mixing.
##' @seealso \link{sim_single_pd_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mulower <- 0.1
##' muupper <- 200
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- c(500,5000)
##' r <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 0
##' n_sim <- 2000
##' compare_mixing_2(mulower, muupper, sigma , alpha_in, k, l, r, distribution, UDL, n_sim)
##' @export
compare_mixing_2 <-  function(mulower, muupper, sigma , alpha_in, k, l, r, distribution, UDL, n_sim){
  mixing_scheme <- NULL
  prob.detection <- NULL
  f_spri <- function(l, k, distribution) {
    sprintf("mixing plan (k = %.0f, l = %.0f, %s)", k, l, distribution)
  }
  mu <- seq(mulower, muupper, 0.1)
  # sim.sum3 <- matrix(NA, nrow = length(mu), ncol = length(l))
  # # set.seed(1, kind = "L'Ecuyer-CMRG")
  # for (i in 1:nrow(sim.sum3)) {
  #   for (j in 1:ncol(sim.sum3)) {
  #     sim.sum3[i,j] <-  sim_single_pd_stages(mu[i], sigma , alpha_in, k, l[j], r, distribution, UDL, n_sim)[l[j]] # If we want to use the probability of detection at the end of mixing, please use this.
  #   }
  # }

  sim.sum3 <- matrix(NA, nrow = length(mu), ncol = 2)
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  for (i in 1:nrow(sim.sum3)) {
    sim.sum3[i,1] <-  sim_single_pd_end_stages(mu[i], sigma , alpha_in, k, l[1], r, distribution, UDL, n_sim)
    sim.sum3[i,2] <-  sim_single_pd_end_stages(mu[i], sigma , alpha_in, k, l[2], r, distribution, UDL, n_sim)
  }
  # MA <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
  result <- data.frame(mu, sim.sum3)
  # result <- data.frame(mu, MA(sim.sum3,50))
  colnames(result) <- c("mu", f_spri(l, k, distribution))
  melten.Prob <- reshape2::melt(result, id = "mu", variable.name = "mixing_scheme", value.name = "prob.detection")
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(prob.detection, group = mixing_scheme, colour = mixing_scheme)) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = prob.detection)) +
    # ggplot2::geom_smooth(stat = "smooth",  method = 'gam', formula = y ~ s(x, bs = "cs"), mapping = ggplot2::aes(x = mu, y = prob.detection), se = FALSE) +
    # ggplot2::ylim(0,1) +
    ggplot2::ylab(expression("Prob.detection at the end of the mixing"~ (P[d[l]]))) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("Mean concentration (" ~ mu*~")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) +
    # ggplot2::ggtitle(label = f_spr(n_sim))+
    ggthemes::scale_colour_colorblind()
  # plot1
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(plot1)
}
