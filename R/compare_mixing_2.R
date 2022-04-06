##' This function provides a graphical display to compare mixing plans using a different number of revolutions based on the estimated average probability of detection after mixing is ended.
##' @param mulower the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muupper the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sigma the standard deviation of the colony-forming units (CFUs) in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in each of the revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return graphical display compares mixing plans using a different number of revolutions based on the estimated average probability of detection after mixing is ended.
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
##' l <- c(50,5000)
##' rate <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 0
##' n_sim <- 2000
##' compare_mixing_2(mulower, muupper, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim)
##' @export
compare_mixing_2 <-  function(mulower, muupper, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim){
  mixing_scheme <- NULL
  prob.detection <- NULL
  f_spri <- function(l, k, distribution) {
    sprintf("mixing plan (k = %.0f, l = %.0f, %s)", k, l, distribution)
  }
  # f_spr <- function(n_sim ) {
  #   sprintf("Simulation results (no.simulations = %.0f)", n_sim)
  # }
  mu <- seq(mulower, muupper, 0.1)
  # stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = length(mu), ncol = length(l))
  # sim.sum3 <- matrix(NA, nrow = length(mu), ncol = 1)
  set.seed(1, kind = "L'Ecuyer-CMRG")
  for (i in 1:nrow(sim.sum3)) {
    for (j in 1:ncol(sim.sum3)) {
      # sim.sum3[i,j] <-  sim_single_pd_stages(mu[i], sigma , alpha_in, k, l[j], rate, distribution, UDL, n_sim)[l[j]] # If we want to use the probability of detection at the end of mixing, please use this.
      sim.sum3[i,j] <-  mean(sim_single_pd_stages(mu[i], sigma , alpha_in, k, l[j], rate, distribution, UDL, n_sim))
    }
  }
  result <- data.frame(mu, sim.sum3)
  colnames(result) <- c("mu", f_spri(l, k, distribution))
  melten.Prob <- reshape2::melt(result, id = "mu", variable.name = "mixing_scheme", value.name = "prob.detection")
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(prob.detection, group = mixing_scheme, colour = mixing_scheme)) +
    # ggplot2::geom_line(ggplot2::aes(x = mu, y = prob.detection)) +
    ggplot2::geom_smooth(stat = "smooth",  method = 'gam', formula = y ~ s(x, bs = "cs"), mapping = ggplot2::aes(x = mu, y = prob.detection), se = FALSE) +
    ggplot2::ylim(0,1) +
    ggplot2::ylab(expression(" Average prob.detection after"~~l~~ "revolutions"~~ (bar(P[d[l]])))) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("Mean concentration (" ~ mu*~")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) +
    # ggplot2::ggtitle(label = f_spr(n_sim))+
    ggthemes::scale_colour_colorblind()
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(plot1)
}
