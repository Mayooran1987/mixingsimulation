##' This function provides a graphical display to compare mixing plans based on the estimated cumulative moving average (CMA) of detection probability at each revolution of the mixing process.
##' @title Graphical comparison of mixing plans based on the estimated cumulative moving average (CMA) of detection probability at each revolution.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units (CFUs) in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit, which depends on the type of microorganisms and testing regulations.
##' @param n_sim number of simulations
##' @return graphical display of estimated cumulative moving average (CMA) of detection probability at each revolution in the mixing.
##' @seealso \link{sim_single_pd_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(10,30,60)
##' l <- 2500
##' r <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' UDL <- 0
##' n_sim <- 2000
##' compare_mixing_1(mu,sigma , alpha_in, k, l, r, distribution, UDL, n_sim)
##' @export
compare_mixing_1 <-  function(mu,sigma , alpha_in, k, l, r, distribution, UDL, n_sim){
  mixing_scheme <- NULL
  prob.detection <- NULL
  f_spri <- function(k, distribution) {
    sprintf("mixing plan (k = %.0f, %s)", k, distribution)
  }
  cummean <- function(x){cumsum(x)/seq_along(x)} # to get a cumulative average calculation
  # f_spr <- function(n_sim ) {
  #   sprintf("Simulation results (no.simulations = %.0f)", n_sim)
  # }
  # mu <- seq(mulower, muupper, 0.1)
  # stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = l, ncol = length(k))
  set.seed(1, kind = "L'Ecuyer-CMRG")
  # for (i in 1:nrow(sim.sum3)) {
  for (j in 1:ncol(sim.sum3)) {
    #set.seed(1000+j)
    # sim.sum3[,j] <-  sim_single_pd_stages(mu, sigma , alpha_in, k[j], l, r, distribution, UDL, n_sim)
    sim.sum3[,j] <-  cummean(sim_single_pd_stages(mu, sigma , alpha_in, k[j], l, r, distribution, UDL, n_sim))
  }
  # }
  stages <- 1:l
  # p_d <- sim_single_pd_stages(mu, sigma , alpha_in, k, l, r, distribution, UDL, n_sim)
  # p_a <- (1-p_d)^k
  # result <- data.frame(stages,p_d,p_a)
  result <- data.frame(stages, sim.sum3)
  colnames(result) <- c("stages", f_spri(k, distribution))
  melten.Prob <- reshape2::melt(result, id = "stages", variable.name = "mixing_scheme", value.name = "prob.detection")
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(prob.detection, group = mixing_scheme, colour = mixing_scheme)) +
    ggplot2::geom_line(ggplot2::aes(x = stages, y = prob.detection)) +
    # ggplot2::stat_smooth(geom =  "smooth",  method = "gam", mapping = ggplot2::aes(x = stages, y = prob.detection), se = TRUE, n = 1000) +
    ggplot2::ylim(NA,1) +
    ggplot2::ylab(expression("Cumulative moving average of prob. detection"~~~(~bar(P[d[l]])))) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("No. of revolutions")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) +
    # ggplot2::ggtitle(label = f_spr(n_sim))+
    ggthemes::scale_colour_colorblind()
  # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(plot1)
}

