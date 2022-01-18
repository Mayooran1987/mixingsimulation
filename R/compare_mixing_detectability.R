##' This function provides a graphical display to compare mixing plans based on the estimated detectability in the mixing process using different mixing parameters such as revolutions, type of distribution and number of primary samples.
##' @title Graphical comparison of mixing plans based on estimated detectability in the mixing process.
##' @param mulower the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muupper the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sigma the standard deviation of the colony-forming units (CFUs) in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in each of the revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UDL the upper decision limit of the expected total CFUs, which can be found from a stabilising point when the mean is about standard deviation
##' @param n_sim number of simulations
##' @return Estimates the detectability-based graphical display for a comparison of mixing plans.
##' @details Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by a contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i},
##'  and let \eqn{l} be the number of revolutions in the mixing. This function provides a graphical display to compare mixing plans based on the estimated average
##'  detectability in the mixing process with different input variables, such as number of revolutions, type of distribution and number of primary samples.
##'
##' The detectability is given by the following formula,
##'
##' \deqn{detectability = 1- {(1-p_d)}^k ;}
##'
##'where \eqn{p_d} is the average probability of detection after \eqn{l} number of revolutions in each mixing plan. The probability of detection at every stage of the mixing process can be estimated by employing function \link{sim_single_pd_stages}. However, if we want to estimate detectability values in each revolution of the mixing process of each plan, we have to utilise function \link{sim_single_detectability}. We can flexibly change the mixing parameters of this function, which depends on what purpose the comparison is needed for.
##' @seealso \link{sim_single_detectability}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' \item McCallum DA (2005) A conceptual guide to detection probability for point counts and other count-based survey methods. USDA Forest Service General Technical Report PSW-GTR-191, \href{https://www.fs.usda.gov/treesearch/pubs/32063}{pp 754-761}.
##' }
##' @examples
##' \dontrun{
##' mulower <- 0
##' muupper <- 200
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(30,30)
##' l <- c(500,25000)
##' rate <- 0.01
##' distribution <- c("Poisson lognormal-Type B","Poisson lognormal-Type B")
##' UDL <- 138
##' n_sim <- 20000
##' compare_mixing_detectability(mulower,muupper,sigma,alpha_in,k,l,rate,distribution,UDL,n_sim)}
##' @export
compare_mixing_detectability <-  function(mulower, muupper, sigma , alpha_in, k, l, rate, distribution, UDL, n_sim){
  mixing_scheme <- NULL
  detectability <- NULL
  f_spri <- function(l, k, distribution) {
    sprintf("mixing plan (k = %.0f, l = %.0f, %s)", k, l, distribution)
    }
  # f_spr <- function(n_sim ) {
  #   sprintf("Simulation results (no.simulations = %.0f)", n_sim)
  # }
  mu <- seq(mulower, muupper, 0.01)
  # stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = length(mu), ncol = length(distribution))
  for (i in 1:nrow(sim.sum3)) {
    for (j in 1:ncol(sim.sum3)) {
      sim.sum3[i,j] <-  sim_single_detectability(mu[i], sigma , alpha_in, k[j], l[j], rate, distribution[j], UDL, n_sim)
    }
    }
  result <- data.frame(mu, sim.sum3)
  colnames(result) <- c("mu", f_spri(l, k, distribution))
  melten.Prob <- reshape2::melt(result, id = "mu", variable.name = "mixing_scheme", value.name = "detectability")
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(detectability, group = mixing_scheme, colour = mixing_scheme)) +
    # ggplot2::geom_line(ggplot2::aes(x = mu, y = detectability))+
    ggplot2::geom_smooth(stat = "smooth",  method = 'gam', formula = y ~ s(x, bs = "cs"), mapping = ggplot2::aes(x = mu, y = detectability), se = FALSE) +
    ggplot2::ylab(expression("Detectability")) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("Mean concentration (" ~ mu*~")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) +
    # ggplot2::ggtitle(label = f_spr(n_sim))+
    ggthemes::scale_colour_colorblind()
    # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(plot1)
}

