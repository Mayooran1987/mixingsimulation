##' This function provides a graphical display to compare mixing plans based on estimated average prevalence in the mixing process using different mixing parameters such as revolutions, type of distribution and number of primary samples.
##' @title Graphical comparison of mixing plans based on estimated average prevalence in the mixing process.
##' @param mulower the lower value of the mean concentration (\eqn{\mu}), which is desired to use in the graphical display's x-axis.
##' @param muupper the upper value of the mean concentration (\eqn{\mu}), which is desired to use in the graphical display's x-axis.
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param UL the upper limit value of the expected total CFU, which can be found from stabilising point when the mean is about standard deviation
##' @param n_sim number of simulations
##' @return estimated average prevalence based graphical display to the comparison of mixing plans.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and \eqn{l} be the number of revolutions in the mixing process.
##' This function provides a graphical display to compares mixing plans based on estimated average prevalence in the mixing process with different input variables such as number of revolutions, type of distribution and number of primary samples.
##'
##' The prevalence is given by the following formula,
##'
##' \deqn{Prevalence = 1- {(1-p_d)}^l ;}
##'
##' where \eqn{p_d} is the probability of detection in each stage of the mixing process, for the comparison purpose, we have applied average prevalence after \eqn{l} number of revolutions in each mixing plan.
##' The probability of detection at every stage of the mixing process can be estimated by employing function \link{sim_single_pd_stages}.
##'
##' However, if we want to estimate prevalence values in each revolution of the mixing process of each plan, we have to utilise function \link{sim_single_prevalence_stages}.
##' We can flexibly change the mixing parameters of this function which is depending on what purpose of the comparison is needed.
##' @seealso \link{sim_single_prevalence_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mulower <- 50
##' muupper <- 200
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(30,30)
##' l <- c(500,25000)
##' rate <- 0.01
##' distribution <- c("Poisson lognormal-Type B","Poisson lognormal-Type B")
##' UL <- 138
##' n_sim <- 2000
##' compare_mixing_prevalence(mulower, muupper, sigma, alpha_in, k, l, rate, distribution, UL, n_sim)
##' @export
compare_mixing_prevalence <-  function(mulower, muupper, sigma , alpha_in, k, l, rate, distribution, UL, n_sim){
  mixing_scheme <- NULL
  Prevalence <- NULL
  f_spri <- function(l, k, distribution) {
    sprintf("mixing plan (k = %.0f,l = %.0f, %s)", k, l, distribution)
  }
  f_spr <- function(n_sim ) {
    sprintf("Simulation results (no.simulations = %.0f)", n_sim)
  }
    mu <- seq(mulower, muupper, 0.1)
    # stages <- 1:l
    sim.sum3 <- matrix(NA, nrow = length(mu), ncol = length(distribution))
    # for(i in 1:length(mu)){
    #   sim.sum3[i,] <-  sim_multiple_prevalence_stages(mu[i], sigma , alpha_in, k, l, rate, distribution, UL, n_sim)
    # }
    for(i in 1:nrow(sim.sum3)){
      for(j in 1:ncol(sim.sum3)){
        sim.sum3[i,j] <-  mean(sim_single_prevalence_stages(mu[i], sigma , alpha_in, k[j], l[j], rate, distribution[j], UL, n_sim))
      }
    }
    result <- data.frame(mu, sim.sum3)
    colnames(result) <- c("mu", f_spri(l, k, distribution))
    melten.Prob <- reshape2::melt(result, id = "mu", variable.name = "mixing_scheme", value.name = "Prevalence")
    plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(Prevalence, group = mixing_scheme, colour = mixing_scheme)) +
      ggplot2::geom_line(ggplot2::aes(x = mu, y = Prevalence))+
      # theme_classic()
      ggplot2::ylab(expression("Prevalence"))+
      ggplot2::theme_classic()+ ggplot2::xlab(expression("Mean concentration"))+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.20,0.75))+
        ggplot2::ggtitle(label = f_spr(n_sim))+ ggthemes::scale_colour_colorblind()+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    # cat("Calculation took", proc.time()[1], "seconds.\n")
    return(plot1)
  }

