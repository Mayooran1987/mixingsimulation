##' This function calculates the resulting expected total number of colony-forming units in the mixed sample in the multiple mixing plans at the single stage of the mixing process.
##' @title The expected total number of colony-forming units in the mixed sample in the multiple mixing schemes at the single stage of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return total number of colony forming units in the multiple mixing scheme
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by contribution of \eqn{k} primary samples mixing and \eqn{N' = \sum N_i}. This function provides the simulated resulting of the expected total number of colony-forming units in the mixed sample in the multiple mixing plans at the single stage of the mixing process.
##' To more details, please refer the details section of  \link{compare_mixing_3}.
##' @seealso \link{sim_single}, \link{compare_mixing_3}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' set.seed(1350)
##' sigma <- 0.8
##' alpha <- c(0.1,5)
##' k <- c(30,30)
##' distribution <-  c("Poisson lognormal-Type B","Poisson lognormal-Type B")
##' n_sim <- 20000
##' f_spr <- function(n_sim) {
##'   sprintf("Simulation results (no.simulations = %.0f)", n_sim)
##' }
##' f_spri <- function(alpha, distribution) {
##'   sprintf("mixing plan (alpha = %.1f, %s)", alpha, distribution)
##' }
##' mu <- seq(100, 200, 0.1)
##' sim.sum3 <- matrix(NA, nrow = length(mu), ncol = length(distribution))
##' for(i in 1:nrow(sim.sum3)){
##'   sim.sum3[i,] <-  colMeans(sim_multiple(mu[i], sigma, alpha, k, distribution, n_sim))
##' }
##' result <- data.frame(mu, sim.sum3)
##' colnames(result) <- c("mu", f_spri(alpha, distribution))
##' melten.Prob <- reshape2::melt(result, id = "mu", variable.name = "mixing_scheme",
##'                               value.name = "Total_CFU")
##' plot_example <-
##' ggplot2::ggplot(melten.Prob, ggplot2::aes(Total_CFU, group = mixing_scheme,colour = mixing_scheme))+
##'   ggplot2::geom_line(stat="density",ggplot2::aes(x = Total_CFU))+
##'   ggplot2::ylab(expression("pmf"))+
##'   ggplot2::theme_classic()+ ggplot2::xlab(expression("Total number of CFU in the mixed sample"))+
##'   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.75))+
##'   ggthemes::scale_colour_colorblind()
##'   plot_example
##' @export
sim_multiple <- function(mu, sigma, alpha, k, distribution, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
    }
  if (length(alpha) != length(k)) {
    warning("length of alpha and length of k are must be equal")
    } else {
      sim.sum1 <- matrix(NA, nrow = n_sim, ncol = length(k))
      for (j in 1:length(k)) {
        sim.sum1[,j] <-  sim_single(mu, sigma, alpha[j], k[j], distribution[j], n_sim, summary = 3)
      }
      }
  result <- data.frame(sim.sum1)
  colnames(result) <- f_spri(mu, k, alpha, distribution)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
  }

