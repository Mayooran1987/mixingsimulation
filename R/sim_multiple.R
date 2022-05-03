##' This function calculates the resulting expected total number of colony-forming units in the mixed sample in the multiple mixing plans at the single stage of the mixing process.
##' @title The expected total number of colony-forming units in the mixed sample in the multiple mixing schemes at the single stage of the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions / primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param summary if we need to get all simulated \eqn{N'}, use \code{summary = 2}; otherwise, if we use \code{summary = 1}, the function provides the mean value of the simulated \eqn{N'}.
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
##' mu <- 100
##' sigma <- 0.8
##' alpha <- c(0.1,5)
##' k <- c(30,30)
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 2000
##' f_spri <- function(alpha, distribution) {
##' sprintf("mixing plan (alpha = %.1f, %s)", alpha, distribution)
##' }
##' sim.sum3 <- sim_multiple(mu, sigma, alpha, k, distribution, summary = 2, n_sim)
##' result <- data.frame(1:n_sim, sim.sum3)
##' colnames(result) <- c("n_sim", f_spri(alpha, distribution))
##' melten.Prob <- reshape2::melt(result, id = "n_sim", variable.name = "mixing_scheme",
##'                               value.name = "Total_CFU")
##' plot_example <-
##' ggplot2::ggplot(melten.Prob, ggplot2::aes(Total_CFU, group = mixing_scheme,colour = mixing_scheme))+
##' ggplot2::geom_line(stat="density",ggplot2::aes(x = Total_CFU))+
##' ggplot2::ylab(expression("pmf"))+
##' ggplot2::theme_classic()+ ggplot2::xlab(expression("Total number of CFU in the mixed sample"))+
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.75))+
##' ggthemes::scale_colour_colorblind()
##' plot_example
##' @export
sim_multiple <- function(mu, sigma, alpha, k, distribution, summary, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
    }
  if (length(alpha) != length(k)) {
    warning("length of alpha and length of k are must be equal")
    } else if (summary == 1) {
      sim.sum1 <- matrix(NA, nrow = 1, ncol = length(k))
      for (j in 1:length(k)) {
        sim.sum1[,j] <-  sim_single(mu, sigma, alpha[j], k[j], distribution, n_sim, summary = 1)
      }
    } else if (summary == 2) {
      sim.sum1 <- matrix(NA, nrow = n_sim, ncol = length(alpha))
      for (j in 1:length(alpha)) {
        sim.sum1[,j] <- sim_single(mu, sigma, alpha[j], k[j], distribution = "Poisson lognormal-Type B", n_sim, summary = 2)
      }
    } else {
      print("please include the summary value, which depends on your expected output")
    }
  result <- data.frame(sim.sum1)
  # result <- sim.sum1
  colnames(result) <- f_spri(mu, k, alpha, distribution)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}

