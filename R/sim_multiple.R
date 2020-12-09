##' This function calculates the resulting total number of colony forming units in the mixed sample in the multiple mixing plans. (to be finished later on)
##' @title The total number of colony-forming units in the mixed sample by simulation result (in the multiple mixing plan).
##' @param mu the average number of colony-forming units in the mixed sample, which is in logarithmic scale if we use a lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the log10 scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"}
##' @param n_sim number of simulations
##' @return total number of colony forming units in the multiple mixing scheme
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}. To more details, please refer the details section of  \link{compare_mixing}. (to be finished later on)
##' @seealso \link{sim_single}, \link{compare_mixing}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- c(100,100)
##' sigma <- c(0.8,0.8)
##' alpha <- c(0.1,10)
##' k <- c(10,10)
##' distribution <-  c("Lognormal-Type B","Lognormal-Type B")
##' n_sim <- 20000
##' n_sim_df <-data.frame(n_simulations = c(1:n_sim))
##' Prob_df <- cbind.data.frame(n_sim_df,sim_multiple(mu, sigma, alpha, k, distribution, n_sim))
##' melten.Prob <- reshape2::melt(Prob_df, id = "n_simulations", variable.name = "mixing_scheme",
##' value.name = "Total_CFU")
##' plot_example <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_point(ggplot2::aes(x = n_simulations,
##' y = Total_CFU, group = mixing_scheme, colour = mixing_scheme))+
##' ggplot2::xlab(expression("Number of simulations"))+ ggplot2::theme_classic()+
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),legend.position = c(0.70,0.90),
##' legend.box.background = ggplot2::element_rect(),legend.box.margin = ggplot2::margin(1,1,1,1))+
##' ggplot2::ylab(expression("Total number of CFU"))+
##' ggthemes::scale_colour_colorblind()
##' print(plot_example)
##' @export
sim_multiple <- function(mu, sigma, alpha, k, distribution, n_sim){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, alpha = %.1f, %s)", mu, k, alpha, distribution)
  }
  if (length(mu)!=length(k)) {
    warning("length of mu and length of k are must be equal")
  } else {
  sim.sum1 <- matrix(NA, nrow = n_sim, ncol = length(k))
  for(j in 1:length(k)){
    sim.sum1[,j] <-  sim_single(mu[j], sigma[j], alpha[j], k[j], distribution[j], n_sim)
  }
   }
  result <- data.frame(sim.sum1)
  colnames(result) <- f_spri(mu, k, alpha, distribution)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}

