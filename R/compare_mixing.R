##' This function provides the graphical displays for a different set of mixing parameters for comparison purpose of mixing schemes.
##' @title The graphical comparison between different mixing schemes with varying parameters of mixing by simulation results.
##' @param mu the average number of colony-forming units in the mixed sample, which is in logarithmic scale if we use a lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the log10 scale (default value 0.8)
##' @param alpha concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"}
##' @param n_sim number of simulations
##' @return graphical comparison between different mixing schemes
##' @details {Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i} and \eqn{N_i} be the number of colony-forming units
##' in the \eqn{i^{th}} primary sample; where \eqn{i = 1,2,....k}. Following \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}, contribution weight of contamination by each primary sample can be defined by the random variable \eqn{w_i}
##' which is possible to be following either uniform distribution with paramater \eqn{1/k} or joint distribution of \eqn{w_1,w_2,\cdots w_k} follows Dirichlet distribution with concentration parameter \eqn{alpha}.
##' From the previous literature, Dirichlet distribution can be formulated by beta or gamma algorithm which are revealed the same results; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}. This function is developed by based on
##' gamma algorithm, it is formulated  the following steps.
##' \deqn{w_i=\frac{x_i}{\sum_{j=1}^{k}{x_j}}~~~~ \forall i = 1,2,\cdots k}; where \eqn{x_i} follows \eqn{gamma (alpha,1)} and also \eqn{\sum w_i} must be equal to one.
##' \itemize{
##' \item Case 1 (Poisson-Type A): \eqn{N_i} follows \eqn{Poisson(\mu/k)}
##' \item Case 2 (Poisson-Type B): \eqn{N_i} follows \eqn{Poisson(\mu*w_i)}
##' \item Case 3 (Lognormal-Type A): \eqn{N_i} follows \eqn{Binomial(M_i,1/k)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' \item Case 4 (Lognormal-Type B): \eqn{N_i}  follows \eqn{Binomial(M_i,w_i)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' }}
##' For this package development, we have employed the notations 'Type-A' and 'Type-B'  to indicate the type of distributions, which are applied in the previous literature as 'fair' and 'beta', respectively; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##' @seealso \link{sim_multiple}, \link{sim_single}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- c(100,100,100)
##' sigma <- c(0.8,0.8,0.8)
##' alpha <- c(0.1,1,10)
##' k <- c(30,30,30)
##' distribution <-  c("Lognormal-Type B","Lognormal-Type B","Lognormal-Type B")
##' n_sim <- 20000
##' compare_mixing (mu, sigma, alpha, k, distribution, n_sim )
##' @export
compare_mixing <- function(mu, sigma, alpha, k, distribution, n_sim){
  Total_CFU <- NULL
  mixing_scheme <- NULL
  f_spr <- function(n_sim) {
    sprintf("Simulation results (no.simulations = %.0f)", n_sim)
  }
  n_sim_df <-data.frame(n_simulations = c(1:n_sim))
  Prob_df <- cbind.data.frame(n_sim_df, sim_multiple(mu, sigma, alpha, k, distribution, n_sim))
  melten.Prob <- reshape2::melt(Prob_df, id = "n_simulations", variable.name = "mixing_scheme", value.name = "Total_CFU")
  # ggplot2::ggplot(df, aes(Total_CFU)) + stat_ecdf(geom = "point")+ theme_classic()
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(Total_CFU, group = mixing_scheme, colour = mixing_scheme)) + ggplot2::stat_ecdf(geom = "step")+ggplot2::ylab(expression("Cumulative probability of N'"))+ ggplot2::theme_classic()+ ggplot2::xlab(expression("N' after mixing"))+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.70,0.25)) + ggplot2::ggtitle(label = f_spr(n_sim))+ ggthemes::scale_colour_colorblind()
  # + ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::ggtitle(label = f_spr(mu, sigma, alpha, k, n_sim))
  # print(plot1)
  return(plot1)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
}



