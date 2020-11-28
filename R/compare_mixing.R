##' This function provides the graphical displays for a different set of mixing plans for comparison purpose. (to be finished later on)
##' @title The graphical comparison between different mixing schemes with varying parameters of mixing by simulation results.
##' @param n_iter number of iterations
##' @param mu average number of colony-forming units in a primary sample which is in logarithmic scale if we use a lognormal distribution
##' @param sigma log standard deviation of the colony-forming units in a primary sample
##' @param b concentration parameter
##' @param k number of small portions/ primary samples
##' @param distribution what suitable distribution type we have employed for simulation such as 'Poisson-Type A' or 'Poisson-Type B' or 'Lognormal-Type A' or 'Lognormal-Type B'
##' @return graphical comparison between different mixing schemes
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary sample and \eqn{N' = \sum(N_i)}  and \eqn{i} be the number of colony-forming units in the \eqn{i^{th}} primary sample; where \eqn{i = 1,2,....k} and \eqn{y_i=x_i/\sum(x_i) = q_i/Q}; where \eqn{x_i} follows \eqn{gamma (b,1)} (to be finished later on)
##' \itemize{
##' \item Case 1 (Poisson-Type A): \eqn{N_i} follows \eqn{Poisson(\mu/k)}
##' \item Case 2 (Poisson-Type B): \eqn{N_i} follows \eqn{Poisson(\mu*y_i)}
##' \item Case 3 (Lognormal-Type A): \eqn{N_i} follows \eqn{Binomial(M_i,1/k)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' \item Case 4 (Lognormal-Type B): \eqn{N_i}  follows \eqn{Binomial(M_i,y_i)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' }
##' @seealso \link{sim_multiple}, \link{sim_single}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' n_iter <- 2000000
##' mu <- c(log(100),log(100),log(100))
##' sigma <- c(0.8,0.8,0.8)
##' b <- c(0.1,1,10)
##' k <- c(30,30,30)
##' distribution <-  c("Lognormal-Type B","Lognormal-Type B","Lognormal-Type B")
##' compare_mixing (n_iter, mu, sigma, b, k, distribution )
##' @export
compare_mixing <- function(n_iter, mu, sigma, b, k, distribution){
  Total_CFU <- NULL
  mixing_scheme <- NULL
  f_spr <- function(n_iter) {
    sprintf("Simulation results (no.iterations = %.0f)", n_iter)
  }
  n_iterations <-data.frame(n_iterations = c(1:n_iter))
  Prob_df <- cbind.data.frame(n_iterations, sim_multiple( n_iter, mu, sigma, b, k, distribution))
  melten.Prob <- reshape2::melt(Prob_df, id = "n_iterations", variable.name = "mixing_scheme", value.name = "Total_CFU")
  # ggplot2::ggplot(df, aes(Total_CFU)) + stat_ecdf(geom = "point")+ theme_classic()
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(Total_CFU, group = mixing_scheme, colour = mixing_scheme)) + ggplot2::stat_ecdf(geom = "step")+ggplot2::ylab(expression("Cumulative probability"))+ ggplot2::theme_classic()+ ggplot2::xlab(expression("N' after mixing"))+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.25)) + ggplot2::ggtitle(label = f_spr(n_iter))+ ggthemes::scale_colour_colorblind()
  # + ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::ggtitle(label = f_spr(mu, sigma, b, k, n_iter))
  # print(plot1)
  return(plot1)
  # cat("Calculation took", proc.time()[1], "seconds.\n")
}



