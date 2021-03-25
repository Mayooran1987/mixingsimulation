##' This function provides the graphical displays for a different set of mixing parameters for comparison purpose of mixing schemes with multiple stages of mixing.
##' @title The graphical comparison between different mixing schemes by the simulation results in the mixing process's multiple stages with varying mixing parameters.
##' @param mu the average number of colony-forming units in the mixed sample, which is in logarithmic scale if we use a lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions/ primary samples
##' @param l number of revolutions/stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"}
##' @param n_sim number of simulations
##' @return graphical comparison between different mixing schemes.
##' @details {Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i} and \eqn{N_i} be the number of colony-forming units
##' in the \eqn{i^{th}} primary sample; where \eqn{i = 1,2,....k}.
##'
##'
##' For this package development, we have employed the notations 'Type-A' and 'Type-B'  to indicate the type of distributions, which are applied in the previous literature as 'fair' and 'beta', respectively; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##'
##' Following \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}, the contribution weight of contamination by each primary sample can be defined by the random variable \eqn{w_i}
##' which is possible to be following either uniform distribution with parameter \eqn{1/k} or joint distribution of \eqn{w_1,w_2,\cdots w_k} follows Dirichlet distribution with concentration parameter \eqn{\alpha}.
##' From the previous literature, Dirichlet distribution can be formulated by beta or gamma algorithm which are revealed the same results; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##'
##' This function is developed based on the beta algorithm, and the following steps formulate it.
##' \deqn{w_i = x_i {\prod_{j=1}^{i-1}{1-x_j}}~~~~ \forall i = 2,3,\cdots k} and \eqn{w_1=x_1};
##'
##' where \eqn{x_i} follows \eqn{Beta (\alpha,\alpha(k-i))} and also \eqn{\sum w_i} must be equal to one.
##' \itemize{
##' \item Case 1 (Poisson-Type A): \eqn{N_i} follows \eqn{Poisson(\mu/k)}
##' \item Case 2 (Poisson-Type B): \eqn{N_i} follows \eqn{Poisson(\mu*w_i)}
##' \item Case 3 (Lognormal-Type A): \eqn{N_i} follows \eqn{Binomial(M_i,1/k)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' \item Case 4 (Lognormal-Type B): \eqn{N_i}  follows \eqn{Binomial(M_i,w_i)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' }}
##'
##'
##' The powder mixing process can be defined as breaking clusters stage by stage. Usually, it will be occurring systematically in the standard powder mixtures. For this package development, we assume that mixing parameters also systematically changing with a fixed rate in each stage of the mixing.
##' The mixing parameter can be defined as revolutions instead of the mixing stage in general. Due to the lack of theoretical results to the dependent random variable sum's distribution, we have chosen simulation techniques for this modelling.
##'
##' Let \eqn{l} be the number of stages or revolution of the mixture, also we assumed a fixed concentration parameter value at the initial phase of the mixing process. Based on the literature in this area, the concentration parameter can be assumed that increasing at every stage of the mixing, which is possible to be systematically. Therefore this function exhibits the graphical display with different quantities of primary sample mixing as a large unit.
##' @seealso  \link{sim_single}, \link{sim_single_stages}, \link{sim_multiple_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 1
##' k <- c(10,30,50)
##' l <- 1000
##' rate <- 0.01
##' distribution <- c("Lognormal-Type B","Lognormal-Type B","Lognormal-Type B")
##' n_sim <- 200000
##' compare_mixing_stages(mu, sigma, alpha_in, k, l, rate, distribution, n_sim)
##' @export
compare_mixing_stages <- function(mu, sigma, alpha_in, k, l, rate, distribution, n_sim){
  Total_CFU <- NULL
  mixing_scheme <- NULL
  f_spri <- function(mu, k, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, %s)", mu, k, distribution)
  }
  f_spr <- function(l) {
    sprintf("Simulation results (no.revolutions = %.0f)", l)
  }
  stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = l, ncol = length(distribution))
  for(j in 1:length(distribution)){
    sim.sum3[,j] <-  sim_single_stages(mu, sigma , alpha_in, k[j], l, rate, distribution[j], n_sim)
    }
  result <- data.frame(stages, sim.sum3)
  colnames(result) <- c("stages", f_spri(mu, k, distribution))
  # return(result)
  melten.Prob <- reshape2::melt(result, id = "stages", variable.name = "mixing_scheme", value.name = "Total_CFU")
  # ggplot2::ggplot(df, aes(Total_CFU)) + stat_ecdf(geom = "point")+ theme_classic()
  plot1 <- ggplot2::ggplot(melten.Prob, ggplot2::aes(Total_CFU, group = mixing_scheme, colour = mixing_scheme)) + ggplot2::stat_ecdf(geom = "step") + ggplot2::ylab(expression("Cumulative probability of N'"))+
    ggplot2::theme_classic()+ ggplot2::xlab(expression("Total CFU after mixing (N')"))+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.70,0.25))+
    ggplot2::ggtitle(label = f_spr(l))+ ggthemes::scale_colour_colorblind()
  # + ggplot2::theme(plot.title = element_text(hjust = 0.5))
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(plot1)
  }
