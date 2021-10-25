##' This function gives a simulated number of CFU after each stage of the mixing process.
##' @title The total number of colony-forming units in the mixed sample by the simulation results in the single mixing plan with \eqn{l}  number of stages.
##' @param l the maximum number of stages in the mixing process
##' @param mu the average number of colony-forming units (\eqn{\mu}) in the mixed sample, which is in logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param summary if we need to get all simulated \eqn{N'}, use \code{summary = 3}; otherwise, if we use \code{summary = 1} or \code{summary = 2}, the function provides the mean value of the simulated \eqn{N'} or generated CFUs in each primary sample, respectively ( default \code{summary = 1}).
##' @param n_sim number of simulations
##' @return average number of colony forming units in the single mixing plan with \eqn{l} number of stages.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by contribution of \eqn{k} primary samples mixing, \eqn{N' = \sum N_i} and \eqn{l} be the number of stages in the mixing process.
##' This function provides simulated number of CFU after each stages of the mixing process. To more details, please refer the details section of  \link{compare_mixing_stages}.
##' @seealso \link{sim_single}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- 25000
##' rate <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 20000
##' no.revolutions <-c(1:l)
##' Prob_df <-
##' data.frame(no.revolutions,sim_single_stages(mu,sigma,alpha_in,k,l,rate,distribution,n_sim))
##' colnames(Prob_df) <- c("no.revolutions","CFU")
##' cummean <- function(x){cumsum(x)/seq_along(x)}
##' cum_mean <- cummean(Prob_df[,2])
##' plot_example <- ggplot2::ggplot(Prob_df) +
##' ggplot2::geom_line(ggplot2::aes(x = no.revolutions, y = CFU))+
##' ggplot2::geom_line( ggplot2::aes(x = no.revolutions, y = cum_mean), color = "red", size = .75)+
##' ggplot2::xlab(expression("Number of revolutions"))+
##' ggplot2::ylab(expression("Expected total number of CFU"))+
##' ggplot2::theme_classic()+
##' ggplot2::ggtitle(label = "Expected total number of CFU versus number of revolutions")+
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
##' ggthemes::scale_colour_colorblind()
##' print(plot_example)
##' @export
sim_single_stages <- function(mu, sigma , alpha_in, k, l, rate, distribution, n_sim, summary = 1){
  f_spri <- function(mu, k, alpha, distribution) {
    sprintf("mixing plan (mu = %.1f, k = %.0f, l = %.1f, %s)", mu, k, l, distribution)
  }
  stages <- 1:l
  alpha <- matrix(NA, nrow =1 , ncol = l)
  for(j in 1:l){
    if (j==1) {
      alpha[,j] <- alpha_in
    } else {
      alpha[,j] <- alpha[,j-1]+ rate
    }
  }
  if (summary == 1){
    sim.sum1 <- matrix(NA, nrow = l, ncol = 1)
    for(j in 1:length(alpha)){
      sim.sum1[j,] <- sim_single(mu, sigma, alpha[1,j], k, distribution, n_sim, summary = 1)
    }
    results <- sim.sum1
    colnames(results) <- f_spri(mu, k, l, distribution)
  } else if (summary == 2){
    sim.sum1 <- matrix(NA, nrow = l, ncol = k)
    for(j in 1:length(alpha)){
      sim.sum1[j,] <- sim_single(mu, sigma, alpha[,j], k, distribution, n_sim, summary = 2)
    }
    results <- sim.sum1
  } else if (summary == 3){
    sim.sum1 <- matrix(NA, nrow = l, ncol = n_sim)
    # result <- round(sum(apply(sim_single(mu, sigma , alpha , k, distribution, n_sim, summary = 3), 1, sum)))
    for(j in 1:length(alpha)){
      sim.sum1[j,] <- sim_single(mu, sigma, alpha[,j], k, distribution, n_sim, summary = 3)
    }
    results <- sim.sum1
  } else {
    print("please include the summary value, which depends on your expected output")
  }
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(results)

}

