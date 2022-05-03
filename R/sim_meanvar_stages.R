##' This function provides the mean and variance of the expected number of CFUs at each mixing stage.
##' @title Mean and variance of the expected number of CFUs at each mixing stage.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return Mean and variance changes at each mixing stage.
##' @details Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}.
##' This function produces a graphical display of the mean and variance changes at each mixing stage. It is helpful to identify the optimal number of revolutions of the mixture, which is a point of mixing that initiates Poisson-like homogeneity.
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- 2500
##' r <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 2000
##' result1 <- sim_meanvar_stages(mu, sigma , alpha_in, k, l, r, distribution, n_sim)
##' melten.Prob <- reshape2::melt(result1, id = "Revolutions",
##' variable.name = "summary", value.name = "Value")
##' plot_example <- ggplot2::ggplot(melten.Prob, ggplot2::aes(x = Revolutions,
##' y = Value, group = summary, colour = summary)) +
##' ggplot2::geom_line(ggplot2::aes(x = Revolutions, y = Value)) +
##' ggplot2::ylab(expression("No. of colony forming units")) +
##' ggplot2::theme_classic() + ggplot2::xlab(expression("No. of revolutions")) +
##' ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
##' legend.position = c(0.75,0.25),legend.title = ggplot2::element_blank()) +
##' ggthemes::scale_colour_colorblind()
##' plot_example
##' @export
sim_meanvar_stages <- function(mu, sigma , alpha_in, k, l, r, distribution, n_sim){
  "Revolutions" <- NULL
  "summary" <- NULL
  "Value" <- NULL
  # f_spri <- function(mu, k, alpha, distribution) {
  #   sprintf("mixing plan (mu = %.1f, k = %.0f, l = %.1f, %s)", mu, k, l, distribution)
  # }
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  stages <- 1:l
  alpha <- matrix(NA, nrow = l , ncol = 1)
  for (j in 1:l) {
    if (j == 1) {
      alpha[j,] <- alpha_in
    } else {
      alpha[j,] <- alpha[j - 1,] + r
    }
  }
  sim.sum1 <- matrix(NA, nrow = l, ncol = 2)
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  for (j in 1:l) {
    sim.sum1[j,] <- sim_meanvar(mu, sigma , alpha[j,1] , k, distribution, n_sim)
  }
  result <- data.frame(stages,sim.sum1)
  colnames(result) <- c("Revolutions", "Mean","Variance")
  # cat("Calculation took", proc.time()[1], "seconds.\n")
  return(result)
}

