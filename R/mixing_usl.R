##' This function provides a graphical display to find out upper stabilizing limit (USL). Based on the cumulative mean at each stage of the mixing,
##' the cumulative mean persists in the same value after a particular number of stages. Therefore, the upper stabilizing limit (USL) is defined as
##' a stabilizing point in the graphical display of the expected total number of CFUs versus the number of revolutions when the cumulative mean remains stable.
##' @title Graphical display to find out upper stabilizing limit (USL).
##' @param l the maximum number of stages in the mixing process
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param rate concentration parameter changing rate in the each revolutions
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return graphical display to find out upper stabilizing limit (USL).
##' @seealso \link{sim_single_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- 30
##' l <- 250
##' rate <- 0.01
##' distribution <-  "Poisson lognormal-Type B"
##' n_sim <- 20
##' mixing_usl(mu, sigma, alpha_in, k, l, rate, distribution, n_sim)
##' @export
mixing_usl <-  function(mu, sigma, alpha_in, k, l, rate, distribution, n_sim){
  no.revolutions <- NULL
  CFU <- NULL
  no.revolutions <- c(1:l)
  Prob_df <-
    data.frame(no.revolutions,sim_single_stages(mu,sigma,alpha_in,k,l,rate,distribution,n_sim))
  colnames(Prob_df) <- c("no.revolutions","CFU")
  cummean <- function(x){cumsum(x)/seq_along(x)}
  cum_mean <- round(cummean(Prob_df[,2]))
  #DescTools::Mode(cum_mean)
  Prob_df1 <- data.frame(Prob_df,cum_mean)
  plot1 <- ggplot2::ggplot(Prob_df1) +
    ggplot2::geom_line(ggplot2::aes(x = no.revolutions, y = CFU)) +
    ggplot2::geom_line( ggplot2::aes(x = no.revolutions, y = cum_mean), color = "red", size = .75) +
    ggplot2::xlab(expression("Number of revolutions")) +
    ggplot2::ylab(expression("Expected total number of CFUs")) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    # ggthemes::scale_colour_colorblind()+ggplot2::geom_vline(xintercept=no.revolutions[which.max(Prob_df$CFU)],linetype = "dashed",colour = "red")+
    ggthemes::scale_colour_colorblind() + ggplot2::geom_vline(xintercept = no.revolutions[length(cum_mean) - max(rle(cum_mean == DescTools::Mode(cum_mean))$lengths[rle(cum_mean == DescTools::Mode(cum_mean))$values == TRUE])],linetype = "dashed",colour = "blue") +
    ggplot2::annotate("text", x = no.revolutions[length(cum_mean) - max(rle(cum_mean == DescTools::Mode(cum_mean))$lengths[rle(cum_mean == DescTools::Mode(cum_mean))$values == TRUE])], y = 0,
                      label = sprintf("\n USL = %0.0f",round(no.revolutions[length(cum_mean) - max(rle(cum_mean == DescTools::Mode(cum_mean))$lengths[rle(cum_mean == DescTools::Mode(cum_mean))$values == TRUE])])), size = 3)
  #
  # ggplot2::annotate("text", x=1.2*no.revolutions[which.max(Prob_df$CFU)], y = 0,
  #                   label = sprintf("\n USL = %0.0f",round(no.revolutions[which.max(Prob_df$CFU)])), size=3)
  return(plot1)
}

