##' This function provides a graphical display to compare mixing plans based on the cumulative distribution of expected total CFUs in the mixing process using different mixing parameters, such as type of distribution and number of primary samples.
##' @title Graphical comparison of mixing plans based on cumulative distribution of expected total CFUs in the mixing process.
##' @param mu the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution
##' @param sigma the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)
##' @param alpha_in concentration parameter at the initial stage
##' @param k number of small portions / primary samples
##' @param l number of revolutions / stages
##' @param r the rate of the concentration parameter changes at each mixing stage
##' @param distribution what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}
##' @param n_sim number of simulations
##' @return Graphical comparison between different mixing schemes.
##' @details {Let \eqn{N'} be the number of CFUs in the mixed sample, which is produced by the mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}  and let \eqn{N_i} be the number of CFUs in the \eqn{i^{th}} primary sample;
##'where \eqn{i = 1,2,....k}.
##'
##'For this package development, we have employed the notations 'Type-A' and 'Type-B' to indicate the type of distributions, which are applied in the previous literature as 'fair' and 'beta', respectively;
##'see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##'
##'
##'Following \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}, the contribution weight of contamination by each primary sample can be defined by the random variable \eqn{w_i},
##'which is possible to be followed by either uniform distribution with parameter \eqn{1/k} or the joint distribution of \eqn{w_1,w_2,\cdots w_k} follows a Dirichlet distribution with concentration parameter \eqn{\alpha}.
##'From the previous literature, a Dirichlet distribution can be formulated by beta or gamma algorithms, which have revealed the same results; see \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{Nauta (2005)}.
##'
##' This function is developed based on the beta algorithm and the following steps formulate it.
##' \deqn{w_i = x_i {\prod_{j=1}^{i-1}{1-x_j}}~~~~ \forall i = 2,3,\cdots k,~~~~~ w_1=x_1};
##'
##' where \eqn{x_i} follows \eqn{Beta (\alpha,\alpha(k-i))} and also \eqn{\sum w_i} must be equal to one.
##' \itemize{
##' \item Case 1 (Poisson-Type A): \eqn{N_i} follows \eqn{Poisson(\mu/k)}
##' \item Case 2 (Poisson-Type B): \eqn{N_i} follows \eqn{Poisson(\mu w_i)}
##' \item Case 3 (Lognormal-Type A): \eqn{N_i} follows \eqn{Binomial(M_i,1/k)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' \item Case 4 (Lognormal-Type B): \eqn{N_i}  follows \eqn{Binomial(M_i,w_i)}; where  \eqn{M_i} follows \eqn{Lognormal(\mu, \sigma)}
##' \item Case 5 (Poisson lognormal-Type A): \eqn{N_i} follows \eqn{Binomial(M_i,1/k)}; where  \eqn{M_i}
##'
##' follows \eqn{Poisson lognormal (\mu, \sigma)}
##' \item Case 6 (Poisson lognormal-Type B): \eqn{N_i}  follows \eqn{Binomial(M_i,w_i)}; where  \eqn{M_i}
##'
##' follows \eqn{Poisson lognormal(\mu, \sigma)}
##' }}
##'
##'The powder-mixing process can be defined as breaking clusters stage-by-stage. Usually, it occurs systematically in the standard powder mixtures. For this package development, we assume that mixing parameters also systematically change with a fixed rate at each stage of the mixing.
##'The mixing parameter can be defined as revolutions instead of the mixing stage in general. Due to the lack of theoretical results for the dependent random variable sum's distribution, we have chosen simulation techniques for this modelling.
##'
##'Let \eqn{l} be the number of stages or revolutions of the mixture, and we also assumed a fixed concentration parameter value at the initial phase of the mixing process. Based on the literature in this area, the concentration parameter can be assumed to increase at every stage of the mixing, which is possible to do systematically.
##'
##' Therefore, this function exhibits the graphical display with different quantities of primary sample mixing as a large unit.
##' @seealso  \link{sim_single}, \link{sim_single_stages}, \link{sim_multiple_stages}
##' @references
##' \itemize{
##' \item Nauta, M.J., 2005. Microbiological risk assessment models for partitioning and mixing during food handling. International Journal of Food Microbiology 100, \href{https://doi.org/10.1016/j.ijfoodmicro.2004.10.027}{311-322}.
##' }
##' @examples
##' mu <- 100
##' sigma <- 0.8
##' alpha_in <- 0.01
##' k <- c(10,30,60)
##' r <- 0.01
##' distribution <- c("Poisson lognormal-Type B","Poisson lognormal-Type B","Poisson lognormal-Type B")
##' n_sim <- 20000
##' plot1 <- compare_mixing_3(mu, sigma, alpha_in, k , l = 50,r, distribution,n_sim) +
##' ggplot2::theme(legend.text = ggplot2::element_text(size = 7.5),
##' legend.title = ggplot2::element_text(size = 7.5),
##' legend.key.size = ggplot2::unit(4, 'mm')) + ggplot2::xlim(0,300)
##' plot2 <- compare_mixing_3(mu, sigma, alpha_in, k , l = 500, r, distribution, n_sim) +
##' ggplot2::theme(legend.title = ggplot2::element_text(size = 7.5),
##' legend.key.size = ggplot2::unit(4, 'mm')) + ggplot2::xlim(0,300)
##' plot3 <- compare_mixing_3(mu, sigma, alpha_in, k, l = 5000,r , distribution,n_sim) +
##' ggplot2::theme(legend.text = ggplot2::element_text(size = 7.5),
##' legend.title = ggplot2::element_text(size = 7.5),
##' legend.key.size = ggplot2::unit(4, 'mm')) + ggplot2::xlim(0,300)
##' plot4 <- compare_mixing_3(mu, sigma , alpha_in , k , l = 25000, r , distribution ,n_sim) +
##' ggplot2::theme(legend.text = ggplot2::element_text(size = 7.5),
##' legend.title = ggplot2::element_text(size = 7.5),
##' legend.key.size = ggplot2::unit(4, 'mm')) + ggplot2::xlim(0,300)
##' gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)
##' @export
compare_mixing_3 <- function(mu, sigma, alpha_in, k, l, r, distribution, n_sim){
  "Total_CFU" <- NULL
  mixing_scheme <- NULL
  x <- NULL
  cd <- NULL
  set.seed(1, kind = "L'Ecuyer-CMRG")
  f_spri <- function(l, k, distribution) {
    sprintf("mixing plan (l = %.0f, k = %.0f, %s)",l, k, distribution)
  }
  # f_spr <- function(l) {
  #   sprintf("Simulation results (no.revolutions = %.0f)", l)
  # }
  stages <- 1:l
  sim.sum3 <- matrix(NA, nrow = l, ncol = length(distribution))
  # set.seed(1, kind = "L'Ecuyer-CMRG")
  for (j in 1:length(distribution)) {
    sim.sum3[,j] <-  sim_single_stages(mu, sigma , alpha_in, k[j], l, r, distribution[j], n_sim, summary = 1)
  }
  result <- data.frame(stages, sim.sum3)
  colnames(result) <- c("stages", f_spri(l, k, distribution))
  # return(result)
  melten.Prob <- reshape2::melt(result, id = "stages", variable.name = "mixing_scheme", value.name = "Total_CFU")
  # ggplot2::ggplot(df, aes(Total_CFU)) + stat_ecdf(geom = "point")+ theme_classic()
  dat <- cbind.data.frame(melten.Prob[,3],melten.Prob[,2])
  colnames(dat) <- c("x", "mixing_scheme")
  # Split the data by group and calculate the smoothed cumulative density for each group
  `%>%` <- magrittr::`%>%`
  dens <-  split(dat, dat$mixing_scheme) %>%
    purrr::map_df(function(d) {
      dens <-  stats::density(d$x, adjust = 0.1, from = min(dat$x) ,
                     to = max(dat$x))
      data.frame(x = dens$x, y = dens$y, cd = cumsum(dens$y)/sum(dens$y), mixing_scheme = d$mixing_scheme[1])
    })
  plot1 <- ggplot2::ggplot() +
     ggplot2::geom_line(data = dens, ggplot2::aes(x, cd, colour = mixing_scheme)) +
    # ggplot2::stat_smooth(data = dens,size = 0.5, method = 'gam', formula = y ~ s(x, bs = "cs") , ggplot2::aes(x, cd, colour = mixing_scheme),se = FALSE) +
     # ggplot2::ylim(0,1) +
    ggplot2::ylab(expression("P(Total CFUs after mixing " <= " N')")) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("Total CFUs after mixing (N')")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.70,0.25)) +
    # ggplot2::ggtitle(label = f_spr(l))+
    # ggplot2::scale_fill_discrete(name="mixing scheme")+
    ggthemes::scale_colour_colorblind()
  return(plot1)
}


