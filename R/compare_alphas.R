##' This function provides a graphical display to investigate the impact of concentration parameter values in the Dirichlet distribution.
##' @title Graphical display to investigate the impact of concentration parameter values in the Dirichlet distribution.
##' @param n number of weights which is generated from Dirichlet distribution
##' @param alpha concentration parameter is given in vector form; each element represents each group.
##' @return graphical display of different concentration parameters of Dirichlet distribution
##' @details We assumed that the concentration parameters are identical in all weights, it means that \eqn{\boldsymbol{\alpha} =[\alpha, \alpha, \cdots \alpha]}. The probability mass
##' function of the weights is given by: \deqn{\displaystyle{f(w_1,w_2, \cdots w_n;\boldsymbol{\alpha} )=\frac{\displaystyle{\Gamma(n\alpha)}} {{\Gamma(\alpha)}^n} \prod_{i=1}^{n}{w_i^{\alpha-1}}}}
##' where \eqn{\boldsymbol{\alpha} =[\alpha, \alpha, \cdots \alpha]} is the vector of concentration parameter, and \eqn{w_i} is the stochastic weights which sum to one.
##' @examples
##' n <- 10000
##' alpha <- c(5,10,20)
##' compare_alphas(n,alpha)
##' @export
compare_alphas <- function(n,alpha){
  concentration_parameter <- NULL
  values <- NULL
  rdirichlet <- function(n,alpha1){
    x <-  matrix(NA, nrow = 1, ncol = n) # If we want to apply a beta algorithm to generate Dirichlet distribution's random numbers.
    for (j in 1:n) {
      x[,j] <- stats::rbeta(1,alpha1, alpha1*(n - j))
    }
    w <-  matrix(NA, nrow = 1, ncol = n)
    for (j in 2:n) {
      w[,1] <- x[,1]
      w[,j] <- x[,j] %*% prod(1 - x[,1:(j - 1)])
    }
    #sum(w)
    return(w)
  }

  f_spri <- function(alpha) {
    sprintf("alpha = %.1f", alpha)
  }
  sim.sum1 <- matrix(NA, nrow = n, ncol = length(alpha))
  for (j in 1:length(alpha)) {
    sim.sum1[,j] <- rdirichlet(n,alpha[j])
  }

  result <- data.frame(1:n, sim.sum1)
  colnames(result) <- c("weight", f_spri(alpha))
  melten.Prob <- reshape2::melt(result, id = "weight", variable.name = "concentration_parameter",
                                value.name = "values")
  plot_example <-
    ggplot2::ggplot(melten.Prob, ggplot2::aes(values, group = concentration_parameter,colour = concentration_parameter)) +
    ggplot2::geom_line(stat = "density",ggplot2::aes(x = values)) +
    ggplot2::ylab(expression("probability density")) +
    ggplot2::theme_classic() + ggplot2::xlab(expression("weights")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75,0.75)) +
    ggthemes::scale_colour_colorblind()
  return(plot_example)
}


