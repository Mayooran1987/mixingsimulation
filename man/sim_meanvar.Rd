% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_meanvar.R
\name{sim_meanvar}
\alias{sim_meanvar}
\title{Mean and variance of the expected number of CFUs in the single mixing stage.}
\usage{
sim_meanvar(mu, sigma, alpha, k, distribution, n_sim)
}
\arguments{
\item{mu}{the average number of CFUs (\eqn{\mu}) in the mixed sample, which is in a logarithmic scale if we use a Lognormal / Poisson lognormal distribution}

\item{sigma}{the standard deviation of the colony-forming units in the mixed sample on the logarithmic scale (default value 0.8)}

\item{alpha}{concentration parameter}

\item{k}{number of small portions / primary samples}

\item{distribution}{what suitable distribution type we have employed for simulation such as \code{"Poisson-Type A"} or \code{"Poisson-Type B"} or \code{"Lognormal-Type A"} or \code{"Lognormal-Type B"} or \code{"Poisson lognormal-Type A"} or \code{"Poisson lognormal-Type B"}}

\item{n_sim}{number of simulations}
}
\value{
Mean and variance changes in the single mixing stage.
}
\description{
This function provides the mean and variance of the expected number of CFUs in the single mixing stage.
}
\details{
Let \eqn{N'} be the number of colony-forming units in the mixed sample which is produced by mixing of \eqn{k} primary samples and \eqn{N' = \sum N_i}.
This function produces a graphical display of the mean and variance changes at each mixing stage. It is helpful to identify the optimal number of revolutions of the mixture, which is a point of mixing that initiates Poisson-like homogeneity.
}
\examples{
mu <- 100
sigma <- 0.8
alpha <- 0.1
k <- 30
distribution <-  "Poisson lognormal-Type B"
n_sim <- 2000
sim_meanvar(mu, sigma , alpha , k, distribution, n_sim)
}
