% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s.m.mcmc.R
\name{s.m.mcmc}
\alias{s.m.mcmc}
\title{MCMC Ensemble sampler with the stretch move (emcee)}
\usage{
s.m.mcmc(f, max.iter, n.walkers, n.dim, init.range, ...)
}
\arguments{
\item{f}{a funtion to sample from}

\item{max.iter}{the maximum number of function evaluations}

\item{n.walkers}{the number of walkers (ensemble size)}

\item{n.dim}{the number of parameters to sample}

\item{init.range}{a matrix(nrow=n.dim, ncol=2) defying the initial range for all the parameters, every row of the matrix should contain the lower and the upper limit}

\item{...}{all additional agruments of f}
}
\value{
List containing: \code{$samples[n.walkers,chain.length,n.dim] and $log.p[n.walkers,chain.length]}
}
\description{
Markov Chain Monte Carlo sampler: using the stretch move (implementation of the Goodman and Ware emcee)
}
