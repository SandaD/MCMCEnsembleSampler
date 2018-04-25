## -------------------------------------------------------
##
## File: main.r
##
## April 24, 2018 -- Andreas Scheidegger
## andreas.scheidegger@eawag.ch
## -------------------------------------------------------

##' @import stats
##' @import coda
NULL



##' MCMCEnsembleSampler
##'
##' This package implements a particle Monte Carlo Markov Chain sampler
##' with two different ways of creating proposals. 
##' @name MCMCEnsembleSampler
##' @author Sanda Kern-Dejanic
##' @docType package
NULL

##' MCMC ensemble sampler
##'
##' Ensemble Markov Chain Monte Carlo sampler with different strategies to generate proposals.
##' Either the \emph{stretch move} as proposed by Goodman and Weare (2010),
##' or a \emph{differential evolution jump move} similar to Braak and Vrugt (2008).
##' 
##' @param f function that returns a value proportional to the log probability
##' density to sample from. 
##' @param lower.inits vector specifying for each parameters the lower value the initial distibution
##' @param upper.inits vector specifying for each parameters the upper value the initial distibution
##' @param max.iter maximum number of function evaluations
##' @param n.walkers number of walkers (ensemble size)
##' @param method method for proposal generation, either \code{"stretch"}, or
##' \code{"differential.evolution"}. The first letter is sufficient.
##' @param coda logical. Should the samples be returned as \code{link[coda]{mcmc.list}} object?
##' @param ... further arguments passed to \code{f}
##' @return if \code{coda==FALSE} a list with
##' \itemize{
##'  \item{samples }{A three dimensional array of samples with dimensions \code{walker} x \code{generation} x \code{parameter}}
##'  \item{log.p }{A matrix with the log density evaluate for aeach walker at each generation.}
##' }
##' if \code{coda==TRUE} a list with
##' \itemize{
##'  \item{samples }{A object of class \code{link[coda]{mcmc.list}} containing all samples.}
##'  \item{log.p }{A matrix with the log density evaluate for aeach walker at each generation.}
##' }
##' @references
##' Goodman, J. and Weare, J. (2010) Ensemble samplers with affine invariance. Communications in Applied Mathematics and Computational Science, 5(1), 65-80.
##' 
##' Braak, C. J. F. ter and Vrugt, J. A. (2008) Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing, 18(4), 435-446.
##' @examples
##' ## a log-pdf to sample from
##' p.log <- function(x) {
##'     B <- 0.03                              # controls 'bananacity'
##'     -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
##' }
##' 
##' ## use stretch move
##' res1 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
##'                      max.iter=3000, n.walkers=10, method="s")
##' str(res1)
##' 
##' 
##' ## use stretch move, return samples as 'coda' object
##' res2 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
##'                      max.iter=3000, n.walkers=10, method="s", coda=TRUE)
##' 
##' summary(res2$samples)
##' plot(res2$samples)
##' 
##' 
##' ## use different evolution move, return samples as 'coda' object
##' res3 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
##'                      max.iter=3000, n.walkers=10, method="d", coda=TRUE)
##' 
##' summary(res3$samples)
##' plot(res3$samples)
##' @export
MCMCEnsemble <- function(f, lower.inits, upper.inits,
                         max.iter, n.walkers=10*length(lower.inits),
                         method=c("stretch", "differential.evolution"),
                         coda=FALSE, ...){

    if(length(lower.inits) != length(upper.inits)) {
        stop("The length of 'lower.inits' and 'lower.inits' is must be identical!")   
    }
    
    n.dim <- length(lower.inits)
    init.range <- cbind(lower.inits, upper.inits)

    ## run mcmc
    method <- match.arg(method)
    if(method == "differential.evolution"){
        print(paste("Using differential evolution move with", n.walkers, "walkers."))
        res <- d.e.mcmc(f, max.iter, n.walkers, n.dim, init.range, ...)
    }
    if(method == "stretch") {
        print(paste("Using stretch move with", n.walkers, "walkers."))
        res <- s.m.mcmc(f, max.iter, n.walkers, n.dim, init.range, ...)
    }

    ## add names
    if(is.null(names(lower.inits))){
        pnames <- paste0("para_", 1:n.dim)
    } else {
        pnames <- names(lower.inits)
    }
    dimnames(res$samples) <- list(paste0("walker_", 1:n.walkers),
                                  paste0("generation_", 1:dim(res$samples)[2]),
                                  pnames)

    dimnames(res$log.p) <- list(paste0("walker_", 1:n.walkers),
                                paste0("generation_", 1:dim(res$samples)[2]))
    
    ## convert to coda object
    if(coda){
        ll <- list()
        for(w in 1:n.walkers){
            ll[[w]] <- coda::as.mcmc(res$samples[w,,])
        }
        res <- list(samples=as.mcmc.list(ll), log.p=res$log.p)
    }
    
    res

}
