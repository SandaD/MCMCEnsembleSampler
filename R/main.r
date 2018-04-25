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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param f 
##' @param lower.inits 
##' @param upper.inits 
##' @param max.iter 
##' @param n.walkers 
##' @param method 
##' @param coda 
##' @param ... 
##' @return 
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
        res <- d.m.mcmc(f, max.iter, n.walkers, n.dim, init.range, ...)
    }
    if(method == "stretch") {
        print(paste("Using stretch move with", n.walkers, "walkers."))
        res <- s.m.mcmc(f, max.iter, n.walkers, n.dim, init.range, ...)
    }

    ## add names
    dimnames(res$samples) <- list(paste0("walker_", 1:n.walkers),
                                  paste0("generation_", 1:dim(res$samples)[2]),
                                  names(lower.inits))

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
