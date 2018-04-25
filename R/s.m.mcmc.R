## ...........................
## MCMC ensemble sampler
## ................
## Sanda Dejanic - R implementation of Goodman and Ware's emcee
## .............................


##' MCMC Ensemble sampler with the stretch move (emcee)
##'
##' Markov Chain Monte Carlo sampler: using the stretch move (implementation of the Goodman and Ware emcee)
##'
##' @param f a funtion to sample from
##' @param max.iter the maximum number of function evaluations
##' @param n.walkers the number of walkers (ensemble size)
##' @param n.dim the number of parameters to sample
##' @param init.range a matrix(nrow=n.dim, ncol=2) defying the initial range for all the parameters, every row of the matrix should contain the lower and the upper limit
##' @param ... all additional agruments of f
##'
##' @return List containing: \code{$samples[n.walkers,chain.length,n.dim] and $log.p[n.walkers,chain.length]}
##' @export
s.m.mcmc = function(f, max.iter, n.walkers, n.dim, init.range, ...) {


    ## initial values


    chain.length = max.iter%/%n.walkers

    log.p = matrix(NA, nrow=n.walkers, ncol=chain.length)
    log.p.old = rep(NA, n.walkers)
    ensemble.old = matrix(NA, nrow=n.walkers, ncol=n.dim)
    ensemble.new = matrix(NA, nrow=n.walkers, ncol=n.dim)
    samples = array(NA, dim=c(n.walkers,chain.length,n.dim))
    mcmc.object = array(NA, dim=c(n.walkers,chain.length,n.dim+1))


    for(k in 1:n.walkers) {
        for(g in 1:n.dim){
            ensemble.old[k,g] = runif(1, init.range[g,1], init.range[g,2])
        }
        log.p.old[k] = f(ensemble.old[k,], ...)
    }


    log.p[,1]=log.p.old
    samples[ , 1, ] = ensemble.old


    ## the loop

    for (l in 2:chain.length){
        for (n in 1:n.walkers){


            z = ((runif(1)+1)^2)/2
            a = sample((1:n.walkers)[-n], 1)
            par.active = ensemble.old[a,]


            ensemble.new[n,] = par.active + z*(ensemble.old[n,] - par.active)


            log.p.new = f(ensemble.new[n,], ...)
            if(!is.finite(log.p.new)){acc=0}
            else {acc = z^(n.dim-1)*exp(log.p.new - log.p.old[n])}

            test = runif(1)

            if (acc > test) {

                samples[n,l,] = ensemble.new[n,]
                ensemble.old[n,] = ensemble.new[n,]
                log.p[n,l] = log.p.new
                log.p.old[n] = log.p.new


            } else {

                samples[n,l,] = ensemble.old[n,]
                log.p[n,l] = log.p.old[n]

            }
        }
    }

    mcmc.list = list(samples=samples, log.p=log.p)

    return(mcmc.list)

}
