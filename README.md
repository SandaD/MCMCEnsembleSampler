# MCMC Ensemble Sampler

Ensemble Markov Chain Monte Carlo samplers with different strategies to generate proposals.
Either the _stretch move_ as proposed by Goodman and Weare (2010),
or a _differential evolution jump move_ (similar to Braak and Vrugt,
2008) is used.


## Installation

1. Install the R package `devtools`
2. Then run from R
```R
library(devtools)
install_github("SandaD/MCMCEnsembleSampler")
```

## Usage

```R
library(MCMCEnsembleSampler)

## a log-pdf to sample from
p.log <- function(x) {
    B <- 0.03                              # controls 'bananacity'
    -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

## use stretch move
res1 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
                     max.iter=3000, n.walkers=10, method="s")
str(res1)


## use stretch move, return samples as 'coda' object
res2 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
                     max.iter=3000, n.walkers=10, method="s", coda=TRUE)

summary(res2$samples)
plot(res2$samples)


## use different evolution move, return samples as 'coda' object
res3 <- MCMCEnsemble(p.log, lower.inits=c(a=0, b=0), upper.inits=c(a=1, b=1),
                     max.iter=3000, n.walkers=10, method="d", coda=TRUE)

summary(res3$samples)
plot(res3$samples)
```

## References

Goodman, J. and Weare, J. (2010) Ensemble samplers with affine invariance. Communications in Applied Mathematics and Computational Science, 5(1), 65–80.

Braak, C. J. F. ter and Vrugt, J. A. (2008) Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing, 18(4), 435–446.
