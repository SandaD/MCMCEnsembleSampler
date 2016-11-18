\name{d.e.mcmc}
\alias{d.e.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To be added
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
d.e.mcmc(f, max.iter, n.walkers, n.dim, init.range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{n.walkers}{
%%     ~~Describe \code{n.walkers} here~~
}
  \item{n.dim}{
%%     ~~Describe \code{n.dim} here~~
}
  \item{init.range}{
%%     ~~Describe \code{init.range} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, max.iter, n.walkers, n.dim, init.range) 
{
    chain.length = max.iter\%/\%n.walkers
    sum.log.p.d.e <<- rep(0, chain.length)
    log.p = matrix(NA, nrow = n.walkers, ncol = chain.length)
    log.p.old = rep(NA, n.walkers)
    ensemble.old = matrix(NA, nrow = n.walkers, ncol = n.dim)
    ensemble.new = matrix(NA, nrow = n.walkers, ncol = n.dim)
    x.chain = array(NA, dim = c(n.walkers, chain.length, n.dim))
    for (k in 1:n.walkers) {
        for (g in 1:n.dim) {
            ensemble.old[k, g] = runif(1, init.range[1], init.range[2])
        }
        log.p.old[k] = f(ensemble.old[k, ])
    }
    log.p[, 1] = log.p.old
    sum.log.p.d.e[1] <<- sum(log.p.old[1:n.walkers])/n.walkers
    x.chain[, 1, ] = ensemble.old
    for (l in 2:chain.length) {
        for (n in 1:n.walkers) {
            z = 2.38/sqrt(2 * n.dim)
            if (l\%\%10 == 0) {
                z = 1
            }
            a = sample((1:n.walkers)[-n], 1)
            b = sample((1:n.walkers)[-c(n, a)], 1)
            par.active.1 = ensemble.old[a, ]
            par.active.2 = ensemble.old[b, ]
            ensemble.new[n, ] = ensemble.old[n, ] + z * (par.active.1 - 
                par.active.2)
            log.p.new = f(ensemble.new[n, ])
            acc = exp(log.p.new - log.p.old[n])
            test = runif(1)
            if (acc > test) {
                x.chain[n, l, ] = ensemble.new[n, ]
                ensemble.old[n, ] = ensemble.new[n, ]
                log.p[n, l] = log.p.new
                log.p.old[n] = log.p.new
            }
            else {
                x.chain[n, l, ] = ensemble.old[n, ]
                log.p[n, l] = log.p.old[n]
            }
            sum.log.p.d.e[l] <<- sum.log.p.d.e[l] + log.p[n, 
                l]/n.walkers
        }
    }
    return(x.chain)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
