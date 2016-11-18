s.m.mcmc <-
function(f, max.iter, n.walkers, n.dim, init.range) {
  
  
  # initial values
  
  
  chain.length = max.iter%/%n.walkers
  sum.log.p.s.m <<- rep(0, chain.length)
  
  log.p = matrix(NA,nrow=n.walkers,ncol=chain.length)
  log.p.old = rep(NA,n.walkers)
  ensemble.old = matrix(NA, nrow=n.walkers, ncol=n.dim)
  ensemble.new = matrix(NA, nrow=n.walkers, ncol=n.dim)
  x.chain = array(NA, dim=c(n.walkers,chain.length,n.dim))
  
  
  for(k in 1:n.walkers) {
    for(g in 1:n.dim){
      ensemble.old[k,g] = runif(1, init.range[1], init.range[2])
    }
    log.p.old[k] = f(ensemble.old[k,])
  }
  
  log.p[,1]=log.p.old
  sum.log.p.s.m[1] <<- sum(log.p.old[1:n.walkers])/n.walkers
  
  x.chain[ , 1, ] = ensemble.old
  
  
  # the loop
  
  for (l in 2:chain.length){
    for (n in 1:n.walkers){
      
      
      z = ((runif(1)+1)^2)/2
      a = sample((1:n.walkers)[-n], 1)
      par.active = ensemble.old[a,]  
      
      
      ensemble.new[n,] = par.active + z^(n.dim-1)*(ensemble.old[n,] - par.active)
      
      
      log.p.new = f(ensemble.new[n,])
      acc = z*exp(log.p.new - log.p.old[n])
      test = runif(1)
      
      if (acc > test) { 
        
        x.chain[n,l,] = ensemble.new[n,]
        ensemble.old[n,] = ensemble.new[n,]
        log.p[n,l] = log.p.new
        log.p.old[n] = log.p.new
        
        
      } else {
        
        x.chain[n,l,] = ensemble.old[n,]
        log.p[n,l] = log.p.old[n]
        
      }
      
      
      sum.log.p.s.m[l] <<- sum.log.p.s.m[l]+log.p[n,l]/n.walkers
    }
  }
  
  return(x.chain)
  
}
