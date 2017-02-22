# ......................................................o
# .......... Test functions for testing samplers.........
# ................. Sanda Dejanic ......................o

library(mvtnorm)

# ................o
# 6 Hump Camel Back
# ................o

log.post.camel = function(xx){
 
  x1 = xx[1]
  x2 = xx[2]
  
  term1 = (4-2.1*x1^2+(x1^4)/3) * x1^2
  term2 = x1*x2
  term3 = (-4+4*x2^2) * x2^2
  
  y = term1 + term2 + term3
  return(-y)
}


# ................o
# 3 Hump Camel Back
# ................o

log.post.camel.3 = function(xx){
 
  x1 = xx[1]
  x2 = xx[2]
  
  term1 = -1.05*x1^4 + (x1^6)/6 + 2*x1^2
  term2 = x1*x2
  term3 = x2^2
  
  y = term1 + term2 + term3
  return(-y)
}


# .........o
# Rosenbrock
# .........o

log.post.rosenbrock = function(xx){
  
  x1 = xx[1]
  x2 = xx[2]
  
  y = - (1-x1)^2 - 100*(x2-x1^2)^2
  return(y)
}


# .........o
# Mild Rosenbrock
# .........o

log.post.mild.rosenbrock = function(xx){
  
  y = - (1-xx[1])^2 - 10*(xx[2]-xx[1]^2)^2
  return(y)
}


# .........o
# 20d Rosenbrock
# .........o


log.post.20d.rosenbrock = function(xx){
  
  y = 0
  
  for (dim in 1:19){
  
  w = - (1-xx[dim])^2 - 100*(xx[dim+1]-xx[dim]^2)^2
  y = y + w
  
  }
  
  return(y)
}


# .........o
# 20d Mild Rosenbrock
# .........o

log.post.20d.mild.rosenbrock = function(xx){
  
  y = 0
  
  for (dim in 1:19){
    
    w = - (1-xx[dim])^2 - 10*(xx[dim+1]-xx[dim]^2)^2
    y = y + w
    
  }
  
  return(y)
}



# .............................o
# comet
# ....................................o

comet = function(xx){
  
  x1 = xx[1]
  x2 = xx[2]
  
  first = 1/2*exp(-x1^2/4) 
  second = 1/2*exp(-x2^2/4) + 1/3*exp(-(6-x2)^2/30)         
  y = -(first * second)
  return(y)
  
}

# ....................................o
# Raindrop
# .......................o

raindrop = function(xx){
  
  x1 = xx[1]
  x2 = xx[2]
  
  first = 1/2*exp(-x1^2/4) 
  second = 1/2*exp(-x2^2/4) + 1/2*exp(-(3-x2)^2/25)         
  y = -(first * second)
  return(y)
  
}


# .....
# multi variate normal
# .....

log.post.normal = function(xx){
  y = dmvnorm(xx, mean = rep(0, length(xx)), sigma = diag(length(xx)), log = T)  
  return(y)
}


# .........
# mv bi normal
# ............

log.post.bi.normal = function(xx){
  
  d1 = dmvnorm(xx, sigma = diag(length(xx))) 
  d2 = dmvnorm(xx, rep(4, length(xx)), sigma = diag(length(xx)))
  md = log(0.65*d1 + 0.35*d2)
  return(md)
}


# .........
# mv multi modal normal - not the same as in julia because of sigma
# ............


log.post.funny.normal = function(xx){
  d1 = 0.6*dmvnorm(xx, sigma = 2*diag(length(xx))) 
  d2 = 0.26*dmvnorm(xx, rep(9, length(xx)), sigma = diag(length(xx)))
  d3 = 0.11*dmvnorm(xx, rep(16, length(xx)), sigma = 0.5*diag(length(xx))) 
  d4 = 0.03*dmvnorm(xx, rep(19, length(xx)), sigma = 0.2*diag(length(xx)))
  d = log(d1 + d2 + d3 + d4) 
  return(d)
}

# ...........o
# image plot
# ................o

# xp=matrix(NA, nrow=40, ncol=2)
# x=matrix(NA, nrow=40, ncol=2)
# xp[,1] = seq(0.1,40)
# xp[,2] = seq(0.1,40)
# x[,1] = xp[,1] - 10
# x[,2] = xp[,2] - 10
# 
# for.image = matrix(apply(expand.grid(x[,1], x[,2]), 1,  log.post.R), nrow=40)
# image(x[,1], x[,2], exp(for.image), col=cm.colors(60))


# .................................
# 1d examples for visualisation funny normal
# ..................

multi.modal.normal = function(xx){
  y = 0.6 * dnorm(xx, mean=0, 2, log=F) + 0.26 * dnorm(xx, mean=9, 1, log=F) + 0.11 * dnorm(xx, mean=16, 0.5, log=F) + 0.03 * dnorm(xx, mean=19, 0.2, log=F)
  return(y)
}


multi.modal.normal.log = function(xx){
y = dnorm(xx, mean=0, 2, log=T) + dnorm(xx, mean=9, 1, log=T)/3 + dnorm(xx, mean=16, 0.5, log=T)/9 + dnorm(xx, mean=19, 0.2, log=T)/27
  return(y)
}


bi.modal.normal = function(xx){
  
  y = 0.65 * dnorm(xx) + 0.35 * dnorm(xx, mean=4) 
  w = log(y)
  return(w)
}


# plot
x = seq(-7, 22, 0.1)
y = multi.modal.normal(x)
#plot(x,y)

# ..............................................o
  
 

