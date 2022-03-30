### generate data for sims ###

dgm2c <- function(i, n, F1 = .05, F0 = .1, tau=1, alpha = 1/2, cbeta = 1, maxtau = 2){
  
  set.seed(i)
  n <- n 
  id <- c(1:n)
  z <- rnorm(n, 0, 1) #confounder
  pa <- plogis(log(0.4/(1-0.4)) - log(4) * 0 + log(4) * z) # probability oftreatment
  a <- rbinom(n, 1, pa)
  h1 <- exp(-6.22  + log(3)*z)
  h0 <- exp(-4.77 + log(3)*z)
  ts1<-rweibull(n, shape=alpha, scale=1/h1) # time to event, a=1
  ts0<-rweibull(n, shape=alpha, scale=1/h0) # time to event, a=0
  tc <- rexp(n) * exp(log(3) + log(2) * a + log(cbeta) * z) 
  
  # potential outcomes
  delta1 <- 1 
  t1 <- ts1 
  delta0 <- 1 
  t0 <- ts0 
  
  # observed outcomes
  ts <- ifelse(a == 1, ts1, ts0)
  t <- pmin(tc, ts, maxtau) 
  delta <- as.numeric((ts <= pmin(tc, maxtau)))
  
  cen <- 1 - delta
  
  dat_a <- data.frame(id, z, a, t, delta, cen, ts1, ts0, delta1, delta0)
  
  return(dat_a)
}

