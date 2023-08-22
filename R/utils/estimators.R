### estimators for breslow paper ####

require(splines)
require(tidyverse)
require(survival)
require(resample)


#' a function to get the risk difference from a set of curves for time maxtau
#' 
getrd <- function(riskfxn, maxtau){
  rd <- riskfxn %>% 
    group_by(a,t) %>%
    summarize(risk = mean(risk)) %>% 
    pivot_wider(names_from = a,  names_prefix = "risk_", values_from = risk) %>% #names_prefix = "risk_",
    add_row(t = 0, risk_1 = as.numeric(0), risk_0 = as.numeric(0)) %>% 
    arrange(t) %>% 
    fill(risk_1, risk_0) %>% 
    mutate(rd = risk_1 - risk_0) 
  
  rd_at_tau <- vector()
  for(r in 1:length(maxtau)){
    rd_at_tau[r]  <- rd %>% 
      slice(which.min(abs(maxtau[r] - replace(t, t>maxtau[r], Inf)))) %>% 
        pull(rd)
  }

  return(rd_at_tau)
}

#' a function to get the standard error of the risk difference from a set of curves for time tau
#' 
get_delta_se <- function(riskfxn, maxtau){
  rdse <- riskfxn %>% 
    # group_by(a) %>%
    pivot_wider(id_cols = c(t), names_from = a,  values_from = c(std.err, risk)) %>% 
    add_row(t = 0, risk_1 = as.numeric(0), risk_0 = as.numeric(0), std.err_1 = 0, std.err_0 = 0) %>% 
    arrange(t) %>% 
    fill(risk_1, risk_0, std.err_1, std.err_0) %>% 
    mutate(rd = risk_1 - risk_0, rdse = sqrt(std.err_1^2 + std.err_0^2))
  
  rdse_at_tau <- vector()
  for(r in 1:length(maxtau)){
    rdse_at_tau[r]  <- rdse %>% 
      slice(which.min(abs(maxtau[r]-t))) %>% 
      pull(rdse)
  }
  
return(rdse_at_tau)
}



#' Tue KM (simulations only)
#' @param data a dataset with variables t (time), delta (event indicator), and a (tx indicator)
#' @return a list with RD (element 1) and risk functions (element 2)


gettruth <- function(data){
  data$deltat <- 1
  km1 <- survfit(Surv(ts1, deltat) ~ 1, data = data)
  km_true1 <- data.frame(summary(km1)[c("surv", "time", "std.err")]) %>% mutate(a = 1)
  km0 <- survfit(Surv(ts0, deltat) ~ 1, data = data)
  km_true0 <- data.frame(summary(km0)[c("surv", "time", "std.err")]) %>% mutate(a = 0)
  km_true <- rbind(km_true1, km_true0) %>% mutate(risk = 1 - surv, method = 1, t = time) %>% select(t, risk, a, method)
  
  return( tibble(km_true) %>% arrange(a, t))

}



##' a function to fit crude estimator
##' 
##' 
crude <- function(data){#}, full = "no"){
  #fit KM to get risks
  kmc <- survfit(Surv(t, delta) ~ a, data = data)
  km_crude <- data.frame(summary(kmc)[c("strata", "surv", "time", "std.err")]) %>% 
    mutate(a = as.numeric(substr(strata, 3, 3)), risk = 1 - surv, method = 2, t = time)%>% 
    select(a, t, risk, std.err)
  
  return(km_crude)

}

##' a function to fit the semiparmetric g comp estimator

getspgcomp <- function(data, covs){#}, full = "no"){
  # under a = 1
  sp1 <- tibble(breslow_intervention(data, "t", "delta", covs, intvar = a, intval = 1, strat_cox = F)) %>% 
                  mutate(a = 1, risk = 1 - s) %>% 
                  select(a, risk, t) #%>% 
                  #rename(t = times))
  # under a = 0
  sp0 <- tibble(breslow_intervention(data, "t", "delta", covs, intvar = a, intval = 0, strat_cox = F)) %>% 
                  mutate(a = 0, risk = 1 - s) %>% 
                  select(a, risk, t) #%>% 
                  #rename(t = times))
  
  return(tibble(rbind(sp0, sp1)))
}
  


##' a function to fit pooled logit gcomp a la hernan
##' @param multiplier quantity to multiply by t to obtain intervals (if t in years, could be 12 or 52. if t in days, could be 1/7 or 1/30 or 1/365)
getplgcomp <- function(dh_dat, covs, multiplier, maxtau, knots = 3){
  
  dh_dat <- dh_dat %>% mutate(id = row_number())

  # expand dataset
  full_dat <- dh_dat %>% 
    mutate(intervals = ceiling(maxtau*multiplier))%>% 
    uncount(intervals, .id = "k")%>%
    mutate( r = ifelse(k <= ceiling(t*multiplier), 1, 0), 
            deltaind = ifelse(k>=ceiling(t*multiplier), delta, 0),
            k = k/multiplier )


  # fit outcome model
  fit.mod <- glm(as.formula( paste("deltaind ~ a + k + ns(k,", knots, ") +", paste(covs, collapse = "+"))), 
                 data = full_dat %>% filter(r == 1), family = "binomial"(link = "logit"))


  # compute predicted probability after setting noart = 1
  full_dat$yhat_1 <- predict(fit.mod, newdata = full_dat %>% mutate(a = 1), type = "response")
  full_dat$yhat_0 <- predict(fit.mod, newdata = full_dat %>% mutate(a = 0), type = "response")
  
  risks <- full_dat %>% 
    group_by(id) %>% 
    mutate(risk1_i = 1 - cumprod(1 - yhat_1), 
           risk0_i = 1 - cumprod(1 - yhat_0)) %>% 
    group_by(k) %>% 
    summarize(risk1 = mean(risk1_i), risk0 = mean(risk0_i)) %>% 
    mutate(t = k) %>% 
    select(t, risk1, risk0) 
  
  risks <- risks %>% 
    pivot_longer(cols = -t, names_to = "a", names_prefix = "risk", values_to = "risk") %>% 
    mutate(a = as.numeric(a))
  
  return(risks)

}


##' a function to bootstrap an arbitrary function to get the RD
##' @param B number of bootstrap samples
##' @param data
##' @param taus vector of timepoints of interest
##' @param function

bootstrap <- function(B, bootdata, taus, rdfunc, ...){
  if(B == 0) return(0)
  if(B>0){
    bsrd <- matrix(NaN, nrow = B, ncol = length(taus))
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      riski <- rdfunc(dati, ...)
      for(k in 1:length(taus)){
        bsrd[i,k] <- getrd(riski, taus[k])
      }
    }
  }
  return(bsrd)
}




##' a function to bootstrap an arbitrary function to get the RD
##' @param B number of bootstrap samples
##' @param data
##' @param taus vector of timepoints of interest
##' @param function

bootstrap_risk <- function(B, data, teval, rdfunc, ...){
  if(B == 0) return(0)
  if(B>0){
    risks <- matrix(0, nrow = B, ncol = 2)
    datbi<-samp.bootstrap(nrow(data), B)
    for(i in 1:B){
      dati <- data[datbi[,i],]
      riski <- rdfunc(dati, ...)
      risks[i,] <- riski %>% 
        group_by(a) %>% 
        filter(t<=teval) %>% 
        arrange(a, t) %>% 
        filter(t == last(t)) %>% 
        pull(risk)
    }
  }
  return(risks)
}

