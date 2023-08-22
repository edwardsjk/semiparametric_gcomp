### Breslow Utils v2.0 ###


##' a function to estimate betahat from cox model
##' @param data a dataset
##' @param model a model 
##' @return betahat 

get_betahat <- function(data, model){
  coxmod <- coxph(as.formula(model),
                  data = data, method = "breslow")
  betahat <- coxmod$coef
  return(betahat)
}

##' a function to compute baseline hazard function
##' @param event_times a vector of event times
##' @param events number of events at each event time
##' @param covars a matrix of covariates
##' @param betahat output from getbetahat() above
##' @return cumulative baseline hazard function

get_baseline_hazard <- function(event_times, events, covars, betahat, output = "unique"){
  event_times <- unname(unlist(event_times))
  if(output == "unique") unique.times <-(unique(event_times))
  if(output == "all") unique.times <-((event_times))
  x <- as.matrix(covars)
  y <- events #unname(unlist(events))
  longy <- matrix(rep(unlist(y), length(unique.times)), nrow= length(unique.times), byrow = T)
  r <- (outer(event_times, event_times, function(a, b) as.integer(a-b<=0 )))
  v <- (outer(unique.times, event_times, function(a, b) as.integer(a-b>=0)))
  den <-  t(exp(x %*% (betahat)) ) %*% t(r)
  H0 <-   (v * longy) %*% t(1/den) 
  return(H0)
}


##' a function to predict survival
##' @param baseline_hazard cumulative baseline hazard function
##' @param x matrix of covariates
##' @param betahat vector of cox model coefficients
##' @return survival function

get_survival_breslow <- function(baseline_hazard, x, betahat){
  Hi <- baseline_hazard %*% t(exp(x %*% betahat)) #exp(x %*% (betahat)) %*% t(baseline_hazard)
  si <- exp(-Hi)
  s <- rowMeans(si) #rowSums(si)/ncol(si)
  return(s)
}


##' wrapper function to fit breslow estimator
##' @param data a dataset
##' @param time_var character string naming the time variable
##' @param event_var character string denoting the event indicator
##' @param covar matrix of covariates
##' @param intvar intervention variable
##' @param intval value to set intervention variable to
##' @param strat_cox indicator of whether to fit stratified cox models
##' @return a dataframe with each event time and survival under exposure intval

breslow_intervention <- function(data, time_var, event_var, covar, intvar, intval, strat_cox = F){

  intvarc <- deparse( substitute(intvar) )

  # part 1: fix cox model
  if(strat_cox == T){
    mod <- paste0("Surv(", paste(time_var, event_var, sep = ","), ") ~ ", 
                  paste(covar, collapse = "+"), 
                  paste0("+ strata(", intvarc, ")"))
    coxdat <- data
  }
  if(strat_cox == F){
    mod <- paste0("Surv(", paste(time_var, event_var, sep = ","), ") ~ ", 
                  paste(covar, collapse = "+"))
    coxdat <- data %>% filter({{intvar}} == intval)
  }
  betahat <- get_betahat(coxdat, mod)
  
  #part 2: fit bh function
  bhdat <- data %>% filter({{intvar}} == intval) #subset to those with A=a
  bh <- get_baseline_hazard(bhdat[,time_var], bhdat[,event_var], as.matrix(bhdat[, covar]), betahat)
  bhtimes <- unique(bhdat[,time_var])
  
  # part 3: predict S for all in dataset and summarize
  
  allx <- as.matrix(data[, covar])
  s <- get_survival_breslow(bh, allx, betahat)
  res <- data.frame(bhtimes, s) %>% arrange(bhtimes) %>% rename(t = bhtimes)
  return(res)
}

