### utils to run sims ###

##' a function to run all estimators
fn <- function(i, n, B, F1 = .15, F0 = .3, tau=3, alpha = 1/2, cbeta = 1, taus, 
               estimators = c("crude",  "gcomp_sp", "gcomp_dh", "gcomp_dh6")){

  dat <- dgm2c(i, n, F1 = F1, F0 = F0, tau=tau, alpha = alpha, cbeta = cbeta)
  
  rd_ests <- vector()
  ses <- vector()
  ### 2) Crude ----
  if("crude" %in% estimators){
    crude_ests <- crude(dat) 
    rd_crude <- getrd(crude_ests, taus)
    se_crude <- get_delta_se(crude_ests, taus)
    rd_ests <- c(rd_ests, rd_crude)
    ses <- c(ses, se_crude)
  }
  

 
  ### 3) Semiparametric g-comp ----
  if("gcomp_sp" %in% estimators){
    risks_sp <- getspgcomp_beta(dat, "z")
    rd_sp <- getrd(risks_sp, taus)
    if(B == 0) se_sp <- rep(0, length(taus))
    if(B>0){
      bootrds <- bootstrap(B, dat, taus = taus, getspgcomp_beta, covs = "z")
      se_sp <- apply(bootrds, 2, sd)
    }
    rd_ests <- c(rd_ests, rd_sp)
    ses <- c(ses, se_sp)
  }

  
  
  ### 4) Pooled logit g-comp ----
  if("gcomp_dh" %in% estimators){
    risks_dh <- getplgcomp(dat, "z", 12, maxtau = max(taus)) 
    rd_dh <- getrd(risks_dh, taus)
    if(B == 0) se_dh <- rep(0, length(taus))
    if(B>0){
      bootrds_dh <- bootstrap(B, dat, taus = taus, getplgcomp, covs = "z", 12, maxtau = max(taus))
      se_dh <- apply(bootrds_dh, 2, sd)
    }
    rd_ests <- c(rd_ests, rd_dh)
    ses <- c(ses, se_dh)
  }
  
  ### 5) Pooled logit g-comp ----
  if("gcomp_dh6" %in% estimators){
    risks_dh6 <- getplgcomp(dat, "z", 6, maxtau = max(taus)) 
    rd_dh6 <- getrd(risks_dh6, taus)
    if(B == 0) se_dh6 <- rep(0, length(taus))
    if(B>0){
      bootrds_dh6 <- bootstrap(B, dat, taus = taus, getplgcomp, covs = "z", 6, maxtau = max(taus))
      se_dh6 <- apply(bootrds_dh6, 2, sd)
    }
    rd_ests <- c(rd_ests, rd_dh6)
    ses <- c(ses, se_dh6)
  }
  
  ### 6) Pooled logit g-comp, 5 knots ----
  if("gcomp_dh_5k" %in% estimators){
    risks_dh_5k <- getplgcomp(dat, "z", 12, maxtau = max(taus), knots = 5) 
    rd_dh_5k <- getrd(risks_dh_5k, taus)
    if(B == 0) se_dh_5k <- rep(0, length(taus))
    if(B>0){
      bootrds_dh_5k <- bootstrap(B, dat, taus = taus, getplgcomp, covs = "z", 12, maxtau = max(taus), knots = 5)
      se_dh_5k <- apply(bootrds_dh_5k, 2, sd)
    }
    rd_ests <- c(rd_ests, rd_dh_5k)
    ses <- c(ses, se_dh_5k)
  }
  
  
  ### 6) Pooled logit g-comp, 6 intervals, 5 knots ----
  if("gcomp_dh6_5k" %in% estimators){
    risks_dh6_5k <- getplgcomp(dat, "z", 6, maxtau = max(taus), knots = 5) 
    rd_dh6_5k <- getrd(risks_dh6_5k, taus)
    if(B == 0) se_dh6_5k <- rep(0, length(taus))
    if(B>0){
      bootrds_dh6_5k <- bootstrap(B, dat, taus = taus, getplgcomp, covs = "z", 6, maxtau = max(taus), knots = 5)
      se_dh6_5k <- apply(bootrds_dh6_5k, 2, sd)
    }
    rd_ests <- c(rd_ests, rd_dh6_5k)
    ses <- c(ses, se_dh6_5k)
  }
  
  rds <- data.frame(rd = rd_ests,
                    se = ses,
                    tau = rep(taus, length(estimators)),
                    method =  rep(estimators,  each = length(taus)))
  print(i)
  return(rds)
}

##' a function to summarize simulation results
summarize_sims <- function(allsims, true, outchr = "t1"){
  res1 <- allsims %>% 
    # group_by(tau) %>% 
    mutate(truth = rep(rep(true, length(unique(allsims$method))), nrow(allsims)/(length(taus)*length(unique(allsims$method)))),
           # ifelse(tau == taus[1], true_1, 
           #               ifelse(tau == taus[2], true_2, true_3)), 
           bias = rd - truth,
           cover = (truth <= (rd + 1.96*se) & truth >= (rd - 1.96*se))) %>% 
    group_by(tau, method) %>% 
    summarize(rd_est = mean(rd, na.rm = T), 
              ese = sd(rd, na.rm = T), 
              ase = mean(se, na.rm = T), 
              cover = mean(cover), 
              mc_err = sd(bias)/sqrt(nrow(allsims)/(length(taus)*length(unique(allsims$method)))), #5 = number of estimators
              bias = mean(bias)) %>% 
    mutate(rmse = sqrt(bias^2 + ese^2)) #consider using ase

  
  # make table for paper
  res1_tab <- res1 %>% filter(method!="true") %>% 
    select(tau, method, bias, ese, ase, rmse, cover)
  
  #save all results 
  write.csv(res1, 
            file = paste0("output/res_dgm", outchr, as.character(Sys.Date()),".csv"))
  
  #save again for paper
  write.csv(res1_tab, 
            file =paste0("output/", outchr, ".csv"))
  return(res1_tab)
}



