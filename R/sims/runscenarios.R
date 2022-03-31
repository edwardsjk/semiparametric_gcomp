# Code to runs sims for all scenarios #


library(foreach)
library(doParallel)
library(survival)
library(tidyverse)
library(resample)
library(splines)
source("R/utils/breslow_utilsv2.R")
source("R/sims/simutils.R")
source("R/sims/gendata.R")
source("R/utils/estimators.R")

scenario <- 1 #update to change scenario
# DGM CONTROL ----
cen <-  ifelse(scenario %in% c(1, 2), 1, 2)
## CONTROL N -----
N <- ifelse(scenario %in% c(1, 3), 1000, 5000)



# common inputs ----
J   <- 1000 
B   <- 500
F0  <- 0.1
F1  <- 0.05
alpha <- 0.5
taus <- c(0.5, 1, 2)
tau <- 1



start <- Sys.time()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
allres <- foreach(k = 1:J, .combine = "rbind",
.packages = c("survival", "tidyverse", "splines", "resample")) %dopar%
  fn(k, N, B, F1 = F1, F0 = F0, tau = tau, alpha = alpha,
    cbeta = cen, taus = taus,
     estimators = c("crude", "gcomp_sp", "gcomp_dh"))
stopCluster(cl)
end <- Sys.time()

elapsed <- end - start
elapsed

write.csv(allres, 
file = paste0("output/", scenario, "_", Sys.Date(), ".csv"))