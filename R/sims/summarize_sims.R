### summarize simulation results ####
library(tidyverse)
source("R/sims/simutils.R")
source("R/sims/gendata.R")
source("R/utils/estimators.R")
taus <- c(0.5, 1, 2)
F0  <- 0.1
F1  <- 0.05
alpha <- 0.5
taus <- c(0.5, 1, 2)
tau <- 1

# read in simulation results
 allres1 <- read.csv(file = "output/Tallres_1_2021-10-27.csv")
 allres3 <- read.csv(file = "output/Tallres_3_2021-10-27.csv")
 allres2 <- read.csv(file = "output/Tallres_2_2021-10-28.csv")
 allres4 <- read.csv(file = "output/Tallres_4_2021-10-28.csv")

#get truth
dat <- dgm2c(213, 1e6, F1 = F1, F0 = F0, tau=tau, alpha = alpha, cbeta = 1, maxtau = 2) # large sample to get truth
true1 <- vector()
true0 <- vector()
for(i in 1:length(taus)){
  true1[i] <- mean(dat$ts1<=taus[i])
  true0[i] <- mean(dat$ts0<=taus[i])
}
true <-  true1-true0

# summarize results
summarize_sims(allres1, true, out = "t1a_1000")
summarize_sims(allres2, true, out = "t1a_2000")
summarize_sims(allres3, true, out = "t1b_1000")
summarize_sims(allres4, true, out = "t1b_2000")
