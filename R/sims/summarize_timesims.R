## summarize time results ###
library(tidyverse)
source("R/sims/simutils.R")
source("R/sims/gendata.R")
source("R/utils/estimators.R")
source("R/utils/breslow_utilsv2.R")
#taus <- c(0.5, 1, 2)
F0  <- 0.1
F1  <- 0.05
alpha <- 0.5
taus <- seq(0, 2, by = 1/52) # evaluate risk each week
tau <- 1

timeres <- read.csv("output/timesres1all2021-10-27.csv")

# get truth

dat <- dgm2c(213, 1e6, F1 = F1, F0 = F0, tau=tau, alpha = alpha, cbeta = 1) # large sample to get truth
true1 <- vector()
true0 <- vector()
for(i in 1:length(taus)){
  true1[i] <- mean(dat$ts1<=taus[i])
  true0[i] <- mean(dat$ts0<=taus[i])
}
true2 <-  true1-true0


timeres_sum <- summarize_sims(timeres, true2, out = "timeres")

timeres_sum$methlab <-  rep(c("Crude", "Discrete time (12)", "Discrete time (6)", "Semiparametric" ), # "Discrete hazards (12, 3k)","Discrete hazards (12, 5k)", "Discrete hazards (6, 3k)", "Discrete hazards (6, 5k)","Semiparametric"
                        length(taus))
# make plot
library(ggthemr)
ggthemr("fresh") #, layout = "clean"
pl <- ggplot()+
  geom_line(aes(x = tau, y = bias*100,  color = factor(methlab), group = factor(methlab)), data= timeres_sum)+
  geom_point(aes(x = tau, y = bias*100, color = factor(methlab), group = factor(methlab), shape = factor(methlab)), size = .75, data= timeres_sum)+
  labs(x = "Time (years)", y = "Bias", color = "Method") +
  scale_color_viridis_d(labels = c("Crude", "Discrete time (12)", "Discrete time (6)", "Semiparametric"), name = "key")+
  scale_shape_discrete(labels = c("Crude", "Discrete time (12)", "Discrete time (6)", "Semiparametric"), name = "key")+
  #scale_color_manual(values = swatch()[2:7])+
  theme(legend.position = c(0.25, 0.85), legend.title = element_blank(),
        legend.key.size = unit(.85, 'lines'),
        legend.background = element_rect(fill="transparent", color = "transparent"), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA, color = NA))
pl

ggsave(pl, file = "text/figs/fig1_time.png", width = 10, height = 7, units = "cm")

