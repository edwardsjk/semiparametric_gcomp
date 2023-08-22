## Code to replicate results ###

# set up environment ----
setwd("/Users/jkedwar/BoxSync/Research/HIVmeths/BreslowGComp/public")
source("R/utils/breslow_utils.R")
source("R/utils/estimators.R")
B <- 1000
tau <- 1
# read in data
dat_b <- read.csv("data/exampledat.csv")

# analysis ----

## crude ----

crude_gcomp <- crude(dat_b)
crude_rd <- getrd(crude_gcomp, 1)
crude_se <- get_delta_se(crude_gcomp, maxtau = 1)
crude_rd_lcl <- crude_rd - 1.96 * crude_se
crude_rd_ucl <- crude_rd + 1.96 * crude_se
crude_risks <- crude_gcomp %>% 
  group_by(a) %>% 
  filter(t == last(t)) %>% 
  mutate(method = "Crude") %>% 
  select(-std.err)
crude_risks

## discrete gcomp, monthly ----

discrete_gcomp <- getplgcomp(dat_b, c("z", "black", "idu", "age", "male", "hispanic"), 12, maxtau = 1, knots = 4)
discrete_rd <- getrd(discrete_gcomp, 1)
discrete_risk <- discrete_gcomp %>% 
  mutate(method = "Discrete time, monthly")
discrete_risk_tau <- discrete_risk %>% 
  group_by(a) %>% 
  filter(t == last(t))
discrete_risk_tau

# get standard error and 95% CI
boot_rds <- bootstrap(B, dat_b, taus = 1, getplgcomp, covs = c("z", "black", "idu", "age", "male", "hispanic"), multiplier = 12, maxtau = 1, knots = 4)
discrete_se <- sd(boot_rds)
discrete_rd_lcl <- discrete_rd - 1.96 * discrete_se
discrete_rd_ucl <- discrete_rd + 1.96 * discrete_se

# print results
print(paste0(discrete_rd, " (",discrete_rd_lcl, ", ", discrete_rd_ucl, ")"))


## discrete gcomp, weekly ----

discrete_gcomp_wk <- getplgcomp(dat_b, c("z", "black", "idu", "age", "male", "hispanic"), 52, maxtau = 1, knots = 4)
discrete_rd_wk <- getrd(discrete_gcomp_wk, 1)
discrete_risk_wk <- discrete_gcomp_wk %>% 
  mutate(method = "Discrete time, weekly")
discrete_risk_tau_wk <- discrete_risk_wk %>% 
  group_by(a) %>% 
  filter(t == last(t))
discrete_risk_tau_wk

# get standard error and 95% CI
boot_rds <- bootstrap(B, dat_b, taus = 1, getplgcomp, covs = c("z", "black", "idu", "age", "male", "hispanic"), multiplier = 52, maxtau = 1, knots = 4)
discrete_se_wk <- sd(boot_rds)
discrete_rd_lcl_wk <- discrete_rd_wk - 1.96 * discrete_se_wk
discrete_rd_ucl_wk <- discrete_rd_wk + 1.96 * discrete_se_wk

# print results
print(paste0(discrete_rd_wk, " (",discrete_rd_lcl_wk, ", ", discrete_rd_ucl_wk, ")"))

## semiparametric g comp ----

breslow_gcomp <- getspgcomp(dat_b, c("z", "black", "idu", "age", "male", "hispanic"))
breslow_rd <- getrd(breslow_gcomp, 1)
breslow_risk <- breslow_gcomp %>% 
  mutate(method = "Semiparametric")
breslow_risk %>% filter(a==0 & t == last(t))
breslow_risk_tau <- breslow_risk %>% 
  arrange(a, t) %>% 
  group_by(a) %>% 
  filter(t == last(t)) %>% 
  distinct(t, .keep_all = T)
breslow_risk_tau

# get standard error and 95% CI
boot_rds <- bootstrap(B, dat_b, taus = 1, getspgcomp_beta, covs = c("z", "black", "idu", "age", "male", "hispanic"))
breslow_se <- sd(boot_rds)
breslow_rd_lcl <- breslow_rd - 1.96 * breslow_se
breslow_rd_ucl <- breslow_rd + 1.96 * breslow_se

# print results
print(paste0(breslow_rd, " (",breslow_rd_lcl, ", ", breslow_rd_ucl, ")"))