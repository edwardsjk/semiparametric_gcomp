## Script to induce confounding in ACTG Data ##


# read in data
base <- "data/actg320.23nov16.dat"
dat_a <- read_fwf(base, fwf_cols(id = c(1,6), male = 8, black = 10, hispanic = 12, idu = 14, 
                                 a = 16, delta = 18, cen = 20, r = 22, age = c(24, 25), 
                                 karnof = c(27,29), days = c(31,33), cd4 = c(35,37), stop = c(39,41) 
)
) %>% 
  mutate(z = cd4/100, t = days/365.25)

#induce confounding
set.seed(5)
tmp1 <- dat_a %>% filter(a == 0 | z<=1)
tmp2 <- dat_a %>% filter(a == 1 & z > 1) %>% sample_frac(0.25)
dat_b <- bind_rows(tmp1, tmp2)


dat_b <- dat_b %>% 
  mutate(agecat = ifelse(age<30, 0, 
                         ifelse(age>=30 & age < 40, 1, 
                                ifelse(age>=40 & age < 50, 2, 
                                       ifelse(age >= 50, 3, NA)))), 
         agecat2 = as.numeric(age>40), 
         lowcd4 = as.numeric(z<1))


write.csv(dat_b, file = "../data/exampledat.csv")