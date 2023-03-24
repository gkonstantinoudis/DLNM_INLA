
# Created 26.09.2022

library(INLA)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

# DLNM with INLA

setwd("E:/Postdoc Imperial/Projects/Temperature and respi/Italy_tutorial/data/")

data <- readRDS("dataItaly")
head(data)


# Create indexes
data$year <- as.numeric(year(data$date))
data$id.space = as.numeric(as.factor(data$SIGLA))
data$id.year <- data$year - 2010
data$id.week <- week(data$date)


data %>%
  rename(lag0 = mean.temp) %>% 
  group_by(id.year, NOME_PROVINCIA) %>% 
  arrange(date) %>% 
  mutate(lag1 = lag(lag0),
         lag2 = lag(lag1),
         lag3 = lag(lag2), 
         temperature = (lag0 + lag1 + lag2 + lag3)/3) -> data

data <- data[!is.na(data$temperature),]
data$id.tmp <- inla.group(data$temperature, n = 50, method = "cut", idx.only = TRUE)

# Run INLA to obtain smart starting points to be used later
formula = 
  deaths ~ 1 + offset(log(expected)) + hol + 
  # f(id.week, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) +
  f(id.tmp, model='rw2', hyper=hyper.iid2, constr = TRUE, scale.model = TRUE) +
  f(id.year, model='iid', hyper=hyper.iid, constr = TRUE) + 
  f(id.space, model='bym2', graph="W.adj", scale.model = TRUE, constr = TRUE, hyper = hyper.bym)

# INLA SET UP
# priors
hyper.bym <- list(theta1 = list('PCprior', c(1, 0.01)), theta2 = list('PCprior', c(0.5, 0.5)))
hyper.iid <- list(theta = list(prior="pc.prec", param=c(1, 0.01)))
hyper.iid2 <- list(theta = list(prior="pc.prec", param=c(0.1, 0.01)))
# Under Poisson uses default set up
control.family=inla.set.control.family.default()


m1 = inla(formula,
         data=data,
         family="Poisson",
         control.family=control.family,
         verbose = TRUE,
         num.threads = round(parallel::detectCores()*.5),
         control.compute=list(config = TRUE),
         # control.mode(theta = c(5.4191981, 8.2504967, 5.5434304, -0.6974143)),
         control.predictor=list(link = 1))


summary(m1)
m1$mode$theta
plot(m1$summary.random$id.tmp$`0.5quant`, type = "o")






# What if I omit the lag dimension for now and focus on time/year


# data$id.tmp.period <- paste0(data$id.period, data$id.tmp) %>% as.factor() %>% as.numeric()
data$id.tmp.period <- interaction(data$id.tmp, data$id.year) %>% as.numeric()


P <- table(data$id.tmp) %>% length()
K <- table(data$id.year) %>% length()
### unstructured
Q_lag_unstruc <- diag(K)
### structured
Q_tmp_struc <- INLA:::inla.rw(n = P, order = 2, scale.model = TRUE, sparse = TRUE)
Q_tmp_period <- Q_tmp_struc %x% Q_lag_unstruc

eigenQ <- eigen(Q_tmp_period)
ids <- which(eigenQ$values < 1e-10)
cMat <- t(eigenQ$vectors[,ids])


formula = 
  deaths ~ 1 + offset(log(expected)) + hol + 
  # f(id.week, model='rw1', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) +
  f(id.year, model='iid', hyper=hyper.iid, constr = TRUE) + 
  f(id.space, model='bym2', graph="W.adj", scale.model = TRUE, constr = TRUE, hyper = hyper.bym) +
  f(id.tmp.period, model = "generic0", Cmatrix = Q_tmp_period, 
    extraconstr = list(A = cMat, e = rep(0, nrow(cMat))), rankdef = nrow(cMat), hyper = hyper.iid2)


m2 = inla(formula,
          data=data,
          family="Poisson",
          control.family=control.family,
          verbose = TRUE,
          num.threads = round(parallel::detectCores()*.5),
          control.compute=list(config = TRUE),
          control.mode(theta = c(5.4191981, 8.2504967, 5.5434304, -0.6974143)),
          control.predictor=list(link = 1))



summary(m2)
data$fit.m2 = m2$summary.random$id.tmp.period$`0.5quant`[data$id.tmp.period]
data$fit.m1 = m1$summary.random$id.tmp$`0.5quant`[data$id.tmp]


plot(m2$summary.random$id.tmp$`0.5quant`)
points(m2$summary.random$id.tmp$`0.5quant`[49:100], col = "red")




data %>% 
  ggplot() + geom_point(aes(x=temperature, y = fit.m2, col = id.year)) +
  geom_line(aes(x=temperature, y = fit.m2, col = id.year))  + 
  geom_line(aes(x=temperature, y = fit.m1))


data %>% filter(id.year == 8) %>% 
  ggplot() +
  geom_line(aes(x=temperature, y = fit.m2))  + 
  geom_line(aes(x=temperature, y = fit.m1))

##
## A better plot

data[!duplicated(data$id.tmp.period),] %>% View()

plot(m2$summary.random$id.tmp.period$`0.5quant`)

# # Now I will try a simple dlnm-like here
# 
# # first calculate lags
# data %>%
#   rename(lag0 = id.tmp) %>% 
#   group_by(id.year) %>% 
#   arrange(date) %>% 
#   mutate(lag1 = lag(lag0),
#          lag2 = lag(lag1),
#          lag3 = lag(lag2)) -> data_2
# 
# data_2$pm25 <- data_2$no2 <- data_2$o3 <- NULL
# 
# 
# data_long <- gather(data_2, lag, temperature, lag0:lag3, factor_key=TRUE)
# data_long <- data_long[complete.cases(data_long$temperature),]
# data_long$id_lag_tmp <- paste0(data_long$temperature, data_long$lag) %>% as.factor() %>% as.numeric()
# P <-  table(data_long$temperature) %>% length()
# ### unstructured
# Q_lag_unstruc <- diag(4)
# ### structured
# Q_tmp_struc <- INLA:::inla.rw(n = P, order = 2, scale.model = TRUE, sparse = TRUE)
# Q_tmp_lag <- Q_tmp_struc %x% Q_lag_unstruc
# 
# eigenQ <- eigen(Q_tmp_lag)
# ids <- which(eigenQ$values < 1e-10)
# cMat <- t(eigenQ$vectors[,ids])
# 
# 
# formula = 
#   deaths ~ 1 + offset(log(expected)) + hol + 
#   f(id.year, model='iid', hyper=hyper.iid, constr = TRUE) + 
#   f(id.space, model='bym2', graph="W.adj", scale.model = TRUE, constr = TRUE, hyper = hyper.bym) +
#   f(id_lag_tmp, model = "generic0", Cmatrix = Q_tmp_lag, 
#     extraconstr = list(A = cMat, e = rep(0, nrow(cMat))), rankdef = nrow(cMat), hyper = hyper.iid)
# 
# 
# m2 = inla(formula,
#          data=data_long,
#          family="Poisson",
#          control.family=control.family,
#          verbose = TRUE,
#          num.threads = round(parallel::detectCores()*.5),
#          control.compute=list(config = TRUE),
#          control.mode(theta = c(5.4191981, 8.2504967, 5.5434304, -0.6974143)),
#          control.predictor=list(link = 1))
# 
# summary(m2)
# data_2$fit.m2 = m2$summary.random$id_lag_tmp$`0.5quant`[data_2$id_lag_tmp]
# data_2 %>% 
#   mutate(id.lag = as.factor(lag)) %>% 
#   filter(!duplicated(id_lag_tmp)) %>% 
#   ggplot() + geom_point(aes(x=temperature, y = fit.m2, col = id.lag)) 
# 
# 

