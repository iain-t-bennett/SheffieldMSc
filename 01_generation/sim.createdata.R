# Purpose: Create main simulations 
# Functions: sim.GenStudy
# Inputs: none
# Outputs: data/sims/sim.specs.rData
#          data/sims/sim.speclists.rData


date()
rm(list = ls())

##### Globals and libraries 
library(dplyr)
library(survival)
library(parallel)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

######### run it

set.seed(1984)

sim1.grid.df = expand.grid(rep = 1:1000, 
                           rho = c( -.4, -.2 , 0, .2, .4, .6, .8, 1),
                           p.switch = c(0, 0.2, 0.4, 0.6), 
                           scenario = c(1:3), 
                           prop.pd.int = c(0,0.1,0.2,0.3,0.4,0.5),
                           sim = 1
)


sim2.grid.df = expand.grid(rep = 1:1000, 
                      rho = c(-.4, -.2 , 0, .2,  .4, .6, .8, 1),
                      p.switch = c(0.2, 0.4, 0.6), 
                      prop.pd.int = 0,
                      scenario = c(1:3), 
                      sim = 2
)
sim3.grid.df = expand.grid(rep = 1:1000, 
                           rho = c(-.4, -.2 , 0, .2,  .4, .6, .8, 1),
                           p.switch = c(0.4, 0.6), 
                           scenario = c(1:4), 
                           prop.pd.int = 0,
                           sim = 3
)

sim4.grid.df = expand.grid(rep = 1:1000, 
                           rho = c( -.4, -.2 , 0, .2, .4, .6, .8, 1),
                           p.switch = c(0.4, 0.6), 
                           scenario = c(1:3), 
                           prop.pd.int = 0,
                           sim = 4
)

# common treatment effect and sustained effect - vary time when switch is allowed
sim1.spec.te <- rbind(data_frame(scenario = 1, beta_1a = log(0.7)),
                      data_frame(scenario = 2, beta_1a = log(0.9)),
                      data_frame(scenario = 3, beta_1a = log(1)))
sim1.spec.te <- mutate(sim1.spec.te, beta_1b = beta_1a,  beta_2a = beta_1a,  beta_2b = beta_1a)



# common treatment effect and sustained effect
sim2.spec.te <- rbind(data_frame(scenario = 1, beta_1a = log(0.7)),
                      data_frame(scenario = 2, beta_1a = log(0.9)),
                      data_frame(scenario = 3, beta_1a = log(1)))
sim2.spec.te <- mutate(sim2.spec.te, beta_1b = beta_1a,  beta_2a = beta_1a,  beta_2b = beta_1a)

# common treatment effect and no sustained effect
sim3.spec.te <- rbind(data_frame(scenario = 1, beta_1a = log(0.01), beta_1b = log(1)),# circa 0.7
                      data_frame(scenario = 2, beta_1a = log(0.4),  beta_1b = log(1)), # circa 0.9
                      data_frame(scenario = 3, beta_1a = log(0.5),  beta_1b = log(0.8)),# circa 0.7
                      data_frame(scenario = 4, beta_1a = log(0.8),  beta_1b = log(0.95))) # circa 0.9
sim3.spec.te <- mutate(sim3.spec.te, beta_2a = beta_1a, beta_2b = beta_1b)

# reduced/increased effect
sim4.spec.te <- rbind(data_frame(scenario = 1, beta_1a = log(0.7), beta_2a = log(0.8)),
                      data_frame(scenario = 2, beta_1a = log(0.7), beta_2a = log(0.9)),
                      data_frame(scenario = 3, beta_1a = log(0.7), beta_2a = log(1.25)))
sim4.spec.te <- mutate(sim4.spec.te, beta_1b = beta_1a, beta_2b = beta_2a)


sim1.spec <- left_join(sim1.grid.df, sim1.spec.te, by = "scenario")
sim2.spec <- left_join(sim2.grid.df, sim2.spec.te, by = "scenario")
sim3.spec <- left_join(sim3.grid.df, sim3.spec.te, by = "scenario")
sim4.spec <- left_join(sim4.grid.df, sim4.spec.te, by = "scenario")

spec.df <- rbind(sim1.spec, sim2.spec, sim3.spec, sim4.spec) 
spec.df <- cbind(spec.df, seed = round(runif(nrow(spec.df),0,1e9)), simid = 1:nrow(spec.df))
spec.df <- transmute(spec.df, sim, scenario, simid, rho, p.switch, beta_1a, beta_1b, beta_2a, beta_2b, prop.pd.int, seed)

# get 'true values'
u.spec.df <- transmute(ungroup(slice(group_by(spec.df, sim, scenario, rho, beta_1a, beta_1b),1)),sim, scenario, beta_1a, beta_1b, rho)
true.hr <- rep(NA, nrow(u.spec.df))
for (i in 1:nrow(u.spec.df)){
  prms <- as.list(u.spec.df[i,3:5])
  if (prms$beta_1a==prms$beta_1b){
    true.hr[i]<-exp(prms$beta_1a)
  } else{
    prms$seed = 1234
    prms$arm.n = 500000
    prms$simid = i
    prms$beta_2a = log(1)
    prms$beta_2b = log(1)
    prms$p.switch = 0.2
    prms$prop.pd.int = 0
    bigdf <- do.call(sim.GenStudy,prms)
    cx <- coxph(Surv(os.o.t, os.o.e)~x.trt, data = bigdf)
    true.hr[i]<-exp(cx$coefficients)[1]
  }
}
u.spec.df <- transmute(u.spec.df, sim, scenario, rho, true.hr = true.hr)
spec.df <- left_join(spec.df, u.spec.df, by=c("sim","scenario","rho"))

save(spec.df, file = "data/sims/sim.specs.rData")


######### save the lists 
rm(list=ls())

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)


load("data/sims/sim.specs.rData")

convspec2list <- function(df){
  df <- transmute(df, seed, simid,p.switch, prop.pd.int, rho, beta_1a, beta_1b, beta_2a, beta_2b)
  rc <- list()
  for (i in 1:nrow(df)){
    rc[[i]] <- as.list(df[i,])
  }
  return(rc)
}

sim1.spec.list <-convspec2list(filter(spec.df,sim==1))
sim2.spec.list <-convspec2list(filter(spec.df,sim==2))
sim3.spec.list <-convspec2list(filter(spec.df,sim==3))
sim4.spec.list <-convspec2list(filter(spec.df,sim==4))

save(sim1.spec.list, sim2.spec.list, sim3.spec.list,sim4.spec.list, file = "data/sims/sim.speclists.rData")
