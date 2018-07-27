# Purpose: Analyse the pilot study
# Functions: sim.GenStudy
# Inputs: data/sims/pilot.sim.specs.rData
#         data/sims/pilot.sim.speclists.rData
# Outputs: data/raw/sim1.res.rData

date()
rm(list = ls())

##### Globals and libraries 
library(dplyr)
library(survival)
library(parallel)
library(pbapply)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

######### load the specs

load(file = "data/sims/pilot.sim.speclists.rData")
load(file = "data/sims/pilot.sim.specs.rData")

# gen sim1 - 

fsim1 <- function(parms){
    df <- do.call(sim.GenStudy, parms)
    cx <- coxph(Surv(os.t,os.e)~x.trt, data = df)
    p.val <- 1-pchisq(cx$score,df=1)
    est <- exp(coef(cx)[1])
    ci  <- exp(confint(cx))
    a.switch <- sum(df$x.switch)/250
    rc  <- c(df$simid[1],est,ci,p.val,a.switch)
    return(rc)
}

n.batch <- length(sim.pilot.spec.list)/10
for (batch in 1:10){
  rc.list <- pblapply(sim.pilot.spec.list[((batch-1)*n.batch+1):(batch*n.batch)],fsim1)
  rc.df   <- as.data.frame(matrix(unlist(rc.list),nrow = length(rc.list), byrow=TRUE))
  names(rc.df) <- c("simid", "hr","cil","ciu","pval","a.switch")
  assign(paste("rc",batch,sep="_"), rc.df)
}

rc <- rbind(rc_1, rc_2, rc_3, rc_4, rc_5, rc_6, rc_7, rc_8, rc_9, rc_10)

############# combine with the sim info

sim1.res.df <- left_join(filter(spec.df, sim==1),rc, by ="simid")

############# save

save(sim1.res.df, file = "data/raw/sim1.res.rData")


