# Purpose: Analyse the main study using simple methods
#          ITT
#          Per Protocol
#          Censoring at Switch
#          Treatment as a time varying covariate
# Functions: sim.GenStudy
# Inputs: data/sims/sim.specs.rData
#         data/sims/sim.speclists.rData
# Outputs: data/raw/sim.simple.res.rData


date()
rm(list = ls())

##### Globals and libraries 
library(dplyr)
library(survival)
library(parallel)

cores <- 40

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

######### load the specs

load(file = "data/sims/sim.speclists.rData")
load(file = "data/sims/sim.specs.rData")

#derive ITT, TVC etc

fsim.simple <- function(parms){
    df <- do.call(sim.GenStudy, parms)
    #itt
    cox.itt   <- coxph(Surv(os.t, os.e) ~ x.trt, data = df)
    
    # per-protocol 
    cox.pp    <- coxph(Surv(os.t, os.e) ~ x.trt, data = filter(df, x.switch!=1))
    
    # Censoring at Xover
    df <- mutate(df, cx.os.t = ifelse(x.switch, t.switch, os.t), cx.os.e = ifelse(x.switch, 0, os.e))
    cox.cx   <- coxph(Surv(cx.os.t, cx.os.e) ~ x.trt, data = df)
    
    # time varying covariate
    df.tvc <- rbind(transmute(filter(df,x.switch==0), pid = patid, t.start = 0,        t.stop = os.t,     tvc.e = os.e, ontrt1 = x.trt, ontrt2=0, ontrt=x.trt),
                    transmute(filter(df,x.switch==1), pid = patid, t.start = 0,        t.stop = t.switch, tvc.e = 0,    ontrt1 = 0,     ontrt2=0, ontrt=0),
                    transmute(filter(df,x.switch==1), pid = patid, t.start = t.switch, t.stop = os.t,     tvc.e = os.e, ontrt1 = 0,     ontrt2=1, ontrt=1))
    
    cox.tvc <- coxph(Surv(t.start, t.stop, tvc.e) ~ ontrt, data = df.tvc)
    cox.tvc2 <- coxph(Surv(t.start, t.stop, tvc.e) ~ ontrt1 + ontrt2, data = df.tvc)
    
    # get est and confint
    frep <- function(cx,sfx){
      est <- exp(coef(cx)[1])
      ci  <- exp(confint(cx))[1,]
      rc <- c(est,ci)
      names(rc) <- paste(c("hr","cil","ciu"),sfx,sep =".")
      return(rc)
    }
    
    rc <- c(simid=df$simid[1],
            frep(cox.itt,"itt"), 
            frep(cox.pp,"pp"),
            frep(cox.cx,"cens"),
            frep(cox.tvc,"tvc"),
            frep(cox.tvc2,"tvc2")
            )
    
    return(rc)
}

fsim.simple.lst <- function(lst){ 
  rc.list <- mclapply(lst,fsim.simple, mc.cores = cores)
  rc.df   <- as.data.frame(matrix(unlist(rc.list),nrow = length(rc.list), byrow=TRUE))
  names(rc.df) <- names(rc.list[[1]])
  return(rc.df)
}

sim2.simple.res.df <-fsim.simple.lst(sim2.spec.list)
sim3.simple.res.df <-fsim.simple.lst(sim3.spec.list)
sim4.simple.res.df <-fsim.simple.lst(sim4.spec.list)

sim.simple.res.df <- rbind(sim2.simple.res.df, sim3.simple.res.df, sim4.simple.res.df)

############# save

sim.simple.res.df <- merge(spec.df, sim.simple.res.df, by ="simid")

save(sim.simple.res.df, file = "data/raw/sim.simple.res.rData")



