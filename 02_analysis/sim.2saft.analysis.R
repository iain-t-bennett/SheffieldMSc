# Purpose: Analyse the main study using following methods
#          Two-stage AFT without recensoring
# Functions: sim.GenStudy
# Inputs: data/sims/sim.speclists.rData
# Outputs: data/raw/sim.2saft.res.rData
# Note: run on AWS then results downloaded

date()
rm(list = ls())

##### Globals and libraries 
library(dplyr)
library(survival)
library(parallel)


##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

######### load the specs

load(file = "data/sims/sim.speclists.rData")

#derive 2saft

parms <- sim2.spec.list[[1]]

fsim.2saft <- function(parms){
    df <- do.call(sim.GenStudy, parms)
  
    #2stage aft
    
    prog.df <- filter(df,x.trt ==0 & pfs.t < os.t)
    df.2saft <- mutate(prog.df, pps.t = os.t - pfs.t)
    
    # get adjusted estimates of survival
    null.mod <- survreg(Surv(pps.t,os.e)~x.switch, data = df.2saft)
    corr.mod <- survreg(Surv(pps.t,os.e)~x.switch + pfs.t, data = df.2saft)
    
    # derive corrected times using the adjusted estimates
    df.2saft.res <- mutate(df,
                           pps.t = os.t - pfs.t,
                           null.coef = coef(null.mod)[2], 
                           corr.coef = coef(corr.mod)[2])
    
    fdf <- mutate(df.2saft.res,
                           corr.os.nmod = ifelse(x.switch==1 & pps.t >0, pfs.t + pps.t * exp(-null.coef),os.t),
                           corr.os.cmod = ifelse(x.switch==1 & pps.t >0, pfs.t + pps.t * exp(-corr.coef),os.t)
                  
                  #         t.cens.nmod  = t.censor * exp(-null.coef),
                  #        cens.os.nmod   = pmin(corr.os.nmod, t.cens.nmod),
                  #        cens.os.nmod.e = ifelse(cens.os.nmod==t.cens.nmod, 0, os.e),
                  #        t.cens.cmod  = t.censor * exp(-corr.coef),
                  #         cens.os.cmod   = pmin(corr.os.cmod, t.cens.cmod),
                  #         cens.os.cmod.e = ifelse(cens.os.cmod==t.cens.cmod, 0, os.e)
    )
    
    x.trt <- fdf$x.trt
    
    # get est and confint
    frep <- function(t,e,sfx){
      
      if (sum(e[x.trt==0])>0 & sum(e[x.trt==1])>0){
        cx <-  coxph(Surv(t, e) ~ x.trt)
        est <- exp(coef(cx)[1])
        ci  <- exp(confint(cx))[1,]
        nev <- cx$nevent
        rc <- c(est,ci,nev,1)
      } else{
        rc <- c(NA,NA,NA,NA,0)
      }
      names(rc) <- paste(c("hr","cil","ciu","nev","conv"),sfx,sep =".")
      return(rc)
    }
    
    rc <- c(simid=df$simid[1],
            frep(t = fdf$corr.os.nmod, e = fdf$os.e,           sfx = "2saft.null"),
            #frep(t = fdf$cens.os.nmod, e = fdf$cens.os.nmod.e, sfx = "2saft.null.cens"),
            frep(t = fdf$corr.os.cmod, e = fdf$os.e,           sfx = "2saft.corr")
            #frep(t = fdf$cens.os.cmod, e = fdf$cens.os.cmod.e, sfx = "2saft.corr.cens")
            )
    
    return(rc)
}

fsim.2saft.lst <- function(lst){ 
  rc.list <- mclapply(lst,fsim.2saft, mc.cores=2)
  rc.df   <- as.data.frame(matrix(unlist(rc.list),nrow = length(rc.list), byrow=TRUE))
  names(rc.df) <- names(rc.list[[1]])
  return(rc.df)
}

sim2.2saft.res.df <-fsim.2saft.lst(sim2.spec.list)
sim3.2saft.res.df <-fsim.2saft.lst(sim3.spec.list)
sim4.2saft.res.df <-fsim.2saft.lst(sim4.spec.list)

sim.2saft.res.df <- rbind(sim2.2saft.res.df, sim3.2saft.res.df, sim4.2saft.res.df)

############# save


save(sim.2saft.res.df, file = "data/raw/sim.2saft.res.rData")




