# Purpose: Analyse the main study using following methods
#          On treatment RPSFT
#          - log rank test
#          - wilcoxon test
#          Treatment Group RPSFT
#          - log rank test
#          - wilcoxon test
# Functions: sim.GenStudy
#            RPSFT.2pass
#            RPSFT
#            RPSFT.cox
#            RPSFT.latent
#            RPSFT.trypsi
# Inputs: data/sims/sim.specs.rData
#         data/sims/sim.speclists.rData
# Outputs: data/raw/sim_sim2_nn_resrpsft.rData (nn = 1,2,..,10)
#          data/raw/sim_sim3_nn_resrpsft.rData (nn = 1,2,..,10)
#          data/raw/sim_sim4_nn_resrpsft.rData (nn = 1,2,..,10)
# Note: run on AWS then results downloaded

date()
rm(list = ls())

##### Globals and libraries 
library(dplyr)
library(survival)
library(parallel)
#library(pbapply)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

cores <- 40

######### load the specs

load(file = "data/sims/sim.speclists.rData")
load(file = "data/sims/sim.specs.rData")

#do the rpsft

fsim <- function(parms){
    df <- do.call(sim.GenStudy, parms)
    
    # TREATMENT GROUP RPSFT
    rpsft.input.TG <- transmute(df, 
                                event.time = os.t,
                                trt.ind = x.trt, 
                                censor.ind = os.e, 
                                t.start = ifelse(x.switch|x.trt, pmax(0,t.switch,na.rm = TRUE), 0),
                                t.stop = ifelse(x.switch|x.trt, os.t                          , 0),
                                t.on  = t.stop - t.start, 
                                t.off = os.t - (t.stop - t.start),
                                cutofftime = t.censor
    )
    
    # log rank
    rpsft.TG.LR        <- RPSFT.2pass(rpsft.input.TG, Grho=0)
    # wilcoxon
    rpsft.TG.W     <- RPSFT.2pass(rpsft.input.TG, Grho=1)

    # ON TREATMENT RPSFT
    rpsft.input.OT <- transmute(df, 
                                event.time = os.t,
                                trt.ind = x.trt, 
                                censor.ind = os.e, 
                                t.on  = ifelse(x.trt, pfs.t, ifelse(x.switch, pmin(pfs.t, t.censor-pfs.t), 0)),
                                t.off = os.t - t.on,
                                cutofftime = t.censor
    )
    
    # log rank
    rpsft.OT.LR        <- RPSFT.2pass(rpsft.input.OT, Grho=0)
    
    # wilcoxon
    rpsft.OT.W     <- RPSFT.2pass(rpsft.input.OT, Grho=1)
    
    tg.lr.conv <- rpsft.TG.LR$psi.unique
    tg.lr.rng  <- max(rpsft.TG.LR$psi.found) - min(rpsft.TG.LR$psi.found) 
    tg.lr.psi  <- rpsft.TG.LR$psi.chosen

    tg.w.conv <- rpsft.TG.W$psi.unique
    tg.w.rng  <- max(rpsft.TG.W$psi.found) - min(rpsft.TG.W$psi.found) 
    tg.w.psi  <- rpsft.TG.W$psi.chosen
    
    ot.lr.conv <- rpsft.OT.LR$psi.unique
    ot.lr.rng  <- max(rpsft.OT.LR$psi.found) - min(rpsft.OT.LR$psi.found) 
    ot.lr.psi  <- rpsft.OT.LR$psi.chosen
    
    ot.w.conv <- rpsft.OT.W$psi.unique
    ot.w.rng  <- max(rpsft.OT.W$psi.found) - min(rpsft.OT.W$psi.found) 
    ot.w.psi  <- rpsft.OT.W$psi.chosen
  
    # get cox models if it converged
    if (-Inf < tg.lr.psi & tg.lr.psi < Inf & !is.na(tg.lr.psi)){
      cox.rpsft.TG.LR.TU <- RPSFT.cox(rpsft.input = rpsft.input.TG, 
                                      rpsft.output = rpsft.TG.LR, 
                                      Grho = 0,
                                      use.latent.only = FALSE)
      
      tg.lr.hr    <- exp(coef(cox.rpsft.TG.LR.TU$cox.rpsft))
      tg.lr.hr.cil <- exp(confint(cox.rpsft.TG.LR.TU$cox.rpsft))[1]
      tg.lr.hr.ciu <- exp(confint(cox.rpsft.TG.LR.TU$cox.rpsft))[2]
    } else {
      tg.lr.hr <- NA
      tg.lr.hr.cil <- NA
      tg.lr.hr.ciu <- NA
    }
    
    if (-Inf < ot.lr.psi & ot.lr.psi < Inf & !is.na(ot.lr.psi)){
      
      cox.rpsft.OT.LR.V  <- RPSFT.cox(rpsft.input = rpsft.input.OT, 
                                      rpsft.output = rpsft.OT.LR, 
                                      Grho = 0,
                                      use.latent.only = TRUE)
      
      cox.rpsft.OT.LR.TU <- RPSFT.cox(rpsft.input  = rpsft.input.OT, 
                                      rpsft.output = rpsft.OT.LR, 
                                      Grho = 0,
                                      use.latent.only = FALSE)
      
      ot.lr.v.hr     <- exp(coef(cox.rpsft.OT.LR.V$cox.rpsft))
      ot.lr.v.hr.cil <- exp(confint(cox.rpsft.OT.LR.V$cox.rpsft))[1]
      ot.lr.v.hr.ciu <- exp(confint(cox.rpsft.OT.LR.V$cox.rpsft))[2]
      
      ot.lr.tu.hr     <- exp(coef(cox.rpsft.OT.LR.TU$cox.rpsft))
      ot.lr.tu.hr.cil <- exp(confint(cox.rpsft.OT.LR.TU$cox.rpsft))[1]
      ot.lr.tu.hr.ciu <- exp(confint(cox.rpsft.OT.LR.TU$cox.rpsft))[2]
      
    } else {
      ot.lr.v.hr <- NA
      ot.lr.v.hr.cil <- NA
      ot.lr.v.hr.ciu <- NA
      ot.lr.tu.hr <- NA
      ot.lr.tu.hr.cil <- NA
      ot.lr.tu.hr.ciu <- NA
    }
    
    
    # get cox models if it converged
    if (-Inf < tg.w.psi & tg.w.psi < Inf & !is.na(tg.w.psi)){
      cox.rpsft.TG.W.TU <- RPSFT.cox(rpsft.input = rpsft.input.TG, 
                                     rpsft.output = rpsft.TG.W, 
                                     Grho = 1,
                                     use.latent.only = FALSE)
      
      tg.w.hr    <- exp(coef(cox.rpsft.TG.W.TU$cox.rpsft))
      tg.w.hr.cil <- exp(confint(cox.rpsft.TG.W.TU$cox.rpsft))[1]
      tg.w.hr.ciu <- exp(confint(cox.rpsft.TG.W.TU$cox.rpsft))[2]
    } else {
      tg.w.hr <- NA
      tg.w.hr.cil <- NA
      tg.w.hr.ciu <- NA
    }
    
    if (-Inf < ot.w.psi & ot.w.psi < Inf & !is.na(ot.w.psi)){
      
      cox.rpsft.OT.W.V  <- RPSFT.cox(rpsft.input = rpsft.input.OT, 
                                      rpsft.output = rpsft.OT.W, 
                                      Grho = 1,
                                      use.latent.only = TRUE)
      
      cox.rpsft.OT.W.TU <- RPSFT.cox(rpsft.input  = rpsft.input.OT, 
                                      rpsft.output = rpsft.OT.W, 
                                      Grho = 1,
                                      use.latent.only = FALSE)
      
      ot.w.v.hr     <- exp(coef(cox.rpsft.OT.W.V$cox.rpsft))
      ot.w.v.hr.cil <- exp(confint(cox.rpsft.OT.W.V$cox.rpsft))[1]
      ot.w.v.hr.ciu <- exp(confint(cox.rpsft.OT.W.V$cox.rpsft))[2]
      
      ot.w.tu.hr     <- exp(coef(cox.rpsft.OT.W.TU$cox.rpsft))
      ot.w.tu.hr.cil <- exp(confint(cox.rpsft.OT.W.TU$cox.rpsft))[1]
      ot.w.tu.hr.ciu <- exp(confint(cox.rpsft.OT.W.TU$cox.rpsft))[2]
      
    } else {
      ot.w.v.hr <- NA
      ot.w.v.hr.cil <- NA
      ot.w.v.hr.ciu <- NA
      ot.w.tu.hr <- NA
      ot.w.tu.hr.cil <- NA
      ot.w.tu.hr.ciu <- NA
    }
    
    rc <- c(simid = df$simid[1],
            tg.lr.conv = tg.lr.conv,
            tg.lr.rng = tg.lr.rng,
            tg.lr.psi = tg.lr.psi,
            
            tg.lr.hr = tg.lr.hr,
            tg.lr.hr.cil = tg.lr.hr.cil,
            tg.lr.hr.ciu = tg.lr.hr.ciu,
            
            tg.w.conv = tg.w.conv,
            tg.w.rng = tg.w.rng,
            tg.w.psi = tg.w.psi,
            
            tg.w.hr = tg.w.hr,
            tg.w.hr.cil = tg.w.hr.cil,
            tg.w.hr.ciu = tg.w.hr.ciu,
            
            ot.lr.conv = ot.lr.conv,
            ot.lr.rng = ot.lr.rng,
            ot.lr.psi = ot.lr.psi,
            
            ot.lr.tu.hr = ot.lr.tu.hr,
            ot.lr.tu.hr.cil = ot.lr.tu.hr.cil,
            ot.lr.tu.hr.ciu = ot.lr.tu.hr.ciu,
            
            ot.lr.v.hr = ot.lr.v.hr,
            ot.lr.v.hr.cil = ot.lr.v.hr.cil,
            ot.lr.v.hr.ciu = ot.lr.v.hr.ciu,
            
            ot.w.conv = ot.w.conv,
            ot.w.rng = ot.w.rng,
            ot.w.psi = ot.w.psi,
            
            ot.w.tu.hr = ot.w.tu.hr,
            ot.w.tu.hr.cil = ot.w.tu.hr.cil,
            ot.w.tu.hr.ciu = ot.w.tu.hr.ciu,
            
            ot.w.v.hr = ot.w.v.hr,
            ot.w.v.hr.cil = ot.w.v.hr.cil,
            ot.w.v.hr.ciu = ot.w.v.hr.ciu)
            
    return(rc)
}

fsim.lst <- function(lst,pfx){ 
  n.batch <- length(lst)/10
  for (batch in 1:10){
    rc.list <- mclapply(lst[((batch-1)*n.batch+1):(batch*n.batch)],fsim, mc.cores = cores)
    rc.df   <- as.data.frame(matrix(unlist(rc.list),nrow = length(rc.list), byrow=TRUE))
    names(rc.df) <- names(rc.list[[1]])
    assign(paste(pfx,batch,sep="_"), rc.df)
    save(list=paste(pfx,batch,sep="_"), file = paste("data/raw/sim",pfx,batch,"resrpsft.rData", sep = "_"))
  }
}
  
fsim.lst(sim2.spec.list,pfx = "sim2")
fsim.lst(sim3.spec.list,pfx = "sim3")
fsim.lst(sim4.spec.list,pfx = "sim4")
