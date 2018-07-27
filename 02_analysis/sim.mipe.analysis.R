# Purpose: Analyse the main study using following methods
#          MIPE
# Functions: sim.GenStudy
#            MIPE.cox
#            IPE_BW
# Inputs: data/sims/sim.specs.rData
#         data/sims/sim.speclists.rData
# Outputs: data/raw/sim_sim2_nn_resmipe.rData (nn = 1,2,..,10)
#          data/raw/sim_sim3_nn_resmipe.rData (nn = 1,2,..,10)
#          data/raw/sim_sim4_nn_resmipe.rData (nn = 1,2,..,10)
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
load(file = "data/sime/sim.specs.rData")

#do the rpsft

fsim <- function(parms){
    df <- do.call(sim.GenStudy, parms)
    
    # IPE
    ipe.input <- transmute(df, 
                           event.time = os.t, 
                           censor.ind = os.e, 
                           trt.ind = x.trt, 
                           switch.time = t.switch, 
                           switch.ind = x.switch, 
                           cutofftime = t.censor
    )
    this.mipe <- IPE_BW(ipe.input, max.iter = 100)
    this.mipe.cx <- MIPE.cox(ipe.input, this.mipe)
    # conv stat
    rc <- c(df$simid[1],
            exp(coef(this.mipe.cx$cox.ipe)),
            exp(confint(this.mipe.cx$cox.ipe)),
            this.mipe.cx$phi.unique,
            this.mipe.cx$weib.hr)
    names(rc) <- c("simid", "hr.mipe", "cil.mipe", "ciu.mipe", "conv.mipe", "hr2.mipe")
    # get weibull es
    rc

    return(rc)
}

fsim.lst <- function(lst,pfx){ 
  n.batch <- length(lst)/10
  for (batch in 1:10){
    rc.list <- mclapply(lst[((batch-1)*n.batch+1):(batch*n.batch)],fsim, mc.cores = cores)
    rc.df   <- as.data.frame(matrix(unlist(rc.list),nrow = length(rc.list), byrow=TRUE))
    names(rc.df) <- names(rc.list[[1]])
    assign(paste(pfx,batch,sep="_"), rc.df)
    save(list=paste(pfx,batch,sep="_"), file = paste("data/raw/sim",pfx,batch,"resmipe.rData", sep = "_"))
  }
}
  
fsim.lst(sim2.spec.list,pfx = "sim2")
fsim.lst(sim3.spec.list,pfx = "sim3")
fsim.lst(sim4.spec.list,pfx = "sim4")
