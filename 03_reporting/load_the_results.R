# Purpose: Combines all the following analysis results in a single data frame for more analysis
#          ITT
#          Per Protocol
#          Censoring at Switch
#          Treatment as a time varying covariate
#
#          On treatment RPSFT
#          - log rank test
#          - wilcoxon test
#          Treatment Group RPSFT
#          - log rank test
#          - wilcoxon test
#         
#          MIPE
# 
#          Two-stage AFT without recensoring
#
#          Two-stage AFT with recensoring
# Functions: n/a
# Inputs: 
#   data/sims/sim.specs.rData
#   data/raw/sim.simple.res.rData
#   data/raw/sim_sim2_nn_resrpsft.rData (nn = 1,2,..,10)
#   data/raw/sim_sim3_nn_resrpsft.rData (nn = 1,2,..,10)
#   data/raw/sim_sim4_nn_resrpsft.rData (nn = 1,2,..,10)
#   data/raw/sim_sim2_nn_resmipe.rData (nn = 1,2,..,10)
#   data/raw/sim_sim3_nn_resmipe.rData (nn = 1,2,..,10)
#   data/raw/sim_sim4_nn_resmipe.rData (nn = 1,2,..,10)
#   data/raw/sim.2saft.res.rData
#   data/raw/sim.2saft.res.cens.rData
# Outputs: 
#   data/derived/res_all_sim234.rData

rm(list = ls())

# load the specs
load("data/sims/sim.specs.rData")
# load the simple results
load("data/raw/sim.simple.res.rData")
# load the 2saft results
load("data/raw/sim.2saft.res.rData")
# load the 2saft results censoring
load("data/raw/sim.2saft.res.cens.rData")

# modify names of censoring 
names(sim.2saft.res.cens.df) <- c("simid", "hr.2saft.null.c", "cil.2saft.null.c", "ciu.2saft.null.c", "nev.2saft.null.c", "conv.2saft.null.c",
                                  "hr.2saft.corr.c", "cil.2saft.corr.c",  "ciu.2saft.corr.c", "nev.2saft.corr.c",  "conv.2saft.corr.c")



load("data/derived/res_all_sim234.rData")
ores.df <- res.df
ores.table <- res.table

# load all the RPSFT and MIPE results (batches 1-10)


res.df <- NULL
for (i in 2:4){
  for (j in 1:10){
    mipe.fname <- paste("data/raw/sim_sim",i,"_",j,"_resmipe.rData", sep ="")
    rpsft.fname <- paste("data/raw/sim_sim",i,"_",j,"_resrpsft.rData", sep ="")
    dname <- paste("sim",i,"_",j, sep ="")
    load(mipe.fname)
    mipe.res <- get(dname)
    load(rpsft.fname)
    rpsft.res <- get(dname)
    assign(dname, merge(mipe.res, rpsft.res, by = "simid"))
    if (is.null(res.df)){
      res.df <- get(dname)
    } else {
      res.df <- rbind(res.df, get(dname))
    }
  }
}

res0.df <- merge(sim.simple.res.df,res.df, by = "simid")

res1a.df <- merge(res0.df, sim.2saft.res.df, by = "simid")
res1.df  <- merge(res1a.df, sim.2saft.res.cens.df, by = "simid")

# transpose

require(dplyr)

itt.df    <- transmute(res1.df, simid, type = "ITT",            hr = hr.itt,              cil = cil.itt,         ciu = ciu.itt,         conv = 1,          rng = NA)
pp.df     <- transmute(res1.df, simid, type = "PP-EX",         hr = hr.pp,               cil = cil.pp,          ciu = ciu.pp,          conv = 1,          rng = NA)
cens.df   <- transmute(res1.df, simid, type = "PP-CENS",        hr = hr.cens,             cil = cil.cens,        ciu = ciu.cens,        conv = 1,          rng = NA)
tvc.df    <- transmute(res1.df, simid, type = "TVC",            hr = hr.tvc,              cil = cil.tvc,         ciu = ciu.tvc,         conv = 1,          rng = NA)
tvc2.df   <- transmute(res1.df, simid, type = "TVC2",           hr = hr.tvc2,             cil = cil.tvc2,        ciu = ciu.tvc2,        conv = 1,          rng = NA)

mipe.df   <- transmute(res1.df, simid, type = "MIPE",           hr = hr.mipe,             cil = cil.mipe,        ciu = ciu.mipe,        conv = conv.mipe,  rng = NA)
mipe2.df  <- transmute(res1.df, simid, type = "MIPE-WEIB",      hr = hr2.mipe,            cil = NA,              ciu = NA,              conv = conv.mipe,  rng = NA)

tglr.df   <- transmute(res1.df, simid, type = "RPSFT-TG-LR",    hr = tg.lr.hr.trt.ind,    cil = tg.lr.hr.cil,    ciu = tg.lr.hr.ciu,    conv = tg.lr.conv, rng = tg.lr.rng)
tgw.df    <- transmute(res1.df, simid, type = "RPSFT-TG-W",     hr = tg.w.hr.trt.ind,     cil = tg.w.hr.cil,     ciu = tg.w.hr.ciu,     conv = tg.w.conv,  rng = tg.w.rng)

otlrv.df  <- transmute(res1.df, simid, type = "RPSFT-OT-LR-V",  hr = ot.lr.v.hr.trt.ind,  cil = ot.lr.v.hr.cil,  ciu = ot.lr.v.hr.ciu,  conv = ot.lr.conv, rng = ot.lr.rng)
otlrtu.df <- transmute(res1.df, simid, type = "RPSFT-OT-LR-TU", hr = ot.lr.tu.hr.trt.ind, cil = ot.lr.tu.hr.cil, ciu = ot.lr.tu.hr.ciu, conv = ot.lr.conv, rng = ot.lr.rng)

otwv.df   <- transmute(res1.df, simid, type = "RPSFT-OT-W-V",   hr = ot.w.v.hr.trt.ind,   cil = ot.w.v.hr.cil,   ciu = ot.w.v.hr.ciu,   conv = ot.lr.conv, rng = ot.lr.rng)
otwtu.df  <- transmute(res1.df, simid, type = "RPSFT-OT-W-TU",  hr = ot.w.tu.hr.trt.ind,  cil = ot.w.tu.hr.cil,  ciu = ot.w.tu.hr.ciu,  conv = ot.lr.conv, rng = ot.lr.rng)

saft1.df  <- transmute(res1.df, simid, type = "2SAFT-NULL",     hr = hr.2saft.null.c,       cil = cil.2saft.null.c,  ciu = ciu.2saft.null.c,  conv = conv.2saft.null.c, rng = NA)
saft2.df  <- transmute(res1.df, simid, type = "2SAFT",          hr = hr.2saft.corr.c,       cil = cil.2saft.corr.c,  ciu = ciu.2saft.corr.c,  conv = conv.2saft.corr.c, rng = NA)

saft3.df  <- transmute(res1.df, simid, type = "2SAFT-NC-NULL",  hr = hr.2saft.null,       cil = cil.2saft.null,  ciu = ciu.2saft.null,  conv = conv.2saft.null, rng = NA)
saft4.df  <- transmute(res1.df, simid, type = "2SAFT-NC",       hr = hr.2saft.corr,       cil = cil.2saft.corr,  ciu = ciu.2saft.corr,  conv = conv.2saft.corr, rng = NA)

res2.df <- rbind(itt.df, pp.df, cens.df, tvc.df, tvc2.df, mipe.df,  mipe2.df, tglr.df, tgw.df, otlrv.df, otlrtu.df, otwv.df, otwtu.df, saft1.df, saft2.df, saft3.df, saft4.df)

res3.df <- merge(res2.df, filter(spec.df, sim !=1), by = "simid")

# get some stats

res4.df <- mutate(res3.df, 
                  bad   = ifelse(hr < 0.001 | hr > 1000 | is.na(hr), 1, 0),
                  hr    = ifelse(bad ==1, NA, hr),
                  cil   = ifelse(bad ==1, NA, cil),
                  ciu   = ifelse(bad ==1, NA, ciu),
                  conv  = ifelse(bad ==1, 0,  conv),
                  trend = as.numeric(hr < 1),
                  inci  = as.numeric(true.hr > cil & true.hr < ciu)
                  )

# derive version including all for RPSFT and MIPE

res5.df <- filter(res4.df, type %in% c("RPSFT-TG-LR",
                                       "RPSFT-TG-W",
                                       "RPSFT-OT-LR-TU",
                                       "RPSFT-OT-LR-V",
                                       "RPSFT-OT-W-V", 
                                       "RPSFT-OT-W-TU",  
                                       "MIPE"))
                                       
                                       
res6.df <- mutate(res5.df,                      type = paste(type,"ALL",sep="-"), conv = !is.na(hr))
res7.df <- mutate(filter(res5.df,type!="MIPE"), type = paste(type,"TOL",sep="-"), conv = !is.na(hr) & rng < 0.2)

res8.df <- rbind(res4.df, res6.df, res7.df)



types <- c(
  "ITT",
  "PP-CENS",            "PP-EX",
  "TVC"     ,           "TVC2"  ,
  
  "RPSFT-TG-LR" ,       "RPSFT-TG-LR-ALL" ,   "RPSFT-TG-LR-TOL"    ,
  "RPSFT-TG-W"   ,      "RPSFT-TG-W-ALL"  ,   "RPSFT-TG-W-TOL"     ,
  
  "RPSFT-OT-LR-TU"    , "RPSFT-OT-LR-TU-ALL", "RPSFT-OT-LR-TU-TOL" ,
  "RPSFT-OT-LR-V"      ,"RPSFT-OT-LR-V-ALL" , "RPSFT-OT-LR-V-TOL"  ,
  "RPSFT-OT-W-TU",      "RPSFT-OT-W-TU-ALL" , "RPSFT-OT-W-TU-TOL"  ,
  "RPSFT-OT-W-V"  ,      "RPSFT-OT-W-V-ALL",   "RPSFT-OT-W-V-TOL" ,
  
  "2SAFT"     ,         "2SAFT-NULL" ,
  "2SAFT-NC"  ,         "2SAFT-NC-NULL" ,
  "MIPE"     ,          "MIPE-ALL"     ,      "MIPE-WEIB"         
)


res9.df <- mutate(ungroup(res8.df), type = factor(type, levels = types, ordered = TRUE))



res.df1 <- ungroup(res9.df)

res.table1 <- summarise(group_by(filter(res.df1, conv ==1), sim, scenario, p.switch, rho, true.hr,type), 
                       mnhr = mean(hr,na.rm = TRUE), 
                       cov   = mean(inci, na.rm = TRUE),
                       nconv = sum(conv),
                       SEmean = sd(hr, na.rm = TRUE)
                       )

res.table2 <- mutate(res.table1,
                     bias = mnhr - true.hr, 
                     pbias = bias/true.hr,
                     pconv  = nconv / 1000,
                     mse = (bias)^2 + SEmean^2 
                     )


res.table3 <- ungroup(res.table2)


# define labels

ids <- unique(transmute(res.df1, sim, scenario, e1a = exp(beta_1a),  e1b = exp(beta_1b),  e2a = exp(beta_2a), e2b = exp(beta_2b)))
ids <- mutate(ids, 
              newsim = ifelse(sim==4 & scenario ==3,  0, (sim!=2) + 1), 
              newscenario = as.numeric(newsim !=0))
ids <- mutate(ids, 
              newscenario = ifelse(sim ==  2, scenario, newscenario),
              lbl1 = ifelse(sim == 2, paste("True Hazard Ratio: ", e1a, sep = ""), "NA"))
ids <- mutate(ids, newscenario = ifelse(sim == 3 & scenario == 3, 1 ,newscenario))
ids <- mutate(ids, newscenario = ifelse(sim == 3 & scenario == 4, 2 ,newscenario))
ids <- mutate(ids, newscenario = ifelse(sim == 3 & scenario == 1, 3 ,newscenario))
ids <- mutate(ids, newscenario = ifelse(sim == 3 & scenario == 2, 4 ,newscenario))
ids <- mutate(ids, newscenario = ifelse(sim == 4, scenario+4, newscenario))

ids <- mutate(ids, lbl1 = ifelse(newsim == 2, paste("Scenario", newscenario)  , lbl1))

res.table4 <- merge(res.table3, ids, by = c("sim","scenario"))
res.df2 <- merge(res.df1, ids, by = c("sim","scenario"))

res.table5 <- mutate(res.table4, lbl = paste(lbl1, "\nProportion switch: ", p.switch*100,"%", sep = ""))
res.df3    <- mutate(res.df2,    lbl = paste(lbl1, "\nProportion switch: ", p.switch*100,"%", sep = ""))

ids2 <- unique(transmute(filter(res.table5, newsim !=0), newsim, newscenario, p.switch, lbl))
labs <- ids2$lbl[order(ids2$newsim,ids2$newscenario,ids2$p.switch)]

res.table6 <- mutate(res.table5, lbl = factor(lbl, labs, ordered = TRUE))
res.df4    <- mutate(res.df3,    lbl = factor(lbl, labs, ordered = TRUE))


res.table <- res.table6
res.df    <- res.df4

# save it
save(res.df, res.table, file = "data/derived/res_all_sim234.rData")


