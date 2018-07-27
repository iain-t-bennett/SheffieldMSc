# Purpose: Creates tex summary tables for body and appendix for sim 3
# Functions: n/a
# Inputs: 
#   data/derived/res_all_sim234.rData
# Outputs: 
#   TeX/tables/t_sim3app_guide.tex
#   TeX/tables/t_sim3app_full.tex
#   TeX/tables/t_sim3res.tex

rm(list=ls())
date()

##### load the results
load(file = "data/derived/res_all_sim234.rData")


require(dplyr)
require(xtable)

rep.table <- filter(mutate(res.table, sim = newsim, scenario = newscenario, 
                           pbias = pbias*100, 
                           cov = cov*100,
                           pconv = pconv * 100
                           ), rho >= 0, sim ==2)


###################################
# make a table


guide.df <- data.frame(expand.grid(rh=seq(from=0, to=1, by =0.2), sc = 1:6))
guide.df <- mutate(guide.df, t = 1:nrow(guide.df))
guide.df <- mutate(guide.df, tab = paste("Table \\ref{T:app.sim3.res.",t,"}", sep =""))
guide.df<- mutate(guide.df, desc = paste("Results for Study 2 Scenario ",sc," and correlation $\\rho =", rh,"$",sep="") )

guid2.df <- select(guide.df, tab, desc)
colnames(guid2.df) <- c("Table", "Title")

guide.xt <- xtable(guid2.df)
print(guide.xt,
      include.rownames = FALSE, 
      hline.after = c(-1,-1,0,which(guide.df$rh==1)), 
      caption.placement = "top", 
      sanitize.text.function = function(x) {x},
      file = "TeX/tables/t_sim3app_guide.tex",
      append = FALSE,
      floating = FALSE
)


mtorep <- c(
"ITT",                "PP-EX",              "PP-CENS",            "TVC",    "TVC2",            "MIPE",           "MIPE-WEIB",
"RPSFT-TG-LR",        "RPSFT-OT-LR-TU",     "RPSFT-OT-LR-V",     "RPSFT-TG-W",         "RPSFT-OT-W-TU",      "RPSFT-OT-W-V",
"RPSFT-TG-W-TOL",         "RPSFT-OT-W-TU-TOL",      "RPSFT-OT-W-V-TOL",
"MIPE-ALL",           "2SAFT", "2SAFT-NULL","2SAFT-NC", "2SAFT-NC-NULL",
"RPSFT-TG-LR-TOL",    "RPSFT-OT-LR-TU-TOL", "RPSFT-OT-LR-V-TOL")


make.sim3.res.xt <- function(sel.rho, sel.sc, tnum){

  rt <- filter(rep.table, 
               abs(rho - sel.rho) <=0.01, 
               abs(scenario - sel.sc) <=0.01,
               type %in% mtorep
               )
  
  rt <- arrange(rt, p.switch, type)
  
  true.hr <- round(rt$true.hr[1],2)
  
  rt1 <- transmute(rt, psw = paste(p.switch*100,"\\%",sep =""), type, mnHR=mnhr, pbias, mse, cov, pconv)
  
  rn2 <- "\\\\ $(\\%)$ & & $\\exp(\\hat{\\beta})$ &  (\\%) & & $(\\%)$ & $(\\%)$"
  
  colnames(rt1) <- c("Switch ", "Method", "Mean Est.", "Bias", "MSE", " Coverage ", paste("Convergence",rn2))

  cap <- paste("Simulation Study 2 Scenario ",sel.sc," Results with $\\rho =", sel.rho,"$ (True HR = ",true.hr,") \\label{T:app.sim3.res.",tnum,"}",sep="")

  xt <- xtable(rt1, caption = cap)
  digits(xt) <- c(0,0,0,2,2,2,1,1)
  print(xt,
        include.rownames = FALSE, 
        hline.after = c(-1,-1,0,which(rt$type %in% c("TVC2", "RPSFT-OT-W-V-TOL", "RPSFT-OT-LR-V-TOL", "MIPE-WEIB"))), 
        scalebox = 0.7, 
        caption.placement = "top", 
        sanitize.text.function = function(x) {x},
        file = fname,
        append = TRUE
        )
}


fname <-  "TeX/tables/t_sim3app_full.tex"

write("%the rest of the tables", file = fname, append = "FALSE")

for (i in 1:nrow(guide.df)){
  make.sim3.res.xt(sel.rho = guide.df$rh[i], sel.sc = guide.df$sc[i], tnum = i)
  write("\\clearpage", file = fname, append = "TRUE")
}

#################################################################################
## get min and max bias per method
mtosum <- c( "ITT", "MIPE-ALL", "MIPE-WEIB",  
             "RPSFT-TG-LR", "RPSFT-OT-LR-TU-TOL", 
             "RPSFT-OT-LR-V-TOL","2SAFT", "PP-CENS", "TVC", "TVC2")

rt <- mutate(rep.table, scs = ifelse(scenario %in% c(1,2), "1 and 2", ifelse(scenario %in% c(3,4), "3 and 4", "5 and 6")))


r2 <- summarise(group_by(filter(rt, type %in% mtosum), scs, type), 
                minb = min(pbias), maxb= max(pbias),
                minc = min(cov), maxc = max(cov), 
                minh = min(mnhr), maxh = max(mnhr), 
                minm = min(mse), maxm = max(mse),
                mincv = min(pconv), maxcv = max(pconv)
                )
r2 <- mutate(r2, 
             ranb = paste(format(minb, nsmall = 2, digits = 0, scientific = FALSE),", ",format(maxb, nsmall = 2, digits = 0, scientific = FALSE), sep = ""),
             ranc = paste(format(minc, nsmall = 1, digits = 0, scientific = FALSE),", ",format(maxc, nsmall = 1, digits = 0, scientific = FALSE), sep = ""),
             ranh = paste(format(minh, nsmall = 2, digits = 0, scientific = FALSE),", ",format(maxh, nsmall = 2, digits = 0, scientific = FALSE), sep = ""),
             ranm = paste(format(minm, nsmall = 2, digits = 0, scientific = FALSE),", ",format(maxm, nsmall = 2, digits = 0, scientific = FALSE), sep = ""),
             rancv = paste(format(mincv, nsmall = 1, digits = 0, scientific = FALSE),", ",format(maxcv, nsmall = 1, digits = 0, scientific = FALSE), sep = "")
                          )

r2$ranc[r2$type %in% c("MIPE-WEIB")] <- ""
r2.df <- select(r2, scs, type, ranb,  ranm, ranc, rancv)


t2 <- "\\\\  &  & $(\\%)$ &   & $(\\%)$ & $(\\%)$"
t3 <- "\\\\  &  & min, max & min, max & min, max & min, max"
colnames(r2.df) <- c("Scenarios", "Method", "Bias ", "MSE" ,"Coverage", paste("Convergence ",t2,t3))


r2.xt <- xtable(r2.df)
align(r2.xt) <- c("l", "l", "r", "r", "r", "r", "r")
print(r2.xt,
      include.rownames = FALSE, 
      hline.after = c(-1,-1,0,which(r2$type %in% c("MIPE-WEIB"))), 
      scalebox = 0.9, 
      floating = FALSE,
      caption.placement = "top", 
      sanitize.text.function = function(x) {x},
      file = "TeX/tables/t_sim3res.tex",
      append = FALSE
)
