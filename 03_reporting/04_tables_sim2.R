# Purpose: Creates tex summary tables for body and appendix for sim 2
# Functions: n/a
# Inputs: 
#   data/derived/res_all_sim234.rData
# Outputs: 
#   TeX/tables/t_sim2app_guide.tex
#   TeX/tables/t_sim2app_full.tex
#   TeX/tables/t_sim2res.tex

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
                           ), rho >= 0, sim ==1)


###################################
# make a table


guide.df <- data.frame(expand.grid(rh=seq(from=0, to=1, by =0.2), hr = c(0.7,0.9,1)))
guide.df <- mutate(guide.df, t = 1:nrow(guide.df))
guide.df <- mutate(guide.df, tab = paste("Table \\ref{T:app.sim2.res.",t,"}", sep =""))
guide.df<- mutate(guide.df, desc = paste("Results for study 1 with true HR $\\exp(\\beta_{1a}) =",hr,"$ and correlation $\\rho =", rh,"$",sep="") )

guid2.df <- select(guide.df, tab, desc)
colnames(guid2.df) <- c("Table", "Title")

guide.xt <- xtable(guid2.df)
print(guide.xt,
      include.rownames = FALSE, 
      hline.after = c(-1,-1,0,which(guide.df$rh==1)), 
      size = "small", 
      caption.placement = "top", 
      sanitize.text.function = function(x) {x},
      file = "TeX/tables/t_sim2app_guide.tex",
      append = FALSE,
      floating = FALSE
)


mtorep <- c(
"ITT",                "PP-EX",              "PP-CENS",            "TVC",    "TVC2",            "MIPE",           "MIPE-WEIB",
"RPSFT-TG-LR",        "RPSFT-OT-LR-TU",     "RPSFT-OT-LR-V",     "RPSFT-TG-W",         "RPSFT-OT-W-TU",      "RPSFT-OT-W-V",
"RPSFT-TG-W-ALL",         "RPSFT-OT-W-TU-ALL",      "RPSFT-OT-W-V-ALL",
"MIPE-ALL",           "2SAFT", "2SAFT-NULL","2SAFT-NC", "2SAFT-NC-NULL",
"RPSFT-TG-LR-ALL",    "RPSFT-OT-LR-TU-ALL", "RPSFT-OT-LR-V-ALL")

make.sim2.res.xt <- function(sel.rho, sel.HR, tnum){

  rt <- filter(rep.table, 
               abs(rho - sel.rho) <=0.01, 
               abs(true.hr - sel.HR) <=0.01,
               type %in% mtorep
               )
  
  rt <- arrange(rt, p.switch, type)
  
  rt1 <- transmute(rt, psw = paste(p.switch*100,"\\%",sep =""), type, mnHR=mnhr, pbias, mse, cov, pconv)
  
  rn2 <- "\\\\ $(\\%)$ & & $\\exp(\\hat{\\beta})$ &  (\\%) & & $(\\%)$ & $(\\%)$"
  
  colnames(rt1) <- c("Switch ", "Method", "Mean Est.", "Bias", "MSE", " Coverage ", paste("Convergence",rn2))

  cap <- paste("Simulation Study 1 Results with True HR =",sel.HR,", $\\rho =", sel.rho,"$ \\label{T:app.sim2.res.",tnum,"}",sep="")

  xt <- xtable(rt1, caption = cap)
  digits(xt) <- c(0,0,0,2,2,2,1,1)
  print(xt,
        include.rownames = FALSE, 
        hline.after = c(-1,-1,0,which(rt$type %in% c("TVC2", "RPSFT-OT-W-V-ALL", "RPSFT-OT-LR-V-ALL", "MIPE-WEIB"))), 
        scalebox = 0.6, 
        caption.placement = "top", 
        sanitize.text.function = function(x) {x},
        file = fname,
        append = TRUE
        )
}


fname <-  "TeX/tables/t_sim2app_full.tex"

write("%the rest of the tables", file = fname, append = "FALSE")

for (i in 1:nrow(guide.df)){
  make.sim2.res.xt(sel.rho = guide.df$rh[i], sel.HR = guide.df$hr[i], tnum = i)
  write("\\clearpage", file = fname, append = "TRUE")
}

#################################################################################
## get min and max bias per method
mtosum <- c( "ITT", "PP-EX", "PP-CENS", "TVC", "TVC2", "MIPE", "MIPE-WEIB",  
             "RPSFT-TG-LR", "RPSFT-OT-LR-TU", "RPSFT-OT-LR-V","MIPE-ALL",  "RPSFT-TG-LR-ALL", "RPSFT-OT-LR-TU-ALL", 
             "RPSFT-OT-LR-V-ALL","2SAFT","2SAFT-NULL","2SAFT-NC","2SAFT-NC-NULL")

r2 <- summarise(group_by(filter(rep.table, type %in% mtosum), type), 
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
r2.df <- select(r2, type, ranb,  ranm, ranc, rancv)


t2 <- "\\\\   & $(\\%)$ &   & $(\\%)$ & $(\\%)$"
t3 <- "\\\\   & min, max & min, max & min, max & min, max"
colnames(r2.df) <- c("Method", "Bias ", "MSE" ,"Coverage", paste("Convergence ",t2,t3))

                     

r2.xt <- xtable(r2.df)
align(r2.xt) <- c("l", "l", "r", "r", "r", "r")
print(r2.xt,
      include.rownames = FALSE, 
      hline.after = c(-1,-1,0,which(r2$type %in% c("TVC2","RPSFT-TG-LR-ALL", "RPSFT-OT-LR-V-ALL", "2SAFT-NC-NULL", "MIPE-WEIB"))), 
      scalebox = 0.9, 
      floating = FALSE,
      caption.placement = "top", 
      sanitize.text.function = function(x) {x},
      file = "TeX/tables/t_sim2res.tex",
      append = FALSE
)
