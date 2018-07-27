# Purpose: Plot of correlation between PFS OS dependent on correlation between TTP and OS
# Functions: 
#   sim.GenStudy
#   mscTheme
# Inputs: 
#   n/a
# Outputs: 
#   TeX/images/chap_simdesign/PFSOScor.png


##### create the summary tables

rm(list=ls())
date()

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

require(dplyr)
require(survival)
require(ggplot2)

###### get km est at 6m for pfs

f3  <- function(rho){
  nt = 200
  df.list <- lapply(floor(runif(nt,0,100000)),sim.GenStudy,arm.n = 500, rho = rho, prop.pd.int = 0, simid =99,p.switch = 0, beta_1a = log(1),beta_1b = log(1),beta_2a = log(1), beta_2b=log(1))
  f4 <- function(t.df){
    sf.pfs <-survfit(Surv(pfs.t,pfs.e)~1,data = t.df)
    sf.os <-survfit(Surv(os.t, os.e)~1,data = t.df)
    km.pfs <- sf.pfs$surv[min(which(sf.pfs$time-0.5>0))-1]
    km.os  <- sf.os$surv[min(which(sf.os$time-1>0))-1]
    rc=c(km.pfs, km.os)
  }
  rc <- lapply(df.list,f4)
  rc2 <- as.data.frame(matrix(unlist(rc),nrow = nt, byrow=TRUE))
  names(rc2)<- c("pfs","os")
  rc3 <- mutate(rc2,rho=rho)
  return(rc3)
}

rho <- seq(from = 0, to = 1, by = 0.1)
corr.l <- lapply(rho,f3)

all.df <- rbind(corr.l[[1]],corr.l[[2]],corr.l[[3]],corr.l[[4]],corr.l[[5]],corr.l[[6]],corr.l[[7]],corr.l[[8]],corr.l[[9]],corr.l[[10]],corr.l[[11]])


ggplot(all.df) + geom_point(aes(x=pfs,y=os))+facet_wrap(~rho, scales = "free", nrow = 3) + xlab("KM estimate of PFS at 6 months") + ylab("KM estimate of OS at 12 months")




plot.df <- summarise(group_by(all.df,rho), c= cor(pfs,os))

p1 <- ggplot(plot.df, aes(x=rho,y=c)) + 
  geom_hline(yintercept = 0, size = 1) +
  geom_vline(xintercept = 0, size = 1) +
  geom_point(size = 3) + geom_line(linetype = 2, size = 1) +
  mscTheme() + 
  xlab("Correlation between TTP and OS") + ylab("Estimated correlation between PFS and OS") + 
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.line.x = element_blank()) +
  theme(axis.line.y = element_blank()) 

ggsave(filename = "TeX/images/chap_simdesign/PFSOScor.png",plot = p1, width = 8, height = 8)

