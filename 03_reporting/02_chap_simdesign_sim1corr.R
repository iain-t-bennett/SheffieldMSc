# Purpose: Generates the a plot showing correlations for simulation 1
# Functions: 
#   mscTheme
# Inputs: 
#   n/a
# Outputs: 
#   TeX/images/chap_simdesign/sim1corr.png

require("dplyr")
require("ggplot2")

rm(list = ls())

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

os.lambda = 0.3
os.gamma = 1.2
pfs.lambda = 2
pfs.gamma = 1.5
rho = c(-.8,0,.8)

n.arm <-500

u <- seq(from = 0, to= 1, by = 1/n.arm)
u<-u[order(u)]
  
ur <- runif(length(u))

os   <- (-(log(u) / (os.lambda)))^(1/os.gamma)
ttp  <- (-(log(u) / (pfs.lambda)))^(1/pfs.lambda)
  
  pfs <- matrix(nrow = length(os), ncol = length(rho))
  ttp2 <- pfs
  os <- pmin(os,3)
  ttp <- pmin(ttp,3)
  
  for (i in 1:length(rho)){
    r <- rho[i]
    u_corr <-  pnorm(r* qnorm(u,0,1) + ((1-r^2)^0.5)*qnorm(ur,0,1),0,1)
    ttp2[,i] <- (-(log(u_corr) / (pfs.lambda)))^(1/pfs.gamma)
    pfs[,i] <- pmin(os,ttp2[,i])
    
  }
  
  rhol <- paste("rho = ",rho)
  
  df <- data.frame(os,ttp = ttp,ttp2=ttp2[,1],pfs=pfs[,1],id=1:length(os)/length(os), rho=rho[1])
  
  df2 <- data.frame(id=1:length(os)/length(os),pfsx=pfs[order(pfs[,1],decreasing = TRUE),1],rho=rho[1])
  
  for (i in 2:length(rho)){
    df <- rbind(df,
                data.frame(os,ttp,ttp2= ttp2[,i],pfs=pfs[,i],id=1:length(os)/length(os),rho=rho[i])
    )
    
    df2<- rbind(df2, 
                data.frame(id=1:length(os)/length(os),pfsx=pfs[order(pfs[,i],decreasing = TRUE),i],rho=rho[i])
    )
  }
  
  df$rho <- factor(df$rho, levels = rho, labels = paste("Correlation:",rho))
  df2$rho <- factor(df2$rho, levels = rho, labels = paste("Correlation:",rho))
  
  cols <- c("red","blue","black")
  
  
  p1 <- ggplot()  + 
    geom_ribbon(data = df,aes(x=id, ymin = 0, ymax = pfs,fill="PFS"), alpha = 0.7) +
    geom_ribbon(data = df,aes(x=id, ymin = pfs, ymax = ttp2,fill="TTP"), alpha = 0.7) +
    geom_ribbon(data = df,aes(x=id, ymin = pfs, ymax = os,fill="OS"), alpha = 0.3) +
    geom_line(data = df2, aes(x=id,y=pfsx,color ="PFS"), size = 1.3) +
    geom_line(data = df, aes(x=id,y=os,color = "OS"),size = 1.3) +
    geom_line(data = df, aes(x=id,y=ttp,color = "TTP"), size = 1) +
    mscTheme() +
    facet_wrap(~rho, ncol = 1) + 
    theme(legend.position = "bottom",
          axis.line.y = element_blank(),
          axis.line.x = element_blank()
          ) + 
    coord_flip(ylim = c(0,3)) +
    scale_x_continuous(labels = scales::percent) +
    geom_vline(xintercept = 0, size = 1) + 
    geom_hline(yintercept = 0, size = 1) +
    
    xlab("Survival") + 
    ylab("Time") +
    scale_colour_manual(name="Endpoint",values=cols) +
    scale_fill_manual(name="Endpoint",values=cols)

p1
  
ggsave("TeX/images/chap_simdesign/sim1corr.png", width = 12, height = 18, plot = p1)
