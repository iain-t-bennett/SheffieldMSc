# Purpose: Creates plots for body and appendix of simulation results
# Functions: 
#   mscTheme
# Inputs: 
#   data/derived/res_all_sim234.rData
# Outputs: 
#   TeX/images/chap_sim2/simple_bias.png
#   TeX/images/chap_sim2/2saft_bias.png
#   TeX/images/chap_sim2/complex_bias.png
#   TeX/images/chap_sim2/simple_cov.png
#   TeX/images/chap_sim2/complex_cov.png
#   TeX/images/chap_sim2/complex_conv.png
#
#   TeX/images/app_sim2res/complexA_cov.png
#   TeX/images/app_sim2res/complexA_bias.png
# 
#   TeX/images/chap_sim3/conv_cdf.png
#   TeX/images/chap_sim3/simple_bias.png
#   TeX/images/chap_sim3/comp_bias12.png
#   TeX/images/chap_sim3/comp_bias34.png
#   TeX/images/chap_sim3/comp_bias56.png
#   TeX/images/chap_sim3/2saft_bias.png
#   TeX/images/chap_sim3/simple_cov.png
#   TeX/images/chap_sim3/comp_cov.png
#   TeX/images/chap_sim3/comp_conv.png
#   
#   TeX/images/app_allres/simp_bias4.png
#   TeX/images/app_allres/rpsft_bias4.png
#   TeX/images/app_allres/saft_bias4.png
#   TeX/images/app_allres/mipe_bias4.png
#

rm(list=ls())
date()

##### load the results
load(file = "data/derived/res_all_sim234.rData")

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

require(dplyr)
require(ggplot2)

plot.df <- filter(mutate(res.table, sim = newsim, scenario = newscenario), rho >= 0)

# plots of mean value
###################################################################???

makebiasplot <- function(rt, nc, lnc = 3) {
  
  g <- guide_legend(ncol = lnc, byrow = TRUE)
  
  rc.p <- ggplot(data = rt) + mscTheme() +
    geom_point(aes(
      x = rho,
      y = pbias,
      color = type,
      shape = type
    ),size = 4
    ) +
    geom_line(aes(x = rho, y = pbias, color = type, linetype = type), size = 1) +
    facet_wrap( ~ lbl, scales = "free", ncol = nc) +
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = 0)  +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    ylab("Bias (%)") + xlab(expression(paste("Correlation (", rho, ")"))) +
    theme(legend.key.size  = unit(5, "lines")) +
    guides(colour = g, linetype = g, shape = g)
  
  return(rc.p)
}

makecovrgplot <- function(rt, nc, lnc = 3) {
  
  g <- guide_legend(ncol = lnc, byrow = TRUE)
  
  
  rc.p <- ggplot(data = rt) + mscTheme() +
    geom_point(aes(
      x = rho,
      y = cov,
      color = type,
      shape = type
    ),size = 4, position = position_dodge(width = 0.08)
    ) +
    geom_line(aes(x = rho, 
                  y = cov, 
                  color = type,
                  linetype = type), 
              size = 1, 
              position = position_dodge(width = 0.08)
    ) +
    facet_wrap( ~ lbl, scales = "free", ncol = nc) +
    geom_hline(yintercept = 0.95, linetype = 2, size = 1)  +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) + 
    scale_y_continuous(labels = scales::percent) +
    ylab("Coverage (%)") + 
    xlab(expression(paste("Correlation (", rho, ")"))) +
    theme(legend.key.size  = unit(5, "lines")) +
    guides(colour = g, linetype = g, shape = g)
  
  
  return(rc.p)
}

makeconvplot <- function(rt, nc, lnc = 3) {

g <- guide_legend(ncol = lnc, byrow = TRUE)

rc <- ggplot(data = rt) + 
  mscTheme() +
  geom_point(aes(
    x = rho,
    y = pconv,
    color = type,
    shape = type
  ),size = 4
  ) +
  geom_line(aes(x = rho,
                y = pconv, 
                color = type,
                linetype = type), 
            size = 1
  ) +
  facet_wrap( ~ lbl, scales = "free", ncol = nc) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) + 
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of simulations where method converged (%)") + 
  xlab(expression(paste("Correlation (", rho, ")"))) +
  theme(legend.key.size  = unit(5, "lines")) +
  guides(colour = g, linetype = g, shape = g)


return(rc)
}

########3
## simulation 2

s1.simp.types <- c("ITT", "PP-CENS", "PP-EX", "TVC", "TVC2")
s1.comp.types <- c("MIPE-ALL", "MIPE-WEIB", "RPSFT-TG-LR", "RPSFT-OT-LR-V-ALL", "RPSFT-OT-LR-TU-ALL")
s1.2saf.types <- c("2SAFT", "2SAFT-NULL", "2SAFT-NC", "2SAFT-NC-NULL")

p.s1.simple    <- makebiasplot(filter(plot.df, sim == 1, type %in% s1.simp.types), nc = 3)
p.s1.complex   <- makebiasplot(filter(plot.df, sim == 1, type %in% s1.comp.types), nc = 3)
p.s1.2saft     <- makebiasplot(filter(plot.df, sim == 1, type %in% s1.2saf.types), nc = 3, lnc =2)
p.s1.complexA  <- makebiasplot(filter(plot.df, sim == 1, type %in% c("MIPE", "MIPE-ALL", "RPSFT-OT-LR-V","RPSFT-OT-LR-V-ALL", "RPSFT-OT-LR-TU","RPSFT-OT-LR-TU-ALL")), nc = 3)

p.s1.simpleC   <- makecovrgplot(filter(plot.df, sim == 1, type %in% s1.simp.types), nc = 3)
p.s1.complexC  <- makecovrgplot(filter(plot.df, sim == 1, type %in% c("MIPE-ALL", "2SAFT" ,"2SAFT-NC", "RPSFT-TG-LR", "RPSFT-OT-LR-V-ALL", "RPSFT-OT-LR-TU-ALL")), nc = 3)
p.s1.complexCA <- makecovrgplot(filter(plot.df, sim == 1, type %in% c("RPSFT-OT-LR-V","RPSFT-OT-LR-V-ALL", "RPSFT-OT-LR-TU","RPSFT-OT-LR-TU-ALL")), nc = 3)

dfX <- filter(plot.df, 
              type %in% c("MIPE", "MIPE-ALL", "RPSFT-TG-LR", "RPSFT-OT-LR-TU", "RPSFT-OT-LR-TU-ALL"),
              sim==1)

dfX <- mutate(dfX, 
             type = ifelse(type == "MIPE", "MIPE", 
                           ifelse(type == "MIPE-ALL", "MIPE-ALL", 
                                  ifelse(type == "RPSFT-TG-LR", "RPSFT-TG-LR", 
                                         ifelse(type == "RPSFT-OT-LR-TU", "RPSFT-OT-LR", "RPSFT-OT-LR-ALL")))))

dfX <- mutate(dfX, type = factor(type, levels = c("RPSFT-TG-LR", "RPSFT-OT-LR", "RPSFT-OT-LR-ALL","MIPE","MIPE-ALL"), ordered = TRUE))

p.s1.complexCV <- makeconvplot(dfX, nc = 3)


ggsave(filename = "TeX/images/chap_sim2/simple_bias.png",
       plot = p.s1.simple + coord_cartesian(ylim = c(-0.5,0.5)),
       height = 20, width = 14
)

ggsave(filename = "TeX/images/chap_sim2/2saft_bias.png",
       plot = p.s1.2saft + coord_cartesian(ylim = c(-.1,.1)),
       height = 20, width = 14
)

ggsave(filename = "TeX/images/chap_sim2/complex_bias.png",
       plot = p.s1.complex + coord_cartesian(ylim = c(-0.5,0.5)),
       height = 20, width = 14
)

ggsave(filename = "TeX/images/app_sim2res/complexA_bias.png",
       plot = p.s1.complexA + coord_cartesian(ylim = c(-0.5,0.5)),
       height = 20, width = 14
)


ggsave(filename = "TeX/images/chap_sim2/simple_cov.png",
       plot = p.s1.simpleC + coord_cartesian(ylim = c(0,1)),
       height = 20, width = 14
)

ggsave(filename = "TeX/images/chap_sim2/complex_cov.png",
       plot = p.s1.complexC + coord_cartesian(ylim = c(0.45,1)),
       height = 20, width = 14
       
)

ggsave(filename = "TeX/images/app_sim2res/complexA_cov.png",
       plot = p.s1.complexCA + coord_cartesian(ylim = c(0.45,1)),
       height = 20, width = 14
)

ggsave(filename = "TeX/images/chap_sim2/complex_conv.png",
       plot = p.s1.complexCV + coord_cartesian(ylim = c(0,1)),
       height = 20, width = 14
)

##########################################################
# sim3 - convergence

p.df <- filter(res.df, type == "RPSFT-OT-LR-TU", newsim==2, rho>=0)

g <- guide_legend(ncol = 6, byrow = TRUE, title = expression(paste("Correlation (", rho, ")")))

p.df <- mutate(
  p.df,
  rng      = ifelse(is.finite(rng),rng,5)
)

p.cv <- ggplot(p.df, aes(rng, colour = as.factor(rho), linetype = as.factor(rho))) + 
  stat_ecdf(size = 1) + 
  mscTheme() +
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_blank()
  )+
  geom_vline(xintercept = 0, size = 1) +
  facet_wrap(~lbl, scales = "free", ncol = 3) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") + 
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of simulations where method converged (%)") + 
  xlab("Tolerance of G-estimation convergence") +
  theme(legend.key.size  = unit(5, "lines")) +
  geom_hline(yintercept = 1, size = 1, linetype = 2) +
  guides(colour = g, linetype = g) +
  theme(legend.title = element_text())

ggsave(filename = "TeX/images/chap_sim3/conv_cdf.png",
       plot = p.cv + coord_cartesian(ylim = c(0,1), xlim = c(0,0.5)),
       height = 20, width = 14
)


#######################################################################
# sim 3 

s2.simp.types <- c("ITT", "PP-CENS", "PP-EX", "TVC", "TVC2")

s2.comp.types <- c("MIPE-ALL", "MIPE-WEIB", "RPSFT-TG-LR", "RPSFT-OT-LR-V-TOL", "RPSFT-OT-LR-TU-TOL")
s2.2saf.types <- c("2SAFT", "2SAFT-NULL", "2SAFT-NC", "2SAFT-NC-NULL")

s2.cmp.types <- c("MIPE-ALL", "MIPE-WEIB", "RPSFT-TG-LR", "RPSFT-OT-LR-TU-TOL","RPSFT-OT-LR-V-TOL", "2SAFT")


p.s2.sall.simple    <- makebiasplot(filter(plot.df, sim == 2, scenario <= 6, type %in% s2.simp.types), nc = 3)

ggsave(filename = "TeX/images/chap_sim3/simple_bias.png",
       plot = p.s2.sall.simple + coord_cartesian(ylim = c(-0.5,0.5)),
       height = 20, width = 14
)


p.s2.s12.cplx <- makebiasplot(filter(plot.df, sim == 2, scenario %in% c(1,2), type %in% s2.cmp.types), nc = 2, lnc = 3) 

p.s2.s34.cplx <- makebiasplot(filter(plot.df, sim == 2, scenario %in% c(3,4), type %in% s2.cmp.types), nc = 2, lnc = 3)

p.s2.s56.cplx <- makebiasplot(filter(plot.df, sim == 2, scenario %in% c(5,6), type %in% s2.cmp.types), nc = 2, lnc = 3)


ggsave(filename = "TeX/images/chap_sim3/comp_bias12.png",
       plot = p.s2.s12.cplx + coord_cartesian(ylim = c(-0.5,0.5)),
       height = 14, width = 14
)

ggsave(filename = "TeX/images/chap_sim3/comp_bias34.png",
       plot = p.s2.s34.cplx + coord_cartesian(ylim = c(-.5,.5)),
       height = 14, width = 14
)

ggsave(filename = "TeX/images/chap_sim3/comp_bias56.png",
       plot = p.s2.s56.cplx + coord_cartesian(ylim = c(-.5,.5)),
       height = 14, width = 14
)




p.s2.s34.saft    <- makebiasplot(filter(plot.df, sim == 2, scenario %in% c(3,4), type %in% s2.2saf.types), nc = 2, lnc = 4)

ggsave(filename = "TeX/images/chap_sim3/2saft_bias.png",
       plot = p.s2.s34.saft + coord_cartesian(ylim = c(-0.3,0.1)),
       height = 14, width = 14
)



#########################################
## coverage

p.s2.sall.simpleC <- makecovrgplot(filter(plot.df, sim == 2, scenario <= 6, type %in% s2.simp.types), nc = 3)

ggsave(filename = "TeX/images/chap_sim3/simple_cov.png",
       plot = p.s2.sall.simpleC + coord_cartesian(ylim = c(0,1)),
       height = 20, width = 14
)

p.s2.sall.cplxC <- makecovrgplot(filter(plot.df, sim == 2, scenario <= 6, type %in% c("MIPE-ALL", "RPSFT-TG-LR", "RPSFT-OT-LR-TU-TOL","RPSFT-OT-LR-V-TOL", "2SAFT")), nc = 3)

ggsave(filename = "TeX/images/chap_sim3/comp_cov.png",
       plot = p.s2.sall.cplxC + coord_cartesian(ylim = c(0,1)),
       height = 20, width = 14
)

#########################################
## convergence

s2.dfX <- filter(plot.df, 
              type %in% c("MIPE", "MIPE-ALL", "RPSFT-TG-LR", "RPSFT-OT-LR-TU", "RPSFT-OT-LR-TU-TOL"),
              sim==2)

s2.dfX <- mutate(s2.dfX, 
              type = ifelse(type == "MIPE", "MIPE", 
                            ifelse(type == "MIPE-ALL", "MIPE-ALL", 
                                   ifelse(type == "RPSFT-TG-LR", "RPSFT-TG-LR", 
                                          ifelse(type == "RPSFT-OT-LR-TU", "RPSFT-OT-LR", "RPSFT-OT-LR-TOL")))))

s2.dfX <- mutate(s2.dfX, type = factor(type, levels = c("RPSFT-TG-LR", "RPSFT-OT-LR", "RPSFT-OT-LR-TOL","MIPE","MIPE-ALL"), ordered = TRUE))

p.s2.complexCV <- makeconvplot(s2.dfX, nc = 3)

ggsave(filename = "TeX/images/chap_sim3/comp_conv.png",
       plot = p.s2.complexCV + coord_cartesian(ylim = c(0,1)),
       height = 20, width = 14
)


require(reshape2)

plot2.df <- dcast(filter(plot.df, newsim %in% c(1,2)), newsim+newscenario+p.switch+rho+true.hr~type, value.var = "pbias")
plot2.df <- mutate(plot2.df, lblp = ifelse(newsim == 1, "Study 2",
                                           ifelse(newscenario %in% c(1,2), "Study 3 Scenario 1 and 2",
                                                  ifelse(newscenario %in% c(3,4), "Study 3 Scenario 3 and 4", "Study 3 Scenario 5 and 6"))))

x <- "RPSFT-TG-LR"
y <- "RPSFT-OT-LR-TU"

scplot <- function(x,y){

  rc <- ggplot(plot2.df) + 
    geom_point(aes(x=get(x, plot2.df), y=get(y, plot2.df), color = lblp, shape = lblp), size = 3) +
    xlab(paste("Bias for", x,"(%)")) +
    ylab(paste("Bias for", y,"(%)")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) 

  return(rc)
    
}


mlist <- c("RPSFT-TG-LR","RPSFT-TG-W","RPSFT-OT-LR-TU-TOL","RPSFT-OT-LR-V-TOL")


plist <- NULL
for (xi in 1:(length(mlist))){
  for (yi in 1:(length(mlist))){
    if (xi==yi){
      px <- ggplot() + theme_void() +  annotate("text", x = 0, y = 0, label = mlist[xi])
    } 
    else{
      px <- scplot(mlist[xi], mlist[yi]) + mscTheme()
    }
    if (is.null(plist)){
      plist <- list(px)
    } else{
      plist <- c(plist, list(px))
    }
  }
}

require(GGally)

p1 <- ggmatrix(plist, 
               nrow = length(mlist), 
               ncol = length(mlist)
               )

p1

p1 <- scplot("RPSFT-TG-LR","RPSFT-OT-LR-TU-TOL")
p2 <- scplot("RPSFT-TG-LR","RPSFT-OT-LR-V-TOL")
p3 <- scplot("RPSFT-TG-LR","RPSFT-TG-W")
p4 <- scplot("RPSFT-TG-LR","MIPE-ALL")


p2 <- ggplot(plot2.df) + geom_point(aes(x=`RPSFT-TG-LR`, y=`RPSFT-OT-LR-TU`, color = as.factor(sim)))  
p3 <- ggplot(plot2.df) + geom_point(aes(x=`RPSFT-TG-LR`, y=`RPSFT-OT-LR-V`, color = as.factor(sim)))  





####################### plots by sim

plot.df2 <- filter(mutate(plot.df, basis.scenario = as.numeric(lbl)), !is.na(basis.scenario))

f.bp <- function(df){
  rc <- ggplot(data = df) + geom_bar(aes(x=type, y=pbias, fill = type), stat = "identity") + 
  facet_grid(~basis.scenario) + mscTheme() + scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = 0, size = 1)  +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) +
  ylab("Bias (%)") + xlab("") +
  theme(axis.text.x  = element_blank()) +
  theme(axis.line.x = element_blank()) 
  
  return(rc)
  }

f.pb2 <- function(sp){

p1 <- f.bp(filter(sp, rho == 0)) + ggtitle(expression(paste("Correlation (",rho,")=0"))) 
p2 <- f.bp(filter(sp, rho == 0.2)) + ggtitle(expression(paste("Correlation (",rho,")=0.2")))
p3 <- f.bp(filter(sp, rho == 0.4)) + ggtitle(expression(paste("Correlation (",rho,")=0.4")))
p4 <- f.bp(filter(sp, rho == 0.6)) + ggtitle(expression(paste("Correlation (",rho,")=0.6")))
p5 <- f.bp(filter(sp, rho == 0.8)) + ggtitle(expression(paste("Correlation (",rho,")=0.8")))
p6 <- f.bp(filter(sp, rho == 1.0)) + ggtitle(expression(paste("Correlation (",rho,")=1.0")))  

require(gridExtra)

rc1 <- grid.arrange(p1 +theme(legend.position="none"),p2 +theme(legend.position="none"),p3, ncol = 1)
rc2 <- grid.arrange(p4 +theme(legend.position="none"),p5 +theme(legend.position="none"),p6, ncol = 1)
rc3 <- grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),p3 +theme(legend.position="none"), p4+theme(legend.position="none"),p5+theme(legend.position="none"),p6, ncol = 1)
rc4 <- grid.arrange(p3 +theme(legend.position="none"), p4 +theme(legend.position="none"),p5, ncol = 1)
rc <- list(p1 = rc1, p2 = rc2, p3 = rc3, p4 = rc4)

return(rc)
}

simpl.bar <- f.pb2(filter(plot.df2, type %in%  c("ITT", "PP-CENS", "PP-EX", "TVC", "TVC2")))
rpsft.bar <- f.pb2(filter(plot.df2, type %in%  c("RPSFT-TG-LR", "RPSFT-TG-W", "RPSFT-OT-LR-TU-TOL", "RPSFT-OT-W-TU-TOL", "RPSFT-OT-LR-V-TOL", "RPSFT-OT-W-V-TOL")))
mipe.bar  <- f.pb2(filter(plot.df2, type %in%  c("RPSFT-TG-LR", "RPSFT-OT-LR-TU-TOL", "MIPE", "MIPE-ALL", "MIPE-WEIB")))
saft.bar  <- f.pb2(filter(plot.df2, type %in%  c("2SAFT", "2SAFT-NC", "2SAFT-NULL", "2SAFT-NC-NULL")))

ggsave(filename = "TeX/images/app_allres/simp_bias4.png",
       plot = simpl.bar$p4 ,
       height = 14, width = 20
)

ggsave(filename = "TeX/images/app_allres/rpsft_bias4.png",
       plot = rpsft.bar$p4 ,
       height = 14, width = 20
)

ggsave(filename = "TeX/images/app_allres/saft_bias4.png",
       plot = saft.bar$p4 ,
       height = 14, width = 20
)

ggsave(filename = "TeX/images/app_allres/mipe_bias4.png",
       plot = mipe.bar$p4 ,
       height = 14, width = 20
)

