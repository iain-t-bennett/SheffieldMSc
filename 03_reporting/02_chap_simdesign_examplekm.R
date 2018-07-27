# Purpose: Example KM plots of different simulations
# Functions: 
#   sim.GenStudy
#   mscTheme
# Inputs: 
#   n/a
# Outputs: 
#   TeX/images/chap_simdesign/examplekm.png


require(survival)
require(ggplot2)
require(dplyr)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)


fkm <- function(parms){
  
  df <- do.call(sim.GenStudy, parms)

  
  f <- function(sf, ep, armcd){
    df <- data_frame(t = c(0,sf$time),
                     s = c(1,sf$surv),
                     c = c(0, as.numeric(sf$n.censor>=1)),
                     ep = ep,
                     armcd = armcd)
    return(df)
  }

  os.df.0 <- f(survfit(Surv(os.t,os.e)~1, data = filter(df, x.trt == 0)), ep = "OS", armcd = "Control")
  os.df.1 <- f(survfit(Surv(os.t,os.e)~1, data = filter(df, x.trt == 1)), ep = "OS", armcd = "Experimental")
  pfs.df.0 <- f(survfit(Surv(pfs.t,pfs.e)~1, data = filter(df, x.trt == 0)), ep = "PFS", armcd = "Control")
  pfs.df.1 <- f(survfit(Surv(pfs.t,pfs.e)~1, data = filter(df, x.trt == 1)), ep = "PFS", armcd = "Experimental")

  km.df <- rbind(os.df.0, os.df.1, pfs.df.0, pfs.df.1)

  return(km.df)

}


parms1 <- list(seed        = 122,
              simid       = 1,
              p.switch    = 0,
              prop.pd.int = 0,
              rho         = 0.6,
              beta_1a     = log(0.7),
              beta_1b     = log(0.7),
              beta_2a     = log(0.7),
              beta_2b     = log(0.7)
)

parms2 <- parms1
parms2$p.switch = 0.6

parms3 <- parms1
parms3$seed = 12
parms3$beta_1a = parms3$beta_2a = log(0.01)
parms3$beta_1b = parms3$beta_2b = log(1)
parms4 <- parms3
parms4$p.switch = 0.6

k1 <- mutate(fkm(parms1), sc = 1)
k2 <- mutate(fkm(parms2), sc = 2)
k3 <- mutate(fkm(parms3), sc = 3)
k4 <- mutate(fkm(parms4), sc = 4)

km.df <- rbind(k1,k2,k3,k4)

km.df$scenario <- factor(km.df$sc, 
                         levels =1:4, 
                         labels = c("Example Scenario 1 (No switch)", 
                                    "Example Scenario 1 (60% switch)", 
                                    "Example Scenario 2 (No switch)", 
                                    "Example Scenario 2 (60% switch)"), 
                         ordered = TRUE)

p1 <- ggplot(data = filter(km.df)) + 
  mscTheme() +
  geom_step(aes(x = t, y = s, linetype = as.factor(ep), color = armcd), direction = "hv", size = 2) +
  geom_point(aes(x= t, y= s,  shape = "Censored"), data = filter(km.df,c == 1), size = 3) +
  xlab("Time (years)")  +
  ylab("Survival") +
  geom_hline(yintercept = 0, size = 1) +
  geom_vline(xintercept = 0, size = 1) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(axis.ticks.length = unit(0,"cm")) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~scenario, scales = "free") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size  = unit(5, "lines")
  ) 

p1

ggsave(filename = "TeX/images/chap_simdesign/examplekm.png",
       plot = p1,
       height = 20, width = 14
)






