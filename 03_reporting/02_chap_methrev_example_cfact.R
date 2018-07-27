# Purpose: Example plot showing counterfactuals 
# Functions: 
#   sim.GenStudy
#   mscTheme
# Inputs: 
#   data/sims/pilot.sim.speclists.rData
# Outputs: 
#   TeX/images/chap_methrev/example_cfact.png
# Note:
#   original png created from earlier simulation which is lost

rm(list=ls())
date()

##### Globals and libraries 
library(dplyr)
library(survival)
library(ggplot2)
require(ggthemes)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

load("data/sims/sim.speclists.rData")
load("data/sims/sim.specs.rData")
filter(spec.df, sim == 1, rho == 0.6, true.hr == 0.7) %>%
  summarise(min(true.hr))
df <- do.call(sim.GenStudy, sim2.spec.list[[5001]]) 

rpsft.input.OT <- transmute(df, 
                            event.time = os.t,
                            trt.ind = x.trt, 
                            censor.ind = os.e, 
                            t.on  = ifelse(x.trt, pfs.t, ifelse(x.switch, pmin(pfs.t, t.censor-pfs.t), 0)),
                            t.off = os.t - t.on,
                            cutofftime = t.censor
                            )

this.rpsft.OT <- RPSFT.2pass(rpsft.input.OT)

rp.TU <- RPSFT.cox(rpsft.input = rpsft.input.OT, 
                   rpsft.output = this.rpsft.OT, 
                   Grho = 0,
                   use.latent.only = FALSE)

rp.VV <- RPSFT.cox(rpsft.input = rpsft.input.OT, 
                   rpsft.output = this.rpsft.OT, 
                   Grho = 0,
                   use.latent.only = TRUE)


tim.df <- rbind(
  data.frame(typ = "T", trt = df$x.trt, tim = rp.TU$cfact.time, cns = rp.TU$cfact.censor.ind),
  data.frame(typ = "V", trt = df$x.trt, tim = rp.VV$cfact.time, cns = rp.VV$cfact.censor.ind)
)

rpsft.output <- this.rpsft.OT
rpsft.input  <- rpsft.input.OT


sf.E.T <- survfit(Surv(tim,cns)~1, data = filter(tim.df, typ == "T", trt == 1))
sf.C.T <- survfit(Surv(tim,cns)~1, data = filter(tim.df, typ == "T", trt == 0))
sf.E.V <- survfit(Surv(tim,cns)~1, data = filter(tim.df, typ == "V", trt == 1))
sf.C.O <- survfit(Surv(os.t,os.e)~1, data = filter(df, x.trt == 0))
sf.E.O <- survfit(Surv(os.t,os.e)~1, data = filter(df, x.trt == 1))

kmdf <- function(sf){
  data.frame(t = c(0,sf$time),
             s = c(1,sf$surv),
             c = c(0, as.numeric(sf$n.censor>=1)),
             stringsAsFactors = TRUE)
}

d1 <- mutate(kmdf(sf.E.O), type = "Observed",       trt = "Experimental", lbl = "Observed")
d2 <- mutate(kmdf(sf.C.O), type = "Observed",       trt = "Control",      lbl = "Observed")
d3 <- mutate(kmdf(sf.E.T), type = "Observed",       trt = "Experimental", lbl = "Observed vs Latent")
d4 <- mutate(kmdf(sf.C.T), type = "Counterfactual", trt = "Control",      lbl = "Observed vs Latent")
d5 <- mutate(kmdf(sf.E.V), type = "Counterfactual", trt = "Experimental", lbl = "Counterfactual vs Latent")
d6 <- mutate(kmdf(sf.C.T), type = "Counterfactual", trt = "Control",      lbl = "Counterfactual vs Latent")


plot.df <- rbind(d1,d2,d3,d4,d5,d6)
plot.df$type <- factor(plot.df$type, levels = c("Observed", "Counterfactual"), ordered = TRUE)
plot.df$lbl <- factor(plot.df$lbl, levels = c("Observed", "Observed vs Latent",  "Counterfactual vs Latent" ), ordered = TRUE)

p1 <- ggplot(data = plot.df) + 
  mscTheme(base_size = 18) +
  geom_step(aes(x = t, y = s, linetype = type, color = trt), size = 1, direction = "hv") + 
  geom_point(aes(x= t, y= s,  shape = "Censored"), data = filter(plot.df,c == 1), size = 2) +
  xlab("Overall Survival Time (years)")  +
  ylab("Survival") +
  coord_cartesian(ylim = c(0,1), xlim = c(0,3)) +
  scale_color_manual(values = c("Red","Blue"), name = " ", breaks=c("Experimental","Control")) +
  #scale_linetype_manual(name = " ", values = c(2,1)) +
  scale_shape(name = " ", solid = FALSE) +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  ) +
  geom_vline(xintercept = 0, size = 1) + 
  geom_hline(yintercept = 0, size = 1) +
  facet_wrap(~lbl, scales = "free")  + 
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position="bottom", legend.direction="horizontal", legend.box = "horizontal") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + theme(axis.ticks.length = unit(0,"cm")) 
p1

ggsave(plot = p1, filename = "TeX/images/chap_methrev/example_cfact.png", width = 12, height = 7)





