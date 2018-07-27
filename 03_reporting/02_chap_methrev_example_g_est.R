# Purpose: Example plot showing counterfactuals 
# Functions: 
#   sim.GenStudy
#   mscTheme
#   RPSFT
# Inputs: 
#   data/sims/pilot.sim.speclists.rData
# Outputs: 
#   TeX/images/chap_methrev/example_g_est.png
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

load("data/sims/pilot.sim.speclists.rData")

df <- do.call(sim.GenStudy, sim.pilot.spec.list[[20000]])

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

attach(rpsft.TG.LR)
df1 <- data.frame(psi = psi.tried, z = z, pass = paste("Pass:", pass), test = "log-rank")
psi.sel <- psi.chosen
detach(rpsft.TG.LR)

plot.df <- df1

#### make a plot

plot.df$pass[round(plot.df$psi,2) == -0.4] <- "Pass: 1"
plot.df$pass[round(plot.df$psi,2) == -0.3] <- "Pass: 1"
plot.df$pass[round(plot.df$psi,2) == -0.2] <- "Pass: 1"
plot.df$pass[round(plot.df$psi,2) == -0] <- "Pass: 1"
plot.df$pass[round(plot.df$psi,2) == 0.1] <- "Pass: 1"


p1 <- ggplot(data = plot.df) +
  mscTheme() +
  theme(
    legend.title = element_text( colour = "black",face = "bold")
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1 ) +
  geom_hline(yintercept = 0, color = "black", size = 1 ) +
  coord_cartesian(xlim = c(-0.75, 0.25), ylim = c(-4,4)) +
  geom_point(aes(x = psi, y = z, shape = pass, color = pass), size = 4) + 
  xlab(expression(psi)) + 
  ylab(expression(paste("log-rank Z(",psi,")"))) +
  theme(legend.position = "bottom") +
  scale_color_discrete(name = "Search step (decreasing grid size with second pass)") +
  scale_shape_discrete(name = "Search step (decreasing grid size with second pass)") +
  geom_vline(xintercept = psi.sel, color = "black" , linetype = 2, size = 1) +
  annotate("text", x=psi.sel + 0.1, y  = 1, label = paste("Estimate:", psi.sel), size = 8)

p1

ggsave(plot = p1, filename = "TeX/images/chap_methrev/example_g_est.png", width = 14, height = 10)




