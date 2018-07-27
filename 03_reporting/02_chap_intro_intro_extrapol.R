# Purpose: Example KM plot showing extrapolation
# Functions: 
#   sim.GenStudy
#   mscTheme
# Inputs: 
#   data/sims/pilot.sim.speclists.rData
# Outputs: 
#   TeX/images/chap_intro/intro_extrapol.png


rm(list=ls())
date()

##### Globals and libraries 
library(dplyr)
library(ggplot2)
library(survival)

##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

#### load an example dataset
load(file = "data/sims/pilot.sim.speclists.rData")
example.df <- do.call(sim.GenStudy, sim.pilot.spec.list[[4]]) %>%
  filter(x.trt == 0)

surv.control <- Surv(example.df$os.o.t, example.df$os.o.e)

fit1 <- survreg(surv.control~1)
intercept <- fit1$coefficients[1]
scale     <- fit1$scale

w.lambda <- exp(-intercept / scale)
w.gamma  <- 1 / fit1$scale

sf <- survfit(surv.control~1)

km.df <- data.frame(t = c(0,sf$time),
                    s = c(1,sf$surv),
                    c = c(0, as.numeric(sf$n.censor>=1)),
                    armcd  = "Control",
                    type   = "Kaplan Meier",
                    emean = NA,
                    elife = NA,
                    w.lambda= NA,
                    w.gamma = NA,
                    stringsAsFactors = TRUE)

exhr0 <- data.frame(t = seq(from = 0, to = 8, by = 0.01))
exhr0 <- mutate(exhr0, 
                    w.lambda = exp(-intercept / scale),
                    w.gamma  = 1 /scale,
                    c  = 0,
                    armcd = "Control",
                    type = "Weibull Extrapolation"
                    )

exhr1 <- transmute(exhr0, t, c, type, w.gamma,
                    w.lambda = w.lambda*0.8,
                    armcd = "Experimental (HR = 0.8)"
                    )

exhr2 <- transmute(exhr0, t, c, type, w.gamma,
                    w.lambda = w.lambda*0.7,
                    armcd = "Experimental (HR = 0.7)"
)


ex.df <- rbind(exhr0, exhr2, exhr1)

ex.df <- mutate(ex.df, 
                s= exp(-w.lambda*t^w.gamma),
                emean = gamma(1+1/w.gamma) / (w.lambda^(1/w.gamma)),
                elife = s*0.01
                )

plot.df <- rbind(km.df, ex.df)


plt2.df <- select(filter(plot.df, type != "Kaplan Meier"), armcd, t, s, emean)
plt2.df <- merge(plt2.df, transmute(filter(plt2.df, armcd == "Control"), t, s0 = s), by = "t")
plt2.df <- merge(plt2.df, transmute(filter(plt2.df, armcd == "Experimental (HR = 0.8)"), t, s1 = s), by = "t")
plt2.df <- merge(plt2.df, transmute(filter(plt2.df, armcd == "Experimental (HR = 0.7)"), t, s2 = s), by = "t")

plt2.df$ym <- 0
plt2.df$ym[plt2.df$armcd=="Experimental (HR = 0.8)"] <- plt2.df$s0[plt2.df$armcd=="Experimental (HR = 0.8)"]
plt2.df$ym[plt2.df$armcd=="Experimental (HR = 0.7)"] <- plt2.df$s0[plt2.df$armcd=="Experimental (HR = 0.7)"]

plt2.df <- mutate(plt2.df, armcd2 = paste(armcd, " - mean life = ", round(emean,2)))

summarise(group_by(plt2.df, armcd), unique(emean))


p1 <- ggplot(data = plot.df) + 
  geom_ribbon(data = plt2.df, aes(x = t, ymin = ym, ymax = s, fill = armcd), alpha = 0.2) +
  geom_step(aes(x = t, y = s, linetype = as.factor(type), color = armcd), direction = "hv") + 
  geom_point(aes(x= t, y= s,  shape = "Censored"), data = filter(plot.df,c == 1)) +
  xlab("Overall Survival Time (years)")  +
  ylab("Survival") +
  coord_cartesian(ylim = c(0,1), xlim = c(0,6)) +
  scale_fill_manual(values = c("Red","Blue", "Green"), name = " ", breaks=c("Experimental (HR = 0.7)","Experimental (HR = 0.8)","Control")) +
  scale_color_manual(values = c("Red","Blue", "Green"), name = " ", breaks=c("Experimental (HR = 0.7)","Experimental (HR = 0.8)","Control")) +
  scale_linetype_manual(name = " ", values = c(1,2)) + 
  scale_shape(name = " ", solid = FALSE) +
  
  theme(legend.position="right")  + 
  annotate("text", x= 1.2, y=0.15, label ="Control\n Extrapolated Mean \n  2.77 years") +
  
  annotate("segment", x = 3.2, xend = 3.5, y = 0.36, yend = 0.50, colour = "black") +
  annotate("text", x= 3.5, y=0.60, label ="Additional vs Control\n Extrapolated Mean \n  0.59 years") +
  
  annotate("segment", x = 4.6, xend = 5, y = 0.29, yend = 0.50, colour = "black") +
  annotate("text", x= 5, y=0.60, label ="Additional vs Control\n Extrapolated Mean \n  1.00 years") 

p1 <- p1 + 
  mscTheme(base_size = 12) + 
  scale_y_continuous(labels = scales::percent) +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  ) +
  geom_vline(xintercept = 0, size = 1) + 
  geom_hline(yintercept = 0, size = 1) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + theme(axis.ticks.length = unit(0,"cm")) 
  

ggsave(plot = p1, filename = "TeX/images/chap_intro/intro_extrapol.png", width = 10, height = 6)

