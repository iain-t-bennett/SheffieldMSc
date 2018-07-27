# Purpose: Example plot showing counterfactuals 
# Functions: 
#   mscTheme
# Inputs: 
#   data/raw/sim1.res.rData
# Outputs: 
#   TeX/images/chap_simdesign/p_switch.png

rm(list=ls())
date()

##### Globals and libraries 
library(dplyr)
library(survival)
library(ggplot2)


##### functions
sapply(list.files(pattern="[.]R$", path="functions", full.names=TRUE), source)

######### load the results

load(file = "data/raw/sim1.res.rData")

######### summarise by scenario

my_theme <-  theme_classic(base_size = 24)+
  theme(
    axis.title = element_text( colour = "black",face = "bold"),
    axis.text.x = element_text( colour = "black",face = "bold"),
    axis.text.y = element_text( colour = "black",face = "bold"),
    legend.text = element_text( colour = "black",face = "bold"),
    
    strip.text.x = element_text( colour = "black",face = "bold"),
    strip.background = element_blank(),
    legend.title =  element_text( colour = "black",face = "bold"),
    axis.ticks.length = unit(0,"cm"),
    plot.margin = unit(c(0,0,0,0), "cm"),
    
) 


sim1.res.df <- mutate(sim1.res.df, cov = as.numeric(true.hr > cil & true.hr < ciu), sig = as.numeric(pval <= 0.05), trend = as.numeric(hr <1))

res.table <- summarise(group_by(sim1.res.df, sim, scenario, p.switch, prop.pd.int, rho, true.hr),hr = mean(hr), cov = mean(cov), trnd = mean(trend), sig = mean(sig), a.switch = mean(a.switch))
res.table <- mutate(res.table, bias = hr-true.hr, pbias = bias/true.hr)



# bias


res.table <- mutate(res.table, 
                    lbl = paste("True HR: ", true.hr,"\nTarget Switch: ",p.switch*100,"%", sep =""),
                    proppfs =  paste(prop.pd.int *100,"%", sep ="")
                    )




p.pswitch <- ggplot(res.table, aes(x=p.switch, y=a.switch,color=proppfs, shape = proppfs)) +
  mscTheme() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size  = unit(5, "lines"),
    axis.line.y = element_line(size = 1, colour = "black"),
    axis.line.x = element_line(size = 1, colour = "black")
  ) + 
  geom_jitter(width = 0.07, size = 4) +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.2,0.4,0.6)) +
  scale_x_continuous(labels = scales::percent, breaks = c(0,0.2,0.4,0.6)) +
  ylab("Actual Switch (%)") + xlab("Target Switch (%)") +
  guides(color = guide_legend(title = "Interim PFS events"),
         shape = guide_legend(title = "Interim PFS events"))

ggsave(filename = "TeX/images/chap_simdesign/p_switch.png",
       plot = p.pswitch,
       height = 8, width = 10
)

  



