# Purpose: Creates plots for body and appendix of simulation results
# Functions: 
#   mscTheme
# Inputs: 
#   data/derived/res_all_sim234.rData
# Outputs: 
#   TeX/images/chap_sim3/truehr.png
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
## hr for scenarios in study 2

p.hr <- ggplot(data = filter(plot.df, sim == 2, scenario <=4)) + 
  mscTheme() + 
  geom_point(aes(x=rho, y=true.hr, color = lbl1, shape = lbl1), size = 3) + 
  geom_line(aes(x=rho, y=true.hr, color = lbl1, linetype = lbl1), size = 1) +
  ylab("True Hazard Ratio") +
  xlab(expression(paste("Correlation (", rho, ")"))) +
  geom_hline(yintercept = c(0.7,0.9), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(0.5,1)) +
  theme(legend.key.size  = unit(5, "lines")) 

ggsave(filename = "TeX/images/chap_sim3/truehr.png",
       plot = p.hr,
       height = 7, width = 9
)

