
# to run in the terminal: Rscript script_name.R 

library(ggplot2)
library(latex2exp)

setwd("~/25Al+d_2014/R-scripts")

g1 <- c(4.50478e+01,5.92260e+00, 6.57838e-02)
g2 <- c(3.73793e+01,6.32535e+00,8.73756e-02)
g3 <- c(4.29107e+01,6.69900e+00,1.44333e-01)
g4 <- c(1.67444e+01,7.34902e+00,2.41747e-01)

fit.calc <- function (x) {
  g1[1]*exp(-0.5*((x-g1[2])/g1[3])^2) + g2[1]*exp(-0.5*((x-g2[2])/g2[3])^2) + g3[1]*exp(-0.5*((x-g3[2])/g3[3])^2) + g4[1]*exp(-0.5*((x-g4[2])/g4[3])^2)
}

# load data into a data.frame
df <- read.csv("ex_energy.csv")
sim <- read.csv("~/25Al+d_2014/dev/InvKinCalc/data/efficieny_data.txt",header = FALSE, sep=" ")

ggplot() + # plot histogram
  theme_classic() + 
  geom_histogram(data = df, aes(x=ex_energy), bins=120, fill="gray") +
  xlab("Excitation Energy (MeV)") + 
  ylab("Counts") +
  geom_vline(xintercept = 5.514, linetype="dashed", color = "black") + 

  # plot a top and bottom axis
  scale_x_continuous(name = TeX("$^{26}Si$ Excitation Energy (MeV)"), limits = c(4.514,9.514),
                     sec.axis = sec_axis(trans = ~.-5.514, breaks = c(0,0.5,1,1.5,2,2.5), name = TeX("Proton Resonance $E_{cm}$ (MeV)" ))) + 
  
  # plot the simulation line on top
  geom_line(data = sim, aes(x=V1,y=V2*50), color='blue') + 
  scale_y_continuous(name = TeX("Counts"), sec.axis = sec_axis(~./50, breaks = c(0,0.5,1.0), name = TeX('Detection Efficiency'))) + 
  theme(axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"), 
        axis.line.y.right = element_line(color = "blue"), 
        axis.ticks.y.right = element_line(color = "blue")) + 
  coord_cartesian(xlim=c(5.1,8.2)) + 

  # plot the fit function on top
  stat_function(fun=fit.calc, n = 1000)

# save the figure at the necessary size
ggsave("ex_energy.png", width = 8.6, height = 8.6, units = "cm")
ggsave("ex_energy.eps", width = 8.6, height = 8.6, units = "cm")

  
  