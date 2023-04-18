
# I have three or four functions, I need to make into a data.frame
# 

library(ggplot2)
library(latex2exp)

library(tidyverse)

rate <- function (t, wg, er){
	3.7318E10 * wg * t^(-1.5) / (0.9805806) * exp(-11.605 *er / t)
}
rate_error <- function (t, wg, er, dwg, der){
    sqrt( (rate(t,wg,er)*dwg/wg)^2 + (rate(t,wg,er)*der*(-11.605 / t))^2 )
}
# print(f(0.2,2.60E-15,0.1614))

r1_e<-0.413 # MeV
r1_eerror<-0.0008
r1_wg<-3.27E-08 # MeV

r2_e<-0.3762
r2_eerror<-0.0010
r2_wg<-2.40E-10

r3_e<-0.1614
r3_eerror<-0.0015
r3_wg<-2.60E-15

t <- seq(0.1,0.6, length.out = 1000)
r1 <- rate(t, r1_wg, r1_e)
r1_er <- rate_error(t, r1_wg, r1_e, 0.25*r1_wg, r1_eerror)
r2 <- rate(t, r2_wg, r2_e)
r2_er <- rate_error(t, r2_wg, r2_e, 0.3*r2_wg, r2_eerror)
r3 <- rate(t, r3_wg, r3_e)
r3_er <- rate_error(t, r3_wg, r3_e, 0.3*r3_wg, r3_eerror)
# rt <- r1 + r2 + r3

data <- data.frame(t, r1, r2, r3, r1_e, r2_e, r3_e)
head(data)

df <- pivot_longer(data, cols = c('r1', 'r2', 'r3','r1_e', 'r2_e', 'r3_e'), 
    names_to = c('state','energy'), 
    values_to = "rate",
    names_sep = "_",)
# df <- pivot_longer(data, cols = c('er_0p', 'er_1p', 'er_3p'))

# df <- mutate(df, error_plus = value + rate_error(t, wg_0p, er_0p, 0.3*wg_0p, er_0p_error))
head(df)

# Data preparation
# df <- data %>%
#   select(t, r1, r2, r3, rt) %>%
#   gather(key = "variable", value = "value", -t)
# head(df)

## plotting complete data.frame at once
# ggplot(df, aes(x = t, y = value)) + 
#     geom_line(aes(color = variable, linetype = variable)) + 
#     geom_ribbon(alpha=0.5) + 
#     theme_classic() +   
#     theme(legend.position = c(0.8, 0.3)) +
#     xlab("Temperature (GK)") + 
#     ylab(parse(text = TeX('$N_{A} < \\sigma  \\nu > (cm^{-3} s^{-1} mol^{-1} )$'))) + 
#     labs(color = "",linetype = "") +  
#     scale_y_log10() 

## plotting each line  
# ggplot(data = df, aes(x = t, y = r1, ymin=r1-er1, ymax=r1+er1, linetype=1), color = "black") + geom_line() +
# geom_line(data = df, aes(x = t, y = r2, ymin=r2-er2, ymax=r2+er2, linetype=1), color = "black")  +

#     # geom_line(aes(y = r2, ymin=r2-er2, ymax=r2+er2, linetype=1), color = "black") +
#     # geom_line(aes(y = r3, ymin=r3-er1, ymax=r3+er3, linetype=1), color = "black") +
#     # geom_line(aes(y = rt), color = "black") +
#     geom_ribbon(alpha=0.5) + 
#     theme_classic() +
#     scale_y_continuous(trans='log10')   


## plotting just the functions 
# ggplot(data.frame(t = c(0.1, 0.6)), aes(x = t)) + 
#     geom_function(fun = f, n = 1000, args = list(wg = 2.60E-15, er = 0.1614)) + 
#     geom_function(fun = f, n = 1000, args = list(wg = 2.40E-10, er = 0.3672)) + 
#     geom_function(fun = f, n = 1000, args = list(wg = 1.08E-08, er = 0.4138)) + 
#     theme_classic() +
#     scale_y_continuous(trans='log10')   



# ggsave("rrate.png", width = 8.6, height = 8.6, units = "cm")
# ggsave("ex_energy.eps", width = 8.6, height = 8.6, units = "cm")