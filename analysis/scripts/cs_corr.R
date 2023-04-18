library(ggplot2)

x1 <- c(0.05,0.2)
y1 <- c(2.84E+00,1.13E+01)

x2 <- c(0.05,0.14,0.2)
y2 <- c(5.41E+00,8.6512E+00,1.08E+01)

ggplot() + 
  geom_line(aes(x=x1,y=y1)) + 
  geom_line(aes(x=x2,y=y2)) + 
  xlab('C2S') + 
  ylab('Cross Section') + 
  geom_hline(yintercept=8, linetype="dashed", color = "red") + 
  theme_classic()

ggsave(file='./cs_corr.eps', width=8.6, height=6.4,units = 'cm')
