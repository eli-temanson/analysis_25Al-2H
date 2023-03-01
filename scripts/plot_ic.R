
# resources: 
#   https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/plotmath.html

library(latex2exp)
library(ggplot2)

df <- read.csv("out.csv")

# 2d histogram with default option
ggplot(df, aes(x=ic_E+ic_dE, y=ic_dE) ) +
  theme_classic() + 
  # theme(panel.background = element_rect(colour = "black", size=1)) +
  geom_bin2d(bins=254, show.legend = FALSE) +
  scale_fill_continuous(type = "viridis") +
  xlab("Ion Chamber E (arb)") + 
  ylab("Ion Chamber dE (arb)") + 
  annotate("text", x=1550, y=530, label=TeX("$^{25}Al^{+13}$")) + 
  annotate("text", x=1100, y=560, label=TeX("$^{24}Mg^{+11}$")) + 
  annotate("text", x=1450, y=400, label=TeX("$^{24}Mg^{+12}$"))
  # theme_bw()

ggsave("ic_ede.eps", width = 8.6, height = 6.45, units = "cm")
ggsave("ic_ede.png", width = 8.6, height = 6.45, units = "cm")
