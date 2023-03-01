
library(ggplot2)
library(ggrepel) 
library(latex2exp)

data = read.csv("C2S.csv")
# data

# ggplot(data, aes(x=C2S_l.0, y=C2S_l.2, color=Exp, label=Label)) +
#   theme_classic() +
#   theme(legend.position = c(0.8,0.94),
#         legend.text=element_text(size=7)) + 
#   labs(x = TeX("$C^{2}S_{l=0}$"), y = TeX("$C^{2}S_{l=2}$"), colour = "") + 
#   geom_point() + 
#   geom_errorbar(aes(xmin=C2S_l.0-error_l.0, xmax=C2S_l.0+error_l.0), width=0.02) + 
#   geom_errorbar(aes(ymin=C2S_l.2-error_l.2, ymax=C2S_l.2+error_l.2), width=0.02) + 
#   geom_text_repel(aes(segment.shape = 0),
#                   size = 3, 
#                   seed = 42,
#                   force_pull = 1, # force of attraction to the point
#                   # nudge_x = 0,
#                   # nudge_y = 0.01,
#                   segment.curvature = -1e-20,
#                   segment.angle = 20,
#                   box.padding = 1.3, # padding around the text label
#                   max.overlaps = Inf,
#                   segment.linetype = 2,
#                   min.segment.length = 0, # draw all lines
#                   color = "black") + 
#                   # point.padding = 0.0, # padding lines around data
#                   # force = 1, # force of repulsion between text
#                   # point.padding = 0.1, # additional padding 
#                   # arrow = arrow(length = unit(0.015,"npc"), 
#                   #               ends = "first", 
#                   #               type = "closed"),
#   coord_cartesian(xlim = c(0,0.8), ylim = c(0, 0.8))


ggplot(data, aes(x=C2S_l.0, y=C2S_l.2, shape=Label, color=Exp)) +
  geom_point(size=2) + 
  scale_shape_manual(values = c(0,1,2,15,16,17)) + 
  theme_classic() +
  theme(legend.position = c(0.76, 0.6),
        legend.text=element_text(size=7),
        legend.spacing.y = unit(0.01, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE)) + 
  labs(x = TeX("$C^{2}S_{l=0}$"), y = TeX("$C^{2}S_{l=2}$"), colour = "", shape="") + 
  geom_errorbar(aes(xmin=C2S_l.0-error_l.0, xmax=C2S_l.0+error_l.0), width=0.02) + 
  geom_errorbar(aes(ymin=C2S_l.2-error_l.2, ymax=C2S_l.2+error_l.2), width=0.02) + 
  coord_cartesian(xlim = c(0.0,0.5), ylim = c(0.0, 0.9))


ggsave("C2S.eps", width = 8.6, height = 8.6, units = "cm")
ggsave("C2S.png", width = 8.6, height = 8.6, units = "cm")





