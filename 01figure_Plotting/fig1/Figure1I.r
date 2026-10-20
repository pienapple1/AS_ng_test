###figure1 I
###R4

library(tidyverse)
library(ggplot2)
data <- readRDS("1g.RDS")
cols <- c('W0'= '#FFE3F0','W2'='#EFA5B9','W6'='#E4D7FF','W10'='#ABA0FB','W20'='#736DC2','W20_N'="#A1CFFA")
ggplot(data, aes(x = thickness_mm, y = CV,color=week,fill=week)) +geom_errorbarh(aes(xmin = Min_thickness, xmax = Max_thickness,
                     color = week), height = 0.01, linewidth = 0.2) +geom_errorbar(aes(ymin = Min_cv, ymax = Max_cv,color = week), 
                width = 0.01, linewidth = 0.2) +
  geom_boxplot(size=4)+ylim(0.2,0.6)+xlim(0.08,0.25)+
  scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
  theme_classic()+theme(axis.text.x = element_text(size=15,color = "black"),
         axis.text.y = element_text(size=15,colour = "black"))+coord_flip()