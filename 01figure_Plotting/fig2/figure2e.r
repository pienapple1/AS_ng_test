####2E
###Author Jinpei Lin
###Times 2024.09
###R4
####Left
library(ggplot2)
m <- readRDS("2e.RDS")
cols <- c("Macrophage"='#df65b0',"Modulated_SMC"="#A1CFFA",
           "SMC"='#1f78b4')
ggplot(m,aes(x=gsva,y = percent,fill=celltype,color=celltype))+theme_classic ()+
labs(subtitle = "")+geom_smooth(method = "loess",span=1,se = T,level=0.95,alpha=0.15)+scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
labs(x="Plaque_score",y="Proportion")+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))
####Rgiht
library(ggtern)
library(RColorBrewer)
library(grid)
library(scales)


ggtern(data=m,aes(x=SMC,y=Macrophage,z=Modulated_SMC))+geom_mask()+
  geom_point(aes(size=gsva_p_all,color=stage_ud),alpha=0.8)+
  scale_colour_manual(values = mycolor)+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  labs(title = "Ternary stage cell")+theme_rgbw()+theme_custom(col.T = "#df65b0",col.L ="#1f78b4" ,col.R = "#A1CFFA",
  tern.panel.background = "white")