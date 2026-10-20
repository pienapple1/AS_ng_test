#####figure2C GC LC plaque score by addmodulescore
###R4

library(ggplot2)
GC <-readRDS("GC.RDS")
LC <-readRDS("LC.RDS")
cols <- c('W0'= '#FFE3F0','W2'='#EFA5B9','W6'='#E4D7FF','W10'='#ABA0FB','W20'='#736DC2','W20_N'="#A1CFFA")
ggplot(GC,aes(x=column1,y = plaque_score,fill=week,color=week))+theme_classic ()+
labs(subtitle = "GC")+geom_smooth(method = "gam",span=1,se = T,level=0.90,alpha=0.15,sep = "~~~~")+scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
labs(x="digitization_Column",y="plaue_score")
ggplot(LC,aes(x=column1,y = plaque_score,fill=week,color=week))+theme_classic ()+
labs(subtitle = "LC")+geom_smooth(method = "gam",span=1,se = T,level=0.90,alpha=0.15,sep = "~~~~")+scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
labs(x="digitization_Column",y="plaue_score")