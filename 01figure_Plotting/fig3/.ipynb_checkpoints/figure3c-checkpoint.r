####3B gene expr by layer
###Author Jinpei Lin
###Times 2024.09
###R4
library(ggplot2)
mer <- readRDS("3b.RDS")
order <- c('W20_N','W0','W2','W6','W10','W20')
week <- order
mer$week <- factor(mer$week,levels = order)
cols <- c('W0'= '#FFE3F0','W2'='#EFA5B9','W6'='#E4D7FF','W10'='#ABA0FB','W20'='#736DC2','W20_N'="#A1CFFA")
for(gene in c("Pi16","Col15a1","Col1a1","Fn1")){
    s=0
    for(i in week){
        s=s+1
        w <- subset(mer,subset = week==i)
        if(s==1){
            m <- w
        }
        else{
            m <- rbind(m,w)
        }
    }
    m$week <- factor(m$week,level=week)
    p <- ggplot(m,aes(x=layer1,y = get(gene),fill=week,color=week))+theme_classic ()+
    labs(subtitle = gene)+geom_smooth(method = "gam",span=1,se = T,level=0.90,alpha=0.15)+scale_fill_manual(values =cols)+scale_color_manual(values = cols)+
    labs(x="digitization_layer",y="expression")+theme_bw()+theme(panel.grid = element_blank())+
    theme_classic()+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))
    print(p)
}