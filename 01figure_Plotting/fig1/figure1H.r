###figure1f_1
data <- readRDS("1H.RDS")
summ <-summarySE(data, measurevar="percent", groupvars=c("week","celltype"))
order <- c("W20_N","W0","W2","W6","W10","W20")
summ$week <- factor(summ$week,levels = order)
# Standard error of the mean
mycolor <- c("Mp"='#df65b0',"MSMC"="#A1CFFA",
                                 "SMC"='#1f78b4',"EC"='#238b45',"Fibro"="#8b4a4b",
                                 "Tcell"='#ce1256')
ggplot(summ, aes(x=week, y=percent,color=celltype,group=celltype))+
geom_errorbar(aes(ymin=percent-se, ymax=percent+se), width=.1) +
scale_color_manual(values = mycolor)+labs(x="",y="")+
geom_line(cex=0.8) +geom_point()+theme_classic()+
theme(axis.text.x = element_text(size=15,color = "black"),axis.text.y = element_text(size=15,color = "black"))