
###R4
####figure3 B plauqe score and gene expr

library(ggplot2)
library(reshape2)

sp1 <- readRDS("satge_bulk_celltype_filter.RDS")#######/00data_Processing/1.8pseudobulk_of_celltype_area.ipynb
exp <- t(as.data.frame(sp1@assays$RNA@data))
mer <- merge(sp1@meta.data,exp,all = TRUE,by = 0)
rownames(mer) <- mer$Row.names
m <- mer[,-1]

###SMC
m1 <- subset(m,subset = cell=="SMC")
gene<-c("Myl9","Acta2","Myh11","Mylk","Tagln")
m1 <- subset(m1,select = c('gsva_p_all','ssgsea_scale','ssgsea',gene))
m2 <- melt(m1,id.vars = c('gsva_p_all','ssgsea_scale','ssgsea'))
colnames(m2)<- c('gsva_p_all','ssgsea_scale','ssgsea',"gene","exp")
p1 <- ggplot(m2,aes(x=gsva_p_all,y = exp,color=gene))+theme_classic ()+scale_color_manual(values = colorRampPalette(c("#97D1A0","#ABDAEC","#6A9ACE","#F55665FF"))(5))+labs(subtitle = "SMC")+geom_smooth(method = "loess",span=1,se = TRUE,level=0.9,fill='#f1eef6')+
labs(x="plaue_score",y="expr")+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))
p1
###MSMC
m1 <- subset(m,subset = cell=="Modulated_SMC")
gene<-c("Sox9","Ly6a","Myh11","Vcam1","Ibsp")
m1 <- subset(m1,select = c('gsva_p_all','ssgsea_scale','ssgsea',gene))
m2 <- melt(m1,id.vars = c('gsva_p_all','ssgsea_scale','ssgsea'))
colnames(smc2)<- c('gsva_p_all','ssgsea_scale','ssgsea',"gene","exp")
p2 <- ggplot(smc2,aes(x=gsva_p_all,y = exp,color=gene))+theme_classic ()+scale_color_manual(values = colorRampPalette(c("#97D1A0","#ABDAEC","#6A9ACE","#F55665FF"))(5))+labs(subtitle = "MSMC")+geom_smooth(method = "loess",span=1,se = TRUE,level=0.9,fill='#f1eef6')+
labs(x="plaue_score",y="expr")+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))
p2
####FIbs
m1 <- subset(m,subset = cell=="Fibroblast")
gene<-c("Pi16","Col1a1","Col4a1","Dcn","Col15a1")
m1 <- subset(m1,select = c('gsva_p_all','ssgsea_scale','ssgsea',gene))
m2 <- melt(m1,id.vars = c('gsva_p_all','ssgsea_scale','ssgsea'))
colnames(m2)<- c('gsva_p_all','ssgsea_scale','ssgsea',"gene","exp")
p3 <- ggplot(m2,aes(x=gsva_p_all,y = exp,color=gene))+theme_classic ()+scale_color_manual(values = colorRampPalette(c("#97D1A0","#ABDAEC","#6A9ACE","#F55665FF"))(5))+labs(subtitle = "Fibro")+geom_smooth(method = "loess",span=1,se = TRUE,level=0.9,fill='#f1eef6')+
labs(x="plaue_score",y="expr")+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))
p3
####EC
m1 <- subset(m,subset = cell=="EC")
gene<-c("Ptprb","Vwf","Fn1","Ece1","Mmp12")
m1 <- subset(m1,select = c('gsva_p_all','ssgsea_scale','ssgsea',gene))
m2 <- melt(m1,id.vars = c('gsva_p_all','ssgsea_scale','ssgsea'))
colnames(m2)<- c('gsva_p_all','ssgsea_scale','ssgsea',"gene","exp")
p4 <- ggplot(m2,aes(x=gsva_p_all,y = exp,color=gene))+theme_classic ()+scale_color_manual(values = colorRampPalette(c("#97D1A0","#ABDAEC","#6A9ACE","#F55665FF"))(5))+labs(subtitle = "EC")+geom_smooth(method = "loess",span=1,se = TRUE,level=0.9,fill='#f1eef6')+
labs(x="plaue_score",y="expr")+theme(axis.text = element_text (size = 20))+#调整坐标轴字体大小
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+#调整xlab字体大小
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+#调整ylab字体大小
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
    theme(legend.title = element_text(size=20),  # 调整图例标题的大小
          legend.text = element_text(size=15))
p4