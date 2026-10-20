####Tfs target use pyscenic

library(psych)

library(pheatmap)
library(stringr)
exp <- readRDS("scenic_select.RDS")
type <- readRDS("scenic_type.RDS")
mer <- merge(type,t(exp),all = TRUE,by = 0)
rownames(mer) <- mer$Row.names
mer <- mer[,-1]
mer$stage_cell <- str_replace(mer$stage_cell,"MSM","MSMC")
unique(mer$stage_cell)
stage_cell <- c('stage1_FIB','stage2_FIB','stage3_FIB','stage4_FIB',"stage1_SMC","stage2_SMC","stage3_SMC","stage4_SMC",
                "stage1_MSMC","stage2_MSMC","stage3_MSMC","stage4_MSMC","stage1_MAC","stage2_MAC","stage3_MAC","stage4_MAC",
                "stage1_EC","stage2_EC","stage3_EC","stage4_EC")
tf <- rownames(exp)
data <- matrix(rnorm(460), nrow = 20, ncol = 23)

df <- as.data.frame(data)

colnames(df) <- tf
rownames(df) <- stage_cell
Pvalue <- df
mean_df <- df
for(i in unique(mer$stage_cell)){
    mer1 <- subset(mer,subset = stage_cell==i)
    for(tfs in tf){
        mer$type<-"others"
        mer[which(mer$stage_cell==i),]$type <- i
        p=t.test(mer[,tfs]~mer$type)$p.value
        means <- mean(mer1[,tfs])
        Pvalue[which(rownames(Pvalue)==i),][,tfs] <- p
        mean_df[which(rownames(mean_df)==i),][,tfs] <- means
    }
}
saveRDS(Pvalue,"4a_p.RDS")
saveRDS(mean_df,"4a_tfs_score.RDS")
#cols；
mycol<-colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200)

Pvalue[Pvalue>=0 & Pvalue < 0.001] <- "***"
Pvalue[Pvalue>=0.001 & Pvalue < 0.01] <- "**"
Pvalue[Pvalue>=0.01 & Pvalue < 0.05] <- "*"
Pvalue[Pvalue>=0.05 & Pvalue <= 1] <- ""
pheatmap(t(mean_df),scale = "none",
         border_color ="white",
         number_color="white",
         cluster_rows=T,angle_col = 45,
         cluster_cols=F,color = mycol,display_numbers = t(Pvalue),
         show_rownames=T,show_colnames=TRUE)