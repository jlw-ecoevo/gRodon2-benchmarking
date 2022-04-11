## JLW 2022 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

# Functions --------------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Data --------------------------------------------------------------------

setwd("~/gRodon2-benchmarking/Data/")
load("zou_synthmetaF10_data.rda")
xmF <- xm %>% 
  mutate(Genome=gsub("_.*","",Genome))
load("zou_synthmeta_data.rda")
xm <- xm %>% 
  mutate(Genome=gsub("_.*","",Genome))
xm <- merge.easy(xm,xmF,key="Genome") %>%
  subset(nHE.y>=100 & nHE.x>=100)

# Caluclate errors -------------------------------------------------------------

xm$mix.d.gr <- xm$mix.d.mt.i
xm$mix.d.gr[xm$mix.Consistency<0.6] <- xm$mix.d[xm$mix.Consistency<0.6]
xm$mix.d.gr.depth <- xm$mix.d.mt.depth.i
xm$mix.d.gr.depth[xm$mix.Consistency.depth<0.6] <- xm$mix.d.depth[xm$mix.Consistency.depth<0.6]

xm$d.gr <- xm$d.mt.i
xm$d.gr[xm$Consistency<0.6] <- xm$d[xm$Consistency<0.6]
xm$d.gr.depth <- xm$d.mt.depth.i
xm$d.gr.depth[xm$Consistency.depth<0.6] <- xm$d.depth[xm$Consistency.depth<0.6]

xm$frag.d.gr <- xm$frag.d.mt.i
xm$frag.d.gr[xm$Consistency.y<0.6] <- xm$frag.d[xm$Consistency.y<0.6]

xm$dgr_avg <- xm$dmti_avg
xm$dgr_avg[xm$mix.Consistency<0.6] <- xm$dm_avg[xm$mix.Consistency<0.6]

xm <- xm %>%
  mutate(dgrd_err=(log10(d.gr.depth)-log10(df_avg))^2,
         fdgr_err=(log10(frag.d.gr)-log10(df_avg))^2,
         mdgr_err=(log10(mix.d.gr.depth)-log10(df_avg))^2)

# Plot -------------------------------------------------------------------------

global_size <- 18

p1 <- ggplot(data=xm,aes(x=df_avg,y=d.gr.depth,fill="Synthetic Metagenome (Assembled)")) +
  geom_point(pch=21,size=3) + 
  geom_point(data=xm,aes(x=df_avg,y=frag.d.gr,fill="Synthetic Metagenome (Reads)"),
             pch=21,size=3) +
  geom_point(data=xm,aes(x=df_avg,y=mix.d.gr.depth,fill="Genome Mixture"),
             pch=21,size=3) +
  geom_smooth(color="black",alpha=.5) + 
  geom_smooth(data=xm,aes(x=df_avg,y=frag.d.gr,fill="Synthetic Metagenome (Reads)"),
              color="black",alpha=.5) +
  geom_smooth(data=xm,aes(x=df_avg,y=mix.d.gr.depth,fill="Genome Mixture"),
              color="black",alpha=.5) +
  geom_abline(slope=1,intercept=0,lty=2) + 
  theme_pubclean(base_size = global_size) +
  scale_x_log10() +
  scale_y_log10() +
  labs(fill="") +
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a")) +
  xlab("Avg. gRodon Full Mode Prediction of Genomes in Mixture") +
  ylab("Prediction on Community")

xm$conclass <- "> 0.6"
xm$conclass[xm$Consistency<0.6] <- "< 0.6"
xms <-  xm %>%
  subset(select=c(conclass,
                  mdgr_err,
                  fdgr_err,
                  dgrd_err))
names(xms) <- c("conclass","Genome Mixture","Synthetic Metagenome (Reads)","Synthetic Metagenome (Assembled)")
xmm <- xms %>%
  melt(id.vars="conclass")
p2 <- ggplot(xmm,aes(x=reorder(variable,value),y=value,fill=variable)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 60, hjust=1))  +
  scale_fill_manual(values=c("#a6cee3","#b2df8a","#1f78b4"))


setwd("~/gRodon2-benchmarking/Figures/")
png(file="Compare_zou_synth_F10.png",width=1000,height=1000)
ggarrange(p1,p2,ncol=2,widths=c(3,1),labels=c("(a)","(b)"))
dev.off()

err.aov <- aov(value~variable,data=xmm)
summary(err.aov)
TukeyHSD(err.aov)

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Zou_frag_CUB_F10.png",width=500,height=500)
ggplot(xm,aes(x=CUBHE.depth,y=CUBHE.y)) + 
  geom_point() + 
  geom_abline(slope=1,intercept=0,lty=2) +
  theme_pubclean() +
  ylab("CUB Ribosomal Proteins (Reads)")  +
  xlab("CUB Ribosomal Proteins (Assemblies)")
dev.off()











