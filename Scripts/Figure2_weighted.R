## JLW 2022 

# Load Packages ----------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)

# Functions --------------------------------------------------------------------

prepXM <- function(xm){
  #Helper function to calculate errors
  xm$mix.d.gr <- xm$mix.d.mt.i
  xm$mix.d.gr[xm$Consistency<0.6] <- xm$mix.d[xm$mix.Consistency<0.6]
  xm$mix.d.gr.depth <- xm$mix.d.mt.depth.i
  xm$mix.d.gr.depth[xm$mix.Consistency.depth<0.6] <- 
    xm$mix.d.depth[xm$mix.Consistency.depth<0.6]
  xm$dgr_err <- (log10(xm$df_avg)-log10(xm$mix.d.gr))^2
  xm$dgr_err_depth <- (log10(xm$df_avg)-log10(xm$mix.d.gr.depth))^2
  return(xm)
}


# Load Data --------------------------------------------------------------------

setwd("~/gRodon2-benchmarking/Data/")
load("zou_synthmix_data.rda")
zou <- prepXM(xm)

setwd("~/gRodon2-benchmarking/Data/")
load("gorg_synthmix_data.rda")
gorg <- prepXM(xm)

setwd("~/gRodon2-benchmarking/Data/")
load("refseq_synthmix_data.rda")
refseq <- prepXM(xm)

# Plot -------------------------------------------------------------------------

global_size <- 12
al <- 0.02

p1refseq <- ggplot(refseq,aes(x=df_avg,y=mix.d.gr,color="mix.d")) + 
  geom_point(alpha=al) +
  geom_point(data=refseq,
             aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
             alpha=al) +
  geom_smooth(aes(color="mix.d"),fill="black",alpha=1) +
  geom_smooth(data=refseq,
              aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
              fill="black",alpha=1) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_pubclean(base_size = global_size) +
  theme(legend.position=c(0.8,0.2)) +
  xlab("Full Mode Prediction on Mixture")+
  ylab("Prediction with gRodon Metagenome Mode v2") +
  ggtitle("RefSeq") +
  labs(color="") +
  scale_color_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                     values=c("#994F00","#006CD1")) + 
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1")) + 
  guides(fill = FALSE) 

refseq$conclass <- "> 0.6"
refseq$conclass[refseq$mix.Consistency.depth<0.6] <- "< 0.6"
refseqs <-  refseq %>%
  subset(select=c(conclass,
                  dgr_err,
                  dgr_err_depth))
names(refseqs) <- c("conclass","No Depth","Depth")
refseqm <- refseqs %>%
  melt(id.vars="conclass")
p2refseq <- ggplot(refseqm,aes(x=conclass,y=value,fill=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. Full Mode of Mixture") +
  labs(fill="") +
  xlab("Consistency") + 
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 60, vjust = 0.5)) +
  labs(color="") +
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1")) +
  stat_compare_means(method="wilcox.test",label.x=1.2,label = "p.signif")



p1zou <- ggplot(zou,aes(x=df_avg,y=mix.d.gr,color="mix.d")) + 
  geom_point(alpha=al) +
  geom_point(data=zou,
             aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
             alpha=al) +
  geom_smooth(aes(color="mix.d"),fill="black",alpha=1) +
  geom_smooth(data=zou,
              aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
              fill="black",alpha=1) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_pubclean(base_size = global_size) +
  theme(legend.position="none") +
  xlab("Full Mode Prediction on Mixture")+
  ylab("Prediction with gRodon Metagenome Mode v2") +
  ggtitle("Human Gut") +
  labs(color="") +
  scale_color_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                     values=c("#994F00","#006CD1")) + 
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1")) + 
  guides(fill = FALSE) 

zou$conclass <- "> 0.6"
zou$conclass[zou$mix.Consistency.depth<0.6] <- "< 0.6"
zous <-  zou %>%
  subset(select=c(conclass,
                  dgr_err,
                  dgr_err_depth))
names(zous) <- c("conclass","No Depth","Depth")
zoum <- zous %>%
  melt(id.vars="conclass")
p2zou <- ggplot(zoum,aes(x=conclass,y=value,fill=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. Full Mode of Mixture") +
  labs(fill="") +
  xlab("Consistency") + 
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 60, vjust = 0.5)) +
  labs(color="") +
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1")) +
  stat_compare_means(method="wilcox.test",label.x=1.2,label = "p.signif")



p1gorg <- ggplot(gorg,aes(x=df_avg,y=mix.d.gr,color="mix.d")) + 
  geom_point(alpha=al) +
  geom_point(data=gorg,
             aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
             alpha=al) +
  geom_smooth(aes(color="mix.d"),fill="black",alpha=1) +
  geom_smooth(data=gorg,
              aes(x=df_avg,y=mix.d.gr.depth,color="mix.d.depth"),
              fill="black",alpha=1) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_x_log10() +
  scale_y_log10(limits=c(1e-1,1e2)) +
  theme_pubclean(base_size = global_size) +
  theme(legend.position="none") +
  xlab("Full Mode Prediction on Mixture")+
  ylab("Prediction with gRodon Metagenome Mode v2") +
  ggtitle("Ocean Surface") +
  labs(color="") +
  scale_color_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                     values=c("#994F00","#006CD1")) + 
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1")) + 
  guides(fill = FALSE) 

gorg$conclass <- "> 0.6"
gorg$conclass[gorg$mix.Consistency.depth<0.6] <- "< 0.6"
gorgs <-  gorg %>%
  subset(select=c(conclass,
                  dgr_err,
                  dgr_err_depth))
names(gorgs) <- c("conclass","No Depth","Depth")
gorgm <- gorgs %>%
  melt(id.vars="conclass")
p2gorg <- ggplot(gorgm,aes(x=conclass,y=value,fill=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. Full Mode of Mixture") +
  labs(fill="") +
  xlab("Consistency") + 
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 60, vjust = 0.5)) +
  labs(color="") +
  scale_fill_manual(labels=c("Unweighted Prediction","Weighted Prediction"),
                    values=c("#994F00","#006CD1"))  +
  stat_compare_means(method="wilcox.test",label.x=1.2,label = "p.signif")


setwd("~/gRodon2-benchmarking/Figures/")
png(file="Compare_all_10depths.png",width=570,height=1000)
ggarrange(ggarrange(p1refseq,p2refseq,widths=c(2,1),ncol=2,labels=c("(a)","(b)")),
          ggarrange(p1gorg,p2gorg,widths=c(2,1),ncol=2,labels=c("(c)","(d)")),
          ggarrange(p1zou,p2zou,widths=c(2,1),ncol=2,labels=c("(e)","(f)")),
          nrow=3)
dev.off()









