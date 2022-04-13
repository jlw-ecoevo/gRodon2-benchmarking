## JLW 2022 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

# Functions --------------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  #aligns dataframes by specified column
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  #boxcox transform/bakctransform w/ user-specified lambda
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

prepXM <- function(xm,xmi,x1,x2,x3){
  #helper function to merge each mixture dataset, calculate errors
  xm$Genome <- xm$Genome %>% 
    gsub(pattern=".fna.gz$",replace="") %>%
    gsub(pattern=".fna$",replace="")
  xm <- xm %>% mutate(Genome=as.character(Genome),
                      Genomes=as.character(Genomes))
  xmi$Genome <- xmi$Genome %>% 
    gsub(pattern=".fna.gz$",replace="") %>%
    gsub(pattern=".fna$",replace="")
  xmi <- xmi %>% mutate(Genome=as.character(Genome))
  xm <- merge.easy(xm,xmi,key="Genome")
  
  for(i in 1:nrow(x1)){
    genomes <- xm$Genomes[xm$Genome==x1$Genome[i]] %>% 
      strsplit(split=",") %>% 
      unlist() %>%
      gsub(pattern=".fna.gz$",replace="") %>%
      gsub(pattern=".fna$",replace="")
    #x1$gc_avg[i] <- mean(x2$GC[x2$Genome %in% genomes])
    x1$dgpm_avg[i] <- mean(x2$d[x2$Genome %in% genomes])
    x1$dgp_avg[i] <- mean(x3$d[x3$Genome %in% genomes])
  }
  xm <- merge.easy(xm,x1,key="Genome") #%>%
  # sample_n(100)
  
  xm$d.gr <- xm$d.mt.i
  xm$d.gr[xm$Consistency<0.6] <- xm$d[xm$Consistency<0.6]
  xm$dgr_avg <- xm$dmt_avg_i
  xm$dgr_avg[xm$Consistency<0.6] <- xm$dm_avg[xm$Consistency<0.6]
  
  xm <- xm %>% 
    mutate(d_err=(log10(d)-log10(df_avg))^2,
           dmngct_err=(log10(d.mngct)-log10(df_avg))^2,
           dmt_err=(log10(d.mt)-log10(df_avg))^2,
           dmngcti_err=(log10(d.mngct.i)-log10(df_avg))^2,
           dmti_err=(log10(d.mt.i)-log10(df_avg))^2,
           dgr_err=(log10(d.gr)-log10(df_avg))^2,
           dgp_err=(log10(d.gp)-log10(df_avg))^2) %>%
    mutate(d_err_self=(log10(d)-log10(dm_avg))^2,
           dmngct_err_self=(log10(d.mngct)-log10(dmngct_avg))^2,
           dmt_err_self=(log10(d.mt)-log10(dmt_avg))^2,
           dmngcti_err_self=(log10(d.mngct.i)-log10(dmngct_avg_i))^2,
           dmti_err_self=(log10(d.mt.i)-log10(dmt_avg_i))^2,
           dgr_err_self=(log10(d.gr)-log10(dgr_avg))^2,
           dgp_err_self=(log10(d.gp)-log10(dgpm_avg))^2)
  xm <- xm %>% subset(!is.na(CUB.i)) %>% subset(df_avg<100)
}

# Load Data --------------------------------------------------------------------

setwd("~/gRodon2-benchmarking/Data/")
load("zou_simmeta10i_data.rda")
xmi <- xm
load("zou_simmeta_data.rda")
x1 <- read.delim("growthpred_zou_10.tbl",head=F)
names(x1) <- c("Genome","d.gp")
x1$GC <- NULL
x2 <- read.delim("growthpred_zou_genomes_meta.tbl",head=F)
names(x2) <- c("Genome","d")
x3 <- read.delim("growthpred_zou_genomes.tbl",head=F)
names(x3) <- c("Genome","d")
zou <- prepXM(xm,xmi,x1,x2,x3)

setwd("~/gRodon2-benchmarking/Data/")
load("gorg_simmeta10i_data.rda")
xmi <- xm
load("gorg_simmeta_data.rda")
x1 <- read.delim("growthpred_gc_gorg.tbl",head=F)
names(x1) <- c("Genome","d.gp")
x1$GC <- NULL
x2 <- read.delim("growthpred_gc_gorg_genomes_meta.tbl",head=F)
names(x2) <- c("Genome","d")
x3 <- read.delim("growthpred_gc_gorg_genomes.tbl",head=F)
names(x3) <- c("Genome","d")
gorg <- prepXM(xm,xmi,x1,x2,x3)

setwd("~/gRodon2-benchmarking/Data/")
load("refseq_simmeta10i_data.rda")
xmi <- xm
load("refseq_simmeta_data.rda")
x1 <- read.delim("growthpred_refseq10.tbl",head=F)
names(x1) <- c("Genome","d.gp")
x1$GC <- NULL
x2 <- read.delim("growthpred_gc_refseq_genomes_meta_GM.tbl",head=F)
names(x2) <- c("Genome","d")
x3 <- read.delim("growthpred_gc_refseq_genomes_G.tbl",head=F)
names(x3) <- c("Genome","d")
refseq <- prepXM(xm,xmi,x1,x2,x3)

# Plot -------------------------------------------------------------------------

#Plotting parameters
global_size <- 12
al <- 0.01

#Merge error info for all mixtures for first panel
xm1 <- refseq %>% subset(select=c(Consistency,dgr_err_self,d_err_self,d.gr,df_avg,dgr_avg,d,d.gp,GC))
xm2 <- gorg %>% subset(select=c(Consistency,dgr_err_self,d_err_self,d.gr,df_avg,dgr_avg,d,d.gp,GC))
xm3 <- zou %>% subset(select=c(Consistency,dgr_err_self,d_err_self,d.gr,df_avg,dgr_avg,d,d.gp,GC))
xm <- rbind(xm1,xm2,xm3)


pERRself <- ggplot(xm,aes(x=Consistency,dgr_err_self,color="MMv2")) + 
  geom_point(alpha=0.01) +
  geom_point(data=xm,aes(x=Consistency,d_err_self,color="MMv1"),alpha=0.01) +
  geom_smooth(fill="black",
              alpha=1) + 
  geom_smooth(data=xm,aes(x=Consistency,d_err_self,color="MMv1"),
              fill="black",
              alpha=1) + 
  scale_y_log10()+
  geom_vline(xintercept=0.6,lty=2) +
  labs(color="") +
  ylab("MSE Predicting Avg. Metagenome Mode of Mixture") +
  scale_y_log10() +
  theme_pubclean(base_size = global_size+6) +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "bottom")


pGCrefseq <- ggplot(data=refseq,
                    aes(x=GC,y=d,color="MMv1")) +
  geom_point(alpha=al) +
  geom_point(data=refseq,
             aes(x=GC,y=d.gr,color="MMv2"),
             alpha=al) +
  theme_pubclean(base_size = global_size) +
  geom_smooth(fill="black",
              alpha=1) +
  geom_smooth(data=refseq,
              aes(x=GC,y=d.gr,color="MMv2"),
              fill="black",
              alpha=1) +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_brewer(palette="Set1") + 
  ggtitle("RefSeq")+
  theme(legend.position = "none")

pGCzou <- ggplot(data=zou,
                 aes(x=GC,y=d,color="MMv1")) +
  geom_point(alpha=al) +
  geom_point(data=zou,
             aes(x=GC,y=d.gr,color="MMv2"),
             alpha=al) +
  theme_pubclean(base_size = global_size) +
  geom_smooth(fill="black",
              alpha=1) +
  geom_smooth(data=zou,
              aes(x=GC,y=d.gr,color="MMv2"),
              fill="black",
              alpha=1) +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_brewer(palette="Set1") + 
  ggtitle("Human Gut")+
  theme(legend.position = "none")+
  labs(color="")

pGCgorg <- ggplot(data=gorg,
                  aes(x=GC,y=d,color="MMv1")) +
  geom_point(alpha=al) +
  geom_point(data=gorg,
             aes(x=GC,y=d.gr,color="MMv2"),
             alpha=al) +
  theme_pubclean(base_size = global_size) +
  geom_smooth(fill="black",
              alpha=1) +
  geom_smooth(data=gorg,
              aes(x=GC,y=d.gr,color="MMv2"),
              fill="black",
              alpha=1) +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_brewer(palette="Set1") + 
  ggtitle("Ocean Surface")+
  theme(legend.position = "none")

#Reformat for error boxplots
refseq$conclass <- "> 0.6"
refseq$conclass[refseq$Consistency<0.6] <- "< 0.6"
refseqs <-  refseq %>%
  subset(select=c(conclass,
                  d_err,
                  dmti_err,
                  dgp_err,
                  dgr_err))
names(refseqs) <- c("conclass","MMv1","MMBC","Growthpred","MMv2")
refseqm <- refseqs %>%
  melt(id.vars="conclass")
pERRrefseq <- ggplot(refseqm,aes(y=value,fill=variable,x=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position="right") +
  ggtitle("RefSeq")+ 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  stat_compare_means(method="anova",label.x=1.2) +       
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "MMv2", hide.ns = F, label.y=1)  

#Reformat for error boxplots
zou$conclass <- "> 0.6"
zou$conclass[zou$Consistency<0.6] <- "< 0.6"
zous <-  zou %>%
  subset(select=c(conclass,
                  d_err,
                  dmti_err,
                  dgp_err,
                  dgr_err))
names(zous) <- c("conclass","MMv1","MMBC","Growthpred","MMv2")
zoum <- zous %>%
  melt(id.vars="conclass")
pERRzou <- ggplot(zoum,aes(y=value,fill=variable,x=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position="right") +
  ggtitle("Human Gut")+ 
  theme(legend.position = c(0.5,0.2),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  stat_compare_means(method="anova",label.x=1.2) +       
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "MMv2", hide.ns = F, label.y=1)  

#Reformat for error boxplots
gorg$conclass <- "> 0.6"
gorg$conclass[gorg$Consistency<0.6] <- "< 0.6"
gorgs <-  gorg %>%
  subset(select=c(conclass,
                  d_err,
                  dmti_err,
                  dgp_err,
                  dgr_err))
names(gorgs) <- c("conclass","MMv1","MMBC","Growthpred","MMv2")
gorgm <- gorgs %>%
  melt(id.vars="conclass")
pERRgorg <- ggplot(gorgm,aes(y=value,fill=variable,x=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position="right") +
  ggtitle("Ocean Surface")+ 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  stat_compare_means(method="anova",label.x=1.2) +       
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "MMv2", hide.ns = F, label.y=1)  


setwd("~/gRodon2-benchmarking/Figures/")
png(file="Compare_all_10i.png",width=1000,height=1000)
ggarrange(pERRself,
          ggarrange(ggarrange(pERRrefseq,
                              pERRgorg,
                              pERRzou,
                              ncol=3,
                              labels=c("(b)","(c)","(d)")),
                    ggarrange(pGCrefseq,
                              pGCgorg,
                              pGCzou,
                              ncol=3,
                              labels=c("(e)","(f)","(g)")),
                    nrow=2),
          ncol=2,
          widths=c(1,1),
          labels=c("(a)",""),
          hjust=0)
dev.off()

setwd("~/gRodon2-benchmarking/Figures/")
pdf(file="Compare_all_10i.pdf",width=14,height=14)
ggarrange(pERRself,
          ggarrange(ggarrange(pERRrefseq,
                              pERRgorg,
                              pERRzou,
                              ncol=3,
                              labels=c("(b)","(c)","(d)")),
                    ggarrange(pGCrefseq,
                              pGCgorg,
                              pGCzou,
                              ncol=3,
                              labels=c("(e)","(f)","(g)")),
                    nrow=2),
          ncol=2,
          widths=c(1,1),
          labels=c("(a)",""),
          hjust=0)
dev.off()

refseq.aov <- aov(value~variable,data=refseqm)
summary(refseq.aov)
TukeyHSD(refseq.aov)

zou.aov <- aov(value~variable,data=zoum)
summary(zou.aov)
TukeyHSD(zou.aov)

gorg.aov <- aov(value~variable,data=gorgm)
summary(gorg.aov)
TukeyHSD(gorg.aov)

# Supplemental - performance at low consistency --------------------------------

pCC3 <- ggplot(gorg,aes(x=df_avg,y=d.gr,fill=factor(Consistency<0.6))) +
  geom_point(alpha=0.5,size=3,pch=21) +
  geom_point(data=gorg %>% subset(Consistency<0.6),
             aes(x=df_avg,y=d.gr),
             alpha=1,size=3,pch=21,fill="red") +
  scale_fill_manual(values=c("white","red"),labels=c(">0.6", "<0.6")) +
  labs(fill="Consistency") +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_pubclean() +
  ylab("Prediction on Mixture") +
  xlab("Avg. Full Mode Prediction for Genomes in Mix.") +
  ggtitle("Marine Surface")

pCC1 <- ggplot(zou,aes(x=df_avg,y=d.gr,fill=factor(Consistency<0.6))) +
  geom_point(alpha=0.5,size=3,pch=21) +
  geom_point(data=zou %>% subset(Consistency<0.6),
             aes(x=df_avg,y=d.gr),
             alpha=1,size=3,pch=21,fill="red") +
  scale_fill_manual(values=c("white","red"),labels=c(">0.6", "<0.6")) +
  labs(fill="Consistency") +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_pubclean() +
  ylab("Prediction on Mixture") +
  xlab("Avg. Full Mode Prediction for Genomes in Mix.") +
  ggtitle("Human Gut")



pCC2 <- ggplot(zoum%>%subset(variable="MMv2"),aes(x=conclass,y=value)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position="right") +
  # ggtitle("Human Gut")+ 
  theme(legend.position = c(0.5,0.2)) +
  xlab("Consistency")

pCC4 <- ggplot(gorgm%>%subset(variable="MMv2"),aes(x=conclass,y=value)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Genomes in Mixture") +
  labs(fill="") +
  theme(legend.position="right") +
  # ggtitle("Marine Surface")+ 
  theme(legend.position = c(0.5,0.2)) +
  xlab("Consistency")

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Compare_consistency_10i.png",width=800,height=800)
ggarrange(ggarrange(pCC1,pCC2,ncol=2,widths=c(3,1),labels=c("(a)","(b)")),
          ggarrange(pCC3,pCC4,ncol=2,widths=c(3,1),labels=c("(c)","(d)")),
          nrow=2)
dev.off()





# Supplemental - constructing gRodon mmv2 --------------------------------------


setwd("~/gRodon2-benchmarking/Data/")
load("refseq_simmeta2i_data.rda")
xmi <- xm
load("refseq_simmeta2_data.rda")
x1 <- read.delim("growthpred_gc_refseq.tbl",head=F)
names(x1) <- c("Genome","d.gp")
x1$GC <- NULL
x2 <- read.delim("growthpred_gc_refseq_genomes_meta_GM.tbl",head=F)
names(x2) <- c("Genome","d")
x3 <- read.delim("growthpred_gc_refseq_genomes_G.tbl",head=F)
names(x3) <- c("Genome","d")
refseq <- prepXM(xm,xmi,x1,x2,x3)
xm <- refseq %>% mutate(di_err=(log10(d.i)-log10(df_avg))^2)

p1 <- ggplot(refseq,aes(x=Consistency,dmti_err_self,color="MMBC")) + 
  geom_point(alpha=0.01) +
  geom_point(data=refseq,aes(x=Consistency,d_err_self,color="MMv1"),
             alpha=0.01) +
  geom_point(data=refseq,aes(x=Consistency,dgr_err_self,color="MMv2"),
             alpha=0.01) +
  geom_smooth(fill="black",
              alpha=1,
              method="loess") + 
  geom_smooth(data=refseq,aes(x=Consistency,d_err_self,color="MMv1"),
              fill="black",
              alpha=1,
              method="loess") + 
  geom_smooth(data=refseq,aes(x=Consistency,dgr_err_self,color="MMv2"),
              fill="black",
              alpha=1,
              method="loess") + 
  scale_y_log10()+
  geom_vline(xintercept=0.6,lty=2) +
  labs(color="") +
  ylab("MSE Predicting Avg. Metagenome Mode of Mixture") +
  scale_y_log10() +
  theme_pubclean(base_size = global_size+6) +
  scale_color_manual(values=c("#377eb8","#e41a1c","#4daf4a")) +
  theme(legend.position = c(0.8,0.1))

refseq$conclass <- "> 0.6"
refseq$conclass[refseq$Consistency<0.6] <- "< 0.6"
refseqs <-  refseq %>%
  subset(select=c(conclass,
                  d_err,
                  dmti_err,
                  dgr_err))
names(refseqs) <- c("conclass","MMv1","MMBC","MMv2")
refseqm <- refseqs %>%
  melt(id.vars="conclass")
p4 <- ggplot(refseqm,aes(x=conclass,y=value,fill=variable)) +
  geom_boxplot() +
  scale_y_log10(limits=c(1e-11,1e2)) +
  theme_pubclean(base_size = global_size) +
  ylab("MSE Predicting Avg. gRodon Full Mode of Mixture") +
  labs(fill="") +
  theme(legend.position=c(0.5,0.1)) +
  scale_fill_manual(values=c("#e41a1c","#377eb8","#4daf4a")) +
  theme()


plot_lims <- c(1e-2,1e2)

p2 <- ggplot(data=xm,
               aes(x=dm_avg,y=d,fill="> 0.7")) +
  geom_point(pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.7),
             aes(x=dm_avg,y=d,fill="< 0.7"),
             pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.6),
             aes(x=dm_avg,y=d,fill="< 0.6"),
             pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.5),
             aes(x=dm_avg,y=d,fill="< 0.5"),
             pch=21,size=3) +
  scale_y_log10(limits=plot_lims) +
  scale_x_log10(limits=plot_lims) +
  theme_pubclean(base_size = global_size) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_fill_manual(values=c("red","orange","yellow","white")) +
  ylab("Prediction on Mixture") +
  xlab("Avg. Prediction for Genomes in Mix.") +
  labs(fill="Consistency") +
  theme(legend.position=c(0.25,0.8)) +
  ggtitle("MMv1")
p3 <- ggplot(data=xm,
               aes(x=dmt_avg_i,y=d.mt.i,fill="> 0.7")) +
  geom_point(pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.7),
             aes(x=dmt_avg_i,y=d.mt.i,fill="< 0.7"),
             pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.6),
             aes(x=dmt_avg_i,y=d.mt.i,fill="< 0.6"),
             pch=21,size=3) +
  geom_point(data=xm %>% subset(Consistency<0.5),
             aes(x=dmt_avg_i,y=d.mt.i,fill="< 0.5"),
             pch=21,size=3) +
  scale_y_log10(limits=plot_lims) +
  scale_x_log10(limits=plot_lims) +
  theme_pubclean(base_size = global_size) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_fill_manual(values=c("red","orange","yellow","white")) +
  ylab("Prediction on Mixture") +
  xlab("Avg. Prediction for Genomes in Mix.") +
  labs(fill="Consistency") +
  theme(legend.position="none") +
  ggtitle("MMBC")

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Constructing_v2_2i.png",width=1000,height=1000)
ggarrange(p1,
          ggarrange(p2,p3,nrow=2,
                    labels=c("(b)","(c)")),
          p4,
          ncol=3,
          widths=c(2,2,1),
          labels=c("(a)","","(d)"))
dev.off()

global_size <- 10
al <- 0.2
errlim <- c(1e-5,1e1)

pGCd <- ggplot(data=xm,
                 aes(x=GC,y=d,color=d_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMv1")+
  theme(legend.position = "none")+
  labs(color="")
pGCdi <- ggplot(data=xm,
               aes(x=GC,y=d.i,color=di_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMv1, Single-Gene")+
  theme(legend.position = "none")+
  labs(color="")
pGCdmti <- ggplot(data=xm,
                     aes(x=GC,y=d.mt.i,color=dmti_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("%GC") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMBC")+
  theme(legend.position = c(0.8,0.2))+
  labs(color="Squared Error")


setwd("~/gRodon2-benchmarking/Figures/")
png(file="Constructing_v2_2i_GC.png",width=750,height=500)
ggarrange(pGCd,pGCdi,pGCdmti,
          ncol=3,
          labels=c("(a)","(b)","(c)"))
dev.off()


pCd <- ggplot(data=xm,
               aes(x=Consistency,y=d,color=d_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("Consistency") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMv1")+
  theme(legend.position = "none")+
  labs(color="")
pCdi <- ggplot(data=xm,
                aes(x=Consistency,y=d.i,color=di_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("Consistency") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMv1, Single-Gene")+
  theme(legend.position = "none")+
  labs(color="")
pCdmti <- ggplot(data=xm,
                  aes(x=Consistency,y=d.mt.i,color=dmti_err)) +
  geom_point(alpha=al) +
  theme_pubclean(base_size = global_size) +
  # geom_smooth(fill="black",
  #             alpha=1,color="gray") +
  ylab("Prediction on Mixture") +
  scale_y_log10(limits=c(1e-3,1e2)) +
  xlab("Consistency") + 
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  ggtitle("MMBC")+
  theme(legend.position = c(0.8,0.2))+
  labs(color="Squared Error")


setwd("~/gRodon2-benchmarking/Figures/")
png(file="Constructing_v2_2i_Consistency.png",width=750,height=500)
ggarrange(pCd,pCdi,pCdmti,
          ncol=3,
          labels=c("(a)","(b)","(c)"))
dev.off()


global_size <- 10
al <- 0.05
errlim <- c(1e-5,1e1)


pCGCab <- ggplot(data=xm,
       aes(y=Consistency,x=abs(0.5-GC),color=d_err)) +
  geom_point(alpha=al,size=2) +
  theme_pubclean(base_size = global_size) +
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  theme(legend.position = "none")+
  labs(color="Squared Error") +
  xlab("|GC-0.5|")


pCGC <- ggplot(data=xm,
       aes(y=Consistency,x=GC,color=d_err)) +
  geom_point(alpha=al,size=2) +
  theme_pubclean(base_size = global_size) +
  scale_color_viridis_c(limits = errlim, oob = scales::squish,trans="log") +
  theme(legend.position = c(0.1,0.8))+
  labs(color="Squared Error")


pGCaberr <- ggplot(refseq,aes(x=abs(0.5-GC),dmti_err_self,color="MMBC")) + 
  geom_point(alpha=0.01) +
  geom_point(data=refseq,aes(x=abs(0.5-GC),d_err_self,color="MMv1"),
             alpha=0.01) +
  geom_point(data=refseq,aes(x=abs(0.5-GC),dgr_err_self,color="MMv2"),
             alpha=0.01) +
  geom_smooth(fill="black",
              alpha=1,
              method="loess") + 
  geom_smooth(data=refseq,aes(x=abs(0.5-GC),d_err_self,color="MMv1"),
              fill="black",
              alpha=1,
              method="loess") + 
  geom_smooth(data=refseq,aes(x=abs(0.5-GC),dgr_err_self,color="MMv2"),
              fill="black",
              alpha=1,
              method="loess") + 
  scale_y_log10()+
  labs(color="") +
  ylab("MSE Predicting Avg. Metagenome Mode of Mixture") +
  scale_y_log10() +
  theme_pubclean(base_size = global_size+6) +
  scale_color_manual(values=c("#377eb8","#e41a1c","#4daf4a")) +
  theme(legend.position = c(0.8,0.1)) +
  xlab("|GC-0.5|")

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Constructing_v2_2i_GC-Consistency.png",width=800,height=800)
ggarrange(ggarrange(pCGC,
                    pCGCab,
                    nrow=2,
                    labels=c("(a)","(b)")),
          pGCaberr,
          ncol=2,
          labels=c("","(c)"))
dev.off()