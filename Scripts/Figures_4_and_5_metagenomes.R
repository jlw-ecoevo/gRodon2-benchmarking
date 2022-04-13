## JLW 2022 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)

# Functions --------------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

prepXM <- function(xm,xmi,x1,x2,x3){
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

#BIOGEOTRACES
load("CodonStatistics_BIOGEOTRACES.RData")
met_df <- read.csv("GEOTRACES_metadata.csv")
met_df$Assembly <- paste0("metaeuk",
                          met_df$NCBI.SRA.Accession.for.assembled.metagenome.contigs,
                     "_",
                     met_df$Sample.name,
                     ".codon.fas")
bio_df <- merge.easy(bio_df,met_df,key="Assembly")

#HMASM
load("CodonStatistics_HMP.Rdata")
hmp_df$SRS.ID <- gsub(".ffn","",hmp_df$Assembly)
hmpmet_df <- read.csv("HMASM.csv")
hmp_df <- merge.easy(hmp_df,hmpmet_df,key="SRS.ID")

# BioGEOTRACES -----------------------------------------------------------------

pBIO1 <- ggplot(data=bio_df%>%subset(Depth..m.<100),aes(x=ConsistencyHE,fill="BIOGEOTRACES Surface")) +
  geom_density(alpha=0.5,color="white") +
  geom_density(data=gorg,aes(x=Consistency,fill="Marine Surface Mixtures"),alpha=0.5,color="white") +
  geom_density(data=refseq,aes(x=Consistency,fill="RefSeq Mixtures"),alpha=0.5,color="white") +
  scale_fill_manual(values=c("#40B0A6","#E1BE6A","gray")) +
  theme_pubclean() +
  geom_vline(xintercept=0.6,lty=2) +
  labs(fill="") +
  xlab("Consistency")

pBIO2 <- ggplot(bio_df,aes(x=Depth..m.,y=d.2)) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_log10() + 
  geom_point() +
  geom_smooth(color="black") +
  theme_pubclean() +
  coord_flip() +
  ylab("Predicted Average Community-Wide Min. Doubling Time (Hours)") +
  xlab("Depth (Meters)")

pBIO3 <- ggplot(data=bio_df,aes(x=d.2,fill="MMv2")) +
  geom_density(alpha=0.5,color="white") +
  geom_density(aes(x=d,fill="MMv1"),alpha=0.5,color="white") +
  xlab("Predicted Average Community-Wide Min. Doubling Time (Hours)") +
  scale_fill_manual(values=c("#FEFE62","#D35FB7")) +
  scale_x_log10() +
  labs(fill="") +
  theme_pubclean() 

setwd("~/gRodon2-benchmarking/Figures/")
png(file="BIOGEOTRACES.png",width=800,height=800)
ggarrange(pBIO1,
          ggarrange(pBIO2,
                   pBIO3,
                   ncol=2,
                   labels=c("(b)","(c)")),
          nrow=2,
          labels=c("(a)",""))
dev.off()

# HMASM ------------------------------------------------------------------------

pHMP1 <- ggplot(data=hmp_df%>%subset(Body.Site=="stool"),aes(x=ConsistencyHE,fill="HMP Stool")) +
  geom_density(alpha=0.5,color="white") +
  geom_density(data=zou,aes(x=Consistency,fill="Human Gut Mixtures"),alpha=0.5,color="white") +
  geom_density(data=refseq,aes(x=Consistency,fill="RefSeq Mixtures"),alpha=0.5,color="white") +
  scale_fill_manual(values=c("#40B0A6","#E1BE6A","gray")) +
  theme_pubclean() +
  geom_vline(xintercept=0.6,lty=2) +
  labs(fill="") +
  xlab("Consistency")

hmp_df$Body.Site <- gsub("_"," ",hmp_df$Body.Site)
sites_keep <- names(table(hmp_df$Body.Site))[table(hmp_df$Body.Site)>100]
pHMP2 <- ggplot(hmp_df %>% subset(Body.Site %in% sites_keep),
       aes(x=reorder(Body.Site,d.2),y=d.2,fill=Body.Site)) +
  geom_boxplot() +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=60,hjust=1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("Body Site") +
  ylab("Predicted Avg Community Min. DT (Hours)") +
  scale_y_log10() +
  stat_compare_means(method="anova",label.x=1.2,label.y=1.1) +       
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "stool", hide.ns = F, label.y=0.8)  
pHMP2


pHMP3 <- ggplot(data=hmp_df,aes(x=d.2,fill="MMv2")) +
  geom_density(alpha=0.5,color="white") +
  geom_density(aes(x=d,fill="MMv1"),alpha=0.5,color="white") +
  xlab("Predicted Average Community-Wide Min. Doubling Time (Hours)") +
  scale_fill_manual(values=c("#FEFE62","#D35FB7")) +
  labs(fill="") +
  scale_x_log10() +
  theme_pubclean() 


setwd("~/gRodon2-benchmarking/Figures/")
png(file="HMP.png",width=800,height=800)
ggarrange(pHMP1,
          ggarrange(pHMP2,
                    pHMP3,
                    ncol=2,
                    labels=c("(b)","(c)"),
                    hjust=0),
          nrow=2,
          labels=c("(a)",""),
          hjust=0)
dev.off()


# Stats ------------------------------------------------------------------------

x <- hmp_df %>% subset(Body.Site %in% sites_keep)
hmp.aov <- aov(log10(d.2)~Body.Site,data=x)
summary(hmp.aov)
TukeyHSD(hmp.aov)

wilcox.test(hmp_df$ConsistencyHE[hmp_df$Body.Site=="stool"],zou$Consistency)
wilcox.test(bio_df$ConsistencyHE.2[bio_df$Depth..m.<100],gorg$Consistency)

wilcox.test(hmp_df$d.2,bio_df$d.2)
mean(hmp_df$d.2)
mean(bio_df$d.2)

min(bio_df$d.2)
max(hmp_df$d.2)

table(hmp_df$d.2<min(bio_df$d.2))/sum(table(hmp_df$d.2<min(bio_df$d.2)))

