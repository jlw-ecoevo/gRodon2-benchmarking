
## JLW 2022 

# Load Packages ----------------------------------------------------------------


library(Biostrings)
library(coRdon)
library(robustbase)

# Functions --------------------------------------------------------------------


truncateGene <- function(gene,len_truncate){
  DNAString(substr(as.character(gene),1,len_truncate)) %>%
    DNAStringSet() %>%
    return()
}

CUBTruncateGenome <- function(genome,len_truncate){
  genes_tl <- genome %>% lapply(truncateGene,len_truncate) 
  names(genes_tl) <- NULL
  genes_t <- do.call(c,genes_tl)
  MILC(codonTable(genes_t)) %>% return()
}

doTruncate <- function(genes,truncate_lengths){
  genes <- genes[width(genes)>=max(truncate_lengths)]
  cubmat <- numeric()
  for(i in truncate_lengths){
    cubmat <- cbind(cubmat,
                    CUBTruncateGenome(genes,i))
  }
  return(cubmat)
}

plotCUBs <- function(cubmat,truncate_lengths,main=""){
  al <- 0.1
  ul <- cubmat[,length(truncate_lengths)] >= mean(cubmat[,length(truncate_lengths)])
  plot(truncate_lengths,
       cubmat[1,],type="l",
       ylim=c(0.1,1.5),
       col=rgb(ul[1],0,1-ul[1],al),
       xlab="Gene Length",
       ylab="CUB (MILC)",
       main=main)
  for(i in 2:nrow(cubmat)){
    lines(truncate_lengths,cubmat[i,],col=rgb(ul[i],0,1-ul[i],al))
  }
  lines(truncate_lengths,colMedians(cubmat),lwd=2)
  CI <- apply(X=cubmat,MARGIN=2,FUN=quantile,c(0.025,0.975))
  lines(truncate_lengths,CI[1,],lwd=1)
  lines(truncate_lengths,CI[2,],lwd=1)
  abline(v=120,lty=2)
  abline(v=240,lty=2)
  abline(h=median(cubmat[,ncol(cubmat)]),lty=2)
}

plotCUBHEs <- function(cubmat,truncate_lengths,main=""){
  al <- 0.5
  ul <- cubmat[,length(truncate_lengths)] >= mean(cubmat[,length(truncate_lengths)])
  plot(truncate_lengths,
       cubmat[1,],type="l",
       ylim=c(0.1,1.5),
       xlab="Gene Length",
       ylab="CUB Ribosomal Proteins (MILC)",
       main=main,
       col=rgb(0,0,0,al),
       lwd=2)
  for(i in 2:nrow(cubmat)){
    lines(truncate_lengths,cubmat[i,],
          col=rgb(0,0,0,al),
          lwd=2)
  }
  lines(truncate_lengths,colMedians(cubmat),lwd=5)
  abline(v=120,lty=2)
  abline(v=240,lty=2)
  abline(h=median(cubmat[,ncol(cubmat)]),lty=2)
}

# Load Data --------------------------------------------------------------------


setwd("~/gRodon2-benchmarking/Data/")
genesB <- readDNAStringSet("GCF_000009045.1_ASM904v1_cds_from_genomic.fna.gz")
genesPM <- readDNAStringSet("GCF_000007925.1_ASM792v1_cds_from_genomic.fna.gz")
genesS <- readDNAStringSet("GCF_000012345.1_ASM1234v1_cds_from_genomic.fna.gz")
genesPC <- readDNAStringSet("GCF_020735445.1_ASM2073544v1_cds_from_genomic.fna.gz")

# Truncate ---------------------------------------------------------------------

truncate_lengths <- c(30,60,90,120,150,180,210,240,270,300,
                      330,360,390,420,450,480,510)
cubmatB <- doTruncate(genesB,truncate_lengths)
cubmatPM <- doTruncate(genesPM,truncate_lengths)
cubmatS <- doTruncate(genesS,truncate_lengths)
cubmatPC <- doTruncate(genesPC,truncate_lengths)

# Plot ------------------------------------------------------------------------

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Truncate_genes.png",width=800,height=800)
par(mfrow=c(2,2))
plotCUBs(cubmatB,truncate_lengths,main="Bacillus subtilis str. 168")
plotCUBs(cubmatPC,truncate_lengths,main="Prevotella copri DSM18205")
plotCUBs(cubmatPM,truncate_lengths,main="Prochlorococcus marinus CCMP1375")
plotCUBs(cubmatS,truncate_lengths,main="Pelagibacter ubique HTCC1062")
par(mfrow=c(1,1))
dev.off()

setwd("~/gRodon2-benchmarking/Figures/")
png(file="Truncate_genes_HE.png",width=800,height=800)
par(mfrow=c(2,2))
heB <- grepl("ribosomal protein",names(genesB[width(genesB)>=510]),ignore.case=T)
cubhematB <- cubmatB[heB,]
plotCUBHEs(cubhematB,truncate_lengths,main="Bacillus subtilis str. 168")
lines(truncate_lengths,colMedians(cubmatB[,]),type="l",col="red",lwd=3)

hePC <- grepl("ribosomal protein",names(genesPC[width(genesPC)>=510]),ignore.case=T)
cubhematPC <- cubmatPC[hePC,]
plotCUBHEs(cubhematPC,truncate_lengths,main="Prevotella copri DSM18205")
lines(truncate_lengths,colMedians(cubmatPC[,]),type="l",col="red",lwd=3)

hePM <- grepl("ribosomal protein",names(genesPM[width(genesPM)>=510]),ignore.case=T)
cubhematPM <- cubmatPM[hePM,]
plotCUBHEs(cubhematPM,truncate_lengths,main="Prochlorococcus marinus CCMP1375")
lines(truncate_lengths,colMedians(cubmatPM[,]),type="l",col="red",lwd=3)

heS <- grepl("ribosomal protein",names(genesS[width(genesS)>=510]),ignore.case=T)
cubhematS <- cubmatS[heS,]
plotCUBHEs(cubhematS,truncate_lengths,main="Pelagibacter ubique HTCC1062")
lines(truncate_lengths,colMedians(cubmatS[,]),type="l",col="red",lwd=3)
par(mfrow=c(1,1))
dev.off()


