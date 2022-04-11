library(gRodon)
library(dplyr)
library(Biostrings)

setwd("/media/ink/pair_metagenomes/genomes/genomes_to_sample/")
genomes <- readLines("genome.list")
#load("gc_model.rda")

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}



getGenomeDat <- function(i,genomes){
  genome <- genomes[i]
  genes <- readDNAStringSet(genome)
  highly_expressed <-  grepl("ribosomal protein",names(genes),ignore.case = T)
  ribo_genes <- genes[highly_expressed] 
  writeXStringSet(ribo_genes,paste0(genome,".ribo.fasta"))  
  growth <- predictGrowth(genes,highly_expressed,mode="full",training_set="madin")
  growthp <- predictGrowth(genes,highly_expressed,mode="partial",training_set="madin")
  growthm <- predictGrowth(genes,highly_expressed,mode="metagenome",training_set="madin")
  growthmt <- predictGrowth(genes,highly_expressed,mode="meta_testing")
  growthmngct <- predictGrowth(genes,highly_expressed,mode="meta_nogc_testing")
  growthmt.i <- predictGrowth(genes,highly_expressed,mode="meta_testing",bg="individual")
  growthmngct.i <- predictGrowth(genes,highly_expressed,mode="meta_nogc_testing",bg="individual")
  growth$gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
  #growth$GCdiv <- abs(growth$gc-0.5)
  #growth$dgc <- boxcoxTransform(predict(mnew_milcgc,growth),
  #                          lambda_milcgc,
  #                          back_transform = TRUE)
  return(data.frame(Genome=genome,
                                GC=growth$gc,
                                CUBHE=growth$CUBHE,
                                CUB=growth$CUB,
                                ConsistencyHE=growth$ConsistencyHE,
                                nHE=growth$nHE,
                                d.f=growth$d,
                                d.p=growthp$d,
                                d.m=growthm$d,
                                d.mt=growthmt$d,
                                d.mngct=growthmngct$d,
                                d.mt.i=growthmt.i$d, 
                                d.mngct.i=growthmngct.i$d))
}

tryG <- function(i,genomes){
  try(getGenomeDat(i,genomes))
}

genome_list <- mclapply(1:length(genomes), 
                     tryG,
                     genomes = genomes,
                     mc.cores = 60)

genome_df <- do.call("rbind",genome_list) %>% 
  as.data.frame(stringsAsFactors=F)

rownames(genome_df) <- genome_df$Genome

save(genome_df,file="refseq_genome_data.rda")
