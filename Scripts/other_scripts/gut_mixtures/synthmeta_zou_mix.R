library(gRodon)
library(Biostrings)
library(coRdon)
library(dplyr)
library(parallel)
library(matrixStats)

setwd("/home/jake/zou/Assemblies/")
load("zou_genome_data.rda")

makeSim <- function(i,genome_df,n_sim){
  print(i)
  sample_i <- sample(1:length(genome_df$Genome),n_sim)
  abs_abundances <- rlnorm(n_sim)
  rel_abundances <- abs_abundances/sum(abs_abundances)
  
  genes <- readDNAStringSet(as.character(genome_df$Genome[sample_i][1]))
  abundance_mix <- numeric(length(genes))+rel_abundances[1]
  for(i in 2:length(sample_i)){
    genes_i <- readDNAStringSet(as.character(genome_df$Genome[sample_i][i]))
    abundance_mix <- c(abundance_mix,numeric(length(genes_i))+rel_abundances[i])
    genes <- c(genes,genes_i)
  }
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T) 
  
  
  print(table(abundance_mix))
  print(length(abundance_mix))
  print(length(genes))
  print("b")
  
  growth_mix <- predictGrowth(genes,highly_expressed,mode="metagenome",training_set="madin")
  growthmti_mix <- predictGrowth(genes,highly_expressed,mode="meta_testing",bg="individual") 
  growth.depth_mix <- predictGrowth(genes,highly_expressed,mode="metagenome",training_set="madin",depth_of_coverage = abundance_mix)
  growthmt.depthi_mix <- predictGrowth(genes,highly_expressed,mode="meta_testing",depth_of_coverage = abundance_mix,bg="individual") 

  gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))

  return(data.frame(GC=gc,
                    mix.d=growth_mix$d,
                    mix.d.mt.i=growthmti_mix$d,
                    mix.d.depth=growth.depth_mix$d,
                    mix.d.mt.depth.i=growthmt.depthi_mix$d, 
                    Genomes=paste(genome_df$Genome[sample_i],collapse=","), 
                    mix.Consistency=growth_mix$ConsistencyHE,
                    mix.Consistency.depth=growth.depth_mix$ConsistencyHE,
                    gc_avg=sum(genome_df[abundance$V1,"GC"]*abundance$V2),
                    nHE_mean=sum(genome_df[abundance$V1,"nHE"]*abundance$V2),
                    df_avg=weightedMedian(genome_df[abundance$V1,"d.f"],w=abundance$V2),
                    dp_avg=weightedMedian(genome_df[abundance$V1,"d.p"],w=abundance$V2), 
                    dm_avg=weightedMedian(genome_df[abundance$V1,"d.m"],w=abundance$V2),
                    dmt_avg=weightedMedian(genome_df[abundance$V1,"d.mt"],w=abundance$V2),
                    dmti_avg=weightedMedian(genome_df[abundance$V1,"d.mt.i"],w=abundance$V2), 
                    dm_avg_flat=weightedMedian(genome_df[abundance$V1,"d.m"],w=(numeric(10)+0.1)), 
                    dmt_avg_flat=weightedMedian(genome_df[abundance$V1,"d.mt"],w=(numeric(10)+0.1)),
                    dmti_avg_flat=weightedMedian(genome_df[abundance$V1,"d.mt.i"],w=(numeric(10)+0.1)), 
                    dmngct_avg_flat=weightedMedian(genome_df[abundance$V1,"d.mngct"],w=(numeric(10)+0.1)), 
                    dmngct_avg=weightedMedian(genome_df[abundance$V1,"d.mngct"],w=abundance$V2)))
}

n <- 1000
n_sim <- 10
trySim <- function(i,genome_df,n_sim){
  try(makeSim(i,genome_df,n_sim))
}
pair_list <- mclapply(1:n, 
         trySim,
         genome_df = genome_df,
         n_sim = n_sim,
         mc.cores = 20)

xm <- do.call("rbind",pair_list) %>% 
  as.data.frame(stringsAsFactors=F)

save(xm,file="zou_synthmix_data.rda")

#warnings()
