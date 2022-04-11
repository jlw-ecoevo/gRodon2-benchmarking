library(gRodon)
library(Biostrings)
library(coRdon)
library(dplyr)
library(parallel)
library(matrixStats)

setwd("/media/ink/gorg/Genes10plus/")
load("gorg_genome_data.rda")

makeSim <- function(i,genome_df,n_sim){
  print(i)
  sample_i <- sample(1:length(genome_df$Genome),n_sim)
  abs_abundances <- rlnorm(n_sim)
  rel_abundances <- abs_abundances/sum(abs_abundances)
  
  genes <- readDNAStringSet(as.character(genome_df$Genome[sample_i][1]))
  ribo <- readLines(paste0(as.character(genome_df$Genome[sample_i][1]),".ribo")) %>% gsub(pattern="^>",replace="")
  #genes <- genes[!grepl("hypothetical",names(genes),ignore.case = T)]
  abundance_mix <- numeric(length(genes))+rel_abundances[1]
  for(i in 2:length(sample_i)){
    genes_i <- readDNAStringSet(as.character(genome_df$Genome[sample_i][i]))
    ribo_i <- readLines(paste0(as.character(genome_df$Genome[sample_i][1]),".ribo")) %>% gsub(pattern="^>",replace="") 
    #genes_i <- genes_i[!grepl("hypothetical",names(genes_i),ignore.case = T)]
    abundance_mix <- c(abundance_mix,numeric(length(genes_i))+rel_abundances[i])
    genes <- c(genes,genes_i)
    ribo <- c(ribo,ribo_i)
  }
  highly_expressed <- names(genes) %in% ribo 
  
  
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
                    gc_avg=sum(genome_df[sample_i,"GC"]*rel_abundances),
                    nHE_mean=sum(genome_df[sample_i,"nHE"]*rel_abundances),
                    df_var=var(genome_df[sample_i,"d.f"]),
                    dm_var=var(genome_df[sample_i,"d.m"]),
                    dmti_var=var(genome_df[sample_i,"d.mt.i"]),
                    df_avg=weightedMedian(genome_df[sample_i,"d.f"],w=rel_abundances),
                    dm_avg=weightedMedian(genome_df[sample_i,"d.m"],w=rel_abundances),
                    dmti_avg=weightedMedian(genome_df[sample_i,"d.mt.i"],w=rel_abundances), 
                    df_mean=weightedMean(genome_df[sample_i,"d.f"],w=rel_abundances),
                    dm_mean=weightedMean(genome_df[sample_i,"d.m"],w=rel_abundances),
                    dmti_mean=weightedMean(genome_df[sample_i,"d.mt.i"],w=rel_abundances),
                    dm_avg_flat=weightedMedian(genome_df[sample_i,"d.m"],w=(numeric(n_sim)+1/n_sim)), 
                    dmti_avg_flat=weightedMedian(genome_df[sample_i,"d.mt.i"],w=(numeric(n_sim)+1/n_sim))))
}

n <- 10000
n_sim <- 10
trySim <- function(i,genome_df,n_sim){
  try(makeSim(i,genome_df,n_sim))
}
pair_list <- mclapply(1:n, 
         trySim,
         genome_df = genome_df,
         n_sim = n_sim,
         mc.cores = 10)

xm <- do.call("rbind",pair_list) %>% 
  as.data.frame(stringsAsFactors=F)

save(xm,file="gorg_synthmix_data.rda")

#warnings()
