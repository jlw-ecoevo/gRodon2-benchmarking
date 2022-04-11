library(gRodon)
library(Biostrings)
library(coRdon)
library(dplyr)
library(parallel)

setwd("/media/ink/HMP/CDS/")
genomes <- readLines("../ffn_files.txt")

makeSim <- function(i,out_folder){  

  genes <- readDNAStringSet(i)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T) 
  bamcov <- read.delim(paste0("/media/ink/HMP/mapped/",gsub(".ffn",".cov",i)),head=F,stringsAsFactors=F)   


  growth <- predictGrowth(genes,highly_expressed,mode="metagenome_v1",training_set="madin",temperature=37,depth_of_coverage = bamcov$V7)
  growth2 <- predictGrowth(genes,highly_expressed,mode="metagenome_v2",temperature=37,depth_of_coverage = bamcov$V7)
  names(growth2) <- paste0(names(growth2),".2")
  gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
  
  out <- c(list(assembly=i,
                GC=gc,
                growth,
                growth2))

  save(out,file=paste0("/media/ink/HMP/Out/",basename(i),".gRodon.rda"))
  return(out)
}

trySim <- function(i){
  try(makeSim(i))
}
pair_list <- mclapply(genomes, 
         trySim,
         mc.cores = 40)

#xm <- do.call("rbind",pair_list) %>% 
#  as.data.frame(stringsAsFactors=F)
#
#save(xm,file="CodonStatistics_HMP.rda")
