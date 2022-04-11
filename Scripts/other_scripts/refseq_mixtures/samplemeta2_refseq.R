library(gRodon)
library(Biostrings)
library(coRdon)
library(dplyr)
library(parallel)

setwd("/media/ink/pair_metagenomes/genomes/genomes_to_sample/")
genomes <- readLines("genome.list")
load("refseq_genome_data.rda")

makeSim <- function(i,genome_df,samp_size){  
  print(i)
  genome <- paste0("sim2_",i,".fna")
  sample_i <- sample(1:nrow(genome_df),samp_size)
  system(paste0("cat ",
                paste(genome_df$Genome[sample_i],collapse=" "),
                " > ",
                "sim2_",i,".fna"))
  #system(paste0("cat ",
  #              paste(paste0(genome_df$Genome[sample_i],".ribo"),collapse=" "),
  #              " > ",
  #              "sim_",i,".fna.ribo"))
  #ribo <- readLines(paste0(genome,".ribo")) %>% gsub(pattern="^>",replace="")
  genes <- readDNAStringSet(genome)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T) 
  #highly_expressed <- names(genes) %in% ribo
  ribo_genes <- genes[highly_expressed]
  writeXStringSet(ribo_genes,paste0(genome,".ribo.fasta")) 
  milc <- MILC(codonTable(genes),
               subsets = list(HE = highly_expressed),
               id_or_name2 = "11")
  #CUBdist_HE <- milc[highly_expressed, 2] %>% mean()
  #print(CUBdist_HE)
  growth <- predictGrowth(genes,highly_expressed,mode="metagenome",training_set="madin")
  growthmt <- predictGrowth(genes,highly_expressed,mode="meta_testing") 
  growthmngct <- predictGrowth(genes,highly_expressed,mode="meta_nogc_testing") 
  #growth.t <- predictGrowth(genes,highly_expressed,mode="metagenome",training_set="madin",temperature=mean(c(genome_df$OGT[sample_i]))) 
  gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
  return(data.frame(GC=gc,
                    d=growth$d,
                    d.mt=growthmt$d,
                    d.mngct=growthmngct$d,
                    #d.t=growth.t$d,
                    #OGT=mean(c(genome_df$OGT[sample_i])),
                    Genome=genome,
                    Genomes=paste(genome_df$Genome[sample_i],collapse=","), 
                    CUBHE=growth$CUBHE,
                    CUB=growth$CUB,
                    Consistency=growth$ConsistencyHE,
                    #CUBdist_HE=CUBdist_HE,
                    nHE=growth$nHE,
                    gc_avg=mean(c(genome_df$GC[sample_i])),
                    gc_var=var(c(genome_df$GC[sample_i])),
                    nHE_mean=mean(c(genome_df$nHE[sample_i])),
                    nHE_var=var(c(genome_df$nHE[sample_i])),
                    df_avg=median(c(genome_df$d.f[sample_i])),
                    dp_avg=median(c(genome_df$d.p[sample_i])), 
                    dm_avg=median(c(genome_df$d.m[sample_i])),
                    dmt_avg=median(c(genome_df$d.mt[sample_i])),
                    dmngct_avg=median(c(genome_df$d.mngct[sample_i])),
                    #dft_avg=median(c(genome_df$d.f.t[sample_i])),
                    #dpt_avg=median(c(genome_df$d.p.t[sample_i])),
                    #dmt_avg=median(c(genome_df$d.m.t[sample_i])),
                    #dgct_avg=median(c(genome_df$d.gc.t[sample_i])), 
                    df_var=var(c(genome_df$d.f[sample_i])),
                    dp_var=var(c(genome_df$d.p[sample_i])),
                    dm_var=var(c(genome_df$d.m[sample_i]))))
}

n <- 10000
samp_size <- 2
trySim <- function(i,genome_df,samp_size){
  try(makeSim(i,genome_df,samp_size))
}
pair_list <- mclapply(1:n, 
         trySim,
         genome_df = genome_df,
         samp_size = samp_size,
         mc.cores = 60)

xm <- do.call("rbind",pair_list) %>% 
  as.data.frame(stringsAsFactors=F)

save(xm,file="refseq_simmeta2_data.rda")
