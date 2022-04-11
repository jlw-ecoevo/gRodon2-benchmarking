library(parallel)
library(gRodon)
library(dplyr)
library(Biostrings)

doAssembly <- function(path_assembly,out_folder,metadata){
  
  srr <- path_assembly %>% 
    gsub(pattern="metaeuk",replace="") %>% 
    gsub(pattern=".codon.fas",replace="")
  ffn_path <- srr %>% gsub(pattern="euk_a",replace="prok_a") %>% paste0(.,".ffn")
  cov_path <- ffn_path %>% gsub(pattern="prok_annot",replace="mapped") %>% gsub(pattern=".ffn",replace=".cov")
  bamcov <- read.delim(cov_path,head=F,stringsAsFactors=F) 

  #print(gsub("_.*","",basename(srr)))
  #print(basename(srr))
  #print(ffn_path)
  print(path_assembly)
  
  genes <- readDNAStringSet(ffn_path)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case=T)
  gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
  print(sum(highly_expressed))
  #stop()
  
  genes_euk <- readDNAStringSet(path_assembly)
  riboprot <- path_assembly %>% 
    gsub(pattern=".codon.fas",replace=".riboprot") %>%
    readLines()
  highly_expressed_euk <- names(genes_euk) %in% riboprot
  
  #print(sum(highly_expressed_euk))
  #stop()

  temperature <- metadata$Temperature[metadata$NCBI.SRA.Accession.for.assembled.metagenome.contigs==gsub("_.*","",basename(srr))]
  

  #print(temperature)
  if(sum(highly_expressed)>1 & sum(highly_expressed_euk)>1 & !is.na(temperature)){
    
    pgm <- predictGrowth(genes, 
                          highly_expressed, 
                          mode = "metagenome_v1",
                          training_set = "madin", 
                          temperature = temperature,
                          depth_of_coverage = bamcov$V7)
    pgm2 <- predictGrowth(genes,
                         highly_expressed,
                         mode = "metagenome_v2",
                         temperature = temperature,
                         depth_of_coverage = bamcov$V7)
    names(pgm2) <- paste0(names(pgm2),".2")
    #pge <- predictGrowth(genes_euk, 
    #                     highly_expressed_euk,
    #                     mode = "eukaryote",
    #                     temperature = temperature)
    #names(pge) <- paste0(names(pge),".eukaryotes")
    
    out <- c(list(assembly = path_assembly,
		  GC=gc,
                  pgm,
                  pgm2))                   
                  #pgmngct,
                  #pgmngcti))
                  #pge))
    
    save(out,file=paste0(out_folder,basename(path_assembly),".gRodon.rda"))
    
    return(out)
    
  } else {
    return(NULL)
  }
}


tryDoAssembly <- function(path_assembly,out_folder,metadata){
  try(doAssembly(path_assembly,out_folder,metadata))
}

setwd("/media/ink/BIOGEOTRACES")
load("BIOGEOTRACES_temp.RData")

system("mkdir /media/ink/BIOGEOTRACES/Out")

system('cd /media/ink/BIOGEOTRACES; find -type f -name "*.riboprot" | awk \'gsub(".riboprot",".codon.fas")\' | awk \'{gsub("\\\\./",""); print}\' > euk_files.txt')

gene_files <- paste0("/media/ink/BIOGEOTRACES/",
                     readLines("/media/ink/BIOGEOTRACES/euk_files.txt"))

out_folder <- "/media/ink/BIOGEOTRACES/Out/"


mclapply(gene_files,
       tryDoAssembly,
       out_folder = out_folder,
       metadata=BIOGEOTRACES_temp,
       mc.cores = 40)

