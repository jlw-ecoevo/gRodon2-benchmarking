

library(dplyr)
library(data.table)
library(Biostrings)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}


system('cd /media/ink/BIOGEOTRACES/Out; find -type f -name "*.rda" | awk \'{gsub("./",""); print}\' > ../BIOGEOTRACES_out_files.txt')
out_files <- paste0("/media/ink/BIOGEOTRACES/Out/",
                     readLines("/media/ink/BIOGEOTRACES/BIOGEOTRACES_out_files.txt"))
mag_list <- list()
for(i in 1:length(out_files)){
  load(out_files[i])
  mag_list[[i]] <- c(assembly=out[[1]],out[[2]],out[[3]],out[[4]])
}
print(mag_list[[1]])
bio_df <- do.call("rbind",mag_list) %>% as.data.frame(stringsAsFactors=F) %>%
  mutate_all(unlist)
bio_df$Assembly <- basename(bio_df$assembly) %>% 
  gsub(pattern = "[.]trinity.*", replace = "")
bio_df$assembly <- NULL


setwd("/media/ink/BIOGEOTRACES/")
save(bio_df,file="CodonStatistics_BIOGEOTRACES.RData")
