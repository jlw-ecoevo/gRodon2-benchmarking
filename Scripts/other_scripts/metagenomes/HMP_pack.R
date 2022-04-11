

library(dplyr)
library(data.table)
library(Biostrings)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}


system('cd /media/ink/HMP/Out; find -type f -name "*.rda" | awk \'{gsub("./",""); print}\' > ../HMP_out_files.txt')
out_files <- paste0("/media/ink/HMP/Out/",
                     readLines("/media/ink/HMP/HMP_out_files.txt"))
mag_list <- list()
for(i in 1:length(out_files)){
  load(out_files[i])
  mag_list[[i]] <- c(assembly=out[[1]],out[[2]],out[[3]],out[[4]])
}
print(mag_list[[1]])
hmp_df <- do.call("rbind",mag_list) %>% as.data.frame(stringsAsFactors=F) %>%
  mutate_all(unlist)
hmp_df$Assembly <- basename(hmp_df$assembly) %>% 
  gsub(pattern = "[.]trinity.*", replace = "")
hmp_df$assembly <- NULL


setwd("/media/ink/HMP/")
save(hmp_df,file="CodonStatistics_HMP.RData")
