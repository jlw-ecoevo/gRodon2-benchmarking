x <- readLines("genome.list")

nsim <- 100
df <- data.frame()
for(i in 1:nsim){
  samp <- sample(x,10)
  df <- rbind(df,data.frame(
    sim=paste0("SIM",i),
    G1=samp[1],
    G2=samp[2],
    G3=samp[3],
    G4=samp[4],
    G5=samp[5],
    G6=samp[6],
    G7=samp[7],
    G8=samp[8],
    G9=samp[9],
    G10=samp[10]))
}

write.table(df,
            file = "synthetic_metagenomes.tbl",
            row.names = F,
            col.names = F,
            quote = F)