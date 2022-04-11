#!/bin/sh

cd /media/ink/zou/Assemblies/

Rscript ../genome_sample_for_synthetic.R 

while read SIMNAME G1 G2 G3 G4 G5 G6 G7 G8 G9 G10; do
  iss generate --n_reads 100M  --cpus 60 --model hiseq --output $SIMNAME --coverage lognormal  --draft $G1 $G2 $G3 $G4 $G5 $G6 $G7 $G8 $G9 $G10
done <synthetic_metagenomes.tbl
