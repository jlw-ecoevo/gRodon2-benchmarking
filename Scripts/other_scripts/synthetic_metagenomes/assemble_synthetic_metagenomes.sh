#!/bin/sh

cd /media/ink/zou/Assemblies/

while read SIMNAME G1 G2 G3 G4 G5 G6 G7 G8 G9 G10; do
  if [ ! -d "${SIMNAME}_assembled" ]; then
    ~/fastp -A -w 16 -i ${SIMNAME}_R1.fastq.gz -I ${SIMNAME}_R2.fastq.gz -o ${SIMNAME}_R1.clean.fastq -O ${SIMNAME}_R2.clean.fastq
    #spades.py -o ${SIMNAME}_assembled --meta --pe1-1 ${SIMNAME}_R1.clean.fastq.gz --pe1-2 ${SIMNAME}_R2.clean.fastq.gz --threads 60
    megahit -1 ${SIMNAME}_R1.clean.fastq -2 ${SIMNAME}_R2.clean.fastq -o ${SIMNAME}_megahit -t 60
  fi
done <synthetic_metagenomes.tbl
