#!/bin/sh

cd /media/ink/zou/Assemblies/

while read SIMNAME G1 G2 G3 G4 G5 G6 G7 G8 G9 G10; do
  bwa index ${SIMNAME}_prokka/${SIMNAME}_prokka.ffn
  bwa mem -t 60 ${SIMNAME}_prokka/${SIMNAME}_prokka.ffn ${SIMNAME}_R1.clean.fastq ${SIMNAME}_R1.clean.fastq | samtools view -bS -  > ${SIMNAME}.bam
  samtools sort ${SIMNAME}.bam -o ${SIMNAME}.sort.bam
  ~/bamcov/bamcov -H ${SIMNAME}.sort.bam > ${SIMNAME}.cov
done <synthetic_metagenomes.tbl
