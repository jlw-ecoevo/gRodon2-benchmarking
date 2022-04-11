#!/bin/sh

cd /media/ink/zou/Assemblies/

while read SIMNAME G1 G2 G3 G4 G5 G6 G7 G8 G9 G10; do
  #if [ ! -d "${SIMNAME}_assembled" ]; then
  #  ~/fastp -w 16 -i ${SIMNAME}_R1.fastq -I ${SIMNAME}_R2.fastq -o ${SIMNAME}_R1.clean.fastq -O ${SIMNAME}_R2.clean.fastq
  #  spades.py -o ${SIMNAME}_assembled --meta --pe1-1 ${SIMNAME}_R1.clean.fastq --pe1-2 ${SIMNAME}_R2.clean.fastq --threads 30
  #fi
  seqtk seq -A ${SIMNAME}_R1.clean.fastq ${SIMNAME}_R2.clean.fastq > ${SIMNAME}_reads.fasta
  FragGeneScanRs --seq-file-name ${SIMNAME}_reads.fasta --thread-num 10 --training-file illumina_5 --output-prefix ${SIMNAME}_fraggenescan
  #blastn -query ${SIMNAME}_fraggenescan.ffn -db ~/riboprot_db/riboprot_growthpred_data -out ${SIMNAME}_ribo -outfmt 6 -num_threads 10
done <synthetic_metagenomes.tbl
