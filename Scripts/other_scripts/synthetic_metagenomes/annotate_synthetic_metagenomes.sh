#!/bin/sh

cd /media/ink/zou/Assemblies/

while read SIMNAME G1 G2 G3 G4 G5 G6 G7 G8 G9 G10; do
  prokka ${SIMNAME}_megahit/final.contigs.fa --norrna --notrna --metagenome --cpus 3 --centre X --compliant --prefix ${SIMNAME}_prokka
done <synthetic_metagenomes.tbl
