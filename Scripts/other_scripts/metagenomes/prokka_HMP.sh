#!/bin/sh

cd /media/ink/HMP/Assemblies

find -type f -name "*.fa" | awk 'gsub("./","")' > gene_list.txt

while read FILE; do
  NAME=${FILE%.scaffolds.fa}
  prokka ${FILE} --norrna --notrna --metagenome --cpus 40 --centre X --compliant --prefix $NAME
done <gene_list.txt
