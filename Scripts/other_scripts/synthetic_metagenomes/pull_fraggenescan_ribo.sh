#!/bin/sh

cd /media/ink/zou/Assemblies

ls *_ribo > ribo.files

while read FILE; do
  awk '$11<1e-5' ${FILE} | awk '$3>99' | awk '{print $1}' > ${FILE}.genes
done <ribo.files
