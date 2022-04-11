#!/bin/sh

cd /media/ink/zou/Assemblies/

ls *_fraggenescan.ffn > fgs.files

while read FILE; do
  seqkit sample -p 0.1 -o ${FILE%.ffn}_sample10.ffn ${FILE}
done <fgs.files
