#/bin/sh

cd /media/ink/HMP/Reads/

mkdir ../mapped

while read SRS; do
  echo ${SRS}
  tar -xvf ${SRS}.tar.bz2
  ~/fastp -w 16 -i ${SRS}/${SRS}.denovo_duplicates_marked.trimmed.1.fastq -I ${SRS}/${SRS}.denovo_duplicates_marked.trimmed.1.fastq -o ${SRS}/${SRS}.clean.1.fastq -O ${SRS}/${SRS}.clean.2.fastq
  bwa index ../CDS/${SRS}.ffn
  bwa mem -t 16 ../CDS/${SRS}.ffn ${SRS}/${SRS}.clean.1.fastq ${SRS}/${SRS}.clean.2.fastq | samtools view -bS - | samtools sort -  > ../mapped/${SRS}.bam
  ~/bamcov/bamcov -H ../mapped/${SRS}.bam > ../mapped/${SRS}.cov
  rm ${SRS}/${SRS}.clean.1.fastq 
  rm ${SRS}/${SRS}.clean.2.fastq 
  rm ../mapped/${SRA}.bam
  rm -r ${SRS}/
done <SRS_todo.txt
