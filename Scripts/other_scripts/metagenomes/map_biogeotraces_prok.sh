#/bin/sh

cd /media/ink/BIOGEOTRACES/

mkdir mapped

while read SAMP SRR SRZ; do
  echo ${SRR}
  echo ${SRZ}
  echo ${SAMP}
  fasterq-dump sra/${SRR}.sra -O /media/ink/BIOGEOTRACES/sra/
  ~/fastp -w 16 -i sra/${SRR}_1.fastq -I sra/${SRR}_2.fastq -o sra/${SRR}_1.clean.fastq -O sra/${SRR}_2.clean.fastq
  bwa index prok_annot/${SRZ}_${SAMP}.ffn
  bwa mem -t 16 prok_annot/${SRZ}_${SAMP}.ffn sra/${SRR}_1.clean.fastq sra/${SRR}_2.clean.fastq | samtools view -bS - | samtools sort -  > mapped/${SRZ}_${SAMP}.bam
  ~/bamcov/bamcov -H mapped/${SRZ}_${SAMP}.bam > mapped/${SRZ}_${SAMP}.cov
  rm sra/${SRR}_1.fastq
  rm sra/${SRR}_2.fastq
  rm sra/${SRR}_1.clean.fastq
  rm sra/${SRR}_2.clean.fastq
  rm mapped/${SRZ}_${SAMP}.bam
done <BIOGEOTRACES_srr_to_srz_todo.tab
