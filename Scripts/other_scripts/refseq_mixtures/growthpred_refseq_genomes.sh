#!/bin/sh

cd /media/ink/pair_metagenomes/

# ----------
# GrowthPred
# ----------

GrowthPred_CONTAINER="GP_G"
GrowthPred_DIR="REFSEQ_GrowthPred_Genomes"
NUCL_DIR="genomes/genomes_to_sample"
RIBO_DIR="genomes/genomes_to_sample"
mkdir -p ${GrowthPred_DIR}

#:<<'COMMENT'
# initialize a container in background
docker run -t -d \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $PWD:/mnt \
           --name ${GrowthPred_CONTAINER} shengwei/growthpred:latest

# without heg
for cds_file in `ls ${NUCL_DIR}/*_genomic.fna`; do
    cds_filename=$(basename $cds_file)
    cds_filestem=${cds_filename%.fna}
    echo ${cds_file}
    echo $RIBO_DIR/${cds_filename}.ribo.fasta


    sudo docker inspect --format="{{.State.Running}}" $GrowthPred_CONTAINER
    if [ $? -eq 0 ];
    then
         echo "existing"
    else
         echo "missing"
         docker run -t -d \
                    -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
                    -u $(id -u $USER):$(id -g $USER) \
                    -v $PWD:/mnt \
                    --name ${GrowthPred_CONTAINER} shengwei/growthpred:latest
    fi

    FILE=${GrowthPred_DIR}/${cds_filestem}_input.fa
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        python ./filter_non_3_nucl_cds.py --out_folder ${GrowthPred_DIR} --prefix ${cds_filestem}_input ${cds_file}
        python ./filter_non_3_nucl_cds.py --out_folder ${GrowthPred_DIR} --prefix ${cds_filestem}_ribo $RIBO_DIR/${cds_filename}.ribo.fasta
    fi


    FILE=${GrowthPred_DIR}/${cds_filestem}.results
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        cmd="/growthpred/growthpred-v1.08.py -d ${GrowthPred_DIR} -g ${cds_filestem}_input.fa -f ${cds_filestem}_ribo.fa \
                       -s -S -c 0 -o ${GrowthPred_DIR}/${cds_filestem} 2>&1 | tee -a ${GrowthPred_DIR}/growthpred.log"
        echo $cmd
        docker exec ${GrowthPred_CONTAINER} /bin/bash -c "$cmd" 2>&1 | tee -a ${GrowthPred_DIR}/growthpred_run.log
    fi

    echo "DONE"
done

# kill container
docker container kill ${GrowthPred_CONTAINER}
#COMMENT

cd ${GrowthPred_DIR}

grep "Predicted minimum generation time:" * | awk 'gsub(".results:Predicted minimum generation time:  ","\t")' | awk 'gsub(" hours.*","")' > ../growthpred_refseq_G.tbl
grep "Nucleotide frequencies:" * | awk 'gsub("Nucleotide frequencies: ","")' | awk '{print $7+$11}' > ../gp_gc_G.txt
paste ../growthpred_refseq_G.tbl ../gp_gc_G.txt > ../growthpred_gc_refseq_genomes_G.tbl

