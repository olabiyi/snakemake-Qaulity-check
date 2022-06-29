#!/usr/bin/env bash


module purge; module load Anaconda3/2020.07; source activate /scratch/user/obayomi/.conda/envs/bioinfo



parallel -j 10 \
 '[ -d 01.trim/{}/ ] || mkdir 01.trim/{}/ && trim_galore --paired -o 01.trim/{}/ --fastqc  ../01.raw_data/{}/{}_R1.fastq.gz ../01.raw_data/{}/{}_R2.fastq.gz' \
  ::: ${SAMPLES[*]} && \
  multiqc --interactive -f 01.trim/ -o 02.QC/
