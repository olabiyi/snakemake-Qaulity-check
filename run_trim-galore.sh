#!/usr/bin/env bash

 source activate /home/jeffbrady/anaconda3/envs/bioinfo

SAMPLES=($(ls -1 01.raw_data/ | grep -Ev "MANIFEST|seq" - |sort -V))

parallel -j 10 \
 '[ -d 06.cutadapt_trim/{}/ ] || mkdir -p 06.cutadapt_trim/{}/ && trim_galore --paired -o 06.cutadapt_trim/{}/ --fastqc  04.Trim_primers/{}/{}_R1.fastq.gz 04.Trim_primers/{}/{}_R2.fastq.gz' \
  ::: ${SAMPLES[*]} && \
  multiqc --interactive -f 06.cutadapt_trim/ -o 06.cutadapt_trim/
