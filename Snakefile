from os import path, getcwd

# Run the pipeline like so:
# snakemake -pr --cores 10 --keep-going --rerun-incomplete --restart-times 3
# snakemake -s Snakefile --rulegraph |dot -Tpng > rulegraph.png
configfile: "config/config.yaml"

onsuccess:
    print("Workflow completed without any error")

onerror:
    print("An error occurred")


RULES=["Count_seqs_pre_trim", "QC_pre_trim", "SummarizeQC_pre_trim", "Trimmomatic_trim_adaptors",
       "Trim_primers", "Trim_galore_trim_adaptors", "SummarizeQC_post_adaptor_trim",
       "Count_seqs_post_trim", "Copy_clean_read", "Copy_seq_reports"]

rule all:
    input:
        "02.Count_seqs_pre_trim/seqs_stat.txt",
        "04.SummarizeQC_pre_trim/multiqc_report.html",
        "08.SummarizeQC_post_adaptor_trim/multiqc_report.html",
        "09.Count_seqs_post_trim/seqs_stat.txt",
         expand(["10.Export/{sample}/{sample}_R1.fastq.gz", "01.raw_data/{sample}/{sample}_R2.fastq.gz"], sample=config['samples']),
         "10.Export/post_trim_multiqc_report.html"


rule Make_logs_directories:
    output:
        directory("logs/Trim_galore_trim_adaptors/"),
        directory("logs/QC_pre_trim/"),
        directory("logs/SummarizeQC_post_adaptor_trim/"),
        directory("logs/Trim_primers/"),
        directory("logs/SummarizeQC_pre_trim/"),
        directory("logs/Trimmomatic_trim_adaptors/")
    threads: 1
    shell:
        """
         [ -d logs/ ] || mkdir -p logs/
         cd logs/
         for RULE in {RULES}; do
          [ -d ${{RULE}}/ ] || mkdir -p ${{RULE}}/
         done
        """


rule Count_seqs_pre_trim:
    input: expand(["01.raw_data/{sample}/{sample}_R1.fastq.gz", "01.raw_data/{sample}/{sample}_R2.fastq.gz"], sample=config['samples'])
    output: "02.Count_seqs_pre_trim/seqs_stat.txt"
    params:
        program=config['programs_path']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        {params.program} stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """



rule QC_pre_trim:
    input:
        rules.Make_logs_directories.output,
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz"
    output:
        forward_html="03.QC_pre_trim/{sample}/{sample}_R1_fastqc.html",
        reverse_html="03.QC_pre_trim/{sample}/{sample}_R2_fastqc.html"
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=1
    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 1
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {params.threads} {input.forward} {input.rev} > {log} 2>&1"


rule SummarizeQC_pre_trim:
    input:
        forward_html=expand("03.QC_pre_trim/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("03.QC_pre_trim/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "04.SummarizeQC_pre_trim/multiqc_report.html"
    log: "logs/SummarizeQC_pre_trim/multiqc.log"
    params:
        program=config['programs_path']['multiqc'],
        in_dir=lambda w, input: path.dirname(input[0]).split("/")[0],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} \
              --interactive \
              -f {params.in_dir} \
              -o {params.out_dir}  > {log} 2>&1
        """



adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']
rule Trimmomatic_trim_adaptors:
    input:
        rules.Make_logs_directories.output,
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz",
        log_dirs=rules.Make_logs_directories.output
    output:
        r1="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        r2="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.unpaired.fastq.gz"
    log:
        "logs/Trimmomatic_trim_adaptors/{sample}/{sample}.log"
    params:
        program=config['programs_path']['trimmomatic'],
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length)
    threads: 1
    resources:
        mem_mb=1024
    shell:
        "{params.program} PE "
        "-threads {threads} "
        "{input.forward} {input.rev} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} > {log} 2>&1 "




# Trim Primers using cutadapt
rule Trim_primers:
    input:
        rules.Make_logs_directories.output,
        forward_reads=rules.Trimmomatic_trim_adaptors.output.r1,
        rev_reads=rules.Trimmomatic_trim_adaptors.output.r2
    output:
        forward_reads="06.Trim_primers/{sample}/{sample}_R1.fastq.gz",
        rev_reads="06.Trim_primers/{sample}/{sample}_R2.fastq.gz"
    log: "logs/Trim_primers/{sample}/{sample}.log"
    params:
        program=config['programs_path']['cutadapt'],
        forward_primer=config['parameters']['cutadapt']['forward_primer'],
        rev_primer=config['parameters']['cutadapt']['reverse_primer'],
        minimum_length=config['parameters']['cutadapt']['minimum_length'],
        quality_cutoff=config['parameters']['cutadapt']['quality_cutoff'],
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        # -a TGGAATTCTCGGGTGCCAAGG sequence here is the small RNA 3' adaptor
        # https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/RNA/Small-RNA/TruSeqSmallRNA.htm
        # -A CTGTCTCTTATACAC sequece here is the nextera transposa sequence 
        {params.program} \
              -g '{params.forward_primer}' \
              -G '{params.rev_primer}' \
              -o {output.forward_reads} \
              -p {output.rev_reads} \
              --minimum-length  {params.minimum_length} \
              --quality-cutoff  {params.quality_cutoff} \
             {input.forward_reads} {input.rev_reads} > {log} 2>&1

        """

# Trim adaptors and quality check using cutadapt and factqc with Trim galore
rule Trim_galore_trim_adaptors:
    input: 
        rules.Make_logs_directories.output,
        forward=rules.Trim_primers.output.forward_reads,
        rev=rules.Trim_primers.output.rev_reads
    output:
        forward_reads="07.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        rev_reads="07.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        forward_html="07.Trim_galore_trim_adaptors/{sample}/{sample}_R1_fastqc.html",
        rev_html="07.Trim_galore_trim_adaptors/{sample}/{sample}_R2_fastqc.html"
    log: "logs/Trim_galore_trim_adaptors/{sample}/{sample}.log"
    threads: 1
    params:
        program=config['programs_path']['trim_galore'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    shell:
        """ 
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
           -o {params.out_dir} \
           --fastqc  \
           --paired {input.forward} {input.rev} > {log} 2>&1

         #Rename the files
         #Fastq files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1.fq.gz {params.out_dir}/{wildcards.sample}_R1.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_val_2.fq.gz {params.out_dir}/{wildcards.sample}_R2.fastq.gz
         #HTML files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.html {params.out_dir}/{wildcards.sample}_R1_fastqc.html
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.html {params.out_dir}/{wildcards.sample}_R2_fastqc.html
         #Zip files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.zip {params.out_dir}/{wildcards.sample}_R1_fastqc.zip
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.zip {params.out_dir}/{wildcards.sample}_R2_fastqc.zip
        """


# Aggregate and summarize the quality check using multiqc
rule SummarizeQC_post_adaptor_trim:
    input:
        forward_html=expand("07.Trim_galore_trim_adaptors/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("07.Trim_galore_trim_adaptors/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "08.SummarizeQC_post_adaptor_trim/multiqc_report.html"
    log: "logs/SummarizeQC_post_adaptor_trim/multiqc.log"
    params:
        program=config['programs_path']['multiqc'],
        in_dir=lambda w, input: path.dirname(input[0]).split("/")[0],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} \
              --interactive \
              -f {params.in_dir} \
              -o {params.out_dir}  > {log} 2>&1
        """




rule Count_seqs_post_trim:
    input: expand(["07.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz", "07.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz"], sample=config['samples'])
    output: "09.Count_seqs_post_trim/seqs_stat.txt"
    params:
        program=config['programs_path']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        {params.program} stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """

# Copy the clean reads and the associated pre and post statistics to a
# new folder for ease of exporting
rule Copy_clean_reads:
    input: 
        forward_reads=rules.Trim_galore_trim_adaptors.output.forward_reads,
        rev_reads=rules.Trim_galore_trim_adaptors.output.rev_reads
    output:
        forward_reads="10.Export/{sample}/{sample}_R1.fastq.gz",
        rev_reads="10.Export/{sample}/{sample}_R2.fastq.gz"
    params:
        sample_dir=lambda w, input: path.dirname(input.forward_reads),
        out_dir=lambda w, output: path.dirname(output[0]).split("/")[0]
    shell:
        """
        [ -d  {params.out_dir}/{wildcards.sample}/ ] || mkdir {params.out_dir}/{wildcards.sample}/
        # copy clean reads for export
        cp  {input.forward_reads} {params.out_dir}/{wildcards.sample}/
        cp  {input.rev_reads} {params.out_dir}/{wildcards.sample}/
        """
rule Copy_seq_reports:
    input:
        pre_count=rules.Count_seqs_pre_trim.output,
        post_count=rules.Count_seqs_post_trim.output,
        pre_summary=rules.SummarizeQC_pre_trim.output,
        post_summary=rules.SummarizeQC_post_adaptor_trim.output
    output:
        pre_count="10.Export/pre_trim_seqs_stat.txt",
        pre_summary="10.Export/pre_trim_multiqc_report.html",
        post_count="10.Export/post_trim_seqs_stat.txt",
        post_summary="10.Export/post_trim_multiqc_report.html"
    params:
        out_dir=lambda w, output: path.dirname(output[0])
    shell:
        """
        [ -d {params.out_dir} ] || mkdir {params.out_dir}
        # copy general sequence statistics
        ## Before trimming
        cp {input.pre_count} {params.out_dir}/ && \
          mv {params.out_dir}/seqs_stat.txt {output.pre_count}
        cp {input.pre_summary} {params.out_dir}/ && \
          mv {params.out_dir}/multiqc_report.html {output.pre_summary}

        ## After trimming
        cp {input.post_count} {params.out_dir}/ && \
          mv {params.out_dir}/seqs_stat.txt {output.post_count}
        cp {input.post_summary} {params.out_dir}/ && \
          mv {params.out_dir}/multiqc_report.html {output.post_summary}

        """



