#!/usr/bin/env python
configfile: "./config/config.json"

"""
This will be a snakemake file to deal with the download operation
"""

onstart: 
    "python3 ./bin/get_list_run_ids.py"

with open(config["logs_analysis"]+"runs_to_download.txt", "r") as inp:
    samples = inp.read().splitlines() 

rule all:
    input: 
        expand(config["fastq_dir"]+"{sample}_1.fastq.gz", sample=samples),
        expand(config["fastq_dir"]+"{sample}_2.fastq.gz", sample=samples)

rule download_fastq:
    output: config["fastq_dir"]+"{sample}_1.fastq.gz", config["fastq_dir"]+"{sample}_2.fastq.gz"
    shell:
        """
        SRAPATH=`srapath {wildcards.sample}`
        wget -O {config[fastq_dir]}/{wildcards.sample} "$SRAPATH"
        fastq-dump --split-files --gzip {config[fastq_dir]}/{wildcards.sample} -O {config[fastq_dir]}
        """

