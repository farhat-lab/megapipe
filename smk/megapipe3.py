configfile: "./config/config_pipeline.json"

#I read the file with the of the samples to analyze to get the list of samples
with open(config["logs_analysis"]+"isolates_to_analyze.txt", "r") as inp:
    samples = inp.read().splitlines()
print(samples)
#print(expand(config["temp_dir"]+ "{sample}/{sample}-combined_1.fastq", sample=samples))

rule all:

    input:
        config["path_ref_genome"]+".bwt",
        expand("results/{sample}/bam/{sample}.duprem.bam.bai", sample=samples),
        expand("results/{sample}/fast-lineage-caller/{sample}.lineage", sample=samples)

rule generate_idx_ref:

    input: config["path_ref_genome"]

    output: config["path_ref_genome"]+".bwt"

    log: config["logs_analysis"]+"indexed_ref.txt"

    shell: "bwa index {input} 2> {log}"

rule trim_with_trimmomatic:

    input:
        read1 = config["fastq_dir"] + "{sample}_1.fastq.gz",
        read2 = config["fastq_dir"] + "{sample}_2.fastq.gz",

    output:
        trimmed1 = config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq",
        trimmed2 = confid["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired = config["temp_dir"] + "{sample}/{sample}_1_trimmed.unpaired.fastq",
        r2_unpaired = config["temp_dir"] + "{sample}/{sample}_2_trimmed.unpaired.fastq",

    log:
        config["logs_analysis"]+"{sample}/trimmomatic.txt"

    params:
        # list of trimmers (see manual) http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
        trimmer=[“ILLUMINACLIP:‘References/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta’:2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        # optional parameters
        # extra=” ”

    threads: 8

    wrapper:
        “0.38.0/bio/trimmomatic/pe”

rule assemble_with_spades:

    input:
        trimmed1 = config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq",
        trimmed2 = confid["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq",
    
    output:
        assembled1 = config["temp_dir"]+ "{sample}/{sample}-assembled_1.fastq",
        assembled2 = config["temp_dir"]+ "{sample}/{sample}-assembled_2.fastq"
    
    log:
        config["logs_analysis"]+"{sample}/spades.txt"

    threads: 8

    resources:
        mem_mem=250000,
        time=60 * 24,

    wrapper:
        "v0.86.0/bio/spades/metaspades"
    

rule align_to_ref:

    input:
        assembled1 = config["temp_dir"]+ "{sample}/{sample}-assembled_1.fastq",
        assembled2 = config["temp_dir"]+ "{sample}/{sample}-assembled_2.fastq"

    output:
        temp(config["temp_dir"]+ "{sample}/{sample}.sam")

    log:
        config["logs_analysis"]+"{sample}/align_to_ref.txt"

    shell:
        """
            minimap2 -ax sr {config[path_ref_genome]} {input.comb_fq1} {input.comb_fq2} > {output} 2> {log}
        """

rule sort_convert_tobam:

    input:
        config["temp_dir"]+ "{sample}/{sample}.sam"

    output:
        temp(config["temp_dir"]+ "{sample}/{sample}.bam")

    log:
        config["logs_analysis"]+"{sample}/sort_convert_tobam.txt"

    shell:
        """
        java -Xmx16G -jar {config[path_to_picard]} SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate > {log} 2>&1
        """

rule remove_duplicates:

    input:
        config["temp_dir"]+ "{sample}/{sample}.bam"

    output:
        outfile="results/{sample}/bam/{sample}.duprem.bam", metrics=config["temp_dir"]+"results/{sample}/bam/{sample}.metrics"

    log:
        config["logs_analysis"]+"{sample}/duprem.txt"

    shell:
        "java -Xmx32G -jar {config[path_to_picard]} MarkDuplicates I={input} O={output.outfile} REMOVE_DUPLICATES=true M={output.metrics} ASSUME_SORT_ORDER=coordinate > {log} 2>&1"

rule calculate_depth:

    input:
        "results/{sample}/bam/{sample}.duprem.bam"

    output:
        "results/{sample}/depth/{sample}.depth.gz", "results/{sample}/depth/{sample}_depth_OK"

    log:
        config["logs_analysis"]+"{sample}/calc_depth.txt"

    shell:
        """
        samtools depth -a {input} > results/{wildcards.sample}/depth/{wildcards.sample}.depth 2> {log}
        ./bin/check_depth.py results/{wildcards.sample}/depth/{wildcards.sample}.depth >> {log} 2>&1
        gzip results/{wildcards.sample}/depth/{wildcards.sample}.depth >> {log} 2>&1
        """

rule indexing_bam:

    input:
        bam="results/{sample}/bam/{sample}.duprem.bam", depth_ok="results/{sample}/depth/{sample}_depth_OK"

    output:
        "results/{sample}/bam/{sample}.duprem.bam.bai"

    log:
        config["logs_analysis"]+"{sample}/indexing_bam.txt"

    shell:
        "samtools index {input.bam} > {log} 2>&1"

rule variant_calling:

    input:
        bam="results/{sample}/bam/{sample}.duprem.bam", depth_ok="results/{sample}/depth/{sample}_depth_OK", bai="results/{sample}/bam/{sample}.duprem.bam.bai"

    output:
        "results/{sample}/pilon/{sample}.vcf", "results/{sample}/pilon/{sample}.fasta", "results/{sample}/pilon/{sample}_full.vcf.gz"

    log:
        config["logs_analysis"]+"{sample}/variant_calling.txt"

    shell:
        """
        pilon -Xmx32G --genome {config[path_ref_genome]} --bam {input.bam} --output results/{wildcards.sample}/pilon/{wildcards.sample} --variant > {log} 2>&1
        mv results/{wildcards.sample}/pilon/{wildcards.sample}.vcf results/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf >> {log} 2>&1
        ./bin/vcf_cutter.py results/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf results/{wildcards.sample}/pilon/{wildcards.sample}.vcf >> {log} 2>&1
        gzip results/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf >> {log} 2>&1
        """

rule lineage_calling:

    input:
        "results/{sample}/pilon/{sample}.vcf"

    output:
        "results/{sample}/fast-lineage-caller/{sample}.lineage"

    log:
        config["logs_analysis"]+"{sample}/lineage_calling.txt"
        #Ideally could have everything we install into SW put into the megapipe repo and downloaded concurrently
    shell:
        "fast-lineage-caller-vcf {input} /home/maf7596/sw/fast-lineage-caller-vcf/snp_schemes/coll.tsv > {output} 2> {log}"

rule rem_temp_files:

    shell:
        """
        cd {config[temp_dir]}
        printf "current directory: ";pwd
        rm -rf -I *
        """

rule rem_results:

    shell:
        """
        cd results
        printf "current directory: ";pwd
        rm -rf -I *
        """
