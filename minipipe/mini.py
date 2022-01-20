# _______
# This is a snakemake file/script for processing of Paired End Illumina WGS generated from Mtb.
### Maximillian Marin (mgmarin@g.harvard.edu)
# Define PATH to the reference genome to be used:
refGenome_FA_PATH = config[“RefGenome_FA_PATH”]
#refGenome_GFF_PATH = config[“RefGenome_GFF_PATH”]
# Kraken2_DB_PATH = config[“Kraken2_DB_PATH”]
# Define PATH of OUTPUT directory
output_Dir = config[“output_dir”]
################################
import pandas as pd
# Read in data regarding input
sampleInfo_DF = pd.read_csv( config[“inputSampleData_TSV”], sep=‘\t’)
# Save a list of SRA/ENA “Run” Accessions
input_SampleIDs_SRA_RunAcc = list( sampleInfo_DF[“SRA_RunID”].values )
input_SampleIDs_WithIllumina = input_SampleIDs_SRA_RunAcc
##### Import Snakemake rules from SMK files #####
# include: “./2_Illumina_Mtb_Mtb_Analysis.smk”
rule all:
    input:
        expand(output_Dir + “/{sampleID_WiIll}/Minimap2_Ill_SPAdesAssembly_AlignTo_H37rv/{sampleID_WiIll}_mm2_Ill_SPAdes_AssemblyToH37rv.vcf”, sampleID_WiIll = input_SampleIDs_WithIllumina),
rule download_FQ_FromSRA_RunID:
    output:
         fq1_unzipped = output_Dir + “/{sampleID_WiIll}/FASTQs/{sampleID_WiIll}_1.fastq”,
         fq2_unzipped = output_Dir + “/{sampleID_WiIll}/FASTQs/{sampleID_WiIll}_2.fastq”,
    params:
        target_DownloadDir = output_Dir + “/{sampleID_WiIll}/FASTQs/”
    conda:
        “CondaEnvs/sratools_2_10_7_Conda.yml”
    shell:
        “fastq-dump --split-files --outdir {params.target_DownloadDir} {wildcards.sampleID_WiIll}\n”
# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_PE_Trimming:
    input:
        r1 = output_Dir + “/{sampleID_WiIll}/FASTQs/{sampleID_WiIll}_1.fastq”,
        r2 = output_Dir + “/{sampleID_WiIll}/FASTQs/{sampleID_WiIll}_2.fastq”,
    output:
        r1 = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq”,
        r2 = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq”,
        # reads where trimming entirely removed the mate
        r1_unpaired = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.unpaired.fastq”,
        r2_unpaired = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.unpaired.fastq”,
    log:
        output_Dir + “/logs/trimmomatic/{sampleID_WiIll}.log”
    params:
        # list of trimmers (see manual)
        trimmer=[“ILLUMINACLIP:‘References/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta’:2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        # optional parameters
        # extra=” ”
    threads: 8
    wrapper:
        “0.38.0/bio/trimmomatic/pe”
# Adding SPADes assembly with UNIcycler to analysis
### A) Assembly with SPAdes through Unicycler
rule unicycler_SPAdes_Assemble_IlluminaWGS:
    input:
        fq1_trimmed = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq”,
        fq2_trimmed = output_Dir + “/{sampleID_WiIll}/FASTQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq”,
        DnaA_Seq_fa = “References/DnaA_MTb_H37Rv_dna.fasta”   #H37rv_DnaA_FA_PATH,
    output:
        assembly_GFA = output_Dir + “/{sampleID_WiIll}/Unicycler_SPAdesAssembly/assembly.gfa”,
        assembly_fa = output_Dir + “/{sampleID_WiIll}/Unicycler_SPAdesAssembly/assembly.fasta”,
    threads: 8
    params:
        Unicycler_OutputDir_PATH = output_Dir + “/{sampleID_WiIll}/Unicycler_SPAdesAssembly/”
    shell:
        “unicycler --vcf -t {threads}  --start_genes {input.DnaA_Seq_fa} ”
        ” -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} ”
        ” -o {params.Unicycler_OutputDir_PATH} ”
        # ” ‑‑mode conservative ” # This version of unicycler doesn’t support the --mode arguement
### B) Call variants from assemblies of isolates using Minimap2 (To H37rv)
rule Minimap2_IlluminaWGS_SPAdesAssembly_AlignTo_H37rv:
    input:
        Ill_SPAdes_Assembly_fa = output_Dir + “/{sampleID_WiIll}/Unicycler_SPAdesAssembly/assembly.fasta”,
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_Ill_SPAdes_To_H37rv_SAM = output_Dir + “/{sampleID_WiIll}/Minimap2_Ill_SPAdesAssembly_AlignTo_H37rv/{sampleID_WiIll}_mm2_Ill_SPAdes_AssemblyToH37rv.sam”,
        MM2_Ill_SPAdes_To_H37rv_BAM = output_Dir + “/{sampleID_WiIll}/Minimap2_Ill_SPAdesAssembly_AlignTo_H37rv/{sampleID_WiIll}_mm2_Ill_SPAdes_AssemblyToH37rv.bam”,
        MM2_Ill_SPAdes_To_H37rv_bai = output_Dir + “/{sampleID_WiIll}/Minimap2_Ill_SPAdesAssembly_AlignTo_H37rv/{sampleID_WiIll}_mm2_Ill_SPAdes_AssemblyToH37rv.bam.bai”,
        MM2_Ill_SPAdes_To_H37rv_VCF = output_Dir + “/{sampleID_WiIll}/Minimap2_Ill_SPAdesAssembly_AlignTo_H37rv/{sampleID_WiIll}_mm2_Ill_SPAdes_AssemblyToH37rv.vcf”,
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
    shell:
        “minimap2 -ax asm10 --cs {input.H37rv_FA} {input.Ill_SPAdes_Assembly_fa} > {output.MM2_Ill_SPAdes_To_H37rv_SAM} \n”
        “samtools view -bS {output.MM2_Ill_SPAdes_To_H37rv_SAM} | samtools sort - > {output.MM2_Ill_SPAdes_To_H37rv_BAM} \n”
        “samtools index {output.MM2_Ill_SPAdes_To_H37rv_BAM} \n”
        “minimap2 -cx asm10 --cs {input.H37rv_FA} {input.Ill_SPAdes_Assembly_fa} | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIll} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_Ill_SPAdes_To_H37rv_VCF} \n”