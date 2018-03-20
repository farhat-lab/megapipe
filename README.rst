Workflow
========

Installation
############

Before being able to use megapipe, you need to install the following softwares on your computer (or use the correponding modules if you are working on o2):

 * bwa (we use it to map the reads to the reference genome)
 * prinseq (we use it to trim the reads)
 * picard (we use it to remove the duplicated reads)
 * pilon (we use it to determine the differences between the query genome and the reference genome and generate .vcf files)
 * spades (we use it to generate an assembly of the query genome from the trimmed reads)
 * samtools (we use it to calculate the depth of our sequencing data and other file format conversions)
 * vrtTools (it is a requirement of the Fast lineage caller) -- https://github.com/farhat-lab/vrt-tools
 * Fast lineage caller (we use it to do the lineage calls) -- https://github.com/farhat-lab/fast-lineage-caller

You need to get the tar.gz present in the dist directory. Then you can install the package: 
::
 pip install <path_to_targz>

Now you need to create and set up the .megapipe.json configuration file

 * This file is needed to tell megapipe where some of the programs used by the pipeline are located (these programs are .jar files that in theory can be located everywhere, so it is difficult to autodetect them). 
 * This file is a hidden file located in your home directory (~/.megapipe.json). You will able to list it if you use the -a option with the command 'ls'.

Here are the steps to create and edit this file. 
::
    cd
    > .megapipe.json
    nano .megapipe.json 

Now you can paste the following code (json format):
::
 {
 "kraken_db":"/n/data1/hms/dbmi/farhat/bin/kraken/tbdb/",
 "picard":"/n/app/picard/2.8.0/bin/picard-2.8.0.jar",
 "pilon":"/n/data1/hms/dbmi/farhat/bin/pilon/pilon-1.22.jar",
 "qualimap":"/home/lf61/mfarhat/bin/qualimap_v2.2.1/qualimap",
 "prinseq":"/home/lf61/mfarhat/bin/prinseq-lite-0.20.4/prinseq-lite.pl",
 "spades":"/home/lf61/sw/spades/3.10.1/bin/spades.py",
 "fasta_ref":"/home/lf61/lf61/repos/megapipe/RefGen/TBRefGen.fasta",
 "lineage_snp_db":"/home/lf61/lf61/repos/fast-lineage-caller/example/db_snps.tsv"
 }


Of course you need to change the paths of each program with the paths. Save the file and exit.

Downloading metadata of public strains
##################################
In order to download the metadata from runs that are available on NCBI you can run the follwing command:
::
 #megapipe-download-metadata-from-ncbi.py <txt_with_run_ids> <tag_dest_files>
 megapipe-download-metadata-from-ncbi.py toDownload.txt dataNCBI

Two files will be generated: one (<tag>.txt with the raw text data and one <tag>_tab.txt with the data in tabular format).

Note: on o2 you may need to load the perl module in order to use this script since eutils, the tools that retrieve the data from NCBI are written in operl. 
::
 module load gcc/6.2.0 perl/5.24.0

In a future version of megapipe we will evaluate the possibility of using a python module instead of the standalone version of eutils.

Create or modify the table for strain identification
###############################################
In order to create a new table for strain identification, you can run the following command:
::
 #megapipe-create-table-identification-strains.py <table_metadata> <table_identification_strains>
 megapipe-create-table-identification-strains.py dataNCBI_tab.txt dataNCBI_table_identification_strains.txt

Notes: 

 * if you create a brand new table, please start tracking the changes with git. So that if something goes wrong you have the chance to go back.
 * you are supposed to create a this table starting from public data. If you want to start from your own data, please change this script.

In order to add new strains to an existing table, you can run the following command:
::
 #megapipe-modify-table-identification-strains.py <table_identification_strains> <table_metadata>
 megapipe-modify-table-identification-strains.py dataNCBI_table_identification_strains.txt new_metadata.txt
Note: I am adding again public data.

In order to add internal strains to the table, you can use the same command:
::
 #megapipe-modify-table-identification-strains.py <table_identification_strains> <table_metadata>
 megapipe-modify-table-identification-strains.py dataNCBI_table_identification_strains.txt new_metadata2.txt

However, plese take into account that internal strains MUST have a public_xref set to "" and MUST have a column "internal_fastq_files" that tells megapipe where to retrive the fastq files. Here is an example of a table for internal strains:
::
 internal_xref   internal_fastq_files
 01-R0902        run1:/home/lf61/mfarhat/fastq_db/pools/01-R0902.1.fastq.gz,/home/lf61/mfarhat/fastq_db/pools/01-R0902.2.fastq.gz

Each sequencing run included into "internal_fastq_files" should have the following format:
::
 <run_name>:<fastq1>,<fastq2>
If there are multiple runs, the synthax becomes the following:
::
 <run_nameA>:<fastq1>,<fastq2>;<run_nameB>:<fastq1>,<fastq2>

Downloading data for public strains (NCBI)
######################################
Retrieving ids of the runs for the public strains:
::
 megapipe-retrieve-runIDs-from-table.py <table> <dir_results> <file_output>
 megapipe-retrieve-runIDs-from-table.py dataNCBI_table_identification_strains.txt results/ runsToDownload.txt

Notes:

 * if you do not have a directory with some results, just create a new directory
 * the script checks the <dir_results> to see if you already analyzed some of the strains. If there is a directory that matches the public_xref of one of the strains, the script will not put the corresponding runs into the output file
 * are you worried about the internal strains? You should have already set the internal_fastq_files variable for these runs, right (see above)? If that's the case, you are all set!

In order to download fastq files from NCBI you can use two utilities:

 * megapipe-download-fastq-from-ncbi.py
 * megapipe-download-fastq-from-ncbi-HT-o2.py

Use "megapipe-download-fastq-from-ncbi.py" when you have a few fastq files to download (5 or less) or you need to dowload the runs sequentially (num_of_threads=1 in this case). First you need to have a text file with the run ids you want to download. For instance:
::
 SRR023455
 SRR023480
 SRR026444

In order to download the runs, open an interactive session and choose the number of cores you need and the amount of memory (10G should be fine):
::
 srun -n 3 -t 0-6:00 --pty -p interactive --mem=10G /bin/bash

Then run the script:
::
 # synthax: megapipe-download-fastq-from-ncbi.py <txt_file_with_run_ids> <dest_directory> <num_of_threads>
 megapipe-download-fastq-from-ncbi.py runsToDownload.txt fastq 3

Note: it takes 45m to download three runs. 

Use "megapipe-download-fastq-from-ncbi-HT-o2.py" if you need to download quickly multiple sequencing runs from NCBI.
First you need to have a text file with the run ids you want to download. For instance:
::
 SRR023455
 SRR023480
 SRR026444

Now you can run the script:
::
 # synthax: megapipe-download-fastq-from-ncbi-HT-o2.py <txt_file_with_run_ids> <dest_directory> <directory_log_files>
 megapipe-download-fastq-from-ncbi-HT-o2.py runsToDownload.txt /n/scratch2/lf61/fastq logs

In order to check if the download finished or not, please use the "squeue" command:
::
 #squeue|grep <your_username>
 squeue|grep lf61

Generating all genomic data
#########################

Create a directory where you want to store your data (if you did not do it before)
::
 mkdir results
 cd results

Run the pipeline
Here is the general synthax  of the command:
::
 megapipe-launch.py <table_identification_strains> <fastq_dir> <output_dir> <scratch_dir> <jobs_to_launch>  

For instance here is an example that show how to launch a megapipe analysis for all the genomes of a dataset of the RESEQTB project:


Note: 
If you are running megapipe on o2, megapipe has to use some modules in order to work properly. Please create a loadmodules.txt file in the directory where you will execute megapipe-launch.py:
::
 gcc/6.2.0|perl/5.24.0|picard/2.8.0


**GOOD LUCK for your analyses!**

**NOTE: remember to clean the scratch from time to time!** 

Version history
===============


Todo
==== 
* v2.0
    * parse the vcf file so that we have a shorter version of it.
    * add a directory to store the pilon results
    * everything should happen in the scratch. Just save the final results on the results directory
    * add lineage calling
        * modify vrtTools so that they work with python 3
    * log the versions of the programs that megapipe uses (important when we want to write papers)
    * improve the output that goes into the grid engine output file

Misc
====

How to deal with pip
##################
How to pack the module:
::
 python setup.py sdist

How to install the module:
::
 pip install megapipe-0.1.0.tar.gz

How to remove the module:
::
 pip uninstall megapipe

How to use the gridmanager module
##############################
Here is an example:
::
 from gridmanager import gridpuppeteer as gp
 a=gp.GridEngine()
 a.generate_script("prova.sh","short","12:00","prova.out","10M","wget http://poisson.phc.unipi.it/~freschi/img/luca.jpg")
 a.launch_job("prova.sh")


