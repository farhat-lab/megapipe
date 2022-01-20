## Installation
1. Go to the directory where you want to install megapipe and type:

`git clone https://github.com/farhat-lab/megapipe.git`

A directory called `megapipe` will be created!

2. Create a conda environment and activate it. Install the mjority of the software packages with conda

````
conda create --name megapipe3 python=3.8
conda activate megapipe3
conda install snakemake pilon samtools bwa sra-tools picard kraken minimap spades trimmomatic pigz
````

**Note**: remember to check if you have the bioconda channel in your `.condarc` file. Otherwise you can add it and reorder the priority of your channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Install the remaining software:

```
# INSTALLING metatools_ncbi
pip install metatools_ncbi

# INSTALLING fast-lineage-caller
mkdir -p sw/
cd sw/
git clone https://github.com/farhat-lab/fast-lineage-caller-vcf.git
cd fast-lineage-caller-vcf/
python setup.py sdist
pip install dist/fast_lineage_caller_vcf-0.1.1.tar.gz
cd ..

# INSTALL AND CONFIGURE gcc
# This is needed in order to compile fastQValidator
# You might need to modify this a little if you are using MacOS
conda install gcc_linux-64 gxx_linux-64
ln -s -f `which x86_64-conda_cos6-linux-gnu-gcc` ${CONDA_PREFIX}/bin/gcc
ln -s -f `which x86_64-conda_cos6-linux-gnu-ar` ${CONDA_PREFIX}/bin/ar
ln -s -f `which x86_64-conda_cos6-linux-gnu-g++` ${CONDA_PREFIX}/bin/g++


# INSTALLING fastQValidator
git clone https://github.com/statgen/fastQValidator.git
git clone https://github.com/statgen/libStatGen.git
cd sw/libStatGen/
make
cd ../fastQValidator/
make all
```

## Important directories

They all must be present!

- `bin/` -- all the scripts used by the pipelines are located here.
- `config/` -- configuration files for the grid engine and the pipeline are stored here.
- `data/` -- files like the reference sequence, the kraken database and the metadata (`data/metadata/[name_run]`) should be here.
- `fq/` -- `fastq` files go here.
- `logs/` -- if something goes wrong this is the directory that will help you understanding why. Also stores the `README` files of all runs.
- `smk/` -- snakemake files are here.
- `temp/` -- intermediate files are stored here.
- `results/` -- your results will be generated here.

## Preparing the run

### with your own Fastq files

1. decide the name of your run. In this case it will be `tutorial`
2. clean the `fq/`, `results/` and `temp/` directories. **Be careful when you use `rm`!**
3. create a directory inside the `log/` directory with the name of the run. In our case:

```
 mkdir -p logs/tutorial/
```

3. create a `README` file inside such directory where you will report all the commands you execute, plus your notes (Everything went smoothly? Did you have any problems? How did you overcome them?) to ensure the run is reproducible. 

```
cd logs/tutorial
touch README_tutorial.md
# edit this file
```

4. Copy the `fastq` files inside the directory `fq/`. 

If you want you can also save the `fastq` files on `fastq_db` (`/n/data1/hms/dbmi/farhat/fastq_db/`) or another directory and create soft links inside the `fq/` directory. Here is an example:

```
for i in `ls /n/data1/hms/dbmi/farhat/Roger/mmpR_BDQ_mutant_project/eis_promoter_mutant_strains_from_Milan/`; 
do
ln -s /n/data1/hms/dbmi/farhat/Roger/mmpR_BDQ_mutant_project/eis_promoter_mutant_strains_from_Milan/${i} fq/${i};
done
```

**Important note**: the files should be in the format `<sequencing_run_id>_[1|2].fastq.gz`.

5. create the `isolates_to_analyze.txt` and the `summary_runs.tsv` files inside your `logs/[run_directory]/`. The former is a list of the ids of the isolates, the latter is used to merge multiple pairs of `fastq` files (when the same strain has been sequenced multiple times).

For instance, let's say our isolates are called `IT1070` and `IT184`. 

The `logs/tutorial/isolates_to_analyze.txt` file will look like this:

```
IT1070
IT184 
```

and the `logs/tutorial/summary_runs.tsv` file will look like this (**it is tab separated!**): <mark>case with multiple run_ids</mark>

```
biosample	run_id                                                                       IT1070	IT1070                                                                      IT184	IT184 
```

6. edit the configuration file `config/config_pipeline.json`. Usually you only need to change the values of `logs_analysis` and `summary_runs`.

In our case this file will look like this:

```
{
    "logs_analysis": "logs/tutorial/",
    "summary_runs": "logs/tutorial/summary_runs.tsv",
    "path_ref_genome": "/home/lf61/lf61/repos/megapipe/RefGen/TBRefGen.fasta",
    "temp_dir":"/home/lf61/lf61/repos/megapipe/megapipe_snakemake/temp/",
    "fastq_dir":"/home/lf61/lf61/repos/megapipe/megapipe_snakemake/fq/",
    "path_to_prinseq":"/home/lf61/mfarhat/bin/prinseq-lite-0.20.4/prinseq-lite.pl",
    "path_to_kraken_db":"/n/data1/hms/dbmi/farhat/bin/kraken/tbdb/",
    "bin_fastqvalidator":"/home/lf61/mfarhat/bin/fastQValidator/bin/fastQValidator",
    "path_to_picard":"/home/lf61/mfarhat/bin/picard/picard/build/libs/picard.jar"
}
```

Now you are ready to execute the run!

### Fastq files / IDs from NCBI

1. decide the name of your run. In this case it will be `tutorial_ncbi`
2. clean the `fq/`, `results/` and `temp/` directories. **Be careful when you use `rm`!**
3. create a directory inside the `log/` directory with the name of the run. In our case:

```
 mkdir -p logs/tutorial_ncbi/
```

3. create a `README` file inside such directory where you will report all the commands you execute, plus your notes (Everything went smoothly? Did you have any problems? How did you overcome them?) to ensure the run is reproducible. 

```
cd logs/tutorial_ncbi
touch README_tutorial_ncbi.md
# edit this file
```

5. Let's say we want to analyze the strains having the following `ERR2510808` and `ERR2199868`. These are NCBI IDs which are not BioSamples (which are in the form `SAM*`, e.g. `SAMN09090639`) -- they are sequencing runs.

We therefore need to convert the IDs to BioSamples.

We create a `.txt` file with the list of the IDs (we named it `logs/tutorial_ncbi/ids.txt`). In our case will look like this:

```
ERR2510808
ERR2199868                                                               
```

 Now we use `metatools_ncbi` to convert the IDs to BioSamples:

```
 metatools_convert logs/tutorial_ncbi/ids.txt logs/tutorial_ncbi/ids_biosample_assoc.txt
```

The output file (`logs/tutorial_ncbi/ids_biosample_assoc.txt`) looks like this. The first column lists the IDs we queried and the second one lists the corresponding BioSamples:

```
ERR2510808     SAMEA1102833
ERR2199868     SAMEA104394480 
```

Given that multiple sequencing runs could be associated to the same BioSample, we may want to check if more sequencing runs are available (**Note**: if this is not relevant to you, you can directly generate the file `summary_runs.tsv` -- see below). First of all we create store our list of BioSamples inside the `log/tutorial_ncbi` directory

```
cat logs/tutorial_ncbi/ids_biosample_assoc.txt | cut -f 2 > logs/tutorial_ncbi/isolates_to_analyze.txt
```

We now download the metadata of such BioSamples:

```
mkdir -p data/metadata/tutorial_ncbi/
metatools_download runs logs/tutorial_ncbi/isolates_to_analyze.txt data/metadata/tutorial_ncbi/
```

Now we can edit the variables in the snakemake configuration file (`config/config_pipeline.json`). Usually you only need to change the values of `logs_analysis` and `summary_runs`.

We generate the file `summary_runs.tsv`:

```
./bin/generate_summary_runs.py data/metadata/tutorial_ncbi/ ./logs/tutorial_ncbi/summary_runs.tsv
```

The file looks like this:

```
biosample	run_id
SAMEA1102833    ERR2510808
SAMEA104394480  ERR2199868  
```

Now I run the script `get_list_run_ids`. It will read the list of the sequencing runs to download from NCBI:

```
./bin/get_list_run_ids.py
```

and generate the file ` logs/tutorial_ncbi/runs_to_download.txt`.

 Finally, we can download the sequencing runs from NCBI:

```
snakemake \
-s smk/download_fastq_from_ncbi.py \
-j 100 \
--cluster-config config/config_slurm_fastq_downloader.json \
--cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -p {cluster.partition}" --keep-going                                                         
```


## Executing a run

```
sbatch \
-t 3-5:00 \
-p medium \
--mem=1G \
-o wrap.log \
-e wrap.log \
--wrap='snakemake -s smk/megapipe3.py -j 5000 --cluster-config config/config_slurm_pipeline.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -p {cluster.partition}" --keep-going --latency-wait 180'
```

You can check if the pipeline is running using `squeue -u <o2_user_id>`:

```
squeue -u lf61
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
27162555     short snakejob     lf61  R       0:32      1 compute-a-16-170
27162556     short snakejob     lf61  R       0:32      1 compute-e-16-188
27162528    medium     wrap     lf61  R       0:53      1 compute-a-16-161  
```

If you want to check what is going on in real time you can have a look at the snakemake log:

```
watch tail -n 35 wrap.log
```

## Checking how things went

You can use `./bin/get_report.py` to get an overview of what went well and what did not go well:

```                                                                     ./bin/get_report.py
* SUCCEDED: 1 (50.00%);FAILED: 1 (50.00%)
* Steps where megapipe stopped:                                           - combine_runs.txt
    ['63']
```

## Moving the results to rollingDB

```
cd results/
rsync -avzP . /n/data1/hms/dbmi/farhat/rollingDB/genomic_data/
```