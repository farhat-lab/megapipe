## Important directories

They all must be present!

- `bin/` -- all the scripts of the pipeline and the snakemake file are located here
- `config/` -- configuration files for the grid engine and the pipeline are stored here. 
- `fq/` -- `fastq` go here
- `logs/` -- if something goes wrong this is the directory that will help you understanding why
- `temp/` -- intermediate files are stored here
- `results/` -- your results will be generated here

## Preparing the run

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

4. The `fastq` files of the isolates you want to analyze are from NCBI or we sequenced them/our collaborators sent us the files?
   - If they are on NCBI, go to step 5
   - if we sequenced them/our collaborators sent us the files, go to step <mark>6</mark>

5. <mark>[expand this]</mark>                                                                                  

6. Copy the `fastq` files inside the directory `fq/`. 

If you want you can also save the `fastq` files on `fastq_db` (`/n/data1/hms/dbmi/farhat/fastq_db/`) or another directory and create soft links inside the `fq/` directory. Here is an example:

```
for i in `ls /n/data1/hms/dbmi/farhat/Roger/mmpR_BDQ_mutant_project/eis_promoter_mutant_strains_from_Milan/`; 
do
ln -s /n/data1/hms/dbmi/farhat/Roger/mmpR_BDQ_mutant_project/eis_promoter_mutant_strains_from_Milan/${i} fq/${i};
done
```

**Important note**: the files should be in the format `<sequencing_run_id>_[1|2].fastq.gz`.

7. create the `isolate_to_analyze.txt` and the `summary_runs.tsv` files inside your `logs/[run_directory]/`. The former is a list of the ids of the isolates, the latter is used to merge multiple pairs of `fastq` files (when the same strain has been sequenced multiple times).

For instance, let's say our isolates are called `IT1070` and `IT184`. 

The `logs/tutorial/isolates_to_analyze.txt` file will look like this:

```
IT1070
IT184 
```

and the `logs/tutorial/summary_runs.tsv` file will look like this (**it is tab separated!**):

```
biosample	run_id                                                                       IT1070	IT1070                                                                             IT184	IT184 
```

8. edit the configuration file `config/config.json`. Usually you only need to change the values of `logs_analysis` and `summary_runs`.

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

Now you are ready to go!

## Executing a run

```
sbatch -t 3-5:00 -p medium --mem=5G -o wrap.log -e wrap.log --wrap='snakemake -s bin/megapipe3.py -j 100 --cluster-config config/cluster_test.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -p {cluster.partition}" --keep-going --latency-wait 180'
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

You can use `/bin/get_report.py` to get an overview of what went well and what did not go well:
```
./bin/get_report.py
* SUCCEDED: 0 (0.00%);FAILED: 14 (100.00%)                                                * Steps where megapipe stopped:
    - combine_runs.txt
['63']                                                                         
```













## Todo



If I have a list of biosamples I can download the metadata about the sequencing runs runs. 

```
metatools_download runs isolates_to_analyze.txt ../../metadata/sources/metatools_ncbi/
```

SAMN02231188 this is an isolate that should work
I build the table with the with the associations <biosample> <run id>:

```
./bin/generate_summary_runs.py metadata/sources/metatools_ncbi/ metadata/summary_tables/test2.tsv
```

I get the list of run ids to download from the table (Please edit the variables _logs\_analysis_ and _summary_runs_ on _config/config.json_).

```
./bin/get_list_run_ids.py
```

I download the runs:

```
snakemake -s bin/download_fastq_from_ncbi.py -j 100 --cluster-config config/cluster_fastq_downloader.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -p {cluster.partition}" 
```

I set up the configuration file for megapipe. 

I run megapipe:

```
snakemake -s bin/megapipe3.py -j 100 --cluster-config config/cluster_test.json --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -p {cluster.partition}"
```

​	

## 