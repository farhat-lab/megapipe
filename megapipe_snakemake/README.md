## Generate the genomic data

I first need to download metatools_ncbi (done). I can download the metadata for the biosamples I care about. I should be able to go from/to a list of NCBI IDs to biosamples.

```
missing functionality
```

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