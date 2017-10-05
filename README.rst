### Version

### Todo
 * build the package
 * reference genome -- do something intelligent -- check if it is present and if not generate it -- the megapipe-launch could check this
 * rename the scripts so that we know that we are talking about megapipe
 * module to run the calculations on slurm // if I can do it so tht it can run the things on bsub -- even better
 * user should install its own dependencies -- remove absolute paths


### Workflow



(0) if megapipe is not in your path, you should go into the directory where the megapipe binaries are located (with the "cd" command) and run the script ./add-mp-to-path.sh. In my case:
```
cd lf61/work-horse/megapipe/bin/
./add-mp-to-path.sh

```
Then, in order to be able to have access to the megapipe commands, you should reload your .bashrc configuration file

```
. ~/.bashrc
```
Theoretically, now you should be ready to use megapipe (if all the paths and permission are OK).


(1) load the python3 module and create a directory where you want to store your data
```
module load dev/python/3.4.2
mkdir mp_out
cd mp_out
```

(2) copy the files of the reference genome -- they are needed by the pipeline -- they are located on the RefGen directory on inside the megapipe directory
```
cp -r ~/lf61/work-horse/megapipe/RefGen .
```

(3) Generate the list of the files to analyze. 

Here is the general synthax  of the command:
```
generate-acclist.py <dir_fastq> <tag_pair_end_fastq1> <extension_fastq1> > <output_file>
```
Here is an example of how to run the command:

```
generate-acclist.py ../fastq_db/reseqtb/IS-1001/ _1 _1.fastq.gz > acclist2.0
```

The output file of generate-acclist.py is a tab separated value file (acclist file or accession list file) that contains the tag of the genome and the absolute paths of the 2 fastq files. The tag is automatically generated from the fastq file names. Here is one line of a sample acclist file:
```
00R0025 ../fastq_db/pools/00-R0025.1.fastq.gz   ../fastq_db/pools/00-R0025.2.fastq.gz

```


(4) Run the pipeline
Here is the general synthax  of the command:
```
megapipelaunch.py <acclist_file> <output_dir> <scratch_dir> <first_genome> <last_genome> 0
megapipelaunch.py <acclist_file> <output_dir> <scratch_dir> <first_genome> <last_genome> 1

```
You need two commands because the first one generates the scripts that will be run on orchestra (runmode 0), while the second command actually launches the jobs (runmode 1).
Here is an example of how to run the command in real life:

```
megapipelaunch.py acclist2.0 test6 /n/scratch2/lf61/mp/ 1 1 0
megapipelaunch.py acclist2.0 test6 /n/scratch2/lf61/mp/ 1 1 1
```
acclist2.0 is the accession list file you generated in the step (3); "test6" is your output directory where the results will be stored; "/n/scratch2/lf61/mp/" is a directory in the scratch that will contain your partial results (please be sure that this directory exists and it is writable).

In the example I proposed, the first and the last genomes variables are both set at 1, so it means that megapipe will read only the first entry of your acclist file and analyze it. If you want to analyze more entries or you want to skip some genomes, you can change these parameters. 

For instance here is an example that show how to launch a megapipe analysis for all the genomes of a dataset of the RESEQTB project:
```
megapipelaunch.py 00-metadata/reseqtb-RESEQTB_Dec16-acclist.tsv 01-mp_out/ /n/scratch2/lf61/mp/ 1 2432 1
```

**GOOD LUCK for your analyses!**

**NOTE: remember to clean the scratch from time to time!** 




```
from gridmanager import gridpuppeteer as gp
a=gp.GridEngine()
a.generate_script("prova.sh","short","12:00","prova.out","10M","wget http://poisson.phc.unipi.it/~freschi/img/luca.jpg")
a.launch_job("prova.sh")
```
