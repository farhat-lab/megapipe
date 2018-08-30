#!/usr/bin/env python

import sys
import os
from os import system, listdir, path, environ
import re
import subprocess as sbp
import json
import pkg_resources
import argparse
from itertools import islice
import shutil

def get_path(*args, **kw):
    """
    Get's the given path and creates it if it doesn't exist, will also construct
    a filename ending and return for quick use.
    """
    my_path = os.path.join(*args)
    if not os.path.isdir(my_path):
        os.makedirs(my_path)
    return os.path.join(my_path, kw['fn']) if 'fn' in kw else my_path

def detect_cfg_file():
    if path.isfile(environ['HOME']+"/.megapipe.json"):
        with open(environ['HOME']+"/.megapipe.json") as json_file:
            data = json.load(json_file)
            #If you update the pipeline, do not forget to update this check
            if(("picard" in data) and ("pilon" in data) and ("prinseq" in data) and ("kraken_db" in data) and ("qualimap" in data) and ("spades" in data) and ("fasta_ref" in data) and ("lineage_snp_db" in data)):
                return(data)
            else:
                print("::I could not find some required variables in the configuration file")
    else:
        print("::I do not have a configuration file. Please check the documentation to understand how to create one!")
        sys.exit()

def detect_weird_read_names(fastq):
    with open(fastq,"r") as inf:
        fline=inf.readline()
        decision=False
        if re.search("#0/[1-9]\n$",fline):
            decision=True
    return decision

def valFQ(file1,file2):
    print("::validate fastq files")
    try:
        out1 = sbp.check_output(["fastQValidator", "--file", file1])
        out2 = sbp.check_output(["fastQValidator", "--file", file2])
        m1 = re.search("FASTQ_SUCCESS", out1)
        m2 = re.search("FASTQ_SUCCESS", out2)
        if(m1 == None or m2 == None):
            return False
        return True
    except:
        return False

def write_msg(file_log,msg):
    if path.exists(file_log):
        mode = 'a'
    else:
        mode = 'w'
    with open(file_log,mode) as outf:
        if type(msg)=="byte":
            outf.write(msg.decode("ascii")+"\n")
        else:
            outf.write(msg+"\n")


def checkDRs(depths):
    """The purpose of this method is to check that at least 95% of base positions in the drug resistance conferring regions of the genome are covered 10x. This method was heavily inspired by and based upon Dr. Farhat's script "qc_report.py" in the work-horse directory"""
    okaybases = 0
    for d in depths:
        if(d >= 10):
            okaybases += 1
    okayp = okaybases / len(depths)
    return okayp

def medianLengthReadsAndReadCount(fastq):
    db_lengths=[]
    with open(fastq,"r") as fq:
        while(True):
            try:
                current_lines = islice(fq, 4)
                current_data=list(current_lines)
                if(len(current_data)==0):
                    break
            except:
                break
            current_length=len(current_data[1])
            db_lengths.append(current_length)
    # I calculate the median length
    num_reads=len(db_lengths)
    #print("* Number of reads: {}".format(num_reads))
    # I calculate the median length
    mod=num_reads % 2
    if mod==1:
        median_length=db_lengths[int(((num_reads-1)/2)-1)]
    else:
        median_length=(db_lengths[int(((num_reads)/2)-1)]+db_lengths[int(((num_reads/2)))])/2
    bp_covered=median_length*num_reads
    return(median_length,num_reads,bp_covered)

def remove_temporary_files(directory_to_remove):
    shutil.rmtree(directory_to_remove)

########
##MAIN##
########

#idea for the parameters: ::usage: {} <tag> <table_identification_strains> <dir_fastq> <output_dir> <scratch_dir> <o2_logs> <n_runs_to_consider[all|num]> <optional_flags[--skip_assembly|--keep_tmp_files]>

parser = argparse.ArgumentParser()
parser.add_argument("tag", type=str,
                    help="Tag for the isolates (usually -> biosample)")
parser.add_argument("table_identification_strains", type=str,
                    help="Path to the table identification strains")
parser.add_argument("dir_fastq", type=str,
                    help="Path to the directory containing the fastq files")
parser.add_argument("output_dir", type=str,
                    help="Path to the directory where the results will be written")
parser.add_argument("scratch_dir", type=str,
                    help="Path to the directory where the temporary files will be written")
parser.add_argument("logs_dir", type=str,
                    help="Path to the directory where the slurm scripts and logs will be written")
parser.add_argument("n_runs_to_consider", type=int,
                    help="Maximum number of sequencing runs to consider for the analyses")
parser.add_argument("-k", "--keep_tmp", action="store_true",
                    help="keeps the temporary files")
parser.add_argument("-s", "--skip_assembly", action="store_true",
                    help="skips the assembly ")
args = parser.parse_args()

tag=args.tag
table=args.table_identification_strains
dir_fastq=args.dir_fastq
out_dir=args.output_dir
scratch_dir=args.scratch_dir
#this directory stores the logs from o2 (output and error)
log_dir=args.logs_dir
runs_to_consider=args.n_runs_to_consider
keep_tmp=args.keep_tmp
skip_assembly=args.skip_assembly

#I read the configuration file
data_json=detect_cfg_file()

if path.exists(path.join(out_dir, tag)):
    print(":: I found some results for {}. Skipping!".format(tag))
    sys.exit()

# I create new directories (if they do not exist)
get_path(out_dir, tag)
get_path(scratch_dir, tag)
get_path(log_dir)

# the log file will the the following
file_log = path.join(out_dir, tag, "{}.out".format(tag))

# I log the version of megapipe I am using
write_msg(file_log,":: I am using megapipe v{}".format(pkg_resources.get_distribution("megapipe").version))

#I check if the reference is ok. If so, I index it (if needed).
write_msg(file_log,":: Indexing the reference")
if path.isfile(data_json["fasta_ref"]):
    if path.isfile(data_json["fasta_ref"]+".bwt"):
        write_msg(file_log,"  * I detected the index of the reference. I do not generate a new index")
    else:
        write_msg(file_log,"  * I did not detect the index of the reference. I will generate it")
        cmd="bwa index "+data_json["fasta_ref"]
        system(cmd)

# I get the data about my strain
write_msg(file_log,":: I am retrieving the data about the runs to analyze from {}".format(table))
runs_to_analyze=[]
with open(table,"r") as inp:
    fields=inp.readline().rstrip("\n").split("\t")
    for line in inp:
        if line =="\n":
            continue
        entry=line.rstrip("\n").split("\t")
        try:
            idx_public_xref=fields.index("public_xref")
            public_xref=entry[idx_public_xref]
            has_public_xref=True
        except:
            has_public_xref=False
        try:
            idx_internal_xref=fields.index("internal_xref")
            internal_xref=entry[idx_internal_xref]
            has_internal_xref=True
        except:
            has_internal_xref=False
        if not (has_public_xref or has_internal_xref):
            print("  [ERROR] I did not find neither the public_xref nor the internal_xref")
            sys.exit()
        if "public_xref" in vars():
            if public_xref==tag:
                # I get the runs
                idx_runs=fields.index("runs")
                runs=entry[idx_runs].split(",")
                # I reformat the information
                for run in runs:
                    current_run_formatted = run + ":" + path.join(dir_fastq, run+"_1.fastq.gz,") + path.join(dir_fastq, run+"_2.fastq.gz")
                    runs_to_analyze.append(current_run_formatted)
        if "internal_xref" in vars():
            if internal_xref==tag:
                # I get the runs
                idx_runs=fields.index("internal_fastq_files")
                internal_runs=[entry[idx_runs]]
                runs_to_analyze.extend(internal_runs)

# Now foreach run I unzip it, I check that the fastq files are ok.

#I unzip the fastq files in a directory on scratch2
write_msg(file_log,":: Checking the runs")

# I will use this dictionary to  the runs (0=everything is fine, 1=there is a problem, so the run will be excluded)
data_runs={}

for current_run in runs_to_analyze:

    # I get the info about this run
    general_info=current_run.split(":")
    run=general_info[0]
    write_msg(file_log,"  * analysing run {}".format(run))
    fastq_files=general_info[1].split(",")

    # the flag is 0 for now.
    data_runs[run]={}
    data_runs[run]["flag"]={}
    # Explanation for the flag: zero is good, one means that the run is not good
    data_runs[run]["flag"]=0
    data_runs[run]["fq1"]={}
    data_runs[run]["fq1"]=fastq_files[0]
    data_runs[run]["fq2"]={}
    data_runs[run]["fq2"]=fastq_files[1]
    data_runs[run]["bp_coverage"]={}
    data_runs[run]["fq2"]=fastq_files[1]

    #I rename the variables fqf1 and fqf2 since the fastq I will use are the ones on the scratch, so that I do not touch the original ones
    fqf1 = path.join(scratch_dir, tag, general_info[0] + "_1.fastq")
    fqf2 = path.join(scratch_dir, tag, general_info[0] + "_2.fastq")

    gz_pattern=re.compile("\.gz$")
    #if the file has .gz extension
    if re.search(gz_pattern,fastq_files[0]):
        # I unzip the fastq files of this run
        #I unzip the fastq files in a directory on scratch2
        write_msg(file_log,"    - Unzipping fastq files (on {0})".format(scratch_dir))
        cmd=["zcat",fastq_files[0]]
        try:
            with open("{0}/{1}/{2}_1.fastq".format(scratch_dir,tag, general_info[0]),"w") as fqf:
                sbp.call(cmd,stdout=fqf)
        except:
            write_msg(file_log,"      + [WARNING] I found a problem while unzipping {}".format(general_info[0]))
            data_runs[run]["flag"]=1
            continue
        cmd=["zcat",fastq_files[1]]

        try:
            with open("{0}/{1}/{2}_2.fastq".format(scratch_dir,tag, general_info[0]),"w") as fqf:
                sbp.call(cmd,stdout=fqf)
        except:
            write_msg(file_log,"      + [WARNING] I found a problem while unzipping {}".format(general_info[0]))
            data_runs[run]["flag"]=1
            continue
    else:
        # If the files are already fastq files...
        write_msg(file_log,"    - The fastq files are already unzipped. I copy the fastq files")
        cmd=["cp",fastq_files[0],fqf1]
        sbp.call(cmd)
        cmd=["cp",fastq_files[1],fqf2]
        sbp.call(cmd)

    #I check the fastq files
    write_msg(file_log,"    - Validating fastq files")
    try:
        valFQ(fqf1,fqf2)
        write_msg(file_log,"      + OK!")
    except:
        write_msg(file_log,"      + [WARNING] Fastq files are NOT valid ")
        data_runs[run]["flag"]=1
        continue

    # I check the names of the reads
    write_msg(file_log,"    - Checking the names of the reads")
    test_names1=detect_weird_read_names(fqf1)
    test_names2=detect_weird_read_names(fqf2)
    if((test_names1==True) or (test_names2==True)):
        write_msg(file_log,"      + I found some weird names. I am fixing them!")
        cmd="megapipe-correct-names-reads.py {}".format(fqf1)
        system(cmd)
        cmd="megapipe-correct-names-reads.py {}".format(fqf2)
        system(cmd)
    else:
        write_msg(file_log,"      + OK!")

    # I trim the reads with Prinseq -- this should happen in the scratch in order to save space
    write_msg(file_log,"    - Trimming reads with Prinseq")
    prinseq = get_path(out_dir, tag, "prinseq")
    path_to_prinseq=data_json["prinseq"]
    cmd=["perl", path_to_prinseq,
        "-fastq",fqf1,
        "-fastq2",fqf2,
        "-out_format", "3",
        "-out_good", path.join(scratch_dir, tag, run + "-trimmed"),
        "-out_bad", "null",
        "-log", path.join(prinseq, "{}-prinseq.log".format(run)),
        "-min_qual_mean", "20",
        "-verbose"]
    print(" ".join(cmd))
    sbp.call(cmd)
    write_msg(file_log,"      + Please see the Prinseq report in {0}/{1}/prinseq/{2}-prinseq.log".format(out_dir,tag, run))
    # I check if the trimming went well.
    trflstem1 = run + "-trimmed_1.fastq"
    trflstem2 = run + "-trimmed_2.fastq"
    trfl1 = path.join(scratch_dir, tag, trflstem1)
    trfl2 = path.join(scratch_dir, tag, trflstem2)
    files = listdir(path.join(scratch_dir, tag))
    if(not(trflstem1 in files and trflstem2 in files)):
        write_msg(file_log, "      + [WARNING] Prinseq failed. Please have a look at the log file")
        data_runs[run]["flag"]=1
        continue

    # I classify the reads with Kraken
    write_msg(file_log,"    - Classifying reads with Kraken")
    kraken = get_path(out_dir, tag, "kraken")
    path_to_krakendb=data_json["kraken_db"]
    cmd=["kraken",
        "--fastq-input", trfl1,
        "--output", path.join(kraken, trflstem1 + ".krkn"),
        "--db", path_to_krakendb]
    print(" ".join(cmd))
    e=sbp.check_output(cmd,stderr=sbp.STDOUT)
    mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
    tbperc1 = float(mat.groups()[1]) / 100
    cmd=["kraken",
        "--fastq-input", trfl2,
        "--output", path.join(kraken, trflstem2 + ".krkn"),
        "--db", path_to_krakendb]
    print(" ".join(cmd))
    e=sbp.check_output(cmd,stderr=sbp.STDOUT)
    mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
    tbperc2 = float(mat.groups()[1]) / 100
    cmd="gzip {0}/{1}/kraken/{2}.krkn".format(out_dir, tag, trflstem1)
    system(cmd)
    cmd="gzip {0}/{1}/kraken/{2}.krkn".format(out_dir, tag, trflstem2)
    system(cmd)
    write_msg(file_log,"      + Please see the Kraken report in {0}/{1}/kraken/{2}.krkn.gz and {0}/{1}/kraken/{3}.krkn.gz".format(out_dir, tag, trflstem1, trflstem2))
    # We still need to check that more than 90% of the sequences are from mycobacterium tuberculosis by making sense of output from Kraken
    # write doen in the log tbperc1 and 2
    if(tbperc1 < 0.9):
        write_msg(file_log,"      + [WARNING] Less than 90% of reads in the first fastq file belonged to Mycobacterium tuberculosis")
        data_runs[run]["flag"]=1
        continue
    if(tbperc2 < 0.9):
        write_msg(file_log,"      + [WARNING] Less than 90% of reads in the second fastq file belonged to Mycobacterium tuberculosis")
        data_runs[run]["flag"]=1
        continue
    # I determine how good is this run
    write_msg(file_log,"    - I count the reads and calculate the median read length")
    qual_metrics=medianLengthReadsAndReadCount(trfl1)
    write_msg(file_log,"      + Num of reads: {}; Median read length: {} bp".format(qual_metrics[1],qual_metrics[0]))
    data_runs[run]["bp_coverage"]=qual_metrics[2]


# Now I can combine the runs that succeeded
fq_comb1 = path.join(scratch_dir, tag, tag + "-combined_1.fastq")
fq_comb2 = path.join(scratch_dir, tag, tag + "-combined_2.fastq")

#I check that at least one run is OK. I define a variable to count the good runs
runs_ok=0

#I create a dictionary: Run => bp_coverage
ranking_runs={}

for run in data_runs:
    if(data_runs[run]["flag"]==0):
        runs_ok=runs_ok+1
        ranking_runs[run]=data_runs[run]["bp_coverage"]

write_msg(file_log,":: Summaryzing the situation of the sequencing runs")
# If there are no good runs:
if runs_ok == 0:
    write_msg(file_log,"  * [ERROR] I found some problems in each one of the runs you provided. The analysis stops here!")
    sys.exit()
else:
    write_msg(file_log,"  * {} runs are OK ({})".format(runs_ok,",".join(ranking_runs.keys())))

selected_runs=[]
# I find the best runs
if(runs_to_consider>=len(ranking_runs)):
    runs_to_consider=len(ranking_runs)
else:
    try:
        runs_to_consider=int(runs_to_consider)
    except:
        write_msg(file_log,"  * [ERROR] I have some problems to understand how many runs I should consider for the mapping and the assembly!")
        sys.exit()
if(len(ranking_runs)>runs_to_consider):
    bp_coverages=sorted(ranking_runs.values(),reverse=True)
    first_two=set(bp_coverages[0:runs_to_consider])
    for run in ranking_runs:
        if ranking_runs[run] in first_two:
            selected_runs.append(run)
else:
    selected_runs=ranking_runs.keys()

write_msg(file_log,"  * Selected runs: {}".format(",".join(selected_runs)))

# I reinitialize the fastq files.
cmd="> {}".format(fq_comb1)
system(cmd)
cmd="> {}".format(fq_comb2)
system(cmd)

for run in selected_runs:
    trflstem1 = path.join(scratch_dir, tag, run + "-trimmed_1.fastq")
    trflstem2 = path.join(scratch_dir, tag, run + "-trimmed_2.fastq")
    cmd="cat {0} >> {1}".format(trflstem1,fq_comb1)
    system(cmd)
    cmd="cat {0} >> {1}".format(trflstem2,fq_comb2)
    system(cmd)


# Aligning reads to the reference
write_msg(file_log,":: Aligning reads with bwa")
samfile = path.join(scratch_dir, tag, "{}.sam".format(tag))
#piper = popen("bwa mem -M -R '@RG\tID:<unknown>\tSM:<unknown>\tPL:<unknown>\tLB:<unknown>\tPU:<unknown>' RefGen/TBRefGen.fasta {0} {1} > {2}".format(sctrfl1, sctrfl2, samfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
try:
    cmd=["bwa","mem","-M", data_json["fasta_ref"], fq_comb1, fq_comb2]
    print(" ".join(cmd))
    with open(samfile,"w") as samf:
        sbp.call(cmd,stdout=samf)
    write_msg(file_log,"  * OK!")
except:
    write_msg(file_log,"  * [ERROR] I had some problems with bwa! Check the logs from the grid engine for more information")
    sys.exit()

# Sorting and converting to bam
write_msg(file_log,":: Sorting sam file with Picard")
bamfile = path.join(scratch_dir, tag, "{}.sorted.bam".format(tag))
path_to_picard=data_json["picard"]
try:
    cmd=["java","-Xmx16G",
    "-jar", path_to_picard, "SortSam",
    "INPUT={}".format(samfile),
    "OUTPUT={}".format(bamfile),
    "SORT_ORDER=coordinate"
    ]
    print(" ".join(cmd))
    sbp.call(cmd)
    write_msg(file_log,"  * OK!")
except:
    write_msg(file_log,"  * [ERROR] I had some problems with picard!")
    sys.exit()

# Quality control
write_msg(file_log,":: Evaluating mapping with Qualimap")
path_to_qualimap=data_json["qualimap"]
system(path_to_qualimap + " bamqc -bam {0} --outfile {1}.pdf --outformat PDF".format(bamfile, tag))
qualimap = get_path(out_dir, tag, "qualimap")
system("mv {0}/* {1}".format(bamfile.replace(".bam", "_stats"), qualimap))
write_msg(file_log,"  * Please see the Qualimap report in {}".format(qualimap))

# Removing duplicates
write_msg(file_log,":: Removing duplicates from bam file with Picard")
out_bam = get_path(out_dir, tag, "bam")
drbamfile = bamfile.replace(".bam", ".duprem.bam")
try:
    cmd=["java","-Xmx32G",
        "-jar",path_to_picard,
        "MarkDuplicates",
        "I={}".format(bamfile),
        "O={}".format(drbamfile),
        "REMOVE_DUPLICATES=true",
        "M={}".format(drbamfile[:-4]+".metrics"),
        "ASSUME_SORT_ORDER=coordinate"
]
    sbp.call(cmd) # help from http://seqanswers.com/forums/showthread.php?t=13192
    sbp.call(["cp", drbamfile, path.join(out_bam,drbamfile[drbamfile.rindex("/") + 1:])])
    write_msg(file_log,"  * OK!")
    write_msg(file_log,"  * Please see the Picard MarkDuplicates report in {0}/{1}".format(scratch_dir,drbamfile[drbamfile.rindex("/") + 1:-4]))
except:
    write_msg(file_log,"  * [ERROR] I had some problems with picard!")
    sys.exit()

write_msg(file_log,":: Calculating depth")
out = sbp.check_output(["samtools", "depth", "-a", drbamfile])
dmat = re.findall("(.+\t)([0-9]+)(\n)", out.decode("ascii"))
depths = [float(t[1]) for t in dmat]

#gencov = sum(depths) / len(depths)
#outputer.write("Genome coverage: {}\n".format(gencov))
gencovprop = checkDRs(depths)
write_msg(file_log,"  * The percent of H37Rv bases that have a coverage of at least 10x is {}.".format(str(gencovprop * 100)))
if(gencovprop < 0.95): # The threshold is 95%
    write_msg(file_log,"    - [ERROR] The percent of H37Rv bases that have a coverage of at least 10x is less than 95%.")
    sys.exit()
fileout_depth = get_path(out_dir, tag, "depth", fn=tag + ".depth")
write_msg(fileout_depth,out.decode("ascii"))
# I compress the depth output file
cmd="gzip {}".format(fileout_depth)
system(cmd)
refcov = 0
for d in depths:
    if(d > 0):
        refcov += 1
write_msg(file_log,"  * Percent of reference genome covered: {}".format(refcov / len(depths)))

write_msg(file_log,":: Indexing {}".format(drbamfile))
sbp.call(["samtools", "index", drbamfile])

# I call the variants with pilon
write_msg(file_log,":: Variant calling with Pilon")
out_pi = get_path(out_dir, tag, "pilon")
get_path(scratch_dir, tag, "pilon")
path_to_pilon=data_json["pilon"]
try:
    # update all subprocess calls so that they look like this.
    out_pilon = os.path.join(scratch_dir, tag, "pilon", tag)
    cmd=["java", '-Xmx32G', '-jar', path_to_pilon,
       '--genome', data_json["fasta_ref"],
       '--bam', drbamfile,
       '--output', out_pilon,
       '--variant']
    print(' '.join(cmd))
    sbp.call(cmd)
    # I reduce the size of the pilon output
    cmd=["megapipe-vcf-cutter.py", out_pilon+".vcf", path.join(out_pi, tag+".vcf")]
    sbp.call(cmd)
    # I copy the data
    ## I copy and compress the original pilon vcf
    system("cp {0}.vcf {1}/{2}_full.vcf".format(out_pilon, out_pi, tag))
    system("gzip {0}/{1}_full.vcf".format(out_pi,tag))
    ## I copy the fasta file
    system("cp {0}.fasta {1}".format(out_pilon, out_pi))
    write_msg(file_log,"  * OK")

except:
    write_msg(file_log,"  * [ERROR] I had some problems with pilon!")
    sys.exit()

# I calculate the lineage using the fast-lineage-caller
path_to_lineage_snp_db=data_json["lineage_snp_db"]
scratch_flc = get_path(scratch_dir, tag, "fast-lineage-caller")
out_flc = get_path(out_dir, tag, "fast-lineage-caller")

try:
    cmd=["vrtTools-vcf2vrt.py", path.join(out_dir, tag, "pilon", tag+".vcf"), path.join(scratch_flc, tag + ".vrt"), "1"]
    sbp.call(cmd)
    cmd=["FastLineageCaller-assign2lineage.py", path_to_lineage_snp_db, path.join(scratch_flc, tag + ".vrt")]
    with open(path.join(out_flc, tag+".lineage"), "w") as lin:
        sbp.call(cmd,stdout=lin)
except:
    write_msg(file_log,"  * [ERROR] I had some problems with the lineage caller!")
    sys.exit()


if skip_assembly:
    write_msg(file_log,":: I skip the genome assembly!")
    if not keep_tmp:
        write_msg(file_log,":: I delete the temporary files")
        remove_temporary_files(os.path.join(scratch_dir, tag))
    sys.exit()

# I generate the assembly with spades
write_msg(file_log,":: Generating the assembly with Spades")
out_sp = get_path(out_dir, tag, "spades")
scratch_sp = get_path(scratch_dir, tag, "spades")
#-t (treads); -m (memory, in Gb)
path_to_spades=data_json["spades"]
try:
    cmd=["python2", path_to_spades,
        '-t','1',
        '-m','30',
        '--careful',
        '--pe1-1',fq_comb1,
        '--pe1-2',fq_comb2,
        '-o', path.join(scratch_dir, tag, "spades")]
    print(' '.join(cmd))
    sbp.call(cmd)
    # I copy the assembly data and the logs
    sbp.call(["cp", path.join(scratch_sp, "scaffolds.fasta"), path.join(out_sp, tag+"_"+"scaffolds.fasta")])
    sbp.call(["cp", path.join(scratch_sp, "contigs.fasta"), path.join(out_sp, tag+"_"+"contigs.fasta")])
    sbp.call(["cp", path.join(scratch_sp, "spades.log"), path.join(out_sp, tag+"_"+"spades.log")])
    write_msg(file_log,"  * OK")
except:
    write_msg(file_log,"  * [ERROR] I had some problems with spades!")
    sys.exit()

if not keep_tmp:
    write_msg(file_log,":: I delete the temporary files")
    remove_temporary_files(os.path.join(scratch_dir, tag))
