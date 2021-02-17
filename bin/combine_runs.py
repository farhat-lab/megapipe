#!/usr/bin/env python
import argparse
import json
import pandas as pd
import sys
import os
import re
import subprocess as sbp
from itertools import islice

parser = argparse.ArgumentParser()
parser.add_argument("biosample", type=str,
                    help="Biosample we are interested in")
parser.add_argument("summary_runs", type=str,
                    help="tsv containing the associations biosample => run_ids")
parser.add_argument("fastq_dir", type=str,
                    help="directory containing the fastq files")
parser.add_argument("temp_dir", type=str,
                    help="directory containing the temporary files of the results")
parser.add_argument("num_runs_to_consider", type=int,
                    help="tsv containing the associations biosample => run_ids")
args = parser.parse_args()

# I define the functions needed to process the samples

def detect_weird_read_names(fastq):
    with open(fastq,"r") as inf:
        fline=inf.readline()
        decision=False
        if re.search("#0/[1-9]\n$",fline):
            decision=True
    return decision

def valFQ(file1,file2):
    print("- validate fastq files")
    try:
        out1 = sbp.check_output([config["bin_fastqvalidator"], "--file", file1])
        out2 = sbp.check_output([config["bin_fastqvalidator"], "--file", file2])
        m1 = re.search("FASTQ_SUCCESS", out1)
        m2 = re.search("FASTQ_SUCCESS", out2)
        if(m1 == None or m2 == None):
            return False
        return True
    except:
        return False

def medianLengthReadsAndReadCount(fastq):
    print(fastq)
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
    print("* Number of reads: {}".format(num_reads))
    # I calculate the median length
    mod=num_reads % 2
    if mod==1:
        median_length=db_lengths[int(((num_reads-1)/2)-1)]
    else:
        median_length=(db_lengths[int(((num_reads)/2)-1)]+db_lengths[int(((num_reads/2)))])/2
    bp_covered=median_length*num_reads
    return(median_length,num_reads,bp_covered)


# I load the configuration file
## It is needed to know where kraken and prinseq are
import json
with open('./config/config.json') as inp:
    config = json.load(inp)

try:
    os.makedirs(config["temp_dir"] + args.biosample)
except:
    if(os.path.isdir(config["temp_dir"] + args.biosample)):
        pass
    else:
        raise Exception("- [ERROR] I cannot create the directory {}".format(config["temp_dir"]))

# I get the data about my strain
print("- I am retrieving the data about the runs to analyze from {}".format(args.summary_runs))
tab=pd.read_csv(args.summary_runs, sep="\t")
tab_sel=tab.loc[tab["biosample"]==args.biosample]
## If I do not find anything, the analysis finishes here
if tab_sel.shape[0] == 0:
    raise Exception("- [ERROR] Sorry! I did not find the biosample you provided in the summary runs file.") 
## If I find multiple lines corresponding to one biosample, something went wrong too
elif tab_sel.shape[0] > 1:
    raise Exception("- [ERROR] There are multiple lines in the summary_runs file that match the biosample you provided.") 
else:
    runs=tab_sel.iloc[0]["run_id"].split(",")
    runs_formatted=[]
    for run in runs:
        current_run_formatted = run + ":" + os.path.join(args.fastq_dir, run+"_1.fastq.gz,") + os.path.join(args.fastq_dir, run+"_2.fastq.gz")
        runs_formatted.append(current_run_formatted)

    # Now foreach run I unzip it, I check that the fastq files are ok.
    print("- Checking the runs")
    # I will use this dictionary to store the info about the runs (0=everything is fine, 1=there is a problem, so the run will be excluded)
    data_runs={}
    for current_run in runs_formatted:
        # I get the info about this run
        general_info=current_run.split(":")
        run=general_info[0]
        print("  - analysing run {}".format(run))
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
        fqf1 = os.path.join(args.temp_dir, args.biosample, run + "_1.fastq")
        fqf2 = os.path.join(args.temp_dir, args.biosample, run + "_2.fastq")

        #if the file has .gz extension
        gz_pattern=re.compile("\.gz$")
        if re.search(gz_pattern,fastq_files[0]):
            # I unzip the fastq files of this run
            #I unzip the fastq files in a directory on scratch2
            print("    - Unzipping fastq files ({0} => {1})".format(args.fastq_dir, args.temp_dir + args.biosample))
            cmd=["zcat",fastq_files[0]]
            try:
                with open(fqf1,"w") as fqf:
                    sbp.call(cmd,stdout=fqf)
            except:
                print("      - [WARNING] I found a problem while unzipping {}".format(fastq_files[0]))
                data_runs[run]["flag"]=1
                continue
            cmd=["zcat",fastq_files[1]]
            try:
                with open(fqf2,"w") as fqf:
                    sbp.call(cmd,stdout=fqf)
            except:
                print("      - [WARNING] I found a problem while unzipping {}".format(fastq_files[1]))
                data_runs[run]["flag"]=1
                continue
        else:
            # If the files are already fastq files...
            print("    - The fastq files are already unzipped. I copy the fastq files")
            cmd=["ln",fastq_files[0],fqf1]
            sbp.call(cmd)
            cmd=["ln",fastq_files[1],fqf2]
            sbp.call(cmd)

        #I check the fastq files
        print("    - Validating fastq files")
        try:
            valFQ(fqf1,fqf2)
            print("      - OK!")
        except:
            print("      - [WARNING] Fastq files are NOT valid ")
            data_runs[run]["flag"]=1
            continue

        # I check the names of the reads
        print("    - Checking the names of the reads")
        test_names1=detect_weird_read_names(fqf1)
        test_names2=detect_weird_read_names(fqf2)
        if((test_names1==True) or (test_names2==True)):
            print("      - I found some weird names. I am fixing them!")
            cmd="./bin/correct_names_reads.py {}".format(fqf1)
            os.system(cmd)
            cmd="./bin/correct_names_reads.py {}".format(fqf2)
            os.system(cmd)
        else:
            print("      - OK!")

        # I trim the reads with Prinseq -- this should happen in the scratch in order to save space
        print("    - Trimming reads with Prinseq")
        dir_prinseq_results = "results/" + args.biosample + "/prinseq/"
        try:
            os.makedirs(dir_prinseq_results)
        except:
            if(os.path.isdir(dir_prinseq_results)):
                pass
            else:
                raise Exception("- [ERROR] I cannot create the directory {}".format(dir_prinseq_results))

        path_to_prinseq=config["path_to_prinseq"]
        cmd=["perl", path_to_prinseq,
            "-fastq",fqf1,
            "-fastq2",fqf2,
            "-out_format", "3",
            "-out_good", os.path.join(args.temp_dir, args.biosample, run + "-trimmed"),
            "-out_bad", "null",
            "-log", os.path.join(dir_prinseq_results, "{}-prinseq.log".format(run)),
            "-min_qual_mean", "20",
            "-verbose"]
        print(" ".join(cmd))
        sbp.call(cmd)
        print("        - Please see the Prinseq report in {0}/{1}/prinseq/{2}-prinseq.log".format("results", args.biosample, run))
        # I check if the trimming went well.
        trflstem1 = run + "-trimmed_1.fastq"
        trflstem2 = run + "-trimmed_2.fastq"
        trfl1 = os.path.join(args.temp_dir, args.biosample, trflstem1)
        trfl2 = os.path.join(args.temp_dir, args.biosample, trflstem2)
        files = os.listdir(os.path.join(args.temp_dir, args.biosample))
        if(not(trflstem1 in files and trflstem2 in files)):
            print("      - [WARNING] Prinseq failed. Please have a look at the log file.")
            data_runs[run]["flag"]=1
            continue

        # I classify the reads with Kraken
        print("    - Classifying reads with Kraken")
        dir_kraken_results = "results/" + args.biosample + "/kraken/"
        try:
            os.makedirs(dir_kraken_results)
        except:
            if(os.path.isdir(dir_kraken_results)):
                pass
            else:
                raise Exception("- [ERROR] I cannot create the directory {}".format(dir_kraken_results))
        path_to_krakendb=config["path_to_kraken_db"]
        cmd=["kraken",
            "--fastq-input", trfl1,
            "--output", os.path.join(dir_kraken_results, trflstem1 + ".krkn"),
            "--db", path_to_krakendb]
        print(" ".join(cmd))
        e=sbp.check_output(cmd,stderr=sbp.STDOUT)
        mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
        tbperc1 = float(mat.groups()[1]) / 100
        cmd=["kraken",
            "--fastq-input", trfl2,
            "--output", os.path.join(dir_kraken_results, trflstem2 + ".krkn"),
            "--db", path_to_krakendb]
        print(" ".join(cmd))
        e=sbp.check_output(cmd,stderr=sbp.STDOUT)
        mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
        tbperc2 = float(mat.groups()[1]) / 100
        cmd="gzip {0}/{1}/kraken/{2}.krkn".format("results", args.biosample, trflstem1)
        os.system(cmd)
        cmd="gzip {0}/{1}/kraken/{2}.krkn".format("results", args.biosample, trflstem2)
        os.system(cmd)
        print("      - Please see the Kraken report in {0}/{1}/kraken/{2}.krkn.gz and {0}/{1}/kraken/{3}.krkn.gz".format("results", args.biosample, trflstem1, trflstem2))
        # We still need to check that more than 90% of the sequences are from mycobacterium tuberculosis by making sense of output from Kraken
        # write doen in the log tbperc1 and 2
        if(tbperc1 < 0.9):
            print("      - [WARNING] Less than 90% of reads in the first fastq file belonged to Mycobacterium tuberculosis")
            data_runs[run]["flag"]=1
            continue
        if(tbperc2 < 0.9):
            print("      - [WARNING] Less than 90% of reads in the second fastq file belonged to Mycobacterium tuberculosis")
            data_runs[run]["flag"]=1
            continue
        # I determine how good is this run
        print("    - I count the reads and calculate the median read length")
        qual_metrics=medianLengthReadsAndReadCount(trfl1)
        print("    - Num of reads: {}; Median read length: {} bp".format(qual_metrics[1],qual_metrics[0]))
        data_runs[run]["bp_coverage"]=qual_metrics[2]


    # Now I can combine the runs that succeeded
    fq_comb1 = os.path.join(args.temp_dir, args.biosample, args.biosample + "-combined_1.fastq")
    fq_comb2 = os.path.join(args.temp_dir, args.biosample, args.biosample + "-combined_2.fastq")

    #I check that at least one run is OK. I define a variable to count the good runs
    runs_ok=0

    #I create a dictionary: Run => bp_coverage
    ranking_runs={}
    for run in data_runs:
        if(data_runs[run]["flag"]==0):
            runs_ok=runs_ok+1
            ranking_runs[run]=data_runs[run]["bp_coverage"]

    print("- Summaryzing the situation of the sequencing runs")
    # If there are no good runs:
    if runs_ok == 0:
        raise Exception("    - [ERROR] I found some problems in each one of the runs you provided. The analysis stops here!")
    else:
        print("    - {} runs are OK ({})".format(runs_ok,",".join(ranking_runs.keys())))

    selected_runs=[]
    # I find the best runs
    runs_to_consider=args.num_runs_to_consider
    if(runs_to_consider>=len(ranking_runs)):
        runs_to_consider=len(ranking_runs)
    else:
        try:
            runs_to_consider=int(runs_to_consider)
        except:
            raise Exception("    - [ERROR] I have some problems to understand how many runs I should consider for the mapping and the assembly!")
            sys.exit()
    if(len(ranking_runs)>runs_to_consider):
        bp_coverages=sorted(ranking_runs.values(),reverse=True)
        first_two=set(bp_coverages[0:runs_to_consider])
        for run in ranking_runs:
            if ranking_runs[run] in first_two:
                selected_runs.append(run)
    else:
        selected_runs=ranking_runs.keys()

    print("    - Selected runs: {}".format(",".join(selected_runs)))

    # I reinitialize the fastq files.
    cmd="> {}".format(fq_comb1)
    os.system(cmd)
    cmd="> {}".format(fq_comb2)
    os.system(cmd)

    for run in selected_runs:
        trflstem1 = os.path.join(args.temp_dir, args.biosample, run + "-trimmed_1.fastq")
        trflstem2 = os.path.join(args.temp_dir, args.biosample, run + "-trimmed_2.fastq")
        cmd="cat {0} >> {1}".format(trflstem1,fq_comb1)
        os.system(cmd)
        cmd="cat {0} >> {1}".format(trflstem2,fq_comb2)
        os.system(cmd)
