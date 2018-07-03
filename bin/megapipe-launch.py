#!/usr/bin/env python

import os
from os import listdir
import argparse
from gridmanager import gridpuppeteer as gp

parser = argparse.ArgumentParser()

parser.add_argument("table_identification_strains", type=str,
                    help="Path to the table identification strains")
parser.add_argument("fastq_dir", type=str,
                    help="Path to the directory containing the fastq files")
parser.add_argument("output_dir", type=str,
                    help="Path to the directory where the results will be written")
parser.add_argument("scratch_dir", type=str,
                    help="Path to the directory where the temporary files will be written")
parser.add_argument("dir_logs", type=str,
                    help="Path to the directory where the slurm scripts and logs will be written")
parser.add_argument("num_jobs_to_launch", type=int,
                    help="Number of jobs to launch")
parser.add_argument("-k", "--keep_tmp", action="store_true",
                    help="keeps the temporary files")
parser.add_argument("-s", "--skip_assembly", action="store_true",
                    help="skips the assembly ")
args = parser.parse_args()


#I initialize the variables
table = args.table_identification_strains
fastq_dir = args.fastq_dir
dir_results = args.output_dir
scratch_dir = args.scratch_dir
num_jobs_to_launch = args.num_jobs_to_launch
#directory with the scripts for slurm and the logs
dir_logs = args.dir_logs

grid_obj = gp.GridEngine()
if not os.path.exists(dir_logs):
    os.makedirs(dir_logs+"/commands")

# I get the data about the strains that have already been analyzed
try: 
    already_done = listdir(dir_results+"/")
except:
    already_done = []
    pass
# I determine which strains I have to analyze
to_analyze = []
col_public_xref=True
col_internal_xref=True
with open(table, "r") as inp:
    fields = inp.readline().rstrip("\n").split("\t")
    try:
        idx_public_xref = fields.index("public_xref")
    except:
        col_public_xref=False
    try:
        idx_internal_xref = fields.index("internal_xref")
    except:
        col_internal_xref=False
    for line in inp:
        if line == "\n":
            continue
        entry = line.rstrip("\n").split("\t")
        if col_public_xref :
            public_xref = entry[idx_public_xref]
            if((public_xref != "") and (public_xref not in already_done)):
                to_analyze.append(public_xref)
        if col_internal_xref :
            internal_xref=entry[idx_internal_xref]
            if((internal_xref != "") and (internal_xref not in already_done)):
                to_analyze.append(internal_xref)
        if ((not "idx_public_xref" in vars()) and (not "idx_internal_xref" in vars())):
                print("[ERROR] An entry does not have neither a public_xref nor an internal_xref")
                print(line.rstrip("\n"))

num_jobs_to_launch = int(num_jobs_to_launch)
if(num_jobs_to_launch > len(to_analyze)): 
    num_jobs_to_launch = len(to_analyze)
for i in range(0, num_jobs_to_launch):
    script_name = dir_logs+"/commands/"+"job_{}.sh".format(to_analyze[i])
    out_file = dir_logs +"/{0}.job".format(to_analyze[i])
    cmd = "megapipe-core.py {} {} {} {} {} {} 100".format(to_analyze[i], table, fastq_dir, dir_results, scratch_dir, dir_logs)
    if args.skip_assembly:
        cmd+=" --skip_assembly"
    if args.keep_tmp:
        cmd+=" --keep_tmp"
    grid_obj.generate_script(script_name, "short", "12:00:00", out_file, "35G", cmd+"\n")
for i in range(0, num_jobs_to_launch):
    script_name = dir_logs + "/commands/" + "job_{}.sh".format(to_analyze[i])
    grid_obj.launch_job(script_name)


