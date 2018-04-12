#!/usr/bin/env python

import os
from os import listdir
import sys
from sys import argv
from gridmanager import gridpuppeteer as gp

if(len(argv) != 7):
    print("::usage: {} <table_identification_strains> <fastq_dir> <output_dir> <scratch_dir> <jobs_to_launch> <dir_logs>".format(argv[0]))
    sys.exit()

#I initialize the variables
table=argv[1]
fastq_dir=argv[2]
dir_results=argv[3]
scratch_dir=argv[4]
num_jobs_to_launch=argv[5]
#directory with the scripts for slurm and the logs
dir_logs=argv[6]

grid_obj=gp.GridEngine()
if not os.path.exists(dir_logs):
    os.makedirs(dir_logs+"/commands")

# I get the data about the strains that have already been analyzed
already_done=listdir(dir_results+"/")
# I determine which strains I have to analyze
to_analyze=[]
with open(table,"r") as inp:
    fields=inp.readline().rstrip("\n").split("\t")
    if "public_xref" in fields:
        idx_public_xref=fields.index("public_xref")
    if "internal_xref" in fields:
        idx_internal_xref=fields.index("internal_xref")
    for line in inp:
        if line =="\n":
            continue
        entry=line.rstrip("\n").split("\t")
        if("idx_public_xref" in vars()):
            public_xref=entry[idx_public_xref]
            if((public_xref!="") and (public_xref not in already_done)):
                to_analyze.append(public_xref)
        if("idx_internal_xref" in vars()):
            internal_xref=entry[idx_internal_xref]
            if((internal_xref!="") and (internal_xref not in already_done)):
                to_analyze.append(internal_xref)
        if ((not "idx_public_xref" in vars()) and (not "idx_internal_xref" in vars()))
                print("[ERROR] An entry does not have neither a public_xref nor an internal_xref")
                print(line.rstrip("\n"))

num_jobs_to_launch=int(num_jobs_to_launch)
if(num_jobs_to_launch > len(to_analyze)): 
    num_jobs_to_launch=len(to_analyze)
for i in range(0, num_jobs_to_launch):
        script_name=dir_logs+"/commands/"+"job_{}.sh".format(to_analyze[i])
        out_file=dir_logs +"/{0}.job".format(to_analyze[i])
        cmd="megapipe-core.py {} {} {} {} {} {} all\n".format(to_analyze[i], table, fastq_dir, dir_results, scratch_dir,dir_logs)
        grid_obj.generate_script(script_name, "short", "12:00:00", out_file, "35G", cmd)
for i in range(0, num_jobs_to_launch):
        script_name=dir_logs + "/commands/" + "job_{}.sh".format(to_analyze[i])
        grid_obj.launch_job(script_name)


