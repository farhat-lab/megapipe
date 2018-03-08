#!/usr/bin/env python
# I import the required libraries
import sys
from os import path,makedirs,chmod
import subprocess as sbp
from glob import glob
#I check that all parameters are there
if(len(sys.argv) != 4):
    print("::usage: {} <filein> <dest_dir> <dir_logs_grid_engine>".format(sys.argv[0]))
    sys.exit()
fileIn=sys.argv[1]
destDir=sys.argv[2]
dirLogsFromGridEngine=sys.argv[3]
if not path.exists(destDir):
    print("[INFO] I create the directory that will store the fastq files %s" % destDir)
    makedirs(destDir)
if not path.exists(dirLogsFromGridEngine):
    print("[INFO] I create the directory that will store the grid engine logs %s" % dirLogsFromGridEngine)
    makedirs(dirLogsFromGridEngine)
print("[INFO] I generate the .sh files that will be submitted to o2")
with open(fileIn,"r") as inp:
    for line in inp:
        idNCBI=line.rstrip("\n")
        if idNCBI=="":
            continue
        with open(dirLogsFromGridEngine+"/"+idNCBI+"_cmds.sh","w") as outf:
            outf.write("#!/bin/bash\n")
            outf.write("path=`srapath \'{}\'`\n".format(idNCBI))
            outf.write("fastq-dump --split-files --gzip $path -O {}\n".format(destDir))
print("[INFO] I submit the jobs to o2")
for script in glob(dirLogsFromGridEngine+"/"+"*_cmds.sh"):
    cmd="sbatch -p short -n 1 -t 10:00 --mem 50M -o {0}.out -e {0}.err {0}".format(script)
    print("Command: "+cmd)
    sbp.call(cmd,shell=True)
print("Bye!")

