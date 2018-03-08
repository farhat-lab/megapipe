#!/usr/bin/env python

import sys
from os import path,makedirs
import subprocess as sbp

if(len(sys.argv) != 4):
	print("::usage: {} <filein> <dest_dir> <threads>".format(sys.argv[0]))
	sys.exit()

fileIn=sys.argv[1]
destDir=sys.argv[2]
numThreads=sys.argv[3]



if not path.exists(destDir):
	print("::I create the directory %s" % destDir)
	makedirs(destDir)

print("::I generate the job file for parallel")
with open("megapipe-temp-job-file.txt","w") as outf:
	with open(fileIn,"r") as inp:
		for line in inp:
			idNCBI=line.rstrip("\n")
			if idNCBI=="":
				continue
			cmd="srapath %s" % idNCBI
			pathID=sbp.check_output(cmd,shell=True)
			pathID=pathID.decode('utf-8').rstrip("\n")
			cmd="fastq-dump --split-files --gzip %s -O %s" % (pathID,destDir)
			outf.write(cmd+"\n")

print("::I download the fastq files using parallel -- this might take some time! (Threads=%s)" % numThreads)
cmd="parallel -n %s < megapipe-temp-job-file.txt" % numThreads

sbp.call(cmd,shell=True)
sbp.call("rm -rf megapipe-temp-job-file.txt",shell=True)
print("Bye!")

