#!/usr/bin/env python
# Akshith Doddi 07/10/2017 (original script), additions by Luca Freschi
# --------------------------------------------------------------------------------------------
# The purpose of this script is to launch or megapipe.py on a large number of fastq files.

import os
import sys
from sys import argv


if(len(argv) != 7):
	print("Usage: {} <file_in> <out_dir> <scratch_dir> <startindex> <endindex> <runmode>".format(argv[0]))
	sys.exit()


file_in=argv[1]
out_dir=argv[2]
scratch_dir=argv[3]
starti = int(argv[4]) - 1
endi = int(argv[5]) - 1

inputer = open(file_in)
info = inputer.read().split("\n")
inputer.close()


if(argv[-1] == "0"):
	
	if not os.path.exists(out_dir+"/commands"):
		os.makedirs(out_dir+"/commands")




	for i in range(starti, endi + 1):
		details=info[i].split("\t")
		outputer = open(out_dir+"/commands/"+"job_{}.sh".format(details[0]), "w")
		outputer.write("#!/bin/bash\n")
		outputer.write("#BSUB -q short\n")
		outputer.write("#BSUB -W 12:00\n")
		outputer.write("#BSUB -o "+ out_dir +"/{0}/{1}.job\n".format(details[0],details[0]))
		outputer.write("#BSUB -R rusage[mem=35000]\n")
		outputer.write("module load dev/python/3.4.2\n")
		outputer.write("module load dev/perl/5.18.1\n")
		outputer.write("module load seq/samtools/1.3\n")
		outputer.write("module load seq/bwa/0.7.8\n")
		outputer.write("module load dev/java/jdk1.8\n")
		outputer.write("megapipe.py {} {} {} {} {}\n".format(details[0], details[1], details[2], scratch_dir, out_dir))
		outputer.close()


elif(argv[-1] == "1"):
	for i in range(starti, endi + 1):
		details=info[i].split("\t")
		os.system("bsub <" + out_dir + "/commands/" + "job_{}.sh".format(details[0]))

else:
	print("Invalid runmode: " + argv[3])


