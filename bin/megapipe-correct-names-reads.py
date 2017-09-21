#!/usr/bin/env python

import sys
import re
from itertools import islice
import os

if(len(sys.argv) != 2):
	print("::usage: {} <fastq_file>".format(sys.argv[0]))
	sys.exit()

fastq_file=sys.argv[1]

print("--cleaning file "+fastq_file)
with open(fastq_file,"r") as inf:
	with open(fastq_file+"_clean","w") as outf:
		while True:
			next_n_lines = list(islice(inf, 4))
			if not next_n_lines:
				break
			first_line=next_n_lines[0]
			if re.search("#0/[1-9]\n$",first_line):
				arr=first_line.split("#")
				next_n_lines[0]="#".join(arr[0:-1])
				next_n_lines[0]=next_n_lines[0]+"\n"
			#I write down these 4 lines
			for line in next_n_lines:
				outf.write(line)

cmd="mv "+fastq_file+"_clean "+fastq_file
#print(cmd)
os.system(cmd)







