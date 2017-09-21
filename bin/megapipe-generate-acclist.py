#!/usr/bin/env python

import glob
import sys

if(len(sys.argv) != 4):
	print("::usage: {} <dir_fastq> <tag_pe1> <ext>".format(sys.argv[0]))
	sys.exit()

fastq_dir=sys.argv[1]
tag_pe1=sys.argv[2]
ext=sys.argv[3]

l_files=glob.glob(fastq_dir+"/*"+ext)

for item in l_files:
	dirs=item.split("/")
	current=dirs[-1]
	tag=current.replace(ext,"").replace("_","").replace("-","")
	fastq1=item
	tag_pe2=tag_pe1.replace("1","2")
	fastq2=item.replace(tag_pe1,tag_pe2)
	final_list=[tag,fastq1,fastq2]
	print("\t".join(final_list))


