#!/usr/bin/env python

import glob
import sys
import os.path

if(len(sys.argv) != 2):
	print("::usage: {} <results_dir> ".format(sys.argv[0]))
	sys.exit()

res_dir=sys.argv[1]

l_dirs=glob.glob(res_dir+"/*")



print("#tag\tvcf_present\tassembly_present")

for d in sorted(l_dirs):

	vcf="NO"
	assembly="NO"

	d2=d.split("/")[-1]
	#I check if the vcf is present
	if os.path.isfile(res_dir+"/"+d2+"/"+d2+".sorted.duprem.vcf"):
		vcf="yes"

	#I check if the assembly is present
	if os.path.isfile(res_dir+"/"+d2+"/spades/scaffolds.fasta"):
		assembly="yes"

	print("\t".join([d2,vcf,assembly]))


