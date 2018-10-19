#!/usr/bin/env python3

import sys
import os.path
from os import scandir
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("results_dir", type=str,
					help="directory with the megapipe results")
parser.add_argument("tsv_out", type=str,
					help="output file (TSV)")
args = parser.parse_args()

d_isolates={}
l_dirs = os.scandir(args.results_dir)
with open(args.tsv_out,"w") as outf:
	outf.write("xref\tlineage\tvcf\tassembly\tbam\n")
	for d in l_dirs:
		vcf="no"
		assembly="no"
		lineage="no"
		bam="no"
		#I get the lineage
		if os.path.isfile(args.results_dir+"/"+d.name+"/fast-lineage-caller/"+d.name+".lineage"):
			with open(args.results_dir+"/"+d.name+"/fast-lineage-caller/"+d.name+".lineage","r") as inp:
				line=inp.readline()
				entry=line.rstrip("\n").split(" ")
				lineage=entry[1]
		#I check if the vcf is present
		if os.path.isfile(args.results_dir+"/"+d.name+"/pilon/"+d.name+".vcf"):
			vcf="YES"
		#I check if the bam is present
		if os.path.isfile(args.results_dir+"/"+d.name+"/bam/"+d.name+".sorted.duprem.bam"):
			bam="YES"
		#I check if the assembly is present
		if os.path.isfile(args.results_dir+"/"+d.name+"/spades/"+d.name+"_scaffolds.fasta"):
			assembly="YES"
		d_isolates[d.name]="\t".join([d.name, lineage, vcf, assembly, bam])+"\n"
	for isolate in sorted(d_isolates):
		outf.write(d_isolates[isolate])
