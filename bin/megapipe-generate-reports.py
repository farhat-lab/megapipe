#!/usr/bin/env python

import glob
import sys
import os.path

if(len(sys.argv) != 2):
    print(":: usage: {} <results_dir> ".format(sys.argv[0]))
    sys.exit()

res_dir=sys.argv[1]

l_dirs=glob.glob(res_dir+"/*")

print("xref\tlineage\tvcf\tassembly")

for d in sorted(l_dirs):
    vcf="no"
    assembly="no"
    lineage="no"
    d2=d.split("/")[-1]
    #I get the lineage
    if os.path.isfile(res_dir+"/"+d2+"/fast-lineage-caller/"+d2+".lineage"):
        with open(res_dir+"/"+d2+"/fast-lineage-caller/"+d2+".lineage","r") as inp:
            line=inp.readline()
            entry=line.rstrip("\n").split(" ")
            lineage=entry[1]
    #I check if the vcf is present
    if os.path.isfile(res_dir+"/"+d2+"/pilon/"+d2+".vcf"):
        vcf="YES"
    #I check if the assembly is present
    if os.path.isfile(res_dir+"/"+d2+"/spades/"+d2+"_scaffolds.fasta"):
        assembly="YES"
    print("\t".join([d2,lineage,vcf,assembly]))


