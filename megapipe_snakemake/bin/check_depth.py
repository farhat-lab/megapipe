#!/usr/bin/env python3

import subprocess as abp
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("depth_file", type=str,
                    help="bam file duplicated reads removed")
args = parser.parse_args()


def checkDRs(depths):
    """The purpose of this method is to check that at least 95% of base positions in the drug resistance conferring regions of the genome are covered 10x. This method was heavily inspired by and based upon Dr. Farhat's script "qc_report.py" in the work-horse directory"""
    okaybases = 0
    for d in depths:
        if(d >= 10):
            okaybases += 1
    okayp = okaybases / len(depths)
    return okayp


# I detect the tag
base=os.path.basename(args.depth_file)
tag=base.replace(".depth","")

depths=[]
print(":: Calculating depth")
with open(args.depth_file,"r") as inp:
    out = inp.read()
    dmat = re.findall("(.+\t)([0-9]+)(\n)", out)
    depths = [float(t[1]) for t in dmat]

#gencov = sum(depths) / len(depths)
#outputer.write("Genome coverage: {}\n".format(gencov))
gencovprop = checkDRs(depths)
print("  * The percent of H37Rv bases that have a coverage of at least 10x is {}.".format(str(gencovprop * 100)))
if(gencovprop < 0.95): # The threshold is 95%
    print("    - [ERROR] The percent of H37Rv bases that have a coverage of at least 10x is less than 95%.")
else:
    with open("results/{0}/depth/{0}_depth_OK".format(tag),"w") as outf:
        outf.write("OK\n")

refcov = 0
for d in depths:
    if(d > 0):
        refcov += 1
print("  * Percent of reference genome covered: {}".format(refcov / len(depths)))

