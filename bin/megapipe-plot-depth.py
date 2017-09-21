#!/usr/bin/env python

import subprocess as sp
import sys
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as mpl

if(len(sys.argv) != 2):
	print("::usage: {} <bam> ".format(sys.argv[0]))
	sys.exit()


threshold_coverage=10

bam=sys.argv[1]

cmd="samtools depth -a {}".format(bam)

out = sp.getoutput(cmd)

x=[]
y=[]

n_bases_grater_thr=0

for entry in out.split("\n"):
	fragments=entry.split("\t")
	x.append(fragments[1])
	y.append(float(fragments[2]))

	if int(fragments[2]) >= threshold_coverage:
		n_bases_grater_thr=n_bases_grater_thr+1

perc_bases_grater_thr=float(n_bases_grater_thr)/float(len(y))*100


print("::percentage of bases that have a coverage > {0}: {1}".format(threshold_coverage,perc_bases_grater_thr))
print("::generating figure -- coverage for each base...")
f=mpl.figure()
mpl.plot(x,y)
f.savefig(bam+"_coverage-genomes.pdf")
print("::generating figure -- histogram coverage...")
f2=mpl.figure()
mpl.hist(y,bins=10)
f2.savefig(bam+"_hist.pdf")

print("::Bye!")

