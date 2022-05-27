#!/usr/bin/env python

import re
import sys
import subprocess as sp

if len(sys.argv) != 3:
    print("::usage: %s <file_in.vcf> <file_out.vcf> " % sys.argv[0])
    sys.exit()

print("--reducing size of vcf file %s" % sys.argv[1])

with open(sys.argv[2],"w") as outf:
    with open(sys.argv[1],"r") as inp:
        for line in inp:

            #skip the comment lines
            if line.startswith("#"):
                outf.write(line)
                continue
            data=line.rstrip("\n").split("\t")
            #if ALT is "." and the REF has only one base -> skip it
            if ((len(data[3])==1) and (data[4]==".")):
                continue
            else:
                outf.write(line)


