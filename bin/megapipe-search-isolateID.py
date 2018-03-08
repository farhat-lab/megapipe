#!/usr/bin/env python

import sys
import os

if(len(sys.argv) != 4):
	print("::usage: {} <dir_results> <txt_filein> <tag_outfiles>".format(sys.argv[0]))
	sys.exit()

dirRes=sys.argv[1]
fileIDs=sys.argv[2]
tag=sys.argv[3]

with open(tag+"_found.txt","w") as outf1:
	with open(tag+"_notFound.txt","w") as outf2:
		with open(fileIDs,"r") as inp:
			for line in inp:
				idIsol=line.rstrip("\n")
				if os.path.isdir(dirRes+"/"+idIsol):
					outf1.write(line)
				else:
					outf2.write(line)
