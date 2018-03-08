#!/usr/bin/env python

import os
import sys
from sys import argv
from gridmanager import gridpuppeteer as gp

if(len(argv) != 7):
	print("::usage: {} <file_in> <out_dir> <scratch_dir> <startindex> <endindex> <runmode>".format(argv[0]))
	sys.exit()

#I initialize the variables
file_in=argv[1]
out_dir=argv[2]
scratch_dir=argv[3]
starti = int(argv[4]) - 1
endi = int(argv[5]) - 1

inputer = open(file_in)
info = inputer.read().split("\n")
inputer.close()

grid_obj=gp.GridEngine()
#If runmode==0 (I generate the scripts)
if(argv[-1] == "0"):
	if not os.path.exists(out_dir+"/commands"):
		os.makedirs(out_dir+"/commands")
	for i in range(starti, endi + 1):
		details=info[i].split("\t")
		script_name=out_dir+"/commands/"+"job_{}.sh".format(details[0])
		out_file=out_dir +"/{0}/{0}.job".format(details[0])
		cmd="megapipe-core.py {} {} {} {} {}\n".format(details[0], details[1], details[2], scratch_dir, out_dir)
		grid_obj.generate_script(script_name, "short", "12:00", out_file, "35000", cmd)
#If runmode==1 (I launch the scripts)
elif(argv[-1] == "1"):
	for i in range(starti, endi + 1):
		details=info[i].split("\t")
		script_name=out_dir + "/commands/" + "job_{}.sh".format(details[0])
		grid_obj.launch_job(script_name)
else:
	print("Invalid runmode: " + argv[3])


