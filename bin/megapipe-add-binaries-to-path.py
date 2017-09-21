#!/usr/bin/env python

import os
from os.path import expanduser

home = expanduser("~")


with open(home+"/.bashrc","a") as outf:
	print("Adding MEGAPIPE binaries to your path (will work on Linux only)!")
	outf.write("#MEGAPIPE\n")
	cwd = os.getcwd()
	outf.write("PATH={}:$PATH\n".format(cwd))
