#!/usr/bin/env python

import subprocess as sbp
import shutil as shu
from os import environ, path
import json

def checkInstalled(humanReadableName,cmd):
    if shu.which(cmd) is None:
        print("--%s (NO) -- please install it!" % humanReadableName)
    else:
        print("--%s (ok)" % humanReadableName)


def checkPath(name,p):
    if not path.isfile(p):
        print("--you did not provide a good path for %s. Please check the configuration file \"~/.megapipe.json\"" % ())

def checkMegapipeCfg():
    if path.isfile(environ['HOME']+"/.megapipe.json"):
        with open(environ['HOME']+"/.megapipe.json","r") as json_file:
            data = json.load(json_file)
            #If you update the pipeline, do not forget to update this check
            progToCheck=["picard","pilon","prinseq","qualimap","spades"]
            for prog in progToCheck:
                checkPath(prog,data[prog])
    else:
        print("::I do not have a configuration file (~/.megapipe.json). Please check the documentation to understand how to create one!")
        sys.exit()



#dependencies: kraken, fastQValidator,
#to add pilon, prinseq, picard, spades -- cfg file



checkInstalled("Kraken","kraken")
checkInstalled("FastQValidator","fastQValidator")
checkInstalled("BWA","bwa")
checkInstalled("Samtools","samtools")
checkMegapipeCfg()


