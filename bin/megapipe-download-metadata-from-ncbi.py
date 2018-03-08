#!/usr/bin/env python

import subprocess as sp
from bs4 import BeautifulSoup
import sys

def getData(cmd):
    results=sp.check_output(cmd,shell=True)
    return(results)

def getEntryFromRunIDng(run_id):
    bioproject=None
    biosample=None
    submitter=None
    collectResults=[]
    collectResults.append("query_run:{}".format(run_id))
    results_query1=getData("esearch -db sra -query \"{0}\" | efetch -format docsum -mode xml".format(run_id))
    soup=BeautifulSoup(results_query1,"lxml")
    try:
        bioproject=soup.html.body.documentsummaryset.documentsummary.bioproject.string
        collectResults.append("bioproject:{}".format(bioproject))
    except:
        print("[INFO] {}|I have found some problems getting the bioproject".format(run_id))
    try:
        biosample=soup.html.body.documentsummaryset.documentsummary.biosample.string
        collectResults.append("biosample:{}".format(biosample))
    except:
        print("[ERR] {}|I have found some problems getting the biosample".format(run_id))
        return()
    results_query2=getData("esearch -db sra -query \"{0}\" | efetch -format docsum -mode xml".format(biosample))
    soupRuns=BeautifulSoup(results_query2,"lxml")
    try:
        submitter=soupRuns.documentsummaryset.documentsummary.submitter['center_name']
        collectResults.append("submitter:{}".format(submitter))
    except:
        print("[INFO] {}|I have found some problems getting the bioproject".format(run_id))
    runs=soupRuns.find_all('run')
    listOfRuns=[]
    for r in range(0,len(runs)):
        run=runs[r]['acc']
        listOfRuns.append(run)
    collectResults.append("runs:{}".format(",".join(listOfRuns)))
    #I get the attributes associated to this biosample
    results_query3=getData("esearch -db biosample -query \"{0}\" | efetch -format docsum -mode xml".format(biosample))
    soupBiosample=BeautifulSoup(results_query3,"lxml")
    try:
        identifiers=soupBiosample.html.body.identifiers.string.replace(" ","").split(";")
        collectResults.append(identifiers[1])
        collectResults.append(identifiers[2])
    except:
        print("[INFO] {}|I have found some problems getting the attributes of the biosample".format(run_id))
        pass
    attributes=soupBiosample.find_all('attribute')
    for entry in attributes:
        try:
            attribute_name=entry['harmonized_name']
            attribute_value=entry.string
            collectResults.append("{0}:{1}".format(attribute_name,attribute_value))
        except:
            pass
    return("|".join(collectResults))

def getDataNCBIWriteFile(fileIn,fileOut):
    failedRuns=[]
    with open(fileOut+".txt","w") as outf:
        with open(fileIn,"r") as inp:
            for line in inp:
                entry=line.rstrip("\n")
                if entry=="":
                    continue
                print("[INFO] Downloading NCBI data. Query_run: {}".format(entry),)
                try:
                    result=getEntryFromRunIDng(entry)
                    outf.write(str(result)+"\n")
                except:
                    failedRuns.append(entry)
    if len(failedRuns)>0:
        print("[INFO] failed query_run(s): {}".format(",".join(failedRuns)))

def convertDataToTable(fileIn,fileOut):
    fields={}
    db={}
    with open(fileIn,"r") as inp:
        with open(fileOut,"w") as outf:
            for line in inp:
                entry=line.rstrip("\n").split("|")
                key=entry[0].split(":")[1]
                if key not in db:
                    db[key]={}
                for element in entry:
                    current=element.split(":")
                    if current[0] not in db[key]:
                        db[key][current[0]]={}
                        current_string=":".join(current[1:])
                        db[key][current[0]]=current_string
                    if current[0] not in fields:
                        fields[current[0]]={}
            # I write down the results
            outf.write("\t".join(sorted(fields))+"\n")
            for key in db:
                newRow=[]
                for field in sorted(fields.keys()):
                    if db[key][field]:
                        newRow.append(db[key][field])
                    else:
                        newRow.append("")
                outf.write("\t".join(newRow)+"\n")

if len(sys.argv) != 3:
	print("::usage: {} <csv_in> <tag_out_file>".format(sys.argv[0]))
	sys.exit()

getDataNCBIWriteFile(sys.argv[1],sys.argv[2])
convertDataToTable(sys.argv[2]+".txt",sys.argv[2]+"_tab.txt")
