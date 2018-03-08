#!/usr/bin/env python

import sys

def generateTableIdenficationFromScratch(fileIn,fileOut):
    fields=[]
    ids={}
    with open(fileIn,"r") as inp:
        with open(fileOut,"w") as outf:
            fields=inp.readline().rstrip("\n").split("\t")
            fields.extend(["public_xref","internal_xref"])
            outf.write("\t".join(fields)+"\n")
            for line in inp:
                entry=line.rstrip("\n").split("\t")
                # I create two new empty fields for "public_xref" and "internal_xref"
                entry.extend(["",""])
                # I get the index of the biosample field
                idx=fields.index("biosample")
                if entry[idx]!="":
                    if entry[idx] in ids:
                        print("[ERROR] Duplicate id detected! Skipping {}".format(entry[idx]))
                        continue
                    entry[-2]=entry[idx]
                    entry[-1]=entry[idx]
                else:
                    idx=fields.index("runs")
                    runs=entry[idx].split(",")
                    if len(runs)==1:
                        if entry[idx] in ids:
                            print("[ERROR] Duplicate id detected! Skipping {}".format(entry[idx]))
                            continue
                        entry[-2]=entry[idx]
                        entry[-1]=entry[idx]
                    else:
                        idx=fields.index("query_run")
                        print("[ERROR] Skipping {}".format(entry[idx]))
                        continue
                outf.write("\t".join(entry)+"\n")
                ids[entry[idx]]={}


if len(sys.argv) != 3:
	print("::usage: {} <table_metadata_NCBI> <table_identification_strains>".format(sys.argv[0]))
	sys.exit()

generateTableIdenficationFromScratch(sys.argv[1],sys.argv[2])

