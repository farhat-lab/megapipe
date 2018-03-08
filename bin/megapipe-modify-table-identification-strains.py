#!/usr/bin/env python

import sys
import subprocess as sbp

def modifyTableIdenfication(table_identification,new_metadata):
    # I get the data from the first table
    fields=[]
    data={}
    with open(table_identification,"r") as inp:
        fields=inp.readline().rstrip("\n").split("\t")
        for line in inp:
            if line =="\n":
                continue
            entry=line.rstrip("\n").split("\t")
            try:
                idx_internal_xref=fields.index("internal_xref")
                internal_xref=entry[idx_internal_xref]
                data[internal_xref]={}
                if internal_xref=="":
                    print("[ERROR] the following entry does not have an internal_xref:")
                    print(line.rstrip("\n"))
                    sys.exit()
            except:
                print("[ERROR] the following entry does not have an internal_xref:")
                print(line.rstrip("\n"))
                sys.exit()
            for idx,item in enumerate(entry):
                data[internal_xref][fields[idx]]={}
                data[internal_xref][fields[idx]]=item
    # I store the fields of the first file
    new_fields=fields
    # Now I read the new data and I update the fields
    with open(new_metadata,"r") as inp2:
        # I add the fields of the second file and then I find the intersection
        fields=inp2.readline().rstrip("\n").split("\t")
        new_fields.extend(fields)
        new_fields=set(new_fields)
        for line in inp2:
            if line =="\n":
                continue
            entry=line.rstrip("\n").split("\t")
            try:
                idx_internal_xref=fields.index("internal_xref")
                internal_xref=entry[idx_internal_xref]
                data[internal_xref]={}
                if internal_xref=="":
                    print("[ERROR] the following entry does not have an internal_xref:")
                    print(line.rstrip("\n"))
                    sys.exit()
            except:
                print("[ERROR] the following entry does not have an internal_xref:")
                print(line.rstrip("\n"))
                sys.exit()
            for idx,item in enumerate(entry):
                if item=="":
                    continue
                try:
                    data[internal_xref][fields[idx]]=item
                except:
                    data[internal_xref][fields[idx]]={}
                    data[internal_xref][fields[idx]]=item
                    pass
            # I set up the public_xref if the biosample is defined
            try:
                idx_biosample=fields.index("biosample")
                idx_public_xref=fields.index("public_xref")
                idx_runs=fields.index("runs")
                public_xref=entry[idx_public_xref]
                biosample=entry[idx_biosample]
                runs=entry[idx_runs]
                if(public_xref=="") and (biosample!=""):
                    data[internal_xref]["public_xref"]={}
                    data[internal_xref]["public_xref"]=biosample
                if(public_xref=="") and (biosample=="") and (runs!=""):
                    data_runs=runs.split(",")
                    if(len(data_runs)==1):
                        data[internal_xref]["public_xref"]={}
                        data[internal_xref]["public_xref"]=runs
            except:
                pass
    with open(table_identification+".mod","w") as outf:
        # I write down the fields
        outf.write("\t".join(new_fields)+"\n")
        # I write down the data
        for key in data:
            newRow=[]
            for field in new_fields:
                if field in data[key]:
                    newRow.append(data[key][field])
                else:
                    newRow.append("")
            outf.write("\t".join(newRow)+"\n")

def executeChanges(fileIn):
    print("[INFO] Ready to commit the changes. You can still review the changes. I created a temporary file: {}".format(fileIn+".mod"))
    response = "none"
    while(response not in ["y","n"]):
        response=input("Should I comit the changes? [y|n]: ")
    if(response=="y"):
        cmd="mv {0}.mod {0}".format(fileIn)
        sbp.check_output(cmd,shell=True)
    else:
        print("[INFO] Aborting!")


if len(sys.argv) != 3:
	print("::usage: {} <table_identification_strains> <table_new_metadata>".format(sys.argv[0]))
	sys.exit()

modifyTableIdenfication(sys.argv[1],sys.argv[2])
executeChanges(sys.argv[1])
