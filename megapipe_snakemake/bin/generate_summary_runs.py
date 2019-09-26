#!/usr/bin/env python

import json
import argparse
import glob

""" for each file I read the info in json format and I get, I get the run IDs and I write down the table """

parser = argparse.ArgumentParser()
parser.add_argument("dir_metadata", type=str,
                    help="directory with the metadata in json format")
parser.add_argument("outfile_summary_runs", type=str,
                    help="output file containing the associations biosample => run_ids")
args = parser.parse_args()

with open(args.outfile_summary_runs, "w") as outf:
    outf.write("biosample\trun_id\n")
    for f in glob.glob(args.dir_metadata + "/*.json"):
        json_content = json.load(open(f))
        biosample=""
        list_of_runs=[]
        for k in json_content:
            list_of_runs.append(json_content[k]['Run'])
            biosample=json_content[k]['BioSample']
        outf.write("{}\t{}\n".format(biosample, ",".join(list_of_runs)))

