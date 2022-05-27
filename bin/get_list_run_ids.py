#!/usr/bin/env python

import json
import pandas as pd

with open('./config/config_pipeline.json') as inp:
    config = json.loads(inp.read())

with open(config["logs_analysis"]+"isolates_to_analyze.txt", "r") as inp:
    samples = inp.read().splitlines()

# I get the data about my strain
print("- I am retrieving the data about the runs to analyze from {}".format(config["summary_runs"]))
tab=pd.read_csv(config["summary_runs"], sep="\t")
tab_sel=tab.loc[tab["biosample"].isin(samples)]
list_run_ids=[val for sublist in tab_sel["run_id"].str.split(",") for val in sublist]

with open(config["logs_analysis"]+"runs_to_download.txt","w") as outf:
    for run_id in list_run_ids:
        outf.write("{}\n".format(run_id))
