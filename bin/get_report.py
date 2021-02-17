#!/usr/bin/env python3

import json
import os

# I read the isolates to analyze file and check the directories that are present in the results folder
## Need to load the json config file
with open('./config/config.json') as json_file:
    config_data = json.load(json_file)

## get the list of isolates
with open(config_data["logs_analysis"] + "isolates_to_analyze.txt") as inp:
    isolates_to_analyze = inp.read().splitlines()

tot_isolates=len(isolates_to_analyze)

# I can check the steps that have a log, and compare them to the expected logs to understand where something went wrong.
## example: TBDM101265. Should be fine. Why it does not have a VCF?
list_all_steps_pipeline_ordered = ["combine_runs.txt", "align_to_ref.txt", "sort_convert_tobam.txt", "duprem.txt", "calc_depth.txt", "indexing_bam.txt", "variant_calling.txt", "lineage_calling.txt"]


steps_megapipe_failed={}
isolates_with_vcf={}
#check if the isolates have the log files
for isolate in isolates_to_analyze:
    #I check if the VCF is there
    if os.path.exists("results/"+isolate+"/pilon/"+isolate+".vcf"):
        isolates_with_vcf[isolate]=1
    else:
        #get the log files
        current_list_files = os.listdir(config_data["logs_analysis"]+"/"+isolate)

        for x, file_name in enumerate(list_all_steps_pipeline_ordered):
            if file_name not in current_list_files:
                if list_all_steps_pipeline_ordered[x-1] not in steps_megapipe_failed:
                    steps_megapipe_failed[list_all_steps_pipeline_ordered[x-1]]={}
                if isolate not in steps_megapipe_failed[list_all_steps_pipeline_ordered[x-1]]:
                    steps_megapipe_failed[list_all_steps_pipeline_ordered[x-1]][isolate]={}
                break


print("* SUCCEDED: {} ({:.2f}%);FAILED: {} ({:.2f}%)".format(len(isolates_with_vcf.keys()), len(isolates_with_vcf.keys())/tot_isolates*100, tot_isolates-len(isolates_with_vcf.keys()),  (tot_isolates-len(isolates_with_vcf.keys()))/tot_isolates*100))
print("* Steps where megapipe stopped:")
for step in steps_megapipe_failed:
    print("  - {}".format(step))
    print(list(steps_megapipe_failed[step].keys()))

