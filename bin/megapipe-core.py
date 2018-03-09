#!/usr/bin/env python

import sys
from os import popen, system, listdir, getcwd, path, environ, makedirs
import re
import subprocess as sbp
import json
from glob import glob
import pkg_resources

def detect_cfg_file():
    if path.isfile(environ['HOME']+"/.megapipe.json"):
        with open(environ['HOME']+"/.megapipe.json") as json_file:
            data = json.load(json_file)
            #If you update the pipeline, do not forget to update this check
            if(("picard" in data) and ("pilon" in data) and ("prinseq" in data) and ("kraken_db" in data) and ("qualimap" in data) and ("spades" in data) and ("fasta_ref" in data)):
                return(data)
            else:
                print("::I could not find some required variables in the configuration file")
    else:
        print("::I do not have a configuration file. Please check the documentation to understand how to create one!")
        sys.exit()

def detect_weird_read_names(fastq):
    with open(fastq,"r") as inf:
        fline=inf.readline()
        decision=False
        if re.search("#0/[1-9]\n$",fline):
            decision=True
    return decision

def valFQ(file1,file2):
    print("::validate fastq files")
    try:
        out1 = sbp.check_output("fastQValidator --file " + file1, shell=True)
        out2 = sbp.check_output("fastQValidator --file " + file2, shell=True)
        m1 = re.search("FASTQ_SUCCESS", out1)
        m2 = re.search("FASTQ_SUCCESS", out2)
        if(m1 == None or m2 == None):
            return False
        return True
    except:
        return False

def write_msg(file_log,msg):
    if path.exists(file_log):
        mode = 'a'
    else:
        mode = 'w'
    with open(file_log,mode) as outf:
        outf.write(msg+"\n")

def checkDRs(depths):
    """The purpose of this method is to check that at least 95% of base positions in the drug resistance conferring regions of the genome are covered 10x. This method was heavily inspired by and based upon Dr. Farhat's script "qc_report.py" in the work-horse directory"""
    okaybases = 0
    for d in depths:
        if(d >= 10):
            okaybases += 1
    okayp = okaybases / len(depths)
    return okayp

########
##MAIN##
########
    
if(len(sys.argv) != 7):
    print("::usage: {} <tag> <table_identification_strains> <dir_fastq> <output_dir> <scratch_dir> <dir_logs>".format(sys.argv[0]))
    sys.exit()

tag=sys.argv[1]
table=sys.argv[2]
dir_fastq=sys.argv[3]
out_dir=sys.argv[4]
scratch_dir=sys.argv[5]
log_dir=sys.argv[6]
#I read the configuration file
data_json=detect_cfg_file()

# the log file will the the following
file_log=log_dir+"/{}.out".format(tag)

# I log the version of megapipe I am using
write_msg(file_log,"[INFO] I am using megapipe v{}".format(pkg_resources.get_distribution("megapipe").version))

# I create new directories (if they do not exist)
if not path.exists(out_dir):
    makedirs(out_dir)
if not path.exists(scratch_dir):
    makedirs(scratch_dir)
if not path.exists(out_dir+"/"+tag+"/"):
    makedirs(out_dir+"/"+tag+"/")
if not path.exists(log_dir):
    makedirs(log_dir)

#I check if the reference is ok. If so, I index it (if needed).
if path.isfile(data_json["fasta_ref"]):
    if path.isfile(data_json["fasta_ref"]+".bwt"):
        write_msg(file_log,"--->I detected the index of the reference. I do not generate a new index")
    else:
        write_msg(file_log,"--->I did not detect the index of the reference. I will generate it")
        cmd="bwa index "+data_json["fasta_ref"]
        system(cmd)

# I get the data about my strain
runs_to_analyze=[]
with open(table,"r") as inp:
    fields=inp.readline().rstrip("\n").split("\t")
    for line in inp:
        if line =="\n":
            continue
        entry=line.rstrip("\n").split("\t")
        try:
            idx_public_xref=fields.index("public_xref")
            public_xref=entry[idx_public_xref]
            has_public_xref=True
        except:
            has_public_xref=False
        try:
            idx_internal_xref=fields.index("internal_xref")
            internal_xref=entry[idx_internal_xref]
            has_internal_xref=True
        except:
            has_internal_xref=False
        if not (has_public_xref or has_internal_xref):
            print("[ERROR] I did not find neither the public_xref nor the internal_xref")
            sys.exit()
        if "public_xref" in vars():
            if public_xref==tag:
                # I get the runs
                idx_runs=fields.index("runs")
                runs=entry[idx_runs].split(",")
                # I reformat the information
                for run in runs:
                    current_run_formatted=run+":"+dir_fastq+"/"+run+"_1.fastq.gz,"+dir_fastq+"/"+run+"_2.fastq.gz"
                    runs_to_analyze.append(current_run_formatted)
        if "internal_xref" in vars():
            if internal_xref==tag:
                # I get the runs
                idx_runs=fields.index("internal_fastq_files")
                internal_runs=list(entry[idx_runs])
                runs_to_analyze.extend(internal_runs)

# Now foreach run I unzip it, I check that the fastq files are ok.

#I unzip the fastq files in a directory on scratch2
write_msg(file_log,"--->Unzipping fastq files (on {})".format(scratch_dir))

# I will use this dictionary to  the runs (0=everything is fine, 1=there is a problem, so the run will be excluded)
data_runs={}

for current_run in runs_to_analyze:

    # I get the info about this run
    general_info=current_run.split(":")
    run=general_info[0]
    fastq_files=general_info[1].split(",")

    # the flag is 0 for now.
    data_runs[run]={}
    data_runs[run]["flag"]={}
    data_runs[run]["fq1"]={}
    data_runs[run]["fq1"]=fastq_files[0]
    data_runs[run]["fq2"]={}
    data_runs[run]["fq2"]=fastq_files[1]

    # I unzip the fastq files of this run
    #I unzip the fastq files in a directory on scratch2
    write_msg(file_log,"----->Unzipping fastq files for run {0} (on {1})".format(general_info[0],scratch_dir))
    cmd="zcat {0} > {1}/{2}_1.fastq".format(fastq_files[0], scratch_dir, general_info[0])
    print(cmd)
    system(cmd)
    cmd="zcat {0} > {1}/{2}_2.fastq".format(fastq_files[1], scratch_dir, general_info[0])
    print(cmd)
    system(cmd)

    #I rename the variables fqf1 and fqf2 since the fastq I will use are the ones on the scratch, so that I do not touch the original ones
    fqf1 = scratch_dir +"/" + general_info[0] + "_1.fastq"
    fqf2 = scratch_dir +"/" + general_info[0] + "_2.fastq"

    #I check the fasta files
    write_msg(file_log,"----->Checking fastq files for run {0}".format(general_info[0]))
    try:
        valFQ(fqf1,fqf2)
        write_msg(file_log,"[INFO] Fastq files are valid for run {}".format(general_info[0]))
    except:
        write_msg(file_log,"[WARNING] Fastq files are NOT valid for run {}".format(general_info[0]))
        data_runs[run]["flag"]=1
        continue

    # I check the names of the reads
    write_msg(file_log,"--->Checking the names of the reads for run {}".format(general_info[0]))
    test_names1=detect_weird_read_names(fqf1)
    test_names2=detect_weird_read_names(fqf2)
    if((test_names1==True) or (test_names2==True)):
        write_msg(file_log,"[INFO] I found some weird names. I am fixing them!")
        cmd="megapipe-check-names.py {}".format(fqf1)
        system(cmd)
        cmd="megapipe-check-names.py {}".format(fqf2)
        system(cmd)
    else:
        write_msg(file_log,"[INFO] No weird names!")

    # I trim the reads with Prinseq -- this should happen in the scratch in order to save space
    write_msg(file_log,"----->Trimming with Prinseq (run {})".format(general_info[0]))
    path_to_prinseq=data_json["prinseq"]
    cmd="perl "+ path_to_prinseq +" -fastq {0} -fastq2 {1} -out_format 3 -out_good {2}/{3}-trimmed -out_bad null -log {4}/{3}-prinseq.log -min_qual_mean 20 -verbose".format(fqf1, fqf2, scratch_dir, run,log_dir)
    print(cmd)
    system(cmd)
    write_msg(file_log,"[INFO] Please see the Prinseq report in {0}/{1}-prinseq.log".format(log_dir, run))
    # I check if the trimming went well.
    trflstem1 = run + "-trimmed_1.fastq"
    trflstem2 = run + "-trimmed_2.fastq"
    trfl1 = scratch_dir + "/" + trflstem1
    trfl2 = scratch_dir + "/" + trflstem2
    files = glob(scratch_dir+"/"+run+"-trimmed*")
    if(not(trflstem1 in files and trflstem2 in files)):
        write_msg(file_log, "Prinseq failed.")
        data_runs[run]["flag"]=1
        continue
    #sctrfl = scratchtrimmedfile // I am aware these two lines are not useful. Please remove them.
    sctrfl1 = scratch_dir +"/" + trflstem1
    sctrfl2 = scratch_dir + "/" + trflstem2

    # I classify the reads with Kraken
    write_msg(file_log,"----->Classifying reads with Kraken")
    path_to_krakendb=data_json["kraken_db"]
    cmd="kraken --fastq-input {0} --output {1} --db {2}".format(sctrfl1, scratch_dir + "/" + trflstem1 + ".krkn", path_to_krakendb)
    e=sbp.check_output(cmd,shell=True)
    mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
    tbperc1 = float(mat.groups()[1]) / 100
    cmd="kraken --fastq-input {0} --output {1} --db {2}".format(sctrfl2, scratch_dir + "/" + trflstem2 + ".krkn", path_to_krakendb)
    e=sbp.check_output(cmd,shell=True)
    mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
    tbperc2 = float(mat.groups()[1]) / 100
    write_msg(file_log,"[INFO] Please see the Kraken report in {0}/{1}.krkn and {0}/{2}.krkn".format(out_dir, trflstem1, trflstem2))
    # We still need to check that more than 90% of the sequences are from mycobacterium tuberculosis by making sense of output from Kraken
    # write doen in the log tbperc1 and 2
    if(tbperc1 < 0.9):
        write_msg(file_log,"Less than 90% of reads in the first fastq file belonged to Mycobacterium tuberculosis")
        data_runs[run]["flag"]=1
        continue
    if(tbperc2 < 0.9):
        write_msg(file_log,"Less than 90% of reads in the second fastq file belonged to Mycobacterium tuberculosis")
        data_runs[run]["flag"]=1
        continue

# Now I can combine the runs that succeeded
fq_comb1=scratch_dir + "/" + tag+"-combined_1.fastq"
fq_comb2=scratch_dir + "/" + tag+"-combined_2.fastq"
for run in data_runs:
    if(data_runs[run]["flag"]==0):
        trflstem1 = scratch_dir+ "/" + run + "-trimmed_1.fastq"
        trflstem2 = scratch_dir+ "/" + run + "-trimmed_2.fastq"
        cmd="cat {} >> {}".format(trflstem1,fq_comb1)
        system(cmd)
        cmd="cat {} >> {}".format(trflstem2,fq_comb2)
        system(cmd)

# Aligning reads to the reference
write_msg(file_log,"--->Aligning reads with bwa")
samfile = scratch_dir + "/{}.sam".format(tag)
#piper = popen("bwa mem -M -R '@RG\tID:<unknown>\tSM:<unknown>\tPL:<unknown>\tLB:<unknown>\tPU:<unknown>' RefGen/TBRefGen.fasta {0} {1} > {2}".format(sctrfl1, sctrfl2, samfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
out=sbp.check_output("bwa mem -M {0} {1} {2} > {3}".format(data_json["fasta_ref"],fq_comb1, fq_comb2, samfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
write_msg(file_log,out)

# Sorting and converting to bam
write_msg(file_log,"--->Sorting sam file with Picard")
bamfile = scratch_dir + "/{}.sorted.bam".format(tag)
path_to_picard=data_json["picard"]
out=sbp.check_output("java -Xmx16G -jar " + path_to_picard + " SortSam INPUT={0} OUTPUT={1} SORT_ORDER=coordinate".format(samfile, bamfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
write_msg(file_log,out)

# Quality control
write_msg(file_log,"--->Evaluating mapping with Qualimap")
path_to_qualimap=data_json["qualimap"]
system(path_to_qualimap + " bamqc -bam {0} --outfile {1}.pdf --outformat PDF".format(bamfile, tag))
qdir = scratch_dir + "/" + bamfile[bamfile.rindex("/") + 1:-4] + "_stats"
system("mv {0} {1}".format(bamfile.replace(".bam", "_stats"), qdir))
write_msg(file_log,"[INFO] Please see the Qualimap report in {}".format(qdir))

# Removing duplicates
write_msg(file_log,"--->Removing duplicates from bam file with Picard")
drbamfile = bamfile.replace(".bam", ".duprem.bam")
out=sbp.check_output("java -Xmx32G -jar " + path_to_picard +" MarkDuplicates I={0} O={1} REMOVE_DUPLICATES=true M={2} ASSUME_SORT_ORDER=coordinate".format(bamfile, drbamfile, drbamfile[:-4])) # help from http://seqanswers.com/forums/showthread.php?t=13192
write_msg(file_log,out)
write_msg(file_log,"[INFO] Please see the Picard MarkDuplicates report in {0}/{1}".format(scratch_dir,drbamfile[drbamfile.rindex("/") + 1:-4]))

write_msg(file_log,"--->Calculating depth")
out = sbp.check_output("samtools depth -a " + drbamfile)
dmat = re.findall("(.+\t)([0-9]+)(\n)", out)
depths = [float(t[1]) for t in dmat]

#gencov = sum(depths) / len(depths)
#outputer.write("Genome coverage: {}\n".format(gencov))
gencovprop = checkDRs(depths)
write_msg(file_log,"[INFO] The percent of H37Rv bases that have a coverage of at least 10x is {}.".format(str(gencovprop * 100)))
if(gencovprop < 0.95): # The threshold is 95%
	write_msg(file_log,"[PROBLEM] The percent of H37Rv bases that have a coverage of at least 10x is less than 95%.")
	sys.exit()
fileout_depth=results_dir+"/"+tag+"/{}.depth".format(tag)
write_msg(fileout_depth,dout)
refcov = 0
for d in depths:
	if(d > 0):
		refcov += 1
write_msg(file_log,"[INFO] Percent of reference genome covered: {}".format(refcov / len(depths)))

write_msg(file_log,"--->Indexing {}".format(drbamfile))
system("samtools index {}".format(drbamfile))

# I call the variants with pilon
write_msg(file_log,"--->Variant calling with Pilon")
path_to_pilon=data_json["pilon"]
out=sbp.check_output("java -Xmx32G -jar " + path_to_pilon +" --genome RefGen/TBRefGen.fasta --bam {0} --output {1}/{2}/{3} --variant".format(drbamfile, scratch_dir, tag, drbamfile[drbamfile.rindex("/") + 1:-4]))
write_msg("{0}/{1}/{1}-pilon.log".format(scratch_dir,tag),out)

# I generate the assembly with spades
write_msg(file_log,"--->Generating the assembly with Spades")
#-t (treads); -m (memory, in Gb)
path_to_spades=data_json["spades"]
sbp.check_output("python2 " + path_to_spades + " -t 1 -m 30 --careful --pe1-1 {0} --pe1-2 {1} -o {2}/{3}/spades".format(fq_comb1,fq_comb2,scratch_dir,tag))

