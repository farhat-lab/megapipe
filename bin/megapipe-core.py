#!/usr/bin/env python

import sys
from os import popen, system, listdir, getcwd, path, environ, makedirs
import re
import subprocess as sbp
import json


def detect_cfg_file():
	if path.isfile(environ['HOME']+"/.megapipe.json"):
		with open(environ['HOME']+"/.megapipe.json") as json_file:
			data = json.load(json_file)
			#If you update the pipeline, do not forget to update this check
			if(("picard" in data) and ("pilon" in data) and ("prinseq" in data) and ("kraken_db" in data) and ("qualimap" in data) and ("spades" in data) and ("fasta_ref" in data)):
				return(data)
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
	out1 = sbp.check_output("fastQValidator --file " + file1, shell=True)
	out2 = sbp.check_output("fastQValidator --file " + file2, shell=True)
	m1 = re.search("FASTQ_SUCCESS", out1)
	m2 = re.search("FASTQ_SUCCESS", out2)
	if(m1 == None or m2 == None):
		return False
	return True


def write_msg(file_log,msg):
	with open(file_log,"a") as outf:
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

if(len(sys.argv) != 6):
	print("::usage: {} <tag> <fastq_file1> <fastq_file2> <scratch_dir> <output_dir>".format(sys.argv[0]))
	sys.exit()

tag=sys.argv[1]
fqf1 = sys.argv[2]
fqf2 = sys.argv[3]
scratch_dir=sys.argv[4]
out_dir=sys.argv[5]

#I read the configuration file
data_json=detect_cfg_file()

file_log=out_dir+"/"+tag+"/{}.out".format(tag)

if not path.exists(out_dir):
	makedirs(out_dir)
if not path.exists(scratch_dir):
	makedirs(scratch_dir)
if not path.exists(out_dir+"/"+tag+"/"):
	makedirs(out_dir+"/"+tag+"/")

#I check if the reference is ok. If so, I index it (if needed).
if path.isfile(data_json["fasta_ref"]):
	if path.isfile(data_json["fasta_ref"]+".bwt"):
		write_msg(file_log,"--->I detected the index of the reference. I do not generate a new index")
	else:
		write_msg(file_log,"--->I did not detect the index of the reference. I will generate it")
		cmd="bwa index "+data_json["fasta_ref"]
		system(cmd)

#I unzip the fastq files in a directory on scratch2
write_msg(file_log,"--->Unzipping fastq files (on scratch2)")
cmd="zcat {0} > {1}/{2}_1.fastq".format(fqf1, scratch_dir, tag)
print(cmd)
system(cmd)
cmd="zcat {0} > {1}/{2}_2.fastq".format(fqf2, scratch_dir, tag)
print(cmd)
system(cmd)

#I rename the variables fqf1 and fqf2 since the fastq I will use are the ones on the scratch, so that I do not touch the original ones
fqf1 = scratch_dir +"/" + tag + "_1.fastq"
fqf2 = scratch_dir +"/" + tag + "_2.fastq"


#I check the fasta files
write_msg(file_log,"--->Checking fastq files")

if(valFQ(fqf1,fqf2)):
	write_msg(file_log,"[INFO] Fastq files are valid!")
else:
	write_msg(file_log,"[PROBLEM] Not valid fastq files.")
	write_msg(file_log,"Exiting now...\n")
	sys.exit()


write_msg(file_log,"--->Checking the names of the reads")

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



write_msg(file_log,"--->Trimming with Prinseq")
path_to_prinseq=data_json["prinseq"]
cmd="perl "+ path_to_prinseq +" -fastq {0} -fastq2 {1} -out_format 3 -out_good {2}/{3}/{3}-trimmed -out_bad null -log {2}/{3}/{3}-trimmed.log -min_qual_mean 20 -verbose".format(fqf1, fqf2, out_dir, tag)
print(cmd)
system(cmd)

write_msg(file_log,"[INFO] Please see the Prinseq report in {0}/{1}/{1}-trimmed.log".format(out_dir, tag))

trflstem1 = tag + "-trimmed_1.fastq"
trflstem2 = tag + "-trimmed_2.fastq"
trfl1 = out_dir + "/" + tag + "/" + trflstem1
trfl2 = out_dir + "/" + tag + "/" + trflstem2

#cmd="cp {0} {1}".format(fqf1,trfl1)
#system(cmd)
#cmd="cp {0} {1}".format(fqf2,trfl2)
#system(cmd)

files = listdir(out_dir+"/"+tag+"/")
if(not(trflstem1 in files and trflstem2 in files)):
	write_msg(file_log, "Prinseq failed.")
	write_msg(file_log, "Exiting now...")
	sys.exit()

write_msg(file_log,"--->Moving trimmed fastq files to scratch")
cmd="mv {0} {1}/".format(trfl1,scratch_dir)
print(cmd)
system(cmd)
cmd="mv {0} {1}/".format(trfl2,scratch_dir)
print(cmd)
system(cmd)


#sctrfl = scratchtrimmedfile
sctrfl1 = scratch_dir +"/" + trflstem1
sctrfl2 = scratch_dir + "/" + trflstem2


write_msg(file_log,"--->Classifying reads with Kraken")
path_to_krakendb=data_json["kraken_db"]
piper = sbp.Popen(["kraken", "--fastq-input", sctrfl1, "--output", out_dir + "/" + tag + "/" + trflstem1 + ".krkn", "--db", path_to_krakendb], stderr = sbp.PIPE, stdout = sbp.PIPE) # help from https://docs.python.org/2/library/subprocess.html#subprocess.Popen and Mr. Martin Owens
(o, e) = piper.communicate()
mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
tbperc1 = float(mat.groups()[1]) / 100

piper = sbp.Popen(["kraken", "--fastq-input", sctrfl2, "--output", out_dir + "/" + tag + "/" + trflstem2 + ".krkn", "--db", path_to_krakendb], stderr = sbp.PIPE, stdout = sbp.PIPE) # help from https://docs.python.org/2/library/subprocess.html#subprocess.Popen and Mr. Martin Owens
(o, e) = piper.communicate()
mat = re.search("( classified \()([0-9]+\.*[0-9]*)(%\))", str(e))
tbperc2 = float(mat.groups()[1]) / 100

write_msg(file_log,"[INFO] Please see the Kraken report in {0}/{1}/{2}.krkn and {0}/{1}/{3}.krkn".format(out_dir, tag, trflstem1, trflstem2))

# We still need to check that more than 90% of the sequences are from mycobacterium tuberculosis by making sense of output from Kraken
if(tbperc1 < 0.9):
	write_msg(file_log,"Less than 90% of reads in the first fastq file belonged to mycobacterium tuberculosis")
	write_msg(file_log,"Exiting now...")
	sys.exit()

if(tbperc2 < 0.9):
	write_msg(file_log,"Less than 90% of reads in the second fastq file belonged to mycobacterium tuberculosis")
	write_msg(file_log,"Exiting now...")
	sys.exit()

write_msg(file_log,"--->Aligning reads with bwa")
samfile = scratch_dir + "/{}.sam".format(tag)
#piper = popen("bwa mem -M -R '@RG\tID:<unknown>\tSM:<unknown>\tPL:<unknown>\tLB:<unknown>\tPU:<unknown>' RefGen/TBRefGen.fasta {0} {1} > {2}".format(sctrfl1, sctrfl2, samfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
piper = popen("bwa mem -M {0} {1} {2} > {3}".format(data_json["fasta_ref"],sctrfl1, sctrfl2, samfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
write_msg(file_log,piper.read())
piper.close()

write_msg(file_log,"--->Sorting sam file with Picard")
bamfile = scratch_dir + "/{}.sorted.bam".format(tag)
path_to_picard=data_json["picard"]
piper = popen("java -Xmx16G -jar " + path_to_picard + " SortSam INPUT={0} OUTPUT={1} SORT_ORDER=coordinate".format(samfile, bamfile)) # help from http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
write_msg(file_log,piper.read())
piper.close()

#piper = popen("samtools view -b -o {0} {1}".format(bamfile, samfile))
#piper.close()
#outputer.write("--->Converting to bam file done.\n")

write_msg(file_log,"--->Evaluating mapping with Qualimap")
path_to_qualimap=data_json["qualimap"]
system(path_to_qualimap + " bamqc -bam {0} --outfile {1}.pdf --outformat PDF".format(bamfile, tag))
qdir = out_dir + "/" + tag+ "/" + bamfile[bamfile.rindex("/") + 1:-4] + "_stats"
system("mv {0} {1}".format(bamfile.replace(".bam", "_stats"), qdir))
write_msg(file_log,"[INFO] Please see the Qualimap report in {}".format(qdir))

write_msg(file_log,"--->Removing duplicates from bam file with Picard")
drbamfile = bamfile.replace(".bam", ".duprem.bam")
piper = popen("java -Xmx32G -jar " + path_to_picard +" MarkDuplicates I={0} O={1} REMOVE_DUPLICATES=true M={2} ASSUME_SORT_ORDER=coordinate".format(bamfile, drbamfile, drbamfile[:-4])) # help from http://seqanswers.com/forums/showthread.php?t=13192
write_msg(file_log,piper.read())
piper.close()
write_msg(file_log,"[INFO] Please see the Picard MarkDuplicates report in {0}/{1}".format(scratch_dir,drbamfile[drbamfile.rindex("/") + 1:-4]))

write_msg(file_log,"--->Calculating depth")
piper = popen("samtools depth -a " + drbamfile)
dout = piper.read()
piper.close()
dmat = re.findall("(.+\t)([0-9]+)(\n)", dout)
depths = [float(t[1]) for t in dmat]

#gencov = sum(depths) / len(depths)
#outputer.write("Genome coverage: {}\n".format(gencov))
gencovprop = checkDRs(depths)
write_msg(file_log,"[INFO] The percent of H37Rv bases that have a coverage of at least 10x is {}.".format(str(gencovprop * 100)))
if(gencovprop < 0.95): # The threshold is 95%
	write_msg(file_log,"[PROBLEM] The percent of H37Rv bases that have a coverage of at least 10x is less than 95%.")
	write_msg(file_log,"Exiting now...\n")
	sys.exit()
fileout_depth=out_dir+"/"+tag+"/{}.depth".format(tag)
write_msg(fileout_depth,dout)
refcov = 0
for d in depths:
	if(d > 0):
		refcov += 1
write_msg(file_log,"[INFO] Percent of reference genome covered: {}".format(refcov / len(depths)))

write_msg(file_log,"--->Indexing {}".format(drbamfile))
system("samtools index {}".format(drbamfile))

write_msg(file_log,"--->Variant calling with Pilon")
path_to_pilon=data_json["pilon"]
piper = popen("java -Xmx32G -jar " + path_to_pilon +" --genome RefGen/TBRefGen.fasta --bam {0} --output {1}/{2}/{3} --variant".format(drbamfile, out_dir, tag, drbamfile[drbamfile.rindex("/") + 1:-4]))
write_msg("{0}/{1}/{1}-pilon.log".format(out_dir,tag),piper.read())
piper.close()


write_msg(file_log,"--->Generating the assembly with Spades")
#-t (treads); -m (memory, in Gb)
path_to_spades=data_json["spades"]
piper = popen("python2 " + path_to_spades + " -t 1 -m 30 --careful --pe1-1 {0} --pe1-2 {1} -o {2}/{3}/spades".format(sctrfl1,sctrfl2,out_dir,tag))
piper.close()





