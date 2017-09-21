import subprocess as sp
import re
import os.path

class GridEngine:
	myGridEngine="unkonwn"
	file_modules="loadmodules.txt"
	
	def __init__(self):
		try:
			out = sp.check_output("which bsub 2>&1", shell=True)
			print(out)
			if bool(re.search("bsub\n$",out)):
				self.myGridEngine="LSF"
		except:
			pass
		try:
			out = sp.check_output("which sbatch 2>&1", shell=True)
			print(out)
			if bool(re.search("sbatch\n$",out)):
				self.myGridEngine="Slurm"
		except:
			pass

	def generate_script(self, sh_file, queue, wall_time, outfile_job, mem_requirements, cmds):

		if(self.myGridEngine=="LSF"):
			with open(sh_file,"w") as outf:
				outf.write("#!/bin/bash\n")
				outf.write("#BSUB -q "+queue+"\n")
				outf.write("#BSUB -W "+wall_time+"\n")
				outf.write("#BSUB -o "+outfile_job+"\n")
				outf.write("#BSUB -R rusage[mem="+mem_requirements+"]\n")
				if os.path.isfile(self.file_modules): 
					with open(self.file_modules, 'r') as inp:
						text = inp.read()	
						modules=re.split(",|\n",text)
						for module in modules:
							if module!="":
								outf.write("module load "+module+"\n")
				outf.write(cmds+"\n")

		elif(self.myGridEngine=="Slurm"):
			with open(sh_file,"w") as outf:
				outf.write("#!/bin/bash\n")
				outf.write("#SBATCH -p "+queue+"\n")
				outf.write("#SBATCH -t "+wall_time+"\n")
				outf.write("#SBATCH -o "+outfile_job+"\n")
				outf.write("#SBATCH --mem "+mem_requirements+"\n")
				if os.path.isfile(self.file_modules): 
					with open(self.file_modules, 'r') as inp:
						text = inp.read()
						modules=re.split(",|\n",text)
						for module in modules:
							if module!="":
								outf.write("module load "+module+"\n")
				outf.write(cmds+"\n")
		else:
			print("I did not recognize your grid engine!\n")

	def launch_job(self, sh_script):
		if(self.myGridEngine=="LSF"):
			cmd="bsub < "+sh_script
			sp.call(cmd,shell=True)
		elif(self.myGridEngine=="Slurm"):
			cmd="sbatch "+sh_script
			sp.call(cmd,shell=True)
		else:
			print("I did not recognize your grid engine!\n")
