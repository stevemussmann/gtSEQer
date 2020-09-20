from __future__ import print_function

from DefaultListOrderedDict import DefaultListOrderedDict

import subprocess
import sys
import multiprocessing
import os.path

class Extractor():
	'Class for pulling genome regions from reference'

	def __init__(self, f, i, e):
		self.targets = f
		self.flank=e
		f=open(i)
		d=f.readlines()
		f.close()
		d=[l.rstrip() for l in d]
		self.indlist = d
		self.coords=list()
		
		#get base of the path where coverage_files will be located
		self.cwd=os.getcwd()
		self.covdir=os.path.join(self.cwd, "coverage_files")
		self.outdir=os.path.join(self.cwd, "extracted_regions")
		if(os.path.isdir(self.outdir) != True):
			os.mkdir(self.outdir)

	def extract(self):
		for chrom, snps in self.targets.items():
			for snp in snps:
				c=list()

				start = int(snp) - self.flank
				end = int(snp) + self.flank

				#make list with coordinates for samtools
				c.append(chrom)
				c.append(":")
				c.append(str(start))
				c.append("-")
				c.append(str(end))

				#concatenate coordinates
				newstr="".join(c)
				self.coords.append(newstr)

	def readFile(self):
		f=open(self.vcf)
		data=f.readlines()
		f.close()
		data = [l.rstrip() for l in data]
		return data

	def makeCommands(self, genome, gz):
		#outfiles=list()
		p=False
		#parallelize here
		self.runCommands(genome, gz)

	def runProgram(self,string):
		try:
			process = subprocess.Popen(string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			output, err = process.communicate()
			print(err)
			if process.returncode != 0:
				print("Non-zero exit status:")
				print(process.returncode)
				raise SystemExit
		except(KeyboardInterrupt, SystemExit):
			raise
		except:
			print("Unexpected error:")
			print(sys.exec_info())
			raise SystemExit

	def runCommands(self, genome, gz):
		for coord in self.coords:
			for ind in self.indlist:
				outfile=coord.replace(":","_")
				outfn=os.path.join(self.outdir,outfile)
				cf=ind + ".cov.txt"
				cfn=os.path.join(self.covdir, cf)
				command = "samtools faidx " + genome + " " + coord + " | bcftools consensus -s " + ind + " -m " + cfn +  " -I " + gz + " >> " + outfn
				#print(command)
				self.runProgram(command)
				
