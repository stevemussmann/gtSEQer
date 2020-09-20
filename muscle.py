from __future__ import print_function

import os.path
import subprocess
import sys

from Bio import AlignIO
from Bio.Align import AlignInfo

class Muscle():
	'Class for running Muscle on fasta files'

	def __init__(self,f):
		self.fas = f
		basename=os.path.basename(f)

		print("Aligning", basename)

		self.cwd=os.getcwd()
		self.aligned=os.path.join(self.cwd, "muscle_aligned")
		self.consensus=os.path.join(self.cwd, "muscle_aligned_consensus")
		if(os.path.isdir(self.aligned) != True):
			os.mkdir(self.aligned)
		if(os.path.isdir(self.consensus) != True):
			os.mkdir(self.consensus)

		self.out=os.path.join(self.aligned, basename)

		command = self.makeCommand()
		print(command)
		self.runProgram(command)
		self.makeConsensus(basename)

	def makeCommand(self):
		string = "muscle -in " + self.fas + " -out " + self.out
		return string

	def makeConsensus(self, bn):
		fn=os.path.join(self.consensus, bn)
		alignment = AlignIO.read(self.out, 'fasta')
		summary_align = AlignInfo.SummaryInfo(alignment)
		con = summary_align.gap_consensus(0.99, 'N')

		#print(con)

		f=open(fn, "w")
		f.write(">")
		f.write(bn)
		f.write("\n")
		f.write(str(con))
		f.write("\n")

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
