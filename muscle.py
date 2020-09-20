from __future__ import print_function

import os.path
import subprocess
import sys

class Muscle():
	'Class for running Muscle on fasta files'

	def __init__(self,f):
		self.fas = f
		basename=os.path.basename(f)
		print("Aligning", basename)

		self.cwd=os.getcwd()
		self.aligned=os.path.join(self.cwd, "muscle_aligned")
		if(os.path.isdir(self.aligned) != True):
			os.mkdir(self.aligned)

		command = self.makeCommand(basename)
		print(command)
		self.runProgram(command)

	def makeCommand(self, f):
		out=os.path.join(self.aligned, f)
		string = "muscle -in " + self.fas + " -out " + out
		return string

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
