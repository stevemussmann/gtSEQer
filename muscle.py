from __future__ import print_function

import os.path

class Muscle():
	'Class for running Muscle on fasta files'

	def __init__(self,f):
		self.fas = f
		print(self.fas)

		self.cwd=os.getcwd()
		self.aligned=os.path.join(self.cwd, "muscle_aligned")
		if(os.path.isdir(self.aligned) != True):
			os.mkdir(self.aligned)

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
