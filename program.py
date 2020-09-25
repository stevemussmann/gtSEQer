from __future__ import print_function

import subprocess
import sys

class Program():
	'Class for executing a system call'

	def __init__(self, s):
		self.string = s

	def runProgram(self):
		try:
			process = subprocess.Popen(self.string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			output, err = process.communicate()
			print(err)
			if process.returncode != 0:
				print("Non-zero exit status:")
				print(process.returncode)
				raise SystemExit
			return output

		except(KeyboardInterrupt, SystemExit):
			raise
		except:
			print("Unexpected error:")
			print(sys.exc_info())
			raise SystemExit
