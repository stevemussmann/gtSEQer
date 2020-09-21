from __future__ import print_function

import os.path

class MakeTestDir():
	'Class for creating dir, checking if dir exists'

	def __init__(self, d):
		self.cwd=os.getcwd()
		self.newdir=os.path.join(self.cwd, d)

	def testDir(self):
		if(os.path.isdir(self.newdir) != True):
			os.mkdir(self.newdir)
		return self.newdir

