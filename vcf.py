from __future__ import print_function

from DefaultListOrderedDict import DefaultListOrderedDict

class VCF():
	'Class for operating on VCF files'

	def __init__(self, f):
		self.vcf = f
		print(self.vcf)
		self.data=self.readFile()

	def getSNPs(self):
		snps=DefaultListOrderedDict()
		for line in self.data:
			if not line.startswith("#"):
				items = line.split()
				snps[items[0]].append(items[1])
		return snps

	def readFile(self):
		f=open(self.vcf)
		data=f.readlines()
		f.close()
		data = [l.rstrip() for l in data]
		return data
