from __future__ import print_function

import os.path
import pandas

class Coverage():
	'Class for identifying genome regions with insufficient coverage'

	def __init__(self, d, l, m):
		self.directory = d
		self.l = l
		self.extension = m

		self.indlist = self.readFile(self.l)

		#make sure output directory exists for writing
		cwd=os.getcwd()
		self.outdir = os.path.join(cwd,"coverage_files")
		if(os.path.isdir(self.outdir) != True ):
			os.mkdir(self.outdir)

	def applyMinCov(self, c, overwrite):
		for f in self.indlist:
			fn=f + ".cov.txt"
			out=os.path.join(self.outdir,fn)	

			#only apply MinCov for files that do not already exist
			if(overwrite==True):
				self.pandasFunctions(f, out, c)
			elif(os.path.isfile(out) != True):
				self.pandasFunctions(f, out, c)

	def readFile(self, fn):
		f=open(fn)
		data=f.readlines()
		f.close()
		data = [l.rstrip() for l in data]
		return data

	def pandasFunctions(self, f, out, c):
		fn = f + "." + self.extension
		fn=os.path.join(self.directory, fn)
		print(fn)
		#read into pandas dataframe
		df=pandas.read_csv(fn, index_col=0,  header=None, sep='\t', compression='gzip')
		#remove values with sequencing depth >= coverage
		df=df[~(df[3]>=c)]

		#keep only index plus columns 1-2
		df=df.iloc[:, : 2]
				
		#convert to 1-based numbering
		df=df+1
		print(df)

		#print df to file
		df.to_csv(out, sep='\t', index=True, header=False)

