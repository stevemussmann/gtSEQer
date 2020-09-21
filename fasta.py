from __future__ import print_function

import os.path
from Bio import SeqIO

from maketestdir import MakeTestDir

class Fasta():
	'Class for operating on fasta files'

	def __init__(self, f):
		self.fasta = f
		basename=os.path.basename(f)
		print("Filtering", basename)
		self.seqs = self.readFasta(self.fasta)
		#print(self.seqs)
		
		#make output path if it doesn't exist
		mtd = MakeTestDir("filtered_sequences")
		self.filt=mtd.testDir()
		#self.cwd=os.getcwd()
		#self.filt=os.path.join(self.cwd, "filtered_sequences")
		#if(os.path.isdir(self.filt) != True):
		#	os.mkdir(self.filt)

		self.printSeqs(basename)
		

	def readFasta(self, f):
		per_n=55
		sequences = dict()
		seqs = SeqIO.parse(open(f), 'fasta')
		for seq_record in seqs:
			sequence = str(seq_record.seq).upper()
			if( (float(sequence.count("N")) / float(len(sequence)) ) * 100 <= per_n):
				if sequence not in sequences:
					sequences[sequence] = seq_record.id
		return sequences

	def printSeqs(self, bn):
		fn=os.path.join(self.filt, bn)
		f=open(fn, "w")
		for seq, name in self.seqs.items():
			if( name.startswith(">") != True ):
				f.write(">")
			f.write(name)
			f.write("\n")
			f.write(seq)
			f.write("\n")
		f.close()
