from __future__ import print_function

from Bio import SeqIO

import os.path

class Genome():
	'Class for operating on fasta files'

	def __init__(self, f):
		self.seqLengths = dict()
		fn="genome_region_sizes.txt"
		if os.path.isfile(fn) != True:
			print("******************************************************")
			print("Chromosome sizes not present. Calculating from genome.")
			print("******************************************************")
			self.seqLengths = self.readFasta(f)
			self.writeChromSizes(fn)
		else:
			print("******************************************************")
			print("Chromosome sizes present. Loading file", fn)
			print("******************************************************")
			self.seqLengths = self.readChromSizes(fn)
		#print(self.seqLengths)
		

	def readFasta(self, f):
		sequences = dict()
		seqs = SeqIO.parse(open(f), 'fasta')
		for seq_record in seqs:
			sequence = str(seq_record.seq).upper()
			seq_len = int(len(sequence))
			sequences[seq_record.id] = seq_len
		return sequences

	def writeChromSizes(self, fn):
		with open(fn, 'w') as fh:
			for key, val in self.seqLengths.items():
				fh.write(key)
				fh.write("\t")
				fh.write(str(val))
				fh.write("\n")

	def readChromSizes(self, fn):
		chroms=dict()
		fh=open(fn)
		data=fh.readlines()
		fh.close()
		data=[l.strip() for l in data]
		for line in data:
			content=line.split()
			chroms[content[0]] = content[1]
		return chroms


