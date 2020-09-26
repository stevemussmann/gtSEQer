from __future__ import print_function

import os.path
import subprocess
import sys

from Bio import AlignIO
from Bio.Align import AlignInfo

from program import Program
from maketestdir import MakeTestDir

class Muscle():
	'Class for running Muscle on fasta files'

	def __init__(self,f, dirname):
		self.fas = f
		basename=os.path.basename(f)

		print("Aligning", basename)

		mtd_align = MakeTestDir(dirname)
		dirnamecon = dirname + "_consensus"
		mtd_consensus = MakeTestDir(dirnamecon)
		self.aligned=mtd_align.testDir()
		self.consensus=mtd_consensus.testDir()

		self.out=os.path.join(self.aligned, basename)

		command = self.makeCommand()
		prog = Program(command)
		prog.runProgram()
		print(command)
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

