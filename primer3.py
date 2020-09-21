from __future__ import print_function

import os.path
import subprocess
import sys

from program import Program
from maketestdir import MakeTestDir

class Primer3():
	'Class for executing primer3'

	def __init__(self, fn):
		# make directories
		mtd_mac = MakeTestDir("muscle_aligned_consensus")
		mac = mtd_mac.testDir()
		self.fas=os.path.join(mac, fn)

		mtd_p3in = MakeTestDir("p3in")
		mtd_p3out = MakeTestDir("p3out")
		self.p3inDir=mtd_p3in.testDir()
		self.p3outDir=mtd_p3out.testDir()
		
		# parse file
		seqinfo = self.parseFile()

		# make p3in file
		self.makeFile(seqinfo["ID"], seqinfo["SEQUENCE"])

		# run primer3_core
		comm = self.makeCommand(fn)
		#print(comm)

		print("Finding primers for", fn)
		prog = Program(comm)
		prog.runProgram()

	def makeFile(self, name, seq):
		fn=name + ".p3in.txt"
		fn=os.path.join(self.p3inDir, fn)
		
		f = open(fn, 'w')
		f.write("SEQUENCE_ID=") #finish this line
		f.write(name)
		f.write("\n")
		f.write("SEQUENCE_TEMPLATE=") #finish
		f.write(seq)
		f.write("\n")
		f.write("SEQUENCE_TARGET=150,3") #add "start coord, length"
		f.write("\n")
		f.write("PRIMER_TASK=pick_detection_primers")
		f.write("\n")
		f.write("PRIMER_PICK_LEFT_PRIMER=1")
		f.write("\n")
		f.write("PRIMER_PICK_INTERNAL_OLIGO=0")
		f.write("\n")
		f.write("PRIMER_PICK_RIGHT_PRIMER=1")
		f.write("\n")
		f.write("PRIMER_NUM_RETURN=1")
		f.write("\n")
		f.write("PRIMER_OPT_SIZE=22")
		f.write("\n")
		f.write("PRIMER_MIN_SIZE=18")
		f.write("\n")
		f.write("PRIMER_MAX_SIZE=24")
		f.write("\n")
		f.write("PRIMER_MAX_NS_ACCEPTED=0")
		f.write("\n")
		f.write("PRIMER_PRODUCT_SIZE_RANGE=75-125")
		f.write("\n")
		f.write("PRIMER_PRODUCT_OPT_SIZE=100")
		f.write("\n")
		f.write("PRIMER_GC_CLAMP=1")
		f.write("\n")
		f.write("PRIMER_MIN_TM=58")
		f.write("\n")
		f.write("PRIMER_OPT_TM=61")
		f.write("\n")
		f.write("PRIMER_MAX_TM=64")
		f.write("\n")
		f.write("PRIMER_PAIR_MAX_DIFF_TM=2")
		f.write("\n")
		f.write("PRIMER_MAX_POLY_X=3")
		f.write("\n")
		f.write("PRIMER_MAX_HAIRPIN_TH=20")
		f.write("\n")
		f.write("PRIMER_MAX_SELF_ANY_TH=20")
		f.write("\n")
		f.write("P3_FILE_FLAG=1")
		f.write("\n")
		f.write("SEQUENCE_INTERNAL_EXCLUDED_REGION=150,3") #add "start coord, length"
		f.write("\n")
		f.write("PRIMER_OPT_GC_PERCENT=50")
		f.write("\n")
		f.write("PRIMER_MIN_GC_PERCENT=40")
		f.write("\n")
		f.write("PRIMER_MAX_GC_PERCENT=60")
		f.write("\n")
		f.write("PRIMER_EXPLAIN_FLAG=1")
		f.write("\n")
		f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/mussmann/local/src/primer3/primer3/src/primer3_config") #add path to thermodynamic
		f.write("\n")
		f.write("=")
		f.write("\n")

		f.close()

	def parseFile(self):
		f=open(self.fas)
		data=f.readlines()
		f.close()
		data=[l.strip() for l in data]

		seqinfo=dict()
		seqinfo["ID"] = data.pop(0).replace(">", "")
		seqinfo["SEQUENCE"] = "".join(data).upper()

		#print(seqinfo["ID"])
		#print(seqinfo["SEQUENCE"])

		return seqinfo

	def makeCommand(self, fn):
		outfn = fn + ".p3out.txt"
		infn = fn + ".p3in.txt"
		out=os.path.join(self.p3outDir, outfn)
		indir=os.path.join(self.p3inDir, infn)
		string = "primer3_core -format_output -output=" + out + " < " + indir
		return string

