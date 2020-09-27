from __future__ import print_function

import os.path
import subprocess
import sys

from program import Program
from maketestdir import MakeTestDir

class Primer3():
	'Class for executing primer3'

	def __init__(self, fn, flank, probeflank, product, productErr, size, sizeErr, temp, tempErr):
		# make directories
		mtd_mac = MakeTestDir("muscle_aligned_consensus")
		mac = mtd_mac.testDir()
		mtd_mapc = MakeTestDir("muscle_aligned_probes_consensus")
		mapc = mtd_mapc.testDir()

		#make names with paths for template and probe sequences
		template=os.path.join(mac, fn)
		probe=os.path.join(mapc, fn)

		# make directories to hold primer3 input and output
		mtd_p3in = MakeTestDir("p3in")
		mtd_p3out = MakeTestDir("p3out")
		self.p3inDir=mtd_p3in.testDir()
		self.p3outDir=mtd_p3out.testDir()
		
		#calculate probe length
		probeLength = (2*probeflank)+1

		# parse sequence files
		seqinfo = self.parseFile(template)
		probeinfo = self.parseFile(probe)

		# find starting position of probe
		index = seqinfo["SEQUENCE"].find(probeinfo["SEQUENCE"])
		print(str(index))

		if(index < 0):
			print("WARNING: Probe not found in sequence", fn)
			print("Attempting to place probe in expected location.")
			print("Manually inspect any locus present in warnlist.txt")
			index=flank-probeflank
			warnfh=open("warnlist.txt", 'a')
			warnfh.write(fn)
			warnfh.write("\n")
			warnfh.close()
		
		# make p3in file
		self.makeFile(seqinfo["ID"], seqinfo["SEQUENCE"], index, probeLength, product, productErr, size, sizeErr, temp, tempErr)

		# run primer3_core
		comm, out = self.makeCommand(fn)
		#print(comm)

		print("Finding primers for", fn)

		prog = Program(comm)
		prog.runProgram()

		self.parseOutput(out, fn)

	def makeFile(self, name, seq, index, probeLength, product, productErr, size, sizeErr, temp, tempErr):
		fn=name + ".p3in.txt"
		fn=os.path.join(self.p3inDir, fn)
	
		# make sure any ambiguities are converted to Ns
		seq=self.convertN(seq)

		f = open(fn, 'w')
		f.write("SEQUENCE_ID=")
		f.write(name) #sequence name
		f.write("\n")
		f.write("SEQUENCE_TEMPLATE=")
		f.write(seq) #sequence
		f.write("\n")
		f.write("SEQUENCE_TARGET=")
		f.write(str(index)) # start coordinate
		f.write(",") # length
		f.write(str(probeLength))
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
		f.write("PRIMER_OPT_SIZE=")
		f.write(str(size))
		f.write("\n")
		f.write("PRIMER_MIN_SIZE=")
		f.write(str(size-sizeErr))
		f.write("\n")
		f.write("PRIMER_MAX_SIZE=")
		f.write(str(size+sizeErr))
		f.write("\n")
		f.write("PRIMER_MAX_NS_ACCEPTED=0")
		f.write("\n")
		f.write("PRIMER_PRODUCT_SIZE_RANGE=")
		f.write(str(product-productErr))
		f.write("-")
		f.write(str(product+productErr))
		f.write("\n")
		f.write("PRIMER_PRODUCT_OPT_SIZE=")
		f.write(str(product))
		f.write("\n")
		f.write("PRIMER_GC_CLAMP=1")
		f.write("\n")
		f.write("PRIMER_MIN_TM=")
		f.write(str(temp-tempErr))
		f.write("\n")
		f.write("PRIMER_OPT_TM=")
		f.write(str(temp))
		f.write("\n")
		f.write("PRIMER_MAX_TM=")
		f.write(str(temp+tempErr))
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
		f.write("SEQUENCE_INTERNAL_EXCLUDED_REGION=")
		f.write(str(index))
		f.write(",") #add "start coord, length"
		f.write(str(probeLength))
		f.write("\n")
		f.write("PRIMER_OPT_GC_PERCENT=50")
		f.write("\n")
		f.write("PRIMER_MIN_GC_PERCENT=40")
		f.write("\n")
		f.write("PRIMER_MAX_GC_PERCENT=60")
		f.write("\n")
		f.write("PRIMER_EXPLAIN_FLAG=1")
		f.write("\n")
		#f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/user/local/src/primer3/primer3-master/src/primer3_config") #add path to thermodynamic
		f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/mussmann/local/src/primer3/primer3/src/primer3_config") #add path to thermodynamic
		f.write("\n")
		f.write("=")
		f.write("\n")

		f.close()

	def convertN(self, seq):
		iupac = ["M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "-"]
		for b in iupac:
			seq=seq.replace(b, "N")
		return seq


	def parseFile(self, fasta):
		f=open(fasta)
		data=f.readlines()
		f.close()
		data=[l.strip() for l in data]

		seqinfo=dict()
		seqinfo["ID"] = data.pop(0).replace(">", "")
		seqinfo["SEQUENCE"] = "".join(data).upper()

		return seqinfo

	def makeCommand(self, fn):
		outfn = fn + ".p3out.txt"
		infn = fn + ".p3in.txt"
		out=os.path.join(self.p3outDir, outfn)
		indir=os.path.join(self.p3inDir, infn)
		string = "primer3_core -format_output -output=" + out + " < " + indir
		return string, out

	def parseOutput(self, f, fn):
		fh=open(f)
		data=fh.readlines()
		fh.close()
		data = [l.strip() for l in data]

		outf = "summary.txt"
		outfh=open(outf, 'a')
		for line in data:
			if(line.startswith("OLIGO") or line.startswith("LEFT") or line.startswith("RIGHT")):
				if(line.startswith("OLIGO")):
					outfh.write("> ")
					outfh.write(fn)
					outfh.write("\n")
				outfh.write(line)
				outfh.write("\n")
		outfh.close()

