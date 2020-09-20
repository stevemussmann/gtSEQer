from __future__ import print_function

class Primer3():
	'Class for executing primer3'

	def __init__(self,fn):
		self.fn=fn
		self.procFile()

	def makeFile(self):
		fn = self.locus + ".p3.txt"
		f = open(fn, 'w')
		f.write("SEQUENCE_ID=") #finish this line
		f.write("SEQUENCE_TEMPLATE=") #finish
		f.write("SEQUENCE_TARGET=") #add "start coord, length"
		f.write("PRIMER_TASK=pick_detection_primers")
		f.write("PRIMER_PICK_LEFT_PRIMER=1")
		f.write("PRIMER_PICK_INTERNAL_OLIGO=0")
		f.write("PRIMER_PICK_RIGHT_PRIMER=1")
		f.write("PRIMER_NUM_RETURN=1")
		f.write("PRIMER_OPT_SIZE=22")
		f.write("PRIMER_MIN_SIZE=18")
		f.write("PRIMER_MAX_SIZE=24")
		f.write("PRIMER_MAX_NS_ACCEPTED=0")
		f.write("PRIMER_PRODUCT_SIZE_RANGE=100-150")
		f.write("PRIMER_GC_CLAMP=2")
		f.write("PRIMER_MIN_TM=57")
		f.write("PRIMER_OPT_TM=60")
		f.write("PRIMER_MAX_TM=65")
		f.write("PRIMER_PAIR_MAX_DIFF_TM=2")
		f.write("PRIMER_MAX_POLY_X=2")
		f.write("PRIMER_MAX_HAIRPIN_TH=20")
		f.write("PRIMER_MAX_SELF_ANY_TH=20")
		f.close("P3_FILE_FLAG=1")
		f.write("SEQUENCE_INTERNAL_EXCLUDED_REGION=") #add "start coord, length"
		f.write("PRIMER_OPT_GC_PERCENT=50")
		f.write("PRIMER_MIN_GC_PERCENT=40")
		f.write("PRIMER_MAX_GC_PERCENT=60")
		f.write("PRIMER_EXPLAIN_FLAG=1")
		f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/user/local/src/primer3/primer3-master/src/primer3_config") #add path to thermodynamic
		f.write("=")

		f.close()

	def procFile(self):
		f=open(self.fn)
		data=f.readlines()
		f.close()
		data=[l.strip() for l in data]

		seqinfo=dict()
		seqinfo["ID"] = data.pop(0).replace(">", "")
		seqinfo["SEQUENCE"] = "".join(data).upper()

		print(seqinfo["ID"])
		print(seqinfo["SEQUENCE"])
		print(len(seqinfo["SEQUENCE"]))

