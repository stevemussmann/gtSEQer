from __future__ import print_function

import argparse
import os.path
import distutils.util

class ComLine():
	'Class for implementing command line options'


	def __init__(self, args):
		parser = argparse.ArgumentParser()
		parser.add_argument("-v", "--vcf",
							dest='vcf',
							required=True,
							help="Specify a vcf file for input that contains target SNPs"
		)
		parser.add_argument("-g", "--genome",
							dest='genome',
							required=True,
							help="Specify the name of a reference genome in fasta format"
		)
		parser.add_argument("-z", "--vcfgz",
							dest='vcfgz',
							nargs='?',
							help="Specify the name of the gzipped vcf (optional)"
		)
		parser.add_argument("-i", "--indlist",
							dest='indlist',
							required=True,
							help="Specify the list of samples to evaluate"
		)
		parser.add_argument("-d", "--directory",
							dest='directory',
							required=True,
							help="Specify the directory containing coverage files."
		)
		parser.add_argument("-m", "--mask",
							dest='mask',
							default="all.genomecov.txt.gz",
							help="Specify the file extension for coverage file names. The beginning of the name should match individuals in the indlist."
		)
		parser.add_argument("-c", "--cover",
							dest='cover',
							type=int,
							default=5,
							help="Specify the minimum sequencing coverage."
		)
		parser.add_argument("-f", "--flank",
							dest='flank',
							type=int,
							default=150,
							help="Specify the amount of flanking sequence you want to retain on either side of the SNP. "
		)
		parser.add_argument("-x", "--overwrite",
							dest='overwrite',
							action='store_true',
							help="Turn on to overwrite previous files"

		)

		self.args = parser.parse_args()

		#sets default value of vcfgz to append "gz" extension to vcf arg
		if self.args.vcfgz is None:
			self.args.vcfgz = self.args.vcf + ".gz"

		#check if files and directories exist
		self.exists( self.args.vcf )
		self.exists( self.args.indlist )
		self.exists( self.args.genome )
		self.dirExists( self.args.directory )


	def exists(self, filename):
		if( os.path.isfile(filename) != True ):
			print(filename, "does not exist")
			print("Exiting program...")
			print("")
			raise SystemExit
	
	def dirExists(self, dirname):
		if( os.path.isdir(dirname) != True ):
			print( dirname, "either does not exist, or is not a directory.")
			print("Exiting program...")
			print("")
			raise SystemExit
