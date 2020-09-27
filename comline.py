from __future__ import print_function

import argparse
import os.path
import distutils.util

class ComLine():
	'Class for implementing command line options'


	def __init__(self, args):
		parser = argparse.ArgumentParser()
		parser.add_argument("-c", "--cover",
							dest='cover',
							type=int,
							default=5,
							help="Specify the minimum sequencing coverage."
		)
		parser.add_argument("-d", "--directory",
							dest='directory',
							required=True,
							help="Specify the directory containing coverage files."
		)
		parser.add_argument("-f", "--flank",
							dest='flank',
							type=int,
							default=150,
							help="Specify the amount of flanking sequence you want to retain on either side of the SNP when extracting loci from the reference genome."
		)
		parser.add_argument("-F", "--probeflank",
							dest='probeflank',
							type=int,
							default=8,
							help="Specify the amount of flanking sequence you want to retain on either side of the SNP whene extracting probes from the reference genome."
		)
		parser.add_argument("-g", "--genome",
							dest='genome',
							required=True,
							help="Specify the name of a reference genome in fasta format"
		)
		parser.add_argument("-i", "--indlist",
							dest='indlist',
							required=True,
							help="Specify the list of samples to evaluate"
		)
		parser.add_argument("-l", "--loci",
							dest='loci',
							required=True,
							help="Specify a file for input that contains target SNP coordinates in the first two columns (Chromosome \\tab Position). This can also be a filtered VCF file."
		)
		parser.add_argument("-m", "--mask",
							dest='mask',
							default="all.genomecov.txt.gz",
							help="Specify the file extension for coverage file names. The beginning of the name should match individuals in the indlist."
		)
		#implement
		parser.add_argument("-p", "--product",
							dest='product',
							type=int,
							default=105,
							help="Specify the optimal PCR product size."
		)
		#implement
		parser.add_argument("-P", "--producterr",
							dest='producterr',
							type=int,
							default=15,
							help="Specify the +/- error for optimal PCR product size."
		)
		#implement
		parser.add_argument("-s", "--size",
							dest='size',
							type=int,
							default=21,
							help="Specify the ideal length for a primer"
		)
		#implement
		parser.add_argument("-S", "--sizeerr",
							dest='sizeerr',
							type=int,
							default=3,
							help="Specify the acceptable +/- error for ideal primer length"
		)
		#implement
		parser.add_argument("-t", "--temp",
							dest='temp',
							type=int,
							default=61,
							help="Specify targeted melting temperature for primers."
		)
		#implement
		parser.add_argument("-T", "--temperr",
							dest='temperr',
							type=int,
							default=3,
							help="Specify the acceptable +/- error for primer melting temperature."
		)
		parser.add_argument("-v", "--vcf",
							dest='vcf',
							required=True,
							help="Specify the name of the gzipped vcf"
		)
		#implement
		parser.add_argument("-x", "--overwrite",
							dest='overwrite',
							action='store_true',
							help="Turn on to overwrite previous files (Not fully implemented yet)."

		)

		self.args = parser.parse_args()

		#check if files and directories exist
		self.exists( self.args.vcf )
		self.exists( self.args.loci )
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
