#!/usr/bin/env python

from comline import ComLine
from coverage import Coverage
from DefaultListOrderedDict import DefaultListOrderedDict
from extractor import Extractor
from fasta import Fasta
from muscle import Muscle
from primer3 import Primer3
from vcf import VCF

import os.path
import sys
from os import listdir


def main():
	cwd = os.getcwd()

	input = ComLine(sys.argv[1:])

	# identify target SNPs and return coordinates in reference genome
	target = VCF(input.args.vcf)
	snps = target.getSNPs()

	# apply coverage cutoff to BEDtools output
	cov = Coverage(input.args.directory, input.args.indlist, input.args.mask)
	cov.applyMinCov(input.args.cover, input.args.overwrite)

	# Extract target plus flanking region from reference genome
	#regions = Extractor(snps, input.args.indlist, input.args.flank)
	#regions.extract()
	#regions.makeCommands(input.args.genome, input.args.vcfgz) #returns list of output files

	## operate on all output Fasta files
	# filter files
	ex=os.path.join(cwd,"extracted_regions")
	print(ex)
	fastaFiles=[f for f in listdir(ex) if os.path.isfile(os.path.join(ex, f))]
	for f in fastaFiles:
		fpath=os.path.join(ex, f)
		fas = Fasta(fpath)

	# use Muscle to align files
	al=os.path.join(cwd, "filtered_sequences")
	filteredFiles=[f for f in listdir(ex) if os.path.isfile(os.path.join(al, f))]
	for f in filteredFiles:
		fpath=os.path.join(al, f)
		filt = Muscle(fpath)

	# Run primer3
	#for f in files:
	#	p3 = Primer3(f)

main()

raise SystemExit
