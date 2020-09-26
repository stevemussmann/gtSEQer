#!/usr/bin/env python

from comline import ComLine
from coverage import Coverage
from DefaultListOrderedDict import DefaultListOrderedDict
from extractor import Extractor
from fasta import Fasta
from genome import Genome
from maketestdir import MakeTestDir
from muscle import Muscle
from primer3 import Primer3
from vcf import VCF

import os.path
import sys
from os import listdir


def main():
	input = ComLine(sys.argv[1:])

	gen = Genome(input.args.genome)

	# identify target SNPs and return coordinates in reference genome
	target = VCF(input.args.vcf)
	snps = target.getSNPs()

	# apply coverage cutoff to BEDtools output
	cov = Coverage(input.args.directory, input.args.indlist, input.args.mask)
	cov.applyMinCov(input.args.cover, input.args.overwrite)

	# Extract target plus flanking region from reference genome
	regions = Extractor(snps, input.args.indlist, input.args.flank, gen.seqLengths)
	regions.extract()
	regions.makeCommands(input.args.genome, input.args.vcfgz) #returns list of output files

	## operate on all output Fasta files
	# filter extracted_regions
	mtdex = MakeTestDir("extracted_regions")
	ex = mtdex.testDir()
	fastaFiles=[f for f in listdir(ex) if os.path.isfile(os.path.join(ex, f))]
	for f in fastaFiles:
		fpath=os.path.join(ex, f)
		fas = Fasta(fpath, "filtered_sequences")

	# filter probes
	mtdProbe = MakeTestDir("extracted_regions_probes")
	probe = mtdProbe.testDir()
	probeFiles = [f for f in listdir(probe) if os.path.isfile(os.path.join(probe, f))]
	for f in probeFiles:
		fpath = os.path.join(probe, f)
		fas = Fasta(fpath, "filtered_probes")

	# use Muscle to align files
	mtdal = MakeTestDir("filtered_sequences")
	al = mtdal.testDir()
	filteredFiles=[f for f in listdir(al) if os.path.isfile(os.path.join(al, f))]
	for f in filteredFiles:
		fpath=os.path.join(al, f)
		filt = Muscle(fpath, "muscle_aligned")

	mtdalp = MakeTestDir("filtered_probes")
	alp = mtdalp.testDir()
	filteredProbes=[f for f in listdir(alp) if os.path.isfile(os.path.join(alp, f))]
	for f in filteredProbes:
		fpath=os.path.join(alp, f)
		filt = Muscle(fpath, "muscle_aligned_probes")

	# Run primer3
	mtdcf = MakeTestDir("muscle_aligned_consensus")
	cf = mtdcf.testDir()
	conFasta=[f for f in listdir(cf) if os.path.isfile(os.path.join(cf, f))]
	for f in conFasta:
		p3 = Primer3(f)

main()

raise SystemExit
