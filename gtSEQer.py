#!/usr/bin/env python

from checkpoint import Checkpoint
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

def runMuscle(dirName, jsonName, muscleDirName):
	mtdal = MakeTestDir(dirName)
	al = mtdal.testDir()
	filteredFiles=[f for f in listdir(al) if os.path.isfile(os.path.join(al, f))]
	muscleTracker=Checkpoint(jsonName)
	muscleDone=muscleTracker.loadJson()
	for f in filteredFiles:
		if f not in muscleDone.keys():
			fpath=os.path.join(al, f)
			filt = Muscle(fpath, muscleDirName)
			muscleDone[f]=1
			muscleTracker.writeJson(muscleDone)


def main():
	input = ComLine(sys.argv[1:])

	gen = Genome(input.args.genome)

	# identify target SNPs and return coordinates in reference genome
	target = VCF(input.args.loci)
	snps = target.getSNPs()

	# apply coverage cutoff to BEDtools output
	cov = Coverage(input.args.directory, input.args.indlist, input.args.mask)
	cov.applyMinCov(input.args.cover, input.args.overwrite)

	# Extract target plus flanking region from reference genome
	regions = Extractor(snps, input.args.indlist, input.args.flank, gen.seqLengths, input.args.probeflank)
	regions.extract()
	regions.makeCommands(input.args.genome, input.args.vcf) #returns list of output files

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

	# use Muscle to align sequences
	runMuscle("filtered_sequences", "muscle_seqs.json", "muscle_aligned")

	# use Muscle to align probes
	runMuscle("filtered_probes", "muscle_probes.json", "muscle_aligned_probes")

	# Run primer3
	if os.path.isfile("summary.txt"):
		os.remove("summary.txt")
	if os.path.isfile("warnlist.txt"):
		os.remove("warnlist.txt")
	mtdcf = MakeTestDir("muscle_aligned_consensus")
	cf = mtdcf.testDir()
	conFasta=[f for f in listdir(cf) if os.path.isfile(os.path.join(cf, f))]
	for f in conFasta:
		p3 = Primer3(f, input.args.flank, input.args.probeflank)

main()

raise SystemExit
