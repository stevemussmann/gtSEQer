from __future__ import print_function

#import json
import os.path
import subprocess
import sys

from checkpoint import Checkpoint
from maketestdir import MakeTestDir
from program import Program

class Extractor():
	'Class for pulling genome regions from reference'

	def __init__(self, f, i, e, l, p):
		# tracks which files have been done
		self.tracker = Checkpoint("extractor.json")
		self.done=self.tracker.loadJson()

		self.targets = f
		self.flank=e
		self.probeFlank=p
		self.genLen=l
		f=open(i)
		d=f.readlines()
		f.close()
		d=[l.rstrip() for l in d]
		self.indlist = d
		self.coords=list()

		self.probeTargets = dict()
		
		#make sure proper directories exist
		cov_ex = MakeTestDir("coverage_files")
		self.covdir=cov_ex.testDir()
		out_ex = MakeTestDir("extracted_regions")
		self.outdir=out_ex.testDir()
		probe_ex = MakeTestDir("extracted_regions_probes")
		self.probedir=probe_ex.testDir()

	def extract(self):
		for chrom, snps in self.targets.items():
			for snp in snps:
				c=list()
				p=list()

				start = int(snp) - self.flank
				if( start < 1 ):
					start=1
				end = int(snp) + self.flank
				if(end > self.genLen[chrom]):
					end = self.genLen[chrom]

				pstart = int(snp) - self.probeFlank
				if( pstart < 1 ):
					pstart = 1
				pend = int(snp) + self.probeFlank
				if(pend > self.genLen[chrom]):
					pend = self.genLen[chrom]

				#make list with coordinates for samtools
				c.append(chrom)
				c.append(":")
				c.append(str(start))
				c.append("-")
				c.append(str(end))

				#concatenate coordinates
				newstr="".join(c)
				self.coords.append(newstr)

				p.append(chrom)
				p.append(":")
				p.append(str(pstart))
				p.append("-")
				p.append(str(pend))

				#concatenate probe coordinates
				pnewstr="".join(p)

				self.probeTargets[newstr] = pnewstr

	def readFile(self):
		f=open(self.vcf)
		data=f.readlines()
		f.close()
		data = [l.rstrip() for l in data]
		return data

	def makeCommands(self, genome, gz):
		p=False
		self.runCommands(genome, gz)

	def runCommands(self, genome, gz):
		for coord in self.coords:
			if coord not in self.done.keys():
				self.done[coord]=list()
			for ind in self.indlist:
				if(ind not in self.done[coord]):
					outfile=coord.replace(":","_")
					outfn=os.path.join(self.outdir,outfile)
					probefn=os.path.join(self.probedir,outfile)
					cf=ind + ".cov.txt"
					cfn=os.path.join(self.covdir, cf)
					print("Extracting", coord, "for", ind)

					# Run command to extract sequence
					command = "samtools faidx " + genome + " " + coord + " | bcftools consensus -s " + ind + " -m " + cfn +  " -I " + gz
					prog = Program(command)
					output = prog.runProgram()
					
					# write sequence to file
					self.parseOutput(output, outfn, ind)

					# Run command to extract probe
					probeCommand = "samtools faidx " + genome + " " + self.probeTargets[coord] + " | bcftools consensus -s " + ind + " -m " + cfn + " -I " + gz
					probeProg = Program(probeCommand)
					probeOutput = probeProg.runProgram()

					# write probe to file
					self.parseOutput(probeOutput, probefn, ind)

					# write to json file
					self.done[coord].append(ind)
					self.tracker.writeJson(self.done)

	def parseOutput(self, o, fn, ind):
		l = o.splitlines()
		l[0] = l[0] + "_" + ind
		with open(fn, 'a') as fh:
			for item in l:
				fh.write(item)
				fh.write("\n")
	
