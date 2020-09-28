# gtSEQer: Primer design tool for GT-seq
A pipeline to identify GT-seq loci from a reference genome, VCF file, and BEDtools genomecov files. 

## Citing gtSEQer
A manuscript will be developed for publication

## Installation & Setup for gtSEQer:

This pipeline was written to be run on Unix based operating systems, such as the various Linux distributions and Mac OS X.  To get started, clone this project to the desired location on your computer.  

### Python Dependencies
* BioPython
* pandas

### Program Dependencies
* samtools
* bcftools
* primer3
* muscle

## Required inputs:

The following options are required:
* (**-l, --loci**): list of target SNPs. This can be a two-colum tab-delimited file with the format of CHROMOSOME \<tab\> POSITION, or you can provide an uncompressed VCF file.
* (**-g, --genome**): reference genome in fasta format. File must be uncompressed.
* (**-v, --vcf**): compressed VCF file containing all detected variants. Before executing the pipeline, this file must be indexed using the command "bcftools index all_variants.vcf"
* (**-i, --indlist**): list of individuals present the VCF input that will be used for GTseq locus development.
* (**-d, --directory**): path to bedtools coverage files for each individual in your list of individuals. These can be generated using the "bedtools genomecov -bga" command.
