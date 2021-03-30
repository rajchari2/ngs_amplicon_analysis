# ngs_amplicon_analysis
Scripts and files to analyze amplicon sequencing data

# input files

Description ("tab") file - a tab-delimited text file which describes parameters of each set of paired end fastq files

Reference fasta file - can either be a genome fasta file or just a fasta file with the amplicon seqeunce (which also is indexed by bwa)

# values in the description file

name of fastq file <tab> reference file <tab> coordinates <tab> amplicon_name <tab> target site <tab> control sample name

if reference file is just the amplicon, use "full" for coordinates
