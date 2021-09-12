# script store all alterations not observed in WT or Plasmid library samples
# updated to not count spurious SNPs
# updated to not count spurious indels too
# updated to not count SNPs at all (only look at insertions and deletions, even if they happen once)
# edited: March 2018

from __future__ import division

import sys
import argparse
import os
import operator
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

def aggregate_files(input_file_list,modality,low_coverage_file,mutation_summary_file):
	if modality=='BE4' or modality=='ABE7.10':
		mutation_summary_file.write('File\tGene\tTotalCount\tPos1\tPos2\tPos3\tPos4\tPos5\tPos6\tPos7\tPos8\tPos9\tPos10\tPos11\tPos12\tPos13\tPos14\tPos15\tPos16\tPos17\tPos18\tPos19\tPos20\n')
	else:
		mutation_summary_file.write('Sample\tDisplay_Name\tTarget_Site\tNHEJ_Mutation_Count\tOOF_Mutation_Count\tTotalCount\tTotal_NHEJ_Mutation_Percentage\tOOF_Mutation_Percentage\tAlterations\n')
	# go through mutation file list
	for infile in input_file_list:
		ifile = open(infile,'r')
		lineCount = 0
		for line in ifile:
			if lineCount==1:
				# split the line
				line = line.rstrip('\r\n')
				parts = line.split('\t')
				if modality=='BE4' or modality=='ABE7.10':
					depth = parts[3]
				else:
					depth = parts[5]
				if int(depth) > 0:
					mutation_summary_file.write(line + '\n')
					if int(depth) < 100:
						head, tail = os.path.split(parts[0])
						# remove the mpileup_tab
						tail = tail.replace('_mpileup.tab','')
						tail = tail.replace('_bwamem_sorted.bam','')
						low_coverage_file.write(tail + '\n')

			lineCount += 1
		ifile.close()

	# close file handles
	mutation_summary_file.close()
	low_coverage_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file_list',nargs='+',required=True)
	parser.add_argument('-o','--mutation_summary_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-l','--low_coverage_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-m','--modality',required=True)
	opts = parser.parse_args(argv)
	aggregate_files (opts.input_file_list,opts.modality,opts.low_coverage_file,opts.mutation_summary_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])