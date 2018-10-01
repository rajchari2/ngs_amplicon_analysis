# script to calculation mutation rate from SAM file
# Use the FLAG tag, MD tag and CIGAR string for calling
# FLAG tag: Column 2 (make sure this is 0)
# Map start: Column 4 (Should be position 1)
# CIGAR: Column 6
# MD tag: Column 13

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

def calculate_NHEJ_mutation_rate (sample_sam_file,control_sam_file,reference_file,target_site,output_file):
	# variables to store information
	sample_read_total = 0
	sample_mutation_count = 0
	control_mutation_dictionary = defaultdict(str)

	# go through control and get all of the reads in the control sample that would be called mutations
	for line in control_sam_file:
		if line.startswith('@')==False:
			# parse the line
			line = line.rstrip('\r\n')
			parts = line.split('\t')
			# only count the read if FLAG=0 and Map Start = 1
			if parts[1]=='0' and parts[3]=='1':
				


def calculate_base_editing_rate (sample_sam_file,control_sam_file,reference_file,target_site,output_file):
	


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s','--sample_sam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_sam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		calculate_base_editing_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.output_file)
	else:
		calculate_NHEJ_mutation_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.output_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])

