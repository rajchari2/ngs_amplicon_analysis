# script to filter SAM alignments 

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

def filter_sam_file(input_file,output_file):
	for line in input_file:
		if line.startswith('@')==False:
			line = line.rstrip('\r\n')
			parts = line.split('\t')
			if 'H' not in parts[5] and 'S' not in parts[5]:
				output_file.write(line + '\n')
		else:
			output_file.write(line)
	input_file.close()
	output_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	filter_sam_file(opts.input_file,opts.output_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])