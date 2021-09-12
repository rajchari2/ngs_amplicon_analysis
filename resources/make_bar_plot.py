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

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

def make_indel_rate_plot(input_file,output_file):
	# read in the data
	x_data_labels = []
	y_data_nhej = []
	y_data_oof = []
	line_count = 0
	mapped_nhej_data = defaultdict(float)
	mapped_oof_data = defaultdict(float)
	for line in input_file:
		if line_count > 0:
			line = line.rstrip('\r\n')
			parts = line.split('\t')
			s_name = parts[1] + ' (N=' + parts[5] + ')'
			mapped_nhej_data[s_name] = float(parts[6])
			mapped_oof_data[s_name] = float(parts[7])
		line_count += 1
	input_file.close()

	# put the data into arrays
	for key in sorted(mapped_nhej_data.keys()):
		x_data_labels.append(key)
		y_data_nhej.append(mapped_nhej_data[key])
		y_data_oof.append(mapped_oof_data[key])

	# number of bars
	num_bars = len(x_data_labels)
	x_locs = np.arange(num_bars)
	w = 0.2

	# add axis labels
	fig, ax = plt.subplots()
	nhej = ax.bar(x_locs,y_data_nhej,width=0.2,color='b')
	oof = ax.bar(x_locs+w,y_data_oof,width=0.2,color='c')
	plt.xlabel('Guide RNA Candidate')
	plt.ylabel('NHEJ mutation rate (%)')
	ax.set_xticks(x_locs+w/2)
	ax.set_xticklabels(x_data_labels,rotation=90)

	# plot legend
	plt.legend([nhej, oof],['Total NHEJ mutation rate (%)', 'Out of frame (OOF) mutation rate (%)'])

	# save to file
	plt.savefig(output_file)

def make_base_editing_plot(input_file,modality,output_file):
	# only going to use position 5 to 12
	# output_file.write('File\tGene\tTotalCount\tPos1\tPos2\tPos3\tPos4\tPos5\tPos6\tPos7\tPos8\tPos9\tPos10\tPos11\tPos12\tPos13\tPos14\tPos15\tPos16\tPos17\tPos18\tPos19\tPos20\n')
	# data structure to hold data
	data_table = defaultdict(list)
	legend_list = []
	# go through each summary
	line_count = 0
	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		if line_count > 0:
			depth = int(parts[2])
			index_to_start = 3
			while index_to_start < len(parts):
				value = int(parts[index_to_start]) / depth
				data_table[parts[0]].append(value)
				index_to_start += 1
		else:
			x_data_labels = parts[3:]
		line_count += 1

	# bar locations
	ind = np.arange(0,3*len(x_data_labels),3)
	width = 0.3
	
	# add axis labels
	fig, ax = plt.subplots(figsize=(20,10))
	ax.set_xticks(ind)
	ax.set_xticklabels(x_data_labels,rotation=90)

	# iterate through data
	iteration = 1
	for key in sorted(data_table):
		bar_set = ax.bar(ind + (iteration * width),data_table[key],width,align="center")
		#print(data_table[key])
		iteration += 1
		# only use the sample name
		head, tail = os.path.split(key)
		# remove the mpileup_tab
		tail = tail.replace('_mpileup.tab','')
		legend_list.append(tail)

	plt.xlabel('Position in spacer')
	if modality=='BE4':
		plt.ylabel('C>T conversion frequency')
	else:
		plt.ylabel('A>G conversion frequency')

	# set the axis tick labels
	ax.legend(legend_list)

	# save to file
	plt.savefig(output_file)

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--mutation_summary_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--mutation_graph_output',required=True)
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		make_base_editing_plot(opts.input_file,opts.modality,opts.output_file)
	else:
		make_indel_rate_plot(opts.mutation_summary_file,opts.mutation_graph_output)
 
if __name__ == '__main__':
	main(sys.argv[1:])
