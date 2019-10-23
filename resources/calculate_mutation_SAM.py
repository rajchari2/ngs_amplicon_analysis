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
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from collections import defaultdict
from collections import Counter
from operator import itemgetter

def process_cigar(cigar_tuple,md_tag,aligned_pairs,query_qualities):
	alterations = []
	pos = 1
	# if the cigar tuple is only 1 entry check if it has a SNP
	if len(cigar_tuple)==1:
		for seq_pos in aligned_pairs:
			pos = int(seq_pos[0]) + 1
			base = seq_pos[2]
			if base.islower() and int(query_qualities[pos-1]) >= 20:
				alt = str(pos) + ':' + str(pos) + ':SNP'
				alterations.append(alt)
	else:
		for entry in cigar_tuple:
			# match
			if int(entry[0])==0:
				pos = pos + int(entry[1])
			elif int(entry[0])==1:
				alt = str(pos) + ':' + str(pos) + ':Ins'
				alterations.append(alt)
			elif int(entry[0])==2:
				end_pos = pos + int(entry[1])
				alt = str(pos) + ':' + str(end_pos) + ':Del'
				pos = end_pos
				alterations.append(alt)
	return alterations

def find_target_indices(target_site,reference_file):
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	gene = ''
	for record in SeqIO.parse(reference_file,'fasta'):
		refDB[str(record.id)] = str(record.seq)
		targetDB[str(record.id)] = target_site
		gene = str(record.id)

	# identify the target site start and end
	targetSeq = targetDB[gene].upper()
	refSeq = refDB[gene].upper()

	# identify where in the sequence it is
	tSeq = Seq(targetSeq)
	targetSeqRC = str(tSeq.reverse_complement())

	# find the sequence
	index = refSeq.find(targetSeq)

	if index==-1:
		rcIndex = refSeq.find(targetSeqRC)
		startIndex = rcIndex
		endIndex = startIndex + len(targetSeqRC)
	else:
		startIndex = index	
		endIndex = startIndex + len(targetSeq)

	return startIndex,endIndex,gene


def calculate_NHEJ_mutation_rate (sample_sam_file,control,reference_file,target_site,output_file,diversity_file):

	# variables to store information
	sample_read_total = 0
	sample_mutation_count = 0
	control_cigar_strings = []

	# control bam file
	control_sam_file = 'processed/' + control + '_bwamem_sorted.bam'

	# get the start and end indexes of the target in the reference
	target_start,target_end,gene = find_target_indices(target_site,reference_file)

	# write header in output file
	output_file.write('Sample\tAmplicon_Name\tTarget_Site\tMutated_Read_Count\tTotal_Read_Count\tNHEJ_Mutation_Rate\n')
	diversity_file.write('Sample\tAmplicon_Name\tTarget_Site\tMutated_Read_Count\tNum_Distinct_Deletions\tNum_Distinct_Insertions\tTotal_Distinct_Events\n')

	# add data structures for insertions and deletions
	distinct_insertions = []
	distinct_deletions = []

	# go through control file to determine "false" mutations
	ctrl_sam = pysam.AlignmentFile(control_sam_file,'rb')
	for read in ctrl_sam.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==0:
			# check if the read has a mutation that involves the target sequence
			md_tag = read.get_tag('MD')
			ctrl_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities)
			unique_entry = md_tag + '_' + read.cigarstring
			# if alteration is non zero
			valid_alteration = False
			if len(ctrl_alterations) > 0:
				for alt in ctrl_alterations:
					[start,end,alt_type] = alt.split(':')
					if ((int(start) >= target_start and int(start) <= target_end) or (int(end) >= target_start and int(end) <= target_end)):
						valid_alteration = True
			# check if valid alteration
			if valid_alteration==True and unique_entry not in control_cigar_strings:
				# add to the control
				control_cigar_strings.append(unique_entry)

	# go through sample BAM file
	samfile = pysam.AlignmentFile(sample_sam_file,"rb")
	for read in samfile.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==0:
			sample_read_total += 1
			md_tag = read.get_tag('MD')
			sample_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities)
			unique_entry = md_tag + '_' + read.cigarstring
			# only count if the tuple is not in the control
			if len(sample_alterations) > 0 and unique_entry not in control_cigar_strings:
				# if alteration is non zero
				valid_alteration = False
				for alt in sample_alterations:
					[start,end,alt_type] = alt.split(':')
					if ((int(start) >= target_start and int(start) <= target_end) or (int(end) >= target_start and int(end) <= target_end)):
						valid_alteration = True
						# add the valid alteration to the tally
						if alt_type=='Ins':
							if alt not in distinct_insertions:
								distinct_insertions.append(alt)
						elif alt_type=='Del':
							if alt not in distinct_deletions:
								distinct_deletions.append(alt)
				if valid_alteration==True:
					sample_mutation_count += 1
	if sample_read_total >= 100:
		rate = (sample_mutation_count / sample_read_total) * 100
	else:
		rate = 'N/A'

	# update sample name
	sample_name = sample_sam_file.name
	sample_name = sample_name.replace('processed/','')
	sample_name = sample_name.replace('_bwamem_sorted.bam','')

	total_distinct = len(distinct_deletions) + len(distinct_insertions)
	output_file.write(sample_name + '\t' + gene + '\t' + target_site + '\t' + str(sample_mutation_count) + '\t' + str(sample_read_total) + '\t' + str(rate) + '\n')
	diversity_file.write(sample_name + '\t' + gene + '\t' + target_site + '\t' + str(sample_mutation_count) + '\t' + str(len(distinct_deletions)) + '\t' + str(len(distinct_insertions)) + '\t' + str(total_distinct) + '\n')

	# close file handles
	samfile.close()
	ctrl_sam.close()
	output_file.close()
	diversity_file.close()

#def calculate_base_editing_rate (sample_sam_file,control_sam_file,reference_file,target_site,output_file):
	


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s','--sample_sam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_sam_file',required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-d','--diversity_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		calculate_base_editing_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.output_file,opts.diversity_file)
	else:
		calculate_NHEJ_mutation_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.output_file,opts.diversity_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])

