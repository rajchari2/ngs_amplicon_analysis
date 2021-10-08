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

def process_cigar(cigar_tuple,md_tag,aligned_pairs,query_qualities,query_alignment_sequence,query_sequence):
	alterations = []
	pos = 1
	# if the cigar tuple is only 1 entry check if it has a SNP
	if len(cigar_tuple)==1:
		for seq_pos in aligned_pairs:
			pos = int(seq_pos[0]) + 1
			base = seq_pos[2]
			if base.islower() and int(query_qualities[pos-1]) >= 20:
				alt = str(pos) + ':' + str(pos) + ':0:' + str(seq_pos[1]+1) + '-' + str(seq_pos[1]+1) + ':SNP-' + base
				alterations.append(alt)
	else:
		for entry in cigar_tuple:
			if int(entry[0])==0:
				pos = pos + int(entry[1])
			elif int(entry[0])==1:
				index = pos-1
				ins_start = aligned_pairs[index-1][1]
				if ins_start==None:
					print(index)
					print(aligned_pairs)
				ins_seq = query_alignment_sequence[index-1:index-1+int(entry[1])]
				alt = str(pos) + ':' + str(pos) + ':' + str(entry[1]) + ':' + str(ins_start+1) + '-' + str(ins_start+1) + ':Ins-' + ins_seq
				alterations.append(alt)
				pos = pos + int(entry[1])
			elif int(entry[0])==2:
				index = pos - 1
				current_tuple = aligned_pairs[index]
				del_start = current_tuple[1]
				del_end = current_tuple[1]
				del_seq = ''
				while index < len(aligned_pairs) and current_tuple[0]==None:
					del_seq = del_seq + current_tuple[2]
					del_end = current_tuple[1]
					index += 1
					current_tuple = aligned_pairs[index]					
				end_pos = pos + int(entry[1])
				# print(end_pos)
				alt =  str(pos) + ':' + str(end_pos) + ':' + str(entry[1]) + ':' + str(del_start+1) + '-' + str(del_end+1) + ':Del-' + del_seq
				pos = end_pos
				alterations.append(alt)
	return alterations

def find_target_indices(target_site,reference_file,coordinates):
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	gene = ''
	if coordinates != 'full':
		chrom,pos = coordinates.split(':')
		start,end = pos.split('-')
	for record in SeqIO.parse(reference_file,'fasta'):
		if coordinates=='full':
			refDB[str(record.id)] = str(record.seq)
			targetDB[str(record.id)] = target_site
			gene = str(record.id)
			break
		elif chrom==str(record.id):
			refDB[str(record.id)] = str(record.seq)[int(start):int(end)]
			targetDB[str(record.id)] = target_site
			gene = str(record.id)
			break

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


def calculate_NHEJ_mutation_rate (sample_sam_file,control,reference_file,target_site,coordinates,display_name,output_file):

	# variables to store information
	sample_read_total = 0
	sample_mutation_count = 0
	sample_oof_count = 0
	control_cigar_strings = []

	# control bam file
	control_sam_file = 'processed/' + control + '_bwamem_sorted.bam'

	# get the start and end indexes of the target in the reference
	target_start,target_end,gene = find_target_indices(target_site,reference_file,coordinates)

	# get reference position
	if coordinates=='full':
		ref_pos = 0
	else:
		chrom,pos = coordinates.split(':')
		start,end = pos.split('-')
		ref_pos = int(start)-1

	# write header in output file
	output_file.write('Sample\tDisplay_Name\tTarget_Site\tMutated_Read_Count\tOOF_Mutated_Read_Count\tTotal_Read_Count\tNHEJ_Mutation_Rate\tOOF_Mutation_rate\tAlterations\n')

	# go through control file to determine "false" mutations
	ctrl_sam = pysam.AlignmentFile(control_sam_file,'rb')
	ctrl_total_count = 0
	ctrl_mut_freq = defaultdict(int)
	for read in ctrl_sam.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
			# check if the read has a mutation that involves the target sequence
			md_tag = read.get_tag('MD')
			ctrl_total_count += 1
			ctrl_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities,read.query_alignment_sequence, read.query_sequence)
			unique_entry = md_tag + '_' + read.cigarstring
			# if alteration is non zero
			valid_alteration = False
			if len(ctrl_alterations) > 0:
				for alt in ctrl_alterations:
					[start,end,indel_len,indel_coords,alt_type] = alt.split(':')
					if ((int(start) >= target_start-3 and int(start) <= target_end+3) or (int(end) >= target_start-3 and int(end) <= target_end+3)):
						valid_alteration = True
			# check if valid alteration
			if valid_alteration==True:
				if unique_entry not in ctrl_mut_freq:
					ctrl_mut_freq[unique_entry] = 1
					control_cigar_strings.append(unique_entry)
				else:
					ctrl_mut_freq[unique_entry] += 1

	# go through sample BAM file
	samfile = pysam.AlignmentFile(sample_sam_file,"rb")
	unique_alts = defaultdict(int)
	for read in samfile.fetch():
		if read.cigarstring != None and 'S' not in read.cigarstring and 'H' not in read.cigarstring and int(read.reference_start)==ref_pos:
			sample_read_total += 1
			md_tag = read.get_tag('MD')
			sample_alterations = process_cigar(read.cigartuples,md_tag,read.get_aligned_pairs(with_seq=True),read.query_qualities,read.query_alignment_sequence, read.query_sequence)
			unique_entry = md_tag + '_' + read.cigarstring
			# only count if the tuple is not in the control at a high frequency
			entry_freq = 0
			if unique_entry in ctrl_mut_freq and ctrl_total_count > 0:
				entry_freq = ctrl_mut_freq[unique_entry] / ctrl_total_count

			if len(sample_alterations) > 0 and entry_freq < 0.1:
				# if alteration is non zero
				valid_alteration = False
				for alt in sample_alterations:
					[start,end,indel_len,indel_coords,alt_type] = alt.split(':')
					if ((int(start) >= target_start-3 and int(start) <= target_end+3) or (int(end) >= target_start-3 and int(end) <= target_end+3)):
						valid_alteration = True
						final_alt = chrom + ':' + indel_coords + '&' + alt_type
						if final_alt not in unique_alts:
							unique_alts[final_alt] = 1
						else:
							unique_alts[final_alt] += 1
				if valid_alteration==True:
					sample_mutation_count += 1
					# see if it's OOF
					if int(indel_len) % 3 != 0:
						sample_oof_count += 1
	if sample_read_total > 0:
		nhej_rate = (sample_mutation_count / sample_read_total) * 100
		oof_rate = (sample_oof_count / sample_read_total) * 100
	else:
		nhej_rate = 'N/A'
		oof_rate = 'N/A'

	# update sample name
	sample_name = sample_sam_file.name
	sample_name = sample_name.replace('processed/','')
	sample_name = sample_name.replace('_bwamem_sorted.bam','')

	# write to output file
	alt_line = []
	sorted_alts = sorted(unique_alts.items(), key=lambda kv: kv[1], reverse=True)
	for u_alt in sorted_alts:
		proportion = round((unique_alts[u_alt[0]] / sample_read_total * 100),2)
		if proportion >= 20:
			alt = u_alt[0] + '&' + str(proportion)
			alt_line.append(alt)

	output_file.write(sample_name + '\t' + display_name + '\t' + target_site + '\t' + str(sample_mutation_count) + '\t' + str(sample_oof_count) + '\t' + str(sample_read_total) + '\t' + str(nhej_rate) + '\t' + str(oof_rate) + '\t' + ';'.join(alt_line) + '\n')

	# close file handles
	samfile.close()
	ctrl_sam.close()
	output_file.close()

#def calculate_base_editing_rate (sample_sam_file,control_sam_file,reference_file,target_site,output_file):

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s','--sample_sam_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_sam_file',required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	parser.add_argument('-i','--coordinates',required=True)
	parser.add_argument('-d','--display_name',required=True)
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		calculate_base_editing_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.coordinates,opts.display_name, opts.output_file)
	else:
		calculate_NHEJ_mutation_rate (opts.sample_sam_file,opts.control_sam_file,opts.reference_file,opts.target_site,opts.coordinates,opts.display_name, opts.output_file)
 
if __name__ == '__main__':
	main(sys.argv[1:])

