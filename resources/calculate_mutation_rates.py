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

# store the scoring table
def read_illumina_table(illumina_scoring_file):
	scoringTable = defaultdict(int)
	for line in illumina_scoring_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		scoringTable[parts[0]] = int(parts[1])
	illumina_scoring_file.close()
	return scoringTable

# complement table
complement = defaultdict(dict)
complement['A'] = 'T'
complement['C'] = 'G'
complement['T'] = 'A'
complement['G'] = 'C'

# function to parse a pileup string
def processPileupString(pileupString,pileupQuality,scoringTable,refBase,modality):
	# set the defaults
	baseList = []
	alterationList = []
	index = 0
	seqIndex = 0
	pileupString = pileupString.upper()
	while index < len(pileupString):
		flagged = 0
		indelStatus = 0
		# case 1: index is a start of a segment
		if pileupString[index]=='^':
			# move two characters to the call
			flagged = 1
			index += 2

		if index < len(pileupString):
			# now check if you have a match or indel
			if flagged==1:
				seqIndex += 1
			else:
				if pileupString[index]=='.':
					baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index]==',':
					baseList.append(complement[refBase])
					seqIndex += 1
				elif pileupString[index]=='*':
					if scoringTable[pileupQuality[seqIndex]] >= 28:
						baseList.append('-')
					else:
						baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index].upper()=='A' or pileupString[index].upper()=='C' or pileupString[index].upper()=='T' or pileupString[index].upper()=='G':
					if scoringTable[pileupQuality[seqIndex]] >= 28:
						baseList.append(pileupString[index])
						alteration = refBase + '>' + pileupString[index].upper() + ':' + str(seqIndex)
						if modality=='ABE7.10' or modality=='BE4':
							alterationList.append(alteration)
					else:
						baseList.append(refBase)
					seqIndex += 1
				elif pileupString[index]=='+':
					# check if the next digit is a number
					if pileupString[index+2].isdigit():
						insSize = int(pileupString[index+1:index+3])
						insertion = pileupString[index:index+3+insSize] + ':' + str(seqIndex-1)
						index +=  insSize + 2
					else:
						insSize = int(pileupString[index+1])
						insertion = pileupString[index:index+2+insSize] + ':' + str(seqIndex-1)
						index += insSize + 1
					alterationList.append(insertion)

				elif pileupString[index]=='-':
					# check if the next digit is a number
					if pileupString[index+2].isdigit():
						delStart = index
						delEnd = index + int(pileupString[index+1:index+3]) + 3
						deletion = pileupString[delStart:delEnd] + ':' + str(seqIndex-1)
						index += int(pileupString[index+1:index+3]) + 2
					else:
						delStart = index
						delEnd = index + int(pileupString[index+1]) + 2
						deletion = pileupString[delStart:delEnd] + ':' + str(seqIndex-1)
						index += int(pileupString[index+1]) + 1
					alterationList.append(deletion)
					indelStatus = 1								
		index += 1
	return alterationList

# function to determine relevant alteration
def determineAlterationStatus(alteration, position, targetSiteStart, targetSiteEnd):
	countAlteration = 0
	if '+' in alteration or '>' in alteration:
		if position >= targetSiteStart and position <= targetSiteEnd:
			countAlteration = 1
	elif '-' in alteration:
		if alteration[2].isdigit():
			size = int(alteration[1:3])
		else:
			size = int(alteration[1:2])
		if (position >= targetSiteStart and position <= targetSiteEnd) or (position + size >= targetSiteStart and position + size <= targetSiteEnd):
			countAlteration = 1
	return countAlteration


def calculate_base_editing_rate (pileup_file,modality,illumina_scoring_file,reference_file,target_site,output_file):
	# first get the scoring table
	scoringTable = read_illumina_table(illumina_scoring_file)
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	for record in SeqIO.parse(reference_file,'fasta'):
		refDB[str(record.id)] = str(record.seq)
		targetDB[str(record.id)] = target_site
	# write the header for the output file
	output_file.write('File\tGene\tTotalCount\tPos1\tPos2\tPos3\tPos4\tPos5\tPos6\tPos7\tPos8\tPos9\tPos10\tPos11\tPos12\tPos13\tPos14\tPos15\tPos16\tPos17\tPos18\tPos19\tPos20\n')
	# check the modality if it's BE4, looking for C and makes it a T. If it's ABE7.10, looking for A and turning it into G
	if modality=='BE4':
		sense = 'C>T'
		antisense = 'G>A'
	else:
		sense = 'A>G'
		antisense = 'T>C'
	# make permanent variables for spacer_start and spacer_end

	# go through the pileup file
	base_change_list = defaultdict(dict)
	avgDepth = defaultdict(list)
	# initialize for the gene
	avgDepth[gene] = []

	for line in pileup_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		gene = parts[0]
		position = int(parts[1])
		refBase = parts[2].upper()
		if gene not in base_change_list:
			base_change_list[gene] = defaultdict(int)
		
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
			spacer_end = rcIndex + 4
			spacer_start = spacer_end + 19
			direction = 'antisense'
		else:
			spacer_start = index + 1
			spacer_end = spacer_start + 19
			direction = 'sense'
		# store alterations for sample 1
		depth = int(parts[3])
		pileupString = parts[4]
		pileupQuality = parts[5]
		alterationSet = processPileupString(pileupString,pileupQuality,scoringTable,refBase,modality)

		# each alteration now holds a position
		for alt in alterationSet:
			altParts = alt.split(':')
			if direction=='sense':
				if altParts[0]==sense:
					if int(parts[1]) not in base_change_list[gene]:
						base_change_list[gene][position] = 1
					else:
						base_change_list[gene][position] += 1
			else:
				if altParts[0]==antisense:
					if int(parts[1]) not in base_change_list[gene]:
						base_change_list[gene][position] = 1
					else:
						base_change_list[gene][position] += 1				

		# adjust the depth
		adjustment = 0
		for c in pileupQuality:
			if scoringTable[c] < 28:
				adjustment += 1

		depth = depth - adjustment
		avgDepth[gene].append(depth)

	# close handles
	pileup_file.close()
	print('Finished: ' + pileup_file.name)

	# go through each gene
	for gene in avgDepth:
		line_to_write = pileup_file.name + '\t' + gene + '\t' + str(max(avgDepth[gene]))
		if gene not in base_change_list:
			line_to_write = line_to_write + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0'
		else:
			if direction=='sense':
				current_pos = spacer_start
				while current_pos <= spacer_end:
					if current_pos in base_change_list[gene]:
						value = base_change_list[gene][current_pos]
					else:
						value = 0
					line_to_write = line_to_write + '\t' + str(value)
					current_pos += 1
			else:
				current_pos = spacer_start
				while current_pos >= spacer_end:
					if current_pos in base_change_list[gene]:
						value = base_change_list[gene][current_pos]
					else:
						value = 0
					line_to_write = line_to_write + '\t' + str(value)
					current_pos -= 1
		output_file.write(line_to_write + '\n')
	output_file.close()


def calculate_rate(pileup_file,control,illumina_scoring_file,reference_file,target_site,output_file,modality):
	# first get the scoring table
	scoringTable = read_illumina_table(illumina_scoring_file)
	
	# store reference sequence
	refDB = defaultdict(str)
	targetDB = defaultdict(str)
	gene = ''
	for record in SeqIO.parse(reference_file,'fasta'):
		refDB[str(record.id)] = str(record.seq)
		targetDB[str(record.id)] = target_site
		gene = str(record.id)
	
	# write headers to the output file
	output_file.write('File\tGene\tMutatedCount\tTotalCount\tMutation_Percentage\n')

	# structure to hold the mutated reads -> gene -> sample -> read# -> [alterationList]
	avgDepth = defaultdict(list)
	mutatedReads = defaultdict(dict)
	controlReads = defaultdict(list)

	# initialize for the gene
	avgDepth[gene] = []
	mutatedReads[gene] = defaultdict(list)
	controlReads[gene] = []

	control_file = 'processed/' + control + '_mpileup.tab'
	control_file_handle = open(control_file,'r')

	# go through control sample
	for line in control_file_handle:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		gene = parts[0]
		position = int(parts[1])
		refBase = parts[2].upper()	

		# store alterations for sample 1
		depth = int(parts[3])
		pileupString = parts[4]
		pileupQuality = parts[5]
		alterationSet = processPileupString(pileupString,pileupQuality,scoringTable,refBase,modality)

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

		# each alteration now holds a position
		for alt in alterationSet:
			altParts = alt.split(':')
			mutEntry = altParts[0] + ':' + parts[1]

			# identify the alteration span, three cases: insertions, deletions and SNPs
			countAlteration = determineAlterationStatus(altParts[0], position, startIndex, endIndex)
			#print 'Start: ' + str(startIndex) + ', ' + 'End index: ' + str(endIndex)
			if countAlteration==1 and '>' not in altParts[0]:
				if mutEntry not in controlReads[gene]:
					controlReads[gene].append(mutEntry)	
	control_file_handle.close()

	# go through the pileup file
	for line in pileup_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		gene = parts[0]
		position = int(parts[1])
		refBase = parts[2].upper()	

		# store alterations for sample 1
		depth = int(parts[3])
		pileupString = parts[4]
		pileupQuality = parts[5]
		alterationSet = processPileupString(pileupString,pileupQuality,scoringTable,refBase,modality)

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

		# each alteration now holds a position
		for alt in alterationSet:
			altParts = alt.split(':')
			mutEntry = altParts[0] + ':' + parts[1]

			# identify the alteration span, three cases: insertions, deletions and SNPs
			countAlteration = determineAlterationStatus(altParts[0], position, startIndex, endIndex)
			#print 'Start: ' + str(startIndex) + ', ' + 'End index: ' + str(endIndex)
			if countAlteration==1 and '>' not in altParts[0] and mutEntry not in controlReads[gene]:
				if altParts[1] not in mutatedReads[gene]:
					mutatedReads[gene][altParts[1]] = []
				mutatedReads[gene][altParts[1]].append(mutEntry)

		# adjust the depth
		adjustment = 0
		for c in pileupQuality:
			if scoringTable[c] < 28:
				adjustment += 1

		depth = depth - adjustment
		avgDepth[gene].append(depth)

	pileup_file.close()
	print('Finished: ' + pileup_file.name)

	# go through all the reads and count whether they are mutated
	for gene in avgDepth:
		mutatedCount = 0
		totalCount = 0
		for readID in mutatedReads[gene]:
			if len(mutatedReads[gene][readID]) > 0:
				mutatedCount += 1
		# see if no data was added
		if len(avgDepth[gene])==0:
			avgDepth[gene].append(0)

		if max(avgDepth[gene]) >= 100:
			rate = (mutatedCount / max(avgDepth[gene])) * 100
		else:
			rate = 'N/A'
		output_file.write(pileup_file.name + '\t' + gene + '\t' + str(mutatedCount) + '\t' + str(max(avgDepth[gene])) + '\t' + str(rate) + '\n')
	output_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-p','--pileup_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-c','--control_file',required=False)
	parser.add_argument('-i','--illumina_scoring_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-r','--reference_file',required=True)
	parser.add_argument('-t','--target_site',required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--output_file',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	if opts.modality=='ABE7.10' or opts.modality=='BE4':
		calculate_base_editing_rate (opts.pileup_file,opts.control_file,opts.modality,opts.illumina_scoring_file,opts.reference_file,opts.target_site,opts.output_file)
	else:
		calculate_rate (opts.pileup_file,opts.control_file,opts.illumina_scoring_file,opts.reference_file,opts.target_site,opts.output_file,opts.modality)
 
if __name__ == '__main__':
	main(sys.argv[1:])


