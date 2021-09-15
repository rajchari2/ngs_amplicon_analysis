# script to make a yaml file from a tab-delimited sample file
# updated 09/18/18 to compute control subtracted rates

# input file: [tab-delimited - filename {tab} display-name {tab} reference file {tab} coordinates {tab} amplicon-name {tab} target site {tab} control_sample
# if coordinates === full -> entire sequence
# output file: yaml file which can run in the pipeline

# import appropriate libraries
import yaml
import sys
import argparse
import re
import glob
import subprocess
from collections import defaultdict

def makeYAML (input_file,project_name,ngs_run,cell_type,modality,output_file,directory_name):
	# distict names
	distinct_names = defaultdict(str)

	# dictionary to hold yaml
	yaml_dict = {}
	yaml_dict['samples'] = []
	yaml_dict['project'] = project_name
	yaml_dict['ngs_run'] = ngs_run
	yaml_dict['modality'] = modality
	yaml_dict['cell_type'] = cell_type

	# collect control files
	ctrl_list = []

	for line in input_file:
		line = line.rstrip('\r\n')
		parts = line.split('\t')

		# variables to store fields
		filename = parts[0]
		display_name = parts[1]
		reference_file = parts[2]
		coordinates = parts[3]
		amplicon_name = parts[4]
		target_site = parts[5]
		control_sample = parts[6]

		# get the size of the amplicon
		chrom,coords = coordinates.split(':')
		start,end = coords.split('-')
		size = int(end) - int(start)
		if size <= 290:
			trim = '150'
		else:
			trim = '250'

		# add to control file list
		if control_sample not in ctrl_list:
			ctrl_list.append(control_sample)
		
		# add to dictionary
		sample_dict = {}
		sample_name = filename + '-' + amplicon_name

		# check if this is distinct
		if sample_name not in distinct_names:
			distinct_names[sample_name] = filename
		else:
			print('Sample name: ' + sample_name + ' not distinct')
			exit()
		sample_dict[sample_name] = {}
		sample_dict[sample_name]['display_name'] = display_name
		sample_dict[sample_name]['reference'] = reference_file
		sample_dict[sample_name]['coordinates'] = coordinates
		sample_dict[sample_name]['amp_name'] = amplicon_name
		sample_dict[sample_name]['target_site'] = target_site
		sample_dict[sample_name]['control_sample'] = control_sample
		sample_dict[sample_name]['trim'] = trim

		# add to yaml dict
		yaml_dict['samples'].append(sample_dict)

	# now copy files
	for sample in distinct_names:
		filename = distinct_names[sample]

		# go through and change the filenames
		path_to_check = directory_name + '/' + filename + '_*.fastq.gz'
		file_list = sorted(glob.glob(path_to_check))

		# and change each file
		cpCommand_R1 = 'cp ' + file_list[0] + ' ' + directory_name + '/' + sample + '_R1.fastq.gz'
		p = subprocess.Popen(cpCommand_R1,shell=True)
		p.communicate()
		cpCommand_R2 = 'cp ' + file_list[1] + ' ' + directory_name + '/' + sample + '_R2.fastq.gz'
		p = subprocess.Popen(cpCommand_R2,shell=True)
		p.communicate()

	# move the control sample files
	for ctrl_sample in ctrl_list:
		path_to_check = directory_name + '/' + ctrl_sample + '*.fastq.gz'
		file_list = sorted(glob.glob(path_to_check))
		# and change each file
		mvCommand_R1 = 'cp ' + file_list[0] + ' ' + directory_name + '/' + ctrl_sample + '_R1.fastq.gz'
		p = subprocess.Popen(mvCommand_R1,shell=True)
		p.communicate()
		mvCommand_R2 = 'cp ' + file_list[1] + ' ' + directory_name + '/' + ctrl_sample + '_R2.fastq.gz'
		p = subprocess.Popen(mvCommand_R2,shell=True)
		p.communicate()		

	# write final yaml file
	with open(output_file, 'w') as yaml_file:
		yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    # close everything
	input_file.close()
	yaml_file.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-p','--project_name',required=True)
	parser.add_argument('-n','--ngs_run',required=True)
	parser.add_argument('-c','--cell_type',required=True)
	parser.add_argument('-d','--data_directory',required=True)
	parser.add_argument('-m','--modality',required=True)
	parser.add_argument('-o','--output_file',required=True)
	opts = parser.parse_args(argv)
	makeYAML(opts.input_file,opts.project_name,opts.ngs_run,opts.cell_type,opts.modality,opts.output_file,opts.data_directory)
 
if __name__ == '__main__':
	main(sys.argv[1:])