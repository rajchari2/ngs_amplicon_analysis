# test

import snakemake
import os
from collections import defaultdict

# get all of the samples
samples = config["samples"]
project_name = config["project"]
cell_type = config["cell_type"]
ngs_run = config["ngs_run"]
modality = config["modality"]
sample_list = []
control_list = []
reference_list = []
sample_to_reference = defaultdict(str)
sample_to_target_site = defaultdict(str)
sample_to_control = defaultdict(str)
sample_to_coordinates = defaultdict(str)
sample_to_display = defaultdict(str)
sample_to_trim = defaultdict(str)
coords_to_amp_name = defaultdict(str)
for sample in samples:
	for sample_name in sample:
		sample_list.append(sample_name)
		ref_file = sample[sample_name]['reference']
		display_name = sample[sample_name]['display_name']
		target_site = sample[sample_name]['target_site']
		control_sample = sample[sample_name]['control_sample']
		coordinates = sample[sample_name]['coordinates']
		amp_name = sample[sample_name]['amp_name']
		trim = sample[sample_name]['trim']
		if ref_file not in reference_list:
			reference_list.append(ref_file)
		sample_to_reference[sample_name] = ref_file
		sample_to_coordinates[sample_name] = coordinates
		sample_to_target_site[sample_name] = target_site
		sample_to_control[sample_name] = control_sample
		coords_to_amp_name[coordinates] = amp_name
		sample_to_display[sample_name] = display_name
		sample_to_trim[sample_name] = trim

		# populate control as well
		if control_sample not in control_list:
			control_list.append(control_sample)
			sample_to_reference[control_sample] = ref_file
			sample_to_target_site[control_sample] = target_site
			sample_to_coordinates[control_sample] = coordinates
			sample_to_trim[control_sample] = trim

# functions to get variables
def get_reference(wildcards):
	return sample_to_reference[wildcards.sample]

def get_target_site(wildcards):
	return sample_to_target_site[wildcards.sample]

def get_control_sample(wildcards):
	return sample_to_control[wildcards.sample]

def get_coordinates(wildcards):
	return sample_to_coordinates[wildcards.sample]

def get_trim(wildcards):
	return sample_to_trim[wildcards.sample]

def get_display_name(wildcards):
	return sample_to_display[wildcards.sample]

rule bwa_index:
	input:
		reference_file = '{reference}'
	output:
		bai = '{reference}.bwt'
	run:
		if os.path.exists(output.bai)==False:
			shell('bwa index {input.reference_file}')

rule build_indexes:
	input:
		ref_list = expand(rules.bwa_index.output.bai,reference=reference_list)
	output:
		index_file_created = 'processed/' + project_name + '_' + ngs_run + '_Indexes_Created.tab'
	run:
		with open(output.index_file_created, "w") as out:
			out.write('References indexed\n')

rule trim_data:
	input:
		r1 = 'data_files/{sample}_R1.fastq.gz',
		r2 = 'data_files/{sample}_R2.fastq.gz'
	params:
		trim = get_trim
	output:
		pr1 = 'data_files/{sample}_trimmed_R1.fastq.gz',
		pr2 = 'data_files/{sample}_trimmed_R2.fastq.gz'
	shell:
		'cutadapt -l {params.trim} -o {output.pr1} {input.r1} && cutadapt -l {params.trim} -o {output.pr2} {input.r2}'

rule flash_merge:
	input:
		r1 = rules.trim_data.output.pr1,
		r2 = rules.trim_data.output.pr2
	params:
		prefix = 'data_files/{sample}'
	output:
		merged = 'data_files/{sample}.extendedFrags.fastq',
	shell:
		'flash -M 100 -o {params.prefix} {input.r1} {input.r2}'

rule mask_seq_qual:
	input:
		flash_merged = rules.flash_merge.output.merged
	output:
		merged_filtered = 'data_files/{sample}.extendedFrags.filtered.fastq'
	shell:
		'python resources/mask_qual.py -i {input.flash_merged} -o {output.merged_filtered} -q 20'

rule bwa_mem:
	input:
		reference_file = get_reference,
		merged_fastq = rules.mask_seq_qual.output.merged_filtered,
	output:
		sam = 'processed/{sample}_bwamem.sam'
	shell:
		'bwa mem {input.reference_file} {input.merged_fastq} > {output.sam}'

rule filter_sam:
	input:
		sam_file = rules.bwa_mem.output.sam,
	params:
		coordinates = get_coordinates,
	output:
		filtered_sam_file = 'processed/{sample}' + '_bwamem_filtered.sam'
	shell:
		'python resources/filter_sam.py -i {input.sam_file} -c {params.coordinates} -o {output.filtered_sam_file}'

rule samtools_view:
	input:
		sam = rules.filter_sam.output.filtered_sam_file,
		reference_file = get_reference
	output:
		bam = 'processed/{sample}'  + '_bwamem.bam'
	shell:
		'samtools view -bt {input.reference_file} -o {output.bam} {input.sam}'

rule samtools_sort:
	input:
		bam = rules.samtools_view.output.bam
	output:
		sorted_bam = 'processed/{sample}' +'_bwamem_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.bam}'

rule bam_index:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam
	output:
		bai = 'processed/{sample}' + '_bwamem_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule calculate_mutation_by_SAM:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam,
		bam_index = rules.bam_index.output.bai,
		reference_file = get_reference
	params:
		target_site = get_target_site,
		cell_type_param = cell_type,
		modality_param = modality,
		control_file = get_control_sample,
		coordinates = get_coordinates,
		amp_name = coords_to_amp_name[coordinates],
		display_name = get_display_name
	output:
		mutation_summary_by_SAM = 'output/{sample}' + '_mutation_summary_by_SAM.tab'
	shell:
		'python resources/calculate_mutation_SAM.py -s {input.sorted_bam} -c {params.control_file} -r {input.reference_file} -t {params.target_site} -m {params.modality_param} -o {output.mutation_summary_by_SAM} -i {params.coordinates} -d {params.display_name}'

rule aggregate_output_SAM:
	input:
		mutation_file_list = sorted(expand(rules.calculate_mutation_by_SAM.output.mutation_summary_by_SAM,sample=sample_list))
	params:
		modality_param = modality
	output:
		mutation_summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary_bySAM.tab',
		low_coverage = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples_bySAM.tab',
	shell:
		'python resources/aggregate_rates.py -i {input.mutation_file_list} -o {output.mutation_summary_file} -m {params.modality_param} -l {output.low_coverage}'

rule graph_output_SAM:
	input:
		mutation_file = rules.aggregate_output_SAM.output.mutation_summary_file,
	params:
		modality_param = modality
	output:
		output_file = 'final_output/' + ngs_run + '_' + project_name + '_mutation_bar_plot_bySAM.png'
	shell:
		'python resources/make_bar_plot.py -i {input.mutation_file} -o {output.output_file} -m {params.modality_param}'

rule build_controls:
	input:
		pileup_list = expand(rules.samtools_sort.output.sorted_bam,sample=control_list),
		bais = expand(rules.bam_index.output.bai,sample=control_list)
	output:
		controls_created = project_name + '_' + ngs_run + '_Controls.tab'
	shell:
		'touch {output.controls_created}'

rule make_all:
	input:
		control_file = rules.build_controls.output.controls_created,
		report_file = rules.graph_output_SAM.output.output_file

rule clean_run:
	params:
		controls_created = project_name + '_' + ngs_run + '_Controls.tab',
		sams = ' '.join(sorted(expand(rules.bwa_mem.output.sam,sample=sample_list))),
		filtered_sams = ' '.join(sorted(expand(rules.filter_sam.output.filtered_sam_file,sample=sample_list))),
		bams = ' '.join(sorted(expand(rules.samtools_view.output.bam,sample=sample_list))),
		sorted_bams = ' '.join(sorted(expand(rules.samtools_sort.output.sorted_bam,sample=sample_list))),
		ctrl_sams = ' '.join(sorted(expand(rules.bwa_mem.output.sam,sample=control_list))),
		ctrl_bams = ' '.join(sorted(expand(rules.samtools_view.output.bam,sample=control_list))),
		ctrl_sorted_bams = ' '.join(sorted(expand(rules.samtools_sort.output.sorted_bam,sample=control_list))),
		bais = ' '.join(sorted(expand(rules.bam_index.output.bai,sample=sample_list))),
		summaries = ' '.join(sorted(expand(rules.calculate_mutation_by_SAM.output.mutation_summary_by_SAM,sample=sample_list))),
		summary_file_SAM = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary_bySAM.tab',
		low_coverage_SAM = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples_bySAM.tab',
		output_file_SAM = 'final_output/' + ngs_run + '_' + project_name + '_mutation_bar_plot_bySAM.png',

	shell:
		'rm -f {params.index_file_created} {params.controls_created} {params.sams} {params.filtered_sams} {params.sorted_bams} {params.bams} {params.summaries} {params.ctrl_sams} {params.ctrl_bams} {params.ctrl_sorted_bams} {params.bais} {params.summary_file_SAM} {params.low_coverage_SAM} {params.output_file_SAM}'

