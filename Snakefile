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
reference_list = []
control_list = []
sample_to_reference = defaultdict(str)
sample_to_target_site = defaultdict(str)
sample_to_control = defaultdict(str)
for sample in samples:
	for sample_name in sample:
		sample_list.append(sample_name)
		ref_file = sample[sample_name]['reference']
		target_site = sample[sample_name]['target_site']
		control_sample = sample[sample_name]['control_sample']
		if ref_file not in reference_list:
			reference_list.append(ref_file)
		sample_to_reference[sample_name] = 'db/' + ref_file
		sample_to_target_site[sample_name] = target_site
		sample_to_control[sample_name] = control_sample
		# populate control as well
		if control_sample not in control_list:
			control_list.append(control_sample)
			sample_to_reference[control_sample] = 'db/' + ref_file

# functions to get variables
def get_reference(wildcards):
	return sample_to_reference[wildcards.sample]

def get_target_site(wildcards):
	return sample_to_target_site[wildcards.sample]

def get_control_sample(wildcards):
	return sample_to_control[wildcards.sample]

rule bwa_index:
	input:
		reference_file = 'db/{reference}'
	output:
		bai = 'db/{reference}.bwt'
	shell:
		'bwa index {input.reference_file}'

rule build_indexes:
	input:
		ref_list = expand(rules.bwa_index.output.bai,reference=reference_list)
	output:
		index_file_created = 'db/' + project_name + '_' + ngs_run + '_Indexes_Created.tab'
	shell:
		'cat db/*.fasta > {output.index_file_created}'

rule flash_merge:
	input:
		r1 = 'data_files/{sample}_R1.fastq.gz',
		r2 = 'data_files/{sample}_R2.fastq.gz',
	params:
		prefix = 'data_files/{sample}'
	output:
		merged = 'data_files/{sample}.extendedFrags.fastq',
	shell:
		'flash -M 100 -o {params.prefix} {input.r1} {input.r2}'

rule bwa_mem:
	input:
		reference_file = get_reference,
		merged_fastq = rules.flash_merge.output.merged
	output:
		sam = 'processed/{sample}_bwamem.sam'
	shell:
		'bwa mem {input.reference_file} {input.merged_fastq} > {output.sam}'

rule samtools_view:
	input:
		sam = rules.bwa_mem.output.sam,
		reference_file = get_reference
	output:
		bam = 'processed/{sample}_bwamem.bam'
	shell:
		'samtools view -bt {input.reference_file} -o {output.bam} {input.sam}'

rule samtools_sort:
	input:
		bam = rules.samtools_view.output.bam
	output:
		sorted_bam = 'processed/{sample}_bwamem_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.bam}'

rule bam_index:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam
	output:
		bai = 'processed/{sample}_bwamem_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule samtools_pileup:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam,
		reference_file = get_reference,
	output:
		pileup = 'processed/{sample}_mpileup.tab'
	shell:
		'samtools mpileup -C50 -f {input.reference_file} -o {output.pileup} {input.sorted_bam}'

rule calculate_mutation_rate:
	input:
		pileup_file = rules.samtools_pileup.output.pileup,
		reference_file = get_reference,
	params:
		illumina_scoring = 'resources/Illumina.Scoring.txt',
		target_site = get_target_site,
		cell_type_param = cell_type,
		modality_param = modality,
		control_file = get_control_sample
	output:
		mutation_summary = 'output/{sample}_mutation_summary.tab',
		diversity_output = 'output/{sample}_mutation_diversity.tab'
	shell:
		'python resources/calculate_mutation_rates.py -p {input.pileup_file} -c {params.control_file} -i {params.illumina_scoring} -r {input.reference_file} -t {params.target_site} -o {output.mutation_summary} -m {params.modality_param} -d {output.diversity_output}'

rule aggregate_output:
	input:
		file_list_mutation = sorted(expand(rules.calculate_mutation_rate.output.mutation_summary,sample=sample_list)),
		file_list_diversity = sorted(expand(rules.calculate_mutation_rate.output.diversity_output,sample=sample_list)),
	params:
		modality_param = modality
	output:
		mutation_summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary.tab',
		diversity_summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_diversity_summary.tab',
		low_coverage = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples.tab'
	shell:
		'python resources/aggregate_rates.py -i {input.file_list_mutation} -d {input.file_list_diversity} -o {output.summary_file} -m {params.modality_param} -l {output.low_coverage} -s {output.diversity_summary_file}'

rule graph_output:
	input:
		summary_file = rules.aggregate_output.output.mutation_summary_file
	params:
		modality_param = modality
	output:
		output_file = 'final_output/' + ngs_run + '_' + project_name + '_mutation_bar_plot.png'
	shell:
		'python resources/make_bar_plot.py --input_file {input.summary_file} --output_file {output.output_file} --modality {params.modality_param}'

rule calculate_mutation_by_SAM:
	input:
		sorted_bam = rules.samtools_sort.output.sorted_bam,
		bam_index = rules.bam_index.output.bai,
		reference_file = get_reference
	params:
		target_site = get_target_site,
		cell_type_param = cell_type,
		modality_param = modality,
		control_file = get_control_sample
	output:
		mutation_summary_by_SAM = 'output/{sample}_mutation_summary_by_SAM.tab',
		diversity_output_by_SAM = 'output/{sample}_mutation_diversity_by_SAM.tab'
	shell:
		'python resources/calculate_mutation_SAM.py -s {input.sorted_bam} -c {params.control_file} -r {input.reference_file} -t {params.target_site} -m {params.modality_param} -o {output.mutation_summary_by_SAM} -d {output.diversity_output_by_SAM}'

rule aggregate_output_SAM:
	input:
		mutation_file_list = sorted(expand(rules.calculate_mutation_by_SAM.output.mutation_summary_by_SAM,sample=sample_list)),
		diversity_file_list = sorted(expand(rules.calculate_mutation_by_SAM.output.diversity_output_by_SAM,sample=sample_list)),
	params:
		modality_param = modality
	output:
		mutation_summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary_bySAM.tab',
		low_coverage = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples_bySAM.tab',
		diversity_summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_diversity_summary_bySAM.tab',
	shell:
		'python resources/aggregate_rates.py -i {input.mutation_file_list} -d {input.diversity_file_list} -o {output.mutation_summary_file} -m {params.modality_param} -l {output.low_coverage} -s {output.diversity_summary_file}'

rule graph_output_SAM:
	input:
		mutation_file = rules.aggregate_output_SAM.output.mutation_summary_file,
		diversity_file = rules.aggregate_output_SAM.output.diversity_summary_file
	params:
		modality_param = modality
	output:
		output_file = 'final_output/' + ngs_run + '_' + project_name + '_mutatiion_bar_plot_bySAM.png',
		diversity_graph = 'final_output/' + ngs_run + '_' + project_name + '_diversity_bar_plot_bySAM.png',
	shell:
		'python resources/make_bar_plot.py -i {input.mutation_file} -o {output.output_file} -m {params.modality_param} -d {input.diversity_file} -g {output.diversity_graph}'


rule build_controls:
	input:
		pileup_list = expand(rules.samtools_pileup.output.pileup,sample=control_list),
		bais = expand(rules.bam_index.output.bai,sample=control_list)
	output:
		controls_created = project_name + '_' + ngs_run + '_Controls.tab'
	shell:
		'touch {output.controls_created}'

# get all the reference files and run bwa index
rule make_all_V1:
	input:
		index_file = rules.build_indexes.output.index_file_created,
		control_file = rules.build_controls.output.controls_created,
		report_file = rules.graph_output.output.output_file

rule make_all_V2:
	input:
		index_file = rules.build_indexes.output.index_file_created,
		control_file = rules.build_controls.output.controls_created,
		report_file = rules.graph_output_SAM.output.output_file

rule clean_run:
	params:
		index_file_created = 'db/' + project_name + '_' + ngs_run + '_Indexes_Created.tab',
		summary_file = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary.tab',
		low_coverage = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples.tab',
		output_file = 'final_output/' + ngs_run + '_' + project_name + '_bar_plot.png',
		controls_created = project_name + '_' + ngs_run + '_Controls.tab',
		sams = ' '.join(sorted(expand(rules.bwa_mem.output.sam,sample=sample_list))),
		bams = ' '.join(sorted(expand(rules.samtools_view.output.bam,sample=sample_list))),
		sorted_bams = ' '.join(sorted(expand(rules.samtools_sort.output.sorted_bam,sample=sample_list))),
		pileups = ' '.join(sorted(expand(rules.samtools_pileup.output.pileup,sample=sample_list))),
		summaries = ' '.join(sorted(expand(rules.calculate_mutation_rate.output.mutation_summary,sample=sample_list))),
		ctrl_sams = ' '.join(sorted(expand(rules.bwa_mem.output.sam,sample=control_list))),
		ctrl_bams = ' '.join(sorted(expand(rules.samtools_view.output.bam,sample=control_list))),
		ctrl_sorted_bams = ' '.join(sorted(expand(rules.samtools_sort.output.sorted_bam,sample=control_list))),
		ctrl_pileups = ' '.join(sorted(expand(rules.samtools_pileup.output.pileup,sample=control_list))),
		bais = ' '.join(sorted(expand(rules.bam_index.output.bai,sample=sample_list))),
		summaries_SAM = ' '.join(sorted(expand(rules.calculate_mutation_by_SAM.output.mutation_summary_by_SAM,sample=sample_list))),
		summary_file_SAM = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_mutation_summary_bySAM.tab',
		low_coverage_SAM = 'final_output/' + ngs_run + '_' + project_name + '_' + cell_type + '_low_coverage_samples_bySAM.tab',
		output_file_SAM = 'final_output/' + ngs_run + '_' + project_name + '_bar_plot_bySAM.png'
	shell:
		'rm -f {params.index_file_created} {params.summary_file} {params.low_coverage} {params.output_file} {params.controls_created} {params.sams} {params.sorted_bams} {params.bams} {params.pileups} {params.summaries} {params.ctrl_sams} {params.ctrl_bams} {params.ctrl_sorted_bams} {params.ctrl_pileups} {params.bais} {params.summaries_SAM} {params.summary_file_SAM} {params.low_coverage_SAM} {params.output_file_SAM}'

