import sys
import os
import re


## Extract variables from the config file
project_name = config['project_name']
workflow_path = config['workflow_path']

input_file = config['input_file']
output_folder = config['output_folder']

vcf_file = config['vcf_file']

memory = config['memory']
threads = config['threads']


## Read all region string in input file into an array
regions = []
with open(input_file, 'r') as reader:
	for line in reader:
		line = str(line).strip("\n").strip("\r").strip("\r\n")
		line = re.sub('\t', '__', line)
		regions.append(line)


## Print variables
print("project_name: ", project_name)
print("workflow_path: ", workflow_path)

print("input_file: ", input_file)
print("output_folder: ", output_folder)

print("vcf_file: ", vcf_file)

print("memory: ", memory)
print("threads: ", threads)

print("regions: ", regions)


## Collect all outputs
rule all:
	input:
		expand(os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.recode.vcf.gz'), region=regions),
		expand(os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.ped'), region=regions),
		expand(os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.map'), region=regions),
		expand(os.path.join(os.path.abspath(output_folder), "Haploview", '{region}.LD'), region=regions)


## Run VCFtools to subset and convert to plink in a snakemake rule
rule vcftools_subset_and_convert_to_plink:
	params:
		in_file = vcf_file,
		chrom_start_end = lambda wildcards: wildcards.region.split("__"),
		output_prefix = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}'),
		out_folder = os.path.join(os.path.abspath(output_folder), "VCFtools")
	output:
		recode_vcf_file = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.recode.vcf.gz'),
		ped_file = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.ped'),
		map_file = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.map')
	resources:
		memory = memory,
		threads = threads
	shell:
		"""
		mkdir -p {params.out_folder};
		vcftools --gzvcf {params.in_file} --chr {params.chrom_start_end[0]} --from-bp {params.chrom_start_end[1]} --to-bp {params.chrom_start_end[2]} --recode --stdout | bgzip > {output.recode_vcf_file};
		vcftools --gzvcf {output.recode_vcf_file} --plink --out {params.output_prefix};
		"""

rule haploview:
	input:
		ped_file = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.ped'),
		map_file = os.path.join(os.path.abspath(output_folder), "VCFtools", '{region}.map')
	params:
		output_prefix = os.path.join(os.path.abspath(output_folder), "Haploview", '{region}'),
		out_folder = os.path.join(os.path.abspath(output_folder), "Haploview")
	output:
		ld_file = os.path.join(os.path.abspath(output_folder), "Haploview", '{region}.LD')
	resources:
		memory = memory,
		threads = threads
	shell:
		"""
		mkdir -p {params.out_folder};
		java -Xmx{resources.memory}g -jar {workflow_path}/tools/Haploview.jar -n -out {params.output_prefix} -pedfile {input.ped_file} -map {input.map_file} -skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05
		"""
