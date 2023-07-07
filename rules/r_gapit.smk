import sys
import os
import re


## Extract variables from the config file
project_name = config['project_name']
workflow_path = config['workflow_path']

input_folder = config['input_folder']
output_folder = config['output_folder']

genotype_hapmap = config['genotype_hapmap']
genotype_data = config['genotype_data']
genotype_map = config['genotype_map']

kinship = config['kinship']
z_matrix = config['z_matrix']
corvariance_matrix = config['corvariance_matrix']

snp_maf = config['snp_maf']
model = config['model']
pca_total = config['pca_total']

ulimit = config['ulimit']
memory = config['memory']
threads = config['threads']


## Extract samples from files in the input folder
samples = [re.sub('(.txt)|(.csv)', '', f) for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, f))]


## Print variables
print("project_name: ", project_name)
print("workflow_path: ", workflow_path)

print("input_folder: ", input_folder)
print("output_folder: ", output_folder)

print("genotype_hapmap: ", genotype_hapmap)
print("genotype_data: ", genotype_data)
print("genotype_map: ", genotype_map)

print("kinship: ", kinship)
print("z_matrix: ", z_matrix)
print("corvariance_matrix: ", corvariance_matrix)

print("snp_maf: ", snp_maf)
print("model: ", model)
print("pca_total: ", pca_total)

print("ulimit: ", ulimit)
print("memory: ", memory)
print("threads: ", threads)

print("samples: ", samples)


## Collect all outputs
rule all:
	input:
		expand(os.path.join(os.path.abspath(output_folder), "GAPIT_log", '{sample}.log'), sample=samples)


## Run GAPIT in a snakemake rule
rule r_gapit:
	input:
		in_file = os.path.join(os.path.abspath(input_folder),'{sample}.txt')
	params:
		genotype_hapmap = genotype_hapmap,
		genotype_data = genotype_data,
		genotype_map = genotype_map,
		kinship = kinship,
		z_matrix = z_matrix,
		corvariance_matrix = corvariance_matrix,
		snp_maf = snp_maf,
		model = model,
		pca_total = pca_total,
		out_folder = os.path.abspath(output_folder),
		ulimit = ulimit
	output:
		out_file = os.path.join(os.path.abspath(output_folder), "GAPIT_log", '{sample}.log')
	resources:
		memory = memory,
		threads = threads
	shell:
		"""
		mkdir -p {params.out_folder};
		Rscript {workflow_path}/scripts/R/GAPIT.R -i {input.in_file} -o {params.out_folder} --genotype_hapmap {params.genotype_hapmap} --genotype_data {params.genotype_data} --genotype_map {params.genotype_map} --kinship {params.kinship} --z_matrix {params.z_matrix} --corvariance_matrix {params.corvariance_matrix} --snp_maf {params.snp_maf} --model {params.model} --pca_total {params.pca_total} > {output.out_file}
		"""
