import sys
import os
import re


## Extract variables from the config file
project_name = config['project_name']
workflow_path = config['workflow_path']

input_folder = config['input_folder']
output_folder = config['output_folder']

chromosome_array = config['chromosome_array']

vcf_folder = config['vcf_folder']
vcf_file_extension = config['vcf_file_extension']

genotype_hapmap_folder = config['genotype_hapmap_folder']
genotype_hapmap_file_extension = config['genotype_hapmap_file_extension']
genotype_data_folder = config['genotype_data_folder']
genotype_data_file_extension = config['genotype_data_file_extension']
genotype_map_folder = config['genotype_map_folder']
genotype_map_file_extension = config['genotype_map_file_extension']

kinship_folder = config['kinship_folder']
kinship_file_extension = config['kinship_file_extension']
covariance_matrix_folder = config['covariance_matrix_folder']
covariance_matrix_file_extension = config['covariance_matrix_file_extension']

snp_maf = config['snp_maf']
model = config['model']
pca_total = config['pca_total']

ulimit = config['ulimit']
memory = config['memory']
threads = config['threads']

## Extract samples from files in the input folder
samples = [
    re.sub('(.txt)|(.csv)','',f) for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder,f))
]

## Print variables
print("project_name: ",project_name)
print("workflow_path: ",workflow_path)

print("input_folder: ",input_folder)
print("output_folder: ",output_folder)

print("chromosome_array: ",chromosome_array)

print("vcf_folder: ",vcf_folder)
print("vcf_file_extension: ",vcf_file_extension)

print("genotype_hapmap_folder: ",genotype_hapmap_folder)
print("genotype_hapmap_file_extension: ",genotype_hapmap_file_extension)
print("genotype_data_folder: ",genotype_data_folder)
print("genotype_data_file_extension: ",genotype_data_file_extension)
print("genotype_map_folder: ",genotype_map_folder)
print("genotype_map_file_extension: ",genotype_map_file_extension)

print("kinship_folder: ",kinship_folder)
print("kinship_file_extension: ",kinship_file_extension)
print("covariance_matrix_folder: ",covariance_matrix_folder)
print("covariance_matrix_file_extension: ",covariance_matrix_file_extension)

print("snp_maf: ",snp_maf)
print("model: ",model)
print("pca_total: ",pca_total)

print("ulimit: ",ulimit)
print("memory: ",memory)
print("threads: ",threads)

print("samples: ",samples)


## Collect all outputs
rule all:
    input:
        expand(os.path.join(os.path.abspath(output_folder),'{chromosome}',"GAPIT_out",'{sample}.log'),chromosome=chromosome_array,sample=samples)


## Run GAPIT in a snakemake rule
rule r_gapit:
    input:
        in_file=os.path.join(os.path.abspath(input_folder),'{sample}.txt'),
        vcf_file=os.path.join(os.path.abspath(vcf_folder),'{chromosome}' + str(vcf_file_extension))
    params:
        sample='{sample}',
        chromosome='{chromosome}',
        genotype_hapmap=lambda wildcards: os.path.join(os.path.abspath(str(genotype_hapmap_folder)),str(wildcards.chromosome) + str(genotype_hapmap_file_extension)) if os.path.exists(os.path.join(os.path.abspath(str(genotype_hapmap_folder)),str(wildcards.chromosome) + str(genotype_hapmap_file_extension))) else 'NULL',
        genotype_data=lambda wildcards: os.path.join(os.path.abspath(str(genotype_data_folder)),str(wildcards.chromosome) + str(genotype_data_file_extension)) if os.path.exists(os.path.join(os.path.abspath(str(genotype_data_folder)),str(wildcards.chromosome) + str(genotype_data_file_extension))) else 'NULL',
        genotype_map=lambda wildcards: os.path.join(os.path.abspath(str(genotype_map_folder)),str(wildcards.chromosome) + str(genotype_map_file_extension)) if os.path.exists(os.path.join(os.path.abspath(str(genotype_map_folder)),str(wildcards.chromosome) + str(genotype_map_file_extension))) else 'NULL',
        kinship=lambda wildcards: os.path.join(os.path.abspath(str(kinship_folder)),str(wildcards.chromosome) + str(kinship_file_extension)) if os.path.exists(os.path.join(os.path.abspath(str(kinship_folder)),str(wildcards.chromosome) + str(kinship_file_extension))) else 'NULL',
        covariance_matrix=lambda wildcards: os.path.join(os.path.abspath(str(covariance_matrix_folder)),str(wildcards.chromosome) + str(covariance_matrix_file_extension)) if os.path.exists(os.path.join(os.path.abspath(str(covariance_matrix_folder)),str(wildcards.chromosome) + str(covariance_matrix_file_extension))) else 'NULL',
        snp_maf=snp_maf,
        model=model,
        pca_total=pca_total,
        out_folder=os.path.join(os.path.abspath(output_folder),'{chromosome}'),
        ulimit=ulimit
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'{chromosome}',"GAPIT_out",'{sample}.log')
    log:
        os.path.join(os.path.abspath(output_folder),'{chromosome}',"GAPIT_log",'{sample}.log')
    threads: threads
    resources:
        memory=memory
    shell:
        """
        mkdir -p {params.out_folder};
        Rscript {workflow_path}/scripts/R/GAPIT.R \
        -i {input.in_file} \
        -o {params.out_folder} \
        --genotype_hapmap {params.genotype_hapmap} \
        --genotype_data {params.genotype_data} \
        --genotype_map {params.genotype_map} \
        --kinship {params.kinship} \
        --covariance_matrix {params.covariance_matrix} \
        --snp_maf {params.snp_maf} \
        --model {params.model} \
        --pca_total {params.pca_total} | tee {log} {output.out_file};
        """
