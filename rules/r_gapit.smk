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
covariance_matrix = config['covariance_matrix']

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

print("genotype_hapmap: ",genotype_hapmap)
print("genotype_data: ",genotype_data)
print("genotype_map: ",genotype_map)

print("kinship: ",kinship)
print("covariance_matrix: ",covariance_matrix)

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
        expand(os.path.join(os.path.abspath(output_folder),"GAPIT_out",'{sample}.log'),sample=samples)


## Run GAPIT in a snakemake rule
rule r_gapit:
    input:
        in_file=os.path.join(os.path.abspath(input_folder),'{sample}.txt')
    params:
        sample='{sample}',
        genotype_hapmap=os.path.abspath(str(genotype_hapmap)) if os.path.exists(os.path.abspath(str(genotype_hapmap))) else 'NULL',
        genotype_data=os.path.abspath(str(genotype_data)) if os.path.exists(os.path.abspath(str(genotype_data))) else 'NULL',
        genotype_map=os.path.abspath(str(genotype_map)) if os.path.exists(os.path.abspath(str(genotype_map))) else 'NULL',
        kinship=os.path.abspath(str(kinship)) if os.path.exists(os.path.abspath(str(kinship))) else 'NULL',
        covariance_matrix=os.path.abspath(str(covariance_matrix)) if os.path.exists(os.path.abspath(str(covariance_matrix))) else 'NULL',
        snp_maf=snp_maf,
        model=model,
        pca_total=pca_total,
        out_folder=os.path.abspath(output_folder),
        ulimit=ulimit
    output:
        out_file=os.path.join(os.path.abspath(output_folder),"GAPIT_out",'{sample}.log')
    log:
        os.path.join(os.path.abspath(output_folder),"GAPIT_log",'{sample}.log')
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
