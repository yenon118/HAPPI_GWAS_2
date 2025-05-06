import sys
import os
import re


## Extract variables from the config file
project_name = config['project_name']
workflow_path = config['workflow_path']

input_folder = config['input_folder']
output_folder = config['output_folder']

feature_column_indexes = config['feature_column_indexes']

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

print("feature_column_indexes: ",feature_column_indexes)

print("ulimit: ",ulimit)
print("memory: ",memory)
print("threads: ",threads)

print("samples: ",samples)


## Collect all outputs
rule all:
    input:
        expand(os.path.join(os.path.abspath(output_folder),'{sample}_BLUP.txt'),sample=samples)


## Run BLUP in a snakemake rule
rule r_blup:
    input:
        in_file=os.path.join(os.path.abspath(input_folder),'{sample}.txt')
    params:
        feature_column_indexes=feature_column_indexes,
        out_folder=os.path.abspath(output_folder),
        ulimit=ulimit
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'{sample}_BLUP.txt')
    log:
        os.path.join(os.path.abspath(output_folder),'{sample}.log')
    threads: threads
    resources:
        memory=memory
    shell:
        """
        mkdir -p {params.out_folder};
        Rscript {workflow_path}/scripts/R/BLUP.R \
        -i {input.in_file} \
        -o {params.out_folder} \
        -e {params.feature_column_indexes} 2>&1 > {log};
        """
