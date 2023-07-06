#!/usr/bin/env python3

import sys
import os
import re
import json
import gzip
import pathlib
import shlex
import subprocess
import argparse

import pandas as pd


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    project_name = args.project_name
    workflow_path = args.workflow_path

    input_folder = args.input_folder
    output_folder = args.output_folder

    vcf_file = args.vcf_file

    gff_file = args.gff_file
    gff_category = args.gff_category
    gff_key = args.gff_key

    genotype_hapmap = args.genotype_hapmap
    genotype_data = args.genotype_data
    genotype_map = args.genotype_map

    kinship = args.kinship
    z_matrix = args.z_matrix
    corvariance_matrix = args.corvariance_matrix

    snp_maf = args.snp_maf
    model = args.model
    pca_total = args.pca_total

    memory = args.memory
    threads = args.threads

    jobs = args.jobs
    latency_wait = args.latency_wait
    cluster = args.cluster

    p_value_filter = args.p_value_filter
    fdr_corrected_p_value_filter = args.fdr_corrected_p_value_filter
    ld_length = args.ld_length

    #######################################################################
    # Check mandatory and at-least-one arguments
    #######################################################################
    # Project name checking
    if project_name is None:
        print('Project name is missing!!!')
        sys.exit(1)
    if project_name == '':
        print('Project name is missing!!!')
        sys.exit(1)

    # Workflow path checking
    if workflow_path is None:
        print('Workflow path is missing!!!')
        sys.exit(1)
    if workflow_path == '':
        print('Workflow path is missing!!!')
        sys.exit(1)
    if not workflow_path.exists():
        print('Workflow path does not exist!!!')
        sys.exit(1)

    # Input folder checking
    if input_folder is None:
        print('Input folder is missing!!!')
        sys.exit(1)
    if input_folder == '':
        print('Input folder is missing!!!')
        sys.exit(1)
    if not input_folder.exists():
        print('Input folder does not exist!!!')
        sys.exit(1)

    # VCF file checking
    if vcf_file is None:
        print('VCF file is missing!!!')
        sys.exit(1)
    if vcf_file == '':
        print('VCF file is missing!!!')
        sys.exit(1)
    if not vcf_file.exists():
        print('VCF file does not exist!!!')
        sys.exit(1)

    # GFF file checking
    if gff_file is None:
        print('GFF file is missing!!!')
        sys.exit(1)
    if gff_file == '':
        print('GFF file is missing!!!')
        sys.exit(1)
    if not gff_file.exists():
        print('GFF file does not exist!!!')
        sys.exit(1)

    # At-least-one arguments checking
    if genotype_hapmap != 'NULL':
        try:
            genotype_hapmap = pathlib.Path(genotype_hapmap)
            if not genotype_hapmap.exists():
                genotype_hapmap = 'NULL'
        except Exception as e:
            genotype_hapmap = 'NULL'
            print(e)
    if genotype_data != 'NULL':
        try:
            genotype_data = pathlib.Path(genotype_data)
            if not genotype_data.exists():
                genotype_data = 'NULL'
        except Exception as e:
            genotype_data = 'NULL'
            print(e)
    if genotype_map != 'NULL':
        try:
            genotype_map = pathlib.Path(genotype_map)
            if not genotype_map.exists():
                genotype_map = 'NULL'
        except Exception as e:
            genotype_map = 'NULL'
            print(e)
    if genotype_hapmap == 'NULL' and genotype_data == 'NULL' and genotype_map == 'NULL':
        print('Genotype hapmap, genotype data, and genotype map are not provided!!!')
        sys.exit(1)

    #######################################################################
    # Check if output parent folder exists
    # If not, create the output parent folder
    #######################################################################
    if not output_folder.exists():
        try:
            output_folder.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not output_folder.exists():
            print("Output folder does not exists!!!")
            sys.exit(1)

    # Create GAPIT output folder
    gapit_output_folder = output_folder.joinpath('GAPIT', model)
    if not gapit_output_folder.exists():
        try:
            gapit_output_folder.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not gapit_output_folder.exists():
            print("GAPIT output folder does not exists!!!")
            sys.exit(1)

    # Create VCFtools_Haploview output folder
    vcftools_haploview_output_folder = output_folder.joinpath('VCFtools_Haploview')
    if not vcftools_haploview_output_folder.exists():
        try:
            vcftools_haploview_output_folder.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not vcftools_haploview_output_folder.exists():
            print("VCFtools_Haploview output folder does not exists!!!")
            sys.exit(1)

    #######################################################################
    # Generate and save GAPIT configuration file
    #######################################################################
    gapit_config_data = {
        "project_name": project_name,
        "workflow_path": str(workflow_path),
        "input_folder": str(input_folder),
        "output_folder": str(gapit_output_folder),
        "genotype_hapmap": str(genotype_hapmap),
        "genotype_data": str(genotype_data),
        "genotype_map": str(genotype_map),
        "kinship": kinship,
        "z_matrix": z_matrix,
        "corvariance_matrix": corvariance_matrix,
        "snp_maf": snp_maf,
        "model": model,
        "pca_total": pca_total,
        "memory": memory,
        "threads": threads
    }

    with open(gapit_output_folder.joinpath("Snakemake_GAPIT_config.json"), 'w', encoding='utf-8') as writer:
        json.dump(gapit_config_data, writer, ensure_ascii=False, indent=4)

    #######################################################################
    # Run GAPIT in Snakemake
    #######################################################################

    # Construct snakemake command array
    snakemake_gapit_command_array = [
        "snakemake",
        '--jobs', str(jobs),
        '--latency-wait', str(latency_wait),
        '-p',
        '--snakefile', str(workflow_path.joinpath('rules', 'r_gapit.smk')),
        '--configfile', str(gapit_output_folder.joinpath('Snakemake_GAPIT_config.json'))
    ]

    if cluster is not None:
        if cluster != '':
            snakemake_gapit_command_array.append('--cluster')
            snakemake_gapit_command_array.append("\""+str(cluster)+"\"")

    # Construct snakemake command string from snakemake command array
    snakemake_gapit_command_string = str(' '.join(snakemake_gapit_command_array))

    print("\n\n{}\n".format(snakemake_gapit_command_string))

    # Split snakemake command string using shlex to prepare for subprocess
    snakemake_gapit_command = shlex.split(snakemake_gapit_command_string)

    try:
        # Run snakemake command using subprocess
        snakemake_gapit_outcome = subprocess.run(
            snakemake_gapit_command,
            capture_output=True,
            check=True,
            universal_newlines=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError("\n\nOutput: \n{} \n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(e.output, e.stdout, e.stderr, e.cmd, e.returncode))
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        with open(gapit_output_folder.joinpath("Snakemake_GAPIT.log"), 'w') as writer:
            writer.write("\n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(snakemake_gapit_outcome.stdout, snakemake_gapit_outcome.stderr, snakemake_gapit_outcome.args, snakemake_gapit_outcome.returncode))
    except Exception as e:
        print(e)
        sys.exit(1)

    #######################################################################
    # Combine GAPIT results
    #######################################################################

    #  Specify GAPIT results file and GAPIT result LD region file
    gapit_gwas_result_file = gapit_output_folder.joinpath("GAPIT.Association.GWAS_Results."+str(model)+".All.txt")
    gapit_gwas_result_ld_region_file = gapit_output_folder.joinpath("GAPIT.Association.GWAS_Results."+str(model)+".Unique_LD_Regions.txt")

    # Data save parameters
    header = True
    mode = 'w'
    # Find GAPIT results files
    if gapit_output_folder.exists():
        gapit_auto_output_folder = gapit_output_folder.joinpath('GAPIT_auto_output')
        if gapit_auto_output_folder.exists():
            gwas_result_files = list(gapit_auto_output_folder.glob('GAPIT.Association.GWAS_Results*.csv'))
            if len(gwas_result_files) > 0:
                for i in range(len(gwas_result_files)):
                    # Read GWAS result file
                    dat = pd.read_csv(filepath_or_buffer = gwas_result_files[i])
                    # Extract the trait from filename
                    trait = re.sub("(GAPIT.Association.GWAS_Results.)|(.csv)", "", str(gwas_result_files[i].name))
                    trait = trait.replace(model, '')
                    trait = trait.lstrip('\.')
                    # Add model, trait, ld_number, ld_start, and ld_end columns
                    dat['Model'] = model
                    dat['Trait'] = trait
                    dat["LD_length"] = ld_length
                    dat["LD_start"] = dat["Pos"] - ld_length
                    dat["LD_end"] = dat["Pos"] + ld_length
                    dat.loc[dat.LD_start < 0, 'LD_start'] = 0
                    # Filter data based on p-value and FDR corrected p-value
                    dat = dat[dat["P.value"] <= p_value_filter]
                    dat = dat[dat["H&B.P.Value"] <= fdr_corrected_p_value_filter]
                    # Write to one output file to append all data
                    dat.to_csv(path_or_buf = gapit_gwas_result_file, sep = '\t', mode = mode, header = header, index = False)
                    # Change parameters to make sure data append
                    if header:
                        header = False
                    if mode == 'w':
                        mode = 'a'
            else:
                print("No GAPIT GWAS result was created in the analysis!!!")
                sys.exit(1)
        else:
            print("No GAPIT auto output folder was created in the analysis!!!")
            sys.exit(1)
    else:
        print("Output folder does not exists!!!")
        sys.exit(1)


    #######################################################################
    # Query LD regions and drop duplicates
    #######################################################################

    # Read GAPIT GWAS result file
    dat = pd.read_table(filepath_or_buffer = gapit_gwas_result_file)
    # Query LD regions
    dat = dat[["Chr", "LD_start", "LD_end"]]
    # Sort by Chr, LD_start, and LD_end
    dat = dat.sort_values(by=['Chr', 'LD_start', 'LD_end'])
    # Drop duplicates
    dat = dat.drop_duplicates()
    # Save unique LD regions (header is not included in the output file so that snakemake part is easier)
    dat.to_csv(path_or_buf = gapit_gwas_result_ld_region_file, sep = '\t', mode = 'w', header = False, index = False)

    #######################################################################
    # Generate and save VCFtools_Haploview configuration file
    #######################################################################
    vcftools_haploview_config_data = {
        "project_name": project_name,
        "workflow_path": str(workflow_path),
        "input_file": str(gapit_gwas_result_ld_region_file),
        "output_folder": str(vcftools_haploview_output_folder),
        "vcf_file": str(vcf_file),
        "memory": memory,
        "threads": threads
    }

    with open(vcftools_haploview_output_folder.joinpath("Snakemake_VCFtools_Haploview_config.json"), 'w', encoding='utf-8') as writer:
        json.dump(vcftools_haploview_config_data, writer, ensure_ascii=False, indent=4)

    #######################################################################
    # Run VCFtools and Haploview in Snakemake
    #######################################################################

    # Construct snakemake command array
    snakemake_vcftools_haploview_command_array = [
        "snakemake",
        '--jobs', str(jobs),
        '--latency-wait', str(latency_wait),
        '-p',
        '--snakefile', str(workflow_path.joinpath('rules', 'vcftools_haploview.smk')),
        '--configfile', str(vcftools_haploview_output_folder.joinpath('Snakemake_VCFtools_Haploview_config.json'))
    ]

    if cluster is not None:
        if cluster != '':
            snakemake_vcftools_haploview_command_array.append('--cluster')
            snakemake_vcftools_haploview_command_array.append("\""+str(cluster)+"\"")

    # Construct snakemake command string from snakemake command array
    snakemake_vcftools_haploview_command_string = str(' '.join(snakemake_vcftools_haploview_command_array))

    print("\n\n{}\n".format(snakemake_vcftools_haploview_command_string))

    # Split snakemake command string using shlex to prepare for subprocess
    snakemake_vcftools_haploview_command = shlex.split(snakemake_vcftools_haploview_command_string)

    try:
        # Run snakemake command using subprocess
        snakemake_vcftools_haploview_outcome = subprocess.run(
            snakemake_vcftools_haploview_command,
            capture_output=True,
            check=True,
            universal_newlines=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError("\n\nOutput: \n{} \n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(e.output, e.stdout, e.stderr, e.cmd, e.returncode))
    except Exception as e:
        print(e)
        sys.exit(1)

    try:
        with open(vcftools_haploview_output_folder.joinpath("Snakemake_VCFtools_Haploview.log"), 'w') as writer:
            writer.write("\n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(snakemake_vcftools_haploview_outcome.stdout, snakemake_vcftools_haploview_outcome.stderr, snakemake_vcftools_haploview_outcome.args, snakemake_vcftools_haploview_outcome.returncode))
    except Exception as e:
        print(e)
        sys.exit(1)

    #######################################################################
    # Extract haploblock regions
    #######################################################################

    #Create haplotype blocks file and write header
    haplotype_blocks_file = vcftools_haploview_output_folder.joinpath("Haplotype_blocks.txt")
    with open(haplotype_blocks_file, 'w') as writer:
        writer.write("Chromosome\tStart\tEnd\tHaploblock_start\tHaploblock_end\n")

    # Find Haploview and VCFtools output folders
    haploview_output_folder = vcftools_haploview_output_folder.joinpath("Haploview")
    vcftools_output_folder = vcftools_haploview_output_folder.joinpath("VCFtools")

    # Find all GABRIELblocks files
    gabrielblocks_files = list(haploview_output_folder.glob("*.GABRIELblocks"))

    # Loop through all GABRIELblocks files, extract chromosome, ld_start, and ld_end, and find corresponding recode VCF file
    # Extract haploblock start and end based on indexes in GABRIELblocks file and positions in recode VCF file
    for i in range(len(gabrielblocks_files)):
        gabrielblocks_file = gabrielblocks_files[i]

        # Extract file prefix from GABRIELblocks file name
        file_prefix = re.sub('.GABRIELblocks', '', str(gabrielblocks_file.name))

        # Extract chromosome, start, and end from GABRIELblocks file prefix
        chrom_start_end = re.split('__', file_prefix)

        chromsome = chrom_start_end[0]
        ld_start = chrom_start_end[1]
        ld_end = chrom_start_end[2]

        # Create recode VCf file path
        recode_vcf_file = vcftools_output_folder.joinpath(file_prefix+".recode.vcf.gz")

        if recode_vcf_file.exists():
            # Create recode VCF chromosome and position arrays
            recode_vcf_chromosome_array = []
            recode_vcf_position_array = []

            # Read recode VCF file and collect chromosome and position data into arrays
            if str(recode_vcf_file).endswith('gz'):
                with gzip.open(recode_vcf_file, 'rt') as reader:
                    header = ""
                    while not header.strip().startswith("#CHROM"):
                        header = reader.readline()
                        header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
                    for line in reader:
                        line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                        recode_vcf_chromosome_array.append(line_array[0])
                        recode_vcf_position_array.append(line_array[1])
            else:
                with open(recode_vcf_file, "r") as reader:
                    header = ""
                    while not header.strip().startswith("#CHROM"):
                        header = reader.readline()
                        header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
                    for line in reader:
                        line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                        recode_vcf_chromosome_array.append(line_array[0])
                        recode_vcf_position_array.append(line_array[1])

        # Create haploblock start and end arrays
        haploblock_start_array = []
        haploblock_end_array = []

        # Read GABRIELblocks file and collect haploblock start and end data into arrays
        with open(gabrielblocks_file, "r") as reader:
            for line in reader:
                if re.search("MARKERS", line) is not None:
                    line = str(re.sub(".*MARKERS:", "", line)).strip("\n").strip("\r").strip("\r\n").strip(" ")
                    line_array = re.split(" ", line)
                    if len(line_array) > 1:
                        haploblock_start_index = int(line_array[0])-1
                        haploblock_end_index = int(line_array[(len(line_array)-1)])-1
                        if (haploblock_start_index > -1) and (haploblock_end_index > -1) and (haploblock_start_index < haploblock_end_index):
                            haploblock_start = recode_vcf_position_array[haploblock_start_index]
                            haploblock_end = recode_vcf_position_array[haploblock_end_index]
                            haploblock_start_array.append(haploblock_start)
                            haploblock_end_array.append(haploblock_end)

        # If haploblock start and end arrays are not empty, write haploblock data to file
        if (len(haploblock_start_array) > 0) and (len(haploblock_end_array) > 0) and (len(haploblock_start_array) == len(haploblock_end_array)):
            with open(haplotype_blocks_file, 'a') as writer:
                for j in range(len(haploblock_start_array)):
                    writer.write("{}\t{}\t{}\t{}\t{}\n".format(chromsome, ld_start, ld_end, haploblock_start_array[j], haploblock_end_array[j]))

    #######################################################################
    # Generate GWAS results file with haploblock data
    #######################################################################

    # Read data in haplotype blocks file into array
    haplotype_blocks_data_array = []
    with(open(haplotype_blocks_file, "r")) as reader:
        header = reader.readline()
        header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
        for line in reader:
            line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
            haplotype_blocks_data_array.append(line_array)

    # Create GWAS results file with haploblock data
    vcftools_haploview_gwas_result_file = vcftools_haploview_output_folder.joinpath("GAPIT.Association.GWAS_Results."+str(model)+".All.txt")

    f_hdl = open(vcftools_haploview_gwas_result_file, "w")

    # Read GAPIT GWAS results file and add haploblock data
    with(open(gapit_gwas_result_file, "r")) as reader:
        header = reader.readline()
        header = str(header).strip("\n").strip("\r").strip("\r\n") + "\tHaploblock_start\tHaploblock_end\n"
        f_hdl.write(header)
        for line in reader:
            line = str(line).strip("\n").strip("\r").strip("\r\n")
            line_array = line.split("\t")
            haploblock_start = ""
            haploblock_end = ""
            for i in range(len(haplotype_blocks_data_array)):
                # Match chromosome
                if (line_array[1] == haplotype_blocks_data_array[i][0]):
                    # Check if marker position is within haplotyle block region
                    if (int(haplotype_blocks_data_array[i][3]) <= int(line_array[2]) <= int(haplotype_blocks_data_array[i][4])):
                        haploblock_start = haplotype_blocks_data_array[i][3]
                        haploblock_end = haplotype_blocks_data_array[i][4]
                        break
            line = line + "\t" + str(haploblock_start) + "\t" + str(haploblock_end) + "\n"
            f_hdl.write(line)

    f_hdl.close()

    #######################################################################
    # Extract GFF regions that enclose haplotype blocks or markers
    #######################################################################

    # Load GFF file into data array
    gff_data_array = []
    with(open(gff_file, "r")) as reader:
        for line in reader:
            line = str(line).strip("\n").strip("\r").strip("\r\n")
            if not line.startswith("#"):
                line_array = line.split("\t")
                if line_array[2] == gff_category:
                    line_array[8] = str(re.sub(';.*', '', re.sub(fr'^.*?{gff_key}', '', line_array[8]))).strip("=").strip(":").strip(";").strip("_").strip("-")
                    gff_data_array.append(line_array)

    # Create GWAS results file with haploblock data and GFF regions
    output_gwas_result_file = output_folder.joinpath("GAPIT.Association.GWAS_Results."+str(model)+".All.txt")

    f_hdl = open(output_gwas_result_file, "w")

    # Read GWAS results file and add GFF regions
    with(open(vcftools_haploview_gwas_result_file, "r")) as reader:
        header = reader.readline()
        header = str(header).strip("\n").strip("\r").strip("\r\n") + "\tGene_start\tGene_end\tGene_ID\tOverlap_strategy\n"
        f_hdl.write(header)
        for line in reader:
            line = str(line).strip("\n").strip("\r").strip("\r\n")
            line_array = line.split("\t")
            gene_start = ""
            gene_end = ""
            gene_id = ""
            overlap_strategy = ""
            for i in range(len(gff_data_array)):
                # Match chromosome
                if (line_array[1] == gff_data_array[i][0]):
                    # Check if haplotype block region exists
                    if (line_array[13] != '') and (line_array[14] != ''):
                        if (int(gff_data_array[i][3]) <= int(line_array[13]) <= int(gff_data_array[i][4])) and (int(gff_data_array[i][3]) <= int(line_array[14]) <= int(gff_data_array[i][4])):
                            gene_start = gff_data_array[i][3]
                            gene_end = gff_data_array[i][4]
                            gene_id = gff_data_array[i][8]
                            overlap_strategy = "Full"
                        elif (overlap_strategy != "Full") and (int(line_array[13]) < int(gff_data_array[i][3])) and (int(gff_data_array[i][3]) <= int(line_array[14]) <= int(gff_data_array[i][4])):
                            gene_start = gff_data_array[i][3]
                            gene_end = gff_data_array[i][4]
                            gene_id = gff_data_array[i][8]
                            overlap_strategy = "Partial"
                        elif (overlap_strategy != "Full") and (int(gff_data_array[i][3]) <= int(line_array[13]) <= int(gff_data_array[i][4])) and  (int(line_array[14]) > int(gff_data_array[i][4])):
                            gene_start = gff_data_array[i][3]
                            gene_end = gff_data_array[i][4]
                            gene_id = gff_data_array[i][8]
                            overlap_strategy = "Partial"
                    # Check if marker position is within GFF region
                    if (overlap_strategy == "") and (int(gff_data_array[i][3]) <= int(line_array[2]) <= int(gff_data_array[i][4])):
                        gene_start = gff_data_array[i][3]
                        gene_end = gff_data_array[i][4]
                        gene_id = gff_data_array[i][8]
                        overlap_strategy = "Marker"

            line = line + "\t" + str(gene_start) + "\t" + str(gene_end) + "\t" + str(gene_id) + "\t" + str(overlap_strategy) + "\n"
            f_hdl.write(line)

    f_hdl.close()


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='HAPPI_GWAS_2', description='HAPPI_GWAS_2')

    parser.add_argument('-p', '--project_name', help='Project name', type=str, required=True)
    parser.add_argument('-w', '--workflow_path', help='Workflow path', type=pathlib.Path, required=True)

    parser.add_argument('-i', '--input_folder', help='Input folder', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder', type=pathlib.Path, required=True)

    parser.add_argument('-v', '--vcf_file', help='VCF file', type=pathlib.Path, required=True)

    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('--gff_category', help='GFF category', type=str, default='gene')
    parser.add_argument('--gff_key', help='GFF key', type=str, default='ID')

    parser.add_argument('--genotype_hapmap', default='NULL', help='Genotype hapmap', type=str)
    parser.add_argument('--genotype_data', default='NULL', help='Genotype data', type=str)
    parser.add_argument('--genotype_map', default='NULL', help='Genotype map', type=str)

    parser.add_argument('--kinship', default='NULL', type=str, help='Kinship matrix file')
    parser.add_argument('--z_matrix', default='NULL', type=str, help='Z matrix file')
    parser.add_argument('--corvariance_matrix', default='NULL', type=str, help='Corvariance matrix file')
    parser.add_argument('--snp_maf', default=0.0, type=float, help='SNP minor allele frequency')
    parser.add_argument('--model', default='MLM', type=str, help='Model')
    parser.add_argument('--pca_total', default=0, type=int, help='Total PCA')

    parser.add_argument('--memory', default=20, type=int, help='Memory')
    parser.add_argument('--threads', default=4, type=int, help='Threads')

    parser.add_argument('--jobs', default=2, type=int, help='Jobs')
    parser.add_argument('--latency_wait', default=60, type=int, help='Latency wait')
    parser.add_argument('--cluster', default='', type=str, help='Cluster parameters')

    parser.add_argument('--p_value_filter', default=1.0, type=float, help='P-value filter')
    parser.add_argument('--fdr_corrected_p_value_filter', default=1.0, type=float, help='FDR corrected p-value filter')
    parser.add_argument('--ld_length', default=10000, type=int, help='LD length')

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
