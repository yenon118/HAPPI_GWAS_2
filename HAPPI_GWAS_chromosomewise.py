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

from statsmodels.stats.multitest import multipletests


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    project_name = args.project_name
    workflow_path = args.workflow_path

    input_folder = args.input_folder
    output_folder = args.output_folder

    chromosome_array = args.chromosome

    vcf_folder = args.vcf_folder
    vcf_file_extension = args.vcf_file_extension

    gff_file = args.gff_file
    gff_category = args.gff_category
    gff_key = args.gff_key

    genotype_hapmap_folder = args.genotype_hapmap_folder
    genotype_hapmap_file_extension = args.genotype_hapmap_file_extension
    genotype_data_folder = args.genotype_data_folder
    genotype_data_file_extension = args.genotype_data_file_extension
    genotype_map_folder = args.genotype_map_folder
    genotype_map_file_extension = args.genotype_map_file_extension

    kinship_folder = args.kinship_folder
    kinship_file_extension = args.kinship_file_extension
    covariance_matrix_folder = args.covariance_matrix_folder
    covariance_matrix_file_extension = args.covariance_matrix_file_extension

    snp_maf = args.snp_maf
    model = args.model
    pca_total = args.pca_total

    ulimit = args.ulimit
    memory = args.memory
    threads = args.threads

    keep_going = args.keep_going
    jobs = args.jobs
    latency_wait = args.latency_wait
    cluster = args.cluster

    p_value_filter = args.p_value_filter
    fdr_corrected_p_value_filter = args.fdr_corrected_p_value_filter
    multipletests_method = args.multipletests_method
    multipletests_p_value_filter = args.multipletests_p_value_filter
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

    # Chromosome array checking
    if chromosome_array is None:
        print('Chromosome array is missing!!!')
        sys.exit(1)
    if not isinstance(chromosome_array, list):
        print('Chromosome array is missing!!!')
        sys.exit(1)
    if len(chromosome_array) < 1:
        print('Chromosome array is missing!!!')
        sys.exit(1)

    # VCF folder checking
    if vcf_folder is None:
        print('VCF folder is missing!!!')
        sys.exit(1)
    if vcf_folder == '':
        print('VCF folder is missing!!!')
        sys.exit(1)
    if not vcf_folder.exists():
        print('VCF folder does not exist!!!')
        sys.exit(1)

    # VCF file extension checking
    if vcf_file_extension is None:
        print('VCF file extension is missing!!!')
        sys.exit(1)
    if vcf_file_extension == '':
        print('VCF file extension is missing!!!')
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
    if genotype_hapmap_folder:
        pass
    else:
        if genotype_data_folder and genotype_map_folder:
            pass
        else:
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
        "workflow_path": str(workflow_path) if workflow_path is not None else None,
        "input_folder": str(input_folder) if input_folder is not None else None,
        "output_folder": str(gapit_output_folder) if gapit_output_folder is not None else None,
        "chromosome_array": chromosome_array,
        "vcf_folder": str(vcf_folder) if vcf_folder is not None else None,
        "vcf_file_extension": vcf_file_extension,
        "genotype_hapmap_folder": str(genotype_hapmap_folder) if genotype_hapmap_folder is not None else None,
        "genotype_hapmap_file_extension": genotype_hapmap_file_extension,
        "genotype_data_folder": str(genotype_data_folder) if genotype_data_folder is not None else None,
        "genotype_data_file_extension": genotype_data_file_extension,
        "genotype_map_folder": str(genotype_map_folder) if genotype_map_folder is not None else None,
        "genotype_map_file_extension": genotype_map_file_extension,
        "kinship_folder": str(kinship_folder) if kinship_folder is not None else None,
        "kinship_file_extension": kinship_file_extension,
        "covariance_matrix_folder": str(covariance_matrix_folder) if covariance_matrix_folder is not None else None,
        "covariance_matrix_file_extension": covariance_matrix_file_extension,
        "snp_maf": snp_maf,
        "model": model,
        "pca_total": pca_total,
        "ulimit": ulimit,
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
        '--printshellcmds',
        '--snakefile', str(workflow_path.joinpath('rules', 'r_gapit_chromosomewise.smk')),
        '--configfile', str(gapit_output_folder.joinpath('Snakemake_GAPIT_config.json'))
    ]

    if keep_going is not None:
        if keep_going:
            snakemake_gapit_command_array.append('--keep-going')

    if cluster is not None:
        if cluster != '':
            snakemake_gapit_command_array.append('--executor cluster-generic --cluster-generic-submit-cmd')
            snakemake_gapit_command_array.append("\"" + str(cluster) + "\"")

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
        with open(gapit_output_folder.joinpath("Snakemake_GAPIT.log"), 'w') as writer:
            writer.write(
                "\n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(
                    snakemake_gapit_outcome.stdout,
                    snakemake_gapit_outcome.stderr,
                    snakemake_gapit_outcome.args,
                    snakemake_gapit_outcome.returncode
                )
            )
    except subprocess.CalledProcessError as e:
        print(
            "\n\nOutput: \n{} \n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(
                e.output, e.stdout, e.stderr, e.cmd, e.returncode
            )
        )
    except Exception as e:
        print(e)

    #######################################################################
    # Combine GAPIT results
    #######################################################################

    #  Specify GAPIT results file and GAPIT result LD region file
    gapit_gwas_result_file = gapit_output_folder.joinpath(
        "GAPIT.Association.GWAS_Results." + str(model) + ".All.txt"
    )
    gapit_gwas_result_ld_region_file = gapit_output_folder.joinpath(
        "GAPIT.Association.GWAS_Results." + str(model) + ".Unique_LD_Regions.txt"
    )

    # Data save parameters
    header = True
    mode = 'w'
    # Find GAPIT results files
    if gapit_output_folder.exists():
        gwas_result_files = list(
            gapit_output_folder.glob(os.path.join('*', 'GAPIT_auto_output', 'GAPIT.Association.GWAS_Results*.csv'))
        )
        gwas_result_files.sort()
        if len(gwas_result_files) > 0:
            for i in range(len(gwas_result_files)):
                # Read GWAS result file
                dat = pd.read_csv(filepath_or_buffer=gwas_result_files[i])
                # Extract the trait from filename
                trait = re.sub(
                    "(GAPIT.Association.GWAS_Results.)|(.csv)",
                    "",
                    str(gwas_result_files[i].name)
                )
                trait = trait.replace(model, '')
                trait = trait.lstrip('\\.')
                # Add model, trait, ld_number, ld_start, and ld_end columns
                dat['Model'] = model
                dat['Trait'] = trait
                dat["LD_length"] = ld_length
                dat["LD_start"] = dat["Pos"] - ld_length
                dat["LD_end"] = dat["Pos"] + ld_length
                dat.loc[dat.LD_start < 0, 'LD_start'] = 0
                # Check if data is not empty
                if dat.shape[0] > 0 and dat.shape[1] > 0:
                    # Write to one output file to append all data
                    dat.to_csv(path_or_buf=gapit_gwas_result_file, sep='\t', mode=mode, header=header, index=False)
                    # Change parameters to make sure data append
                    if header:
                        header = False
                    if mode == 'w':
                        mode = 'a'
        else:
            print("No GAPIT GWAS result was created in the analysis!!!")
            sys.exit(1)
    else:
        print("Output folder does not exists!!!")
        sys.exit(1)

    #######################################################################
    # Filter GWAS results
    #######################################################################

    # Data save parameters
    header = True
    mode = 'w'
    # Check if GAPIT GWAS result file exists
    if gapit_gwas_result_file.exists():
        # Read GAPIT GWAS result file
        dat = pd.read_table(filepath_or_buffer=gapit_gwas_result_file)
        # Check if data is not empty
        if dat.shape[0] > 0 and dat.shape[1] > 0:
            # Write multipletests_method into the data table
            dat["Multipletests_Method"] = multipletests_method
            # Compute for adjusted p-values
            dat['Multipletests_Adjusted_P_Value'] = dat.groupby(['Model', 'Trait'], group_keys=False)[
                ['P.value']].apply(
                lambda g: pd.Series(multipletests(g['P.value'], method=multipletests_method)[1], index=g.index),
                include_groups=False
            )
            # Filter data based on p-value and corrected p-value
            dat = dat[dat["P.value"] <= p_value_filter]
            dat = dat[dat["H&B.P.Value"] <= fdr_corrected_p_value_filter]
            dat = dat[dat["Multipletests_Adjusted_P_Value"] <= multipletests_p_value_filter]
            if dat.shape[0] > 0 and dat.shape[1] > 0:
                dat.to_csv(path_or_buf=gapit_gwas_result_file, sep='\t', mode=mode, header=header, index=False)
            else:
                print("GAPIT GWAS result is empty after filtering!!!")
                sys.exit(1)
    else:
        print("GAPIT GWAS result file does not exists!!!")
        sys.exit(1)

    #######################################################################
    # Query LD regions and drop duplicates
    #######################################################################

    # Check if GAPIT GWAS result file exists
    if gapit_gwas_result_file.exists():
        # Read GAPIT GWAS result file
        dat = pd.read_table(filepath_or_buffer=gapit_gwas_result_file)
        # Check if data is not empty
        if dat.shape[0] > 0 and dat.shape[1] > 0:
            # Query LD regions
            dat = dat[["Chr", "LD_start", "LD_end"]]
            # Sort by Chr, LD_start, and LD_end
            dat = dat.sort_values(by=['Chr', 'LD_start', 'LD_end'])
            # Drop duplicates
            dat = dat.drop_duplicates()
            # Check if data is not empty
            if dat.shape[0] > 0 and dat.shape[1] > 0:
                # Save unique LD regions (header is not included in the output file so that snakemake part is easier)
                dat.to_csv(path_or_buf=gapit_gwas_result_ld_region_file, sep='\t', mode='w', header=False, index=False)
        else:
            print("No GAPIT GWAS result was created in the analysis!!!")
            sys.exit(1)
    else:
        print("GAPIT GWAS result file does not exists!!!")
        sys.exit(1)

    if not gapit_gwas_result_ld_region_file.exists():
        print("GAPIT GWAS result LD region file cannot be created as no LD region found!!!")
        sys.exit(1)

    #######################################################################
    # Generate and save VCFtools_Haploview configuration file
    #######################################################################
    vcftools_haploview_config_data = {
        "project_name": project_name,
        "workflow_path": str(workflow_path) if workflow_path is not None else None,
        "input_file": str(gapit_gwas_result_ld_region_file) if gapit_gwas_result_ld_region_file is not None else None,
        "output_folder": str(
            vcftools_haploview_output_folder
        ) if vcftools_haploview_output_folder is not None else None,
        "vcf_folder": str(vcf_folder) if vcf_folder is not None else None,
        "vcf_file_extension": vcf_file_extension,
        "ulimit": ulimit,
        "memory": memory,
        "threads": threads
    }

    with open(
            vcftools_haploview_output_folder.joinpath("Snakemake_VCFtools_Haploview_config.json"),
            'w',
            encoding='utf-8'
    ) as writer:
        json.dump(vcftools_haploview_config_data, writer, ensure_ascii=False, indent=4)

    #######################################################################
    # Run VCFtools and Haploview in Snakemake
    #######################################################################

    # Construct snakemake command array
    snakemake_vcftools_haploview_command_array = [
        "snakemake",
        '--jobs', str(jobs),
        '--latency-wait', str(latency_wait),
        '--printshellcmds',
        '--snakefile', str(workflow_path.joinpath('rules', 'vcftools_haploview_chromosomewise.smk')),
        '--configfile', str(vcftools_haploview_output_folder.joinpath('Snakemake_VCFtools_Haploview_config.json'))
    ]

    if keep_going is not None:
        if keep_going:
            snakemake_vcftools_haploview_command_array.append('--keep-going')

    if cluster is not None:
        if cluster != '':
            snakemake_vcftools_haploview_command_array.append('--executor cluster-generic --cluster-generic-submit-cmd')
            snakemake_vcftools_haploview_command_array.append("\"" + str(cluster) + "\"")

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
        with open(vcftools_haploview_output_folder.joinpath("Snakemake_VCFtools_Haploview.log"), 'w') as writer:
            writer.write(
                "\n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(
                    snakemake_vcftools_haploview_outcome.stdout,
                    snakemake_vcftools_haploview_outcome.stderr,
                    snakemake_vcftools_haploview_outcome.args,
                    snakemake_vcftools_haploview_outcome.returncode
                )
            )
    except subprocess.CalledProcessError as e:
        print(
            "\n\nOutput: \n{} \n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(
                e.output, e.stdout, e.stderr, e.cmd, e.returncode
            )
        )
    except Exception as e:
        print(e)

    #######################################################################
    # Extract haploblock regions
    #######################################################################

    # Create haplotype blocks file and write header
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
    if len(gabrielblocks_files) > 0:
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
            recode_vcf_file = vcftools_output_folder.joinpath(file_prefix + ".recode.vcf.gz")

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
                            haploblock_start_index = int(float(line_array[0])) - 1
                            haploblock_end_index = int(float(line_array[(len(line_array) - 1)])) - 1
                            if (haploblock_start_index > -1) and (haploblock_end_index > -1) and (
                                    haploblock_start_index < haploblock_end_index):
                                haploblock_start = recode_vcf_position_array[haploblock_start_index]
                                haploblock_end = recode_vcf_position_array[haploblock_end_index]
                                haploblock_start_array.append(haploblock_start)
                                haploblock_end_array.append(haploblock_end)

            # If haploblock start and end arrays are not empty, write haploblock data to file
            if (len(haploblock_start_array) > 0) and (len(haploblock_end_array) > 0) and (
                    len(haploblock_start_array) == len(haploblock_end_array)):
                with open(haplotype_blocks_file, 'a') as writer:
                    for j in range(len(haploblock_start_array)):
                        writer.write(
                            "{}\t{}\t{}\t{}\t{}\n".format(
                                chromsome,
                                ld_start,
                                ld_end,
                                haploblock_start_array[j],
                                haploblock_end_array[j]
                            )
                        )

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
    vcftools_haploview_gwas_result_file = vcftools_haploview_output_folder.joinpath(
        "GAPIT.Association.GWAS_Results." + str(model) + ".All.txt"
    )

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
            if len(haplotype_blocks_data_array) > 0:
                for i in range(len(haplotype_blocks_data_array)):
                    # Match chromosome
                    if line_array[1] == haplotype_blocks_data_array[i][0]:
                        # Check if marker position is within haplotyle block region
                        if int(float(haplotype_blocks_data_array[i][3])) <= int(float(line_array[2])) <= int(
                                float(haplotype_blocks_data_array[i][4])):
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
                    line_array[8] = str(
                        re.sub(';.*', '', re.sub(fr'^.*?{gff_key}', '', line_array[8]))
                    ).strip("=").strip(":").strip(";").strip("_").strip("-")
                    gff_data_array.append(line_array)

    # Create GWAS results file with haploblock data and GFF regions
    output_gwas_result_file = output_folder.joinpath("GAPIT.Association.GWAS_Results." + str(model) + ".All.txt")

    f_hdl = open(output_gwas_result_file, "w")

    # Read GWAS results file and add GFF regions
    with(open(vcftools_haploview_gwas_result_file, "r")) as reader:
        header = reader.readline()
        header = str(header).strip("\n").strip("\r").strip(
            "\r\n") + "\tGene_start\tGene_end\tGene_ID\tOverlap_strategy\n"
        f_hdl.write(header)
        for line in reader:
            line = str(line).strip("\n").strip("\r").strip("\r\n")
            line_array = line.split("\t")
            gene_start = ""
            gene_end = ""
            gene_id = ""
            overlap_strategy = ""
            if len(gff_data_array) > 0:
                for i in range(len(gff_data_array)):
                    # Match chromosome
                    if line_array[1] == gff_data_array[i][0]:
                        # Check if haplotype block region exists
                        if (line_array[15] != '') and (line_array[16] != ''):
                            if (int(float(gff_data_array[i][3])) <= int(float(line_array[15])) <= int(
                                    float(gff_data_array[i][4]))) and (
                                    int(float(gff_data_array[i][3])) <= int(float(line_array[16])) <= int(
                                float(gff_data_array[i][4]))):
                                gene_start = gff_data_array[i][3]
                                gene_end = gff_data_array[i][4]
                                gene_id = gff_data_array[i][8]
                                overlap_strategy = "Full"
                            if (overlap_strategy != "Full") and (
                                    int(float(line_array[15])) < int(float(gff_data_array[i][3]))) and (
                                    int(float(gff_data_array[i][3])) <= int(float(line_array[16])) <= int(
                                float(gff_data_array[i][4]))):
                                gene_start = gff_data_array[i][3]
                                gene_end = gff_data_array[i][4]
                                gene_id = gff_data_array[i][8]
                                overlap_strategy = "Partial"
                            if (overlap_strategy != "Full") and (
                                    int(float(gff_data_array[i][3])) <= int(float(line_array[15])) <= int(
                                float(gff_data_array[i][4]))) and (
                                    int(float(line_array[16])) > int(float(gff_data_array[i][4]))):
                                gene_start = gff_data_array[i][3]
                                gene_end = gff_data_array[i][4]
                                gene_id = gff_data_array[i][8]
                                overlap_strategy = "Partial"
                        # Check if marker position is within GFF region
                        if (overlap_strategy == "") and (
                                int(float(gff_data_array[i][3])) <= int(float(line_array[2])) <= int(
                            float(gff_data_array[i][4]))):
                            gene_start = gff_data_array[i][3]
                            gene_end = gff_data_array[i][4]
                            gene_id = gff_data_array[i][8]
                            overlap_strategy = "Marker"

            line = line + "\t" + str(gene_start) + "\t" + str(gene_end) + "\t" + str(gene_id) + "\t" + str(
                overlap_strategy) + "\n"
            f_hdl.write(line)

    f_hdl.close()


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='HAPPI_GWAS', description='HAPPI_GWAS')

    parser.add_argument('-p', '--project_name', help='Project name', type=str, required=True)
    parser.add_argument('-w', '--workflow_path', help='Workflow path', type=pathlib.Path, required=True)

    parser.add_argument('-i', '--input_folder', help='Input folder', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_folder', help='Output folder', type=pathlib.Path, required=True)

    parser.add_argument('-c', '--chromosome', help='Chromosome', type=str, action='append', required=True)

    parser.add_argument('-v', '--vcf_folder', help='VCF folder', type=pathlib.Path, required=True)
    parser.add_argument('-x', '--vcf_file_extension', help='VCF file extension', type=str, required=True)

    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('--gff_category', help='GFF category', type=str, default='gene')
    parser.add_argument('--gff_key', help='GFF key', type=str, default='ID')

    parser.add_argument(
        '--genotype_hapmap_folder',
        help='Genotype hapmap folder',
        type=pathlib.Path
    )
    parser.add_argument(
        '--genotype_hapmap_file_extension',
        help='Genotype hapmap file extension',
        type=str
    )
    parser.add_argument(
        '--genotype_data_folder',
        help='Genotype data folder',
        type=pathlib.Path
    )
    parser.add_argument(
        '--genotype_data_file_extension',
        help='Genotype data file extension',
        type=str
    )
    parser.add_argument(
        '--genotype_map_folder',
        help='Genotype map folder',
        type=pathlib.Path
    )
    parser.add_argument(
        '--genotype_map_file_extension',
        help='Genotype map file extension',
        type=str
    )

    parser.add_argument(
        '--kinship_folder',
        type=pathlib.Path,
        help='Kinship matrix folder'
    )
    parser.add_argument(
        '--kinship_file_extension',
        type=str,
        help='Kinship matrix file extension'
    )
    parser.add_argument(
        '--covariance_matrix_folder',
        type=pathlib.Path,
        help='Covariance matrix folder'
    )
    parser.add_argument(
        '--covariance_matrix_file_extension',
        type=str,
        help='Covariance matrix file extension'
    )
    parser.add_argument('--snp_maf', default=0.0, type=float, help='SNP minor allele frequency')
    parser.add_argument('--model', default='MLM', type=str, help='Model')
    parser.add_argument('--pca_total', default=0, type=int, help='Total PCA')

    parser.add_argument('--ulimit', default=0, type=int, help='Ulimit')
    parser.add_argument('--memory', default=20, type=int, help='Memory')
    parser.add_argument('--threads', default=4, type=int, help='Threads')

    parser.add_argument('--keep_going', action='store_true', help='Keep going')
    parser.add_argument('--jobs', default=2, type=int, help='Jobs')
    parser.add_argument('--latency_wait', default=180, type=int, help='Latency wait')
    parser.add_argument('--cluster', default='', type=str, help='Cluster parameters')

    parser.add_argument('--p_value_filter', default=1.0, type=float, help='P-value filter')
    parser.add_argument(
        '--fdr_corrected_p_value_filter',
        default=1.0,
        type=float,
        help='FDR corrected p-value filter'
    )
    parser.add_argument(
        '--multipletests_method',
        default='fdr_bh',
        type=str,
        help='Multipletests method'
    )
    parser.add_argument(
        '--multipletests_p_value_filter',
        default=1.0,
        type=float,
        help='Multipletests corrected p-value filter'
    )
    parser.add_argument('--ld_length', default=10000, type=int, help='LD length')

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
