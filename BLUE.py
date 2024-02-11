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

	feature_column_indexes = args.feature_column_indexes

	ulimit = args.ulimit
	memory = args.memory
	threads = args.threads

	keep_going = args.keep_going
	jobs = args.jobs
	latency_wait = args.latency_wait
	cluster = args.cluster

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

	#######################################################################
	# Generate and save BLUE configuration file
	#######################################################################
	blue_config_data = {
		"project_name": str(project_name),
		"workflow_path": str(workflow_path),
		"input_folder": str(input_folder),
		"output_folder": str(output_folder),
		"feature_column_indexes": str(feature_column_indexes),
		"ulimit": ulimit,
		"memory": memory,
		"threads": threads
	}

	with open(output_folder.joinpath("Snakemake_BLUE_config.json"), 'w', encoding='utf-8') as writer:
		json.dump(blue_config_data, writer, ensure_ascii=False, indent=4)

	#######################################################################
	# Run BLUE in Snakemake
	#######################################################################

	# Construct snakemake command array
	snakemake_blue_command_array = [
		"snakemake",
		'--jobs', str(jobs),
		'--latency-wait', str(latency_wait),
		'--printshellcmds',
		'--snakefile', str(workflow_path.joinpath('rules', 'r_blue.smk')),
		'--configfile', str(output_folder.joinpath('Snakemake_BLUE_config.json'))
	]

	if keep_going is not None:
		if keep_going:
			snakemake_blue_command_array.append('--keep-going')

	if cluster is not None:
		if cluster != '':
			snakemake_blue_command_array.append('--executor cluster-generic --cluster-generic-submit-cmd')
			snakemake_blue_command_array.append("\""+str(cluster)+"\"")

	# Construct snakemake command string from snakemake command array
	snakemake_blue_command_string = str(' '.join(snakemake_blue_command_array))

	print("\n\n{}\n".format(snakemake_blue_command_string))

	# Split snakemake command string using shlex to prepare for subprocess
	snakemake_blue_command = shlex.split(snakemake_blue_command_string)

	try:
		# Run snakemake command using subprocess
		snakemake_blue_outcome = subprocess.run(
			snakemake_blue_command,
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
		with open(output_folder.joinpath("Snakemake_BLUE.log"), 'w') as writer:
			writer.write("\n\nStandard output: \n{} \n\nStandard error: \n{} \n\nCommand: \n{} \n\nReturn code: \n{} \n\n".format(snakemake_blue_outcome.stdout, snakemake_blue_outcome.stderr, snakemake_blue_outcome.args, snakemake_blue_outcome.returncode))
	except Exception as e:
		print(e)
		sys.exit(1)


if __name__ == "__main__":
	#######################################################################
	# Parse arguments
	#######################################################################
	parser = argparse.ArgumentParser(prog='BLUE', description='BLUE')

	parser.add_argument('-p', '--project_name', help='Project name', type=str, required=True)
	parser.add_argument('-w', '--workflow_path', help='Workflow path', type=pathlib.Path, required=True)

	parser.add_argument('-i', '--input_folder', help='Input folder', type=pathlib.Path, required=True)
	parser.add_argument('-o', '--output_folder', help='Output folder', type=pathlib.Path, required=True)

	parser.add_argument('-e', '--feature_column_indexes', help='Feature column indexes', type=str, default='1,2')

	parser.add_argument('--ulimit', default=10240, type=int, help='Ulimit')
	parser.add_argument('--memory', default=20, type=int, help='Memory')
	parser.add_argument('--threads', default=4, type=int, help='Threads')

	parser.add_argument('--keep_going', action='store_true', help='Keep going')
	parser.add_argument('--jobs', default=2, type=int, help='Jobs')
	parser.add_argument('--latency_wait', default=180, type=int, help='Latency wait')
	parser.add_argument('--cluster', default='', type=str, help='Cluster parameters')

	args = parser.parse_args()

	#######################################################################
	# Call main function
	#######################################################################
	main(args)
