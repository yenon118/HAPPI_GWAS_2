#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

library(GAPIT)


##################################################
# Constants/Variables
##################################################


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

# Mandatory arguments
parser$add_argument("-i", "--input_file", type="character", help="Input file")
parser$add_argument("-o", "--output_path", type="character", help="Output folder")

# At-least-one arguments
parser$add_argument("--genotype_hapmap", type="character", default=NULL, help="Genotype hapmap file")
parser$add_argument("--genotype_data", type="character", default=NULL, help="Genotype data file")
parser$add_argument("--genotype_map", type="character", default=NULL, help="Genotype map file")

# Optional arguments
parser$add_argument("--kinship", type="character", default=NULL, help="Kinship matrix file")
parser$add_argument("--covariance_matrix", type="character", default=NULL, help="Covariance matrix file")

parser$add_argument("--snp_maf", type="numeric", default=0.0, help="SNP minor allele frequency")
parser$add_argument("--model", type="character", default="MLM", help="Model")
parser$add_argument("--pca_total", type="integer", default=0, help="Total PCA")

# Parse arguments
args <- parser$parse_args()


# Check input file
input_file <- tryCatch({file.path(args$input_file)}, error = function(e) {return(NULL)})

if (is.null(input_file) | identical(input_file, character(0))) {
	print("Invalid input file!!!")
	quit(status=0)
}

# Check output path
output_path <- tryCatch({file.path(args$output_path)}, error = function(e) {return(NULL)})

if (is.null(output_path) | identical(output_path, character(0))) {
	print("Invalid output path!!!")
	quit(status=0)
}

# At-least-one argument assignment
G <- args$genotype_hapmap
GD <- args$genotype_data
GM <- args$genotype_map

if (is.null(G) & is.null(GD) & is.null(GM)) {
	print("Genotype hapmap, genotype data, and genotype map do not exists!!!")
	quit(status=0)
}

# Optional argument assignment
KI <- args$kinship
CV <- args$covariance_matrix
SNP.MAF <- args$snp_maf
model <- args$model
PCA.total <- args$pca_total


##################################################
# Read in input file
##################################################

if (!is.null(input_file)) {
	if (input_file != "NULL") {
		if (endsWith(input_file, "csv")) {
			Y <- tryCatch({
				read.csv(
					file = input_file,
					header = TRUE,
					check.names = FALSE,
					stringsAsFactors = FALSE
				)
			}, error = function(e) {
				print("Input file cannot be read!!!")
				quit(status=0)
			})
		} else {
			Y <- tryCatch({
				read.table(
					file = input_file,
					sep = "\t",
					header = TRUE,
					check.names = FALSE,
					stringsAsFactors = FALSE
				)
			}, error = function(e) {
				print("Input file cannot be read!!!")
				quit(status=0)
			})
		}
	} else {
		print("Input file (Y) is set to NULL!!!")
		quit(status=0)
	}
	if (!is.null(Y)) {
		if (nrow(Y) > 0 & ncol(Y) > 0){
			print("Input file (Y):")
			print(head(Y))
			print(dim(Y))
		}
	} else {
		print("Input file (Y) is set to NULL!!!")
		quit(status=0)
	}
} else {
	print("Input file (Y) is set to NULL!!!")
	quit(status=0)
}


##################################################
# Output folder
##################################################

# Construct output folder
output_path <- file.path(output_path)
gapit_auto_output_path <- file.path(output_path, "GAPIT_auto_output")

# Check and create output folder
if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		print("Output path cannot be created!!!")
		quit(status=0)
	}
}

# Check and create GAPIT auto output folder
if(!dir.exists(gapit_auto_output_path)){
	dir.create(gapit_auto_output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(gapit_auto_output_path)){
		print("Output path cannot be created!!!")
		quit(status=0)
	}
}


##################################################
# Read in at-least-one file
##################################################

if (!is.null(G)) {
	if (G != "NULL") {
		G <- tryCatch({
			read.table(
				file = G,
				sep = "\t",
				header = FALSE,
				comment.char = "",
				check.names = FALSE,
				stringsAsFactors = FALSE
			)
		}, error = function(e) {
			return(NULL)
		})
	} else {
		G <- NULL
	}
	if (!is.null(G)) {
		if (nrow(G) > 0 & ncol(G) > 0){
			print("Genotype hapmap (G):")
			print(head(G))
			print(dim(G))
		}
	} else {
		print("Genotype hapmap (G) is set to NULL!!!")
	}
} else {
	print("Genotype hapmap (G) is set to NULL!!!")
}

if (!is.null(GD)) {
	if (GD != "NULL") {
		GD <- tryCatch({
			read.table(
				file = GD,
				sep = "\t",
				header = TRUE,
				check.names = FALSE,
				stringsAsFactors = FALSE
			)
		}, error = function(e) {
			return(NULL)
		})
	} else {
		GD <- NULL
	}
	if (!is.null(GD)) {
		if (nrow(GD) > 0 & ncol(GD) > 0){
			print("Genotype data (GD):")
			print(head(GD))
			print(dim(GD))
		}
	} else {
		print("Genotype data (GD) is set to NULL!!!")
	}
} else {
	print("Genotype data (GD) is set to NULL!!!")
}

if (!is.null(GM)) {
	if (GM != "NULL") {
		GM <- tryCatch({
			read.table(
				file = GM,
				sep = "\t",
				header = TRUE,
				check.names = FALSE,
				stringsAsFactors = FALSE
			)
		}, error = function(e) {
			return(NULL)
		})
	} else {
		GM <- NULL
	}
	if (!is.null(GM)) {
		if (nrow(GM) > 0 & ncol(GM) > 0){
			print("Genotype map (GM):")
			print(head(GM))
			print(dim(GM))
		}
	} else {
		print("Genotype map (GM) is set to NULL!!!")
	}
} else {
	print("Genotype map (GM) is set to NULL!!!")
}

# Check three of the at-least-one file
if (is.null(G) & is.null(GD) & is.null(GM)) {
	print("Genotype hapmap, genotype data, and genotype map cannot be read!!!")
	quit(status=0)
}


##################################################
# Read in optional file
##################################################

if (!is.null(KI)) {
	if (KI != "NULL") {
		KI <- tryCatch({
			read.table(
				file = KI,
				sep = "\t",
				header = TRUE,
				check.names = FALSE,
				stringsAsFactors = FALSE
			)
		}, error = function(e) {
			return(NULL)
		})
	} else {
		KI <- NULL
	}
	if (!is.null(KI)) {
		if (nrow(KI) > 0 & ncol(KI) > 0){
			print("Kinship matrix (KI):")
			print(head(KI))
			print(dim(KI))
		}
	} else {
		print("Kinship matrix (KI) is set to NULL!!!")
	}
} else {
	print("Kinship matrix (KI) is set to NULL!!!")
}

if (!is.null(CV)) {
	if (CV != "NULL") {
		CV <- tryCatch({
			read.table(
				file = CV,
				sep = "\t",
				header = TRUE,
				row.names = 1,
				check.names = FALSE,
				stringsAsFactors = FALSE
			)
		}, error = function(e) {
			return(NULL)
		})
	} else {
		CV <- NULL
	}
	if (!is.null(CV)) {
		if (nrow(CV) > 0 & ncol(CV) > 0){
			print("Covariance matrix (CV):")
			print(head(CV))
			print(dim(CV))
		}
	} else {
		print("Covariance matrix (CV) is set to NULL!!!")
	}
} else {
	print("Covariance matrix (CV) is set to NULL!!!")
}


##################################################
# Print oprional arguments
##################################################

print(paste0("SNP.MAF is set to ", SNP.MAF, "."))
print(paste0("model is set to ", model, "."))
print(paste0("PCA.total is set to ", PCA.total, "."))


##################################################
# Run GAPIT
##################################################

current_directory <- getwd()

setwd(gapit_auto_output_path)


print("---------------------------------- GAPIT starts ----------------------------------")

gapit_result <- GAPIT(
	Y = Y,
	G = G,
	GD = GD,
	GM = GM,
	KI = KI,
	CV = CV,
	SNP.MAF = SNP.MAF,
	model = model,
	PCA.total = PCA.total
)

print("---------------------------------- GAPIT ends ----------------------------------")


setwd(current_directory)
