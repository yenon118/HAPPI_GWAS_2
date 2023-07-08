#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)


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
parser$add_argument("-e", "--feature_column_indexes", type="character", default="1,2", help="Feature column indexes")

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


feature_column_indexes <- args$feature_column_indexes


##################################################
# Read in input file
##################################################

if (!is.null(input_file)) {
	if (input_file != "NULL") {
		if (endsWith(input_file, "csv")) {
			dat <- tryCatch({
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
			dat <- tryCatch({
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
		print("Input file is set to NULL!!!")
		quit(status=0)
	}
	if (!is.null(dat)) {
		if (nrow(dat) > 0 & ncol(dat) > 0){
			print("Input file:")
			print(head(dat))
			print(dim(dat))
		}
	} else {
		print("Input file is set to NULL!!!")
		quit(status=0)
	}
} else {
	print("Input file is set to NULL!!!")
	quit(status=0)
}


##################################################
# Output folder
##################################################

# Construct output folder
output_path <- file.path(output_path)

# Check and create output folder
if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		print("Output path cannot be created!!!")
		quit(status=0)
	}
}


##################################################
# Extract feature variables and target variables
##################################################

feature_column_indexes <- strsplit(feature_column_indexes, "[,;]")[[1]]
feature_column_indexes <- as.numeric(feature_column_indexes)
feature_column_indexes <- feature_column_indexes[!is.null(feature_column_indexes)]
feature_column_indexes <- feature_column_indexes[!is.na(feature_column_indexes)]
feature_column_indexes <- feature_column_indexes[is.finite(feature_column_indexes)]

target_column_indexes <- 1:ncol(dat)
target_column_indexes <- target_column_indexes[!(target_column_indexes %in% feature_column_indexes)]


##################################################
# Outlier removal
##################################################

# Create lmer formula
if (length(feature_column_indexes) > 0) {
	termlabels <- rep("", length(feature_column_indexes))
	if (length(feature_column_indexes) > 0) {
		termlabels[1] <- colnames(dat)[feature_column_indexes[1]]
		if (length(feature_column_indexes) > 1) {
			for (i in 2:length(feature_column_indexes)) {
				termlabels[i] <- paste0("(1|", colnames(dat)[feature_column_indexes[i]], ")")
			}
		}
	}
	termlabels <- termlabels[termlabels != ""]
} else {
	print("No feature column is specified!!!")
	quit(status=0)
}

# Calculate threshold
threshold <- qt(1-.05/(2*nrow(dat)), (nrow(dat)-3))

# Find outlier and generate outlier removed dataset and outlier dataset
outlier_removed_dat <- dat
outlier_dat <- dat
for (i in target_column_indexes){
	# Clean miscellaneous data
	dat[is.infinite(dat[,i]),i] = NA
	dat[is.nan(dat[,i]),i] = NA
	dat[is.null(dat[,i]),i] = NA
	dat <- dat[rowSums(is.na(dat)) != ncol(dat),]
	dat <- dat[,colSums(is.na(dat)) != nrow(dat)]

	# Fit lmer model
	lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

	res <- residuals(lme)
	H <- hatvalues(lme)
	sigma <- summary(lme)$sigm
	sres <- as.numeric(res/(sigma*sqrt(1-H)))

	# Generate outlier flags
	flags <- (abs(sres) >= threshold)
	flags[is.na(flags)] <- FALSE

	# Generate outlier removed dataset and outlier dataset
	outlier_removed_dat[flags,i] <- NA
	outlier_dat[!flags,i] <- NA
}

# Update data
dat <- outlier_removed_dat


##################################################
# Box-cox Transformation
##################################################

# Run transformation for each trait
transformed_dat <- dat
lambda_dat <- data.frame(Trait = colnames(dat)[target_column_indexes], Lambda = rep(NA, length(target_column_indexes)))
for(i in target_column_indexes){
	# Clean miscellaneous data
	dat[is.infinite(dat[,i]),i] = NA
	dat[is.nan(dat[,i]),i] = NA
	dat[is.null(dat[,i]),i] = NA
	dat <- dat[rowSums(is.na(dat)) != ncol(dat),]
	dat <- dat[,colSums(is.na(dat)) != nrow(dat)]

	# Create lmer model
	lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

	# Run power transform
	transformed_out <- tryCatch({
		car::powerTransform(lme, family="bcPower", lambda=c(-2, 2))
	}, error = function (e){
		cat(rep("\n", 2))
		print(paste0("Lambda cannot be calculated for ", colnames(dat)[i], "!!!"))
		return(1)
	})

	# Extract lambda
	if (is.list(transformed_out)) {
		lambda <- tryCatch({
			transformed_out$lambda
		},error = function (e) {
			cat(rep("\n", 2))
			print(paste0("Fail to extract lambda for ", colnames(dat)[i], "!!!"))
			return(1)
		})
	} else {
		lambda <- 1
	}

	# Update lambda in lambda table
	lambda_dat[lambda_dat$Trait == colnames(dat)[i],"Lambda"] <- lambda

	# Transform data
	if (lambda != 0) {
		transformed_dat[,i] <- transformed_dat[,i] ^ lambda
	} else {
		transformed_dat[,i] <- log(transformed_dat[,i])
	}
}

# Update data
dat <- transformed_dat


##################################################
# Generate BLUE
##################################################

# Generate BLUE for each trait
blue_dat <- NULL
for(i in target_column_indexes){
	# Clean miscellaneous data
	dat[is.infinite(dat[,i]),i] = NA
	dat[is.nan(dat[,i]),i] = NA
	dat[is.null(dat[,i]),i] = NA
	dat <- dat[rowSums(is.na(dat)) != ncol(dat),]
	dat <- dat[,colSums(is.na(dat)) != nrow(dat)]

	# Create lmer model
	lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

	# estimate BLUE
	blue_model <- lme4::fixef(lme)

	# extract BLUE and add grand mean when only one repetition present
	blue_model_out <- blue_model
	blue_model_out[2:length(blue_model_out)] = blue_model_out[2:length(blue_model_out)] + summary(lme)$coefficients[1]
	blue_model_out_dat <- as.data.frame(blue_model_out, check.names=FALSE, stringsAsFactors=FALSE)
	colnames(blue_model_out_dat) <- colnames(dat)[i]
	blue_model_out_dat <- tibble::rownames_to_column(blue_model_out_dat, var = "Accession")
	blue_model_out_dat[,1] <- sub(colnames(dat)[feature_column_indexes[1]], "", blue_model_out_dat[,1])
	existing_original_accessions <- unique(sort(dat[!is.na(dat[[i]]),colnames(dat)[feature_column_indexes[1]]]))
	blue_model_out_dat[1,1] <- existing_original_accessions[!(existing_original_accessions %in% blue_model_out_dat[,1])][1]

	# Merge BLUE data
	if (is.null(blue_dat)) {
		blue_dat <- blue_model_out_dat
	} else {
		blue_dat <- dplyr::full_join(blue_dat, blue_model_out_dat, by = "Accession")
	}
}

# Sort BLUE data
blue_dat <- arrange(blue_dat, Accession)
blue_dat <- as.data.frame(blue_dat, check.names=FALSE, stringsAsFactors=FALSE)


##################################################
# Write output
##################################################

file_prefix <- gsub("\\.[^.]*$", "", basename(input_file))

write.table(
	x = outlier_removed_dat,
	sep = "\t",
	file = file.path(output_path, paste0(file_prefix, "_outlier_removed.txt")),
	na = "",
	row.names = FALSE,
	quote = FALSE
)

write.table(
	x = outlier_dat,
	sep = "\t",
	file = file.path(output_path, paste0(file_prefix, "_outlier.txt")),
	na = "",
	row.names = FALSE,
	quote = FALSE
)

write.table(
	x = lambda_dat,
	sep = "\t",
	file = file.path(output_path, paste0(file_prefix, "_lambda.txt")),
	na = "",
	row.names = FALSE,
	quote = FALSE
)

write.table(
	x = transformed_dat,
	sep = "\t",
	file = file.path(output_path, paste0(file_prefix, "_box_cox_transformed.txt")),
	na = "",
	row.names = FALSE,
	quote = FALSE
)

write.table(
	x = blue_dat,
	sep = "\t",
	file = file.path(output_path, paste0(file_prefix, "_BLUE.txt")),
	na = "",
	row.names = FALSE,
	quote = FALSE
)
