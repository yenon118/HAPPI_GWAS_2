# HAPPI_GWAS_2

<!-- badges: start -->
<!-- badges: end -->

The HAPPI_GWAS_2 is a pipeline built for genome-wide association study (GWAS).

## Requirements

In order to run the HAPPI_GWAS_2, users need to install Miniconda and prepare the Miniconda environment in their computing systems.

Miniconda can be downloaded from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html).

Installation of the Miniconda is required, and Miniconda environment needs to be activated every time before running the HAPPI_GWAS_2 pipeline.

Write a Conda configuration file (.condarc) before creating a Conda environment:

```
nano ~/.condarc
```

Put the following text into the Conda configuration file (make sure you change *envs_dirs* and *pkgs_dirs*) then save the file:

```
envs_dirs:
	- /new/path/to/miniconda/envs
pkgs_dirs:
	- /new/path/to/miniconda/pkgs
channels:
	- bioconda
	- conda-forge
	- defaults
```

Create a Conda environment named *happigwas*:

```
conda create -n happigwas openjdk=8.0.192 vcftools htslib pandas snakemake \
r-devtools r-biocmanager r-argparse r-dplyr r-tidyr r-tibble r-stringr r-ggplot2 \
r-bh r-mvtnorm r-viridislite r-stringi r-rcpp r-uuid r-nlme r-digest r-matrix r-ape \
r-bigmemory r-genetics r-gplots r-htmltools r-lattice r-magrittr r-lme4 r-mass \
bioconductor-multtest r-plotly r-rcpparmadillo r-rgl r-gridextra r-scatterplot3d \
r-snowfall bioconductor-snpstats r-biganalytics r-biglm r-car r-foreach r-doparallel
```

Activate *happigwas* Conda environment:

```
conda activate happigwas
```

Start R in terminal:

```
R
```

Install required R packages (Do not update any packages if any messages with multiple choices pop-up):

```
install.packages("EMMREML", repos = "https://cloud.r-project.org/")
devtools::install_github('christophergandrud/DataCombine', force=TRUE)
devtools::install_github("SFUStatgen/LDheatmap", force=TRUE)
devtools::install_github("jiabowang/GAPIT", force=TRUE)
```

Quit R:

```
q()
```

## Installation

You can install the HAPPI_GWAS_2 from [Github](https://github.com/yenon118/HAPPI_GWAS_2.git) with:

```
git clone https://github.com/yenon118/HAPPI_GWAS_2.git
```

## Usage

The HAPPI_GWAS_2 pipeline is a command line based pipeline that can be ran on any Linux computing systems. It consists of BLUP.py for best linear unbiased prediction, BLUE.py for best linear unbiased estimation, and HAPPI_GWAS.py for GWAS, haploblock analysis, and condidate gene identification. The command and arguments of each tool are shown as below:

#### BLUP.py

```
usage: python BLUP.py [-h] -p PROJECT_NAME -w WORKFLOW_PATH -i INPUT_FOLDER -o OUTPUT_FOLDER [-e FEATURE_COLUMN_INDEXES]
						[--ulimit ULIMIT] [--memory MEMORY] [--threads THREADS]
						[--keep_going] [--jobs JOBS] [--latency_wait LATENCY_WAIT] [--cluster CLUSTER]

mandatory arguments:
  -p PROJECT_NAME, --project_name PROJECT_NAME
                        Project name
  -w WORKFLOW_PATH, --workflow_path WORKFLOW_PATH
                        Workflow path
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Input folder
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder

optional arguments:
  -h, --help            show this help message and exit
  -e FEATURE_COLUMN_INDEXES, --feature_column_indexes FEATURE_COLUMN_INDEXES
                        Feature column indexes
  --ulimit ULIMIT       Ulimit
  --memory MEMORY       Memory
  --threads THREADS     Threads
  --keep_going          Keep going
  --jobs JOBS           Jobs
  --latency_wait LATENCY_WAIT
                        Latency wait
  --cluster CLUSTER     Cluster parameters
```

#### BLUE.py

```
usage: python BLUE.py [-h] -p PROJECT_NAME -w WORKFLOW_PATH -i INPUT_FOLDER -o OUTPUT_FOLDER [-e FEATURE_COLUMN_INDEXES]
						[--ulimit ULIMIT] [--memory MEMORY] [--threads THREADS]
						[--keep_going] [--jobs JOBS] [--latency_wait LATENCY_WAIT] [--cluster CLUSTER]

mandatory arguments:
  -p PROJECT_NAME, --project_name PROJECT_NAME
                        Project name
  -w WORKFLOW_PATH, --workflow_path WORKFLOW_PATH
                        Workflow path
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Input folder
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder

optional arguments:
  -h, --help            show this help message and exit
  -e FEATURE_COLUMN_INDEXES, --feature_column_indexes FEATURE_COLUMN_INDEXES
                        Feature column indexes
  --ulimit ULIMIT       Ulimit
  --memory MEMORY       Memory
  --threads THREADS     Threads
  --keep_going          Keep going
  --jobs JOBS           Jobs
  --latency_wait LATENCY_WAIT
                        Latency wait
  --cluster CLUSTER     Cluster parameters
```

#### HAPPI_GWAS.py

```
usage: python3 HAPPI_GWAS.py [-h] -p PROJECT_NAME -w WORKFLOW_PATH -i INPUT_FOLDER -o OUTPUT_FOLDER -v VCF_FILE -g GFF_FILE [--gff_category GFF_CATEGORY] [--gff_key GFF_KEY]
								[--genotype_hapmap GENOTYPE_HAPMAP] [--genotype_data GENOTYPE_DATA] [--genotype_map GENOTYPE_MAP]
								[--kinship KINSHIP] [--z_matrix Z_MATRIX] [--corvariance_matrix CORVARIANCE_MATRIX]
								[--snp_maf SNP_MAF] [--model MODEL] [--pca_total PCA_TOTAL]
								[--ulimit ULIMIT] [--memory MEMORY] [--threads THREADS]
								[--keep_going] [--jobs JOBS] [--latency_wait LATENCY_WAIT] [--cluster CLUSTER]
								[--p_value_filter P_VALUE_FILTER] [--fdr_corrected_p_value_filter FDR_CORRECTED_P_VALUE_FILTER] [--ld_length LD_LENGTH]

mandatory arguments:
  -p PROJECT_NAME, --project_name PROJECT_NAME
                        Project name
  -w WORKFLOW_PATH, --workflow_path WORKFLOW_PATH
                        Workflow path
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Input folder
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder
  -v VCF_FILE, --vcf_file VCF_FILE
                        VCF file
  -g GFF_FILE, --gff_file GFF_FILE
                        GFF file

optional arguments:
  -h, --help            show this help message and exit
  --gff_category GFF_CATEGORY
                        GFF category
  --gff_key GFF_KEY     GFF key
  --genotype_hapmap GENOTYPE_HAPMAP
                        Genotype hapmap
  --genotype_data GENOTYPE_DATA
                        Genotype data
  --genotype_map GENOTYPE_MAP
                        Genotype map
  --kinship KINSHIP     Kinship matrix file
  --z_matrix Z_MATRIX   Z matrix file
  --corvariance_matrix CORVARIANCE_MATRIX
                        Corvariance matrix file
  --snp_maf SNP_MAF     SNP minor allele frequency
  --model MODEL         Model
  --pca_total PCA_TOTAL
                        Total PCA
  --ulimit ULIMIT       Ulimit
  --memory MEMORY       Memory
  --threads THREADS     Threads
  --keep_going          Keep going
  --jobs JOBS           Jobs
  --latency_wait LATENCY_WAIT
                        Latency wait
  --cluster CLUSTER     Cluster parameters
  --p_value_filter P_VALUE_FILTER
                        P-value filter
  --fdr_corrected_p_value_filter FDR_CORRECTED_P_VALUE_FILTER
                        FDR corrected p-value filter
  --ld_length LD_LENGTH
                        LD length
```

## Examples

These are a few basic examples which show you how to use the HAPPI_GWAS_2:

#### BLUP.py

```
cd /path/to/HAPPI_GWAS_2

conda activate happigwas

python BLUP.py -p Test -w /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/output/BLUP_Arabidopsis360
```

#### BLUE.py

```
cd /path/to/HAPPI_GWAS_2

conda activate happigwas

python BLUE.py -p Test -w /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/output/BLUE_Arabidopsis360
```

#### HAPPI_GWAS.py

```
cd /path/to/HAPPI_GWAS_2

conda activate happigwas

python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLM \
-v /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--p_value_filter 0.01
```

```
cd /path/to/HAPPI_GWAS_2

conda activate happigwas

python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM \
-v /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--p_value_filter 0.01
```

```
cd /path/to/HAPPI_GWAS_2

conda activate happigwas

python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_FarmCPU \
-v /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--p_value_filter 0.01 \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G --output=log_2023_06_15_r_gapit_\%A-\%a.out"
```

```
cd /path/to/HAPPI_GWAS_2

sbatch test_HAPPI_GWAS_MLMM_Arabidopsis1001_Chr1.sbatch
```

## Remarks

1. The execution time of the HAPPI_GWAS_2 pipeline mainly depends on the size of the data and the available computing resources on the machine.


