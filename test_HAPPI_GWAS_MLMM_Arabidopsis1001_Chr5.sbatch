#!/bin/sh

#SBATCH --account=xulab
#SBATCH --partition=Lewis,BioCompute,hpc5,General
#SBATCH --time=1-21:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=33G
#SBATCH --job-name=2023_07_06_run_HAPPI_GWAS_MLMM_Arabidopsis1001_Chr5
#SBATCH --output=log_2023_07_06_run_HAPPI_GWAS_MLMM_Arabidopsis1001_Chr5_%A-%a.out


cd $(pwd)


source activate happigwas


python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM_Arabidopsis1001_Chr5 \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/vcf/Arabidopsis1135_5.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/gff/Athaliana_TAIR10.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/genotype_hapmap/Chr5.hmp.txt \
--model MLMM \
--p_value_filter 1e-5 \
--memory 90 \
--keep_going \
--jobs 20 \
--latency_wait 300 \
--cluster "sbatch --account=xulab --nodes=1 --ntasks=1 --cpus-per-task=3 --time=1-21:00 --partition=Lewis,BioCompute,hpc5,General --mem=100G"
