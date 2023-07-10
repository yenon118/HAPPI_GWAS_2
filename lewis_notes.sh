
Rscript BLUP.R -i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data/05_22_2019_Arabidopsis_360_BCAA_raw.csv \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/BLUP_Arabidopsis360 \
-e 1,2


Rscript BLUE.R -i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data/05_22_2019_Arabidopsis_360_BCAA_raw.csv \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/BLUE_Arabidopsis360 \
-e 1,2


Rscript GAPIT.R \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/raw_data/mdp_traits.txt \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_GLM/GAPIT/GLM \
--genotype_data /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/genotype_data/mdp_numeric.txt \
--genotype_map /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/genotype_map/mdp_SNP_information.txt \
--model GLM


python BLUP.py -p Test -w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/BLUP_Arabidopsis360


python BLUE.py -p Test -w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/BLUE_Arabidopsis360


python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLM \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--p_value_filter 0.01 \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G"


python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_FarmCPU \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--p_value_filter 0.01 \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G --output=log_2023_06_15_r_gapit_\%A-\%a.out"


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --recode --stdout | bgzip > \
/storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9__11331274__11351274.recode.vcf.gz


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --plink \
--out /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9__11331274__11351274


snakemake -p -n --snakefile rules/vcftools_haploview.smk \
--configfile output/HAPPI_GWAS_MLMM/VCFtools_Haploview/Snakemake_VCFtools_Haploview_config.json


java -jar tools/Haploview.jar \
-n -out output/HAPPI_GWAS_MLMM/VCFtools_Haploview/Haploview/9__11331274__11351274 \
-pedfile output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9__11331274__11351274.ped \
-map output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9__11331274__11351274.map \
-skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05


python3 HAPPI_GWAS.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM_Arabidopsis1001_Chr1 \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/vcf/Arabidopsis1135_1.vcf.gz \
-g /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/gff/Athaliana_TAIR10.gff3 \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/genotype_hapmap/Chr1.hmp.txt \
--model MLMM \
--p_value_filter 1e-5 \
--memory 90 \
--jobs 10 \
--cluster "sbatch --account=xulab --nodes=1 --ntasks=1 --cpus-per-task=3 --time=1-21:00 --partition=Lewis,BioCompute,hpc5,General --mem=100G"
