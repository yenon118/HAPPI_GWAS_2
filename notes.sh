
Rscript GAPIT.R \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data/mdp_traits.txt \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_GLM \
--genotype_data /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_data/mdp_numeric.txt \
--genotype_map /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_map/mdp_SNP_information.txt \
--model GLM


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLM \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--p_value_filter 0.01


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--p_value_filter 0.01 \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G"


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_FarmCPU \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--p_value_filter 0.01 \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G --output=log_2023_06_15_r_gapit_\%A-\%a.out"


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --recode --stdout | bgzip > \
/storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9_11331274_11351274.recode.vcf.gz


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --plink \
--out /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9_11331274_11351274


snakemake -p -n --snakefile rules/vcftools_haploview.smk \
--configfile output/HAPPI_GWAS_MLMM/VCFtools_Haploview/Snakemake_VCFtools_Haploview_config.json


java -jar tools/Haploview.jar \
-n -out output/HAPPI_GWAS_MLMM/VCFtools_Haploview/Haploview/9_11331274_11351274 \
-pedfile output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9_11331274_11351274.ped \
-map output/HAPPI_GWAS_MLMM/VCFtools_Haploview/VCFtools/9_11331274_11351274.map \
-skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05
