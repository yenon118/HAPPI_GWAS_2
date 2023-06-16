
Rscript GAPIT.R \
-i ../../data/mdp_traits.txt \
-o ../../output/GLM \
--genotype_data ../../data/mdp_numeric.txt \
--genotype_map ../../data/mdp_SNP_information.txt \
--model GLM


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_Output \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_Output \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G"


python3 HAPPI_GWAS_2.py \
-p Test \
-w /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2 \
-i /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/raw_data_split \
-o /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_Output \
-v /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--genotype_hapmap /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--cluster "sbatch --account=xulab --cpus-per-task=3 --time=0-02:00 --partition=Lewis,BioCompute,hpc5,General --mem=64G --output=log_2023_06_15_r_gapit_\%A-\%a.out"


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --recode --stdout | bgzip > \
/storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_Output/VCFtools/VCFtools_plink/9_11331274_11351274.recode.vcf.gz


vcftools --gzvcf /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/data/vcf/mdp_genotype_test.vcf.gz \
--chr 9 --from-bp 11331274 --to-bp 11351274 --plink \
--out /storage/htc/joshilab/yenc/projects/2022_07_22_RuthieAngelovici/HAPPI_GWAS_2/output/HAPPI_GWAS_Output/VCFtools/VCFtools_plink/9_11331274_11351274


snakemake -p -n --snakefile rules/vcftools_subset_and_convert_to_plink.smk \
--configfile output/HAPPI_GWAS_Output/VCFtools/Snakemake_VCFtools_config.json


java -jar tools/Haploview.jar \
-n -out output/HAPPI_GWAS_Output/VCFtools/9_11331274_11351274 \
-pedfile output/HAPPI_GWAS_Output/VCFtools/VCFtools_plink/9_11331274_11351274.ped \
-map output/HAPPI_GWAS_Output/VCFtools/VCFtools_plink/9_11331274_11351274.map \
-skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05
