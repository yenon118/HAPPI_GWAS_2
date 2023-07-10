
python BLUP.py -p Test -w /scratch/yenc/projects/HAPPI_GWAS_2 \
-i /scratch/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /scratch/yenc/projects/HAPPI_GWAS_2/output/BLUP_Arabidopsis360


python BLUE.py -p Test -w /scratch/yenc/projects/HAPPI_GWAS_2 \
-i /scratch/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /scratch/yenc/projects/HAPPI_GWAS_2/output/BLUE_Arabidopsis360


python3 HAPPI_GWAS.py \
-p Test \
-w /scratch/yenc/projects/HAPPI_GWAS_2 \
-i /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /scratch/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLM \
-v /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /scratch/yenc/projects/HAPPI_GWAS_2 \
-i /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /scratch/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM \
-v /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /scratch/yenc/projects/HAPPI_GWAS_2 \
-i /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /scratch/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_FarmCPU \
-v /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /scratch/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--p_value_filter 0.01
