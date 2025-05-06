
python BLUP.py -p Test -w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/BLUP_Arabidopsis360


python BLUE.py -p Test -w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis360_example_data/original_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/BLUE_Arabidopsis360


python3 HAPPI_GWAS.py \
-p Test \
-w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLM \
-v /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM \
-v /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model MLMM \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/raw_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_FarmCPU \
-v /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/vcf/mdp_genotype_test.vcf.gz \
-g /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/gff/Zea_mays.AGPv3.26.gff3 \
--genotype_hapmap /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Maize_example_data/genotype_hapmap/mdp_genotype_test.hmp.txt \
--model FarmCPU \
--p_value_filter 0.01


python3 HAPPI_GWAS.py \
-p Test \
-w /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2 \
-i /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/raw_data_split \
-o /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/output/HAPPI_GWAS_MLMM_Arabidopsis1001_Chr1 \
-v /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/vcf/Arabidopsis1135_1.vcf.gz \
-g /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/gff/Athaliana_TAIR10.gff3 \
--genotype_hapmap /project/bilyeu_soybean_genomic_merge/yenc/projects/HAPPI_GWAS_2/data/Arabidopsis1001_example_data/genotype_hapmap/Chr1.hmp.txt \
--model MLMM \
--p_value_filter 1e-5 \
--memory 90 \
--keep_going \
--jobs 20 \
--latency_wait 300
