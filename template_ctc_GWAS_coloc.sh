#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="Coloc"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

# Load Java and other relevant packages for running Nextflow
module load Java/11.0.16

set -f

nextflow_path=[path to Nextflow executable]

eqtl_finemap=[path to parsed eQTL fine-mapping results]
ref=[path to eQTLGen variant reference: 1000G-30x_index.parquet]

ld=[LD panel]
gwas_finemap=[folder with fine-mapped GWAS loci]
gwas_info_file=[File with GWAS information]

gtf=[ENSEMBL .gtf file]
output_folder=[Output folder for colocalisation]

NXF_VER=23.04.4 ${nextflow_path}/nextflow run main.nf \
--finemapped_eqtl ${eqtl_finemap} \
--finemapped_gwas ${gwas_finemap} \
--allele_info ${ref} \
--gwas_info_file ${gwas_info_file} \
--ld_reference ${ld} \
--gtf ${gtf} \
--coloc_threshold 0.8 \
--full_output false \
--mrlink2 true \
--OutputDir ${output_folder} \
-profile singularity,slurm \
-resume
