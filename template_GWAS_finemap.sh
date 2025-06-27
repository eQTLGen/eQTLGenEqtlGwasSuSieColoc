#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="GWAS finemap"

# Load Java and other relevant packages for running Nextflow
module load Java/11.0.16

set -f

nextflow_path=[path to Nextflow executable]

signals=[combined_finemapping_zip]
ref=[path to eQTLGen variant reference: 1000G-30x_index.parquet]

ld=[LD panel]
gwas=[Folder with harmonised GWAS files]
gwas_info_file=[File with GWAS information]
af_file=[eQTLGen allele frequency file]

gtf=[ENSEMBL .gtf file]
output_folder=[Output folder for fine-mapping]

NXF_VER=23.04.4 ${nextflow_path}/nextflow run finemap_gwas.nf \
--finemapped_signals ${signals} \
--allele_info ${ref} \
--gwas_database ${gwas} \
--gwas_info_file ${gwas_info_file} \
--ld_reference ${ld} \
--gtf ${gtf} \
--coloc_threshold 0.8 \
--full_output false \
--OutputDir ${output_folder} \
-profile singularity,slurm \
-resume
