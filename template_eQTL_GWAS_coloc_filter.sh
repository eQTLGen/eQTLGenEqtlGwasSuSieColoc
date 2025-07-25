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

eqtl_gwas_coloc=[path to coloc results]
ref=[path to eQTLGen variant reference: 1000G-30x_index.parquet]

ld=[LD panel]

output_folder=[Output folder for colocalisation]

NXF_VER=23.04.4 ${nextflow_path}/nextflow run main.nf \
--coloc ${eqtl_gwas_coloc} \
--variant-reference ${ref} \
--ld_reference ${ld} \
--OutputDir ${output_folder} \
-profile singularity,slurm \
-resume
