#!/bin/sh

#SBATCH -A mdiop
#SBATCH -J concatenate
#SBATCH -p standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --error="concatenate.error"
#SBATCH --output="concatenate.out"    
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=mdiop@mrc.gm

module load bcftools

Path=$1

cd $Path

vcf_files=`ls -p *.vcf.gz | grep -v / | tr '\n' ' '`
output=basename $(pwd)

bcftools concat -Oz -o ${output} ${vcf_files}
