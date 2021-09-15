#!/usr/bin/bash

#==========================================
# This script splits VCF file by chromosome
#==========================================

module load tabix

VCF=$1
DIRPATH=$2

# Create output directory
mkdir -p $DIRPATH

# get basename and add .gz compression extension
VCFGZ="${VCF##*/}.gz"

#compress vcf
bgzip -c $VCF > ${DIRPATH}/${VCFGZ}

# index compressed vcf
tabix -p vcf ${DIRPATH}/${VCFGZ}

# save all the chromosome names into a file
tabix --list-chroms $VCFGZ > ${DIRPATH}/chromosomes.txt

# make an individual vcf for each chromosome
while IFS= read -r line; do
  tabix $VCFGZ $line > ${DIRPATH}/${line}.vcf;
done < ${DIRPATH}/chromosomes.txt

# NB = Also, you can use that --list-chroms option just to find out the names of your chromosomes.