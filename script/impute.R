
args = (commandArgs(TRUE))

vcf <- args[1]
Dir <- args[2]

library(data.table)
library(tidyverse)
library(statip)

source("vcf_to_fasta_alignment_functions.R")

# vcf <- "/media/Data/Data/Documents_Karim/Fadel/BalancingSelection/results/3D_Structures/raw_data/PEP.vcf.gz"
# Dir <- "/media/Data/Data/Documents_Karim/Fadel/BalancingSelection/results/3D_Structures/Imputed"

if(!dir.exists(Dir)) dir.create(Dir)

cat("Extract the genotype and the read depth from ", gsub("_[\\.aA-zZ]*", "", basename(vcf)))
cat("\n")
Genotypes <- file.path(Dir, 'Genotypes.txt')
AllelicDeph <- file.path(Dir, 'AllelicD.txt')
GTexpression <- '%CHROM\t%POS\t%REF\t%ALT\t%SNPEFF_AMINO_ACID_CHANGE\t%SNPEFF_CODON_CHANGE[\t%GT]\n'
ADexpression <- '%CHROM\t%POS\t%REF\t%ALT\t%SNPEFF_AMINO_ACID_CHANGE\t%SNPEFF_CODON_CHANGE[\t%AD]\n'

system(paste0("bcftools query -l ", vcf, " > " , paste(Dir, "SampleIDs.txt", sep = "/")))
system(sprintf("bcftools query -f'%s' %s > %s", GTexpression, vcf, Genotypes))  
system(sprintf("bcftools query -f'%s' %s > %s", ADexpression, vcf, AllelicDeph))

#======================================
# ---- Phased the final vcf file:
#======================================

genotypeData <- fread(Genotypes, header = FALSE)
allelicDepthData <- fread(AllelicDeph, header = FALSE)

firstColumns <- genotypeData %>% filter(V5 != ".") %>% select(c(1:6)) %>% as_tibble() 
write.table(firstColumns, file.path(Dir, "firstColumns.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

genotypeData <- genotypeData %>% filter(V5 != ".") %>% select(-c(1:6)) %>% as_tibble() 
allelicDepthData <- allelicDepthData %>% filter(V5 != ".") %>% select(-c(1:6)) %>% as_tibble()

genotypeData <- as.matrix(genotypeData)
genotypeData[genotypeData == '0/0'] <- 0
genotypeData[genotypeData == '1/1'] <- 1
genotypeData[genotypeData == './.'] <- '.'
genotypeData <- as.data.frame(genotypeData)

Sys.time()
PhasedData <- genotypePhasing(genotypeData, 
                              allelicDepthData
                              )
Sys.time()
#============================================
#----------- Impute missing data ------------
#============================================
system(paste0('mkdir -p ', Dir, '/Imputation_Folder'))

numberOfSimulation <- 100 ; 
# Name <- paste0(gsub("\\.[aA-zZ]*", "", basename(vcf)), '_ImputedData.txt')
Name <- paste0(Dir, "/", mgsub::mgsub(basename(vcf), c("_[\\.aA-zZ]*", "\\.[aA-zZ]*"), c("", "")), 
               '_ImputedData.txt')

cat("#======== Starting the timer =========\n")
cat("\n")
outputDir <- paste0(Dir, '/Imputation_Folder')
Start <- Sys.time()
imputeMissingGenotypes(PhasedData, firstColumns, numberOfSimulation, outputDir, Name)
cat ("\n")
End <- Sys.time()
cat("Imputation took", round(End - Start, digits = 2), units.difftime(End - Start), "\n")

file.remove(Genotypes, AllelicDeph)
# delete a directory -- must add recursive = TRUE
unlink(outputDir, recursive = TRUE)
