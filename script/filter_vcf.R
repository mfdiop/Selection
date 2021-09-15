
arguments <- commandArgs(trailingOnly = TRUE)

Dir <- arguments[1]

library(data.table)
library(tictoc)
library(tidyverse)
options(scipen=999)

#-------- computing missingness on individuals
computeMissingness = function(dataFrame)
{
    dataFrame = t(dataFrame)
    missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
    for(i in 1:dim(dataFrame)[1])
    {
        m=0
        for(j in 1:dim(dataFrame)[2])
        {
            if(dataFrame[i,j] == './.')
                m=m+1
        }
        missingness[i] = m/(dim(dataFrame)[2])
    }
    return(missingness)
}

#-------- computing missingness on SNPs
computeMissingnessOnSNPs = function(dataFrame)
{
    dataFrame = as.matrix(dataFrame)
    missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
    for(i in 1:dim(dataFrame)[1])
    {
        m=0
        for(j in 1:dim(dataFrame)[2])
        {
            if(dataFrame[i,j] == './.')
                m=m+1
        }
        missingness[i] = m/(dim(dataFrame)[2]) #-1
    }
    return(missingness)
}

# Calculates the amount of missing data in a row or column where missing data is defined as "./."

# missi <- function(x)
# {
#     sum(1*(x == "./."))/ length(x)
# }

## TRIQUAD ##
# Calculates the number of different homozygous calls present in the data (1 - 4). 
# Take heterozygous calls into account. 
# Also does not take into account the frequency of each allele.

triquad <- function(x)
{
    xx <- x[x != "./."];
    res <- 1*(sum(1*(xx=="0/0"))>0) + 1*(sum(1*(xx=="1/1"))>0) + 1*(sum(1*(xx=="0/1"))>0);
    res
}

#------------ Extracting genotypes from the raw data:
setwd(Dir)

vcf <- arguments[2]
maf <- arguments[3]

Filename <- gsub(".vcf.gz", "", basename(vcf))
Genotypes <- paste0('../Data/processed/', Filename, '.Genotypes.txt')

expression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, Genotypes))

#------- Computing the missingness on SNPs and samples
#------- Get the sample names from the VCF file

sampleList <- paste0('../Data/processed/', Filename, '.SampleList.txt')
system(sprintf("bcftools query -l %s > %s", vcf, sampleList))
AfricanSamplesList <- fread(sampleList, header = FALSE)

#---- put the first four columns in a variable
Genotype <- fread(Genotypes, header = FALSE)

first4Column <- subset(Genotype, select=c(1:4))
Genotype <- subset(Genotype, select=-c(1:4))

#===================
## Remove invariants
#==================
varsnp <- apply(Genotype, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, Genotype)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#======================================
#--- computing missingness on SNPs
#======================================

tic();snpMissingness = computeMissingnessOnSNPs(geno);toc()

geno <- cbind(first4Column, geno, snpMissingness)

geno <- geno[which(geno$snpMissingness <= 0.2), ]
first4Column <- subset(geno, select = c(1:4))
geno <- subset(geno, select = -c(1:4, ncol(geno)))

#================================
#=== Compute Missingness on samples
#================================
tic();sampleMissingness = computeMissingness(geno);toc()

geno <- t(geno)
rownames(geno) <- AfricanSamplesList$V1
geno <- as.data.frame(geno)
geno$Missingness <- sampleMissingness

index <- which(geno$Missingness > 0.15)

geno <- subset(geno, select = -c(ncol(geno)))
geno <- as.data.frame(t(geno))

if(!is_empty(index)) geno <- geno[, -index]

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

## Save positions and isolates to keep for downstream analysis
write.table(colnames(geno), paste0("../Data/processed/", Filename, ".samplesTokeep.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(first4Column[,1:2], paste0("../Data/processed/", Filename, ".snpsTokeep.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

#=========================================================
#---- Removing SNPs and Isolates to discard from vcf file
#========================================================
samplesToKeep = paste0("../Data/processed/", Filename, ".samplesTokeep.txt")
snpsToKeep = paste0("../Data/processed/", Filename, ".snpsTokeep.txt")

filtered.vcf <- paste0("../Data/processed/", gsub(".vcf.gz", ".filtered", basename(vcf)))
system(paste0("vcftools --gzvcf ", vcf,
              " --keep ", samplesToKeep, 
              " --positions ", snpsToKeep,
              " --not-chr Pf3D7_API_v3", 
              " --recode --recode-INFO-all --out ", filtered.vcf))

## Rename filtered VCF file
system(paste0("mv ", filtered.vcf, ".recode.vcf", " ", filtered.vcf, ".vcf"))

# Remove temporally data
file.remove(Genotypes, samplesToKeep, snpsToKeep, sampleList)

## Compress and index filtered vcf file
system(paste0("bgzip ", filtered.vcf, ".vcf"))
system(paste0("tabix ", filtered.vcf, ".vcf.gz"))

## Remove sites with low allele frequency
filtered.vcf <- paste0(filtered.vcf, ".vcf.gz")
maf.vcf <- paste0("../Data/reference/", gsub(".filtered", paste0(".maf", maf), basename(filtered.vcf)))

system(paste0("bcftools view -q ", maf, ":minor -Oz -o ", maf.vcf, " ", filtered.vcf))