library(data.table)
library(tictoc)
library(statip)
options(scipen=999)

# Calculates the amount of missing data in a row or column where missing data is defined as a -, N, NA

missi <- function(x)
{
    sum(1*(x == "./."))/ length(x)
}

## TRIQUAD ##

# Calculates the number of different homozygous calls present in the data (1 - 4). 
# Does not take heterozygous calls into account as these can be excluded later via the mixes call. 
# Also does not take into account the frequency of each allele.

triquad <- function(x)
{
    xx <- x[x != "./."];
    res <- 1*(sum(1*(xx=="0/0"))>0) + 1*(sum(1*(xx=="1/1"))>0) + 1*(sum(1*(xx=="0/1"))>0);
    res
}


#------------ Extracting genotypes from the raw data:
setwd("/media/Data/Data/Documents_Karim/Fadel/Alfred/Data/raw")
vcf <- "Farafenni_1990.vcf.gz"
Genotypes = '../filtered/Genotypes.txt'

expression = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, Genotypes))

#------- Computing the missingness on SNPs and samples
#---- put the first four columns in a variable
Genotype = fread(Genotypes, header = FALSE)

first4Column = subset(Genotype, select=c(1:4))
Genotype = subset(Genotype, select=-c(1:4))

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

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp <= 0.8)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#=======================================
# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

## Remove isos with too much missingness
isomissallow <- 0.2
geno <- as.data.frame(geno)
geno <- geno[,which(missiso < isomissallow)]
nrow(geno)
ncol(geno)

############# ADDED
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

# Determine missingness on SNPs/ Isolates
#======================================
tic();misssnp <- apply(geno, 1, missi);toc()
tic();missiso <- apply(geno, 2, missi);toc()
par(mfrow=c(2,1))
plot(missiso, pch=16)
plot(misssnp, pch='.')

#=====================================
# Remove SNPs with too much missingness
#=====================================
keepsnp <- which(misssnp <= 0.2)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepsnp,]

geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)
