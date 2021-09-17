
#======================================
#==========PHASING GENOTYPES===========
#======================================
genotypePhasing <-function (genotypeData, allelicDepth)
{
    phasedData = matrix(9, nrow=dim(genotypeData)[1], ncol=dim(genotypeData)[2])
    for(j in 1:nrow(genotypeData))
    {
        cat("Processing row ", j, " in Phasing process \n")
        k=1
        while(k<=ncol(genotypeData))
        {
            if(genotypeData[j,k] == 0)
                phasedData[j,k] <- 0   
            else if(genotypeData[j,k] == '.')
                phasedData[j,k] <- '.'
            else if(genotypeData[j,k] == 1)
                phasedData[j,k] <- 1
            else if(genotypeData[j,k] == '0/1')
            {
                target <- as.integer(unlist(strsplit(as.character(allelicDepthData[j,k]), ',')))
                if(sum(target)>0)
                {
                    if(target[1]<target[2])
                    {
                        minor <- target[1]
                        maf <- minor/sum(target)
                        phasedData[j,k] <- rbern(1,maf)
                    }
                    else
                    {
                        minor <- target[2] 
                        maf <- minor/sum(target)
                        phasedData[j,k] <- rbern(1,maf)
                    }
                }
                else
                {
                    maf <- target[1]
                    phasedData[j,k] <- rbern(1,maf)
                }
                
            }
            k <- k+1
        }
    }
    return(phasedData)
}

#======================================
#============= IMPUTATION =============
#======================================
imputeMissingGenotypes <- function(PhasedData, firstColumns, numberOfSimulation, outputDir, Name)
{
    maf = vector(mode = "numeric", length = dim(PhasedData)[1])
    PhasedData <- as.data.frame(PhasedData)
    print(paste0("computing the MAF before simulation"), quote = FALSE)
    
    for(i in 1:dim(PhasedData)[1])
    {
        print(paste0("i = ",i), quote = FALSE)
        count1 <- 0
        count2 <- 0
        for(j in 1:dim(PhasedData)[2])
        {
            if(PhasedData[i,j] == '.')
                count1 <- count1+1
            if(PhasedData[i,j] == '1')
                count2 <- count2+1
        }
        validCount = dim(PhasedData)[2] - count1
        count3 = validCount-count2
        if(count2 < count3)
            maf[i] <- count2/validCount
        else if(count2 > count3)
            maf[i] <- count3/validCount
        else
            maf[i] <- count2/validCount
    }
    
    print(paste0("Imputing the missing genotypes"), quote = FALSE)
    s=1
    correlationCoef = vector(mode = "numeric", length = numberOfSimulation)
    while(s <= numberOfSimulation)
    {
        imputedGenotypeMatrix = matrix(-9, nrow = dim(PhasedData)[1], ncol = dim(PhasedData)[2])
        for(j in 1:dim(PhasedData)[1])
        {
            k=1
            while(k<=dim(PhasedData)[2])
            {
                if(PhasedData[j,k] == '0')
                    imputedGenotypeMatrix[j,k] <- 0   
                else if(PhasedData[j,k] == '1')
                    imputedGenotypeMatrix[j,k] <- 1
                else if(PhasedData[j,k] == '.')
                    imputedGenotypeMatrix[j,k] <- rbern(1,maf[j])
                k=k+1
            }
        }
        
        print(paste0("computing the MAF after simulation"), quote = FALSE)
        mafAfter <- vector(mode = "numeric", length = dim(imputedGenotypeMatrix)[1])
        count4 <- rowSums(imputedGenotypeMatrix)
        count5 <- dim(imputedGenotypeMatrix)[2]-count4
        for(l in 1:length(mafAfter))
        {
            if(count4[l] < count5[l])
                mafAfter[l] = count4[l]/dim(imputedGenotypeMatrix)[2]
            else
                mafAfter[l] = count5[l]/dim(imputedGenotypeMatrix)[2]
            
        }
        correlationCoef[s] = cor(maf, mafAfter)
        
        imputedGenotypeMatrix <- cbind(firstColumns, as.data.frame(imputedGenotypeMatrix))
        outFileName <- paste0(outputDir, '/Imputed_', s, ".txt")
        write.table(imputedGenotypeMatrix, outFileName, quote = FALSE, row.names = FALSE, col.names = FALSE)
        print(paste0("simulation = ", s), quote = FALSE)
        s <- s+1
    }
    
    index <- which(correlationCoef == max(correlationCoef))
    cat('--- Save The Imputed File = ', index, '---')
    cat("\n")
    samples <- scan(paste0(dirname(outputDir), '/SampleIDs.txt'), 
                           what = 'character')
    Imputation <- fread(paste0(outputDir, '/Imputed_', index, ".txt"))
    names(Imputation) <- c("CHROM", "POS", "REF", "ALT", "SNPEFF_EFFECT", samples)
    write.table(Imputation, Name, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}


#==========================================================
#============= Convert Genotypes to Haplotype format ======
#==========================================================
genotype_to_haplotype <- function(genotype)
{
    library(tidyverse)
    # genotype <- as.data.frame(genotype)
    
    first4Columns <- genotype %>% select(c(1:6))
    genotype <- genotype %>% select(-c(1:6))
    header <- names(genotype)
    
    hapFormat <- data.frame(matrix(ncol = ncol(genotype), nrow = nrow(genotype)))
    
    for (i in 1:nrow(genotype)) 
    {
        cat(paste0("Processing SNPs ", i, "\n"))
        cat("\n")
        k <- 1
        for (j in 1:ncol(genotype)) 
        {
            if(genotype[i,j]== 0) hapFormat[i,k] <- first4Columns[i,3]
            else if(genotype[i,j]== 1) hapFormat[i,k] <- first4Columns[i,4]
            
            k <- k+1
        }
    }
    
    hapFormat <- as_tibble(hapFormat)
    names(hapFormat) <- header
    return(hapFormat)
}
