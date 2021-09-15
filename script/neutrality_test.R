

require(data.table)
require(PopGenome)

args <- commandArgs(trailingOnly = TRUE)
print(args)
input.vcf <- as.character(args[1])
gtFile <- as.character(args[2])
pathToOutput <- as.character(args[3])

if(!file.exists(input.vcf))
    stop(input.vcf,": No such file or directory")
if(!file.exists(gtFile))
    stop(gtFile,": No such file or directory")
if(!dir.exists(pathToOutput))
    dir.create(pathToOutput)

Filename <- gsub(".maf4.vcf.gz", "", basename(input.vcf))
tempDir <- paste0(pathToOutput, "/", Filename)
system(paste0("mkdir -p ", tempDir))

## Load Plasmodium falciparum GFF file
geneFile <- fread(gtFile, header = FALSE)
inc <- 1

fuli <- NULL

for (i in 1:nrow(geneFile)) 
{
    cat("\n")
    cat("i = ", i)
    cat("\n")
    cat("Extract ", geneFile$V6[i], " from ", input.vcf)
    cat("\n")
    File <- paste0(pathToOutput, "/", Filename, '_tajima.xlsx')
    region <- paste0(geneFile$V1[i], ":", geneFile$V4[i], "-", geneFile$V5[i])
    output <- paste0(tempDir, '/', geneFile$V6[i], ".vcf.gz")
    system(paste0("bcftools view -r ", region," -Oz -o ", output, " ", input.vcf))
    system(paste0("tabix ", output))
    
    ## Check if gene not NULL
    skipNum <- as.numeric(system(paste("zgrep \"##\"", output," | wc -l"), TRUE))
    vcf  <- read.table(output, 
                       skip = skipNum, 
                       header = TRUE, 
                       comment.char = "", 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)
    NbreSNPs = nrow(vcf)
    
    if(NbreSNPs >= 3)
    {
        GENOME.class <- readVCF(output, NbreSNPs, geneFile$V1[i], geneFile$V4[i], geneFile$V5[i], include.unknown = TRUE)
        Stats <- neutrality.stats(GENOME.class, FAST=TRUE)
        
        if(inc == 1)
        {
            geneName <- gsub(".vcf.gz", "", basename(output))
            ff <- as.data.frame(subset(get.neutrality(Stats)[[1]], select=-c(3,6:9)))
            fuli <- cbind(geneName, geneFile$V1[i], rownames(ff), ff)
            
            cat("Save output as ", File)
            system(paste0("touch ", File))
            write.table(fuli, File, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
        }
        else
        {
            geneName <- gsub(".vcf.gz", "", basename(output))
            ff <- as.data.frame(subset(get.neutrality(Stats)[[1]], select=-c(3,6:9)))
            fuli <- cbind(geneName, geneFile$V1[i], rownames(ff), ff)
            
            # Save file
            write.table(fuli, File, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE, sep = '\t')
        }
        
        inc <- inc + 1
        file.remove(output, paste0(output, ".tbi"))
        
    }
    else
    {
        cat('\n')
        file.remove(output, paste0(output, ".tbi"))
    }
}

system(paste0("rm -rf ", tempDir))