#=============================================================================
## This script is to create / reconstruct gene fasta sequence from vcf file
## by phasing mixed calls and imputing missing genotypes
## It will generate the full sequence of the gene for each sample from the vcf
## Dependencies: Imputed data / Output directory
#=============================================================================
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    cat("||======================================\n")
    cat("|| Usage: Rscript vcf_to_fasta_alignment.R \n||\timputed.txt \n||\treference.fasta\n")
    cat("||======================================\n")
    stop("Not enough arguments supplied ...\n")
}

library(tidyverse)
source("vcf_to_fasta_alignment_functions.R")
#============================================
#==== Compare sample sequences to reference
#===========================================
input <- args[1]
input <- gsub("\"", "", input)
ref <- args[2]
ref <- gsub("\"", "", ref)

#==============
# Converting
#==============
cat("Read ", basename(input), "imputed data \n")
Data <- read_tsv(input)
firstColumns <- Data %>% select(c(1:6))
cat("\n")
cat("Convert ", basename(input), " to haplotype format\n")
cat("\n")
hapFormat <- genotype_to_haplotype(Data)

#=====================
# FORMAT THE POSITION
#=====================
firstColumns$AMINO_ACID <- as.numeric(gsub("[aA-zZ-\\*]*", "", firstColumns$AMINO_ACID))

for (i in 1:nrow(firstColumns)) {
    cat("Processing row ", i)
    cat("\n")
    amino.string <- unlist(strsplit(firstColumns$CODON[i], "/"))[1]
    amino.length <- sum(nchar(amino.string))
    amino.matrix <- matrix(nrow = 1, ncol = amino.length)
    
    for (j in 1:amino.length) {
        amino.matrix[1, j] <- substr(amino.string, j, j)
    }
    
    index <- which(grepl("^[[:upper:]]+$", amino.matrix))
    
    if (index == 1)  firstColumns$POS[i] <- firstColumns$AMINO_ACID[i] * 3 - 2
    else if (index == 2)  firstColumns$POS[i] <- firstColumns$AMINO_ACID[i] * 3 - 1
    else  firstColumns$POS[i] <- firstColumns$AMINO_ACID[i] * 3
}

#====================================================
vcf.name <- gsub("_[\\.aA-zZ]*", "", basename(input))

ref.string <- toupper(as.matrix(paste(readLines(ref)[-1], collapse = "")))
seq.length <- sum(nchar(ref.string))
ref.matrix <- matrix(nrow = 1, ncol = seq.length)

for (i in 1:seq.length) {
    ref.matrix[1, i] <- substr(ref.string, i, i)
}

cat("\n")
my.vcf.2.fasta <- matrix()

# Save the reference sequence
my.vcf.2.fasta <- rbind(my.vcf.2.fasta,
                        paste(">", vcf.name),
                        paste(ref.matrix, collapse = ""))

for (i in 1:ncol(hapFormat)) {
    
    seq.matrix <- ref.matrix

    cat("Reading Column ", i, "\n")
    sample_seqs <- hapFormat %>% select(i) %>% pull()

    for (j in 1:nrow(firstColumns)) {
        if(seq.matrix[1, as.numeric(as.character(firstColumns$POS[j]))] != sample_seqs[j])
            seq.matrix[1, as.numeric(as.character(firstColumns$POS[j]))] <- as.character(sample_seqs[j])
        else
            next
    }

    my.vcf.2.fasta <- rbind(my.vcf.2.fasta,
                            paste(">", colnames(hapFormat[i])),
                            paste(seq.matrix, collapse = ""))
    cat("Reading Sample", i, "of", ncol(hapFormat), "\n\r")
}

my.vcf.2.fasta <- as.matrix(my.vcf.2.fasta[-1, ])
cat("\n")
cat("Finished!\n")

write.table(my.vcf.2.fasta,
            file = file.path(dirname(input), paste(vcf.name, ".fasta", sep = "")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("Saved as ", vcf.name, ".fasta in ", dirname(input), sep = "")
cat("\n")
